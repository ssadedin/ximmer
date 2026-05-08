
package ximmer

import com.xlson.groovycsv.PropertyMapper
import gngs.BED
import gngs.FASTA
import gngs.ProgressCounter
import gngs.Regions
import gngs.ToolBase
import gngs.XPos
import gngs.Region
import graxxia.TSV
import groovy.json.JsonSlurper
import groovy.transform.CompileStatic
import groovy.util.logging.Log
import htsjdk.variant.variantcontext.*
import htsjdk.variant.vcf.*

/**
 * Converts the CNV report produced by SummarizeCNVs into VCF format. If the form with copy number
 * information is supplied, includes that in the genotyping column.
 * 
 * <p>Rules for converting multiple variant types to single output allele:</p>
 * <ol>
 *   <li>If variant contains DEL or INV, output as DEL with copy number capped at 1</li>
 *   <li>If variant is pure DUP (no DEL/INV), output as DUP with original copy number</li>
 *   <li>For combined types (e.g. DUP,DEL or DUP,INV), DEL takes precedence</li>
 *   <li>Any variant output as DEL will have copy number capped at 1 regardless of input</li>
 * </ol>
 */
@Log
class TSVtoVCF extends ToolBase {
    
    static final String FILTER_LOW_CALLERS = 'LOW_CALLERS'
    static final String FILTER_FEW_TARGETS = 'FEW_TARGETS'
    
    static void main(String[] args) {
        cli('TSVtoVCF -i <CNV report> -s <sample> [-s <sample2> ...] -o <vcf output>', 'Converts Ximmer CNV report to VCF format', args) {
            i 'CNV report produced by Ximmer SummarizeCNVs', args: 1, required: true, type: File
            r 'Reference genome for computing reference alleles at CNV start positions', args:1, required: true, type:File
            hl 'Additional header lines to add', args:'*', type: String
            s 'Sample to include', args: '*', required: true, type: String
            t 'Target regions BED file for computing target overlap counts', args:1, required: false, type: File
            pass_targets 'Minimum number of overlapping target regions for PASS filter', args:1, required: false, type: Integer
            pass_caller_count 'Minimum number of callers for PASS filter', args:1, required: false, type: Integer
            source 'Optional Source to specify in header line', args:1, required: false, type: String
            o 'VCF output file', args:1, required: true, type: File
        }
    }

    @Override
    public void run() {
        FASTA ref = new FASTA(opts.r.absolutePath)
        List<Map> cnvReport
        if(opts.i.name.endsWith('.tsv'))
            cnvReport = new TSV(opts.i).toListMap()
        else
            cnvReport = new JsonSlurper().parse(opts.i)

        log.info "Converting $opts.i to VCF format"
        log.info "Using reference: $opts.r"
        
        Regions targetRegions = null
        if(opts.t) {
            log.info "Loading target regions from $opts.t"
            targetRegions = new BED(opts.t.absolutePath).load()
            log.info "Loaded ${targetRegions.numberOfRanges} target regions"
        }
        
        Integer passTargets = opts.pass_targets ? opts.pass_targets as Integer : null
        Integer passCallerCount = opts.pass_caller_count ? opts.pass_caller_count as Integer : null
        
        List<String> contigs = cnvReport*.chr.grep { !Region.isMinorContig(it) }.sort()
        
        opts.o.withWriter { w ->
            this.createVCF(w, ref, contigs, cnvReport, targetRegions, passTargets, passCallerCount)
        }
        
        log.info "Wrote $opts.o"
    }

    void createVCF(final Writer w, final FASTA genomeRef, final List<String> contigs, final List<Map> tsv,
                   final Regions targetRegions = null, final Integer passTargets = null, final Integer passCallerCount = null) {
        List<String> samples = opts.ss
        
        Set allHeaders = createHeaders(genomeRef, contigs, passTargets, passCallerCount)

        VCFHeader header = new VCFHeader(allHeaders, samples)

        writeHeader(w, header, opts.ss)
        
        VCFEncoder encoder = new VCFEncoder(header, true, true)
       
        ProgressCounter p = new ProgressCounter(withRate: true, withTime: true, log: log)
        
        List<VariantContext> outputVariants = new ArrayList(10000)
        
        for(Map line in tsv) {
            VariantContext vctx = createVariantFromLine(line, genomeRef, samples, targetRegions, passTargets, passCallerCount)
            if(vctx) {
                outputVariants.add(vctx)
            }
            p.count()
        }
        p.end()
        
        
        List<VariantContext> sortedVariants = outputVariants.sort { XPos.computePos(it.contig, it.start)}
        for(vctx in sortedVariants) {
                encoder.write(w, vctx)
                w.write('\n')
        }

        log.info "Sorting and writing ${outputVariants.size()} output variants ..."
        log.info "Wrote ${sortedVariants.size()} variants to VCF"
    }
    
    /**
     * Creates a VariantContext from a single line of TSV input
     * 
     * @param line The line from the TSV file containing variant information
     * @param genomeRef Reference genome for getting reference bases
     * @param samples List of samples to include in the output
     * @param targetRegions Optional target regions for computing overlap count
     * @param passTargets Minimum number of overlapping targets for PASS (null = no filter)
     * @param passCallerCount Minimum number of callers for PASS (null = no filter)
     * @return VariantContext object if valid, null if variant should be skipped
     */
    VariantContext createVariantFromLine(Map line, FASTA genomeRef, List<String> samples, 
                                         Regions targetRegions = null, Integer passTargets = null, 
                                         Integer passCallerCount = null) {
        String ref = genomeRef.basesAt(line.chr, line.start, line.start+1)[0]

        // Ignore non-primary assembly contigs because they can return blank reference sequence
        if(Region.isMinorContig(line.chr) && !ref.trim())
            return null
            
        // Skip variants where reference base is N
        if(ref.toUpperCase() == 'N') {
            log.warning("Skipping variant at ${line.chr}:${line.start} because reference base is N")
            return null
        }
                                         
        int svLen = (line.end - line.start) * (line.type == 'DEL' ? -1 : 1 )
        Allele refAllele = Allele.create(ref, true)
        
        List<String> types = line.type.tokenize(',')
        
        boolean has_cn_info = line.containsKey('copy_number')
        
        // Check if variant contains INV or DEL
        boolean hasDelOrInv = types.contains('DEL') || types.contains('INV')
        
        // Simplify to single alt allele - DUP only if it's the only type and no INV/DEL, DEL for everything else
        String altType = types.contains('DUP') && !hasDelOrInv ? 'DUP' : 'DEL'
        List<Allele> altAlleles = [Allele.create('<' + altType + '>')]
        
        // For the SVTYPE attribute, use DEL if DEL or INV is present, otherwise use original type
        String svType = hasDelOrInv ? 'DEL' : line.type
        
        // Cap copy number at 1 for DEL variants or combined DUP,DEL
        // Set minimum copy number of 3 for DUP variants
        Integer copyNumber = null
        if (has_cn_info) {
            if (svType == 'DEL' || types.contains('DEL')) {
                // DEL variants: cap at 1
                copyNumber = Math.min(1, (line.copy_number?:0) as int)
            } else if (altType == 'DUP') {
                // DUP variants: minimum of 3
                copyNumber = Math.max(3, (line.copy_number?:0) as int)
            } else {
                // Other variants: use as-is
                copyNumber = (line.copy_number?:0) as int
            }
        }
        
        if(has_cn_info && line.copy_number == null) {
            log.warning("Copy number assigned as null for $line")
        }
        
        Allele firstAllele = refAllele
        if(has_cn_info) {
            if(types[0] == 'DEL' && (line.copy_number?:0) == 0) {
                firstAllele = altAlleles[0]
            }
        }
        
        boolean has_cr_info = line.containsKey('coverage_ratio')
        
        // Default CR values when coverage_ratio not available
        double defaultCR = altType == 'DEL' ? 0.5 : 3.0
        
        List<Allele> alleles = [
            refAllele,
            *altAlleles
        ]
        
        if(!(line.sample in samples)) {
            return null
        }

        List<Genotype> gts = samples.collect {
            if(it == line.sample) {
                def formatFields =  [
                    CR : has_cr_info ? line.coverage_ratio : defaultCR,
                    NC : line.count
                ]
                return GenotypeBuilder.create(it, [firstAllele, altAlleles[0]], formatFields)
            }
            return GenotypeBuilder.create(it, [refAllele, refAllele])
        }
        
        boolean has_combined_qual = line.containsKey('combined_qual')
        double combined_qual = 20
        if(has_combined_qual) {
            combined_qual = line.combined_qual?:0
        }
        else {
            // Calculate assuming Phred scaled values b/w 0 and 100
            // clip at 100 to avoid a single caller dominating the score
            combined_qual = line*.key.grep { it.endsWith('_qual') && line[it.split('_')[0]] == 'TRUE' }
            .collect { line[it].toDouble() }
            .collect { qual ->
                Math.min(100d, Math.max(0d, qual))
            }.sum()
        }
        
        List<String> allCallers = line.findAll { it.key.endsWith('_qual') }.collect { Map.Entry e -> e.key.tokenize('_')[0] }
        
        List<String> calledBy = allCallers.grep { line[it] == "TRUE" }
        
        assert calledBy.size() == line.count

        // Compute target overlap count if target regions are provided
        Integer targetOverlapCount = null
        if(targetRegions != null) {
            Region cnvRegion = new Region(line.chr, line.start..line.end)
            targetOverlapCount = targetRegions.getOverlapRegions(cnvRegion).size()
        }
        
        // Determine filters
        Set<String> filters = computeFilters(line.count as int, targetOverlapCount, passTargets, passCallerCount)

        def builder = new VariantContextBuilder()
                .chr(line.chr)
                .start(line.start)
                .stop(line.end)
                .log10PError(-combined_qual/10)
                .attribute("SVTYPE", svType)
                .attribute("END", line.end)
                .attribute("SVLEN", svLen)
                .attribute("CR", has_cr_info ? line.coverage_ratio : defaultCR)
                .attribute("CN", has_cn_info ? copyNumber : '.')
                .attribute("CALLERS", line.count)
                .attribute("CALLEDBY", calledBy)
                .alleles((Collection)alleles)
                .genotypes(gts)
        
        if(targetOverlapCount != null) {
            builder.attribute("TARGETS", targetOverlapCount)
        }
        
        if(filters.isEmpty()) {
            builder.passFilters()
        } else {
            builder.filters(filters)
        }
        
        return builder.make()
    }
    
    /**
     * Compute the set of filters that should be applied to a variant based on caller count
     * and target overlap count.
     * 
     * @param callerCount Number of callers that called this variant
     * @param targetOverlapCount Number of target regions overlapping this variant (null if not computed)
     * @param passTargets Minimum target overlap count for PASS (null = no filter applied)
     * @param passCallerCount Minimum caller count for PASS (null = no filter applied)
     * @return Set of filter names to apply (empty set means PASS)
     */
    Set<String> computeFilters(int callerCount, Integer targetOverlapCount, Integer passTargets, Integer passCallerCount) {
        Set<String> filters = new LinkedHashSet<>()
        
        if(passCallerCount != null && callerCount < passCallerCount) {
            filters.add(FILTER_LOW_CALLERS)
        }
        
        if(passTargets != null && targetOverlapCount != null && targetOverlapCount < passTargets) {
            filters.add(FILTER_FEW_TARGETS)
        }
        
        return filters
    }

    /**
     * Build the VCF headers to output based on contigs in the file and known output fields
     * 
     * @param genomeRef
     * @param contigs
     * @param passTargets
     * @param passCallerCount
     * @return
     */
    private Set createHeaders(FASTA genomeRef, List contigs, Integer passTargets = null, Integer passCallerCount = null) {
        Map<String, Integer> contigLengths = genomeRef.contigs
        Set contigHeaderLines = contigLengths*.key
            .grep { !Region.isMinorContig(it) }
            .collect { String contig ->
                if(!contigLengths.containsKey(contig))
                    throw new IllegalStateException("Contig $contig was not found in the supplied FASTA file")
                new VCFSimpleHeaderLine("contig", [ ID: contig, length: genomeRef.contigs[contig]])
            } as Set
            
        def formatHeaderLines = [ 
            new VCFFormatHeaderLine('CR', 1, VCFHeaderLineType.Float, "Ratio of expected to observed coverage depth"),
            new VCFFormatHeaderLine('NC', 1, VCFHeaderLineType.Integer, "Count of callers supporting the CNV call")
        ] as Set
            
        Set referenceHeaderLine = 
                 [ new VCFHeaderLine("reference", "GRCh38") ] as Set

        Set sourceHeaderLine = []
        if(opts.source) {
             sourceHeaderLine = [ new VCFHeaderLine("source", opts.source) ] as Set
        }

        Set headerLines = [
            new VCFInfoHeaderLine('SVTYPE', 1, VCFHeaderLineType.String, "Type of structural variant"),
            new VCFInfoHeaderLine('SVLEN', 1, VCFHeaderLineType.Integer, "Difference in length between REF and ALT alleles"),
            new VCFInfoHeaderLine('END', 1, VCFHeaderLineType.Integer, "End position of the variant described in this record"),
            new VCFInfoHeaderLine('CN', 1, VCFHeaderLineType.Integer, "Inferred copy number"),
            new VCFInfoHeaderLine('CR', 1, VCFHeaderLineType.Integer, "Ratio of observed to expected coverage depth over event"),
            new VCFInfoHeaderLine('CALLERS', 1, VCFHeaderLineType.Integer, "Number of callers that identified the event"),
            new VCFInfoHeaderLine('CALLEDBY', 1, VCFHeaderLineType.String, "Comma separated list of callers that identified the event"),
            new VCFInfoHeaderLine('TARGETS', 1, VCFHeaderLineType.Integer, "Number of target regions overlapping the event")
        ] as Set

        Set filterHeaderLines = [] as Set
        if(passCallerCount != null) {
            filterHeaderLines.add(new VCFFilterHeaderLine(FILTER_LOW_CALLERS, "Variant called by fewer than $passCallerCount callers"))
        }
        if(passTargets != null) {
            filterHeaderLines.add(new VCFFilterHeaderLine(FILTER_FEW_TARGETS, "Variant overlaps fewer than $passTargets target regions"))
        }

        Set allHeaders = sourceHeaderLine + referenceHeaderLine + formatHeaderLines + contigHeaderLines + headerLines + filterHeaderLines

        return allHeaders
    }

    @CompileStatic
    private void writeHeader(Writer output, VCFHeader header, List<String> samples) {
        output.write(header.metaDataInInputOrder.collect { '##' + it }.join('\n') + '\n')
        output.write((['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT'] + samples).join('\t') + '\n')
    }

}
