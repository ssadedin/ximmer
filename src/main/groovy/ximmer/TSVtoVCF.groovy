
package ximmer

import com.xlson.groovycsv.PropertyMapper
import gngs.FASTA
import gngs.ProgressCounter
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
    
    static void main(String[] args) {
        cli('TSVtoVCF -i <CNV report> -s <sample> [-s <sample2> ...] -o <vcf output>', 'Converts Ximmer CNV report to VCF format', args) {
            i 'CNV report produced by Ximmer SummarizeCNVs', args: 1, required: true, type: File
            r 'Reference genome for computing reference alleles at CNV start positions', args:1, required: true, type:File
            hl 'Additional header lines to add', args:'*', type: String
            s 'Sample to include', args: '*', required: true, type: String
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
        
        
        List<String> contigs = cnvReport*.chr.grep { !Region.isMinorContig(it) }.sort()
        
        opts.o.withWriter { w ->
            this.createVCF(w, ref, contigs, cnvReport)
        }
        
        log.info "Wrote $opts.o"
    }

    void createVCF(final Writer w, final FASTA genomeRef, final List<String> contigs, final List<Map> tsv) {
        List<String> samples = opts.ss
        
        Set allHeaders = createHeaders(genomeRef, contigs)

        VCFHeader header = new VCFHeader(allHeaders, samples)

        writeHeader(w, header, opts.ss)
        
        VCFEncoder encoder = new VCFEncoder(header, true, true)
       
        ProgressCounter p = new ProgressCounter(withRate: true, withTime: true, log: log)
        
        List<VariantContext> outputVariants = new ArrayList(10000)
        
        for(Map line in tsv) {
            VariantContext vctx = createVariantFromLine(line, genomeRef, samples)
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
     * @return VariantContext object if valid, null if variant should be skipped
     */
    VariantContext createVariantFromLine(Map line, FASTA genomeRef, List<String> samples) {
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
        
        boolean has_cn_info = 'copy_number' in line
        
        // Check if variant contains INV or DEL
        boolean hasDelOrInv = types.contains('DEL') || types.contains('INV')
        
        // Simplify to single alt allele - DUP only if it's the only type and no INV/DEL, DEL for everything else
        String altType = types.contains('DUP') && !hasDelOrInv ? 'DUP' : 'DEL'
        List<Allele> altAlleles = [Allele.create('<' + altType + '>')]
        
        // For the SVTYPE attribute, use DEL if DEL or INV is present, otherwise use original type
        String svType = hasDelOrInv ? 'DEL' : line.type
        
        // Cap copy number at 1 for DEL variants or combined DUP,DEL
        Integer copyNumber = has_cn_info ? (
            (svType == 'DEL' || types.contains('DEL')) ? 
                Math.min(1, line.copy_number as int) : 
                line.copy_number as int
        ) : null
        
        Allele firstAllele = refAllele
        if(has_cn_info) {
            if(types[0] == 'DEL' && line.copy_number == 0) {
                firstAllele = altAlleles[0]
            }
        }
        
        boolean has_cr_info = 'coverage_ratio' in line
        
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
        
        boolean has_combined_qual = ('combined_qual' in line)
        double combined_qual = 20
        if(has_combined_qual) {
            combined_qual = line.combined_qual
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

        return new VariantContextBuilder()
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
                .make()
    }

    /**
     * Build the VCF headers to output based on contigs in the file and known output fields
     * 
     * @param genomeRef
     * @param contigs
     * @return
     */
    private Set createHeaders(FASTA genomeRef, List contigs) {
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
            new VCFInfoHeaderLine('CR', 1, VCFHeaderLineType.Integer, "Ratio of observed to expected coverage depth over event")
        ] as Set

        Set allHeaders = sourceHeaderLine + referenceHeaderLine + formatHeaderLines + contigHeaderLines + headerLines

        return allHeaders
    }

    @CompileStatic
    private void writeHeader(Writer output, VCFHeader header, List<String> samples) {
        output.write(header.metaDataInInputOrder.collect { '##' + it }.join('\n') + '\n')
        output.write((['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT'] + samples).join('\t') + '\n')
    }

}
