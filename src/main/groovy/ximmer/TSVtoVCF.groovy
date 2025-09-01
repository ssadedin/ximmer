
package ximmer

import java.awt.IllegalComponentStateException

import com.xlson.groovycsv.PropertyMapper
import gngs.FASTA
import gngs.ProgressCounter
import gngs.ToolBase
import gngs.XPos
import gngs.Region
import graxxia.TSV
import groovy.transform.CompileStatic
import groovy.util.logging.Log
import htsjdk.variant.variantcontext.*
import htsjdk.variant.vcf.*

/**
 * Converts the CNV report produced by SummarizeCNVs into VCF format. If the form with copy number
 * information is supplied, includes that in the genotyping column.
 */
@Log
class TSVtoVCF extends ToolBase {
    
    static void main(String[] args) {
        cli('TSVtoVCF -i <CNV report> -s <sample> [-s <sample2> ...] -o <vcf output>', 'Converts Ximmer CNV report to VCF format', args) {
            i 'CNV report produced by Ximmer SummarizeCNVs', args: 1, required: true, type: File
            r 'Reference genome for computing reference alleles at CNV start positions', args:1, required: true, type:File
            hl 'Additional header lines to add', args:'*', type: String
            s 'Sample to include', args: '*', required: true, type: String
            o 'VCF output file', args:1, required: true, type: File
        }
    }

    @Override
    public void run() {
        FASTA ref = new FASTA(opts.r.absolutePath)
        TSV cnvReport = new TSV(opts.i)
        log.info "Converting $opts.i to VCF format"
        log.info "Using reference: $opts.r"
        
        
        TSV scanTSV = new TSV(opts.i)
        List<String> contigs = scanTSV*.chr.grep { !Region.isMinorContig(it) }.sort()
        
        opts.o.withWriter { w ->
            this.createVCF(w, ref, contigs, cnvReport)
        }
        
        log.info "Wrote $opts.o"
    }

    void createVCF(final Writer w, final FASTA genomeRef, final List<String> contigs, final TSV tsv) {
        List<String> samples = opts.ss
        
        Set allHeaders = createHeaders(genomeRef, contigs)

        VCFHeader header = new VCFHeader(allHeaders, samples)

        writeHeader(w, header, opts.ss)
        
        VCFEncoder encoder = new VCFEncoder(header, true, true)
       
        ProgressCounter p = new ProgressCounter(withRate: true, withTime: true, log: log)
        
        List<VariantContext> outputVariants = new ArrayList(10000)
        
        for(PropertyMapper line in tsv) {
            VariantContext vctx = createVariantFromLine(line, genomeRef, samples)
            if(vctx) {
                outputVariants.add(vctx)
            }
            p.count()
        }
        p.end()
        
        log.info "Sorting and writing ${outputVariants.size()} output variants ..."
        
    }
    
    /**
     * Creates a VariantContext from a single line of TSV input
     * 
     * @param line The line from the TSV file containing variant information
     * @param genomeRef Reference genome for getting reference bases
     * @param samples List of samples to include in the output
     * @return VariantContext object if valid, null if variant should be skipped
     */
    private VariantContext createVariantFromLine(PropertyMapper line, FASTA genomeRef, List<String> samples) {
        String ref = genomeRef.basesAt(line.chr, line.start, line.start+1)[0]

        // Ignore non-primary assembly contigs because they can return blank reference sequence
        if(Region.isMinorContig(line.chr) && !ref.trim())
            return null
                                         
        int svLen = (line.end - line.start) * (line.type == 'DEL' ? -1 : 1 )
        Allele refAllele = Allele.create(ref, true)
        
        List<String> types = line.type.tokenize(',')
        List<Allele> altAlleles = types.collect { Allele.create('<' + it + '>')}
        
        Allele firstAllele = refAllele
        boolean has_cn_info = 'copy_number' in line.columns
        if(has_cn_info) {
            if(types[0] == 'DEL' && line.copy_number == 0) {
                firstAllele = altAlleles[0]
            }
        }
        
        boolean has_cr_info = 'coverage_ratio' in line.columns
        
        List<Allele> alleles = [
            refAllele,
            *altAlleles
        ]
        
        if(line.sample in samples) {
            List<Genotype> gts = samples.collect {
                if(it == line.sample) {
                    def formatFields =  [
                        CR : has_cr_info ? line.coverage_ratio : null,
                        NC : has_cn_info ? line.count : null
                    ]
                    return GenotypeBuilder.create(it, [firstAllele, altAlleles[0]], formatFields)
                }
                return GenotypeBuilder.create(it, [refAllele, refAllele])
            }
            
            boolean has_combined_qual = ('combined_qual' in line.columns)
            double combined_qual = 20
            if(has_combined_qual) {
                combined_qual = line.combined_qual
            }
            else {
                // Calculate assuming Phred scaled values b/w 0 and 100
                // clip at 100 to avoid a single caller dominating the score
                combined_qual = line.columns*.key.grep { it.endsWith('_qual') && line[it.split('_')[0]] == 'TRUE' }
                .collect { line[it].toDouble() }
                .collect { qual ->
                    Math.min(100d, Math.max(0d, qual))
                }.sum()
            }

            return new VariantContextBuilder()
                    .chr(line.chr)
                    .start(line.start)
                    .stop(line.end)
                    .log10PError(-combined_qual/10)
                    .attribute("SVTYPE", line.type)
                    .attribute("END", line.end)
                    .attribute("SVLEN", svLen)
                    .attribute("CR", has_cr_info ? line.coverage_ratio : '.')
                    .attribute("CN", has_cn_info ? line.copy_number : '.')
                    .attribute("CALLERS", line.count)
                    .alleles((Collection)alleles)
                    .genotypes(gts)
                    .make()
        }
        return null
        
        List<VariantContext> sortedVariants = outputVariants.sort { XPos.computePos(it.contig, it.start)}
        for(vctx in sortedVariants) {
                encoder.write(w, vctx)
                w.write('\n')
        }
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
                    throw new IllegalComponentStateException("Contig $contig was not found in the supplied FASTA file")
                new VCFSimpleHeaderLine("contig", [ ID: contig, length: genomeRef.contigs[contig]])
            } as Set
            
        def formatHeaderLines = [ 
            new VCFFormatHeaderLine('CR', 1, VCFHeaderLineType.Float, "Ratio of expected to observed coverage depth"),
            new VCFFormatHeaderLine('NC', 1, VCFHeaderLineType.Integer, "Count of callers supporting the CNV call")
        ] as Set
            
        Set referenceHeaderLine = 
                 [ new VCFSimpleHeaderLine("reference", "GRCh38", "Reference file") ] as Set
                 

        Set headerLines = [
            new VCFInfoHeaderLine('SVTYPE', 1, VCFHeaderLineType.String, "Type of structural variant"),
            new VCFInfoHeaderLine('SVLEN', 1, VCFHeaderLineType.Integer, "Difference in length between REF and ALT alleles"),
            new VCFInfoHeaderLine('END', 1, VCFHeaderLineType.Integer, "End position of the variant described in this record"),
            new VCFInfoHeaderLine('CN', 1, VCFHeaderLineType.Integer, "Inferred copy number"),
            new VCFInfoHeaderLine('CR', 1, VCFHeaderLineType.Integer, "Ratio of observed to expected coverage depth over event")
        ] as Set

        Set allHeaders = referenceHeaderLine + formatHeaderLines + contigHeaderLines + headerLines

        return allHeaders
    }

    @CompileStatic
    private void writeHeader(Writer output, VCFHeader header, List<String> samples) {
        output.write(header.metaDataInInputOrder.collect { '##' + it }.join('\n') + '\n')
        output.write((['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT'] + samples).join('\t') + '\n')
    }

}
