
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
            
            String ref = genomeRef.basesAt(line.chr, line.start, line.start+1)[0]

            // Ignore non-primary assembly contigs because they can return blank reference sequence
            if(Region.isMinorContig(line.chr) && !ref.trim())
                continue
                                         
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
            
            List<Allele> alleles = [
                refAllele,
                *altAlleles
            ]
            
            Genotype gt = GenotypeBuilder.create(line.sample, alleles)
            if(line.sample in samples) {
                List<Genotype> gts = samples.collect {
                    if(it == line.sample) {
                        return GenotypeBuilder.create(it, [firstAllele, altAlleles[0]])
                    }
                    return GenotypeBuilder.create(it, [refAllele, refAllele])
                }

                VariantContext vctx = 
                    new VariantContextBuilder()
                        .chr(line.chr)
                        .start(line.start)
                        .stop(line.end)
                        .log10PError(has_cn_info ? -line.combined_qual/10 : -2)
                        .attribute("SVTYPE", line.type)
                        .attribute("END", line.end)
                        .attribute("SVLEN", svLen)
                        .attribute("CR", has_cn_info ? line.coverage_ratio : '.')
                        .attribute("CN", has_cn_info ? line.copy_number : '.')
                        .attribute("CALLERS", line.count)
                        .attribute("CALLERS", line.count)
                        .alleles((Collection)alleles)
                        .genotypes(gts)
                        .make()
                        
                outputVariants.add(vctx)
            }
            p.count()
        }
        p.end()
        
        log.info "Sorting and writing ${outputVariants.size()} output variants ..."
        
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
                new VCFSimpleHeaderLine(contig, [ ID: contig, length: genomeRef.contigs[contig]])
            } as Set

        Set headerLines = [
            new VCFInfoHeaderLine('SVTYPE', 1, VCFHeaderLineType.String, "Type of structural variant"),
            new VCFInfoHeaderLine('SVLEN', 1, VCFHeaderLineType.Integer, "Difference in length between REF and ALT alleles"),
            new VCFInfoHeaderLine('END', 1, VCFHeaderLineType.Integer, "End position of the variant described in this record"),
            new VCFInfoHeaderLine('CN', 1, VCFHeaderLineType.Integer, "Inferred copy number"),
            new VCFInfoHeaderLine('CR', 1, VCFHeaderLineType.Integer, "Ratio of observed to expected coverage depth over event")
        ] as Set

        Set allHeaders = contigHeaderLines + headerLines
        return allHeaders
    }

    @CompileStatic
    private void writeHeader(Writer output, VCFHeader header, List<String> samples) {
        output.write(header.metaDataInInputOrder.collect { '##' + it }.join('\n') + '\n')
        output.write((['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT'] + samples).join('\t') + '\n')
    }

}
