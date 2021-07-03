
package ximmer

import com.xlson.groovycsv.PropertyMapper
import gngs.FASTA
import gngs.ProgressCounter
import gngs.ToolBase
import gngs.XPos
import graxxia.TSV
import groovy.transform.CompileStatic
import groovy.util.logging.Log
import htsjdk.variant.variantcontext.*
import htsjdk.variant.vcf.*

@Log
class TSVtoVCF extends ToolBase {
    
    static void main(String[] args) {
        cli('TSVtoVCF -i <CNV report> -s <sample> [-s <sample2> ...] -o <vcf output>', 'Converts Ximmer CNV report to VCF format', args) {
            i 'CNV report produced by Ximmer SummarizeCNVs', args: 1, required: true, type: File
            r 'Reference genome for computing reference alleles at CNV start positions', args:1, required: true, type:File
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
        opts.o.withWriter { w ->
            this.createVCF(w, ref, cnvReport)
        }
    }

    void createVCF(final Writer w, final FASTA genomeRef, final TSV tsv) {
        
        List<String> samples = opts.ss

        List headerLines = [ 
            new VCFInfoHeaderLine('SVTYPE', 1, VCFHeaderLineType.String, "Type of structural variant"),
            new VCFInfoHeaderLine('SVLEN', 1, VCFHeaderLineType.Integer, "Difference in length between REF and ALT alleles"),
            new VCFInfoHeaderLine('END', 1, VCFHeaderLineType.Integer, "End position of the variant described in this record")
        ]

        VCFHeader header = new VCFHeader(headerLines as Set, samples)

        writeHeader(w, header, opts.ss)
        
        VCFEncoder encoder = new VCFEncoder(header, true, true)
       
        ProgressCounter p = new ProgressCounter(withRate: true, withTime: true, log: log)
        
        List<VariantContext> outputVariants = new ArrayList(10000)
        
        for(PropertyMapper line in tsv) {
            String ref = genomeRef.basesAt(line.chr, line.start, line.start+1)[0]
            int svLen = (line.end - line.start) * (line.type == 'DEL' ? -1 : 1 )
            Allele refAllele = Allele.create(ref, true)
            List<Allele> alleles = [
                refAllele,
                *line.type.tokenize(',').collect { Allele.create('<' + it + '>')}
            ]
            
//            Genotype gt = new GenotypeBuilder().alleles(alleles).GQ(line.count)..make()
            Genotype gt = GenotypeBuilder.create(line.sample, alleles)
            
            if(line.sample in samples) {
                List<Genotype> gts = samples.collect {
                    it == line.sample ? gt : GenotypeBuilder.create(it, [refAllele, refAllele])
                }

                VariantContext vctx = 
                    new VariantContextBuilder()
                        .chr(line.chr)
                        .start(line.start)
                        .stop(line.end)
                        .attribute("SVTYPE", line.type)
                        .attribute("END", line.end)
                        .attribute("SVLEN", svLen)
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

    @CompileStatic
    private void writeHeader(Writer output, VCFHeader header, List<String> samples) {
        output.write(header.metaDataInInputOrder.collect { '##' + it }.join('\n') + '\n')
        output.write((['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT'] + samples).join('\t') + '\n')
    }

}
