
package ximmer

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
            new VCFContigHeaderLine("##contig=<ID=chr1,length=248956422,assembly=38>", VCFHeaderVersion.VCF4_2, "contig", 0),
            new VCFContigHeaderLine("##contig=<ID=chr2,length=242193529,assembly=38>", VCFHeaderVersion.VCF4_2, "contig", 1),
            new VCFContigHeaderLine("##contig=<ID=chr3,length=198295559,assembly=38>", VCFHeaderVersion.VCF4_2, "contig", 2),
            new VCFContigHeaderLine("##contig=<ID=chr4,length=190214555,assembly=38>", VCFHeaderVersion.VCF4_2, "contig", 3),
            new VCFContigHeaderLine("##contig=<ID=chr5,length=181538259,assembly=38>", VCFHeaderVersion.VCF4_2, "contig", 4),
            new VCFContigHeaderLine("##contig=<ID=chr6,length=170805979,assembly=38>", VCFHeaderVersion.VCF4_2, "contig", 5),
            new VCFContigHeaderLine("##contig=<ID=chr7,length=159345973,assembly=38>", VCFHeaderVersion.VCF4_2, "contig", 6),
            new VCFContigHeaderLine("##contig=<ID=chr8,length=145138636,assembly=38>", VCFHeaderVersion.VCF4_2, "contig", 7),
            new VCFContigHeaderLine("##contig=<ID=chr9,length=138394717,assembly=38>", VCFHeaderVersion.VCF4_2, "contig", 8),
            new VCFContigHeaderLine("##contig=<ID=chr10,length=133797422,assembly=38>", VCFHeaderVersion.VCF4_2, "contig", 9),
            new VCFContigHeaderLine("##contig=<ID=chr11,length=135086622,assembly=38>", VCFHeaderVersion.VCF4_2, "contig", 10),
            new VCFContigHeaderLine("##contig=<ID=chr12,length=133275309,assembly=38>", VCFHeaderVersion.VCF4_2, "contig", 11),
            new VCFContigHeaderLine("##contig=<ID=chr13,length=114364328,assembly=38>", VCFHeaderVersion.VCF4_2, "contig", 12),
            new VCFContigHeaderLine("##contig=<ID=chr14,length=107043718,assembly=38>", VCFHeaderVersion.VCF4_2, "contig", 13),
            new VCFContigHeaderLine("##contig=<ID=chr15,length=101991189,assembly=38>", VCFHeaderVersion.VCF4_2, "contig", 14),
            new VCFContigHeaderLine("##contig=<ID=chr16,length=90338345,assembly=38>", VCFHeaderVersion.VCF4_2, "contig", 15),
            new VCFContigHeaderLine("##contig=<ID=chr17,length=83257441,assembly=38>", VCFHeaderVersion.VCF4_2, "contig", 16),
            new VCFContigHeaderLine("##contig=<ID=chr18,length=80373285,assembly=38>", VCFHeaderVersion.VCF4_2, "contig", 17),
            new VCFContigHeaderLine("##contig=<ID=chr19,length=58617616,assembly=38>", VCFHeaderVersion.VCF4_2, "contig", 18),
            new VCFContigHeaderLine("##contig=<ID=chr20,length=64444167,assembly=38>", VCFHeaderVersion.VCF4_2, "contig", 19),
            new VCFContigHeaderLine("##contig=<ID=chr21,length=46709983,assembly=38>", VCFHeaderVersion.VCF4_2, "contig", 20),
            new VCFContigHeaderLine("##contig=<ID=chr22,length=50818468,assembly=38>", VCFHeaderVersion.VCF4_2, "contig", 21),
            new VCFContigHeaderLine("##contig=<ID=chrX,length=156040895,assembly=38>", VCFHeaderVersion.VCF4_2, "contig", 22),
            new VCFContigHeaderLine("##contig=<ID=chrY,length=57227415,assembly=38>", VCFHeaderVersion.VCF4_2, "contig", 23),
            new VCFContigHeaderLine("##contig=<ID=chrM,length=16569,assembly=38>", VCFHeaderVersion.VCF4_2, "contig", 24),
            new VCFFormatHeaderLine("CONC_ST", VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "The genotype concordance contingency state"),
            new VCFFormatHeaderLine("CN", 1, VCFHeaderLineType.Integer, "Predicted copy state"),
            new VCFFormatHeaderLine("CNQ", 1, VCFHeaderLineType.Integer, "Read-depth genotype quality"),
            new VCFFormatHeaderLine("EV", VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Classes of evidence supporting final genotype"),
            new VCFFormatHeaderLine("GQ", 1, VCFHeaderLineType.Integer, "Genotype Quality"),
            new VCFFormatHeaderLine("GT", 1, VCFHeaderLineType.String, "Genotype"),
            new VCFFormatHeaderLine("PE_GQ", 1, VCFHeaderLineType.Integer, "Paired-end genotype quality"),
            new VCFFormatHeaderLine("PE_GT", 1, VCFHeaderLineType.Integer, "Paired-end genotype"),
            new VCFFormatHeaderLine("RD_CN", 1, VCFHeaderLineType.Integer, "Predicted copy state"),
            new VCFFormatHeaderLine("RD_GQ", 1, VCFHeaderLineType.Integer, "Read-depth genotype quality"),
            new VCFFormatHeaderLine("SR_GQ", 1, VCFHeaderLineType.Integer, "Split read genotype quality"),
            new VCFFormatHeaderLine("SR_GT", 1, VCFHeaderLineType.Integer, "Split-read genotype"),
            new VCFInfoHeaderLine("ALGORITHMS", VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Source algorithms"),
            new VCFInfoHeaderLine('CALLERS', 1, VCFHeaderLineType.Integer, "MCRI, number of SV callers used"),
            new VCFInfoHeaderLine("CHR2", 1, VCFHeaderLineType.String, "Chromosome for END coordinate"),
            new VCFInfoHeaderLine("CPX_INTERVALS", VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Genomic intervals constituting complex variant."),
            new VCFInfoHeaderLine("CPX_TYPE", 1, VCFHeaderLineType.String, "Class of complex variant."),
            new VCFInfoHeaderLine("END", 1, VCFHeaderLineType.Integer, "End position of the structural variant"),
            new VCFInfoHeaderLine("END2", 1, VCFHeaderLineType.Integer, "Position of breakpoint on CHR2"),
            new VCFInfoHeaderLine("EVIDENCE", VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Classes of random forest support."),
            new VCFInfoHeaderLine("SOURCE", 1, VCFHeaderLineType.String, "Source of inserted sequence."),
            new VCFInfoHeaderLine("STRANDS", 1, VCFHeaderLineType.String, "Breakpoint strandedness [++,+-,-+,--]"),
            new VCFInfoHeaderLine("SVLEN", 1, VCFHeaderLineType.Integer, "SV length"),
            new VCFInfoHeaderLine("SVTYPE", 1, VCFHeaderLineType.String, "Type of structural variant"),
            new VCFInfoHeaderLine("UNRESOLVED_TYPE", 1, VCFHeaderLineType.String, "Class of unresolved variant."),
            new VCFInfoHeaderLine("PREDICTED_BREAKEND_EXONIC", VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Gene(s) for which the SV breakend is predicted to fall in an exon."),
            new VCFInfoHeaderLine("PREDICTED_COPY_GAIN", VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Gene(s) on which the SV is predicted to have a copy-gain effect."),
            new VCFInfoHeaderLine("PREDICTED_DUP_PARTIAL", VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Gene(s) which are partially overlapped by an SV's duplication, but the transcription start site is not duplicated."),
            new VCFInfoHeaderLine("PREDICTED_INTERGENIC", 0, VCFHeaderLineType.Flag, "SV does not overlap any protein-coding genes."),
            new VCFInfoHeaderLine("PREDICTED_INTRAGENIC_EXON_DUP", VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Gene(s) on which the SV is predicted to result in intragenic exonic duplication without breaking any coding sequences."),
            new VCFInfoHeaderLine("PREDICTED_INTRONIC", VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Gene(s) where the SV was found to lie entirely within an intron."),
            new VCFInfoHeaderLine("PREDICTED_INV_SPAN", VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Gene(s) which are entirely spanned by an SV's inversion."),
            new VCFInfoHeaderLine("PREDICTED_LOF", VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Gene(s) on which the SV is predicted to have a loss-of-function effect."),
            new VCFInfoHeaderLine("PREDICTED_MSV_EXON_OVERLAP", VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Gene(s) on which the multiallelic SV would be predicted to have a LOF, INTRAGENIC_EXON_DUP, COPY_GAIN, DUP_PARTIAL, TSS_DUP, or PARTIAL_EXON_DUP annotation if the SV were biallelic."),
            new VCFInfoHeaderLine("PREDICTED_NEAREST_TSS", VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Nearest transcription start site to an intergenic variant."),
            new VCFInfoHeaderLine("PREDICTED_NONCODING_BREAKPOINT", VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Class(es) of noncoding elements disrupted by SV breakpoint."),
            new VCFInfoHeaderLine("PREDICTED_NONCODING_SPAN", VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Class(es) of noncoding elements spanned by SV."),
            new VCFInfoHeaderLine("PREDICTED_PARTIAL_EXON_DUP", VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Gene(s) where the duplication SV has one breakpoint in the coding sequence."),
            new VCFInfoHeaderLine("PREDICTED_PROMOTER", VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Gene(s) for which the SV is predicted to overlap the promoter region."),
            new VCFInfoHeaderLine("PREDICTED_TSS_DUP", VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Gene(s) for which the SV is predicted to duplicate the transcription start site."),
            new VCFInfoHeaderLine("PREDICTED_UTR", VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Gene(s) for which the SV is predicted to disrupt a UTR."),
            new VCFInfoHeaderLine("AN", 1, VCFHeaderLineType.Integer, "Total number of alleles genotyped (for biallelic sites) or individuals with copy-state estimates (for multiallelic sites)."),
            new VCFInfoHeaderLine("AC", VCFHeaderLineCount.A, VCFHeaderLineType.Integer, "Number of non-reference alleles observed (for biallelic sites) or individuals at each copy state (for multiallelic sites)."),
            new VCFInfoHeaderLine("AF", VCFHeaderLineCount.A, VCFHeaderLineType.Float, "Allele frequency (for biallelic sites) or copy-state frequency (for multiallelic sites)."),
            new VCFInfoHeaderLine("N_BI_GENOS", 1, VCFHeaderLineType.Integer, "Total number of individuals with complete genotypes (biallelic sites only)."),
            new VCFInfoHeaderLine("N_HOMREF", 1, VCFHeaderLineType.Integer, "Number of individuals with homozygous reference genotypes (biallelic sites only)."),
            new VCFInfoHeaderLine("N_HET", 1, VCFHeaderLineType.Integer, "Number of individuals with heterozygous genotypes (biallelic sites only)."),
            new VCFInfoHeaderLine("N_HOMALT", 1, VCFHeaderLineType.Integer, "Number of individuals with homozygous alternate genotypes (biallelic sites only)."),
            new VCFInfoHeaderLine("FREQ_HOMREF", 1, VCFHeaderLineType.Float, "Homozygous reference genotype frequency (biallelic sites only)."),
            new VCFInfoHeaderLine("FREQ_HET", 1, VCFHeaderLineType.Float, "Heterozygous genotype frequency (biallelic sites only)."),
            new VCFInfoHeaderLine("FREQ_HOMALT", 1, VCFHeaderLineType.Float, "Homozygous alternate genotype frequency (biallelic sites only)."),
            new VCFInfoHeaderLine("MALE_AN", 1, VCFHeaderLineType.Integer, "Total number of MALE alleles genotyped (for biallelic sites) or MALE individuals with copy-state estimates (for multiallelic sites)."),
            new VCFInfoHeaderLine("MALE_AC", VCFHeaderLineCount.A, VCFHeaderLineType.Integer, "Number of non-reference MALE alleles observed (for biallelic sites) or MALE individuals at each copy state (for multiallelic sites)."),
            new VCFInfoHeaderLine("MALE_AF", VCFHeaderLineCount.A, VCFHeaderLineType.Float, "MALE allele frequency (for biallelic sites) or MALE copy-state frequency (for multiallelic sites)."),
            new VCFInfoHeaderLine("MALE_N_BI_GENOS", 1, VCFHeaderLineType.Integer, "Total number of MALE individuals with complete genotypes (biallelic sites only)."),
            new VCFInfoHeaderLine("MALE_N_HOMREF", 1, VCFHeaderLineType.Integer, "Number of MALE individuals with homozygous reference genotypes (biallelic sites only)."),
            new VCFInfoHeaderLine("MALE_N_HET", 1, VCFHeaderLineType.Integer, "Number of MALE individuals with heterozygous genotypes (biallelic sites only)."),
            new VCFInfoHeaderLine("MALE_N_HOMALT", 1, VCFHeaderLineType.Integer, "Number of MALE individuals with homozygous alternate genotypes (biallelic sites only)."),
            new VCFInfoHeaderLine("MALE_FREQ_HOMREF", 1, VCFHeaderLineType.Float, "MALE homozygous reference genotype frequency (biallelic sites only)."),
            new VCFInfoHeaderLine("MALE_FREQ_HET", 1, VCFHeaderLineType.Float, "MALE heterozygous genotype frequency (biallelic sites only)."),
            new VCFInfoHeaderLine("MALE_FREQ_HOMALT", 1, VCFHeaderLineType.Float, "MALE homozygous alternate genotype frequency (biallelic sites only)."),
            new VCFInfoHeaderLine("FEMALE_AN", 1, VCFHeaderLineType.Integer, "Total number of FEMALE alleles genotyped (for biallelic sites) or FEMALE individuals with copy-state estimates (for multiallelic sites)."),
            new VCFInfoHeaderLine("FEMALE_AC", VCFHeaderLineCount.A, VCFHeaderLineType.Integer, "Number of non-reference FEMALE alleles observed (for biallelic sites) or FEMALE individuals at each copy state (for multiallelic sites)."),
            new VCFInfoHeaderLine("FEMALE_AF", VCFHeaderLineCount.A, VCFHeaderLineType.Float, "FEMALE allele frequency (for biallelic sites) or FEMALE copy-state frequency (for multiallelic sites)."),
            new VCFInfoHeaderLine("FEMALE_N_BI_GENOS", 1, VCFHeaderLineType.Integer, "Total number of FEMALE individuals with complete genotypes (biallelic sites only)."),
            new VCFInfoHeaderLine("FEMALE_N_HOMREF", 1, VCFHeaderLineType.Integer, "Number of FEMALE individuals with homozygous reference genotypes (biallelic sites only)."),
            new VCFInfoHeaderLine("FEMALE_N_HET", 1, VCFHeaderLineType.Integer, "Number of FEMALE individuals with heterozygous genotypes (biallelic sites only)."),
            new VCFInfoHeaderLine("FEMALE_N_HOMALT", 1, VCFHeaderLineType.Integer, "Number of FEMALE individuals with homozygous alternate genotypes (biallelic sites only)."),
            new VCFInfoHeaderLine("FEMALE_FREQ_HOMREF", 1, VCFHeaderLineType.Float, "FEMALE homozygous reference genotype frequency (biallelic sites only)."),
            new VCFInfoHeaderLine("FEMALE_FREQ_HET", 1, VCFHeaderLineType.Float, "FEMALE heterozygous genotype frequency (biallelic sites only)."),
            new VCFInfoHeaderLine("FEMALE_FREQ_HOMALT", 1, VCFHeaderLineType.Float, "FEMALE homozygous alternate genotype frequency (biallelic sites only)."),
            new VCFInfoHeaderLine("gnomAD_V2_SVID", 1, VCFHeaderLineType.String, "Allele frequency (for biallelic sites) or copy-state frequency (for multiallelic sites) of an overlapping event in gnomad."),
            new VCFInfoHeaderLine("gnomAD_V2_AF", 1, VCFHeaderLineType.Float, "Allele frequency (for biallelic sites) or copy-state frequency (for multiallelic sites) of an overlapping event in gnomad."),
            new VCFInfoHeaderLine("gnomAD_V2_AC", 1, VCFHeaderLineType.Float, "Allele frequency (for biallelic sites) or copy-state frequency (for multiallelic sites) of an overlapping event in gnomad."),
            new VCFInfoHeaderLine("gnomAD_V2_AN", 1, VCFHeaderLineType.Float, "Allele frequency (for biallelic sites) or copy-state frequency (for multiallelic sites) of an overlapping event in gnomad."),
            new VCFInfoHeaderLine("StrVCTVRE", 1, VCFHeaderLineType.String, "StrVCTVRE score"),
        ]
        VCFHeader header = new VCFHeader(headerLines as Set, samples)

        writeHeader(w, header, opts.ss)
        
        VCFEncoder encoder = new VCFEncoder(header, true, true)
       
        ProgressCounter p = new ProgressCounter(withRate: true, withTime: true, log: log)
        
        List<VariantContext> outputVariants = new ArrayList(10000)
        
        for(PropertyMapper line in tsv) {

            String ref = genomeRef.basesAt(line.chr, line.start, line.start+1)[0]

            // Ignore non-primary assembly contigs because they can return blank reference sequence
            if(Region.isMinorContig(line.chr) && !ref.trim())
                continue

            line.type.tokenize(',').each { String svType ->
                int svLen = (line.end - line.start) * (svType == 'DEL' ? -1 : 1 )
                Allele refAllele = Allele.create(ref, true)
                List<Allele> alleles = [
                    refAllele,
                    Allele.create('<' + svType + '>')
                ]
                
                Genotype gt = GenotypeBuilder.create(line.sample, alleles)
                
                if(line.sample in samples) {
                    List<Genotype> gts = samples.collect {
                        it == line.sample ? gt : GenotypeBuilder.create(it, [Allele.NO_CALL, Allele.NO_CALL])
                    }

                    // E.g. cxpInterval: DUP_chr1:1499897-1499974
                    String cpxInterval = "${svType}_${line.chr}:${line.start}-${line.end}"

                    VariantContext vctx = 
                        new VariantContextBuilder()
                            .chr(line.chr)
                            .start(line.start)
                            .stop(line.end)
                            .attribute("SVTYPE", svType)
                            .attribute("END", line.end)
                            .attribute("SVLEN", svLen)
                            .attribute("CALLERS", line.count)
                            .attribute("ALGORITHMS", "XIMMER")
                            .attribute("CPX_INTERVALS", Arrays.asList(cpxInterval))
                            .attribute("AC", Arrays.asList(1, 1))
                            .attribute("AN", 1)
                            .attribute("AF", Arrays.asList(1.0, 1.0))
                            .attribute("gnomAD_V2_AC", 1)
                            .attribute("gnomAD_V2_AN", 1)
                            .attribute("gnomAD_V2_AF", 1.0)
                            .attribute("gnomAD_V2_SVID", 1)
                            .attribute("N_HOMALT", 1)
                            .attribute("N_HET", 1)
                            .attribute("StrVCTVRE", '1.0')
                            .alleles((Collection)alleles)
                            .genotypes(gts)
                            .make()
                            
                    outputVariants.add(vctx)
                }
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
