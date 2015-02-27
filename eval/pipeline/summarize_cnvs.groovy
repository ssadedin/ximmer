// vim: expandtab:sw=4:ts=4:cindent
touch_chr = {
   exec """
        echo $chr > $output.txt
   """
}

extract_sample_files = {
    requires sample_info : "Sample meta data object"
    var sample : branch.name
    branch.sample = sample
    branch.sample_info = sample_info
    println "Forwarding files for sample $sample : " + branch.sample_info[sample].files.all
    forward branch.sample_info[sample].files.all
}

extract_cnv_regions = {
    doc """"
        Extracts regions where CNVs were called from a BAM file, based on the 
        universal CNV format
        """

    // Slop gives the additional upstream and downstream flanking
    // flanking region to be included
    var slop : 500 

    transform("bam") to("cnvs.bam") {
        exec """
             set -o pipefail

             JAVA_OPTS="-Xmx1g" $GROOVY -cp tools/groovy-ngs-utils/1.0/groovy-ngs-utils.jar -e 'new RangedData("$input.tsv").load().each { println([it.chr, it.from-$slop, it.to+$slop].join("\\t"))  }' | 
                 $SAMTOOLS view -b -L - $input.bam > $output.bam
        """
    }
}

summarize_cnvs = {

    output.dir="runs/$branch.dir/report"

    var true_cnvs : "true_cnvs.bed",
        create_plots : "T",
        chr : "*"

    println "Summarizing CNVs ..."

    List vcfs = []
    String hasVCFs = sample_info.every { it.value.files.vcf } ? "TRUE" : "FALSE"
    if(hasVCFs == "TRUE")
        vcfs = sample_info.collect{"'${it.value.files.vcf[0]}'"}
    else 
        println "No VCF files present for one or more samples: variant heterozygosity will not be annotated for CNV calls"

    println "Using callers : " + cnv_callers

    def output_prefix = (chr == "*" ? "" : "." + chr)
    produce(batch_name + ".cnv${output_prefix}.summary.tsv", batch_name +"${output_prefix}.plotdata.RData") {

        def result_list = []
        if(callers.contains("xhmm")) 
            result_list << "xhmm = load_xhmm_results('$input.xcnv')"
        if(callers.contains("ec")) 
            result_list << "ec = load_exome_copy_results('$input.exome_copy.cnvs.tsv')"
        if(callers.contains("ex")) 
            result_list << "ex = load_excavator_results('$input.excavator.cnvs.tsv')"
        if(callers.contains("ed")) 
            result_list << "ed = load_exomedepth_results(file.name='$input.exome_depth.cnvs.tsv')"
        if(callers.contains("mops")) 
            result_list << "mops = load_cn_mops_results(file.name='$input.cn_mops_call_cnvs.tsv')"
        if(callers.contains("ang")) 
            result_list << "ang = load_angel_results(file.name='$input.angel.cnvs.bed')"
        if(simulation) 
            result_list << "truth = load_angel_results(file.name='$input.true_cnvs.bed')"

        println "Running R code ..."
        R {"""
            source("./tools/r-utils/cnv_utils.R")

            cnv.samples = c(${sample_names.collect{"'$it'"}.join(",")})

            target.bed = read.bed.ranges("$target_bed")
            cnv.results = list(
                ${result_list.join(",\n")}
            )
            cnv.all = combine_cnv_caller_results(cnv.samples, cnv.results, target.bed, chr="$chr")

            if($hasVCFs) {
                vcf_files = gsub(".vcf\$", ".vcf.bgz", c(${vcfs.join(",")}))

                # Figure out the siez of the chromosome from VCF header (assume all samples the same)
                hdr = scanBcfHeader(vcf_files[[1]])[[1]]
                chr.end = as.integer(hdr$Header$contig['$chr','length'])

                vcfs = sapply(cnv.samples, function(s) {
                    v = readVcf(file=vcf_files[[which(cnv.samples==s)]],
                                genome='hg19', 
                                param=ScanVcfParam(which=GRanges(seqnames='$chr',ranges=IRanges(start=0, end=chr.end))))

                    #v = readVcf(file=vcf_files[[which(cnv.samples==s)]],genome="h19")
                    seqlevels(v,force=T) = hg19.chromosomes
                    return(v)
                }, USE.NAMES=T)

                cnv.all = annotate_variant_counts(vcfs, cnv.all)

                cnv.output = data.frame(as.data.frame(cnv.all),
                    vcf=sapply(cnv.all$sample, function(s) {vcf_files[[which(cnv.samples==s)]]})
                )
            } else {
                cnv.output = as.data.frame(cnv.all)
            }

            chr='$chr'
            if(chr=='*') {
               chr='all'
             }
            write.table(as.data.frame(cnv.output), 
                        file="${output.tsv}",
                        quote=F,
                        sep='\\t',
                        row.names=F
                        )

            save.image(file="$output.RData")
        """}
    }
}

combine_cnv_report = {
    doc "Combine together CNV reports that are created for multiple chromosomes"

    output.dir = "$branch.dir/report"

    exec """
        cat $inputs.summary.tsv | awk '{ if((NR==1) || (\$1!="seqnames")) print \$0; }' > ${output(batch_name+".cnv.summary.tsv")}
    """
}

cnv_report = {

    requires batch_name : "Name of batch for which report is being generated",
	     callers : "Names of CNV callers in the CNV summary file" // see summarize_cnvs

    var inline_js : true,
        cnv_summary : input.cnv.summary.tsv,
        types : ["DEL","DUP"]

    produce("$branch.dir/report/cnv_report.html", "cnv_report.zip") {
        send report('templates/cnv_report.html') to file: output.html
        exec """
            zip -r $output.zip $branch.dir/report -x \\*.RData report
        """
    }
}

