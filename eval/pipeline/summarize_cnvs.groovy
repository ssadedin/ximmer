
touch_chr = {
   exec """
        echo $chr > $output.txt
   """
}

extract_sample_files = {
    requires sample_info : "Sample meta data file"
    var sample : branch.name

    branch.sample = sample

    if(sample_info instanceof String) {
        branch.sample_info = SampleInfo.parse_sample_info(sample_info)
    }
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

             $GROOVY -e 'new RangedData("$input.tsv").load().each { println([it.chr, it.from-$slop, it.to+$slop].join("\\t"))  }' | 
                 $SAMTOOLS view -b -L - $input.bam > $output.bam
        """
    }
}

sample_coverage = {
    branch.sample = branch.name
    branch.bam = all_samples[sample].files.bam[0]
    from(bam) transform("coverage.txt") {
        exec """
          coverageBed -d  -abam ${bam} -b $target_bed > $output
        """

        println "Assigning coverage $output for sample $sample"
        all_samples[sample].files.coverage = output
    }
}

summarize_cnvs = {

    output.dir="runs/$batch_name/report"

    var true_cnvs : "true_cnvs.bed",
        create_plots : "T",
        chr : "*"

    println "Summarizing CNVs ..."

    // branch.callers = [ "xhmm","ex","ed","ang","truth" ]
    // doesn't satisfy requires below :-(
    def sourceFiles = []
    if(callers.contains("ed"))
        sourceFiles << "*.exome_depth.cnvs.tsv"
    if(callers.contains("xhmm"))
        sourceFiles << "*.xhmm_discover.xcnv"
    if(callers.contains("ex"))
        sourceFiles << "${batch_name}.excavator.cnvs.tsv"
    if(callers.contains("ang"))
        sourceFiles << "*.angel.cnvs.bed"
    if(callers.contains("truth"))
        sourceFiles << true_cnvs

    def output_prefix = (chr == "*" ? "" : "." + chr)
    from(sourceFiles) produce(batch_name + ".cnv${output_prefix}.summary.tsv", batch_name +"${output_prefix}.plotdata.RData") {
        R {"""
            source("$DSDSCRIPTS/cnv_utils.R")

            cnv.samples = c(${sample_names.collect{"'$it'"}.join(",")})

            target.bed = read.bed.ranges("$target_bed")

            cnv.results = list(
                xhmm = load_xhmm_results("$input.xcnv"),
                ex = load_excavator_results("$input2.tsv"),
                ed = load_exomedepth_results(file.name="$input1.tsv"),
                ang = load_angel_results(file.name="$input1.bed") """ + 
                    (simulation ? """, truth = load_angel_results(file.name="$input2.bed")""" : "") +
                """ 
            )
            cnv.all = combine_cnv_caller_results(cnv.samples, cnv.results, target.bed, chr="$chr")

            batch.cov = load.normalised.coverage(c(${all_samples.collect { "'$it.key'" }.join(",")}),
                                                 list(${all_samples.collect { "'${it.value.files.coverage[0]}'" }.join(",")}), 
                                                 chr="$chr")

            vcf_files = gsub(".vcf\$", ".vcf.bgz", c(${all_samples.collect{"'${it.value.files.vcf[0]}'"}.join(",")}))

            # Figure out the size of the chromosome from VCF header (assume all samples the same)
            hdr = scanBcfHeader(vcf_files[[1]])[[1]]
            chr.end = as.integer(hdr$Header$contig['$chr','length'])

            vcfs = sapply(cnv.samples, function(s) {
                v = readVcf(file=vcf_files[[which(cnv.samples==s)]],
                            genome='hg19', 
                            param=ScanVcfParam(which=GRanges(seqnames='$chr',ranges=IRanges(start=0, end=chr.end))))

                seqlevels(v,force=T) = hg19.chromosomes
                return(v)
            }, USE.NAMES=T)

            cnv.all = annotate_variant_counts(vcfs, cnv.all)

            cnv.output = data.frame(as.data.frame(cnv.all),
                vcf=sapply(cnv.all$sample, function(s) {vcf_files[[which(cnv.samples==s)]]})
            )

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
            if($create_plots && (length(cnv.all)>0)) {
                for(i in 1:length(cnv.all)) {
                  cat(sprintf("Producing plot for cnv %d (%s)\\n", i, cnv.all[i]$chrpos))
                  png(sprintf('report/cnv_%s_%d.rel.png',chr, i), width=1600, height=800)
                  tryCatch( plot.cnv(cnv.all[i], cnv.results, batch.cov, vcfs, add.title.text=F, cex.cnv.label=1.7),  error = function(e) print(e$message) )
                  dev.off()
                  png(sprintf('report/cnv_%s_%d.abs.png',chr, i), width=1600, height=800)
                  tryCatch(plot.cnv(cnv.all[i], cnv.results, batch.cov, vcfs, add.title.text=F, cex.cnv.label=1.7,scale="absolute"),  error = function(e) print(e$message))
                  dev.off()
                }
            }
        """}
    }
}

combine_cnv_report = {
    doc "Combine together CNV reports that are created for multiple chromosomes"
    // from("*.chr*.summary.tsv") { 
        exec """
            cat $inputs.summary.tsv | awk '{ if((NR==1) || (\$1!="seqnames")) print \$0; }' > ${output("report/" +batch_name+".cnv.summary.tsv")}
        """
    //}
}

cnv_report = {

    requires batch_name : "Name of batch for which report is being generated",
             callers : "Names of CNV callers in the CNV summary file" // see summarize_cnvs

    var inline_js : true,
        cnv_summary : input.cnv.summary.tsv,
        types : ["DEL","DUP"]

    // println "Inlining Javascript? " + inline_js.class.name
    produce("report/cnv_report.html", "cnv_report.zip") {
        send report('templates/cnv_report.html') to file: output.html
        exec """
                zip -r $output.zip report -x \\*.RData report
        """
    }
}


