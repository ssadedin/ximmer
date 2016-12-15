// vim: expandtab:sw=4:ts=4:cindent


touch_chr = {

   branch.chromosome = branch.name

   exec """
        echo $chr > $output.txt
   """
}

extract_sample_files = {
    requires sample_info : "Sample meta data object"
    var sample : branch.name
    branch.sample = sample
    branch.sample_info = sample_info
    println "Processing files for sample $sample : " + branch.sample_info[sample].files.all
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

             JAVA_OPTS="-Xmx1g" $GROOVY -cp $TOOLS/groovy-ngs-utils/1.0/groovy-ngs-utils.jar -e 'new RangedData("$input.tsv").load().each { println([it.chr, it.from-$slop, it.to+$slop].join("\\t"))  }' | 
                 $SAMTOOLS view -b -L - $input.bam > $output.bam
        """
    }
}

plot_cnv_coverage = {
    requires target_bed : "Flattened, sorted BED file describing target regions, with ID column containing gene",
             refgene : "UCSC refGene database (usually named refGene.txt)"

    output.dir="$branch.dir/report"

    var reportSamples : false,
        draw_cnvs : true

    def reportSamplesFlag = reportSamples ? reportSamples.split(",")*.trim().collect { " -sample " + it }.sum()  : ""
    if(!draw_cnvs) {
        println "Skip drawing CNVs because disabled by setting"
        return
    }
    
    uses(threads:2..8) {
        from("cnv_report.tsv") { produce("cnv_${chr}_*.png") {

            def caller_opts = []

            /*
            if('xhmm' in cnv_callers)
               caller_opts << "-xhmm $input.xcnv"

            if('ed' in cnv_callers)
               caller_opts << "-ed $input.exome_depth.cnvs.tsv"

            if('mops' in cnv_callers)
               caller_opts << "-cnmops $input.cnmops.cnvs.tsv"

            if('angelhmm' in cnv_callers) 
               caller_opts << "-angel $input.angelhmm.cnvs.bed"

            if('cfr' in cnv_callers) 
               caller_opts << "-cfr $input.conifer.cnvs.tsv"
               
            */
               
            batch_cnv_results.each { resultsEntry ->
                String caller = resultsEntry.key.tokenize('_')[0]
                String caller_label = resultsEntry.key
                caller_opts << "-$caller $caller_label:$resultsEntry.value"
            }
              
            if(simulation) 
                caller_opts << "-generic truth:$input.true_cnvs.bed"
                 
            exec """
                unset GROOVY_HOME 

                JAVA_OPTS="-Xmx8g -Djava.awt.headless=true" $GROOVY -cp $GNGS_JAR:$XIMMER_SRC $XIMMER_SRC/CNVDiagram.groovy
                    -chr $chromosome
                    -cnvs $input.tsv
                    -targets $input.bed
                    -o ${output.dir+"/cnv.png"} $reportSamplesFlag
                    -t $threads ${caller_opts.join(" ")} ${inputs.vcf.withFlag("-vcf")} ${inputs.bam.withFlag("-bam")}
                    -refseq $refgene
            ""","plot_cnv_coverage"
        }
      }
    }
}

create_cnv_report = {

    var angel_quality_threshold : 8.0f,
        batch_name : false,
        bam_file_path : "http://172.16.56.202/$batch_name/",
        sample_info : false,
        simulation: false,
        imgpath: false

    List vcfs = []
    /*
    if(sample_info) 
        vcfs = sample_info.collect{it.value.files.vcf[0]}
    else
        vcfs = all_samples.collect{it.value.files.vcf[0]}
    */

    output.dir="$branch.dir/report"
    
    def caller_opts = []

    batch_cnv_results.each { resultsEntry ->
            String caller = resultsEntry.key.tokenize('_')[0]
            String caller_label = resultsEntry.key
            caller_opts << "-$caller $caller_label:$resultsEntry.value"
    }
        
       
    produce("cnv_report.html", "cnv_report.tsv") {

        def true_cnvs = ""
        if(simulation) 
            true_cnvs = "-truth $input.true_cnvs.bed"
        
        exec """
            JAVA_OPTS="-Xmx8g -noverify" $GROOVY -cp $GNGS_JAR:$XIMMER_SRC:$XIMMER_SRC/../resources $XIMMER_SRC/SummarizeCNVs.groovy
                -target $target_bed ${caller_opts.join(" ")}
                ${inputs.snpeff.vcf.withFlag("-vcf")} -bampath "$bam_file_path"
                -tsv $output.tsv ${imgpath?"-imgpath "+imgpath.replaceAll('#batch#',batch_name):""}
                -o $output.html ${batch_name ? "-name $batch_name" : ""} ${inputs.bam.withFlag('-bam')}
                -dgv $DGV_CNVS $true_cnvs
        """, "create_cnv_report"
    }
}

zip_cnv_report = {
    produce("cnv_report.zip") {
        exec """
            zip -r $output.zip $branch.dir/report -x \\*.RData report
        """
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

