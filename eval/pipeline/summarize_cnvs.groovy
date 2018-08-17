// vim: expandtab:sw=4:ts=4:cindent


init_chr = {
   branch.chromosome = branch.name
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

             JAVA_OPTS="-Xmx1g -noverify" $GROOVY -cp $TOOLS/groovy-ngs-utils/1.0/groovy-ngs-utils.jar -e 'new RangedData("$input.tsv").load().each { println([it.chr, it.from-$slop, it.to+$slop].join("\\t"))  }' | 
                 $SAMTOOLS view -b -L - jinput.bam > $output.bam
        """
    }
}

plot_cnv_coverage = {
    requires target_bed : "Flattened, sorted BED file describing target regions, with ID column containing gene",
             refgene : "UCSC refGene database (usually named refGene.txt)"

    def chromosome = branch.name
    
    output.dir="$branch.dir/report"

    var reportSamples : false,
        draw_cnvs : true

    def reportSamplesFlag = reportSamples ? reportSamples.split(",")*.trim().collect { " -sample " + it }.sum()  : ""
    if(!draw_cnvs) {
        println "Skip drawing CNVs because disabled by setting"
        return
    }
    
    from("cnv_report.tsv", target_bed) { produce("cnv_${chr}_*.js") {

        def caller_opts = []

        batch_cnv_results.each { resultsEntry ->
            String caller = resultsEntry.key.tokenize('_')[0]
            String caller_label = resultsEntry.key
            caller_opts << "-$caller $caller_label:$resultsEntry.value"
        }
              
        if(simulation)  {
            caller_opts << "-generic truth:$input.true_cnvs.bed"
        }
                 
        exec """
            unset GROOVY_HOME 

            JAVA_OPTS="-Xmx8g -Djava.awt.headless=true -noverify" $GROOVY -cp $GNGS_JAR:$XIMMER_SRC $XIMMER_SRC/CNVDiagram.groovy
                -chr $chromosome
                -cnvs $input.tsv
                -ref $HGFA
                -gatkcov common/xhmm
                -targets $input.bed
                -json -nopng
                -o ${output.dir+"/cnv.png"} $reportSamplesFlag
                -t $threads ${caller_opts.join(" ")} ${inputs.vcf.withFlag("-vcf")} ${inputs.vcf.gz.withFlag("-vcf")} ${inputs.bam.withFlag("-bam")}
                -refseq $refgene
        ""","plot_cnv_coverage"
      }
    }
}

create_cnv_report = {
    
    requires refgene : 'Path to RefGene database',
             DGV_CNVS : 'Path to dgvMerged file downloaded from UCSC'

    var ([ angel_quality_threshold : 8.0f,
        batch_name : false,
        bam_file_path : "http://172.16.56.202/$batch_name/",
        sample_info : false,
        simulation: false,
        imgpath: false,
        genome_build : false,
        sample_id_mask : false,
        gene_filter: '',
        exclude_genes: '',
        minimum_category: false,
        sample_map: false,
        DDD_CNVS: false,
        file_name_prefix : "" ] + 
            batch_cnv_results*.key.collectEntries {  caller_label ->
                [ caller_label + '_quality_filter', false ]
            }
     ) 
    
    String refGeneOpts = ""
    if(genome_build != false) {
        refGeneOpts = "-refgene download -genome $genome_build"
    }
    else {
        refGeneOpts = "-refgene $refgene "
    }
    
    List qualityParams = []

    output.dir="$branch.dir/report"
    
    def caller_opts = []

    batch_cnv_results.each { resultsEntry ->
            String caller = resultsEntry.key.tokenize('_')[0]
            String caller_label = resultsEntry.key
            caller_opts << "-$caller $caller_label:$resultsEntry.value"
    }
    
    def geneFilterOpts = gene_filter ? " -genefilter $gene_filter " : ''
    def excludeGenesOpts = exclude_genes ? " -exgenes $exclude_genes " : ''
    def geneListOpts = genelists.collect { name, f -> 
        "-genelist $name=${file(f).absolutePath}"
    }.join(' ')
    
    def sampleMapParam = sample_map ? "-samplemap $sample_map" : "" 
    
    def minCatOpt = minimum_category ? "-mincat $minimum_category " : ""
    
    def idMaskOpt = sample_id_mask ? "-idmask '$sample_id_mask'" : ""
    
    def dddOpt = DDD_CNVS ? "-ddd $DDD_CNVS" : ""
        
    produce("${file_name_prefix}cnv_report.html", "${file_name_prefix}cnv_report.tsv", "${file_name_prefix}combined_cnvs.json") {

        def true_cnvs = ""
        if(simulation) {
            true_cnvs = "-truth $input.true_cnvs.bed"
        }
        
        exec """
            unset GROOVY_HOME

            JAVA_OPTS="-Xmx12g -noverify" $GROOVY -cp $GNGS_JAR:$XIMMER_SRC:$XIMMER_SRC/../resources:$XIMMER_SRC/../js $XIMMER_SRC/SummarizeCNVs.groovy
                -target $target_bed ${caller_opts.join(" ")} $refGeneOpts
                ${inputs.vcf.withFlag("-vcf")} ${inputs.vcf.gz.withFlag("-vcf")} -bampath "$bam_file_path"
                -tsv $output.tsv -json $output.json ${imgpath?"-imgpath "+imgpath.replaceAll('#batch#',batch_name):""}
                -dgv $DGV_CNVS $dddOpt $true_cnvs $idMaskOpt $geneFilterOpts $excludeGenesOpts $geneListOpts $minCatOpt $sampleMapParam
                ${batch_quality_params.join(" ")} -o $output.html 
                ${batch_name ? "-name $batch_name" : ""} ${inputs.bam.withFlag('-bam')}
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

calc_qc_stats = {
    
    var batch: 'ximmer',
        type: 'qc'
    
    output.dir = "common/$type"
    
    produce(["${batch}_per_base.coverage.tsv.gz", "${batch}.coeffv.js", "${batch}.correlations.js","${batch}.correlations.tsv", "${batch}.cov.js"]) {
        exec """
            set -o pipefail

            unset GROOVY_HOME

            $JAVA -Xmx4g -cp $GROOVY_ALL_JAR:$GNGS_JAR gngs.tools.MultiCov
                    -cvj $output.js
                    -stats 
                    -cv  
                    -corr .
                    -2pass
                    -covo $output.cov.js
                    -co $output.tsv
                    -co $output.correlations.js
                    -bed $target_bed $inputs.bam | gzip -c > $output.gz; 
        """, "calc_qc_stats"
    }
}
