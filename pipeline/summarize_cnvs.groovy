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

    var autoFilterCNVDiagrams : false
    
    def chromosome = branch.name
    
    output.dir="$branch.dir/report"

    var reportSamples : false,
        draw_cnvs : true

    def reportSamplesFlag = reportSamples ? reportSamples.split(",")*.trim().collect { " -sample " + it }.sum()  : ""
    if(!draw_cnvs) {
        println "Skip drawing CNVs because disabled by setting"
        return
    }
    
    def autoFilterOption = ''
    if(autoFilterCNVDiagrams) {
        autoFilterOption = '-autoFilter'
    }
    
    from("cnv_report.tsv", target_bed) { produce("cnv_diagrams_${chromosome}_created.txt") {

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
                -covjs $input.cov.js
                -targets $input.bed $autoFilterOption
                -json -nopng
                -o ${output.dir+"/cnv.png"} $reportSamplesFlag
                -t $threads ${caller_opts.join(" ")} ${inputs.vcf.withFlag("-vcf")} ${inputs.vcf.gz.withFlag("-vcf")} ${inputs.bam.withFlag("-bam")}
                -refseq $refgene 

            ls cnv_${chromosome}_* | wc > $output.txt

        ""","plot_cnv_coverage"
      }
    }
}

create_cnv_report = {
    
    requires refgene : 'Path to RefGene database',
             DGV_CNVS : 'Path to dgvMerged file downloaded from UCSC'

    var ([ 
            angel_quality_threshold : 8.0f,
            batch_name : false,
            bam_file_path : false,
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
            file_name_prefix : "",
            mergeOverlapFraction: 0.4,
            cnvMergeMode: "sharedtargets",
            control_samples : false,
            report_chr : false
        ] + 
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
    
    def reportChrFlag = report_chr ? "-chr $report_chr" : ""
    
    String sampleMapParam = sample_map ? "-samplemap $sample_map" : "" 
    
    String minCatOpt = minimum_category ? "-mincat $minimum_category " : ""
    
    String idMaskOpt = sample_id_mask ? "-idmask '$sample_id_mask'" : ""
    
    String dddOpt = DDD_CNVS ? "-ddd $DDD_CNVS" : ""

    String samplesOption = ""
    if(control_samples) {
        def test_samples = sample_names.grep { !(it in control_samples) }
        samplesOption = "-samples ${test_samples.join(',')}"
    }
    
    def dgvFlag = DGV_CNVS ? "-dgv $DGV_CNVS" : ""
    
    produce("${file_name_prefix}cnv_report.html", "${file_name_prefix}cnv_report.tsv", "${file_name_prefix}combined_cnvs.json") {

        def true_cnvs = ""
        if(simulation) {
            true_cnvs = "-truth $input.true_cnvs.bed"
        }
        
        exec """
            unset GROOVY_HOME

            JAVA_OPTS="-Xmx12g -noverify" $GROOVY -cp $GNGS_JAR:$XIMMER_SRC:$XIMMER_SRC/../resources:$XIMMER_SRC/../js $XIMMER_SRC/SummarizeCNVs.groovy
                -target $target_bed ${caller_opts.join(" ")} $refGeneOpts $reportChrFlag ${inputs.vcf.withFlag("-vcf")} ${inputs.vcf.gz.withFlag("-vcf")}
                -tsv $output.tsv -json $output.json ${imgpath?"-imgpath "+imgpath.replaceAll('#batch#',batch_name):""} -mergefrac $mergeOverlapFraction
                -mergeby $cnvMergeMode $dgvFlag $dddOpt $true_cnvs $idMaskOpt $geneFilterOpts $excludeGenesOpts $geneListOpts $minCatOpt $sampleMapParam $samplesOption
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

calc_target_covs = {

    var coverage_cache_dir : false

    output.dir = "common/rawqc/individual"

    def sample = new gngs.SAM(input.bam.toString()).samples[0]
    def intervalSummaryPath = sample + '.calc_target_covs.sample_interval_summary'
    def statsPath = sample + '.stats.tsv'

    if(coverage_cache_dir) {
        def cachedPath = new File(coverage_cache_dir, intervalSummaryPath) 
        def statsFile = new File(coverage_cache_dir, statsPath) 
        if(cachedPath.exists()) {

            forward(cachedPath.absolutePath, statsFile.absolutePath) 

            sample_to_control_cov_files[sample] = [statsFile.absolutePath, cachedPath.absolutePath]
                /*
            exec """
                ln -s $cachedPath.absolutePath $output.sample_interval_summary

                ln -s $statsFile.absolutePath $output.tsv
            ""","local"
            */
            println "Using cached path $cachedPath for coverage values for $sample"
            return
        }
    }

    produce(statsPath, intervalSummaryPath) {

        exec """
            unset GROOVY_HOME;  

            $JAVA -Xmx${memory}g -cp $GROOVY_ALL_JAR:$GNGS_JAR
                gngs.tools.Cov 
                -L $target_bed
                -o /dev/null 
                -samplesummary $output.stats.tsv
                -intervalsummary $output.sample_interval_summary
                $input.bam
        """, "calc_single_cov"
    }
    
    sample_to_control_cov_files[sample] = [output.stats.tsv, output.sample_interval_summary]
}

forward_all_cov_files = {
    forward(sample_to_control_cov_files*.value.flatten())
}

calc_combined_correlations = {
    
    var batch: 'ximmer',
        type: 'qc'

    output.dir = "common/$type"
    
    uses(threads:1..4) {
        produce([batch +'.combined.correlations.tsv', batch + '.combined.correlations.js', batch + '.combined.cov.js', batch + '.combined.coeffv.js', batch + '.combined.sample_interval_summary']) {
            exec """
                JAVA_OPTS="-Xmx15g -Djava.awt.headless=true -noverify" $GROOVY -cp $GNGS_JAR:$XIMMER_SRC $XIMMER_SRC/ximmer/CalculateCombinedStatistics.groovy
                -corrTSV $output1
                -corrJS $output2
                -covJS $output3
                -coeffvJS $output4
                -stats $output5
                -threads $threads
                $inputs.sample_interval_summary
                $inputs.stats.tsv
            """, "calc_combined_correlations"
        }
    }
}

calc_qc_stats = {
    
    var batch: 'ximmer',
        type: 'qc'
        
    
    output.dir = "common/$type"
    
    def relFlag = type == 'rawqc' ? '' : '-rel'
    
    produce(["${batch}_per_base.coverage.tsv.gz", "${batch}.coeffv.js", "${batch}.correlations.js","${batch}.correlations.tsv", "${batch}.cov.js", "${batch}.merge.sample_interval_summary"]) {
        exec """
            set -o pipefail

            unset GROOVY_HOME

            $JAVA -Xmx${memory}g -cp $GROOVY_ALL_JAR:$GNGS_JAR gngs.tools.MultiCov
                    -cvj $output.js
                    -stats 
                    -cv  
                    -corr . $relFlag
                    -2pass
                    -targetmeans $output.sample_interval_summary
                    -covo $output.cov.js
                    -co $output.tsv
                    -co $output.correlations.js
                    -bed $target_bed $inputs.bam | gzip -c > $output.gz; 
        """, "calc_qc_stats"
    }
}

convert_to_vcf = {
    output.dir = "vcfs"
    
    branch.sample = branch.name
    
    from('local_cnv_report.tsv') produce(sample + '.cnv.vcf') {
        exec """
            unset GROOVY_HOME

            JAVA_OPTS="-Xmx12g -noverify" $GROOVY -cp $GNGS_JAR:$XIMMER_SRC:$XIMMER_SRC/../resources:$XIMMER_SRC/../js $XIMMER_SRC/ximmer/TSVtoVCF.groovy
                -i $input.tsv
                -s $sample
                -r $HGFA
                -o $output.vcf
        ""","convert_to_vcf"
    }
}
