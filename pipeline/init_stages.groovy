import gngs.*

init_excavator = { branch.excavator_batch=batch_name }
init_xhmm = { branch.xhmm_batch_name=batch_name }
init_exome_depth = {  }
init_cn_mops = { }
init_conifer = { }

create_analysable_target = {
    doc """
          Different CNV callers have various different constraints on what they
          can handle. This stage filters the target regions by a variety of common
          rules to ensure that all the CNV callers included can analyse the data. Also 
          supported is an arbitrary BED file for excluding any regions specified by the user.
        """

    requires target_bed : "The target region BED file",
             sample_names : "The names of samples to analyse"
             
             
    var min_target_size : 30,
        exclude_regions : '',
        filter_target_regions : true

        
    if(!filter_target_regions) {
        from(target_bed) {
            println "Filtering of target regions is disabled: CNV callers will be passed the raw BED file: $input.bed"
            branch.analysable_target = input.bed
        }
    }
    else {
            
        def numSamples = sample_names.size()
    
        from(target_bed) filter('analysable') {
    
            groovy """
               import gngs.*
    
               INCLUDE_CHROMOSOMES = "${chromosomes.join(/,/)}".split(",")
    
               targetRegion = new BED("$input.bed", withExtra:true).load()
    
               excludeRegions = "$exclude_regions" ? new BED("$exclude_regions").load() : null
    
               chrs = targetRegion*.chr.unique()
    
               chrCounts = targetRegion.countBy { it.chr }
    
               filteredTargets = targetRegion.grep { 
                   (it.chr in INCLUDE_CHROMOSOMES) && 
                   (chrCounts[it.chr] > $numSamples) && 
                   (it.to - it.from > $min_target_size) && 
                   ((excludeRegions == null) || !it.overlaps(excludeRegions))
               } as Regions
    
               
               filteredTargets.save("$output.bed", sorted:true)
            """,config:"small"
       }
       branch.analysable_target = output.bed
    }
       
   // Overwrite global variable
   branch.analysable_chromosomes = new BED(branch.analysable_target.toString()).load()*.chr.unique().grep { it in INCLUDE_CHROMOSOMES }
   
}

select_controls = {
    
    doc "Selects a subset of controls to use from the supplied controls based on correlation with test samples"
    
    var control_correlation_threshold : 0.9,
        control_samples : false
        
    if(!control_samples) {
        
        println "All bams: $all_bams"
        
        println "sample names: $branch.sample_names"
        
        branch.filtered_bams = all_bams
        filtered_sample_names = branch.sample_names ?: sample_names
        println "No control samples are specified: skipping control selection, samples are: " + filtered_sample_names.join(',')
        return
    }

    println "Control samples are: " + control_samples

    List all_control_samples = control_samples
    List test_samples = branch.sample_names.grep { !(it in control_samples) }
    
    if(!test_samples) {
        fail "All ${branch.sample_names.size()} samples were listed as controls. Please identify one or more samples as a test sample"
    }

    println "Test samples out of $branch.sample_names are: " + test_samples
    
    produce('filtered_controls.txt') {
        exec """
            JAVA_OPTS="-Xmx8g -Djava.awt.headless=true -noverify" $GROOVY -cp $GNGS_JAR:$XIMMER_SRC $XIMMER_SRC/ximmer/FilterControls.groovy
                -corr $input.correlations.js
                -thresh $control_correlation_threshold ${control_samples.collect { '-control ' + it}.join(' ')}
                > $output.txt
        ""","local"
    }
        
    
    List filtered_control_samples = file(output.txt).readLines()*.trim()
    
    List all_remaining_samples = filtered_control_samples + test_samples

    println "All remaining samples are: " + all_remaining_samples.join(',')
    branch.filtered_bams = all_bams.grep {
        new gngs.SAM(it).samples[0] in all_remaining_samples
    }
    
    branch.sample_info = sample_info.grep { it.key in filtered_control_samples }.collectEntries()
    
    branch.sample_names = branch.sample_names.grep { (it in all_remaining_samples) }
    
    filtered_sample_names = branch.sample_names.grep { (it in test_samples) || (it in filtered_control_samples) }

    println "Sample names after filtering are : $branch.sample_names"
    
    def stats_files = branch.sample_names.collect { s -> sample_to_control_cov_files[s] }
    
    forward(filtered_bams + stats_files.flatten())
}


