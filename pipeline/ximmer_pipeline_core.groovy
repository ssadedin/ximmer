
var enable_kmer_normalisation : false

set_sample = {
    branch.sample = branch.name
}

finish = { 
    exec "echo 'Finished'" 
}

init_batch = { 
    
    branch.batch_name = branch.name
    
    // Map of caller label (eg: "xhmm_1") to the result file for that label
    branch.batch_cnv_results = [:]
    branch.batch_quality_params = []
    
    branch.dir = batch_name 
    println "=" * 100
    println "Batch ${batch_name} Analysing ${target_samples.size()} samples:\n\n${target_samples.join(',')}\n"
    println "=" * 100
    println "Chromosomes: " + chromosomes
    println "=" * 100
    
    List<String> param_files = file(batch_name).listFiles().grep { it.name.endsWith('.params.txt') }*.path
    
    if(param_files.isEmpty())
        throw new bpipe.PipelineError("No cnv caller parameter files were found for batch $batch_name: this means no analysis will be done, and probably indicates a configuration error.\n\nPlease check you have specified CNV callers correctly.")

    println "Parameter files for $batch_name are $param_files"
    
    forward param_files
}

init_caller_params = {
    
    requires caller: 'A caller should be passed by the using() function'
    
    branch.params = branch.name
    
    branch.caller_label = (params == "default") ? caller : "${caller}_${params}"
    
    branch.dir = branch.dir + "/$caller_label"
    
    println "Using parameters " + params + " for $caller with output to " + branch.dir
    
    branch.exome_depth_split_chrs = true
    
    def paramsPath = "$batch_name/${caller}.${params}.params.txt"
    load(paramsPath)
    
    var exclude_samples: []
    
    String qualPropertyKey = caller_label + "_quality_filter"
    
    Map qualProps = [:]
    qualProps[qualPropertyKey] = false
    var(qualProps)
    
    if(delegate.getProperty(qualPropertyKey)) {
        def qualFilterValue = getProperty(caller_label + '_quality_filter')
        println "Setting quality filter for $caller to $qualFilterValue"
        batch_quality_params << " -quality ${caller_label}:$qualFilterValue"
    }
    else {
        println "No quality parameters found for $caller under key " + caller_label + '_quality_filter'
    }
    
    if(branch.exclude_samples){
        println "Excluding samples: " + branch.exclude_samples.join(',')
        forward(inputs.bam.grep { bamFile -> !(new gngs.SAM(bamFile).samples[0] in branch.exclude_samples) })
    }
    else {
        println "No samples are excluded for $caller_label: "
    }
}

register_caller_result = {
    
    batch_cnv_results[caller_label] = caller_result
}


init = {
    branch.all_bams = inputs.bam.collect { it.toString() } // clone
}

reset_bams = {
    println "Resetting BAMs to: $filtered_bams"
    forward(filtered_bams)
}

compute_kmer_profiles = segment {
    '%.bam' * [ calc_kmer_profile] + merge_kmer_profiles
}

ximmer_core = segment {
    
    def cnv_reports = []
    if(simulation) {
        cnv_reports << create_cnv_report
    }
    
    cnv_reports << create_cnv_report.using(file_name_prefix:"local_", imgpath: "") 
    
    caller_pipelines = [

       ex  :  (init_excavator + excavator_pipeline),
       
       ed  : (init_exome_depth + exome_depth_pipeline),
       
       xhmm: (init_xhmm + xhmm_pipeline),
       
       cnmops: (init_cn_mops + cn_mops_call_cnvs),
       
       cfr:  (init_conifer + run_conifer),
       
       cdx : codex_pipeline,
       
       dfn: delfin,
       
       savvy: savvy_cnv,
    ]  

    caller_stages = cnv_callers.collect { caller ->
        (caller + '.%.params.txt') * [ init_caller_params.using(caller:caller) + caller_pipelines[caller] + register_caller_result ]
    }
    
    println "There are ${caller_stages.size()} cnv caller stages to run in batches: $batch_dirs"
    
    init + create_analysable_target + compute_kmer_profiles.when { enable_kmer_normalisation } + '%.bam' * [ calc_target_covs ] + forward_all_cov_files + calc_combined_correlations.using(type:'rawqc') + select_controls + [
        calc_combined_correlations.using(type:'qc') + // + calc_qc_stats.using(type:'qc'),
        batch_dirs * [
            init_batch + caller_stages + reset_bams +
                 cnv_reports +
                 chromosomes * [ plot_cnv_coverage ]  
         ]
     ]
} 
