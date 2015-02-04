// vim: ts=4:expandtab:sw=4:cindent
////////////////////////////////////////////////////////////////////
//
// Integrated Exome / Targeted Sequencing CNV Calling Pipeline
//
// Runs multiple CNV callers in parallel, then combines the results
// to a common format and generates an HTML format report summarizing
// the results, including levels of concordance, support from 
// heterozygosity information and other annotation based information.
//
///////////////////////////////////////////////////////////////////

inputs "bam" : "Bam files for analysis"

requires batch_name : """
                         A unique identifier for the run. This will be used to make a directory 
                         that stores all the analysis intermediate files and results. This allows
                         you to run different analyses and compare results from different 
                         settings, etc.
                      """,
         target_bed : """
                         BED file describing target region on which analysis is to be run. Usually 
                         this should be the whole region covered by probes / amplicons, not just 
                         the actual genes that are of interest, and should be unique
                      """,
         dgv_file :   """
                         File downloaded from UCSC, containing population
                         frequencies of structural variants. This is the dgvMerged file / table,
                         usually called something like 'dgvMerged.txt'.
                      """

// Basic configuration
load 'config.groovy'

// A flag indicating whether this run is analysing simulated deletions or not.
// When analysing simulated data, only the X chromosome results are reported,
// and various other reports that are only relevant to simulations are created.
simulation = true

// For simulations we force chr=chrX
if(simulation) {
    chr = "chrX"
    chrs = chr('X', filterInputs:false) 
}
else {
    chrs = chr(1..22,'X',filterInputs:false)
}

sample_info = false

// If sample_info is provided as a file, then use that
if(sample_info) {
    println "Parsing from $sample_info"
    all_samples = SampleInfo.parse_sample_info(sample_info) 
    println "All samples are " + all_samples
    println "Bam files are " + all_samples*.value*.files*.bam
    if(args)
        args += all_samples*.value*.files*.bam.flatten()
    else
        args = all_samples*.value*.files*.bam.flatten()
}
else {
    all_samples = SampleInfo.fromFiles(args)
}

// Can be overridden from command line
sample_info = all_samples
target_samples = all_samples.keySet()
if(target_samples instanceof String)
    run_samples = target_samples.split(",")*.trim()
else
    run_samples = target_samples

sample_names = run_samples

println "Run samples are $run_samples"

load 'excavator.groovy'
load 'xhmm.groovy'
load 'exome_depth.groovy'
load 'cn_mops.groovy'
load 'summarize_cnvs.groovy'

// If not overridden by command line, assume all the callers are to be run
callers = "xhmm,ed,mops,truth"

cnv_callers = callers.split(",") as List

bpipe.Config.userConfig.autoFilter = "false"

set_sample = {
    branch.sample = branch.name
}

finish = { 
    exec "echo 'Finished'" 
}

// Initialization stages - these just set specific output directories for each 
// analysis tool
init_excavator = { branch.dir="$batch_name/excavator"; branch.excavator_batch=batch_name }
init_xhmm = { branch.dir="$batch_name/xhmm"; branch.xhmm_batch_name=batch_name }
init_exome_depth = { branch.dir="$batch_name/exome_depth" }
init_exome_copy = { branch.dir="$batch_name/exome_copy" }
init_cn_mops = { branch.dir="${batch_name}/cn_mops" }
init = { 
    branch.dir = batch_name 
    println "=" * 100
    println "Analysing samples: $all_samples"
    println "=" * 100
}

run {

    caller_stages = [ ]

    if('ex' in cnv_callers) 
            caller_stages << (init_excavator + excavator_pipeline)

    if('ed' in cnv_callers) 
            caller_stages << (init_exome_depth + exome_depth_pipeline)

    if('xhmm' in cnv_callers) 
            caller_stages << (init_xhmm + xhmm_pipeline)

    if('mops' in cnv_callers)
        caller_stages << (init_cn_mops + cn_mops_call_cnvs)

    init + caller_stages + 
         chrs * [ touch_chr + summarize_cnvs ] + combine_cnv_report +
         sample_names * [ extract_sample_files +  extract_cnv_regions ] + 
         cnv_report.using(callers:cnv_callers) 
} 
