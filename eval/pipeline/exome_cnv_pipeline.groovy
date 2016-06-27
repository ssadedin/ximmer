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

requires batches : """
                         Analysis batch identifiers for the run. This will be used to make directories 
                         that store all the analysis intermediate files and results. This allows
                         you to run different analyses and compare results from different 
                         settings, etc.
                      """,
         target_bed : """
                         BED file describing target region on which analysis is to be run. Usually 
                         this should be the whole region covered by probes / amplicons, not just 
                         the actual genes that are of interest, and should be unique
                      """,
         refgene :    """
                         A copy of the UCSC refgene database.
                      """
// Basic configuration
load 'config.groovy'

// A flag indicating whether this run is analysing simulated deletions or not.
// When analysing simulated data, only the X chromosome results are reported,
// and various other reports that are only relevant to simulations are created.
simulation = true

// For simulations we force chr=chrX
/*
if(simulation) {
    chr = "chrX"
    chrs = chr('X', filterInputs:false) 
}
else {
    chrs = chr(1..22,'X',filterInputs:false)
}
*/

println "Parsing samples from " + args

sample_info = SampleInfo.fromFiles(args)

// Can be overridden from command line
target_samples = sample_info.keySet()
if(target_samples instanceof String)
    run_samples = target_samples.split(",")*.trim()
else
    run_samples = target_samples

sample_names = run_samples

// The list of chromosomes to consider for analysis
// Not all of them will be analysed - they will be filtered
// by other checks, such as which of them overlaps the target 
// region, which have enough targets within them to be analysable,
// etc.
INCLUDE_CHROMOSOMES = (["chrX"] + (1..22).collect { "chr" + it })

// Only parallelise over the chromosomes actually in the target bed file
chromosomes = new BED(target_bed).load()*.chr.unique().grep { it in INCLUDE_CHROMOSOMES }

// If provided on command line as option
if(chromosomes instanceof String) {
    chromosomes = chromosomes.split(",") as List
}

load 'excavator.groovy'
load 'xhmm.groovy'
load 'exome_depth.groovy'
load 'cn_mops.groovy'
load 'conifer.groovy'
load 'summarize_cnvs.groovy'
load 'init_stages.groovy'

callers = "xhmm,ed,cnmops,truth"

cnv_callers = callers.split(",") as List

bpipe.Config.userConfig.autoFilter = "false"



set_sample = {
    branch.sample = branch.name
}

finish = { 
    exec "echo 'Finished'" 
}

batch_dirs = batches.tokenize(",")

init = { 
    
    branch.batch_name = branch.name
    
    branch.dir = batch_name 
    println "=" * 100
    println "Batch ${batch_name} Analysing ${sample_info.keySet().size()} samples:\n\n${sample_info.keySet().join('\n')}\n"
    println "=" * 100
    println "Chromosomes: " + chromosomes
    println "=" * 100
    
    load batch_name + '/caller.params.txt'
}


run {
    
    caller_stages = [ ]
    
    common_stages = [ ]

    if('ex' in cnv_callers) 
        caller_stages << (init_excavator + excavator_pipeline)

    if('ed' in cnv_callers) 
        caller_stages << (init_exome_depth + exome_depth_pipeline)

    if('xhmm' in cnv_callers)  {
        caller_stages << (init_xhmm + xhmm_pipeline)
        common_stages << "%.bam" * [ gatk_depth_of_coverage ]
    }

    if('cnmops' in cnv_callers)
        caller_stages << (init_cn_mops + cn_mops_call_cnvs)
        
    if('cfr' in cnv_callers)
        caller_stages << (init_conifer + run_conifer)

    common_stages + 
        batch_dirs * [
            init + create_analysable_target + caller_stages + create_cnv_report +
                 INCLUDE_CHROMOSOMES * [ touch_chr + plot_cnv_coverage ]  +
                 sample_names * [ extract_sample_files ] 
         ]
} 
