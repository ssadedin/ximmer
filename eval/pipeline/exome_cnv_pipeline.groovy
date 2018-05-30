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



println "=" * bpipe.Config.config.columns 

println """
  | |/ /(_)___ ___  ____ ___  ___  _____
  |   // / __ `__ \\/ __ `__ \\/ _ \\/ ___/
 /   |/ / / / / / / / / / / /  __/ /    
/_/|_/_/_/ /_/ /_/_/ /_/ /_/\\___/_/     

Ximmer Integrated CNV Analysis Pipeline 
"""

println "=" * bpipe.Config.config.columns 

// println "Parsing samples from " + args

sample_info = SampleInfo.fromFiles(args)

genelist_ids = this.binding.parameters.grep { it.startsWith('genelist:') }.collect { it.replaceAll('^genelist:','') }
genelists = [:]
for(gl in genelist_ids) {
    genelists[gl] = this['genelist:'+gl]
    if(!file(genelists[gl]).exists())
        throw new IllegalArgumentException("File provided for gene list $gl does not exist: ${genelists[gl]}")
}

if(!genelists.isEmpty())
    println "The following genelists were provided: $genelist_ids"

// Can be overridden from command line
target_samples = sample_info.keySet()
if(target_samples instanceof String)
    run_samples = target_samples.split(",")*.trim()
else
    run_samples = target_samples

sample_names = run_samples

target_regions = new gngs.BED(target_bed).load()

// Sniff the target BED file to see if there is a 'chr' prefix or not
chr_prefix = target_regions[0].chr.startsWith('chr') ? "chr" : "" 

// The list of chromosomes to consider for analysis
// Not all of them will be analysed - they will be filtered
// by other checks, such as which of them overlaps the target 
// region, which have enough targets within them to be analysable,
// etc.
INCLUDE_CHROMOSOMES = ([chr_prefix + "X"] + (1..22).collect { chr_prefix + it })

// Only parallelise over the chromosomes actually in the target bed file
println "Scanning target region ..."
chromosomes = target_regions*.chr.unique().grep { it in INCLUDE_CHROMOSOMES }

// If provided on command line as option
if(chromosomes instanceof String) {
    chromosomes = chromosomes.split(",") as List
}

if(chromosomes.isEmpty())
    throw new RuntimeException("No entries were found in the configured BED file ($target_bed) corresponding to the configured chromosomes ($chromosomes, $INCLUDE_CHROMOSOMES)!")
    
// NOTE: this global variable is modified by the create_analysable_target pipeline stage!
analysable_chromosomes = chromosomes.clone()
    
println "Chromosomes for analysis are: $chromosomes"

load 'excavator.groovy'
load 'xhmm.groovy'
load 'exome_depth.groovy'
load 'cn_mops.groovy'
load 'conifer.groovy'
load 'codex.groovy'
load 'summarize_cnvs.groovy'
load 'init_stages.groovy'

callers = "xhmm,ed,cnmops,truth"

cnv_callers = callers.split(",") as List

bpipe.Config.userConfig.autoFilter = "false"

tmpFile = new File("tmpdata")
tmpFile.mkdirs()
TMPDIR=tmpFile.path

set_sample = {
    branch.sample = branch.name
}

finish = { 
    exec "echo 'Finished'" 
}

batch_dirs = batches.tokenize(",")

init_batch = { 
    
    branch.batch_name = branch.name
    
    // Map of caller label (eg: "xhmm_1") to the result file for that label
    branch.batch_cnv_results = [:]
    branch.batch_quality_params = []
    
    branch.dir = batch_name 
    println "=" * 100
    println "Batch ${batch_name} Analysing ${sample_info.keySet().size()} samples:\n\n${sample_info.keySet().join('\n')}\n"
    println "=" * 100
    println "Chromosomes: " + chromosomes
    println "=" * 100
    
    List<String> param_files = file(batch_name).listFiles().grep { it.name.endsWith('.params.txt') }*.path
    
    println "Parameter files for $batch_name are $param_files"
    
    forward param_files
}

init_caller_params = {
    
    requires caller: 'A caller should be passed by the using() function'
    
    branch.params = branch.name
    
    branch.caller_label = (params == "default") ? caller : "${caller}_${params}"
    
    branch.dir = branch.dir + "/$caller_label"
    
    println "Using parameters " + params + " for $caller with output to " + branch.dir
    
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

init_common = {
    println "Running Common Stages"
}

all_bams = []

init = {
    all_bams = inputs.bam.collect { it } // clone
}

reset_bams = {
    forward(all_bams)
}

run {
    
    common_stages = [ init_common, "%.bam" * [ gatk_depth_of_coverage  ] ]
       
    caller_pipelines = [
       ex  :  (init_excavator + excavator_pipeline),
       
       ed  : (init_exome_depth + exome_depth_pipeline),
       
       xhmm: (init_xhmm + xhmm_pipeline),
       
       cnmops: (init_cn_mops + cn_mops_call_cnvs),
       
       cfr:  (init_conifer + run_conifer),
       
       cdx : codex_pipeline
    ]    

    caller_stages = cnv_callers.collect { caller ->
        (caller + '.%.params.txt') * [ init_caller_params.using(caller:caller) + caller_pipelines[caller] + register_caller_result ]
    }
    
    init + create_analysable_target + common_stages + 
        batch_dirs * [
            init_batch + caller_stages + reset_bams +
                 create_cnv_report + create_cnv_report.using(file_name_prefix:"local_", imgpath: "") +
                 INCLUDE_CHROMOSOMES * [ plot_cnv_coverage ]  +
                 sample_names * [ extract_sample_files ] 
         ]
} 
