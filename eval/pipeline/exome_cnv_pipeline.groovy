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

exome_depth_split_chrs = true

bpipe.Config.userConfig.autoFilter = "false"

tmpFile = new File("tmpdata")
tmpFile.mkdirs()
TMPDIR=tmpFile.path

batch_dirs = batches.tokenize(",")

load 'ximmer_pipeline_core.groovy'

run { ximmer_core } 
