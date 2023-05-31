# Ximmer Simulations

## Introduction

This section explains the steps to add a new CNV caller to Ximmer.

Please note that this topic is more advanced and requires knowledge of
how to edit Bpipe scripts.

## Add Pipeline Stage to Perform Analysis

Create a new file with a name like `<caller name>.groovy`, and place it
in pipeline directory. In this file, define a pipeline stage that performs the 
analysis based on:

 * input BAM files provided as `$input.bam`
 * target regions provided as `$input.bed`

If you want to split the analysis by chromosome, you can do so by referencing the
`$chr` variable. If you do this, you need to add the splitting based on the $chromosomes
variable (see examples in other stages that split by chromosome, for example, exome_depth).

Once you have defined your necessary pipeline stages in their own file, 
ensure that file is loaded in `exome_cnv_pipeline.groovy`, and add the pipeline stage
in the section where CNV callers are defined further down:

```groovy
   caller_pipelines = [
       ex  :  (init_excavator + excavator_pipeline),
       
       ed  : (init_exome_depth + exome_depth_pipeline),
       
       xhmm: (init_xhmm + xhmm_pipeline),
       
       cnmops: (init_cn_mops + cn_mops_call_cnvs),
       
       cfr:  (init_conifer + run_conifer),
       
       cdx : (init_codex + codex_pipeline)

       // Add your stages in this map here
    ]    
```

## Add a New Class for Parsing the Results

This should be a groovy class in src/main/groovy, and it should extend from CNVResults. 
This should not require a lot of code, especially if using the RangedData class to 
do the loading. One thing you must do, however, is nominate a field to use as a quality 
score for ranking CNVs.
