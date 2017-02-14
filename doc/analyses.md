# Ximmer Analyses

## Introduction

Ximmer configuration has two basic parts:

 * how CNVs are simulated ("simulation")
 * how CNVs callers are run ("analyses")

You can choose to do either or both of simulation and analysis. If you just want some 
CNV results from one or more CNV callers, just do analysis. If you just want to simulate 
some CNVs and then take the data to analyse separately, just do the simulation part. If 
you want everything, then do both.

For the analysis part, Ximmer is designed to allow you to easily run many different callers 
with many different settings, or the same caller with many different settings, and compare 
all the results together. Each separate set of settings applied across a selection 
of CNV callers is called an "analysis". 

Whether you're doing simulation or analysis or both, the starting point is to create 
a configuration file that controls the analysis. This file can have a lot of settings,
but in its most minimal form, all it needs is:

 * bam_files setting describing where the BAM files to analyse or simulate from are
 * target_regions setting describing the capture region
 * callers section describing which CNV callers to run (see below)

## Default Analysis

Some CNV callers have a lot of adjustable parameters. Therefore it is inconvenient to
have to set every parameter every time you run them. For this reason, Ximmer has a 
section that allows you to define the default parameters. Then any configuration you make is
simply overriding the defaults for only the parameters you are interested in modifying. Each
separate set of modified parameters is called an "analysis". The set of default parameters 
is the "default analysis". The default analysis is configured in the "callers" section of the
configuration file. An example is shown below:

```
callers {
    xhmm {
        exome_wide_cnv_rate=1e-04
        mean_number_of_targets_in_cnv=3
    }
    exomedepth {
        transition_probability=0.0001
    }
    cnmops {
        prior_impact=5
        min_width=1
    }
    conifer {
        conifer_svd_num=1
    }
}
```

If you don't configure anything else, only the default analysis is what
will be run to find CNVs in the data.

## Customized Analyses

A very common task in using a CNV detection tool is to try different 
setting sto find out what works best on your particular data. To do that you need
to compare between different settings for the same tool. Each group of settings
that you wish to run with is called an "analysis". These are configured in the
analysis section of the configuration file. An example is below:

```
analyses {

    'xhmmtune' {
        xhmm_1 { 
            exome_wide_cnv_rate=1e-02 
        }
        xhmm_2 { 
            exome_wide_cnv_rate=1e-02; 
            xhmm_pve_mean_factor=0.2; 
        }
        xhmm_3 { 
            exome_wide_cnv_rate=1e-04; 
        }
    }
}
```

This example configuration only runs XHMM. The name for the analysis is `xhmmtune`. The use of XHMM 
is *inferred* from the prefix `xhmm_` for the label of each individual block 
within the analyses. The configuration parameters themselves are specified within each 
block and are specific to each caller (see table TODO).

TODO:

| Caller | Parameter |  Description | Example |
|--------|-----------|--------------|---------|
| foo   | bar   | frog | house |


## Full Configuration Example

Below is a very simple, minimal but working configuration which analyses
a set of BAM files to find CNVs using XHMM and ExomeDepth using their
default settings:

```
bam_files="/home/simon/data/*.bam"
target_regions="/home/simon/data/EXOME.bed"
concurrency=20
callers {
    xhmm {}
    exomedepth {}
}
```










