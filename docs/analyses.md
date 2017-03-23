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

## Setting the BAM Files

The most important input to any CNV detection method is the BAM files to analyse. With 
Ximmer these are specified directly in the configuration file for each analysis. The 
specification can be either string value or a list of strings. Each string itself
can be a wildcard pattern to match multiple BAM files.

Example of single wildcard expression:

```
bam_files="/home/simon/data/bam_files/*.bam"
```

Example of multiple wildcard expressions:

```
bam_files=[
    "/home/simon/data/bam_files1/*.bam",
    "/home/simon/data/bam_files2/*.bam"
]
```

## Setting Sample Sex

The sex of each sample is important for two reasons: firstly, one of Ximmer's 
simulation methods utilises the sex of each sample directly in the simulation method.
However sex is important even when just doing analysis because of the differing ploidy 
of the X-chromosome between males and females. 

Ximmer can automatically detect the sex of samples, so specifying it is optional. However
the sex-detection algorithm takes some time to run, and can occasionally be inaccurate 
if data has unexpected characteristics. Therefore it's better to specify the sample 
sexes if you know them. There are two ways to do this. The first way, is to supply 
a [PED file](https://www.broadinstitute.org/haploview/input-file-formats) that specifies 
the sexes. This method is convenient when you already have such a file:

```
ped_file="/home/simon/data/samples.ped"
```

The second way is to explicitly list the males and females:

```
samples {
    males = [
        "SAMPLE_X123",
        "SAMPLE_X542"
        ...
    ]

    females = [
        "SAMPLE_X921",
        "SAMPLE_X291",
        ...
    ]
}
```

Note that in all cases, the samples specified must match the sample ids specified in the BAM 
files supplied.

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
settings to find out what works best on your particular data. To do that you need
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



## Specifying Variants

It can be informative to know when variants such as SNVs and indels overlap CNV calls. 
this is important for two reasons:

 * Heterozygosity and allele balance helps inform about whether a CNV call is accurate
 * Overlapping loss of function variants can form compound heterzygous configurations 
   that result in a complete loss of a gene.

Ximmer can incorporate variant calls for samples into the analysis. You can provide 
these by specifying a list of variant calls under the `variants` attribute in the 
configuration file:

```
 variants="/home/simon/bams/my_project_variants.vcf"
```

The variants attribute can also be specified as a list of VCF files:

```groovy
 variants=[
     "/home/simon/bams/sample1_variants.vcf",
     "/home/simon/bams/sample2_variants.vcf",
     ...
 ]
     
```

## Identity Masking

Sometimes you do not wish to display the full id that is attached to samples in the 
BAM and VCF files in your CNV report. This might be because the sample ids are very 
long and unwieldy, or it could also be because there is sensitive or private information 
in the ids. Ximmer supports a function to mask out portions of the actual sample ids
from being displayed in the report. The function allows the user to specify a regular 
expression which must have a single group (ie. section enclosed in parentheses). Ximmer
will display the section in parentheses only when showing sample ids in the report.

**Important: the full sample id is still accessible internally within the report. This 
feature does not provide security against a malicious user wishing to unmask identities. 
The underlying sample ids are readily accessible by use of Javascript and potentially 
may leak into some parts of the user interface as well.**

Example: Trim a trailing portion from sample ids in the form _SNN:

```
sample_id_mask="(.*)_S[0-9]*"
```

Example: Retain only a trailing SNN portion of the sample id:

```
sample_id_mask=".*(_S[0-9]*)"
```

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












