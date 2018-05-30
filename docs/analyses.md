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

```groovy
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
block and are specific to each caller (see table).

| Caller     | Parameter              | Description                                   | Example / Default |
|------------|------------------------|-----------------------------------------------|-------------------|
| ExomeDepth | transition_probability |                                               | 10e-4             |
|            | expected_cnv_length    |                                               | 50000             |
|            |                        |                                               |                   |
| XHMM       | exome_wide_cnv_rate    |                                               | 10e-8             |
|            | xhmm_pve_mean_factor   |                                               | 0.7               |
|            |                        |                                               |                   |
| Conifer    | conifer_svd_num        |                                               |                   |
|            | conifer_call_threshold |                                               |                   |
|            |                        |                                               |                   |
| cn.MOPs    | prior_impact           | Weighting of prior probability of CNV         | 10                |
|            | min_width              | Min target regions to call a CNV              | 5                 |
|            | lower_threshold        | Affects threshold on coverage for CNV calling | -0.8              |
|            | panel_type             | Sets a range of parameters for panel vs exome | exome (or blank)  |
|            |                        |                                               |                   |
| CODEX      | k_offset               | Adjusts CODEX's preferred k by given amount   | 0                 |
|            | max_k                  | Sets the maximum value of k to be tried       |                   |

## Filtering by quality

If a CNV caller produces many false positives, you may wish to filter out results that 
have a low quality score assigned by the caller. You can specify a quality score threshold
using the `quality_filter` setting for each analysis. Note that how this is interpreted 
is specific to each CNV caller. Ximmer uses a fixed quality metric for each CNV caller
(see publication). 

Example:

```groovy
    cnmops {
        prior_impact=5
        quality_filter=1.5
    }
```


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

Note that each entry in either of these forms can be a Unix style "glob" to match
multiple VCF files.

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

## Excluding Specific Regions from Analysis

For a variety of reasons it is sometimes desirable to exclude some regions of
the genome from analysis, even when they are included in the exome target regions.
Some common reasons can include:

 * to exclude regions where CNV calling is difficult and thus causes
   large numbers of false positives
 * to mask regions where CNV calls might result in incidental findings (for
   example, that violate ethics constraints for research).

**Note**: regions excluded from analysis are _not_ automatically excluded
from simulation or known CNVs provided as true positives. Thus excluding regions
may result in loss of sensitivity in the output. To exclude regions completely
from use by Ximmer, adjust the `target_regions` parameter.

## Excluding Results Overlapping Specific Genes

Although Ximmer can exclude some regions from analysis, often the reason to do this 
is to avoid including specific genes from the results. Ximmer supports this option 
via the `exclude_genes` setting. Set this option to a text file containing one HGNC 
gene symbol per line to exclude CNVs overlapping the specified genes from your 
results. Note that these CNVs will be excluded even if they overlap genes specified 
by the `gene_filter` option (see below).

Example:

```
exclude_genes="/home/ximmer/genes/excluded_genes.txt"
```


## Filtering Results to Specific Genes

Although Ximmer supports interactively filtering to specific genes in the curation 
interface, it may be desirable to hard filter the result set to a gene list. This 
can ensure only specific genes are looked at, such as when ethics or a clinical indication
for testing limits the scope of the investigation.

To set a gene list, add the `gene_filter` configuration attribute, set to a file containing
a list of HGNC gene symbols (one per line). CNVs will be removed from the results unless 
they overlap at least one gene from the provided set.

Example:

```
gene_filter="/home/ximmer/genes/gene_list.txt
```

## Gene Lists

More advanced filtering and ranking of genes can be set up using gene lists. A gene 
list is a list of gene symbols, with each symbol accompanied by a number indicating 
its priority, which is typically in the range 1 - 4. 

You can assign multiple gene lists in a `genelists` section. Each gene list is identified
by symbol which should be a short sequence of upper case letters, which should be assigned 
to a path to a file that defines the gene symbols and priorities (also called "categories"),
separated by tab characters.

Example configuration:

```
genelists {
   CARDIAC='/home/simon/genelists/cardiac_genes.txt'
}
```

An example gene list would look like:

```
DVL1    3
SCN5A   5
```

**NOTE**: the gene list is applied to all genes overlapped by a CNV, and the whole 
CNV is considered to have the rank of its highest ranked gene.

## Filtering Results Based on Genelists

By default, setting a gene list only causes genes to be highlighted and searchable 
by category in the user interface. However you can also completely exclude genes that 
are not of interest from appearing in your results. To do this, set a minimum category
by adding a `filter` section to your genelists:

```
genelists {
   CARDIAC='/home/simon/genelists/cardiac_genes.txt'
   filter {
     miminum_category=1
   }
}
```


## Minimal Complete Configuration Example

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












