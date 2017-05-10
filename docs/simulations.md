# Ximmer Simulations

## Introduction

This section explains how to configure Ximmer to simulate CNVs. Ximmer simulates
CNVs using two alternative mechanisms. 

 * by taking input files (BAM/CRAM format) and downsampling them over
   sequential target regions to create artificial "single copy deletions".
 * by taking female input files (again, CRAM/BAM format) and replacing 
   all the reads aligned to sections of the X chromosome with reads from
   the same segment in a male sample. For this to be valid, normalisation
   needs to be done to make the read counts comparable before swapping the
   reads between samples.

## Specifying Simulation Type

As mentioned above, Ximmer can do two kinds of simulation. You should specify 
in your configuration file which kind you want of the following:

 * replace
 * downsample
 * none

Example:

```
simulation_type='downsample'
```

*Note*: The `none` simulation type is intended for cases where you wish to reanalyse
pre-existing simulated data, or data that has a known set of identified CNVs in it.
For that case, you should be specifying the known CNVs via the `known_cnvs` attribute 
(see below). When you specify `none`, Ximmer will not create new BAM files and instead
will run the analysis directly on the source BAM files.

## Specifying CNV Sizes

One of the most important aspects of CNV detection is what size of event you
are looking to find. If you are after single exon deletions, for example, you
may end up needing quite different settings to those you would have for detecting 
megabase size events. Accordingly, you can set how big you want the CNVs to be
that Ximmer creates. In Ximmer this is controlled by setting the number of target
regions to include in each event. (Note that Ximmer will never begin or end a
simulated CNV inside a target region). The number of regions is specified using a 
range with `..` syntax. For example, to specify that CNVs should be between 5 and 20 
target regions in size, the following setting would be used:

```
regions=5..20
```

## Specifying how many CNVs to Simulate in each Sample

By default Ximmer will only simulate 1 CNV per sample. This minimises the 
potential interference in accurate detection through contamination of the 
normalisation methods by additional CNVs in each sample. However it is less
efficient and means that you need more simulated BAM files than you will if
you simulate more CNVs in each BAM. You can change the number using:

```
deletionsPerSample=10
```

Ximmer will randomly select the number of target regions to include in each 
CNV from this range. Note that these will not always be exactly honored because
sometimes the ranges chosen are expanded if there is continuous read coverage
between the target regions. Since they are sometimes made larger, but never 
smaller, you will generally see that the result is slightly inflated compared to 
the range you select here.

## Specifying Sample Sex

For X-replacement to work, Ximmer needs to known which samples are female 
and which are male. You don't *have* to specify sample sex: if you don't,
Ximmer will guess it by counting the number of reads on the X, Y and non-sex
chromosomes. However this takes some non-trivial computing resources and can 
sometimes be inaccurate, so Ximmer allows you to specify this. There are two ways
to specify it: you can supply a PED file:

```
ped_file="/path/to/your/ped_file.ped"
```

Or in a `samples` section in the configuration:

```groovy
samples {
    males = [
        "SAMPLE_1",
        "SAMPLE_2",
        ...
    ]

    females = [
        "SAMPLE_3",
        "SAMPLE_4",
        ...
    ]
}
```

Note that these sample ids *must* overlap those found in the BAM files that you specify 
as input files. 

*Note*: you can also use this method to make your simulation or analysis operate on
a subset of the samples provided in the BAM files you specify. This can be useful
if you have a single large directory of BAM files but you want to analyse or simulate
from only a portion of them.


## Specifying Number of Runs

Getting access to enough samples to generate good power for determining accuracy
can sometimes be difficult if you do not sequence a large number of samples. To
help with this, Ximmer allows you to increase the number of CNVs simulated in two 
dimensions:

 * Increase the number of CNVs simulated per sample - this can help but is limited
   especially for small targeted panels as you will not want to simulate CNVs too 
   close together, nor to simulate too many in any individual sample as this will
   bias read counts away from normal overall and thus compromise the simulation 
   accuracy

To set the number of CNVs per sample to simulate use the `deletionsPerSample` option:

```
deletionsPerSample=10
```

 * Reuse each sample multiple times. This concept is referred to as having 
   "multiple runs", because you are running the simulation itself multiple times and
   aggregating the results.

To specify multiple runs, set the `runs` parameter in the configuration file to an
integer indicating how many times to use each sample. Eg:

```groovy
runs=5
```

## Specifying Known CNVs

By default Ximmer assumes that every CNV that is detected but not a simulated
CNV is a false positive. This, of course, is not true because samples can 
have real CNVs in them. It is also important to avoid simulating a CNV on top
of another CNV. For this reason, Ximmer allows you to specify known CNVs
in a BED-like tab-separated format consisting of columns:

 * chromosome
 * start position
 * end position
 * sample id

*Note*: in the current code, the type of CNV (duplication vs deletion) is not specified; 
an event will be counted as a true positive if it is detected overlapping the known CNV 
even if the type of event detected is wrong.


## Advanced: Using Per-Run Settings

A more advanced configuration option exists that lets you specify multiple runs
using different source files for each run. In this case, set the `runs` attribute
to a sub-configuration containing specific `known_cnvs` and `bam_files` entries
for each run. Each sub-configuration can be tagged with a specific name to make 
it recognisable in the results:

```groovy
runs {
    "sim0_1" { 
        known_cnvs="/home/simon/bam_files/sim0.1/true_cnvs.bed" 
        bam_files="/home/simon/bam_files/sim0.1/*.bam" 
    }
    "sim0_2" { 
        known_cnvs="/home/simon/bam_files/sim0.2/true_cnvs.bed" 
        bam_files="/home/simon/bam_files/sim0.2/*.bam" 
    }
    "sim3_1" { 
        known_cnvs="/home/simon/bam_files/sim3.1/true_cnvs.bed" 
        bam_files="/home/simon/bam_files/sim3.1/*.bam" 
    }
}
```

## Specifying Analysis Settings

After Ximmer simulates CNVs in the sample data, it then runs the configured CNV callers
on the data to test the sensitivity and specificity of each caller. The settings for
how the CNV callers are run are specificed by two blocks:

 * `callers` which defines which callers are run and the default settings to use
 * `analyses` (optional) which defines groups of settings to use when running each
    CNV caller.
    
See the [Analysis](analysis.md) section for details about how to configure these
sections in detail.

## Using CRAM Format

CRAM format saves a lot of space but requires that a reference sequence be specified. 

This needs to be configured in two places for CRAM format to work:

 * Set the correct FASTA file for the reference in eval/pipeline/config.groovy
 * Set the environment variable XIMMER_REF to the absolute path of the FASTA reference 
   file, eg:

```
export XIMMER_REF=/path/to/your/reference.fasta
```

## Full Configuration Example

```groovy
/**
 * Ximmer Example Configuration
 *
 * This example simulates CNVs in exome data from the Simons Simplex 
 * collection. It then runs 4 CNV callers on the output data in 4 different
 * configurations and creates a report comparing the sensitivity and 
 * specificity of each caller and configuration. The simulation incorporates
 * a set of known CNVs published along with the data into the process.
 */

title="Test Simulation"

bam_files="/home/simon/example/*.recal.bam"

target_regions="/home/simon/example/NIMBLEGENV2.bed"

simulation_type="downsample"

concurrency=80

known_cnvs="/home/simon.sadedin/work/ximmer/eval/k_e_cnvs.tsv"

// Number of separate runs to complete
runs=4

regions=2..10

deletionsPerSample=10

samples {
    females = [
        "SRR1301839",
        "SRR1301908",
        "SRR1301932",
        "SRR1301896",
        "SRR1301260",
        "SRR1301870",
        "SRR1301936",
        "SRR1301912",
        "SRR1301312",
        "SRR1301605",
        "SRR1301928",
        "SRR1301609",
        "SRR1301924",
        "SRR1301554",
        "SRR1301855",
        "SRR1301785",
        "SRR1301893",
        "SRR1301402",
        "SRR1301456"
    ]
}

dgv {
    max_freq = 0.01
    min_study_size = 5
}

callers {
    xhmm {
        exome_wide_cnv_rate=1e-08
        mean_number_of_targets_in_cnv=3
    }
    exomedepth { transition_probability=0.0001 }
    cnmops { prior_impact=10; min_width=2; lower_threshold=-0.8 }
    conifer { conifer_svd_num=1 }
}

analyses {
    '1' {
        xhmm { exome_wide_cnv_rate=1e-04 }
        exome_depth { transition_probability=0.01 }
        conifer { conifer_svd_num=2 }
        cnmops { prior_impact=5 }
    }
    '2' {
        xhmm { exome_wide_cnv_rate=1e-02 }
        exome_depth { transition_probability=0.000001 }
        conifer { conifer_svd_num=4 }
        cnmops { prior_impact=20 }
    }
    '3' {
        xhmm { exome_wide_cnv_rate=1e-04; xhmm_pve_mean_factor=0.9 }
        cnmops { prior_impact=5;  lower_threshold=-0.6 }
    }
    '4' {
        xhmm { exome_wide_cnv_rate=1e-04; xhmm_pve_mean_factor=0.5 }
        cnmops { prior_impact=5;  lower_threshold=-0.4 }
    }
}
```
