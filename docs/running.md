# Running Ximmer 

## Introduction

Ximmer has two separate parts:

 * simulation of CNVs 
 * analysis of data containing CNVs

You can choose to do either or both of simulation and analysis. If you just want some 
CNV results from one or more CNV callers, then just do analysis. If you just 
want to simulate some CNVs and then take the data to analyse separately, 
then just do the simulation part. If you want everything, then do both.

## The Configuration File

Whichever way you are using Ximmer, the starting point is to create 
a configuration file that controls the analysis. This file can have a lot of settings,
but if you are only doing analysis, the most minimal form can simply contain:

 * bam_files setting describing where the BAM files to analyse or simulate from are
 * target_regions setting describing the capture region
 * callers section describing which CNV callers to run (see [Analysis](analyses.md)).

If you are doing simulation, then you need to also specify:

 * simulation_type (one of 'replacement' or 'downsample')
 * regions - a range indicating the number of target regions that should be included in the simulated CNVs

For details on how to set these parameters, see the relevant sections:

 * for analysis, see [Analysis Configuration](analyses.md)
 * for simulation, [Simulation Configuration](simulations.md)

The configuration file can be created in any directory. Typically, you should create 
a directory to run Ximmer in for a data set, and put the configuration file there. 
Ximmer will create a large number of files in this directory, so make sure it has 
plenty of space.

## Running Ximmer

Once you have created a configuration file describing your run (we'll call it `config.groovy`), you 
can start Ximmer. Running Ximmer is similar for both Simulation and Analysis.

### Simulation and Analysis Mode

To run in both modes, use a command such as the following:

```
<ximmer install dir>/bin/ximmer -v -c config.groovy -o results
```

This will run the full Ximmer process, placing the results into the `results` directory. You 
will find the output report in the `results` directory, as an HTML file. The name of 
the file may vary dependending on your configuration, but by default you will find 
it as "analysis.html". This report shows the full details of the CNVs simulated, as well 
as plots indicating the sensitivity and specificity of the different callers tested.

Note that in this mode, all calls that are not either simulated or provided as a 
set of pre-defined true positives are considered as false positives.

### Analysis Mode

To run in Analysis Mode, just add `-nosim` as an argument:

```
<ximmer install dir>/bin/ximmer -v -nosim -c config.groovy -o results
```

After the analysis runs, a report will be produced in the following location:

```
results/run1/analysis/report/local_cnv_report.html
```

This report contains all the CNV calls by every caller used as well as informative 
plots in an interface that supports filtering and annotation of the data.









