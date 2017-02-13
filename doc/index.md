# Ximmer User Guide

## Introduction

Ximmer is a tool designed to help users of targeted high throughput (or "next generation") 
genomic sequencing data (such as exome data) to accurately detect copy number variants
(CNVs). Ximmer is not a copy number detection tool itself. Rather, it is a framework for
running and evaluating other copy number detection tools. It offers three essential features
that users of CNV detection tools need:

 * A suite of pipelines for running a variety of well known CNV detection tools
 * A simulation tool that can create artificial CNVs in sequencing data for 
   the purpose of evaluating performance
 * A visualisation and curation tool that can combine results from multiple 
   CNV detection tools and allow the user to inspect them, along with 
   relevant annotations.

We created Ximmer because although there are very many CNV detection tools,
they can be hard to run and their performance can be highly variable and
hard to estimate. This is why Ximmer builds in simulation: to allow 
a quick and easy estimation of the performance of any tool on any data set.


## Installation and Requirements

To make Ximmer easier to use we have included support to automatically 
download and build a range of tools. You should make sure before starting
that you have at minimum the following requirements:

 * Java 1.7 or higher
 * Python 2.7, preferably the Anaconda installation
 * R 3.2 or higher

Ideally, these should all be directly accessible from your environment. 
If necessary, you can specify custom locations for them in the configuration file.

You should also make sure you have internet access while doing the installation
because Ximmer will try to download some components. It may be necessary to set 
the "http_proxy" environment variable if your network uses a proxy.


### Run Installer

Ximmer includes a simple installer script to help set up and configure
it for basic operation. To get started:

```
git clone git@github.com:ssadedin/ximmer.git
cd ximmer
./bin/install
```

### Set Configuration Parameters

Ximmer needs some basic settings configured before it can be used. Copy the file 
`eval/pipeline/config.groovy.template` to `eval/pipeline/config.groovy`. Then edit 
the file to set the BASE install location 
as the absolute path to the `eval` directory and the absolute path to the indexed
FASTA file for your reference sequence.


## Configuring Ximmer

Ximmer is designed to run on _targeted_ sequencing data (where that term includes 
the most common type of targeted data, _exomes_). That means to run it,
you need three things:

 * Some exome or targeted capture data consisting of at least 5 
   separate samples.
 * A BED file describing the target regions that were captured
 * An indexed human genome reference in FASTA format

