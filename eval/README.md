CNV Evaluation Framework README
===============================

This document describes how to set up and use the CNV evaluation framework
included with Ximmer.

The framework is based on [Bpipe](http://bpipe.org), a data analysis pipeline
tool that is designed to make running complex analysis pipelines easier. 

Install / Setup
----------------

**1. Install Bpipe**

In order to run the evaluation pipeline, you first need to install Bpipe. Please
use Bpipe 0.9.8.6 or later.

**2. Install the tools**

You need to install the following:

  * EXCAVATOR - see instructions [here](tools/excavator/README.md)
  * XHMM - see instructions [here](tools/xhmm/README.md)
  * ExomeDepth - install into your R environment using:

    install.packages(ExomeDepth)

  * cn.MOPS - install into your R environment using:
 
    source("http://bioconductor.org/biocLite.R")
    biocLite("cn.mops")

The pipeline also needs the VariantAnnotation and rootSolve package in your R environment:

    source("http://bioconductor.org/biocLite.R")
    biocLite("VariantAnnotation")
    install.packages("rootSolve")

**3. Provide a BED file of capture region**

CNV callers need to know which parts of the genome to analyze. The providers
of exome capture kits provide a BED file that specifies the regions covered by 
probes which are suitable for this purpose. Note: if the BED file contains overlapping
regions then you should first flatten them.

Create a directory for your BED file under "designs" and copy it to that directory.
You will need to specify this BED file as the design in step 4.

**4. Provide a human genome reference in FASTA format**

This should be an HG19 reference. The FASTA should be indexed with samtools
faidx.

**5. Copy config.groovy.template to config.groovy and edit**

You should read through the entries and edit as appropriate. You will need to 
specify things such as your target region BED file (step 3), and your human genome
reference (step 4).

Running
----------------

The inputs to the pipeline are indexed BAM files. Optionally, you can also
provide a VCF file for each BAM file.

To run the CNV pipeline, execute it using Bpipe as follows:

bpipe run ./pipeline/exome_cnv_pipeline.groovy <bam1>.bam <bam2>.bam ....

