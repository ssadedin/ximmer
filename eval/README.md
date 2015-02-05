CNV Evaluation Framework README
===============================

This document describes how to set up and use the CNV evaluation framework
included with Ximmer.

The framework is based on [Bpipe](http://bpipe.org), a data analysis pipeline
tool that is designed to make running complex analysis pipelines easier. 

Install / Setup
----------------

**1. Install the tools**

You need to install the following:

  * EXCAVATOR - see instructions [here](tools/excavator/README.md)
  * XHMM - see instructions [here](tools/xhmm/README.md)
  * ExomeDepth - install into your R environment using:

    install.packages("ExomeDepth")

  * cn.MOPS - install into your R environment using:
 
    source("http://bioconductor.org/biocLite.R");
    biocLite("cn.mops")

  * samtools - should be in your default path, but you can edit config.groovy to change that
  * R - again, should be in your path

The pipeline also needs some packages in your R environment:

    source("http://bioconductor.org/biocLite.R")
    biocLite("VariantAnnotation")
    biocLite("Rsamtools")
    biocLite("GenomicRanges")
    install.packages("rootSolve")

**2. Provide a BED file of capture region**

CNV callers need to know which parts of the genome to analyze. The providers
of exome capture kits provide a BED file that specifies the regions covered by 
probes which are suitable for this purpose. Note: if the BED file contains overlapping
regions then you should first flatten them.

Create a directory for your BED file under "designs" and copy it to that directory.
You will need to specify this BED file as the design in step 4.

**3. Provide a human genome reference in FASTA format**

This should be an HG19 reference. The FASTA should be indexed with samtools
faidx.

**4. Copy pipeline/config.groovy.template to pipeline/config.groovy and edit**

You should read through the entries and edit as appropriate. You will need to 
specify things such as your target region BED file (step 3), and your human genome
reference (step 4).

Running
----------------

The inputs to the pipeline are indexed BAM files. Optionally, you can also
provide a VCF file for each BAM file. Providing a VCF will add annotation of 
variants in the region of each CNV that is called in the final output.

To run the CNV pipeline on data simulated by Ximmer, you should run it like this:

    ./bpipe run -p true_cnvs=<true cnvs bed file> \
                -p batch_name=<name for your run> \
                -p target_bed=<target region BED file>\
                ./pipeline/exome_cnv_pipeline.groovy <bam1>.bam <bam2>.bam ... <vcf1>.vcf <vcf2.vcf> ....

In this case it will run only on the X chromosome and it will compare the 
CNVs called to the "truth" set provided in true_cnvs.bed.

If you just want an analysis of non-simulated data, you can run it non-simulation mode:

    ./bpipe run -p callers=ed,xhmm,mops \
                -p simulation=false \
                -p batch_name=<name for your run> \
                -p target_bed=<target region BED file>\
                ./pipeline/exome_cnv_pipeline.groovy <bam1>.bam <bam2>.bam ... <vcf1>.vcf <vcf2.vcf> ....

Here we have specified to run only a subset of the CNV callers (ExomeDepth, XHMM and cn.MOPS). You can 
include EXCAVATOR by adding "ex". We find EXCAVATOR can take a long time to run, and does not
add a lot of calls that you don't get from the others, so if you want a quick analysis you can
leave it out.

**NOTE:** In both cases above the VCFs are optional. If you provie them then variant allele balance
will be annotated over CNVs that are called so that you can review the evidence from LOH or distortion
of balance for support of the CNVs called.

What are the results?
---------------------

The results are zipped up into a file called "cnv_report.zip" that will appear in the runs/&lt;batch_name&gt; directory.

There is an HTML file called cnv_report.html that you can open to look at CNVs interactively, and there is a TSV
file that you can open ina program like R or Excel to look at results programmatically. 




