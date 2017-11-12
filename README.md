Ximmer
======


Ximmer is a tool designed to help users of exome and targeted genomic
sequencing data accurately detect and interpret copy number variants (CNVs).
Ximmer is not a copy number detection tool itself. Rather, it is a framework
for running other copy number detection tools and interpreting their results.
It offers three essential features that users of CNV detection tools need:

 * A suite of pipelines for running a variety of well known CNV detection tools
 * A simulation tool that can create artificial CNVs in sequencing data for 
   the purpose of evaluating performance
 * A visualisation and curation tool that can combine results from multiple 
   CNV detection tools and allow the user to inspect them, along with 
   relevant annotations.

All of these are integrated into one streamlined package that you can run
easily on any data set you want to analyse.

We created Ximmer because although there are very many CNV detection tools,
they can be hard to run and their performance can be highly variable and
hard to estimate. This is why Ximmer builds in simulation: to allow 
a quick and easy estimation of the performance of any tool on any data set.

See an online [example report](http://example.ximmer.org) to get an idea what 
Ximmer's output looks like.


Requirements
============

 * Java 1.7 (note: Java 1.8 does not work, unless you upgrade the bundled GATK)
 * Python 2.7, preferably the Anaconda installation
 * R 3.2 or higher
 * 24GB of RAM

or 

 * Use [Docker](https://ssadedin.github.io/ximmer/docker.html)!


Building and Running It
=======================

See the online [documentation](https://ssadedin.github.io/ximmer/) for more details!



