Ximmer
======


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


Requirements
============

 * Java 1.7 (note: Java 1.8 does not work, unless you upgrade the bundled GATK)
 * Python 2.7, preferably the Anaconda installation
 * R 3.2 or higher


Building and Running It
=======================

See [documentation](doc/index.md) for more details!


