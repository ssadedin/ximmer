Ximmer
======

Ximmer is a tool for simulation of single copy deletions in targeted sequencing data 
such as exome data.

Ximmer can simulate deletions by two methods:

  * swapping reads into a female from a male on the X chromosome
  * downsamping the number of reads by 50% over a select number of targets

Requirements
============

All you need is Java 1.6+

Building It
===========

Try this (on Unix-like environments):

    git clone https://github.com/ssadedin/ximmer.git
    cd ximmer
    ./gradlew jar

The result should be build/libs/ximmer.jar. This is a runnable JAR file you can use to launch Ximmer.

Running It
===========

Ximmer is fairly simple to run. Since one of the simulation modes depends on the sex of the samples, 
you currently need to specify the female sample names using the -f option (repeatedly, if you want) and 
the male samples using -m. Deletions will currently only be simulatd in female samples.

To see the help:

    java -jar build/libs/ximmer.jar

A simple example would be:

    java -Xmx1g -jar ximmer.jar -n 1 -f NA12878 -m NA19239 -bam na12878.bam -bam na19239.bam -r target.bed -o deletions.bed 

This will create a single output BAM file that isn named starting with
the female sample name (NA12878) suffixed with the region of the deletion that was 
simulated. An output BED file specifying the deletions would be created as 'deletions.bed'.

Note that you could specify multiple males, in which case the male to use would be selected 
randomly. If you specify multiple females, each female will have a single deletion simulated.

Also note that for a reproducible result, you may like to specify the random seed using the 
-seed option.

Limitations
===========

Ximmer currently assumes that there is only a single sample in each BAM file. 
