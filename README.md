Ximmer
======

Ximmer is a tool for simulation deletions in targeted sequencing data such as exome data.

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

To see the help:

    java -jar build/libs/ximmer.jar

A simple example would be:

    java -Xmx1g -jar ximmer.jar -n 1 -f NA12878 -m NA19239 -bam na12878.bam -bam na19239.bam -r target.bed -o deletions.bed 

This will create a single output BAM file that isn named starting with
the female sample name (NA12878) suffixed with the region of the deletion that was 
simulated. An output BED file specifying the deletions would be created as 'deletions.bed'.

Note that you could specify multiple males, in which case the male to use would be selected randomly. 

Also note that for a reproducible result, you may like to specify the random seed using the 
-seed option.

Limitations
===========

Ximmer currently assumes that there is only a single sample in each BAM file. 
