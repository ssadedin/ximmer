# Ximmer User Guide

## Introduction

Docker makes it especially easy to build and install all the different tools 
and CNV callers that Ximmer uses. Running using Docker requires some small
extra accommodations when invoking Ximmer.  This section describes these.

## Building

Build the docker image using the following docker invocation:

```
git clone git@github.com:ssadedin/ximmer.git
cd ximmer/docker
docker build -t ximmer . 
```

## Running

After building, Ximmer is runnable straight from docker:

```bash
$ docker run ximmer ximmer
**************************************************

#     #  ###  #     #  #     #  #######  ######   
 #   #    #   ##   ##  ##   ##  #        #     #  
  # #     #   # # # #  # # # #  #        #     #  
   #      #   #  #  #  #  #  #  #####    ######   
  # #     #   #     #  #     #  #        #   #    
 #   #    #   #     #  #     #  #        #    #   
#     #  ###  #     #  #     #  #######  #     #  

**************************************************
error: Missing required options: c, o
usage: ximmer -c <config> -o <output_directory>
 -c <arg>      Configuration file
 -nosim        Analyse samples directly (no simulation)
 -o <arg>      Output directory
 -seed <arg>   Random seed
 -simonly      Run only simulation component
 -v            Verbose output
```

However there is a problem: Ximmer cannot see files outside Docker. So you would
need to map the file systems where the files reside into the docker image. Specifically:

 * you have to map the files on your file system to where the Docker
   container can see them. These include:
   - Human genome reference (FASTA)
   - Ximmer analysis configuration
   - Source data files: BAM files, target region BED files, etc.
 * you need to pass through the reference using an environment variable when
   you run Ximmer.

Here is an example of an invocation of Ximmer which does these things:

```bash
docker run
     -e XIMMER_REF=/reference/ucsc.hg19.fasta \  # tell Ximmer the reference FASTA 
     -v /Volumes/reference:/reference \ # map the file system containing the reference
     -v /Volumes/data:/data \  # map the file system containing the data
        ximmer ximmer \ # run the ximmer executable in the ximmer image
     -c /data/config.groovy \
     -v \
     -o /data/output.docker # send output to the /data
```

Also important is that the `config.groovy` file you pass will contain file paths itself.
These must be relative to the locations *inside Docker*, not locations on the computer
outside.

You should also ensure that machine Ximmer is running on has a reasonable amount of 
memory. Ximmer will require at least 24g of memory to run, unless you adjust the internal
settings.
