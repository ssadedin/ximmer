FROM ubuntu:14.04
MAINTAINER Simon Sadedin "simon.sadedin@mcri.edu.au" 
RUN apt-get update; 
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
RUN echo 'deb http://cran.rstudio.com/bin/linux/ubuntu trusty/' >> /etc/apt/sources.list
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9
RUN apt-get update && \
    apt-get install -y software-properties-common python-software-properties && \
    apt-get install -y curl wget &&  \
    apt-get install -y apt-transport-https && apt-get update
RUN apt-get update && apt-get install -y libcurl4-openssl-dev && apt-get install -y libxml2-dev && apt-get install -y libmariadbclient-dev && apt-get install -y r-base
RUN add-apt-repository ppa:openjdk-r/ppa && apt-get update 
RUN apt-get install -y openjdk-8-jre && apt-get install -y libfreetype6-dev pkg-config python-dev python-pip
RUN mkdir -p /usr/local/ximmer
ADD bin /usr/local/ximmer/bin 
ADD src /usr/local/ximmer/src
ADD eval /usr/local/ximmer/eval
RUN cd /usr/local/ximmer && \
    ./bin/install -q && \
    echo 'JAVA="java"' >> /usr/local/ximmer/eval/pipeline/config.groovy && \
    echo 'java { executable="java" }' >> /usr/local/ximmer/eval/pipeline/bpipe.config && \
    cd /usr/local/ximmer; mkdir cache && cd cache && wget 'http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/dgvMerged.txt.gz' && \
    cd /usr/local/ximmer/cache && wget 'http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz' 
ENV PATH="/usr/local/ximmer/bin:${PATH}"
ENV JAVA_HOME="/usr/lib/jvm/java-8-openjdk-amd64/jre"

