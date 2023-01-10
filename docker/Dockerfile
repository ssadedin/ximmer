FROM ubuntu:20.04
MAINTAINER Simon Sadedin "simon.sadedin@mcri.edu.au" 
RUN apt update && apt install -y software-properties-common curl wget apt-transport-https
RUN apt update && add-apt-repository universe && add-apt-repository ppa:openjdk-r/ppa && apt update
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9
RUN add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/'
RUN apt update && apt install -y libcurl4-openssl-dev libxml2-dev libmariadbclient-dev r-base-core
RUN apt install -y openjdk-11-jre libfreetype6-dev pkg-config python2 python2-dev python3 python3-dev python3-pip
RUN echo 'local({r <- getOption("repos"); r["CRAN"] <- "http://cran.r-project.org"; options(repos=r); }); install.packages("BiocManager"); BiocManager::install("ExomeDepth");' | Rscript -
RUN echo 'local({r <- getOption("repos"); r["CRAN"] <- "http://cran.r-project.org"; options(repos=r); }); install.packages("BiocManager"); BiocManager::install("CODEX");' | Rscript -

RUN mkdir -p /usr/local/ximmer
ADD bin /usr/local/ximmer/bin 
ADD src /usr/local/ximmer/src
ADD tools /usr/local/ximmer/tools
ADD pipeline /usr/local/ximmer/pipeline
RUN cd /usr/local/ximmer && \
    ./bin/install -q && \
    echo 'JAVA="java"' >> /usr/local/ximmer/eval/pipeline/config.groovy && \
    echo 'java { executable="java" }' >> /usr/local/ximmer/eval/pipeline/bpipe.config && \
    cd /usr/local/ximmer; mkdir cache && cd cache && wget 'http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/dgvMerged.txt.gz' && \
    cd /usr/local/ximmer/cache && wget 'http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz' 
ENV PATH="/usr/local/ximmer/bin:${PATH}"

