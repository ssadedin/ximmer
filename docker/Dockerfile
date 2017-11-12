FROM ubuntu:14.04
MAINTAINER Simon Sadedin "simon.sadedin@mcri.edu.au" 
RUN apt-get update; 
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
RUN apt-get update && apt-get install -y software-properties-common python-software-properties && apt-get install -y apt-transport-https && apt-get update
RUN curl -s "https://get.sdkman.io" | bash && . "/root/.sdkman/bin/sdkman-init.sh" && sed -i -e 's/sdkman_auto_answer=false/sdkman_auto_answer=true/g' /root/.sdkman/etc/config
RUN apt-get install -y libcurl4-openssl-dev && apt-get install -y libxml2-dev && apt-get install -y r-base
RUN apt-get install -y git
RUN sdk install java 7u141-zulu && source ~/.bashrc
RUN pushd /usr/local && \
    git clone https://github.com/ssadedin/ximmer.git && cd ximmer && ./bin/install -q && \
    echo 'JAVA="/root/.sdkman/candidates/java/current/bin/java"' >> /usr/local/ximmer/eval/pipeline/config.groovy && \
    echo 'java { executable="/root/.sdkman/candidates/java/current/bin/java" }' >> /usr/local/ximmer/eval/pipeline/bpipe.config && \
    pushd /usr/local/ximmer; mkdir cache && cd cache && wget 'http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/dgvMerged.txt.gz' && popd && \
    pushd /usr/local/ximmer/cache && wget 'http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz' && popd && \
    popd
ENV PATH="/usr/local/ximmer/bin:${PATH}"
ENV JAVA_HOME=/root/.sdkman/candidates/java/current/bin/java

