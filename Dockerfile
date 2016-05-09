FROM ubuntu:latest

RUN apt-get update
RUN apt-get install -y libcurl4-openssl-dev libxml2-dev apt-transport-https
RUN echo 'deb https://cran.cnr.berkeley.edu/bin/linux/ubuntu trusty/' >> /etc/apt/sources.list
RUN apt-get update -y --force-yes

RUN apt-get install -y --force-yes r-base-dev

RUN Rscript -e 'source("http://bioconductor.org/biocLite.R"); biocLite("limma",ask=FALSE);'

COPY Dockerfile /opt/

COPY DE_comparison.R /opt/

CMD ["/bin/bash"]

MAINTAINER Jing Zhao jing.zhao@sbgenomics.com