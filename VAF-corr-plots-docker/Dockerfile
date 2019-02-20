FROM ubuntu:18.04
FROM python:3.6
FROM r-base

CMD source ~/.bashrc

MAINTAINER Khushbu Patel <patelk26@email.chop.edu>

# installing Git
RUN apt-get update && apt-get install -y git
 
# clone the most recent version of the pipeline
WORKDIR /home/user/
RUN git clone https://github.com/marislab/create-pptc-pdx-corplots.git


# Install all the software needed to run the pipeline
RUN apt-get install -qqy python3-setuptools python3-docutils python3-flask python3.6 python-pip
RUN pip install --upgrade pip && \
	pip install --no-cache-dir pandas


# Setting working directory
WORKDIR /home/user/VAF-corr-docker


# Executing scripts
RUN [ "/usr/bin/bash", "/home/user/VAF-corr-docker/pipeline.sh" ]

