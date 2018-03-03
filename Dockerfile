# Java 7 is not available on Debian Stretch; use Debian Jessie
FROM debian:jessie
MAINTAINER David J. H. Shih <djh.shih@gmail.com>

ENV DEBIAN_FRONTEND noninteractive

# Install Java Run Time 7, Python 2.7, and 
# other prerequisites for Matlab Common Runtime
RUN apt-get update && apt-get install -y \
    default-jre-headless \
    python \
    wget \
    tar \
    unzip \
    xorg \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*

# Install Matlab Common Runtime 2013a
RUN mkdir -p /tmp/mcr-install && \
    cd /tmp/mcr-install && \
    wget -nv http://ssd.mathworks.com/supportfiles/MCR_Runtime/R2013a/MCR_R2013a_glnxa64_installer.zip && \
    unzip MCR_R2013a_glnxa64_installer.zip && \
    mkdir /opt/mcr && \
    ./install -destinationFolder /opt/mcr -agreeToLicense yes -mode silent && \
    cd / && \
    rm -rf /tmp/mcr-install

# Copy workflow contents
ADD . /opt/ffpe
ENV PATH="${PATH}:/opt/ffpe"
