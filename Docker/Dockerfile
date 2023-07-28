FROM ubuntu:20.04


RUN apt-get update && apt-get install -y software-properties-common && \
    apt-get update && apt-get install -y \
        build-essential \
        cmake \
        git \
        python3 \
        python3-pip \
        libboost-filesystem-dev \
        libboost-regex-dev \
        libboost-system-dev \
        libbz2-dev \
        libcurl3-dev \
        liblzma-dev \
        libpthread-stubs0-dev \
        wget \
	default-jre \
        zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*








WORKDIR /usr/local/src
RUN git clone --recursive https://github.com/MethodsDev/rnaseqc_plus.git && \
      cd rnaseqc_plus/rnaseqc && \
      make && \
      ln -sf /usr/local/src/rnaseqc_plus/rnaseqc/rnaseqc /usr/local/bin/rnaseqc && \
      cp ../plugins/trimmomatic.jar /usr/local/bin && \
      cp ../fastq_read_trimmer/fastq_polyAT_trimmer.py /usr/local/bin


ENV SRC /usr/local/src
ENV BIN /usr/local/bin



########
# Samtools

ENV SAMTOOLS_VERSION=1.9

RUN SAMTOOLS_URL="https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2" && \
   cd $SRC && \
   wget $SAMTOOLS_URL && \
   tar xvf samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
   cd samtools-${SAMTOOLS_VERSION}/htslib-${SAMTOOLS_VERSION} && ./configure && make && make install && \
   cd ../ && ./configure --without-curses && make && make install

##############
## STAR

ENV STAR_VERSION=2.7.8a
RUN STAR_URL="https://github.com/alexdobin/STAR/archive/${STAR_VERSION}.tar.gz" &&\
    wget -P $SRC $STAR_URL &&\
    tar -xvf $SRC/${STAR_VERSION}.tar.gz -C $SRC && \
    mv $SRC/STAR-${STAR_VERSION}/bin/Linux_x86_64_static/STAR /usr/local/bin

