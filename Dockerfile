FROM ubuntu:20.04
LABEL org.opencontainers.image.authors="xhu556@gmail.com"

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get -qq update && apt-get -qq -y install \
    automake \
    build-essential \
    bzip2 \
    curl \
    wget \
    unzip \
    make \
    git \
    perl \
    default-jre \
    libssl-dev \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libz-dev \
    libfreetype6-dev \
    libpng-dev \
    libcurl4-gnutls-dev \
    libdeflate-dev \
    libncurses-dev \
    autogen \
    libtool \
    shtool \
    libblas-dev \
    pkg-config \
    python3 \
    python3-dev \
    python3-pip \
    python3-matplotlib \
    python3-distutils \
    rsync \
    texlive-latex-base \
    tzdata \
    nano \
    apt-utils

## set up tool config and deployment area:
ENV SRC=/opt
ENV BIN=/usr/local/bin
ENV LD_LIBRARY_PATH=/usr/local/lib

#rename python3
RUN ln -sf /bin/python3 /bin/python

#install python packages
#install biopython
RUN pip install biopython

#install Parallel
RUN apt-get -qq -y install parallel

#install multiqc
RUN pip install multiqc

#install pandas
RUN pip install pandas

#apt-get install fastqc
WORKDIR $SRC
RUN wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip && \
    unzip fastqc_v0.12.1.zip && \
    chmod 755 $SRC/FastQC/fastqc && \
    ln -s $SRC/FastQC/fastqc $BIN/.

#install Trimmomatic
WORKDIR $SRC
RUN wget https://github.com/usadellab/Trimmomatic/files/5854859/Trimmomatic-0.39.zip && \
    unzip Trimmomatic-0.39.zip && \
    rm -rf Trimmomatic-0.39.zip

#install quast
WORKDIR $SRC
RUN wget https://github.com/ablab/quast/releases/download/quast_5.2.0/quast-5.2.0.tar.gz && \
    tar -xvzf quast-5.2.0.tar.gz && \
    cd quast-5.2.0 && \
    python setup.py install

#install Cmake
WORKDIR $SRC
RUN wget https://github.com/Kitware/CMake/releases/download/v3.25.2/cmake-3.25.2-linux-x86_64.tar.gz && \
    tar -xvzf cmake-3.25.2-linux-x86_64.tar.gz && \
    cd cmake-3.25.2-linux-x86_64 && \
    ln -s $SRC/cmake-3.25.2-linux-x86_64/bin/* $BIN/.

#install autoconf-2.69, this version is required by Trinity2.15
WORKDIR $SRC
RUN wget ftp://ftp.gnu.org/gnu/autoconf/autoconf-2.69.tar.gz && \
    tar zxvf autoconf-2.69.tar.gz && \
    cd autoconf-2.69 && \
    ./configure --prefix=/usr/ && \
    make && make install


#install jellyfish
RUN wget https://github.com/gmarcais/Jellyfish/releases/download/v2.3.1/jellyfish-2.3.1.tar.gz && \
    tar -xvzf jellyfish-2.3.1.tar.gz && \
    cd jellyfish-2.3.1 && \
    ./configure && \
    make  && \
    make install

#install salmon
RUN wget https://github.com/COMBINE-lab/salmon/releases/download/v1.9.0/salmon-1.9.0_linux_x86_64.tar.gz && \
    tar -xvzf salmon-1.9.0_linux_x86_64.tar.gz 
RUN ln -s $SRC/salmon-1.9.0_linux_x86_64/bin/salmon $BIN/salmon    
    

#Install MUMmer4
WORKDIR $SRC
RUN wget https://github.com/mummer4/mummer/releases/download/v4.0.0rc1/mummer-4.0.0rc1.tar.gz && \
    tar -xvzf mummer-4.0.0rc1.tar.gz && \
    cd mummer-4.0.0rc1 && \
    ./configure && \
    make && make install && \
    ldconfig

#install BBtools
#RUN apt-get -qq -y install bbmap  #don't install Ubuntu BBmap
WORKDIR $SRC
RUN wget https://downloads.sourceforge.net/project/bbmap/BBMap_39.19.tar.gz && \
    tar -xvzf BBMap_39.19.tar.gz
RUN ln -s $SRC/bbmap/*.sh $BIN/.

WORKDIR $SRC
#install Bowtie2
#apt-get -qq -y install bowtie2
RUN wget https://github.com/BenLangmead/bowtie2/releases/download/v2.5.4/bowtie2-2.5.4-linux-x86_64.zip && \
    unzip bowtie2-2.5.4-linux-x86_64.zip  && \
    mv bowtie2-2.5.4-linux-x86_64 bowtie2-2.5.4  && \
    ln -s $SRC/bowtie2-2.5.4/bowtie2* $BIN/.

WORKDIR $SRC
#install bedtools
#RUN apt-get -qq -y install bedtools
RUN mkdir bedtools2 && \
    cd bedtools2  && \
    wget https://github.com/arq5x/bedtools2/releases/download/v2.31.0/bedtools.static  && \
    mv bedtools.static bedtools  && \
    chmod a+x bedtools  && \
    cd .. && \
    ln -s $SRC/bedtools2/bedtools $BIN/bedtools

#install seqtk
WORKDIR $SRC
RUN git clone https://github.com/lh3/seqtk.git  && \
    cd seqtk  && \
    make  && \
    ln -s $SRC/seqtk/seqtk $BIN/seqtk 

#Install BWA
WORKDIR $SRC
RUN git clone https://github.com/lh3/bwa.git && \
    cd bwa; make  && \
    cd ..  && \
    ln -s $SRC/bwa/bwa $BIN/bwa

#Install Krona
WORKDIR $SRC
RUN wget https://github.com/marbl/Krona/releases/download/v2.8.1/KronaTools-2.8.1.tar && \
    tar -xvf KronaTools-2.8.1.tar && \
    cd KronaTools-2.8.1 && \
    ./install.pl
RUN ln -s $SRC/KronaTools-2.8.1/updateAccessions.sh $BIN/updateAccessions.sh
RUN ln -s $SRC/KronaTools-2.8.1/updateTaxonomy.sh $BIN/updateTaxonomy.sh
RUN rm -rf $SRC/KronaTools-2.8.1.tar

#instal NCBI EDirect
WORKDIR $SRC
RUN wget ftp://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/edirect.tar.gz  && \
    gunzip -c edirect.tar.gz | tar xf -  && \
    cd edirect && \
    wget ftp://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/xtract.Linux.gz &&\
    gunzip -f xtract.Linux.gz && \
    chmod +x xtract.Linux && \
    cd .. && \
    find $SRC/edirect -type f -executable -exec ln -s {} $BIN/ \;


#install htslib
WORKDIR $SRC
RUN wget https://github.com/samtools/htslib/releases/download/1.21/htslib-1.21.tar.bz2 && \
    tar -xvjf htslib-1.21.tar.bz2 && \
    cd htslib-1.21 && \
    autoheader  && \
    autoconf -Wno-syntax  && \
    ./configure && \
    make && make install

#install samtools
WORKDIR $SRC
RUN wget https://github.com/samtools/samtools/releases/download/1.21/samtools-1.21.tar.bz2 && \
    tar -xvjf samtools-1.21.tar.bz2 && \
    cd samtools-1.21 && \
    autoheader  && \
    autoconf -Wno-syntax  && \
    ./configure --with-htslib=/usr/local && \
    make  && \
    make install

WORKDIR $SRC
#RUN apt-get -qq -y install bcftools
#install bcftools
RUN wget https://github.com/samtools/bcftools/releases/download/1.21/bcftools-1.21.tar.bz2 && \
    tar -xvjf bcftools-1.21.tar.bz2 && \
    cd bcftools-1.21 && \
    autoheader  && \
    autoconf -Wno-syntax  && \
    ./configure --with-htslib=/usr/local && \
    make  && \
    make install

#Installing NCBI-BLAST+
WORKDIR $SRC
RUN wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.16.0+-x64-linux.tar.gz && \
    tar -xvzf ncbi-blast-2.16.0+-x64-linux.tar.gz 
RUN ln -s $SRC/ncbi-blast-2.16.0+/bin/* $BIN/.

#install diamond
WORKDIR $SRC
RUN wget https://github.com/bbuchfink/diamond/releases/download/v2.1.11/diamond-linux64.tar.gz && \
    tar -xvzf diamond-linux64.tar.gz && \
    mv diamond $BIN/diamond

#Install Kaiju
WORKDIR $SRC
RUN wget https://github.com/bioinformatics-centre/kaiju/archive/refs/tags/v1.10.1.tar.gz && \
    tar -xvzf v1.10.1.tar.gz && \
    cd kaiju-1.10.1/src && \
    make && \
    ln -s $SRC/kaiju-1.9.2/bin/* $BIN/  && \
    cd ..

WORKDIR $SRC
#install recent version Kraken2
#apt-get -qq -y install kraken2 #download library fail error, use git version
RUN wget https://github.com/DerrickWood/kraken2/archive/refs/tags/v2.14.tar.gz && \
    tar -xvzf v2.14.tar.gz && \
    cd kraken2-2.14 && \
    ./install_kraken2.sh $SRC/kraken2.14 && \
    ln -s $SRC/kraken2.14/kraken2* $BIN/  && \
    cd ..

WORKDIR $SRC
#install SPADes
RUN wget https://github.com/ablab/spades/releases/download/v4.1.0/SPAdes-4.1.0-Linux.tar.gz && \
    tar -xzf SPAdes-4.1.0-Linux.tar.gz  && \
    ln -s $SRC/SPAdes-4.1.0-Linux/bin/spades.py $BIN/spades.py

#install Trinity
WORKDIR $SRC
#RUN apt-get -qq -y install rsync  jellyfish salmon
RUN wget https://github.com/trinityrnaseq/trinityrnaseq/releases/download/Trinity-v2.15.1/trinityrnaseq-v2.15.1.FULL.tar.gz && \
    tar -xvzf trinityrnaseq-v2.15.1.FULL.tar.gz && \
    cd trinityrnaseq-v2.15.1  && \
    #sed -i '2 a #include <string>' trinity-plugins/bamsifter/sift_bam_max_cov.cpp && \
    #sed -i 's/HTSCODECS_VERSION_TEXT/HTSCODECS_VERSION/' trinity-plugins/bamsifter/htslib/htscodecs/htscodecs/htscodecs.c && \
    make && make plugins && \
    make install


#install snakemake
WORKDIR $SRC
RUN pip install Cython
RUN pip install setuptools>=61 
RUN pip install versioneer[toml]
RUN pip install snakemake  #==7.32.4
#RUN wget https://github.com/snakemake/snakemake/archive/refs/tags/v8.28.0.tar.gz && \
#    tar -xvzf v8.28.0.tar.gz && \
#    cd snakemake-8.28.0 && \
#    python setup.py install
RUN pip install pulp==2.7.0  #required by snakemake>7.0
#install pysam
RUN pip install pysam

RUN apt-get -qq -y remove git && \
    apt-get -qq -y autoremove && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /var/log/dpkg.log

RUN rm -rf $SRC/mummer-4.0.0rc1
RUN rm -rf $SRC/quast-5.2.0
RUN rm -rf $SRC/htslib-1.21
RUN rm -rf $SRC/samtools-1.21
RUN rm -rf $SRC/bcftools-1.21
RUN rm -rf $SRC/autoconf-2.69
RUN rm -rf $SRC/trinityrnaseq-v2.15.1
RUN rm -rf $SRC/mummer-4.0.0rc1
RUN rm -rf $SRC/kraken2-2.14

# some cleanup
WORKDIR $SRC
RUN rm -r *.tar.gz *.zip *.bz2

# copy PhytoPipe code
# phytopipe folder contains the application
ADD phytopipe $SRC/phytopipe

################################################
## be sure this is last!
COPY Dockerfile $SRC/phytopipe/.
