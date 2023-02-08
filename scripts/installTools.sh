#/usr/bin/env bash
# This script is for installing PhytoPipe required tools on a Ubuntu 20.04 machine

SRC=$1 #/opt

BIN=/usr/local/bin
sudo apt-get update  
sudo apt-get -qq -y install git curl wget unzip make apt-utils default-jdk build-essential libssl-dev zlib1g-dev libbz2-dev liblzma-dev libz-dev libfreetype6-dev libpng-dev autogen libtool shtool 

#sudo apt-get -qq -y install pkg-config 
#Configuring tzdata
#Geographic area: 2
#Time zone: 106

sudo apt-get -qq -y install python3-dev python3-pip python3-versioneer python3-matplotlib
#rename python3
#sudo apt-get update && sudo apt-get -qq -y install python3
#sudo ln -sf /bin/python3 /bin/python

#install Diamond
sudo apt-get -qq -y install diamond-aligner

#install Bowtie2
sudo apt-get -qq -y install bowtie2

#install Parallel
sudo apt-get -qq -y install parallel

#install Kraken2
sudo apt-get -qq -y install kraken2

#install biopython
pip install biopython

#install multiqc
pip install multiqc

#install Cython
pip install Cython

#apt-get install fastqc
cd $SRC
wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip && \
    unzip fastqc_v0.11.9.zip && \
    chmod 755 $SRC/FastQC/fastqc && \
    ln -s $SRC/FastQC/fastqc $BIN/.

#install Trimmomatic
wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip && \
    unzip Trimmomatic-0.39.zip && \

#install quast
wget https://github.com/ablab/quast/releases/download/quast_5.2.0/quast-5.2.0.tar.gz && \
    tar -xvzf quast-5.2.0.tar.gz && \
    cd quast-5.2.0 && \
    python setup.py install
cd ..

#install snakemake
wget https://github.com/snakemake/snakemake/archive/refs/tags/v7.18.0.tar.gz && \
    tar -xvzf v7.18.0.tar.gz && \
    cd snakemake-7.18.0 && \
    python setup.py install
cd  ..

#install Cmake
wget https://github.com/Kitware/CMake/releases/download/v3.25.2/cmake-3.25.2-linux-x86_64.tar.gz && \
    tar -xvzf cmake-3.25.2-linux-x86_64.tar.gz && \
    cd cmake-3.25.2-linux-x86_64
sudo ln -s $SRC/cmake-3.25.2-linux-x86_64/bin/* $BIN/.
cd ..

#install autoconf-2.69, this version is required by Trinity2.15
wget ftp://ftp.gnu.org/gnu/autoconf/autoconf-2.69.tar.gz && \
    tar zxvf autoconf-2.69.tar.gz && \
    cd autoconf-2.69 && \
    ./configure --prefix=/usr/ && \
    make && \
    sudo make install
cd ..

#install htslib
wget https://github.com/samtools/htslib/releases/download/1.16/htslib-1.16.tar.bz2 && \
    tar -xvjf htslib-1.16.tar.bz2 && \
    cd htslib-1.16 && \
    autoheader  && \
    autoconf -Wno-syntax  && \
    ./configure && \
    make && \
    sudo make install
cd ..

#install samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.16.1/samtools-1.16.1.tar.bz2 && \
    tar -xvjf samtools-1.16.1.tar.bz2 && \
    cd samtools-1.16.1 && \
    autoheader  && \
    autoconf -Wno-syntax  && \
    ./configure --with-htslib=/usr/local && \
    make  && \
    sudo make install
cd ..


#RUN apt-get -qq -y install bcftools
#install bcftools
#sudo apt-get -qq -y install bcftools
wget https://github.com/samtools/bcftools/releases/download/1.16/bcftools-1.16.tar.bz2 && \
    tar -xvjf bcftools-1.16.tar.bz2 && \
    cd bcftools-1.16 && \
    autoheader  && \
    autoconf -Wno-syntax  && \
    ./configure --with-htslib=/usr/local && \
    make  && \
    sudo make install
cd ..

#install Trinity
sudo apt-get -qq -y install rsync jellyfish salmon
wget https://github.com/trinityrnaseq/trinityrnaseq/releases/download/Trinity-v2.15.0/trinityrnaseq-v2.15.0.FULL.tar.gz && \
    tar -xvzf trinityrnaseq-v2.15.0.FULL.tar.gz && \
    cd trinityrnaseq-v2.15.0  && \
    sed -i '2 a #include <string>' trinity-plugins/bamsifter/sift_bam_max_cov.cpp && \
    sed -i 's/HTSCODECS_VERSION_TEXT/HTSCODECS_VERSION/' trinity-plugins/bamsifter/htslib/htscodecs/htscodecs/htscodecs.c && \
    make && make plugins && \
    sudo make install
cd ..

#Install Kaiju
git clone https://github.com/bioinformatics-centre/kaiju.git && \
    cd $SRC/kaiju/src && \
    make
sudo cp ../bin/* $BIN/
cd ..

#Install MUMmer4
wget https://github.com/mummer4/mummer/releases/download/v4.0.0rc1/mummer-4.0.0rc1.tar.gz && \
    tar -xvzf mummer-4.0.0rc1.tar.gz && \
    cd mummer-4.0.0rc1 && \
    ./configure && \
    make && \
    sudo make install
ldconfig
cd ..

#Install Krona
wget https://github.com/marbl/Krona/releases/download/v2.8.1/KronaTools-2.8.1.tar && \
    tar -xvf KronaTools-2.8.1.tar && \
    cd KronaTools-2.8.1
sudo ./install.pl
sudo ln -s $SRC/KronaTools-2.8.1/updateAccessions.sh $BIN/updateAccessions.sh
sudo ln -s $SRC/KronaTools-2.8.1/updateTaxonomy.sh $BIN/updateTaxonomy.sh
rm -rf $SRC/KronaTools-2.8.1.tar
cd ..

#Installing NCBI-BLAST+
#sudo apt-get -qq -y install ncbi-blast+  #version > 2.10. v2.9 doesn't work
wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.13.0+-x64-linux.tar.gz && \
    tar -xvzf ncbi-blast-2.13.0+-x64-linux.tar.gz
sudo ln -s $SRC/ncbi-blast-2.13.0+/bin/* $BIN/.

#install BBtools
#RUN apt-get -qq -y install bbmap  #Ubuntu BBmap doesn't work
wget https://downloads.sourceforge.net/project/bbmap/BBMap_39.01.tar.gz && \
    tar -xvzf BBMap_39.01.tar.gz
sudo ln -s $SRC/bbmap/bbmap.sh $BIN/bbmap.sh
sudo ln -s $SRC/bbmap/bbduk.sh $BIN/bbduk.sh
sudo ln -s $SRC/bbmap/clumpify.sh $BIN/clumpify.sh
sudo ln -s $SRC/bbmap/bbsplit.sh $BIN/bbsplit.sh

#install SPADes
wget http://cab.spbu.ru/files/release3.15.5/SPAdes-3.15.5-Linux.tar.gz && \
    tar -xzf SPAdes-3.15.5-Linux.tar.gz
sudo ln -s $SRC/SPAdes-3.15.5-Linux/bin/spades.py $BIN/spades.py

#install bedtools
mkdir bedtools2
cd bedtools2
wget https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools.static.binary
mv bedtools.static.binary bedtools
chmod a+x bedtools
sudo ln -s $SRC/bedtools2/bedtools $BIN/bedtools
cd ..

#install seqtk
git clone https://github.com/lh3/seqtk.git
cd seqtk
make
sudo ln -s $SRC/seqtk/seqtk $BIN/seqtk

wget ftp://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/edirect.tar.gz && \
    gunzip -c edirect.tar.gz | tar xf -
cd $SRC/edirect
./nquire -dwn ftp.ncbi.nlm.nih.gov entrez/entrezdirect xtract.Linux.gz
gunzip -f xtract.Linux.gz
chmod +x xtract.Linux
cd ..
find $SRC/edirect -type f -executable -exec ln -s {} $BIN/ \;

#install pandas
pip install pandas
#install psam
pip install pysam
#install matplotlib
pip install matplotlib

# some cleanup
rm -r *.tar.gz *.zip *.bz2
rm -rf $SRC/mummer-4.0.0rc1
rm -rf $SRC/quast-5.2.0
rm -rf $SRC/htslib-1.13
rm -rf $SRC/snakemake-7.18.0
rm -rf $SRC/autoconf-2.69
rm -rf $SRC/trinityrnaseq-v2.15.0
rm -rf $SRC/mummer-4.0.0rc1
cd ..

sudo apt-get -qq -y autoremove
sudo apt-get clean 
sudo rm -rf /var/lib/apt/lists/* /var/log/dpkg.log


#Check software versions
echo "Python3 version"
python --version
echo "BBtools version"
bbmap.sh version
echo "Java version"
java --version
echo "FastQC version"
fastqc -v 
echo "Trimmomatic version" 
java -jar $SRC/Trimmomatic-0.39/trimmomatic-0.39.jar -version
echo "QUAST version" 
quast.py --version
echo "snakemake version"
snakemake -v
echo "cmake version"
cmake --version
echo "SPAdes version"
spades.py --version
echo "multiqc version"
multiqc --version
echo "Kraken2 version"
kraken2 --version
echo "Kaiju version"
kaiju
echo "MUMmer4 version"
mummer -v
echo "Blast version"
blastn -version
echo "Diamond version"
diamond --version
echo "BWA version"
bwa
echo "Bowtie2 version"
bowtie2 --version
echo "Bcftools version"
bcftools --version
echo "samtools version"
samtools version
echo "bedtools version"
bedtools --version
echo "Seqtk version"
seqtk
echo "Trinity version"
Trinity --version
echo "Krona version"
ls $SRC | grep Krona
