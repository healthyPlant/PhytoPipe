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

sudo apt-get -qq -y install python3-dev python3-pip
#rename python3
#sudo apt-get update && sudo apt-get -qq -y install python3
#sudo ln -sf /bin/python3 /bin/python

#install Parallel
sudo apt-get -qq -y install parallel

#install biopython
sudo pip install biopython

#install multiqc
sudo pip install multiqc

#install Cython
sudo pip install Cython

#apt-get install fastqc
cd $SRC
wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip
unzip fastqc_v0.12.1.zip 
chmod 755 $SRC/FastQC/fastqc
sudo ln -s $SRC/FastQC/fastqc $BIN/fastqc

#install Trimmomatic
wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip 
unzip Trimmomatic-0.39.zip

#install quast
wget https://github.com/ablab/quast/releases/download/quast_5.2.0/quast-5.2.0.tar.gz
tar -xvzf quast-5.2.0.tar.gz 
cd quast-5.2.0
sudo python setup.py install
cd ..


#install Cmake
wget https://github.com/Kitware/CMake/releases/download/v3.25.2/cmake-3.25.2-linux-x86_64.tar.gz
tar -xvzf cmake-3.25.2-linux-x86_64.tar.gz
sudo ln -s $SRC/cmake-3.25.2-linux-x86_64/bin/* $BIN/.

#install autoconf-2.69, this version is required by Trinity2.15
wget ftp://ftp.gnu.org/gnu/autoconf/autoconf-2.69.tar.gz 
tar zxvf autoconf-2.69.tar.gz 
cd autoconf-2.69 
./configure --prefix=/usr
make
sudo make install
cd ..

#install jellyfish
wget https://github.com/gmarcais/Jellyfish/releases/download/v2.3.1/jellyfish-2.3.1.tar.gz 
tar -xvzf jellyfish-2.3.1.tar.gz 
cd jellyfish-2.3.1 
./configure 
make
sudo make install
cd ..

#install salmon
wget https://github.com/COMBINE-lab/salmon/releases/download/v1.9.0/salmon-1.9.0_linux_x86_64.tar.gz && \
tar -xvzf salmon-1.9.0_linux_x86_64.tar.gz 
ln -s $SRC/salmon-1.9.0_linux_x86_64/bin/salmon $BIN/salmon    
    
#Install MUMmer4
wget https://github.com/mummer4/mummer/releases/download/v4.0.0rc1/mummer-4.0.0rc1.tar.gz
tar -xvzf mummer-4.0.0rc1.tar.gz
cd mummer-4.0.0rc1
./configure
make
sudo make install
sudo ldconfig
cd ..

#install BBtools
wget https://downloads.sourceforge.net/project/bbmap/BBMap_39.19.tar.gz
tar -xvzf BBMap_39.19.tar.gz
sudo ln -s $SRC/bbmap/*.sh $BIN/.

#install Bowtie2
#sudo apt-get -qq -y install bowtie2
wget https://github.com/BenLangmead/bowtie2/releases/download/v2.5.4/bowtie2-2.5.4-linux-x86_64.zip
unzip bowtie2-2.5.4-linux-x86_64.zip
mv bowtie2-2.5.4-linux-x86_64 bowtie2-2.5.4
sudo ln -s $SRC/bowtie2-2.5.4/bowtie2* $BIN/.

#install bedtools
mkdir bedtools2
cd bedtools2
wget https://github.com/arq5x/bedtools2/releases/download/v2.31.0/bedtools.static
mv bedtools.static bedtools
chmod a+x bedtools
cd ..
sudo ln -s $SRC/bedtools2/bedtools $BIN/bedtools

#install seqtk
git clone https://github.com/lh3/seqtk.git
cd seqtk
make
sudo cp $SRC/seqtk/seqtk $BIN/seqtk
cd ..

#Install BWA
git clone https://github.com/lh3/bwa.git
cd bwa; make
cd ..
sudo ln -s $SRC/bwa/bwa $BIN/bwa


#Install Krona
wget https://github.com/marbl/Krona/releases/download/v2.8.1/KronaTools-2.8.1.tar
tar -xvf KronaTools-2.8.1.tar
cd KronaTools-2.8.1
sudo ./install.pl  #command kt* to /usr/local/bin
sudo ln -s $SRC/KronaTools-2.8.1/updateAccessions.sh $BIN/updateAccessions.sh
sudo ln -s $SRC/KronaTools-2.8.1/updateTaxonomy.sh $BIN/updateTaxonomy.sh
cd ..

#instal NCBI EDirect
wget ftp://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/edirect.tar.gz
gunzip -c edirect.tar.gz | tar xf -
cd edirect
./nquire -dwn ftp.ncbi.nlm.nih.gov entrez/entrezdirect xtract.Linux.gz
gunzip -f xtract.Linux.gz
chmod +x xtract.Linux
cd ..
sudo find $SRC/edirect -type f -executable -exec ln -s {} $BIN/ \;


#install htslib
wget https://github.com/samtools/htslib/releases/download/1.21/htslib-1.21.tar.bz2
tar -xvjf htslib-1.21.tar.bz2
cd htslib-1.21 
autoheader 
autoconf -Wno-syntax 
./configure 
make
sudo make install
cd ..

#install samtools
wget https://github.com/samtools/samtools/releases/download/1.21/samtools-1.21.tar.bz2
tar -xvjf samtools-1.21.tar.bz2
cd samtools-1.21
autoheader 
autoconf -Wno-syntax 
./configure --with-htslib=/usr/local
make 
sudo make install
cd ..

#install bcftools
#sudo apt-get -qq -y install bcftools
wget https://github.com/samtools/bcftools/releases/download/1.21/bcftools-1.21.tar.bz2
tar -xvjf bcftools-1.21.tar.bz2
cd bcftools-1.21
autoheader 
autoconf -Wno-syntax 
./configure --with-htslib=/usr/local
make 
sudo make install
cd ..

#Installing NCBI-BLAST+
wget  https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.16.0+-x64-linux.tar.gz
tar -xvzf ncbi-blast-2.16.0+-x64-linux.tar.gz
sudo ln -s $SRC/ncbi-blast-2.16.0+/bin/* $BIN/.

#install Diamond
#https://diamond.readthedocs.io/en/latest/Getting-Started/Installation/
#sudo apt-get -qq -y install diamond-aligner
wget https://github.com/bbuchfink/diamond/releases/download/v2.1.11/diamond-linux64.tar.gz
tar -xvzf diamond-linux64.tar.gz
mkdir diamond2.1.11
mv diamond diamond2.1.11/
sudo ln -s $SRC/diamond2.1.11/diamond $BIN/diamond


#Install Kaiju
wget https://github.com/bioinformatics-centre/kaiju/archive/refs/tags/v1.10.1.tar.gz
tar -xvzf v1.10.1.tar.gz
cd kaiju-1.10.1/src
make
ln -s $SRC/kaiju-1.9.2/bin/* $BIN/ 
cd ..

#install recent version Kraken2
#sudo apt-get -qq -y install kraken2 #download library fail
wget https://github.com/DerrickWood/kraken2/archive/refs/tags/v2.14.tar.gz
tar -xvzf v2.14.tar.gz
cd kraken2-2.14
./install_kraken2.sh $SRC/kraken2.14
ln -s $SRC/kraken2.14/kraken2* $BIN/
cd ..

#install SPADes
wget https://github.com/ablab/spades/releases/download/v4.1.0/SPAdes-4.1.0-Linux.tar.gz
tar -xzf SPAdes-4.1.0-Linux.tar.gz
sudo ln -s $SRC/SPAdes-4.1.0-Linux/bin/spades.py $BIN/spades.py


#install Trinity
sudo apt-get -qq -y install rsync jellyfish salmon
wget https://github.com/trinityrnaseq/trinityrnaseq/releases/download/Trinity-v2.15.1/trinityrnaseq-v2.15.1.FULL.tar.gz
tar -xvzf trinityrnaseq-v2.15.1.FULL.tar.gz
cd trinityrnaseq-v2.15.1 
#sed -i '2 a #include <string>' trinity-plugins/bamsifter/sift_bam_max_cov.cpp
#sed -i 's/HTSCODECS_VERSION_TEXT/HTSCODECS_VERSION/' trinity-plugins/bamsifter/htslib/htscodecs/htscodecs/htscodecs.c
make && make plugins
sudo make install
cd ..

#install snakemake
#pip install snakemake
wget https://github.com/snakemake/snakemake/archive/refs/tags/v7.18.0.tar.gz 
tar -xvzf v7.18.0.tar.gz
cd snakemake-7.18.0
sudo python setup.py install
cd  ..

sudo pip install setuptools>=61 
sudo pip install versioneer[toml]
sudo pip install pulp==2.7.0  #required by snakemake>7.0

#install pandas
sudo pip install pandas
#install pysam
sudo pip install pysam
#install matplotlib
sudo pip install matplotlib

# some cleanup
rm -r *.tar.gz *.zip *.bz2
rm -rf $SRC/mummer-4.0.0rc1
rm -rf $SRC/quast-5.2.0
rm -rf $SRC/htslib-1.13
rm -rf $SRC/snakemake-7.18.0
rm -rf $SRC/autoconf-2.69
rm -rf $SRC/trinityrnaseq-v2.15.1
rm -rf $SRC/mummer-4.0.0rc1
rm -rf $SRC/kraken2-2.14
rm -rf $SRC/seqtk
rm -rf $SRC/KronaTools-2.8.1.tar
cd ..

sudo apt-get -qq -y autoremove
sudo apt-get clean 

#Check software versions
echo "Python3 version"
python --version
echo "BBtools version"
bbmap.sh version
echo "Java version"
java -version
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
echo "EDirect version"
esearch -h
