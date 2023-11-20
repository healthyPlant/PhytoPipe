#/usr/bin/env bash
# This script is for installing PhytoPipe required tools on a Ubuntu 20.04 machine

SRC=$1 #/opt

BIN=/usr/local/bin

#************************************
#This is for Ubuntu libary installing,
#For other systems, please find the related libaries.
sudo apt-get update  
sudo apt-get -qq -y install git curl wget unzip make apt-utils default-jdk build-essential libssl-dev zlib1g-dev libbz2-dev liblzma-dev libz-dev libfreetype6-dev libpng-dev autogen libtool shtool cpanm

#sudo apt-get -qq -y install pkg-config 
#Configuring tzdata
#Geographic area: 2
#Time zone: 106

sudo apt-get -qq -y install python3-dev python3-pip python3-versioneer
#rename python3
#sudo apt-get update && sudo apt-get -qq -y install python3
#sudo ln -sf /bin/python3 /bin/python

#********************************************************
#The following codes should be for Linux or Mac.

#install biopython
sudo pip install biopython

#install multiqc
sudo pip install multiqc

#install Cython
sudo pip install Cython

#apt-get install fastqc
cd $SRC
wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip
unzip fastqc_v0.11.9.zip 
chmod 755 $SRC/FastQC/fastqc
sudo ln -s $SRC/FastQC/fastqc $BIN/fastqc

#install Trimmomatic
wget https://github.com/usadellab/Trimmomatic/files/5854859/Trimmomatic-0.39.zip 
unzip Trimmomatic-0.39.zip

#install quast
wget https://github.com/ablab/quast/releases/download/quast_5.2.0/quast-5.2.0.tar.gz
tar -xvzf quast-5.2.0.tar.gz 
cd quast-5.2.0
sudo python setup.py install
cd ..

#install snakemake
wget https://github.com/snakemake/snakemake/archive/refs/tags/v7.18.0.tar.gz 
tar -xvzf v7.18.0.tar.gz
cd snakemake-7.18.0
sudo python setup.py install
cd  ..

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

#install htslib
wget https://github.com/samtools/htslib/releases/download/1.16/htslib-1.16.tar.bz2
tar -xvjf htslib-1.16.tar.bz2
cd htslib-1.16 
autoheader 
autoconf -Wno-syntax 
./configure 
make
sudo make install
cd ..

#install samtools
wget https://github.com/samtools/samtools/releases/download/1.16.1/samtools-1.16.1.tar.bz2
tar -xvjf samtools-1.16.1.tar.bz2
cd samtools-1.16.1
autoheader 
autoconf -Wno-syntax 
./configure --with-htslib=/usr/local
make 
sudo make install
cd ..

#install bcftools
#sudo apt-get -qq -y install bcftools
wget https://github.com/samtools/bcftools/releases/download/1.16/bcftools-1.16.tar.bz2
tar -xvjf bcftools-1.16.tar.bz2
cd bcftools-1.16
autoheader 
autoconf -Wno-syntax 
./configure --with-htslib=/usr/local
make 
sudo make install
cd ..

#install Trinity
#install perl module DB_file
sudo apt-get install libdb-dev
sudo cpanm install DB_File

sudo apt-get -qq -y install rsync jellyfish salmon
wget https://github.com/trinityrnaseq/trinityrnaseq/releases/download/Trinity-v2.15.1/trinityrnaseq-v2.15.1.FULL.tar.gz
tar -xvzf trinityrnaseq-v2.15.1.FULL.tar.gz
cd trinityrnaseq-v2.15.1 
make && make plugins
sudo make install
cd ..

#Install Kaiju
git clone https://github.com/bioinformatics-centre/kaiju.git
cd kaiju/src
make
sudo cp ../bin/* $BIN/
cd ../..

#Install MUMmer4
wget https://github.com/mummer4/mummer/releases/download/v4.0.0rc1/mummer-4.0.0rc1.tar.gz
tar -xvzf mummer-4.0.0rc1.tar.gz
cd mummer-4.0.0rc1
./configure
make
sudo make install
sudo ldconfig
cd ..

#Install Krona
wget https://github.com/marbl/Krona/releases/download/v2.8.1/KronaTools-2.8.1.tar
tar -xvf KronaTools-2.8.1.tar
cd KronaTools-2.8.1
sudo ./install.pl  #command kt* to /usr/local/bin
sudo ln -s $SRC/KronaTools-2.8.1/updateAccessions.sh $BIN/updateAccessions.sh
sudo ln -s $SRC/KronaTools-2.8.1/updateTaxonomy.sh $BIN/updateTaxonomy.sh
cd ..

#Installing NCBI-BLAST+
#sudo apt-get -qq -y install ncbi-blast+  #version > 2.10. v2.9 doesn't work
wget  https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.14.0/ncbi-blast-2.14.0+-x64-linux.tar.gz
tar -xvzf ncbi-blast-2.14.0+-x64-linux.tar.gz
sudo ln -s $SRC/ncbi-blast-2.14.0+/bin/* $BIN/.

#install BBtools
#RUN apt-get -qq -y install bbmap  #Ubuntu BBmap doesn't work
wget https://downloads.sourceforge.net/project/bbmap/BBMap_39.01.tar.gz
tar -xvzf BBMap_39.01.tar.gz
sudo ln -s $SRC/bbmap/*.sh $BIN/.

#install SPADes
wget http://cab.spbu.ru/files/release3.15.5/SPAdes-3.15.5-Linux.tar.gz
tar -xzf SPAdes-3.15.5-Linux.tar.gz
sudo ln -s $SRC/SPAdes-3.15.5-Linux/bin/spades.py $BIN/spades.py

#install bedtools
mkdir bedtools2
cd bedtools2
wget https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools.static.binary
mv bedtools.static.binary bedtools
chmod a+x bedtools
cd ..
sudo ln -s $SRC/bedtools2/bedtools $BIN/bedtools

#install seqtk
git clone https://github.com/lh3/seqtk.git
cd seqtk
make
sudo cp $SRC/seqtk/seqtk $BIN/seqtk
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

#install recent version Kraken2
#sudo apt-get -qq -y install kraken2 #download library fail
git clone https://github.com/DerrickWood/kraken2.git
cd kraken2
./install_kraken2.sh $SRC/kraken2.1.2
sudo cp $SRC/kraken2.1.2/kraken2* $BIN/
cd ..

#install Bowtie2
#sudo apt-get -qq -y install bowtie2
wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.4.2/bowtie2-2.4.2-sra-linux-x86_64.zip
unzip bowtie2-2.4.2-sra-linux-x86_64.zip
mv bowtie2-2.4.2-sra-linux-x86_64 bowtie2-2.4.2
sudo ln -s $SRC/bowtie2-2.4.2/bowtie2* $BIN/.

#Install BWA
git clone https://github.com/lh3/bwa.git
cd bwa; make
cd ..
sudo ln -s $SRC/bwa/bwa $BIN/bwa

#install Diamond
#https://diamond.readthedocs.io/en/latest/Getting-Started/Installation/
#sudo apt-get -qq -y install diamond-aligner
wget https://github.com/bbuchfink/diamond/releases/download/v2.1.1/diamond-linux64.tar.gz
tar -xvzf diamond-linux64.tar.gz
mkdir diamond2.1.1
mv diamond diamond2.1.1/
sudo ln -s $SRC/diamond2.1.1/diamond $BIN/diamond

#install Parallel
#sudo apt-get -qq -y install parallel
wget --no-check-certificate https://ftp.gnu.org/gnu/parallel/parallel-20230222.tar.bz2
tar -xvjf parallel-latest.tar.bz2
cd parallel-20230222
./configure --prefix=/ppq/data1/software/parallel
make
make install
sudo ln -s $SRC/parallel/bin/parallel $BIN/parallel

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
rm -rf $SRC/trinityrnaseq-v2.15.0
rm -rf $SRC/mummer-4.0.0rc1
rm -rf $SRC/kraken2
rm -rf $SRC/seqtk
rm -rf $SRC/KronaTools-2.8.1.tar
rm -rf $SRC/kaiju
rm -rf $SRC/parallel-20230222
cd ..

sudo apt-get -qq -y autoremove
sudo apt-get clean 

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
echo "EDirect version"
esearch -h
