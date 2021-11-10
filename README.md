### Overview

&emsp; &emsp; MHtyper is an end-to-end pipeline for recognized the Forensic microhaplotypes in Nanopore  sequencing data. It
is implemented using Python.

### MHtyper workflow

#### Step1：

&emsp;Sequencing data was filtered by NanoFilt and aligment with minimap2

#### Step2：

&emsp;Files with both BED and pileup format are generated

#### Step3：

&emsp;Phasing with margin ; Correction with isONcorrect; Haplotype analysis by MHtyper ; Integrate the analysis results
to get the final micro haplotype results

### Python environment construction and required software installation

```
   conda create -n MHtyper 
   conda activate MHtyper
   conda config --add channels bioconda 
   conda config --add channels
   conda-forge conda install -y NanoFilt minimap2 samtools bedtools
```
### isONcorrect & Margin installation
&emsp;
            isONcorrect: https://github.com/ksahlin/isONcorrect [1]
            Margin: https://github.com/UCSC-nanopore-cgl/margin

```
# isONcorrect installation

   git clone https://github.com/ksahlin/isONcorrect.git
   cd isONcorrect 
   ./isONcorrect

# Marin installation
# step1:
   sudo apt-get install git make gcc g++ autoconf zlib1g-dev libcurl4-openssl-dev libbz2-dev libhdf5-dev
   wget https://github.com/Kitware/CMake/releases/download/v3.14.4/cmake-3.14.4-Linux-x86_64.sh && sudo mkdir /opt/cmake &&
   sudo sh cmake-3.14.4-Linux-x86_64.sh --prefix=/opt/cmake --skip-license && sudo ln -s /opt/cmake/bin/cmake
   /usr/local/bin/cmake cmake --version

# step2: Check out the repository and submodules:

   git clone https://github.com/UCSC-nanopore-cgl/margin.git
   cd margin git submodule update --init

# step3: Make build directory:

   mkdir build cd build

# step4: Generate Makefile and run:

   cmake .. 
   make ./margin
```

### MHtyper installation
```
   git clone https://github.com/willow2333/MHtyper.git
   cd MHtyper 
   python run.py --h
   
   usage: run.py [-h] [--fastqfiles FASTQFILES] [--reference REFERENCE] [--prefix PREFIX] [--truthvcf TRUTHVCF] [--marginpath MARGINPATH]

    optional arguments:
      -h, --help            show this help message and exit
      --fastqfiles FASTQFILES
                            The input *.fq.gz files.
      --reference REFERENCE
                            The path of your ref.
      --prefix PREFIX       The name of your Sample, default is "Test".
      --truthvcf TRUTHVCF   The truth variant files in your research.
      --marginpath MARGINPATH
                            The setup path of "Margin".
```

###  Illustration
#### 1.Test
```
   cd ./Test
   python ../run.py --fastqfiles test.fq.gz --reference path/hg19.fa --prefix Test --truthvcf truthvcf.txt  --marginpath path/margin
```
#### 2. The sites vcf files needed
&emsp;The snp-sites.txt that contained the information of samples must needed 
#### 3. Output
&emsp;The analysis results of microhaplotypes is in finalphase.txt

### Citation
&emsp;&emsp;1.Sahlin, K., Medvedev, P. Error correction enables use of Oxford Nanopore technology for reference-free transcriptome analysis. Nat Commun 12, 2 (2021). https://doi.org/10.1038/s41467-020-20340-8 Link.


© 2021 by  Yiping Hou (forensic@scu.edu.cn), Zheng Wang (wangzhengtim@scu.edu.cn), Liu Qin (liu.qin@qitantech.com)

