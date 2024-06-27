
# SnpAP for snp Allele Probability

The snpAlleleProbability (SnpAP) is a bioinformatic tool to calculate the likelyhood of snp alleles

# Run
Type this command in command terminal:

`python3 snpAlleleProbability.py`

# Usage: 

  python3 snpAlleleProbability.py [options]
    Options:
    -h, --help                  Show this help message and exit
    -a, --alignment=FILE        Input BAM file for alignment, the bam index file .bai should be present
    -t, --thread=NUMBER         Number of threads to use,default [12]
    -r, --reference=FILE        Reference genome file in FASTA format
    -b, --bed=FILE              1-based position file with regions of interest
    -o, --outprefix=PREFIX      Output prefix for result files
    -c, --chromosomeFile=FILE   Chromosome list, tab separtated text file, e.g. chr1 chr2 chr3
    -q, --quality=INT           The qulaity of mapped reads,default [25]
    version 1.0


# Input

1. SNP target site file, 1-based, tab separated four columns for each variant site: #Chr	Position_1-based	Allele_1	Allele_2
2. alignment in cram or bam format
3. reference sequence, e.g. hg38


# Output
 1.  a tsv file with likelyhood for each targeted site. 
  #Chr	Position_1-based	Allele_1	Allele_2	Likelihood_Hom1:Het:Hom2
  chr21	5034230	C	T	0:30:44
  chr21	5034244	C	T	0:33:89
  chr21	5034245	G	A	.
  chr21	5034257	G	T	0:41:133
  chr22	10572984	C	G	.

2. .vcf the output from snpAD
   

# Dependencies: 
1. python3
2. samtools
3. snpAD an ancient-DNA aware genotyper


# Installation

Download this tool from github via .zip by click on the webpage or the command:

`git clone https://github.com/XuewenWangUGA/snpAP`

# Installation of SnpAD:
The following is the snpAD. The user needs to compile snp and them put the compliled three files in the same directory of snpAlleleProbability.

This version of snpAD 0.3.5 is cloned from https://bioinf.eva.mpg.de/snpAD/ and tested to run our tool. 

-------------

To compile, run make. You'll need:
- g++ with support for C11
- the optimization library nlopt
- the commandline parser lib popt
- libRmath

If you are using Ubuntu, you can install the necessary libraries using the command:

apt-get install libnlopt-dev libpopt-dev r-mathlib

To install to /usr/local run 

make install

If you prefer the binaries in another subdirectory try

PREFIX=/my/dir make install


Running:
--------

To run snpAD, you need to follow three steps:
1) Prepare an Inputfile
2) Run Estimator 
3) Write VCF file

== Preparing the genotype input files ==

The input files with the sequence data can be generated with
Bam2snpAD/Bam2snpAD. The program takes an indexed bam file, and can filter on
the fly for mapping quality and base quality of the sequences in the bam and
for read-groups. The output will consist of one line per genomic position
showing chromosome, position, reference base followed by the information about
each read covering this position. This information per read gives the strand,
the base and an "error profile" number. 

In the simplest case, when you just have one type of library, you can generate
an inputfile that assigns separate error profiles for the first and last 15
positions in the read and one error profile for any position in the interior of
the read by running this commandline:

Bam2snpAD -f chr1.fa -i chr1.bai -r $i -s 15 -o 0 chr1.bam > genotyper_input/chr1.snpAD

The example commandline below generates two outputfiles with an error profile
0-14 for the first 15 positions from the 5' end, then the first 15 positions
from the 3' end, followed by a profile for the middle. The first outputfile
assigns the error profiles 0-14, 15-29, 30 to those three classes of positions,
the second outputfile assigns 31-45, 46-60, 61  to these three classes of
positions (i.e. it adds an offset of 31). 

Bam2snpAD -f chr1.fa -i chr1.bai -r $i -Q 25 -q 30 -s 15 -R Loschbour -o 0 chr1.bam > genotyper_input/out1.snpAD
Bam2snpAD -f chr1.fa -i chr1.bai -r $i -Q 25 -q 30 -s 15 -R Loschbour1,Loschbour2,Loschbour3,Loschbour4 -o 31 chr1.bam > genotyper_input/out2.snpAD

By merging these two files (using Bam2snpAD/snpADjoin out1.snpAD out2.snpAD >
both.snpAD) you generate an input-file that fits the above error profile with
entries for two types of libraries (identified by read-groups Loschbour and
Loschbour[1234], respectively) and 15 positions from either end + middle base.

== Running the Estimator ==

This step estimates the priors for each genotype and the error profile using an
EM algorithm. Here is an example commandline that uses 40 cpus and prints
genotype priors and error profile to the files priors.txt and errors.txt,
respectively. 

snpAD -c 40 -o priors.txt -O errors.txt input.snpAD > log.tab 2> log.err

== Write the VCF file ==

You can then use snpADcall to convert the snpAD inputfile to a VCF file:

snpADCall -N Sample-Name -e errors.txt -p "`cat priors.txt`" input.snpAD > output.VCF


License:
--------

snpAD is free software (GNU General Public License v3). Note that this package includes 
