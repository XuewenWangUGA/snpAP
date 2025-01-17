
# SnpAP for snp Allele Probability

The snpAlleleProbability (SnpAP) is a bioinformatics tool to calculate the likelihood of SNP alleles

# Run
Type this command in a command terminal:

`python3 snpAlleleProbability.py`

# Usage: 

        Usage: python3 snpAlleleProbability.py [options]
    Options:
    -h, --help                  Show this help message and exit
    -a, --alignment=FILE        Input BAM file for alignment, the bam index file .bai should be present
    -t, --thread=NUMBER         Number of threads to use, default [12]
    -r, --reference=FILE        Reference genome file in FASTA format
    -b, --bed=FILE              BED file with regions of interest
    -o, --outprefix=PREFIX      Output prefix for result files
    -c, --chromosomeFile=FILE   Chromosome list, tab separate text file, e.g. chr1 chr2 chr3
    -q, --quality=INT           The quality of mapped reads, default [25]
    -m, --minqual=INT           the minimal quality score of the base for analysis [10]
    -n, --normal                T|F. T will normalize the most likely PL to 0, default [F]    
    version 1.0



# Input

1. SNP target site file, 1-based, a tab-separated four columns for each variant site:
           #Chr	Position_1-based	Allele_1	Allele_2
3. alignment in cram or bam format
4. reference sequence, e.g. hg38


# Output
 1.  a plain text  file in the standard VCF format with likelihood for each targeted site. The smaller the value the higher confidence of the allele.


                ##fileformat=VCFv4.1
                ##FORMAT=<ID=A,Number=2,Type=Integer,Description="Number of A bases on forward and reverse strand">
                ##FORMAT=<ID=C,Number=2,Type=Integer,Description="Number of C bases on forward and reverse strand">
                ##FORMAT=<ID=G,Number=2,Type=Integer,Description="Number of G bases on forward and reverse strand">
                ##FORMAT=<ID=T,Number=2,Type=Integer,Description="Number of T bases on forward and reverse strand">
                ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
                ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
                ##FORMAT=<ID=GQ,Number=1,Type=Float,Description="Genotype Quality (difference between lowest and second lowest PP value)">
                ##FORMAT=<ID=PL,Number=3,Type=Integer,Description="Phred-scaled likelihood for genotypes hom1, het, and hom2">
                #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	gambian_lowcov
                chr21	5034230	.	C	T	.	.	.	GT:PL	0/0:0,30,44
                chr21	5034244	.	C	T	.	.	.	GT:PL	0/0:0,33,89
                chr21	5034245	.	G	A	.	.	.	GT:PP	0/0:.,.,.
                chr21	5034257	.	G	T	.	.	.	GT:PL	0/0:0,41,133
                chr22	10572984	.	C	G	.	.	.	GT:PP	0/0:.,.,.


3. .vcf the output from snpAD
   

# Dependencies: 
1. python3
2. samtools
3. snpAD: an ancient-DNA-aware genotyper


# Installation

Download this tool from GitHub via .zip by clicking on the webpage or the command:

`git clone https://github.com/XuewenWangUGA/snpAP`


# Utility: snpAD_bamFix

The snpAD_bamFix will fix issues that will happen when calling snp with the original snpAD. To run the snpAD_bamFix first and then use the fixed bam for snpAP.

This tool takes .cram or .bam as input.

to run, type this in a command line:
`./snpAD_bamFix.sh`

    
        Usage: ./snpAD_bamFix.sh [options]
        
        Options:
          -h, --help    Display this help message
          -v, --version Display the version number
          -i, --input   an input directory with cram or bam files
          -o, --output  an output directory to store fixed bam files
          -f, --format  the format of input alignment files: cram or bam, default [bam]
          -r, --reference  the file name of the reference fasta, default [hg38]
          -t, --thread  int, the number of parallel computing threads, default [40]
         e.g.: ./snpAD_bamFix.sh -i inputBamDir -o outFixBamDir-f alignment_format -r ref_genomeFasta_file -t threadNumber


# Installation of SnpAD:
The following is the snpAD. The user needs to compile snp and then put the compliled three files in the same directory of snpAlleleProbability.

This version of snpAD 0.3.5 is cloned from https://bioinf.eva.mpg.de/snpAD/ and tested to run our tool. 

-------------

To compile, run make. You'll need:
- g++ with support for C11
- the optimization library nlopt
- the command line parser lib popt
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

--------


# Developement

This tool is under activate development, programmed by Xuewen Wang. 
