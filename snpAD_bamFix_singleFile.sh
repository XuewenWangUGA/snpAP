#!/bin/bash

#SBATCH --job-name="snpAD_bamFix_singleFile"

##SBATCH -p ss

##SBATCH -N 1

#SBATCH --nodes=1

#SBATCH --ntasks=1

##SBATCH -n 40
##SBATCH --mem=500G

#SBATCH --cpus-per-task=1

#SBATCH -t 2-00:00:00

#SBATCH --mail-user=xuewen.wang@unthsc.edu

#SBATCH --mail-type=BEGIN,END,FAIL

#SBATCH -o "%x.%j.o"

#SBATCH -e "%x.%j.e"


## squeue -u xuewen
## squeue -j your_job_id

#scriptdir=/eva/home/xuewen/snpADanalysis/scripts
#working dir
#cd /eva/home/xuewen/snpADanalysis

#  bash snpAD_bamFix_singleFile.sh -i /eva/home/xuewen/snpADanalysis/realData/SC_GMWOF5428715.alt_bwamem_GRCh38DH.20151208.WOLOFF.gambian_lowcov.bam -f bam -o /eva/home/xuewen/snpADanalysis/realData/SC_GMWOF5428715.alt_bwamem_GRCh38DH.20151208.WOLOFF.gambian_lowcov.fixed.bam

function display_help() {
    echo "Usage: $0 [options]"
    echo
    echo "Options:"
    echo "  -h, --help    Display this help message"
    echo "  -v, --version Display the version number"
	echo "  -i, --input   the input name with a cram or bam file"
    echo "  -o, --output  the output file name to store the fixed bam file"
    echo "  -f, --format  the format of input alignment files: cram or bam, default [bam]"
	echo "  -r, --reference  the file name of the reference fasta, default [hg38]. required for cram input only"
	echo "  -t, --thread  int, the number of paralelle computing threads, default [40]"	
    echo "version 1.0, July,5th 2024. "
	echo "Support: Contact Xuewen Wang"
    echo
}

# Check if no arguments were passed
if [ $# -eq 0 ]; then
    echo "No arguments provided"
    display_help
    exit 1
fi


##test command:
# -i /eva/home/xuewen/snpADanalysis/realData/

##inititalization
  #inAlignDir=$1
  inAlignFile=/eva/home/xuewen/snpADanalysis/realData/SC_GMWOF5428715.alt_bwamem_GRCh38DH.20151208.WOLOFF.gambian_lowcov.bam
  #alignFormat=$2  #bam or cram
  alignFormat="bam"
  #refh38=$3  
  refh38=/eva/home/xuewen/snpADanalysis/realData/GRCh38_full_analysis_set_plus_decoy_hla.fa  
  #outdir=$4
  outFile=/eva/home/xuewen/snpADanalysis/realData/SC_GMWOF5428715.alt_bwamem_GRCh38DH.20151208.WOLOFF.gambian_lowcov.fixed.bam
  threadN=40


# Parse command line arguments
while (( "$#" )); do
    case "$1" in
        -h|--help)
            display_help
            exit 0
            ;;
        -v|--version)
            echo "Version 1.0"
            exit 0
            ;;
        -i|--input)
            shift
            if [ -n "$1" ]; then
                inFile="$1"
                echo "Input File is $inFile"
            else
                echo "No input alignment specified"
                display_help
                exit 1
            fi
            ;;
        -o|--output)
            shift
            if [ -n "$1" ]; then
                outFile="$1"
                echo "Output alignment file is $outFile"
            else
                echo "No output alignmeng file specified"
                display_help
                exit 1
            fi
            ;;
        -r|--reference)
            shift
            if [ -n "$1" ]; then
                refh38="$1"
                echo "reference seqeunce is $refh38"
            else
                echo "no reference seqeunce specified"
                #display_help
                #exit 1
            fi
            ;;
        -f|--format)
            shift
            if [ -n "$1" ]; then
                alignFormat="$1"
                echo "alignment format is $alignFormat"
            else
                echo "No alignment format specified"
                display_help
                exit 1
            fi
            ;;
        -t|--thread)
            shift
            if [ -n "$1" ]; then
                threadN="$1"
                echo "thread is $outFile"
            else
                echo "No thread specified"
                display_help
                exit 1
            fi
            ;;			
        *)
            echo "Unknown option: $1"
            display_help
            exit 1
            ;;
    esac
    shift
done
 

  #get list of cram and convert to bam
  alignFiles=""
  if [ $alignFormat == "cram" ] ;then
      echo "cram format is given; convert to bam now"
         cramFile=$inFile
		 cramFile=`basename $cram ".cram"`
		 bamfile=${cramFile}.bam
		 samtools view -@ $threadN -T $refh38 -b $cram > $bamfile
		 samtools index -b  -@ $threadN $bamfile
  
  elif [ $alignFormat == "bam" ];then
      echo "bam format is given"
         bamfile=$inFile
		 if [ ! -f $bamfile ];then
		     samtools index -b  -@ $threadN $bamfile
		 fi
  else
      echo "bam or cram format is needed"
  fi

   
   ## fixation 
   echo "running the fixation of the chr ID in bam head"
 
       bamInput=$bamfile
       bamOut=${bamfile/.bam/}.fixed.bam # .fixed.bam
       samtools view -@ $threadN -H $bamInput >$bamInput.head.tmp
       awk 'BEGIN {OFS="\t"} /^@SQ/ {print $1,$2,$3; next} {print $0}' $bamInput.head.tmp  >$bamInput.head.fixed.tmp
       # @SQ     SN:chr3 LN:198295559 
       samtools view -@ $threadN $bamInput >${bamInput/.bam/}.nohead.sam.tmp
       cat $bamInput.head.fixed.tmp ${bamInput/.bam/}.nohead.sam.tmp >$bamOut.tmp
       samtools view -@ $threadN $bamOut.tmp -b >$bamOut
       samtools index -@ $threadN -b $bamOut
       rm ${bamInput/.bam/}*.tmp
	   printf "fixed bam: $bamOut"

   
   