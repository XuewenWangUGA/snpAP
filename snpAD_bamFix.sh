#!/bin/bash

#SBATCH --job-name="snpAD_bamFix"

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

function display_help() {
    echo "Usage: $0 [options]"
    echo
    echo "Options:"
    echo "  -h, --help    Display this help message"
    echo "  -v, --version Display the version number"
	echo "  -i, --input   the input directory with all cram or bam files"
    echo "  -o, --output  the output directory to store all fixed bam files"
    echo "  -f, --format  the format of input alignment files: cram or bam, default [bam]"
	echo "  -r, --reference  the file name of the reference fasta, default [hg38]"
	echo "  -t, --thread  int, the number of paralelle computing threads, default [40]"	
	echo " e.g.: ./snpAD_bamFix.sh -i inputBamDir -o outFixBamDir-f alignment_format -r ref_genomeFasta_file -t threadNumber"
	echo "programmed by Xuewen Wang"
    echo
}

# Check if no arguments were passed
if [ $# -eq 0 ]; then
    echo "No arguments provided"
    display_help
    exit 1
fi


##inititalization
  #inAlignDir=$1
  inAlignDir=/eva/projects/research/genotyping_assessment/real_data
  #alignFormat=$2  #bam or cram
  alignFormat="bam"
  #refh38=$3  
  refh38=/eva/edatums/reference_materials/reference_genomes/hg38/GRCh38_full_analysis_set_plus_decoy_hla.fa  
  #outdir=$4
  outdir=/eva/home/xuewen/snpADanalysis/wrapperOut
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
                echo "Input directory is $inFile"
            else
                echo "No input directory specified"
                display_help
                exit 1
            fi
            ;;
        -o|--output)
            shift
            if [ -n "$1" ]; then
                outFile="$1"
                echo "Output directory is $outFile"
            else
                echo "No output directory specified"
                display_help
                exit 1
            fi
            ;;
        -r|--reference)
            shift
            if [ -n "$1" ]; then
                refh38="$1"
                echo "reference seqeunce is $outFile"
            else
                echo "no reference seqeunce specified"
                display_help
                exit 1
            fi
            ;;
        -f|--format)
            shift
            if [ -n "$1" ]; then
                alignFormat="$1"
                echo "alignment format is $outFile"
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

  
  #prepare outdir
  if [ ! -d $outdir ];then
      mkdir $outdir
  fi  
  
  #get list of cram and convert to bam
  cd $outdir
  alignFiles=""
  if [ $alignFormat == "cram" ] ;then
      echo "cram format is given; convert to bam now"
      alignFiles=`ls $inAlignDir/*.cram |tr "\n" ' '`
	  for cram in $alignFiles
	  do
		 cramFile=`basename $cram ".cram"`
		 bamfile=${cramFile}.bam
		 samtools view -@ $threadN -T $refh38 -b $cram > $outdir/$bamfile
		 samtools index -b  -@ $threadN $outdir/$bamfile
	  done	  
  elif [ $alignFormat == "bam" ];then
      echo "bam format is given"
      alignFiles=`ls $inAlignDir/*.bam |tr "\n" ' '`
	  echo $alignFiles
	  for bam in $alignFiles
	  do
	     bamlnk=`basename $bam ".bam"`
		 echo $bamlnk
		 ln -s $bam $bamlnk.bam
		 samtools index -b  -@ $threadN $bamlnk.bam
	  done
	  
  else
      echo "bam or cram format is needed"
  fi

 
   
   ## fixation the chr ID in bam head
   cd $outdir
   bamfixedFiles=$(ls *.bam|tr '/\n' ' '|sed 's/ *$//')
   
   for b in $bamfixedFiles
   do 
       bamInput=$b
       bamOut=${b/.bam/}.fixed.bam # .fixed.bam
       samtools view -H $bamInput >tmp.chr.22.bam.head
       awk 'BEGIN {OFS="\t"} /^@SQ/ {print $1,$2,$3; next} {print $0}' tmp.chr.22.bam.head  >tmp.chr.22.bam.head.fixed
       # @SQ     SN:chr3 LN:198295559 
       samtools view $bamInput >tmp.chr.22.nohead.sam
       cat tmp.chr.22.bam.head.fixed tmp.chr.22.nohead.sam >tmp.$bamOut
       samtools view tmp.$bamOut -b >$bamOut
       samtools index -b $bamOut
       rm tmp.*
	   printf "fixed bam: $bamOut"
   done
   
   