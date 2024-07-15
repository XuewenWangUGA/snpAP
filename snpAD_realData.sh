#!/bin/bash

#SBATCH --job-name="snpAD_realData"

##SBATCH -p ss

##SBATCH -N 1

#SBATCH --nodes=1

#SBATCH --ntasks=1

##SBATCH -n 38
##SBATCH --mem=500G

#SBATCH --cpus-per-task=1

#SBATCH -t 2-00:00:00

#SBATCH --mail-user=xuewen.wang@unthsc.edu

#SBATCH --mail-type=BEGIN,END,FAIL

#SBATCH -o "%x.%j.o"

#SBATCH -e "%x.%j.e"


## squeue -u xuewen
## squeue -j your_job_id

scriptdir=/eva/home/xuewen/snpADanalysis/scripts

#working dir
cd /eva/home/xuewen/snpADanalysis

#mkdir realData
cd realData

# which snpAD
# ~/bin/snpAD

#use real data :
#snpAD 0.3.11

  #real data:
  #convert cram to bam   #   
   refh38=/eva/edatums/reference_materials/reference_genomes/hg38/GRCh38_full_analysis_set_plus_decoy_hla.fa

   sourcedir=/eva/projects/research/genotyping_assessment/real_data
   cram="SC_GMWOF5428715.alt_bwamem_GRCh38DH.20151208.WOLOFF.gambian_lowcov.cram"
   #$done cp $sourcedir/$cram ./

   

   
   #works: no chr and shortened ID, not used now. works with chr
   works="refh38nochr="GRCh38_full_analysis_set_plus_decoy_hla_nochr.fa"
   sed -E 's/chr//' $refh38 | sed -E 's/ .*//' >$refh38nochr
   refh38=$refh38nochr"
   
   #shortened ID
   refh38fixID="GRCh38_full_analysis_set_plus_decoy_hla_fixID.fa"
   cat $refh38 | sed -E 's/ .*//' >$refh38fixID
   refh38=$refh38fixID
   
   #cram to bam
   samtools view -@ 40 -T $refh38 -b $cram > ${cram/.cram/}.bam
   bamfile=${cram/.cram/}.bam
   samtools index -b  -@ 40 $bamfile 

   #get chr22 bam
   bam22=${cram/.cram/}.22.bam
   samtools view -@ 40 -b $bamfile chr22 >$bam22
   samtools index -b  -@ 40 $bam22
   
   #get chr2123 bam for multiple chromosomes
    bam2123=${cram/.cram/}.2123.bam
    samtools view -@ 40 -b $bamfile chr21 chr22 >$bam2123
    samtools index -b  -@ 40 $bam2123
   
 
   
   ## with chr and fix ID in bam files including head
   #bamInput="SC_GMWOF5428715.alt_bwamem_GRCh38DH.20151208.WOLOFF.gambian_lowcov.22.bam"
   bamInput=$bam2123
   bamOut=${bam2123/.bam/}.fixed.bam # .2123.fixed.bam
   samtools view -H $bamInput >tmp.chr.22.bam.head
   awk 'BEGIN {OFS="\t"} /^@SQ/ {print $1,$2,$3; next} {print $0}' tmp.chr.22.bam.head  >tmp.chr.22.bam.head.fixed
   # @SQ     SN:chr3 LN:198295559 
   samtools view $bamInput >tmp.chr.22.nohead.sam
   cat tmp.chr.22.bam.head.fixed tmp.chr.22.nohead.sam >tmp.$bamOut
   samtools view tmp.$bamOut -b >$bamOut
   #rm chr.22.fixed.sam
   samtools index -b $bamOut
   rm tmp.*
   
   bam22="chr.22.fixed.bam"
   
   
   #change head line of bam file: works for nochr, changed bam headline, bam alignment lines
   works="awk 'BEGIN {OFS="\t"} /^@SQ/ {print $1,$2,$3; next} {print $0}' nochr.22.bam.head  >nochr.22.bam.head.fixed
   samtools view nochr.22.bam >nochr.22.nohead.sam
   cat nochr.22.bam.head.fixed nochr.22.nohead.sam >nochr.22.fixed.sam
   samtools view nochr.22.fixed.sam -b >nochr.22.fixed.bam
   samtools index -b nochr.22.fixed.bam
   bam22="nochr.22.fixed.bam"
   "
   
   
# (3) Generate input file f:
Usage=' Bam2snpAD [-bAP?] [-i|--bam_index=FILE] [-r|--region=CHR:START-END]
        [-F|--region_file=STRING] [-b|--region_file_is_bedformat] [-A|--region_file_is_array]
        [-R|--readgroups=STRING] [-f|--fasta=FASTA] [-Q|--map_qual=MQ] [-q|--min_qual=Q]
        [--max_qual=Q] [-s|--size=N] [-o|--offset=N] [-P|--nopaired] [--merge_ends] [--version]
        [-?|--help] [--usage] [OPTIONS]* <bamfile>'
 snpADFile="chr22.snpAD"
 snpADFile="chr2122.snpAD"
 
 
 
 ##prepare bed file from vcf
 vcfFile="/eva/edatums/reference_materials/imputation/GnomadPanelV2/Panel/hgdp1kgp_chr22.filtered.SNV_INDEL.phased.shapeit5.unrelateds.nopp.vcf.gz"
 bedfile="chr22_nochr.bed" 
 bcftools view -HG $vcfFile |cut -f1,2 | awk 'BEGIN{ OFS="\t"} {print $1, $2-1, $2}'|sed -E 's/chr//'>$bedfile
 
 bedfile="chr22_chr.bed" 
 bcftools view -HG $vcfFile |cut -f1,2 | awk 'BEGIN{ OFS="\t"} {print $1, $2-1, $2}' >$bedfile
 
 #chr22,chr23
 regionsFile="regions.txt"
 printf "chr21\t1\t46709983\nchr22\t1\t50818468\n" >$regionsFile
 
 Bam2snpAD -F $regionsFile -f $refh38 -Q 25 -i $bamOut.bai $bamOut >$snpADFile
 
 #Bam2snpAD -r chr22 -f $refh38 -Q 25 -b $bedfile -i $bam22.bai $bam22 >$snpADFile # with -b not working
 # Bam2snpAD -r chr22 -f $refh38 -Q 25 -i $bam22.bai $bam22 >$snpADFile  #works with chr, and without bedfile, and fixed headline
 # Bam2snpAD -r 22 -f $refh38 -Q 25 -i $bam22.bai $bam22 >$snpADFile  #works without chr, and without bedfile, and fixed headline, fixed bam 
 # Bam2snpAD -r chr22 -f $refh38 -Q 25 -b $bedfile -i $bamfile.bai $bamfile >$snpADFile
 #-o 31 (giving error matrix numbers starting from 31)
 
# (6) Run parameter estimation (specify cores to use with -c; be aware that this steps can take very long):
 snpAD -o priors.txt -c 40 -O errors.txt $snpADFile >snpAD.log
 
# (7) Produce VCF:
  name="gambian_lowcov"
  outVcf="chr21chr22.snpAD.vcf"
 snpADCall -N $name -e errors.txt -p priors.txt $snpADFile > $outVcf
 

##view vcf

 #bcftools view output.VCF |less



  ##test new compiled : Bam2snpAD_AW snpADjoin_AW: works
  ~/Bam2snpAD_AW -r chr22 -f /eva/edatums/reference_materials/reference_genomes/hg38/GRCh38_full_analysis_set_plus_decoy_hla.fa -Q 25 -i SC_GMWOF5428715.alt_bwamem_GRCh38DH.20151208.WOLOFF.gambian_lowcov.bam.bai SC_GMWOF5428715.alt_bwamem_GRCh38DH.20151208.WOLOFF.gambian_lowcov.bam >snpADFile.Bam2snpAD_AW
  snpAD -o priors.txt -c 40 -O errors.txt snpADFile.Bam2snpAD_AW >snpAD_AW.log
   name="gambian_lowcov"
  outVcf="chr22.snpAD_AW.vcf"
  snpADCall -N $name -e errors.txt -p priors.txt snpADFile.Bam2snpAD_AW > $outVcf
  
  