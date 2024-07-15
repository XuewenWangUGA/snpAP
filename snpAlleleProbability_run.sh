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
 
   d=/eva/home/xuewen/snpADanalysis
   
   python $d/scripts/snpAlleleProbability.py \
   -c $d/realData/chrListAll22.txt \
   -a $d/realData/SC_GMWOF5428715.alt_bwamem_GRCh38DH.20151208.WOLOFF.gambian_lowcov.2122.fixed.bam \
   -t 40 \
   -f $d/realData/GRCh38_full_analysis_set_plus_decoy_hla_fixID.fa \
   -o $d/realData/testPyProb \
   -b $d/realData/hgdp1kgp_autos.filtered.SNV_INDEL.phased.shapeit5.unrelateds.nopp.biallelicsnps.chr2122.tsv \
   -m 10
 
  