#!/usr/bin/python
import sys
import getopt
import os
import subprocess

'''
Version: 1.0, June - 1 - 2024
Author: Xuewen Wang
Email: xwang.kib @ gmail.com
'''
__author__="Xuewen Wang"

script_path = __file__
directory_name = os.path.dirname(script_path)
script_name = os.path.basename(script_path)
separator = os.sep

def print_help():
    help_message =f"""
    Usage: python3 {script_name} [options]
    Options:
    -h, --help                  Show this help message and exit
    -a, --alignment=FILE        Input BAM file for alignment, the bam index file .bai should be present
    -t, --thread=NUMBER         Number of threads to use,default [12]
    -f, --fastaseq=FILE        fasta formatted sequence file
    -b, --bed=FILE              a tabular 1-based position file for target sites: \n\t\t\t\tchrID startPos ref alternative. e.g. chr21   5034230 C       T
    -o, --outprefix=PREFIX      Output prefix for result files 
    -c, --chromosomeFile=FILE   Chromosome list, tab separated text file, one chr per line, e.g. chr1
    -q, --quality=INT           The quality of mapped reads,default [25]
    -m, --minqual=INT           the minimal quality score of the base for analysis [10]
    -n, --normal                T|F. T will normalize the most likely PL to 0, default [F]
    version 1.0
    """
    print(help_message)


#parse arguments
def xpar(argv):
    bed  = ""
    chrFile="chrListAll22.txt"
    prefix="out"
    bam="align.bam"
    refGenome="hg38.fasta"
    threadn="12"
    quality="25"
    normal="F"
    basequal="10"

    if not argv:
        print_help()
        sys.exit(2)
    #e.g.
    # d=/eva/home/xuewen/snpADanalysis
    # python $d/scripts/$0.py -c $d/realData/chrList.txt -a $d/realData/SC_GMWOF5428715.alt_bwamem_GRCh38DH.20151208.WOLOFF.gambian_lowcov.2122.bam -t 40 -f /eva/edatums/reference_materials/reference_genomes/hg38/GRCh38_full_analysis_set_plus_decoy_hla.fa -o $d/realData/testPy -b hgdp1kgp_autos.filtered.SNV_INDEL.phased.shapeit5.unrelateds.nopp.biallelicsnps.tsv

    try:
        opts, args = getopt.getopt(argv, "ha:t:f:b:o:c:q:n:m:", ["alignment=", "thread=","fastaseq=",
                                                          "bed=","outprefix=","chromosomeFile=","quality=","normal=","minqual="])
    except getopt.GetoptError:
        print(f"python {script_name} [options] \n version 1.0 \n for help: python3 {script_name} -h")
        sys.exit(1)
    for opt, arg in opts:
        if opt == '-h':
            print_help()
            sys.exit()
        elif opt in ("-a", "align.bam"):
            bam = arg
        elif opt in ("-t", "--thread"):
            threadn = arg
        elif opt in ("-f", "--fastaseq"):
            refGenome = arg
        elif opt in ("-b","--bed"):
            bed=arg
        elif opt in ("-o", "--outprefix"):
            prefix=arg
        elif opt in ("-c", "--chromosomeFile"):
            chrFile = arg
        elif opt in ("-q", "--quality"):
            quality = arg
        elif opt in ("-n", "--normal"):
            normal= arg
        elif opt in ("-m", "--minqual"):
            basequal=arg
        else:
            print("unknown options")
            print_help()
    return bam, threadn, refGenome, bed, prefix,chrFile,quality, normal, basequal

def parseSnpVcf(inVcf):
    ##get genotype and PL for vcf sites
    vcfDict={}
    with open(inVcf, 'rt') as INVCF:
        # store var likelihood to dict
        # chr21   5030315 .       T       .       28      .       .       GT:DP:A:C:G:T:PP:GQ
        # 0/0:1:0,0:0,0:0,0:0,1:45,46,46,0,77,73,37,79,31,36:31
        # cells=INVCF.read() # whole
        for line in INVCF:
            if not line.strip():
                continue  # skip empty lines

            line = line.strip()
            if (line.startswith("#")):
                continue
            cells = line.split('\t')
            vcfKey = cells[0] + "\t" + cells[1] #chr21   5030315
            vcfLikehood = cells[9].split(':')[-2] #45,46,46,0,77,73,37,79,31,36
            vcfDict.update({vcfKey:vcfLikehood})
    return vcfDict

def getSnpGenotypeLikelihood(likehood):
    aa,cc,gg,tt,ac,ag,at,cg,ct,gt=likehood.split(',')
    likehoodDict={"AA":aa, "CC":cc, "GG":gg, "TT":tt,"AC":ac,"AG":ag, "AT":at,"CG":cg, "CT":ct, "GT":gt,"CA":ac,"GA":ag, "TA":at,"GC":cg, "TC":ct, "TG":gt}
    return likehoodDict

def getNormlikelihoods(likelihoodsList):
    minLh=min(likelihoodsList)
    normlikelihoodsList=[]
    if minLh >0:
        for L in likelihoodsList:
            normlikelihoodsList.append(L-minLh)
        return normlikelihoodsList
    else:
        return likelihoodsList

def creatNewVcfHeader(snpADvcfFileName: str):
    # read in orginal vcf from snpAD output and then add new information to create a new header
    with open(snpADvcfFileName) as svcf:
        vcfHeadlines=[]
        for line in svcf:
            line=line.strip()
            if line.startswith("#"):
                if "FORMAT=<ID=PP" not in line:
                    vcfHeadlines.append(line)
            else:
                break
        #add new infor
        anno='##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred-scaled likelyhood for genotypes hom1, het, and hom2">'
        vcfHeadlines.insert(-1,anno)
        headInfo='\n'.join(vcfHeadlines)

        return  headInfo





if __name__=='__main__':
    #get options
    bam, threadn, refGenome, bed, prefix,chrFile,quality,normal,basequal=xpar(sys.argv[1:])
    cmd=directory_name+separator+"Bam2snpAD"
    Bam2snpADopts_run=[cmd, "-F", chrFile ,"-f", refGenome,  "-Q", quality, "--min_qual", basequal, "-i", bam+".bai", bam]
    print("run: ", " ".join(Bam2snpADopts_run))
    outfileName=prefix+".snpADFile"
    with open(outfileName, 'wt') as of:
        process1=subprocess.run(Bam2snpADopts_run, check=True, stdout=of,  universal_newlines=True) #stdout=subprocess.PIPE
        print(" Step1 done")

    #snpAD -o priors.txt -c 40 -O errors.txt $snpADFile >snpAD.log
    cmd = directory_name + separator + "snpAD"
    snpADopts_run=[cmd, "-o", prefix+".priors.txt", "-c", threadn, "-O", prefix+".errors.txt", prefix+".snpADFile"]
    print("run: ", " ".join(snpADopts_run),"\n")
    process2 = subprocess.run(snpADopts_run, check=True, universal_newlines=True)
    print(" Step2 done")

    #snpADCall -N $name -e errors.txt -p priors.txt $snpADFile > $outVcf
    name="XX"
    outVcf=prefix+".vcf"
    cmd=directory_name+separator+"snpADCall"
    snpADCallopts_run=[cmd,"-N", name, "-e", prefix+".errors.txt", "-p", prefix+".priors.txt", prefix+".snpADFile"]
    print("run: ", " ".join(snpADCallopts_run))
    with open(outVcf, 'wt') as of2:
       process3 = subprocess.run(snpADCallopts_run, check=True, stdout=of2,stderr=subprocess.STDOUT, universal_newlines=True)
       print(" Step3 done")

    #print(f"done snpAD")

    #prepare new header of vcf
    #outVcf="C:/snpADtestData/chr21chr22.snpADtestSmall.vcf"
    #outVcf = "C:\snpADtestData\chr21chr22.snpAD.new.vcf"
    #outVcf ="/eva/home/xuewen/snpADanalysis/realData/chr21chr22.snpAD.new.vcf"
    newHeader=creatNewVcfHeader(outVcf)
    #print(newHeader)


    ##add genotype likelyhood to tsv sites
    print("run: to get genotype likelihoods for target sites in a tsv file")
    insnpADVcf =outVcf

    vcfDict = parseSnpVcf(insnpADVcf)
    #targetSites = "C:\snpADtestData\hgdp1kgp_autos.chr2122.testSmall.tsv"
    #targetSites="/eva/home/xuewen/snpADanalysis/realData/hgdp1kgp_autos.filtered.SNV_INDEL.phased.shapeit5.unrelateds.nopp.biallelicsnps.chr2122.tsv"
    #targetSites = "C:\snpADtestData\hgdp1kgp_autos.filtered.SNV_INDEL.phased.shapeit5.unrelateds.nopp.biallelicsnps.chr2122.tsv"
    targetSites =bed

    
    PLmaxi=1000000
    outAlleleProbResult = prefix+".prob.vcf"
    VCFPROB=open(outAlleleProbResult,"wt")
    VCFPROB.write(newHeader+"\n")
    with open(targetSites,'rt') as TS:
        for site in TS:
            #chr1	14487	G	A
            site=site.strip()
            scells=site.split('\t')
            skey=scells[0]+"\t"+scells[1]
            sgenotypes=[scells[2]+scells[2], scells[2]+scells[3],scells[3]+scells[3]]
            if skey in vcfDict:
                allLikelyhoods=vcfDict[skey] ##45,46,46,0,77,73,37,79,31,36
                glhoods=[]
                for gtype in sgenotypes:
                    glhood=getSnpGenotypeLikelihood(allLikelyhoods)[gtype]
                    if int(glhood) >= PLmaxi:
                        glhood=PLmaxi
                    glhoods.append(glhood)

                # normal PL
                if normal=="T":
                    glhoods=getNormlikelihoods(glhoods).copy()

                plvalues = ",".join(glhoods)
                #add id to 3rd column for scell1to5
                #col6-10:
                col6to10=['.','.','.','GT:PL','0/0:'+plvalues]
                scells.insert(2, '.')
                scells.extend(col6to10) #oneVcfLine
                #print("\t".join(scells))
                VCFPROB.write("\t".join(scells)+"\n")
            else:
                col6to10 = ['.', '.', '.', 'GT:PL', '0/0:0,0,0']
                scells.insert(2, '.')
                scells.extend(col6to10)
                #print("\t".join(scells))
                VCFPROB.write("\t".join(scells)+"\n")

    print(f"VCF file with allele probability: {outAlleleProbResult}")
    print(" All done")



    # to do: Pl int range(done), norm 10,20,30 >>>0,10,20 (done); min genotype Quality (), vcf output(done)

    #vcf output
    # CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  gambian_lowcov
   # chr21    5030315 .T.     28..GT: DP:A: C:G: T:PP: GQ     0 / 0: 1:0, 0: 0, 0: 0, 0: 0, 1: 45, 46, 46, 0, 77, 73, 37, 79, 31, 36: 31
VCFPROB.close()














