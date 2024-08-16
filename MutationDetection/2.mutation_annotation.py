##############################################################################################################
#
#
#
#
##############################################################################################################
import os
from multiprocessing import Pool
import argparse

parser=argparse.ArgumentParser(description='mutation_annotation.py is a mulit processes tool for mutation annotation')
parser.add_argument('-input',type=str,help='Input your cells txt',required=True)
parser.add_argument('-mafDir',type=str,help='Input your output maf file directory',required=True)
parser.add_argument('-vepDir',type=str,help='Input your output vep directory',required=True)
parser.add_argument('-vcfDir',type=str,help='Input your input vcf directory',required=True)

args=parser.parse_args()

input_cells = args.input
mafDir = args.mafDir
vepDir = args.vepDir
vcfDir = args.vcfDir


if not os.path.exists(mafDir):
    os.makedirs(mafDir)

if not os.path.exists(vepDir):
    os.makedirs(vepDir)

with open(input_cells,'r') as cell_file:
    cell_n = cell_file.readlines()
    
cell_n = [i.rstrip() for i in cell_n]

def mutation_annotation(item):
    item = item.rstrip('\n')
    # select PASS mutations
    os.system("./app/bcftools-1.10/bcftools filter -i \'FILTER==\"PASS\"\' \
    -o {vepDir}/{cell}.filter.vcf {vcfDir}/{cell}.filter.vcf".format(vcfDir=vcfDir,vepDir=vepDir,cell=item))
    #  mutation annotation
    os.system("vep --species homo_sapiens --assembly GRCh38 --offline --no_progress --no_stats --sift b --ccds --uniprot --hgvs --symbol \
    --numbers --domains --gene_phenotype --canonical --protein --biotype --tsl --pubmed --variant_class --shift_hgvs 1 --check_existing --total_length \
    --allele_number --no_escape --xref_refseq --failed 1 --vcf --minimal --flag_pick_allele --pick_order canonical,tsl,biotype,rank,ccds,length \
    --dir ./ReferenceGenome/VEP_GRCh38/ \
    --fasta ./ReferenceGenome/VEP_GRCh38/homo_sapiens/102_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa  \
    --input_file {vepDir}/{cell}.filter.vcf \
    --output_file {vepDir}/{cell}.vep.vcf \
    --polyphen b --af --af_1kg --af_esp --regulatory --fork 10 --force_overwrite --buffer_size  10000".format(vepDir=vepDir,cell=item))
    os.system("rm {vepDir}/{cell}.filter.vcf".format(vepDir=vepDir,cell=item))
    # vcf to maf format
    os.system("perl ./app/vcf2maf-1.6.21/vcf2maf.pl --input-vcf {vepDir}/{cell}.vep.vcf \
                --output-maf {mafDir}/{cell}.vep.maf --tumor-id {cell} \
                --vep-path ./app/miniconda3/envs/VEP/bin/ \
                --vep-data ./ReferenceGenome/VEP_GRCh38/ \
                --ref-fasta ./ReferenceGenome/VEP_GRCh38/homo_sapiens/102_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa \
                --ncbi-build GRCh38 --vep-overwrite --inhibit-vep ".format(vepDir=vepDir,mafDir=mafDir,cell=item))

pool = Pool(12)
pool.imap(mutation_annotation, cell_n)
pool.close()
pool.join()