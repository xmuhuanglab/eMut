# %%
import pysam
import collections
from multiprocessing import Pool
import os
import subprocess as sb
import argparse
import pybedtools
from pathlib import Path
from glob import glob
import pandas as pd
script_dir = os.path.dirname(os.path.abspath(__file__))
# %%
parser = argparse.ArgumentParser()
parser.add_argument('-b', '--bam', required=True, help='bam file to call variants')
parser.add_argument('-o', '--outputDir', required=True, help='Directory for output files')
parser.add_argument('-p', '--peaks', required=True, help='BED file for GATK to call variants')
parser.add_argument('--chromsize', required=False, default=f"{script_dir}/Datas/hg38.chrom.sizes.bed", help=f'chromsize (BED) file, default="{script_dir}/Datas/hg38.chrom.sizes.bed"')
parser.add_argument('--reference', required=True, help='fasta file of reference genome')
parser.add_argument('--vepcache', required=True, help='Directory of vep cache')
args = parser.parse_args()


# %%
# file_bam = "/cluster/huanglab/mzhu/projects/aml/summary/eMut/eMut-main/test_run/demo.bam"
# peaks = "/cluster/huanglab/mzhu/projects/aml/summary/eMut/eMut-main/test_run/peaks.bed"
# outputDir = "/cluster/huanglab/mzhu/projects/aml/summary/eMut/eMut-main/test_run/GATK"
# chromsize =  "/cluster/huanglab/mzhu/projects/aml/summary/eMut/eMut-main/eMut/Datas/hg38.chrom.sizes.bed"
# reference = "/cluster2/huanglab/jzhu/ReferenceGenome/GATK_GRCh38/Homo_sapiens_assembly38.fasta"
# out_bam_dir = outputDir + "/bam"
# out_vcf_dir = outputDir + "/vcf"
# mafDir = outputDir + "/maf"
# vepDir = outputDir + "/vep"

# VEP_reference = "/cluster/huanglab/mzhu/reference/VEP_hg38"
# VEP_fasta = "/cluster2/huanglab/jzhu/ReferenceGenome/VEP_GRCh38/homo_sapiens/102_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa"
# %%
file_bam = os.path.abspath(args.bam)
outputDir = os.path.abspath(args.outputDir)
out_bam_dir = outputDir + "/bam"
out_vcf_dir = outputDir + "/vcf"
mafDir = outputDir + "/maf"
vepDir = outputDir + "/vep"
peaks = os.path.abspath(args.peaks)
chromsize = os.path.abspath(args.chromsize)
reference = os.path.abspath(args.reference)
VEP_cache = os.path.abspath(args.vepcache)

# %%
def read_peak(bed):
    with open(bed, 'r') as infile:
        file_peaks = infile.readlines()
    file_peaks = ["\t".join(i.split("\t")[0:3]) + "\n" for i in file_peaks]
    file_peak_tmp = 'peak_tmp'
    with open(file_peak_tmp, 'w') as outfile:
        outfile.writelines(file_peaks)
    a = pybedtools.BedTool('peak_tmp').intersect(chromsize).moveto(bed.split("/")[-1])
    del(a)
    os.remove('peak_tmp')

def write_bam(cell):
    output_tmp = cell + "_tmp"
    ## create header for GATK
    header['RG'] = [{'ID': cell, 'LB': 'library1', 'PL': 'illumina', 'SM': cell, 'PU': 'pl1'}]
    outbam = pysam.AlignmentFile(output_tmp, 'wb', header=header)
    for i in seq[cell]:
        ## create tag (RG:Z:cell) for GATK
        i.set_tag("RG", cell, "Z")
        outbam.write(i)
    outbam.close()
    pysam.sort('-@', '8', output_tmp, '-o', cell + ".bam")
    os.remove(output_tmp)

def run_gatk(bam):
    ## Only call SNPs in peaks
    if not os.path.exists("singlecellvcf"): os.mkdir("singlecellvcf")
    if not os.path.exists("single_filter_vcf"): os.mkdir("single_filter_vcf")
    bed = out_vcf_dir + "/" + peaks.split("/")[-1]
    output_vcf = "singlecellvcf/" + bam.split("/")[-1].rstrip(".bam") + ".vcf"
    output_flt_vcf = "single_filter_vcf/" + bam.split("/")[-1].rstrip(".bam") + ".filter.vcf"
    sb.run('gatk Mutect2 --java-options "-Xmx16g" -L {} -R {} -I {} -O {}'.format(bed, reference, bam, output_vcf), shell=True)
    sb.run('gatk FilterMutectCalls -R {} -V {} -O {}'.format(reference, output_vcf, output_flt_vcf), shell=True)

#####
# %%  Index for inout bam file
if not os.path.exists(file_bam + '.bai'):
    pysam.index('-@', '8', file_bam)

# %% reading cell barcodes in cells
seq = collections.defaultdict(list)
input_bam = pysam.AlignmentFile(file_bam, 'rb', threads = 16)
for i in input_bam:
    cell = i.query_name.split("_")[1]
    seq[cell].append(i)

# %% split bam per cell
if not os.path.exists(out_bam_dir): os.makedirs(out_bam_dir)
os.chdir(out_bam_dir)

# %% read header
input_bam = pysam.AlignmentFile(file_bam, 'rb')
header = input_bam.header.to_dict()

# %% write bam per cell
pool = Pool(5)
pool.imap(write_bam, seq.keys())
pool.close()
pool.join()

# %%
if not os.path.exists(out_vcf_dir): os.makedirs(out_vcf_dir)
os.chdir(out_vcf_dir)
# %% GATK limit
read_peak(peaks)
# %% index single bam
try:
    bam_list = [out_bam_dir+"/"+cell+".bam" for cell in seq.keys()]
except NameError:
    bam_list = os.listdir(out_bam_dir)
    bam_list = [out_bam_dir+"/"+i for i in bam_list if i.endswith("bam")]
# pool = Pool(5)
# pool.imap(pysam.index, bam_list)
# pool.close()
# pool.join()

# # %% run GATK
# print("run gatk" + file_bam)
# pool = Pool(6)
# pool.imap(run_gatk, bam_list)
# pool.close()
# pool.join()

# %%
if not os.path.exists(mafDir):
    os.makedirs(mafDir)

if not os.path.exists(vepDir):
    os.makedirs(vepDir)

# %%
cells = [Path(i).stem.split('.')[0] for i in glob(f"{out_vcf_dir}/single_filter_vcf/*.vcf")]
# %%
def mutation_annotation(cell):
    # select PASS mutations
    os.system(f"bcftools filter -i \'FILTER==\"PASS\"\' \
    -o {vepDir}/{cell}.filter.vcf {out_vcf_dir}/single_filter_vcf/{cell}.filter.vcf")
    #  mutation annotation
    os.system(f"vep --species homo_sapiens --assembly GRCh38 --offline --no_progress --no_stats --sift b --ccds --uniprot --hgvs --symbol \
    --numbers --domains --gene_phenotype --canonical --protein --biotype --tsl --pubmed --variant_class --shift_hgvs 1 --check_existing --total_length \
    --allele_number --no_escape --xref_refseq --failed 1 --vcf --minimal --flag_pick_allele --pick_order canonical,tsl,biotype,rank,ccds,length \
    --dir {VEP_cache} --cache_version 102 \
    --fasta {reference} \
    --input_file {vepDir}/{cell}.filter.vcf \
    --output_file {vepDir}/{cell}.vep.vcf \
    --polyphen b --af --af_1kg --regulatory --fork 10 --force_overwrite --buffer_size 10000")
    os.system(f"rm {vepDir}/{cell}.filter.vcf")
    # vcf to maf format
    os.system(f"vcf2maf.pl --input-vcf {vepDir}/{cell}.vep.vcf \
                --output-maf {mafDir}/{cell}.vep.maf --tumor-id {cell} \
                --vep-data {VEP_cache} \
                --ref-fasta {reference} \
                --ncbi-build GRCh38 --vep-overwrite --inhibit-vep ")
# %%
pool = Pool(12)
pool.imap(mutation_annotation, cells)
pool.close()
pool.join()
# %%
with open(f"{outputDir}/{Path(file_bam).stem}_cells.txt",'w') as f:
    f.writelines([Path(i).stem + '\n' for i in glob("/cluster/huanglab/mzhu/projects/aml/summary/eMut/eMut-main/test_run/GATK/bam/*.bam")])
# %%
maf_list = []
for i in glob(f"{mafDir}/*.maf"):
    maf_list.append(pd.read_table(i, header=1))
df = pd.concat(maf_list).to_csv(f"{outputDir}/{Path(file_bam).stem}.maf", sep='\t', index=False)
# %%
