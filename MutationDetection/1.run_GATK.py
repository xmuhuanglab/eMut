import pysam
import collections
from multiprocessing import Pool
import os
import subprocess as sb
import argparse
import pybedtools

parser = argparse.ArgumentParser()
parser.add_argument('-b', '--bam', required=True, help='bam file to call variants')
parser.add_argument('-o', '--outputDir', required=True, help='Directory for output files')
parser.add_argument('-p', '--peaks', required=True, help='BED file for GATK to call variants')
parser.add_argument('-g', '--genome', required=True, help='chromsize (BED) file')
args = parser.parse_args()


file_bam = args.bam
out_bam_dir = args.outputDir + "/bam"
out_vcf_dir = args.outputDir + "/vcf"

def read_peak(bed):
    with open(bed, 'r') as infile:
        file_peaks = infile.readlines()
    file_peaks = ["\t".join(i.split("\t")[0:3]) + "\n" for i in file_peaks]
    file_peak_tmp = 'peak_tmp'
    with open(file_peak_tmp, 'w') as outfile:
        outfile.writelines(file_peaks)
    a = pybedtools.BedTool('peak_tmp').intersect(args.genome).moveto(bed.split("/")[-1] + ".bed")
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
    bed = out_vcf_dir + "/" + args.peaks.split("/")[-1] + ".bed"
    genome_ref = "./reference/hg38.fa"
    output_vcf = "singlecellvcf/" + bam.split("/")[-1].rstrip(".bam") + ".vcf"
    output_flt_vcf = "single_filter_vcf/" + bam.split("/")[-1].rstrip(".bam") + ".filter.vcf"
    sb.run('gatk Mutect2 --java-options "-Xmx16g" -L {} -R {} -I {} -O {}'.format(bed, genome_ref, bam, output_vcf), shell=True)
    sb.run('gatk FilterMutectCalls -R {} -V {} -O {}'.format(genome_ref, output_vcf, output_flt_vcf), shell=True)

#####
## Index for inout bam file
if not os.path.exists(file_bam + '.bai'):
    pysam.index('-@', '8', file_bam)

## reading reads in cells
seq = collections.defaultdict(list)
input_bam = pysam.AlignmentFile(file_bam, 'rb', threads = 16)
for i in input_bam:
    cell = i.query_name.split("_")[1]
    seq[cell].append(i)

## split bam per cell
if not os.path.exists(out_bam_dir): os.makedirs(out_bam_dir)
os.chdir(out_bam_dir)

## read header
input_bam = pysam.AlignmentFile(file_bam, 'rb')
header = input_bam.header.to_dict()

## write bam per cell
pool = Pool(5)
pool.imap(write_bam, seq.keys())
pool.close()
pool.join()

if not os.path.exists(out_vcf_dir): os.makedirs(out_vcf_dir)
os.chdir(out_vcf_dir)
## GATK limit
read_peak(args.peaks)
## index single bam
try:
    bam_list = [out_bam_dir+"/"+cell+".bam" for cell in seq.keys()]
except NameError:
    bam_list = os.listdir(out_bam_dir)
    bam_list = [out_bam_dir+"/"+i for i in bam_list if i.endswith("bam")]
pool = Pool(5)
pool.imap(pysam.index, bam_list)
pool.close()
pool.join()

## run GATK
print("run gatk" + file_bam)
pool = Pool(6)
pool.imap(run_gatk, bam_list)
pool.close()
pool.join()