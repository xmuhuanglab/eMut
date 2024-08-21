# %%
from pathlib import Path
import os
import argparse
script_dir = os.path.dirname(os.path.abspath(__file__))

parser = argparse.ArgumentParser()
parser.add_argument('-b', '--bam', required=True, help='Input BAM file for variant calling')
parser.add_argument('-o', '--outputDir', required=True, help='Directory to store output files')
parser.add_argument('--monopogen', required=True, default='path to monopogen directory')
parser.add_argument('--reference', required=True, help='Reference genome in FASTA format')
parser.add_argument('--region', required=False, help=f'region list contain which chromosomes to run, default="{script_dir}/Datas/region.lst", to run all chromosomes', default=f"{script_dir}/Datas/region.lst")
parser.add_argument('--variant', required=True, help='Variant panel file')
parser.add_argument('--cells', required=True, help='cell barcodes of cells to run')
args = parser.parse_args()

# %%
monopogen_path = Path(args.monopogen).absolute()
try:
    os.environ['LD_LIBRARY_PATH'] = str(monopogen_path /"apps:") + os.environ['LD_LIBRARY_PATH']
except:
    os.environ['LD_LIBRARY_PATH'] = str(monopogen_path /"apps")

bam = Path(args.bam).absolute()
outputdir = Path(args.outputDir).absolute()
variant_panel = Path(args.variant).absolute()
region = Path(args.region).absolute()
fasta = Path(args.reference).absolute()
cells = Path(args.cells).absolute()
# %%
# bam = Path("../test_run/demo.bam").absolute()
# outputdir = "/cluster/huanglab/mzhu/projects/aml/summary/eMut/eMut-main/test_run/Monopogen"
# variant_panel = "/cluster2/huanglab/jzhu/app/Monopogen/1KG3"
# fasta = "/cluster2/huanglab/jzhu/ReferenceGenome/GATK_GRCh38/Homo_sapiens_assembly38.fasta"

# %%
outputdir = Path(outputdir)
outputdir.mkdir(exist_ok=True, parents=True)
bam_list = outputdir/"bam_list.txt"
with open(bam_list, "w") as f:
    f.write(f"{bam.stem},{bam}")



# %%
os.system(f"{monopogen_path}/src/Monopogen.py preProcess -b {bam_list} \
  -o {outputdir/bam.stem}  -a {monopogen_path}/apps -t 8")
# # %%

os.system(f"{monopogen_path}/src/Monopogen.py germline \
     -a {monopogen_path}/apps -t 5 -r {region} \
     -p {variant_panel}/ \
     -g {fasta} -m 3 -s all \
     -o {outputdir/bam.stem} ")
# %%

os.system(f"python {monopogen_path}/src/Monopogen.py somatic \
    -a {monopogen_path}/apps -r {region} -t 6 \
    -i {outputdir/bam.stem} -l {cells} -s all \
    -g {fasta}")


#os.system(f"python {monopogen_path}/src/Monopogen.py  somatic  \
#     -a {monopogen_path}/apps -r {region} -t 3 \
#     -i {outputdir/bam.stem} -l {cells} -s cellScan \
#     -g {fasta}")

#os.system(f"python {monopogen_path}/src/Monopogen.py  somatic  \
#     -a {monopogen_path}/apps -r {region} -t 5 \
#     -i {outputdir/bam.stem} -l {cells} -s LDrefinement \
#     -g {fasta}")