######
###   for example (AML10)
###################################  Monopogen  ###########################################################################################
########   1K3G download
chr=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" chr22")
url="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05"
for i in ${chr[@]}; do \
axel --num-connections=6 ${url}_${i}.filtered.shapeit2-duohmm-phased.vcf.gz;done 
for i in ${chr[@]}; do \
axel --num-connections=6 ${url}_chr2.filtered.shapeit2-duohmm-phased.vcf.gz.tbi

axel --num-connections=6 ${url}_chr2.filtered.shapeit2-duohmm-phased.vcf.gz

#########    reference genome download
axel --num-connections=6 https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/latest/hg38.fa.gz
axel --num-connections=6 https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/latest/md5sum.txt
md5sum hg38.fa.gz
gzip -d hg38.fa.gz
nohup bwa index hg38.fa -p ./hg38.fa 2>>hg38.fa_index.log &
samtools faidx hg38.fa


######     1) mutation detection
path="./app/Monopogen" # where Monopogen is downloaded
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${path}/apps

sam=AML10
if [ ! -d "./Monopogen/${sam}" ]; then
    mkdir ./Monopogen/${sam}
    echo "已创建文件夹"
fi

python  ${path}/src/Monopogen.py  preProcess -b ./Monopogen/bam.lst \
 -o ./Monopogen/${sam}  -a ${path}/apps -t 8

python  ${path}/src/Monopogen.py  germline  \
    -a   ${path}/apps -t 5   -r  ./Monopogen/region.lst \
    -p  ${path}/1KG3/ \
    -g  ./ReferenceGenome/UCSC/hg38.fa -m 3 -s all  -o ./Monopogen/${sam}


cd ./Monopogen/sampleCell/
python  ${path}/src/Monopogen.py  somatic  \
    -a   ${path}/apps  -r  ./Monopogen/region.lst  -t 6 \
    -i  ./Monopogen/${sam}  -l  ${sam}_readNum.csv   -s featureInfo     \
    -g ./ReferenceGenome/UCSC/hg38.fa

python  ${path}/src/Monopogen.py  somatic  \
    -a   ${path}/apps  -r  ./Monopogen/region.lst  -t 3 -w 50MB \
    -i  ./Monopogen/${sam}  -l ${sam}_readNum.csv   -s cellScan     \
    -g ./ReferenceGenome/UCSC/hg38.fa

python  ${path}/src/Monopogen.py  somatic  \
    -a   ${path}/apps  -r  ./Monopogen/region.lst  -t 5 -w 50MB \
    -i  ./Monopogen/${sam}  -l ${sam}_readNum.csv   -s LDrefinement     \
    -g ./ReferenceGenome/UCSC/hg38.fa;done


######      2) combine SNV for each sample
# eMut.ipynb


######      2) mutations annotation
sort -k1,1 -k2n ./Monopogen/summary/SNV.vcf > ./Monopogen/summary/SNV.sorted.vcf
vep --species homo_sapiens --assembly GRCh38 --offline --no_progress --no_stats --sift b --ccds --uniprot --hgvs --symbol \
    --numbers --domains --gene_phenotype --canonical --protein --biotype --tsl --pubmed --variant_class --shift_hgvs 1 --check_existing --total_length \
    --allele_number --no_escape --xref_refseq --failed 1 --vcf --minimal --flag_pick_allele --pick_order canonical,tsl,biotype,rank,ccds,length \
    --dir ./ReferenceGenome/VEP_GRCh38/ \
    --fasta ./ReferenceGenome/VEP_GRCh38/homo_sapiens/102_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa  \
    --input_file ./Monopogen/summary/SNV.sorted.vcf \
    --output_file ./Monopogen/summary/SNV.vep.vcf \
    --polyphen b --af --af_1kg --af_esp --regulatory --fork 10 --force_overwrite --buffer_size  10000

perl /cluster2/huanglab/jzhu/app/vcf2maf-1.6.21/vcf2maf.pl --input-vcf ./Monopogen/summary/SNV.vep.vcf \
                --output-maf ./Monopogen/summary/SNV.vep.maf \
                --vep-path ./app/miniconda3/envs/VEP/bin/ \
                --vep-data ./ReferenceGenome/VEP_GRCh38/ \
                --ref-fasta ./ReferenceGenome/VEP_GRCh38/homo_sapiens/102_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa \
                --ncbi-build GRCh38 --vep-overwrite --inhibit-vep 

#####################################  GATK Mutect2    ##############################################################################################
#####      1) mutation detection
python ./1.run_GATK.py \
    -b ./AML10_ATAC_dedup.bam \
    -o ./AML10 \
    -p ./peaks.bed \   #  peak location file to save time for mutation detection
    -g ./reference/UCSC/chrom_sizes/hg38.chrom.sizes_flt.bed

#####     2) mutation annotation(VEP)
cd ./AML10/vcf/single_filter_vcf/
for i in `ls *.filter.vcf`; do echo ${i%%.*} >> ./AML10/cells.txt; done
nohup python ./2.mutation_annotation.py \
    -input  ./AML10/cells.txt \
    -mafDir ./AML10/maf \
    -vepDir ./AML10/vep \
    -vcfDir ./AML10/vcf/single_filter_vcf > /dev/null &


