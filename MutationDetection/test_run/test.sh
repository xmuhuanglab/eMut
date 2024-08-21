python ../eMut/run_GATK.py -b demo.bam -p peaks.bed -o GATK \
    --reference ./ReferenceGenome/VEP_GRCh38/homo_sapiens/102_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa \
    --vepcache ./ReferenceGenome/VEP_GRCh38/

python ../eMut/run_monopogen.py -b demo.bam -o Monopogen \
    --variant /cluster2/huanglab/jzhu/app/Monopogen/1KG3 \
    --reference ./ReferenceGenome/GATK_GRCh38/Homo_sapiens_assembly38.fasta \
    --monopogen ./eMut/Monopogen \
    --cells ./test_run/demo_cells.txt \
    --region ./test_run/region.lst

