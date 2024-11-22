#!/bin/bash 

ROOT="/public/home/hpc224712204/software/NanoFilter/NanoFilter"
REF="/public/home/hpc224712204/ref/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna"
BAM="/public/home/hpc224712204/data/hg002_tag/r10_duplex/all_30x.bam"
VCF="/public/home/hpc224712204/bench/hapdup_GRCH38_hg002/r10_duplex/raw_snp_filtered_indel_weight_th_qual15/pepper/PEPPER_VARIANT_FULL.vcf"
OUT_DIR="/public/home/hpc224712204/software/NanoFilter/test_out"

python $ROOT/main.py --reference $REF --bam $BAM --vcf $VCF --out-dir $OUT_DIR -t 48 --filter-all
