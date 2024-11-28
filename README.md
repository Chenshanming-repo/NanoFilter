# NanoFilter
A tool for consistency filtering of SNP and INDEL sites to make them more suitable for genotyping and assembly.

## Singularity
```
singularity pull nanofilter.v1.1.sif docker://chenshanming/nanofilter:v1.1
```
Here is an example to run NanoFilter
```
#!/bin/bash 
REF="/public/home/hpc224712204/ref/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna"
BAM="/public/home/hpc224712204/data/hg002_tag/r10_duplex/all_30x.bam"
VCF="/public/home/hpc224712204/bench/hapdup_GRCH38_hg002/r10_duplex/raw_snp_filtered_indel_weight_th_qual15/pepper/PEPPER_VARIANT_FULL.vcf"
OUT_DIR="/public/home/hpc224712204/software/NanoFilter/test_out"
THREAD="48"
COMMAND="source /opt/conda/bin/activate NanoFilterEnv && python src/main.py --reference $REF --bam $BAM --vcf $VCF --out-dir $OUT_DIR -t $THREAD --filter-all"
singularity exec nanofilter.v1.1.sif bash -c "$COMMAND"
```
