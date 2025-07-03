[![DOI](https://zenodo.org/badge/892613239.svg)](https://doi.org/10.5281/zenodo.15354949)
# NanoFilter

## Overview

This software filters low-consistency variants from a VCF file using a reference genome and BAM alignment. It provides options to filter both SNPs and INDELs based on consistency scores.

## Usage

```bash
python main.py --reference <path_to_reference> \
                 --bam <path_to_bam> \
                 --vcf <path_to_vcf> \
                 --out-dir <output_directory> \
                 [--filter-all | --filter-indel] \
                 [-t <num_threads>] \
                 [-p1 <filter_param1>] \
                 [-p2 <filter_param2>] \
                 [-p3 <filter_param3>] \
                 [-q <quality_threshold>]
```

## Arguments

### Required Arguments

- `--reference <path>`: Path to the reference genome.
- `--bam <path>`: Path to the BAM file containing aligned reads.
- `--vcf <path>`: Path to the raw VCF file.
- `--out-dir <path>`: Path to the output directory.

### Optional Arguments

- `-t, --threads <int>`: Number of parallel threads (default: `10`).
- `--filter-all`: Filter both SNPs and INDELs based on consistency.
- `--filter-indel`: Filter only INDELs while keeping SNPs unchanged.
- `--filter-parameter1, -p1 <float>`: Filtering threshold for step 1 (default: `0.7`).
- `--filter-parameter2, -p2 <float>`: Filtering threshold for step 2 (default: `0.8`).
- `--filter-parameter3, -p3 <float>`: Filtering threshold for step 3 (default: `0.9`).
- `--quality, -q <float>`: Quality filtering threshold for low-coverage regions (default: `15`).

## Output

The script outputs a filtered VCF file in the specified output directory.

## Example

```bash
python src/main.py --reference ref.fasta \
                 --bam sample.bam \
                 --vcf raw_variants.vcf \
                 --out-dir filtered_output \
                 --filter-all \
                 -t 20 \
                 -p1 0.75 -p2 0.85 -p3 0.95 \
                 -q 20
```

## Installation & Dependencies

The script can be run like this:

### 1. Run from Source Code

Dependencies:

- Python 3.x
- `HapCUT2`
- Libraries specified in `environment.yml`

#### Installation Steps:

```bash
git clone <repository_url>
cd <repository_name>
conda env create -f environment.yml
conda activate <environment_name>

# install HapCUT2
git clone https://github.com/vibansal/HapCUT2.git
cd HapCUT2
make
# add HAPCUT2 and extractHAIRS to enviroment path
export PATH=$PATH:{PATH_TO_HAPCUT2}/build

python src/main.py --reference ref.fasta --bam sample.bam --vcf raw_variants.vcf --out-dir filtered_output
```



## License

This script is open-source and available under the MIT License.

## Contact

For issues or questions, please open an issue on the project repository or contact the developer.


