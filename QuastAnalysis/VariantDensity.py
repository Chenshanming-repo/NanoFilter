from Utils import read_contig_info_from_phaseblock_name

def read_vcf_file(vcf_file):

    print(f"read {vcf_file}")

    res = {}

    def process_line(_line):
        if _line.startswith('#'):
            return

        fields = _line.strip().split('\t')
        if len(fields) < 9:  # Need at least 9 columns for genotype info
            return

        chrom = fields[0]
        pos = int(fields[1])
        ref = fields[3]
        alts = fields[4].split(',')
        format_fields = fields[8].split(':')
        sample_info = fields[9].split(':')

        if not chrom in res:
            res[chrom] = []
        # res[chrom].append()


        # Check for heterozygosity (GT field should be 0/1 or 1/0)
        gt_index = format_fields.index('GT') if 'GT' in format_fields else -1
        if gt_index == -1:
            return

        gt = sample_info[gt_index]
        if gt not in ('0/1', '1/0', '0|1', '1|0', '1/2', '1|2', '2/1', '2|1', '0/2', '2/0', '0|2', '2|0'):
            return

        # Determine variant type
        is_snp = True
        ref_len = len(ref)
        for alt in alts:
            if len(alt) != ref_len or ref_len != 1:
                is_snp = False
                break
        res[chrom].append({
            'chrom': chrom,
            'pos': pos,
            'is_snp': is_snp
        })






    try:
        line_no = 0
        if vcf_file.endswith('.gz'):
            import gzip
            with gzip.open(vcf_file, 'rt') as f:
                for line in f:
                    process_line(line)
                    line_no += 1
                    #if line_no % 100000 == 0:
                    #    print(line_no)
        else:
            with open(vcf_file, 'r') as f:
                for line in f:
                    process_line(line)
                    line_no += 1
                    #if line_no % 100000 == 0:
                    #    print(line_no)
    except Exception as e:
        raise ValueError(f"Error reading VCF file: {str(e)}")

    return res

def cal_variant_density_region(vcf_dict, contig, start, end):
    if not contig in vcf_dict:
        return {
        'snp_count': 0,
        'indel_count': 0,
        'total_variants': 0,
        'region_length': end - start + 1,
        'snp_density': 0,
        'indel_density': 0
    }
    variants = vcf_dict[contig]
    snp_count = 0
    indel_count = 0
    region_length = end - start + 1

    for variant in variants:
        if start <= variant['pos'] <= end:
            if variant['is_snp']:
                snp_count += 1
            else:
                indel_count += 1

    # Calculate densities (per kb)
    snp_density = (snp_count / region_length) * 1000 if region_length > 0 else 0
    indel_density = (indel_count / region_length) * 1000 if region_length > 0 else 0

    return {
        'snp_count': snp_count,
        'indel_count': indel_count,
        'total_variants': snp_count + indel_count,
        'region_length': region_length,
        'snp_density': snp_density,
        'indel_density': indel_density
    }

def output_contig_variant_density(vcf_dict):
    total_length = 0
    total_snp = 0
    total_indel = 0
    for contig in vcf_dict:
        contig_start = 1_000_000_000_000
        contig_end = -1
        for variant in vcf_dict[contig]:
            # print(variant["pos"])
            contig_start = min(variant["pos"], contig_start)
            contig_end = max(variant["pos"], contig_end)
        #print(contig_start, contig_end, contig)
        res = cal_variant_density_region(vcf_dict, contig, contig_start, contig_end)
        total_length += res['region_length']
        total_snp += res['snp_count']
        total_indel += res['indel_count']
    #print(total_snp, total_indel, total_length)
    print(f"Total:\nSNP Density: {(total_snp / total_length) * 1000}")
    print(f"Total:\nINDEL Density: {(total_indel / total_length) * 1000}")




def output_variant_density(compare_file, vcf_dict, output_file, window_size=100000):
    total_length = 0
    total_indel = 0
    total_snp = 0
    with open(output_file, 'w') as o:
        o.write("#CONTIG_NAME\tCONTIG_START\tCONTIG_END\tSNP_COUNT\tSNP_DENSITY\tINDEL_COUNT\tINDEL_DENSITY\n")
        with open(compare_file, 'r') as f:
            for line in f:
                if line.startswith("#"):
                    continue
                fields = line.strip().split("\t")
                contig_name = read_contig_info_from_phaseblock_name(fields[1])[1]
                #misassembly_start = int(fields[5]) - min(int(fields[8]), window_size)
                #misassembly_end = int(fields[6]) + min(int(fields[9]), window_size)
                misassembly_start = max(int(fields[5]) - window_size, 0)
                misassembly_end = int(fields[6]) + window_size

                res = cal_variant_density_region(vcf_dict, contig_name, misassembly_start, misassembly_end)
                o.write(f"{contig_name}\t{misassembly_start}\t{misassembly_end}\t{res['snp_count']}\t{res['snp_density']}\t{res['indel_count']}\t{res['indel_density']}\t\n")
                total_indel += res['indel_count']
                total_snp+= res['snp_count']
                total_length += res['region_length']
    print("snp density", (total_snp / total_length) * 1000, compare_file)
    # print("indel density", (total_indel / total_length) * 1000, compare_file)
