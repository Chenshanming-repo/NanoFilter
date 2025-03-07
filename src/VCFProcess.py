import gzip
import os.path
from typing import List, Dict
from VariantRecord import *
"""
    输入vcf绝对路径
    输出VariantRecord数组
"""
def read_clair3_vcf(vcf_fn : str) -> Dict[str, List[VariantRecord]]:
    variant_records_dict = dict[str, List[VariantRecord]]()

    INDEL_num = INDEL_01_num = SNP_num = SNP_01_num = 0

    if not os.path.exists(vcf_fn):
        # logger.error(f"Absolute path '{vcf_fn}' doesn't exists.")
        return variant_records_dict

    if vcf_fn.endswith(".gz"):
        vcf_file = gzip.open(vcf_fn, "rt")
    else:
        vcf_file = open(vcf_fn, "r")
    vcf_record_no = 1

    for line in vcf_file.readlines():
        if type(line) == bytes:
            line = line.decode("utf-8")
        if line.startswith("#"):
            continue
        # CHROM  POS ID  REF ALT QUAL    FILTER  INFO    FORMAT  HG002
        variant_record = parse_clair3_vcf_line(line, vcf_record_no)

        cur_ctg = variant_record.ctg
        gt = variant_record.genotype


        if not cur_ctg in variant_records_dict:
            variant_records_dict[cur_ctg] = []

        variant_records_dict[cur_ctg].append(variant_record)

        if variant_record.is_INDEL():
            INDEL_num += 1
            if gt == "0/1":
                INDEL_01_num += 1
        else:
            SNP_num += 1
            if gt == "0/1":
                SNP_01_num += 1
        vcf_record_no += 1

    # logger.info(f"SNP: {SNP_num} (SNP 0/1: {SNP_01_num}), INDEL: {INDEL_num} (INDEL 0/1: {INDEL_01_num})")

    return variant_records_dict

"""
    输入line为clair3默认输出的vcf中variant记录的一行
        line_no为有效记录的行号
    输出该行对应的
"""
def parse_clair3_vcf_line(line : str, line_no : int) -> VariantRecord:
    if line.startswith("#"):
        # slogger.debug(f"Not a variant record line which starts with #.")
        return VariantRecord(DEFAULT_STRING, DEFAULT_INTEGER, DEFAULT_STRING, DEFAULT_STRING, DEFAULT_INTEGER, DEFAULT_STRING, VariantType.VARIANT_UNKNOWN, DEFAULT_INTEGER)

    ctg, pos, _, REF, ALT, QUAL, FILTER, INFO, FORMAT, SAMPLE = line.split("\t")
    pos = int(pos)
    gt = SAMPLE.split(':')[0]

    if gt == "./.":
        # logger.debug(f"{line}")
        gt = "1/1"

    # task1 calculate alt1/ref and alt2 from REF and ALT
    # task2 variant type (SNP, INDEL)
    allele1 = allele2 = ""
    variant_type = VariantType.VARIANT_UNKNOWN

    # only consider 0/1 case, 0 must in gt
    if "0" in gt and "1" in gt:
        allele1 = REF
        allele2 = ALT
        variant_type = VariantType.VARIANT_SNP if len(allele1) == len(allele2) else VariantType.VARIANT_INDEL
        # if variant_type == VariantType.VARIANT_INDEL:
            # pos += 1

    return VariantRecord(ctg, pos, allele1, allele2, int(float(QUAL)), gt, variant_type, line_no)

