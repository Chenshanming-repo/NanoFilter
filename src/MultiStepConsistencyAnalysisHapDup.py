import gzip
import os.path
from VCFProcess import parse_clair3_vcf_line
import  numpy as np
from merge_fragment_file import read_fragment_file_list


# logger = logging.getLogger()

# def _enable_logging(log_file, debug, overwrite):
#     """
#     Turns on logging, sets debug levels and assigns a log file
#     """
#     log_formatter = logging.Formatter("[%(asctime)s] %(name)s: %(levelname)s: "
#                                       "%(message)s", "%Y-%m-%d %H:%M:%S")
#     console_formatter = logging.Formatter("[%(asctime)s] %(levelname)s: "
#                                           "%(message)s", "%Y-%m-%d %H:%M:%S")
#     console_log = logging.StreamHandler()
#     console_log.setFormatter(console_formatter)
#     if not debug:
#         console_log.setLevel(logging.INFO)
#
#     if overwrite:
#         open(log_file, "w").close()
#     file_handler = logging.FileHandler(log_file, mode="a")
#     file_handler.setFormatter(log_formatter)
#
#     logger.setLevel(logging.DEBUG)
#     logger.addHandler(console_log)
#     logger.addHandler(file_handler)

def read_pepper_vcf(vcf_fn):
    variant_records_list = []

    INDEL_num = INDEL_01_num = SNP_num = SNP_01_num = 0

    if not os.path.exists(vcf_fn):
        # logger.error(f"Absolute path '{vcf_fn}' doesn't exists.")
        return variant_records_list

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


        variant_records_list.append(variant_record)

        if variant_record.is_INDEL():
            INDEL_num += 1
            if gt == "0/1":
                INDEL_01_num += 1
        else:
            SNP_num += 1
            if gt == "0/1":
                SNP_01_num += 1
        vcf_record_no += 1

    print(f"SNP: {SNP_num} (SNP 0/1: {SNP_01_num}), INDEL: {INDEL_num} (INDEL 0/1: {INDEL_01_num})")

    return variant_records_list

def calculate_consistency_mp(vcf_fn, fragment_fn, output_dir):
    pass


def calculate_consistency(vcf_fn, fragment_fn, output_dir):
    # from hapdup.merge_fragment_file import read_fragment_file_list

    # logger.info("Using fragment list.")

    fragment_line_items = read_fragment_file_list(fragment_fn)
    clair3_vcf_records_list = read_pepper_vcf(vcf_fn)
    # clair3_vcf_records_dict = { vcf_record.key : vcf_record for vcf_record in clair3_vcf_records_list }

    consistency_file = open(f"{output_dir}/consistency.results", "w")
    log_file = open(f"{output_dir}/consistency.log", "w")
    log_file.write(f"#vcf_fn: {vcf_fn}.\n")
    log_file.write(f"#fragment_fn: {fragment_fn}.\n")
    log_file.write(f"#out_dir: {output_dir}.\n")


    consistency_dict = {}

    # from common.Timer import Timer
    # timer = Timer()
    n_processed_reads = 0

    fragment_line_num = len(open(fragment_fn).readlines())

    for line_items in fragment_line_items:
        n_processed_reads += 1
        if n_processed_reads % 5000 == 0:
            # logger.info(f"Processed {n_processed_reads} / {fragment_line_num}")
            print(f"Processed {n_processed_reads} / {fragment_line_num}")
        tags = line_items[0].split(" ")
        i = 0

        tags_dict = {}
        qual_dict =  {}

        while i < len(tags):
            block_start_idx = int(tags[i])
            block_tags = tags[i+1]
            j = 0
            while j < len(block_tags):
                tags_dict[f"{clair3_vcf_records_list[block_start_idx+j-1].key}"] = int(block_tags[j])
                qual_dict[f"{clair3_vcf_records_list[block_start_idx+j-1].key}"] = float(clair3_vcf_records_list[block_start_idx+j-1].qual)
                j += 1
            i += 2

        # if "chr2:92508951" in tags_dict or "chr2:92508950" in tags_dict:
        # print(tags_dict)

        for key1 in tags_dict:
            # 只有INDEL计算与SNP一致性 SNP直接过
            for key2 in tags_dict:
                if key1 == key2:
                    continue
                key1_ctg = key1.split(":")[0]
                if not key1_ctg in consistency_dict:
                    consistency_dict[key1_ctg] = {}
                    # logger.info(f"{key1_ctg} not in consistency_dict")
                    # print(f"{key1_ctg} not in consistency_dict")
                if not key1 in consistency_dict[key1_ctg]:
                    consistency_dict[key1_ctg][key1] = {}
                if not key2 in consistency_dict[key1_ctg][key1]:
                    consistency_dict[key1_ctg][key1][key2] = [[0, 0], [0, 0]]
                consistency_dict[key1_ctg][key1][key2][tags_dict[key1]][tags_dict[key2]] += 1

    single_var_num = 0

    log_file.write(f"# {single_var_num} single variant.\n")
    consistency_file.close()
    log_file.close()


    # with open(f"{output_dir}/saved_dictionary.pkl", "wb") as f:
        # pickle.dump(consistency_dict, f)

    return consistency_dict

def filter_by_consistency(consistency_dict_ctg, n_minimum_covered_reads, n_min_variant_position, consistency_threshold, vcf_fn, filtered_keys = None, filter_not_enough_position=False,filter_snp=True, use_avg_consistency=True, only_cal_with_snp=False):
    filtered_keys = set([]) if filtered_keys is None else set(filtered_keys)
    single_snp_num = 0
    single_indel_num = 0
    filtered_indel_num = 0
    filtered_snp_num = 0
    clair3_vcf_records_list = read_pepper_vcf(vcf_fn)
    clair3_vcf_records_dict = { vcf_record.key : vcf_record for vcf_record in clair3_vcf_records_list }

    # logger.info(f"n_minimum_covered_reads: {n_minimum_covered_reads}")
    # logger.info(f"n_min_variant_position: {n_min_variant_position}")
    # logger.info(f"consistency_threshold: {consistency_threshold}")
    # logger.info(f"filter_snp: {filter_snp}")
    # logger.info(f"filter_not_enough_position: {filter_not_enough_position}")
    # logger.info(f"use_avg_consistency: {use_avg_consistency}")

    # with open(consistency_dict_fn, "rb") as f:
    # consistency_dict = pickle.load(f)

    # logger.info("Loaded consistency dictionary.")
    filter_by_qual_num = 0
    consistency_info_list = []
    for consistency_dict_ctg_key in consistency_dict_ctg:
        consistency_dict = consistency_dict_ctg[consistency_dict_ctg_key]
        # consistency_info_list = []
        for key1 in consistency_dict:
            if key1 in filtered_keys:
                continue
            # if not filter_snp and not clair3_vcf_records_dict[key1].is_INDEL():
                #continue
            weight_list = []
            consistency_list = []
            linked_sites = []
            for key2 in consistency_dict[key1]:
                if key2 in filtered_keys:
                    continue
                if only_cal_with_snp and clair3_vcf_records_dict[key2].is_INDEL():
                    continue
                mat = consistency_dict[key1][key2]
                cnt_01_10 = mat[0][1] + mat[1][0]
                cnt_00_11 = mat[0][0] + mat[1][1]
                if cnt_00_11 + cnt_01_10 < n_minimum_covered_reads:
                    continue
                #linked_sites.append(key2)
                consistency = max(cnt_00_11, cnt_01_10) / (cnt_00_11 + cnt_01_10)
                consistency_list.append(consistency)
                linked_sites.append(f"{key2}:{consistency}")
                weight_list.append(0.5 if clair3_vcf_records_dict[key2].is_INDEL() else 1)

            if len(consistency_list) == 0:
                if clair3_vcf_records_dict[key1].is_INDEL():
                    single_indel_num += 1
                else:
                    single_snp_num += 1
            if not filter_snp and not clair3_vcf_records_dict[key1].is_INDEL():
                consistency_key1 = 0
                weight_sum = 0
                for consistency_value, weight in list(zip(consistency_list, weight_list)):
                    consistency_key1 += consistency_value * weight
                    weight_sum += weight
                if len(consistency_list) > 0:
                    consistency_key1 = consistency_key1 /  weight_sum
                    consistency_info_list.append((key1, str(clair3_vcf_records_dict[key1].qual), str(consistency_key1), str(len(consistency_list)),  "PASS", "AUTO", "INDEL" if clair3_vcf_records_dict[key1].is_INDEL() else "SNP", ", ".join(linked_sites)))
                else:
                    consistency_info_list.append((key1, str(clair3_vcf_records_dict[key1].qual), "-", str(len(consistency_list)),  "PASS", "AUTO", "INDEL" if clair3_vcf_records_dict[key1].is_INDEL() else "SNP", ", ".join(linked_sites)))
                
                continue


            if len(consistency_list) >= n_min_variant_position:
                # avg_consistency = np.mean(consistency_list)
                # max_consistency = np.max(consistency_list)
                # consistency_key1 = np.mean(consistency_list) if use_avg_consistency else np.max(consistency_list)
                consistency_key1 = 0
                weight_sum = 0
                for consistency_value, weight in list(zip(consistency_list, weight_list)):
                    consistency_key1 += consistency_value * weight
                    weight_sum += weight
                consistency_key1 = consistency_key1 /  weight_sum

                if consistency_key1 < consistency_threshold:
                    filtered_keys.add(key1)
                    if clair3_vcf_records_dict[key1].is_INDEL():
                        filtered_indel_num += 1
                    else:
                        filtered_snp_num += 1
                    consistency_info_list.append((key1, str(clair3_vcf_records_dict[key1].qual), str(consistency_key1), str(len(consistency_list)),  "FAIL", "CONSISTENCY", "INDEL" if clair3_vcf_records_dict[key1].is_INDEL() else "SNP", ", ".join(linked_sites)))
                else:
                    consistency_info_list.append((key1, str(clair3_vcf_records_dict[key1].qual), str(consistency_key1), str(len(consistency_list)), "PASS", "CONSISTENCY", "INDEL" if clair3_vcf_records_dict[key1].is_INDEL() else "SNP", ", ".join(linked_sites)))
            else:
                if filter_not_enough_position:
                    if float(clair3_vcf_records_dict[key1].qual) < 15:
                        filtered_keys.add(key1)
                        filter_by_qual_num += 1
                        consistency_info_list.append((key1, str(clair3_vcf_records_dict[key1].qual), "-", str(len(consistency_list)), "FAIL", "QUAL", "INDEL" if clair3_vcf_records_dict[key1].is_INDEL() else "SNP", ", ".join(linked_sites)))
                    else:
                        consistency_info_list.append((key1, str(clair3_vcf_records_dict[key1].qual), "-", str(len(consistency_list)), "PASS", "QUAL", "INDEL" if clair3_vcf_records_dict[key1].is_INDEL() else "SNP", ", ".join(linked_sites)))
                else:
                    consistency_info_list.append((key1, str(clair3_vcf_records_dict[key1].qual), "-", str(len(consistency_list)), "PASS", "QUAL", "INDEL" if clair3_vcf_records_dict[key1].is_INDEL() else "SNP", ", ".join(linked_sites))) 
    # logger.info(f"{filter_by_qual_num} variants filtered by qual 10.")

    # logger.info(f"single snp number: {single_snp_num}")
    # logger.info(f"single indel number: {single_indel_num}")
    # logger.info(f"filtered snp number: {filtered_snp_num}")
    # logger.info(f"filtered indel number: {filtered_indel_num}")

    return set(filtered_keys), consistency_info_list

def filter_vcf_file_by_positions(vcf_fn, new_vcf_fn, positions_list):

    positions_set = set(positions_list)

    with open(new_vcf_fn, "w") as f_clair3:
        with open(vcf_fn, "r") as clair3:
            for line in clair3.readlines():
                if line.startswith("#"):
                    f_clair3.write(f"{line.strip()}\n")
                    continue
                items = line.split("\t")
                if not f"{items[0]}:{int(items[1])}" in positions_set:
                    f_clair3.write(f"{line.strip()}\n")

def statistic_filtered_keys(vcf_fn, filtered_keys):
    clair3_vcf_records_list = read_pepper_vcf(vcf_fn)
    clair3_vcf_records_dict = { vcf_record.key : vcf_record for vcf_record in clair3_vcf_records_list }
    filter_indel_num = 0
    filter_snp_num = 0
    for k in filtered_keys:
        if clair3_vcf_records_dict[k].is_INDEL():
            filter_indel_num += 1
        else:
            filter_snp_num += 1

    # logger.info(f"Total filtered snp number: {filter_snp_num}")
    # logger.info(f"Total filtered indel number: {filter_indel_num}")

def consistency_analysis_step_1(consistency_dict, n_minimum_covered_reads, n_min_variant_position, consistency_threshold, vcf_fn, filter_snp=True):
    return filter_by_consistency(consistency_dict, n_minimum_covered_reads, n_min_variant_position, consistency_threshold, vcf_fn, filter_not_enough_position=False,filter_snp=filter_snp, use_avg_consistency=False)
def consistency_analysis_step_2(consistency_dict, n_minimum_covered_reads, n_min_variant_position, consistency_threshold, vcf_fn, filtered_keys, filter_snp=True):
    return filter_by_consistency(consistency_dict, n_minimum_covered_reads, n_min_variant_position, consistency_threshold, vcf_fn, filtered_keys=filtered_keys,filter_not_enough_position=True,filter_snp=filter_snp, use_avg_consistency=True)
def consistency_analysis_step_3(consistency_dict, n_minimum_covered_reads, n_min_variant_position, consistency_threshold, vcf_fn, filtered_keys, filter_snp=True):
    return filter_by_consistency(consistency_dict, n_minimum_covered_reads, n_min_variant_position, consistency_threshold, vcf_fn, filtered_keys=filtered_keys, filter_not_enough_position=True,filter_snp=filter_snp, use_avg_consistency=True)


def log_consistency_info_list(log_file, consistency_info_list):
    with open(log_file, "w") as f:
        for t in consistency_info_list:
            # t[0] = t[0] + "%"
            f.write("\t".join(t))
            f.write("\n")

def filter_main(VCF_FN, FRAGMENT_FN, OUTPUT_VCF, OUT_DIR, FILTER_SNP=True, ctg=None):

    # hapdup_log = os.path.join(OUT_DIR, "filter_main.log" if ctg is None else f"filter_main_{ctg}.log")
    # _enable_logging(hapdup_log, debug=False, overwrite=False)

    N_MIN_COVERED_READS_STEP1 = 5
    N_MIN_VARIANT_POSITION_STEP1 = 1
    MIN_CONSISTENCY_STEP1 = 0.7

    N_MIN_COVERED_READS_STEP2 = 10
    N_MIN_VARIANT_POSITION_STEP2 = 5
    MIN_CONSISTENCY_STEP2 = 0.8

    N_MIN_COVERED_READS_STEP3 = 10
    N_MIN_VARIANT_POSITION_STEP3 = 5
    MIN_CONSISTENCY_STEP3 = 0.9

    # N_MIN_COVERED_READS_STEP4 = 10
    # N_MIN_LR_VARIANT_POSITION_STEP4 = 1
    # MIN_CONSISTENCY_STEP4 = 0.8

    consistency_res_dict = calculate_consistency(VCF_FN, FRAGMENT_FN, OUT_DIR)
    # logger.info("============STEP1=============")
    filtered_keys_step1, consistency_info_list_step1 = consistency_analysis_step_1(consistency_res_dict, N_MIN_COVERED_READS_STEP1, N_MIN_VARIANT_POSITION_STEP1, MIN_CONSISTENCY_STEP1, VCF_FN, filter_snp=FILTER_SNP)
    log_consistency_info_list(f"{OUT_DIR}/consistency.info.step1", consistency_info_list_step1) 
    filter_vcf_file_by_positions(VCF_FN, OUTPUT_VCF + ".step1", filtered_keys_step1)
    # logger.info("============STEP2=============")
    filtered_keys_step2, consistency_info_list_step2 = consistency_analysis_step_2(consistency_res_dict, N_MIN_COVERED_READS_STEP2, N_MIN_VARIANT_POSITION_STEP2, MIN_CONSISTENCY_STEP2, VCF_FN, filtered_keys_step1, filter_snp=FILTER_SNP)
    log_consistency_info_list(f"{OUT_DIR}/consistency.info.step2", consistency_info_list_step2)
    filter_vcf_file_by_positions(VCF_FN, OUTPUT_VCF + ".step2", filtered_keys_step2)
    # logger.info("============STEP3=============")
    filtered_keys_step3, consistency_info_list_step3 = consistency_analysis_step_3(consistency_res_dict, N_MIN_COVERED_READS_STEP3, N_MIN_VARIANT_POSITION_STEP3, MIN_CONSISTENCY_STEP3, VCF_FN, filtered_keys_step2, filter_snp=FILTER_SNP)
    log_consistency_info_list(f"{OUT_DIR}/consistency.info.step3", consistency_info_list_step3)


    filter_vcf_file_by_positions(VCF_FN, OUTPUT_VCF, filtered_keys_step3)
if __name__ == '__main__':
    DIR="/public/home/hpc224712204/software/NanoFilter/test_out/tmp/chr14_KI270724v1_random"
    filter_main(f"{DIR}/PEPPER_VARIANT_FULL.vcf.chr14_KI270724v1_random", f"{DIR}/fragment_file", f"{DIR}/PEPPER_VARIANT_FULL.vcf.chr14_KI270724v1_random.filtered", DIR, FILTER_SNP=False, ctg=None)
    
