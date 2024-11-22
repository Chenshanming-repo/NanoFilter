import os
import subprocess

import pysam

from MultiStepConsistencyAnalysisHapDup import filter_main


def bool_to_int(bool_var):
    return 1 if bool_var else 0

def run_hapcut2_contigs(ctg, ref, bam, vcf, out_dir, only_snvs, rtype, only_extract_hairs):
    EXTRACT_HAIRS=f"extractHAIRS"
    HAPCUT2=f"HAPCUT2"
    FRAGMENT_FILE = os.path.abspath(os.path.join(out_dir, "fragment_file"))
    hapcut2_log = os.path.abspath(os.path.join(out_dir, "hapcut2.log"))
    phase_log = os.path.abspath(os.path.join(out_dir, "phase.log"))

    vcf_ctg = f"{out_dir}/{os.path.split(vcf)[1]}.{ctg}"

    subprocess.check_call(f"bcftools view -o {vcf_ctg} {vcf + '.gz'} {ctg}", shell=True)


    extract_hairs_cmd = [EXTRACT_HAIRS, "--indels", str(bool_to_int(not only_snvs)), "--ont" if rtype == "ont" or rtype == "ont_r10"  else "--pacbio", "1",
                         "--bam", os.path.abspath(bam), "--VCF", vcf_ctg,
                         "--out", os.path.abspath(FRAGMENT_FILE),
                         "--ref", os.path.abspath(ref),
                         "--region", ctg, "2>&1", f"|tee {hapcut2_log}"]

    subprocess.check_call(" ".join(extract_hairs_cmd), shell=True)

    if only_extract_hairs:
        return

    PHASED_VCF = os.path.abspath(os.path.join(out_dir, "haplotype_output_file.phased.VCF"))
    PHASED_VCF_GZ = os.path.abspath(os.path.join(out_dir, "haplotype_output_file.phased.VCF.gz"))
    hapcut2_phase_cmd = [HAPCUT2, "--outvcf", "1", "--fragments", FRAGMENT_FILE, "--VCF", vcf_ctg,
                         "--output", os.path.abspath(os.path.join(out_dir, "haplotype_output_file")), "--verbose 1", "2>&1", f"|tee {phase_log}"]
    print("PHASING COMMAND:", " ".join(hapcut2_phase_cmd))
    subprocess.check_call(" ".join(hapcut2_phase_cmd), shell=True)

    subprocess.check_call(f"bgzip -c {PHASED_VCF} > {PHASED_VCF_GZ}", shell=True)
    subprocess.check_call(f"tabix -p vcf -f {PHASED_VCF_GZ}", shell=True)

    # subprocess.check_call(f"rm {vcf_ctg}", shell=True)

# def vcf_concat(out_dir, contigs, output_vcf):
#     add_head = False
#     with open(output_vcf, "w") as out_vcf:
#         for ctg in contigs:
#             with open(os.path.join(out_dir, ctg, "haplotype_output_file.phased.VCF")) as in_vcf:
#                 for line in in_vcf.readlines():
#                     if line.startswith("#"):
#                         if add_head:
#                             out_vcf.write(line.strip() + "\n")
#                             continue
#                     add_head = True
#                     out_vcf.write(line.strip() + "\n")

def file_check(path):
    if not os.path.isfile(path):
        # logger.error("Missing output: %s", path)
        raise Exception("Missing output: %s", path)

def run_hapcut2_filter_multiprocess(ref, bam, vcf, new_vcf, out_dir, only_snvs, rtype, threads_num, FILTER_SNP):
    from multiprocessing import Pool

    print(f"only extract hairs (filter), threads: {threads_num}")
    p = Pool(threads_num)


    for ctg in pysam.VariantFile(vcf + ".gz").header.contigs:
        out_dir_ctg = os.path.join(out_dir, ctg)
        if not os.path.isdir(out_dir_ctg):
            os.mkdir(out_dir_ctg)
        p.apply_async(run_hapcut2_contigs, args=(ctg, ref, bam, vcf, out_dir_ctg, only_snvs, rtype, True))


        #  run_hapcut2_contigs(ctg, ref, bam, vcf, out_dir_ctg, only_snvs, rtype)

    p.close()
    p.join()

    print(f"filter variant, threads: {threads_num}")

    p = Pool(threads_num)
    for ctg in pysam.VariantFile(vcf + ".gz").header.contigs:
        p.apply_async(filter_main, args=(f"{out_dir}/{ctg}/{os.path.split(vcf)[1]}.{ctg}", os.path.abspath(os.path.join(out_dir, ctg, "fragment_file")), f"{out_dir}/{ctg}/{os.path.split(vcf)[1]}.{ctg}.filtered", os.path.abspath(os.path.join(out_dir, ctg)), FILTER_SNP, ctg), error_callback=on_error)
    p.close()
    p.join()

    add_head = True
    with open(new_vcf, "w") as out_vcf:
        for ctg in pysam.VariantFile(vcf + ".gz").header.contigs:
            if not os.path.exists(os.path.join(out_dir, ctg, f"{os.path.split(vcf)[1]}.{ctg}.filtered")):
                continue
            with open(os.path.join(out_dir, ctg, f"{os.path.split(vcf)[1]}.{ctg}.filtered")) as in_vcf:
                for line in in_vcf.readlines():
                    if line.startswith("#"):
                        if add_head:
                            out_vcf.write(line.strip() + "\n")
                        continue
                    add_head = False
                    out_vcf.write(line.strip() + "\n")
    add_head = True
    with open(new_vcf + ".step1", "w") as out_vcf:
        for ctg in pysam.VariantFile(vcf + ".gz").header.contigs:
            if not os.path.exists(os.path.join(out_dir, ctg, f"{os.path.split(vcf)[1]}.{ctg}.filtered.step1")):
                continue
            with open(os.path.join(out_dir, ctg, f"{os.path.split(vcf)[1]}.{ctg}.filtered.step1")) as in_vcf:
                for line in in_vcf.readlines():
                    if line.startswith("#"):
                        if add_head:
                            out_vcf.write(line.strip() + "\n")
                        continue
                    add_head = False
                    out_vcf.write(line.strip() + "\n")
    
    add_head = True
    with open(new_vcf + ".step2", "w") as out_vcf:
        for ctg in pysam.VariantFile(vcf + ".gz").header.contigs:
            if not os.path.exists(os.path.join(out_dir, ctg, f"{os.path.split(vcf)[1]}.{ctg}.filtered.step2")):
                continue
            with open(os.path.join(out_dir, ctg, f"{os.path.split(vcf)[1]}.{ctg}.filtered.step2")) as in_vcf:
                for line in in_vcf.readlines():
                    if line.startswith("#"):
                        if add_head:
                            out_vcf.write(line.strip() + "\n")
                        continue
                    add_head = False
                    out_vcf.write(line.strip() + "\n")

def on_error(e):
    print(f"Error: {e}")

def run_hapcut2_multiprocess(ref, bam, vcf, new_vcf, out_dir, only_snvs, rtype, threads_num, only_extract_hairs=False):
    from multiprocessing import Pool

    p = Pool(threads_num)

    for ctg in pysam.VariantFile(vcf + ".gz").header.contigs:
        out_dir_ctg = os.path.join(out_dir, ctg)
        if not os.path.isdir(out_dir_ctg):
            os.mkdir(out_dir_ctg)
        p.apply_async(run_hapcut2_contigs, args=(ctg, ref, bam, vcf, out_dir_ctg, only_snvs, rtype, only_extract_hairs))
        #  run_hapcut2_contigs(ctg, ref, bam, vcf, out_dir_ctg, only_snvs, rtype)

    p.close()
    p.join()

    if not only_extract_hairs:

        # create list of input vcf files
        with open(os.path.join(out_dir, "vcf_list.txt"), "w") as f:
            for ctg in pysam.VariantFile(vcf + ".gz").header.contigs:
                f.write(os.path.join(out_dir, ctg, "haplotype_output_file.phased.VCF.gz"))
                f.write("\n")
        vcf_concat_cmd = ["bcftools", "concat", "-f", os.path.join(out_dir, "vcf_list.txt"), ">", new_vcf]
        subprocess.check_call(" ".join(vcf_concat_cmd), shell=True)

        file_check(new_vcf)
    # else:
    #     print("only extract hairs")
    #     fragment_cat_cmd = ["cat", ]
    #     for ctg in pysam.VariantFile(vcf + ".gz").header.contigs:
    #         fragment_cat_cmd.append(os.path.join(out_dir, ctg, "fragment_file"))
    #     fragment_cat_cmd.append(">")
    #     fragment_cat_cmd.append(os.path.join(out_dir, "fragment_file"))
    #     subprocess.check_call(" ".join(fragment_cat_cmd), shell=True)

    # vcf_concat(out_dir, pysam.VariantFile(vcf + ".gz").header.contigs, new_vcf)


