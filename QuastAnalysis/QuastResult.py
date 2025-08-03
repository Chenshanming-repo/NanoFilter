import os
import re
from typing import DefaultDict
from bs4 import BeautifulSoup
import json5
from tqdm import tqdm
from FastaFile import FastaFile
from Utils import read_contig_info_from_phaseblock_name
from VariantDensity import cal_variant_density_region, output_variant_density, read_vcf_file, \
    output_contig_variant_density


class QuastResult:
    def __init__(self, quast_path, assembly_file_1, assembly_file_2):
        self.quast_path = quast_path
        #self.assembly_1 = self.assembly_1 = None
        #self.assembly_file_1 = assembly_file_1
        #self.assembly_file_2 = assembly_file_2

        self.misassemblies_num = 0

        self.chromosomes = []
        for i in range(1, 23):
            self.chromosomes.append(f"chr{i}_MATERNAL")
            self.chromosomes.append(f"chr{i}_PATERNAL")
        self.chromosomes.append("chrX_MATERNAL")
        self.chromosomes.append("chrY_PATERNAL")

        self.phaseblock_in_chromosomes = DefaultDict()
        self.data = DefaultDict()
        for chromosome in self.chromosomes:
            self.data[chromosome] = DefaultDict()
        self.read_html()

        # read assembly, record the order and length of phaseblocks
        print(f"Loading assembly file {assembly_file_1}")
        self.assembly_1 = FastaFile(assembly_file_1)
        print(f"Loading assembly file {assembly_file_2}")
        self.assembly_2 = FastaFile(assembly_file_2)

        self.assembly_phaseblock_1 = DefaultDict()
        self.assembly_phaseblock_2 = DefaultDict()

        for phaseblock_id in self.assembly_1.get_ids():
            hap, contig_name, phaseblock_no = read_contig_info_from_phaseblock_name(phaseblock_id)
            if not contig_name in self.assembly_phaseblock_1:
                self.assembly_phaseblock_1[contig_name] = []
            self.assembly_phaseblock_1[contig_name].append((phaseblock_no, phaseblock_id, len(self.assembly_1.get_sequence(phaseblock_id))))
        for contig_name in self.assembly_phaseblock_1.keys():
            self.assembly_phaseblock_1[contig_name].sort(key=lambda x: x[0])

        for phaseblock_id in self.assembly_2.get_ids():
            hap, contig_name, phaseblock_no = read_contig_info_from_phaseblock_name(phaseblock_id)
            if not contig_name in self.assembly_phaseblock_2:
                self.assembly_phaseblock_2[contig_name] = []
            self.assembly_phaseblock_2[contig_name].append((phaseblock_no, phaseblock_id, len(self.assembly_2.get_sequence(phaseblock_id))))
        for contig_name in self.assembly_phaseblock_2.keys():
            self.assembly_phaseblock_2[contig_name].sort(key=lambda x: x[0])


    def test(self):
        # 是否有同一个contig在一个染色体下有两个hap的phaseblock
        #for chromosome in self.data:
        #    for contig_name in self.data[chromosome]:
        #        phaseblock_no_dict = DefaultDict()
        #        for phaseblock in self.data[chromosome][contig_name]:
        #            phaseblock_no = phaseblock.split("_")[-1]
        #            if phaseblock_no in phaseblock_no_dict.keys():
        #                print(f"{chromosome} {contig_name} {phaseblock} {phaseblock_no_dict[phaseblock_no]}")
        #            else:
        #                phaseblock_no_dict[phaseblock_no] = phaseblock
        # 测试contig region
        print(self.get_phaseblock_region("hap2_contig_1194_phaseblock_0"))
        print(self.get_phaseblock_region("hap2_contig_1194_phaseblock_4"))
        print(self.get_phaseblock_region("hap1_contig_1194_phaseblock_0"))
        print(self.get_phaseblock_region("hap1_contig_1194_phaseblock_3"))

    def read_html(self):
        phaseblock_name_pattern = r'contig_lengths\[[^\]]+\]\["([^"]+)"\]'

        for chromosome in tqdm(self.chromosomes):
            with open(os.path.join(self.quast_path, "icarus_viewers", f"{chromosome}.html"), "r", encoding="utf-8") as html:

                soup = BeautifulSoup(html, "html.parser")
                scripts = soup.find_all("script")
                script_contig = ""
                for script in scripts:
                    if script.text.startswith("var references_by_id = {};"):
                        script_contig = script.text

                lines = script_contig.split("\n")
                line_no = 0
                while line_no < len(lines):
                    line = lines[line_no]
                    if line.startswith("contig_lengths["):
                        match = re.search(phaseblock_name_pattern, line)
                        if match is None:
                            line_no += 1
                            continue
                        phaseblock_name = match.group(1)
                        line_no += 2
                        json_str = "["
                        while line_no < len(lines) and not line == "];":
                            line = lines[line_no]
                            json_str += line
                            line_no += 1
                        jsons = json5.loads(json_str[:-1])
                        self.add_phaseblock(phaseblock_name, jsons)

                    line_no += 1

    def add_phaseblock(self, phaseblock_name, phaseblock_jsons):

        hap, contig_name, phaseblock_no = read_contig_info_from_phaseblock_name(phaseblock_name)
        if len(phaseblock_jsons) == 0:
            return
        map_chromosome = self.chromosomes[int(phaseblock_jsons[0]["chr"])]

        json_results = []
        i = 0
        while i < len(phaseblock_jsons):
            if "contig_type" in phaseblock_jsons[i]:

                self.misassemblies_num += 1

                if i - 1 >= 0 and i + 1 < len(phaseblock_jsons):
                    phaseblock_jsons[i]["start_in_contig"] = max(phaseblock_jsons[i - 1]["start_in_contig"],
                                                                 phaseblock_jsons[i - 1]["end_in_contig"])
                    phaseblock_jsons[i]["end_in_contig"] = min(phaseblock_jsons[i + 1]["start_in_contig"],
                                                               phaseblock_jsons[i + 1]["end_in_contig"])
                    phaseblock_jsons[i]["last_fragment_length"] = abs(phaseblock_jsons[i - 1]["start_in_contig"] - phaseblock_jsons[i - 1]["end_in_contig"])
                    phaseblock_jsons[i]["next_fragment_length"] = abs(phaseblock_jsons[i + 1]["start_in_contig"] - phaseblock_jsons[i + 1]["end_in_contig"])

                else:
                    phaseblock_jsons[i]["start_in_contig"] = max(phaseblock_jsons[i - 1]["start_in_contig"],
                                                                 phaseblock_jsons[i - 1]["end_in_contig"])
                    phaseblock_jsons[i]["end_in_contig"] = phaseblock_jsons[i]["start_in_contig"] + 1
                    phaseblock_jsons[i]["last_fragment_length"] = abs(phaseblock_jsons[i - 1]["start_in_contig"] - phaseblock_jsons[i - 1]["end_in_contig"])
                    phaseblock_jsons[i]["next_fragment_length"] = 0

            json_results.append(phaseblock_jsons[i])
            i += 1


        if len(map_chromosome) == 0:
            print("Unmapped phaseblock: " + phaseblock_name)
            return

        if contig_name not in self.data[map_chromosome]:
            self.data[map_chromosome][contig_name] = DefaultDict()

        self.data[map_chromosome][contig_name][phaseblock_name] = json_results
        self.phaseblock_in_chromosomes[phaseblock_name] = map_chromosome

        # print(contig_name, phaseblock_jsons)

    def get_phaseblock_region(self, phaseblock_name):
        assert phaseblock_name.startswith("hap")
        hap, contig_name, phaseblock_no = read_contig_info_from_phaseblock_name(phaseblock_name)
        assembly_phaseblock_dict = self.assembly_phaseblock_1 if hap == 1 else self.assembly_phaseblock_2
        assembly_phaseblock_list = assembly_phaseblock_dict[contig_name]

        start = 0
        end = 0
        total_length = 0
        for phaseblock in assembly_phaseblock_list:
            total_length += phaseblock[2]
            if phaseblock[0] == phaseblock_no:
                end = total_length
                start = end - phaseblock[2] + 1


        # chrom = self.phaseblock_in_chromosomes[phaseblock_no]

        return contig_name, start, end


    def get_misassembly_fragments_in_region(self, chromosome, contig_name, start, end):
        if not contig_name in self.data[chromosome]:
            return []
        phaseblock_dict = self.data[chromosome][contig_name]
        fragments_results = []
        for phaseblock_name in phaseblock_dict.keys():
            _, _start, _end = self.get_phaseblock_region(phaseblock_name)
            for fragment in phaseblock_dict[phaseblock_name]:
                if not ("contig_type" in fragment and fragment["mstype"] == "real"):
                    continue
                start_in_phaseblock = min(fragment["start_in_contig"], fragment["end_in_contig"])
                end_in_phaseblock = max(fragment["start_in_contig"], fragment["end_in_contig"])
                fragment_start_in_contig = _start + start_in_phaseblock - 1
                fragment_end_in_contig = _start + end_in_phaseblock - 1

                if (fragment_start_in_contig <= start <= fragment_end_in_contig) or (fragment_start_in_contig <= end <= fragment_end_in_contig) or (start <= fragment_start_in_contig and end >= fragment_end_in_contig):
                    fragments_results.append((fragment, fragment_start_in_contig, fragment_end_in_contig, phaseblock_name))

        return fragments_results




    def compare(self, quast_results_compare, out_file):

        num = 0

        with open(out_file, "w", encoding="utf-8") as out:
            out.write(f"#Base quast results: {self.quast_path}\n")
            out.write(f"#compare quast results: {quast_results_compare.quast_path}\n")
            out.write("##CHROM\tBASE_PHASEBLOCK\tBASE_CONTIG_START\tBASE_CONTIG_END\tCOMP_PHASEBLOCK\tCOMP_CONTIG_START\tCOMP_CONTIG_END\tMIS_MSG\tLAST_FRAGMENT_LENGTH\tNEXT_FRAGMENT_LENGTH\n")

            for chromosome in self.chromosomes:
                for contig_name in self.data[chromosome]:
                    self_phaseblock_dict = self.data[chromosome][contig_name]
                    for phaseblock_name in self_phaseblock_dict.keys():
                        phaseblock = self_phaseblock_dict[phaseblock_name]

                        # 只要组装完整的phaseblock
                        has_misassemblies = False
                        for fragment in phaseblock:
                            if "contig_type" in fragment and fragment["contig_type"] == "M":
                                has_misassemblies = True
                        if has_misassemblies:
                            continue

                        _, start, end = self.get_phaseblock_region(phaseblock_name)
                        fragments = quast_results_compare.get_misassembly_fragments_in_region(chromosome, contig_name, start, end)
                        # print(fragments)
                        for fragment in fragments:
                            if "msg" in fragment[0]:
                                num += 1
                                out.write(f'{chromosome}\t{phaseblock_name}\t{start}\t{end}\t{fragment[3]}\t{fragment[1]}\t{fragment[2]}\t{fragment[0]["msg"]}\t{fragment[0]["last_fragment_length"]}\t{fragment[0]["next_fragment_length"]}\n')
        return num



def compare_main_tmp():
    init = True
    data_type="r9"
    quast_results_label = ["rs", "rsri", "rsfi", "fsfi"]
    vcf_file = rf"C:\Users\Xiaotuantuan\Documents\quast_2\{data_type}\PEPPER_VARIANT_FULL.PASS.vcf"

    if init:

        quast_result_rs = QuastResult(rf"C:\Users\Xiaotuantuan\Documents\quast_2\{data_type}\raw_snp\quast_results",
                                      rf"C:\Users\Xiaotuantuan\Documents\quast_2\{data_type}\raw_snp\hapdup_phased_1.fasta",
                                      rf"C:\Users\Xiaotuantuan\Documents\quast_2\{data_type}\raw_snp\hapdup_phased_2.fasta")
        # print(quast_result_rs.misassemblies_num)

        quast_result_rsri = QuastResult(rf"C:\Users\Xiaotuantuan\Documents\quast_2\{data_type}\raw_snp_raw_indel\quast_results",
                                        rf"C:\Users\Xiaotuantuan\Documents\quast_2\{data_type}\raw_snp_raw_indel\hapdup_phased_1.fasta",
                                        rf"C:\Users\Xiaotuantuan\Documents\quast_2\{data_type}\raw_snp_raw_indel\hapdup_phased_2.fasta")
        # print(quast_result_rsri.misassemblies_num)

        quast_result_rsfi = QuastResult(
            rf"C:\Users\Xiaotuantuan\Documents\quast_2\{data_type}\raw_snp_filtered_indel\quast_results",
            rf"C:\Users\Xiaotuantuan\Documents\quast_2\{data_type}\raw_snp_filtered_indel\hapdup_phased_1.fasta",
            rf"C:\Users\Xiaotuantuan\Documents\quast_2\{data_type}\raw_snp_filtered_indel\hapdup_phased_2.fasta")
        # print(quast_result_rsfi.misassemblies_num)

        quast_result_fsfi = QuastResult(
            rf"C:\Users\Xiaotuantuan\Documents\quast_2\{data_type}\filtered_snp_filtered_indel\quast_results",
            rf"C:\Users\Xiaotuantuan\Documents\quast_2\{data_type}\filtered_snp_filtered_indel\hapdup_phased_1.fasta",
            rf"C:\Users\Xiaotuantuan\Documents\quast_2\{data_type}\filtered_snp_filtered_indel\hapdup_phased_2.fasta")
        # print(quast_result_fsfi.misassemblies_num)

        quast_results = [quast_result_rs, quast_result_rsri, quast_result_rsfi, quast_result_fsfi]

        vcf_dict = read_vcf_file(vcf_file)

        for i in range(4):
            for j in range(4):
                if i == j:
                    continue
                out_compare_file = rf"C:\Users\Xiaotuantuan\PycharmProjects\QuastAnalysis\output\{data_type}\{data_type}_{quast_results_label[i]}_comp_{quast_results_label[j]}.txt"
                out_compare_density_file = rf"C:\Users\Xiaotuantuan\PycharmProjects\QuastAnalysis\output\{data_type}\{data_type}_{quast_results_label[i]}_comp_{quast_results_label[j]}.density.txt"

                num = quast_results[i].compare(quast_results[j], out_compare_file)
                output_variant_density(out_compare_file, vcf_dict, out_compare_density_file)

                print(quast_results_label[i], quast_results_label[j], num)

        output_contig_variant_density(vcf_dict)


    # rs_comp_rsri exists,   rsfi_comp_rsri exists, fsfi_comp_rsri exists

    # rs_comp_rsri exists, rs_comp_rsfi exists, rs_comp_fsfi exists
    results_dict = {}
    results_dict_rev = {}
    for label_compared in quast_results_label[1:]:
        with open(
                rf"C:\Users\Xiaotuantuan\PycharmProjects\QuastAnalysis\output\{data_type}\{data_type}_{quast_results_label[0]}_comp_{label_compared}.txt") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                fields = line.strip().split("\t")
                contig_name = read_contig_info_from_phaseblock_name(fields[1])[1]
                if not contig_name in results_dict:
                    results_dict[contig_name] = {"rsri": 0, "rsfi": 0, "fsfi": 0}
                results_dict[contig_name][label_compared] += 1
        with open(
                rf"C:\Users\Xiaotuantuan\PycharmProjects\QuastAnalysis\output\{data_type}\{data_type}_{label_compared}_comp_{quast_results_label[0]}.txt") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                fields = line.strip().split("\t")
                contig_name = read_contig_info_from_phaseblock_name(fields[1])[1]
                if not contig_name in results_dict_rev:
                    results_dict_rev[contig_name] = {"rsri": 0, "rsfi": 0, "fsfi": 0}
                results_dict_rev[contig_name][label_compared] += 1

    for contig_name in results_dict:
        if results_dict[contig_name]["rsri"] > 0 and results_dict[contig_name]["rsfi"] > 0 and \
                results_dict[contig_name]["fsfi"] > 0:
            if not contig_name in results_dict_rev:
                print(contig_name, results_dict[contig_name]["rsri"], results_dict[contig_name]["rsfi"],
                      results_dict[contig_name]["fsfi"])






if __name__ == '__main__':
    compare_main_tmp()

