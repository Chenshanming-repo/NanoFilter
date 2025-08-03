def read_contig_info_from_phaseblock_name(phaseblock_name):
    if not phaseblock_name.startswith("hap"):
        phaseblock_name = f"hap0_{phaseblock_name}"
    parts = phaseblock_name.split("_")
    hap = int(parts[0][-1])
    contig_name = f"{parts[1]}_{parts[2]}"
    phaseblock_no = -1
    if "phaseblock" in phaseblock_name:
        phaseblock_no = int(parts[-1])

    return hap, contig_name, phaseblock_no

if __name__ == '__main__':
    print(read_contig_info_from_phaseblock_name("hap2_contig_1060_phaseblock_0"))
    print(read_contig_info_from_phaseblock_name("hap2_contig_1060"))