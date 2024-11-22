def read_fragment_file(fragment_fn):

    fragment_dict = {}

    with open(fragment_fn) as f:
        for line in f.readlines():
            line_items = line.strip().split(" ")
            # blk_num = int(line_items[0])
            fragment_dict[line_items[1]] = (" ".join(line_items[2:-1]), line_items[-1], int(line_items[0]))
    return fragment_dict

def read_fragment_file_list(fragment_fn):

    fragment_list = []

    with open(fragment_fn) as f:
        for line in f.readlines():
            line_items = line.strip().split(" ")
            # blk_num = int(line_items[0])
            fragment_list.append((" ".join(line_items[2:-1]), line_items[-1], int(line_items[0])))
    return fragment_list

def merge(read_fragment_id, fragment1, fragment2):
    blocks1, qual1, _ = fragment1
    blocks2, qual2, _ = fragment2
    # DEBUG: indel质量改为.
    qual2 = "." * len(qual2)
    blocks1 = blocks1.strip().split(" ")
    blocks2 = blocks2.strip().split(" ")

    block_list = []

    assert len(blocks1) % 2 == 0 and len(blocks2) % 2 == 0

    for i in range(len(blocks1) // 2):
        start_idx, haplotypes = int(blocks1[2 * i]), blocks1[2 * i + 1]
        j = 0
        while j < len(haplotypes):
            block_list.append((start_idx + j, haplotypes[j], qual1[j]))
            j += 1
        qual1 = qual1[j:]
    for i in range(len(blocks2) // 2):
        start_idx, haplotypes = int(blocks2[2 * i]), blocks2[2 * i + 1]
        j = 0
        while j < len(haplotypes):
            block_list.append((start_idx + j, haplotypes[j], qual2[j]))
            j += 1
        qual2 = qual2[j:]


    block_list.sort(key=(lambda x : x[0]))

    blk_num = 0
    new_fragment_blocks = ""
    new_fragment_qual = ""

    if len(block_list) == 0:
        return ""

    last_idx = block_list[0][0]
    cur_haplotypes = f"{block_list[0][1]}"
    cur_start_idx = block_list[0][0]
    new_fragment_qual += block_list[0][2]

    i = 1
    while i < len(block_list):
        if block_list[i][0] == last_idx + 1:
            cur_haplotypes += block_list[i][1]
            if i == len(block_list) - 1:
                new_fragment_blocks += f"{cur_start_idx} {cur_haplotypes} "
                blk_num += 1
        else:
            new_fragment_blocks += f"{cur_start_idx} {cur_haplotypes} "
            blk_num += 1
            cur_haplotypes = f"{block_list[i][1]}"
            cur_start_idx = block_list[i][0]
        last_idx = block_list[i][0]
        new_fragment_qual += block_list[i][2]
        i += 1

    return f"{blk_num} {read_fragment_id} {new_fragment_blocks.strip()} {new_fragment_qual}"



if __name__ == '__main__':
    indel_fragment_fn = "/home/chensm/tmp/fragment_files/indel.fragment"
    snp_fragment_fn = "/home/chensm/tmp/fragment_files/snp.fragment"
    merge_fragment_fn = "/home/chensm/tmp/fragment_files/indel_snp.fragment"

    snp_fragment_dict = read_fragment_file(snp_fragment_fn)
    indel_fragment_dict = read_fragment_file(indel_fragment_fn)
    with open(merge_fragment_fn, "w") as merge_f:
        for read_id in snp_fragment_dict:
            snp_fragment = snp_fragment_dict[read_id]
            if not read_id in indel_fragment_dict:
                # indel_fragment = ""
                merge_f.write(f"{snp_fragment[2]} {read_id} {snp_fragment[0]} {snp_fragment[1]}\n")
                continue

            indel_fragment = indel_fragment_dict[read_id]
            # indel_fragment = indel_fragment_dict[read_id]

            merge_f.write(merge(read_id, snp_fragment, indel_fragment) + "\n")


