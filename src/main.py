import os
import argparse
import shutil
from filter_multiprocess import run_filter_multiprocess
from Utils import float_range
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Filtering out low-consistency variants from VCF")

    parser.add_argument("--reference", dest="reference",
                        metavar="path", required=True,
                        help="path to reference")
    parser.add_argument("--bam", dest="bam", required=True, metavar="path",
                        default=None, help="path to the alignment of reads on the ref in bam format")
    parser.add_argument("--vcf", dest="vcf", required=True, metavar="path",
                        default=None, help="path to the raw vcf")
    parser.add_argument("--out-dir", dest="out_dir", required=True, metavar="path",
                        default=None, help="path to the output directory")

    parser.add_argument("-t", "--threads", dest="threads", type=int,
                        default=10, metavar="int", help="number of parallel threads [10]")


    # additional
    parser.add_argument("--filter-all", dest="filter_all",
                        default=False, action="store_true",
                        help="Use consistency to filter SNP and INDEL")
    parser.add_argument("--filter-indel", dest="filter_indel",
                        default=False, action="store_true",
                        help="Use consistency to filter INDEL, and remain SNP")

    # filter parameter
    parser.add_argument("--filter-parameter1", "-p1", dest="p_step1", type=float_range,
                        default=0.7,
                        help="Filtering parameter of step 1")
    parser.add_argument("--filter-parameter2", "-p2", dest="p_step2", type=float_range,
                        default=0.8,
                        help="Filtering parameter of step 2")
    parser.add_argument("--filter-parameter3", "-p3", dest="p_step3", type=float_range,
                        default=0.9,
                        help="Filtering parameter of step 3")
    parser.add_argument("--quality", "-q", dest="qual", type=float_range,
                        default=15,
                        help="Filtering parameter of quality in low-coverage region where consistency is not convincible")

    args = parser.parse_args()

    print(f"[INFO] Filter parameters: {(args.p_step1, args.p_step2, args.p_step3)}")

    for file_path in [args.reference, args.bam, args.vcf]:
        if not os.path.exists(file_path):
            print(f"[ERROR] {file_path} not exists.")

    if not os.path.isdir(args.out_dir):
        os.mkdir(args.out_dir)

    tmp_dir = os.path.join(args.out_dir, "tmp")
    if not os.path.isdir(tmp_dir):
        os.mkdir(tmp_dir)

    out_vcf = os.path.join(args.out_dir, "filtered.vcf")

    run_filter_multiprocess(args.reference, args.bam, args.vcf,
                            out_vcf, tmp_dir,
                            False, "ont", int(args.threads), args.filter_all, (args.p_step1, args.p_step2, args.p_step3, args.qual))

    try:
        shutil.rmtree(tmp_dir)
        print(f'Temporary folder({tmp_dir}) and its content removed')
    except:
        print(f'Temporary folder({tmp_dir}) not deleted')


