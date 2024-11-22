import os
import argparse
from hapcut2_multiprocess import run_hapcut2_filter_multiprocess

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

    args = parser.parse_args()

    if not os.path.isdir(args.out_dir):
        os.mkdir(args.out_dir)

    tmp_dir = os.path.join(args.out_dir, "tmp")
    if not os.path.isdir(tmp_dir):
        os.mkdir(tmp_dir)

    out_vcf = os.path.join(args.out_dir, "filtered.vcf")

    run_hapcut2_filter_multiprocess(args.reference, args.bam, args.vcf,
                                    out_vcf, tmp_dir,
                                    False, "ont", int(args.threads), FILTER_SNP=args.filter_all)




