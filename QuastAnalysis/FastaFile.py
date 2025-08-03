from typing import DefaultDict

from Bio import SeqIO
class FastaFile:
    def __init__(self, fasta_path):
        self.fasta_path = fasta_path
        self.data = DefaultDict()
        for seq in SeqIO.parse(fasta_path, "fasta"):
            if not seq.id in self.data:
                self.data[seq.id] = seq.seq

    def get_ids(self):
        return self.data.keys()

    def get_sequence(self, seq_id):
        return self.data[seq_id]


if __name__ == '__main__':
    f = FastaFile(r"C:\Users\Xiaotuantuan\Documents\quast_2\r10\raw_snp\hapdup_phased_1.fasta")
    print(f.get_sequence("contig_1_phaseblock_0"))