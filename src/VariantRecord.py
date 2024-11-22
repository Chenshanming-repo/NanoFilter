from enum import Enum

DEFAULT_STRING = "UNKNOWN"
DEFAULT_INTEGER = -1
DEFAULT_FLOAT = -1.0

class VariantType(Enum):
    VARIANT_UNKNOWN = -1
    VARIANT_REFERENCE = 0
    VARIANT_SNP = 1
    VARIANT_INDEL = 2
    VARIANT_METH = 3
    VARIANT_SV = 4

class VariantRecord:
    def __init__(self, ctg : str, pos : int, allele1 : str, allele2 : str , qual : int, genotype : str, variant_type : VariantType, index : int, phased_gt = None):
        self.ctg = ctg
        self.pos = pos
        self.allele1 = allele1
        self.allele2 = allele2
        self.genotype = genotype
        self.qual = qual
        self.index = index
        self.variant_type = variant_type
        self.phased = False
        self.key = f"{ctg}:{pos}"

        if "|" not in self.genotype:
            if not self.genotype in ["0/1", "1/0", "1/2", "2/1", "0/2", "2/0", "1/1"]:
                # print(self.genotype)
                pass
            # assert self.genotype in ["0/1", "1/0", "1/2", "2/1", "0/2", "2/0", "1/1"]
        else:
            self.phased = True

    def get_ctg(self):
        return self.ctg
    def get_pos(self):
        return self.pos


    """
        输出是否需要重比对
    """
    def reaalignable(self) -> bool:
        return self.variant_type in [VariantType.VARIANT_SNP, VariantType.VARIANT_INDEL, VariantType.VARIANT_SV] and self.genotype in ["0/1"]

    def is_INDEL(self) -> bool:
        return self.variant_type == VariantType.VARIANT_INDEL

    def is_phased(self) -> bool:
        return self.phased

    def is_heter(self):
        return "0" in self.genotype and not self.genotype[0] == self.genotype[1]


    # def get_hap_seq




