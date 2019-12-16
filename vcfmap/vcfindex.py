import logging
logging.basicConfig(level=logging.INFO)
import numpy as np
import gzip

"""
Enabling lookup of a numpy array of 0/1s, telling which haplotypes support a given variant
"""


class VcfIndex:
    def __init__(self, haplotypes):
        self._haplotypes = haplotypes
        pass

    def haplotypes_at_variant(self, variant_number):
        return self._haplotypes[variant_number]

    def n_variants(self):
        return len(self._haplotypes)

    def to_file(self, file_name):
        np.save(file_name, self._haplotypes)

    @classmethod
    def from_file(cls, file_name):
        data = np.load(file_name)
        return cls(data)

    @classmethod
    def from_vcf_file(cls, file_name, n_variants):
        variants = None
        n_haplotypes = None

        with gzip.open(file_name) as f:
            for j, line in enumerate(f):
                line = line.decode("utf-8")
                if line.startswith("#"):
                    if line.startswith("#CHROM"):
                        n_haplotypes = 2*(len(line.split()) - 9)
                        logging.info("Making matrix with %d variants and %d haplotypes" % (n_variants, n_haplotypes))
                        variants = np.zeros((n_variants, n_haplotypes))
                    continue

                if j % 1000 == 0:
                    logging.info("%d lines processed" % j)

                line = line.split()
                haplotypes = np.zeros(n_haplotypes)

                for i, genotype in enumerate(line[9:]):
                    haplotype0_id = i * 2
                    haplotype1_id = i * 2 + 1

                    if genotype == "0|1":
                        haplotypes[haplotype0_id] = 0
                        haplotypes[haplotype1_id] = 1
                    elif genotype == "1|1":
                        haplotypes[haplotype0_id] = 1
                        haplotypes[haplotype1_id] = 1
                    elif genotype == "0|0":
                        haplotypes[haplotype0_id] = 0
                        haplotypes[haplotype1_id] = 0
                    elif genotype == "1|0":
                        haplotypes[haplotype0_id] = 1
                        haplotypes[haplotype1_id] = 0

                variants[j] = haplotypes

        return cls(variants)
