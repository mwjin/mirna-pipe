"""
This module contains definitions and functions used in this project frequently.
"""

__all__ = ['chrom_list', 'chr_order']

_CHROM_LIST = ["chr%s" % x for x in range(1, 23)] + ['chrX', 'chrY']
_CHROM_TO_ORDER = {chrom: (i + 1) for i, chrom in enumerate(_CHROM_LIST)}  # the order of chromosomes on the BED files


def chrom_list():
    return _CHROM_LIST


def chr_order(chrom):
    return _CHROM_TO_ORDER[chrom]
