#!/opt/ohpc/pub/bin/anaconda3/envs/minu/bin/python3
"""
For each miRNA, get the target sites by referring 3UTR sequence of genes and
documents these target sites as a bed file format.

Prerequisite
    1. Run 01_gene_filter.py from 'repr-gene' project
    2. Run 01_preprocess.py
"""
from baeklab.genomics.microrna import Mature
from baeklab.genomics.sequence import Seq
from lab.bed.narrow import NarrowPeak
from lab.utils import eprint
from settings import *

import os
import pickle


def main():
    """
    Bootstrap
    """
    # input path settings
    mirna_fa_path = '%s/data/mature_human.fa' % PROJECT_DIR
    mirna_list_path = '%s/data/cons_mature_human_list.txt' % PROJECT_DIR
    refflat_data_path = '/extdata6/Minwoo/projects/repr-gene/results/repr-iso/%s/refFlat_171015_NM.dat' % GENOME_VER

    # output path settings
    mirna_bed_dir = '%s/results/target-sites/%s/bed' % (PROJECT_DIR, GENOME_VER)
    all_mir_bed_path = '%s/high_cons_mir.bed' % mirna_bed_dir  # a file for merged peaks of all miRNA target sites
    os.makedirs(mirna_bed_dir, exist_ok=True)

    make_mirna_target_bed(mirna_bed_dir, mirna_fa_path, mirna_list_path, refflat_data_path)
    merge_target_site(all_mir_bed_path, mirna_bed_dir, mirna_list_path)


def make_mirna_target_bed(result_dir, mirna_fa_path, mirna_list_path, refflat_data_path):
    """
    Make BED files documented miRNA target sites for representative isoforms
    :param result_dir: a result directory
    :param mirna_fa_path: a path of miRNA fasta file from miRBase
    :param mirna_list_path: a path of a list of miRNAs currently focused on
    :param refflat_data_path: a path of the data file (.dat) containing a list of genes (representative isoforms)
    """
    eprint('[LOG] Make bed files for miRNA target sites')
    mirnas = get_mirnas(mirna_fa_path, mirna_list_path)

    with open(refflat_data_path, 'rb') as gene_file:
        genes = pickle.load(gene_file)

    for mirna in mirnas:
        eprint('[LOG] --- miRNA: %s' % mirna.name)
        target_site_peaks = []  # element: a 'NarrowPeak' object

        for gene in genes:
            coords_3utr = get_3utr_coord(gene)
            rna_3utr_seq = Seq(gene.seq_3utr.replace('T', 'U'))
            target_sites = mirna.find_targetsites(rna_3utr_seq, 'cst')

            # convert target site coordinates to peaks
            for target_site in target_sites:
                if target_site.type == '6mer':  # skip the 6mer
                    continue

                target_coords = find_target_coord(target_site, coords_3utr)

                if len(target_coords) == 1:
                    target_name = '%s;%s;%s;%d' % (gene.symbol, gene.id, target_site.type, 0)
                    target_peak = NarrowPeak(gene.chrom, target_coords[0][0], target_coords[0][1], gene.strand,
                                             name=target_name)
                    target_site_peaks.append(target_peak)
                else:
                    for i, coord in enumerate(target_coords):
                        target_name = '%s;%s;%s;%d' % (gene.symbol, gene.id, target_site.type, i + 1)
                        target_peak = NarrowPeak(gene.chrom, target_coords[0][0], target_coords[0][1], gene.strand,
                                                 name=target_name)
                        target_site_peaks.append(target_peak)

        target_site_peaks.sort(key=lambda peak: (peak.chrom[3:], peak.start, peak.end))

        # make bed files for the target sites
        mirna_bed_path = '%s/%s.bed' % (result_dir, mirna.name)

        with open(mirna_bed_path, 'w') as mirna_bed_file:
            for target_site_peak in target_site_peaks:
                print(target_site_peak, file=mirna_bed_file)
    eprint()


def get_mirnas(mirna_fa_path, mirna_list_path):
    """
    Collect miRNAs only in our list of miRNAs and return the miRNA objects as 'Mature' objects
    :param mirna_fa_path: a path of miRNA fasta file from miRBase
    :param mirna_list_path: a path of a list of miRNAs currently focused on
    :return a list of 'Mature' objects
    """
    with open(mirna_list_path, 'r') as mirna_list_file:
        mirna_names = mirna_list_file.read().splitlines()

    with open(mirna_fa_path, 'r') as mirna_file:
        mir_name_to_mirna = {}

        while True:
            mir_header = mirna_file.readline()

            if not mir_header:
                break

            mir_seq = mirna_file.readline().strip()
            mir_name = mir_header.split()[0][1:]
            mirna = Mature(mir_name, mir_seq)
            mir_name_to_mirna[mir_name] = mirna

    mirnas = []

    for mir_name in mirna_names:
        mirnas.append(mir_name_to_mirna[mir_name])

    return mirnas


def gene_filter(genes, gene_candid_ids):
    """
    Remain only genes whose ID is in the gene candidate IDs
    :param genes: a list of 'RefFlat' objects
    :param gene_candid_ids: NM ID of the cancer genes
    :return: a list of 'RefFlat' objects with NM ID in the 'gene_candid_ids'
    """
    genes.sort(key=lambda x: int(x.id[3:]))
    gene_candid_ids.sort(key=lambda x: int(x[3:]))

    gene_cnt = len(genes)
    candid_id_cnt = len(gene_candid_ids)

    gene_idx = 0
    candid_id_idx = 0

    gene_candids = []

    while gene_idx < gene_cnt and candid_id_idx < candid_id_cnt:
        gene_id_num = int(genes[gene_idx].id[3:])
        candid_id_num = int(gene_candid_ids[candid_id_idx][3:])

        if gene_id_num < candid_id_num:
            gene_idx += 1
        elif gene_id_num == candid_id_num:
            gene_candids.append(genes[gene_idx])
            gene_idx += 1
            candid_id_idx += 1
        else:
            candid_id_idx += 1

    return gene_candids


def get_3utr_coord(gene):
    """
    Return the 3'UTR coordinates of the 'RefFlat' object on the genome
    :param gene: a 'RefFlat' object
    :return: a list of (start (0-based), end)
    """
    assert gene.__class__.__name__ == 'RefFlat'

    coords = []

    if gene.strand == '+':
        exon_idx_3utr_start = -1

        for i in range(gene.exon_cnt):
            if gene.exon_starts[i] < gene.cds_end <= gene.exon_ends[i]:
                exon_idx_3utr_start = i
                break

        assert exon_idx_3utr_start != -1

        if gene.cds_end < gene.exon_ends[exon_idx_3utr_start]:
            coords.append((gene.cds_end, gene.exon_ends[exon_idx_3utr_start]))

        for j in range(exon_idx_3utr_start + 1, gene.exon_cnt):
            coords.append((gene.exon_starts[j], gene.exon_ends[j]))

    else:  # gene.strand == '-'
        exon_idx_3utr_end = -1

        for i in range(gene.exon_cnt):
            if gene.exon_starts[i] <= gene.cds_start < gene.exon_ends[i]:
                exon_idx_3utr_end = i
                break

        assert exon_idx_3utr_end != -1

        for j in range(exon_idx_3utr_end):
            coords.append((gene.exon_starts[j], gene.exon_ends[j]))

        if gene.exon_starts[exon_idx_3utr_end] < gene.cds_start:
            coords.append((gene.exon_starts[exon_idx_3utr_end], gene.cds_start))

    return coords


def find_target_coord(target_site, coords_3utr):
    """
    By using the relative position of miRNA target site to the 3'UTR,
    find absolute coordinates of the miRNA target site on the genome and return it
    :param target_site: an object of 'Coor' that represents miRNA target site
    :param coords_3utr: a list of (start, end) (position of the 3'UTR parts of exons)
    :return: a list of (start, end)
    """
    acc_len_list = [0]  # dummy item
    acc_len = 0

    for coord in coords_3utr:
        acc_len += (coord[1] - coord[0])
        acc_len_list.append(acc_len)

    # indexing
    target_start_idx = -1
    target_end_idx = -1
    coord_cnt = len(coords_3utr)

    for i in range(coord_cnt):
        if acc_len_list[i] <= target_site.start < acc_len_list[i + 1]:
            target_start_idx = i

        if acc_len_list[i] < target_site.end <= acc_len_list[i + 1]:
            target_end_idx = i

    assert target_start_idx != -1 and target_end_idx != -1
    assert target_start_idx <= target_end_idx

    # get absolute coordinate
    target_start_coord = target_site.start - acc_len_list[target_start_idx] + coords_3utr[target_start_idx][0]
    target_end_coord = target_site.end - acc_len_list[target_end_idx] + coords_3utr[target_end_idx][0]

    target_coords = []

    if target_start_idx == target_end_idx:
        target_coords.append((target_start_coord, target_end_coord))
    else:
        target_coords.append((target_start_coord, coords_3utr[target_start_idx][1]))

        for i in range(target_start_idx + 1, target_end_idx):
            target_coords.append(coords_3utr[i])

        target_coords.append((coords_3utr[target_end_idx][0], target_end_coord))

    return target_coords


def merge_target_site(merge_bed_path, mirna_bed_dir, mirna_list_path):
    """
    Merge all target sites of miRNAs in our miRNA list and save the merged target sites as a BED file
    :param merge_bed_path: a path of the file for the merged peaks
    :param mirna_bed_dir: a directory BED files for miRNA target sites are saved
    :param mirna_list_path: a path of the file containing out miRNAs' names
    """
    eprint('[LOG] Merge all bed files for miRNA target sites')
    with open(mirna_list_path, 'r') as mirna_list_file:
        mirna_names = mirna_list_file.read().splitlines()

    eprint('[LOG] --- Concatenate all bed files')
    concat_cmd = 'cat'
    sort_cmd = 'sort -k1,1 -k2,2n -k3,3n -V'

    for mirna_name in mirna_names:
        concat_cmd += ' %s/%s.bed' % (mirna_bed_dir, mirna_name)

    concat_bed_path = '%s/temp.bed' % mirna_bed_dir  # it is a temporary file
    cmd = '%s | %s > %s;' % (concat_cmd, sort_cmd, concat_bed_path)
    os.system(cmd)

    eprint('[LOG] --- Merge all target sites using bedtools')
    merge_cmd = 'bedtools merge -s -c 4,5,6,7,8,9,10 -o distinct -i %s > %s;' % (concat_bed_path, merge_bed_path)
    os.system(merge_cmd)
    os.system('rm %s' % concat_bed_path)


if __name__ == '__main__':
    main()
