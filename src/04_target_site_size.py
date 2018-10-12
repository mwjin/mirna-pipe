#!/opt/ohpc/pub/bin/anaconda3/envs/minu/bin/python3
"""
Get the size of RBP binding site peaks for each genic regions

Prerequisite
    1. Run 03_parser.py
"""
import array
import os
import sys
import pickle

from settings import *
from utils import *
from lab.utils import time_stamp, caller_file_and_line, eprint
from lab.gene.anno import genic_region_list, get_repr_anno, parse_anno_val
from lab.genome import genome

genome.set_genome(GENOME_FILE_PATH)


def main():
    # qsub settings
    script = os.path.abspath(__file__)
    queue = 'workq'
    is_test = False

    job_name = 'Minu.Get.miRNA.TS.Size'
    log_dir = '%s/log/%s/%s' % (PROJECT_DIR, job_name, time_stamp())

    if not is_test:
        os.makedirs(log_dir, exist_ok=True)

    # param settings
    chroms = chrom_list()
    cons_score_cutoff = 2.0
    step = 1

    # path settings
    # %s in 'phylop_path_format': a chromosome ID
    phylop_path_format = '/extdata6/Minwoo/data/phyloP/{0}/100way-data/%s.phyloP100way.dat'.format(GENOME_VER)
    mirna_list_path = '%s/data/cons_mature_human_list.txt' % PROJECT_DIR
    target_site_dir = '%s/results/target-sites/%s' % (PROJECT_DIR, GENOME_VER)
    peak_data_dir = '%s/peak-data' % target_site_dir
    peak_size_dir = '%s/peak-size/phyloP-%.1f' % (target_site_dir, cons_score_cutoff)  # a result directory
    all_chr_peak_size_dir = '%s/all' % peak_size_dir

    if not os.path.isdir(peak_data_dir):
        eprint('[ERROR] in %s' % caller_file_and_line())
        sys.exit('\t\'%s\' does not exist. Run 03_parser.py.' % peak_data_dir)

    with open(mirna_list_path, 'r') as mirna_list_file:
        mirnas = mirna_list_file.read().splitlines()

    mirnas.append('high_cons_mir')  # merged miRNA target sites

    if step == 1:  # get peak sizes for binding sites of each RBP on each chromosome
        for mirna in mirnas:
            cmd = ''
            cmd_cnt = 0
            cmd_idx = 1

            for chrom in chroms:
                # miRNA and chromosome-specific path settings
                chr_peak_data_dir = '%s/%s' % (peak_data_dir, chrom)
                chr_peak_size_dir = '%s/%s' % (peak_size_dir, chrom)
                os.makedirs(chr_peak_size_dir, exist_ok=True)

                peak_data_path = '%s/%s.dat' % (chr_peak_data_dir, mirna)
                peak_size_path = '%s/%s.txt' % (chr_peak_size_dir, mirna)
                cons_score_path = phylop_path_format % chrom

                cmd += '%s make_peak_size_stats %s %s %s %s %s %.1f;' % \
                       (script, peak_size_path, peak_data_path, chrom, 'True', cons_score_path, cons_score_cutoff)
                cmd_cnt += 1

                if cmd_cnt == 4:
                    if is_test:
                        print(cmd)
                    else:
                        one_job_name = '%s.%s.%s' % (job_name, mirna, cmd_idx)
                        one_log_path = '%s/%s.txt' % (log_dir, one_job_name)
                        os.system('echo "%s" | qsub -j oe -o %s -q %s -N %s' % (cmd, one_log_path, queue, one_job_name))

                    # reset
                    cmd = ''
                    cmd_cnt = 0
                    cmd_idx += 1

    if step == 2:  # combine results from all chromosomes
        job_name = 'Minu.Combine.All.Chr.Peak.Size'
        log_dir = '%s/log/%s/%s' % (PROJECT_DIR, job_name, time_stamp())

        if not is_test:
            os.makedirs(log_dir, exist_ok=True)

        os.makedirs(all_chr_peak_size_dir, exist_ok=True)  # make a result directory for this step

        for mirna in mirnas:
            cmd = '%s combine_all_chr_stats %s %s %s' % (script, all_chr_peak_size_dir, peak_size_dir, mirna)

            if is_test:
                print(cmd)
            else:
                one_job_name = '%s.%s' % (job_name, mirna)
                one_log_path = '%s/%s.txt' % (log_dir, one_job_name)
                os.system('echo "%s" | qsub -j oe -o %s -q %s -N %s' % (cmd, one_log_path, queue, one_job_name))

    if step == 3:  # concat all peak sizes of all miRNA target sites for each gene-based annotation
        job_name = 'Minu.Concat.Peak.Size.by.Anno'
        log_dir = '%s/log/%s/%s' % (PROJECT_DIR, job_name, time_stamp())

        if not is_test:
            os.makedirs(log_dir, exist_ok=True)

        concat_peak_size_dir = '%s/by-anno' % peak_size_dir
        os.makedirs(concat_peak_size_dir, exist_ok=True)

        cmd = '%s concat_peak_size_stats %s %s %s' % \
              (script, concat_peak_size_dir, all_chr_peak_size_dir, mirna_list_path)

        if is_test:
            print(cmd)
        else:
            one_job_name = job_name
            one_log_path = '%s/%s.txt' % (log_dir, one_job_name)
            os.system('echo "%s" | qsub -j oe -o %s -q %s -N %s' % (cmd, one_log_path, queue, one_job_name))


def make_peak_size_stats(peak_size_path, peak_data_path, chrom, only_repr,
                         cons_score_path=None, cons_score_cutoff=-14.0):
    """
    Get the peak sizes from the peak data and save the result as a file
    :param peak_size_path: a path of the result
    :param peak_data_path: a path of the peak data
    :param chrom: a chromosome ID
    :param only_repr: if True, only consider representative genic region
    :param cons_score_path: a path of a file containing an array of conservation scores in the same chromosome
    :param cons_score_cutoff: a float (-14.0 is the minimum)
    """
    eprint('[LOG] Get the peak size from the peak data')
    only_repr = eval(only_repr)
    cons_score_cutoff = float(cons_score_cutoff)
    genic_regions = genic_region_list()
    chr_size = genome.get_chr_size(chrom)

    if cons_score_path is None:
        chr_cons_scores = array.array('f', [cons_score_cutoff] * chr_size)
    else:
        with open(cons_score_path, 'rb') as infile:
            chr_cons_scores = array.array('f', [])
            chr_cons_scores.fromfile(infile, chr_size)

    if not os.path.isfile(peak_data_path):
        eprint('[LOG] --- \'%s\' does not exist.' % peak_data_path)
        return

    with open(peak_data_path, 'rb') as peak_data_file:
        peaks = pickle.load(peak_data_file)

    peak_size_dict = {genic_region: 0 for genic_region in genic_regions}
    peak_size_dict['UTR'] = 0  # UTR: all UTR (5'UTR or 3'UTR)
    peak_size_dict['all'] = 0  # all genic regions

    for peak in peaks:
        peak_cons_scores = chr_cons_scores[peak.start:peak.end]
        anno_vals = peak.get_anno_vals()

        for i, anno_val in enumerate(anno_vals):
            if peak_cons_scores[i] >= cons_score_cutoff:
                peak_size_dict['all'] += 1

                if only_repr:
                    both_utr, repr_genic_region = get_repr_anno(anno_val)

                    if both_utr:
                        peak_size_dict['5UTR'] += 1
                        peak_size_dict['3UTR'] += 1
                    else:
                        peak_size_dict[repr_genic_region] += 1

                    if repr_genic_region.endswith('UTR'):
                        peak_size_dict['UTR'] += 1
                else:
                    anno_dict = parse_anno_val(anno_val)

                    for genic_region in genic_regions:
                        if anno_dict[genic_region]:
                            peak_size_dict[genic_region] += 1

                            if genic_region.endswith('UTR'):
                                peak_size_dict['UTR'] += 1

    eprint('[LOG] Save the result')
    genic_regions += ['UTR', 'all']

    with open(peak_size_path, 'w') as peak_size_file:
        for genic_region in genic_regions:
            print(genic_region, peak_size_dict[genic_region], sep='\t', file=peak_size_file)


def combine_all_chr_stats(result_dir, peak_size_dir, mirna):
    """
    Combine the peak size stats of all chromosomes and save the result
    :param result_dir: a result directory
    :param peak_size_dir: a directory files for peak size stats are stored
    :param mirna: a name of RBP
    """
    eprint('[LOG] Combine all peak size statistics of %s' % mirna)
    genic_regions = genic_region_list()
    genic_regions += ['UTR', 'all']
    chroms = chrom_list()

    peak_size_dict = {genic_region: 0 for genic_region in genic_regions}

    for chrom in chroms:
        peak_size_path = '%s/%s/%s.txt' % (peak_size_dir, chrom, mirna)

        if not os.path.isfile(peak_size_path):
            eprint('[LOG] --- \'%s\' does not exist.' % peak_size_path)
            continue

        with open(peak_size_path, 'r') as peak_size_file:
            for line in peak_size_file.readlines():
                fields = line.strip().split('\t')
                genic_region = fields[0]
                peak_size = int(fields[1])
                peak_size_dict[genic_region] += peak_size

    eprint('[LOG] Save the result')
    with open('%s/%s.txt' % (result_dir, mirna), 'w') as peak_size_file:
        for genic_region in genic_regions:
            print(genic_region, peak_size_dict[genic_region], sep='\t', file=peak_size_file)


def concat_peak_size_stats(result_dir, peak_size_dir, mirna_list_path):
    """
     Concat all peak sizes from all RBPs save them as a file for each genic region
    :param result_dir: a result directory
    :param peak_size_dir: a directory peak size files for all RBPs are saved
    :param mirna_list_path: a path of a file containing a mirna list
    """
    eprint('[LOG] Concatenate all peak sizes for each genic region')
    with open(mirna_list_path, 'r') as mirna_list_file:
        mirnas = mirna_list_file.read().splitlines()

    mirnas.append('all_merge')

    genic_regions = genic_region_list()
    genic_regions += ['UTR', 'all']

    result_file_dict = {genic_region: open('%s/%s.txt' % (result_dir, genic_region), 'w')
                        for genic_region in genic_regions}

    for mirna in mirnas:
        peak_size_file_path = '%s/%s.txt' % (peak_size_dir, mirna)

        with open(peak_size_file_path, 'r') as peak_size_file:
            for line in peak_size_file.readlines():
                fields = line.strip().split('\t')
                genic_region = fields[0]
                peak_size = int(fields[1])

                print(mirna, peak_size, sep='\t', file=result_file_dict[genic_region])

    for genic_region in genic_regions:
        result_file_dict[genic_region].close()


if __name__ == '__main__':
    if len(sys.argv) == 1:
        main()
    else:
        function_name = sys.argv[1]
        function_parameters = sys.argv[2:]

        if function_name in locals().keys():
            locals()[function_name](*function_parameters)
        else:
            sys.exit('ERROR: function_name=%s, parameters=%s' % (function_name, function_parameters))
