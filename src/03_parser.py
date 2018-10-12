#!/opt/ohpc/pub/bin/anaconda3/envs/minu/bin/python3
"""
Parsing NarrowPeak files of miRNA target sites

Prerequisite
    1. Run 03_anno_genome.py from the project 'repr-gene'
    2. Run 02_target_site.py
"""
from settings import *
from utils import *
from lab.utils import time_stamp, caller_file_and_line, eprint
from lab.bed.mutpeak import MutPeak
from lab.genome import genome

import array
import os
import sys
import pickle

genome.set_genome(GENOME_FILE_PATH)


def main():
    """
    Bootstrap
    """
    # settings for a job scheduler
    script = os.path.abspath(__file__)
    queue = 'workq'
    is_test = False

    job_name = 'Minu.Parse.miRNA.Target'
    log_dir = '%s/log/%s/%s' % (PROJECT_DIR, job_name, time_stamp())

    if not is_test:
        os.makedirs(log_dir, exist_ok=True)

    # path settings
    anno_dir = '/extdata6/Minwoo/projects/repr-gene/results/genome-anno/%s' % GENOME_VER
    mirna_list_path = '%s/data/cons_mature_human_list.txt' % PROJECT_DIR

    target_site_dir = '%s/results/target-sites/%s' % (PROJECT_DIR, GENOME_VER)
    mirna_ts_bed_dir = '%s/bed' % target_site_dir
    mirna_ts_data_dir = '%s/peak-data' % target_site_dir  # a result directory
    
    if not os.path.isdir(anno_dir):
        eprint('[ERROR] in %s' % caller_file_and_line())
        sys.exit('\t\'%s\' does not exist. Run repr-gene/03_anno_genome.py first.' % anno_dir)
    
    if not os.path.isdir(mirna_ts_bed_dir):
        eprint('[ERROR] in %s' % caller_file_and_line())
        sys.exit('\t\'%s\' do not exist. Run 02_target_site.py first.' % mirna_ts_bed_dir)

    os.makedirs(mirna_ts_data_dir, exist_ok=True)
    
    with open(mirna_list_path, 'r') as mirna_list_file:
        mirnas = mirna_list_file.read().splitlines()

    chroms = chrom_list()

    # Target sites of an individual miRNA
    for mirna in mirnas:
        mirna_ts_bed_path = '%s/%s.bed' % (mirna_ts_bed_dir, mirna)
        assert os.path.isfile(mirna_ts_bed_path)

        cmd = ''
        cmd_cnt = 0
        cmd_idx = 1

        for chrom in chroms:
            cmd += '%s parse_peak_file %s %s %s %s;' % \
                   (script, mirna_ts_data_dir, mirna_ts_bed_path, anno_dir, chrom)
            cmd_cnt += 1

            if cmd_cnt == 4:  # one job for 4 chromosomes
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

    # Parsing merged miRNA target sites
    merged_mirna_bed_path = '%s/high_cons_mir.bed' % mirna_ts_bed_dir
    assert os.path.isfile(merged_mirna_bed_path)

    for chrom in chroms:
        cmd = '%s parse_peak_file %s %s %s %s' % \
              (script, mirna_ts_data_dir, merged_mirna_bed_path, anno_dir, chrom)

        if is_test:
            print(cmd)
        else:
            one_job_name = '%s.%s.%s' % (job_name, 'All', chrom)
            one_log_path = '%s/%s.txt' % (log_dir, one_job_name)
            os.system('echo "%s" | qsub -j oe -o %s -q %s -N %s' % (cmd, one_log_path, queue, job_name))


def parse_peak_file(result_dir, bed_file_path, anno_dir, chrom):
    """
    Make MutPeak objects by parsing the BED file and save the objects
    :param result_dir: a directory parsed peak will be saved.
    :param bed_file_path: a file that has a narrow peak bed file format
    :param anno_dir: a directory genome annotation data is saved.
    :param chrom: a chromosome ID
    """
    eprint('[LOG] Parse the peaks from \'%s\'' % bed_file_path)
    eprint('[LOG] ---- Chromosome ID: %s' % chrom)
    peaks = []

    # Parse the NarrowPeak bed file and make 'MutPeak' objects
    with open(bed_file_path, 'r') as bed_file:
        for line in bed_file:
            if line.startswith('%s\t' % chrom):
                peak = MutPeak()
                peak.parse_peak_entry(line.strip())
                peaks.append(peak)

    # parse the genic regions on the peaks
    chr_size = genome.get_chr_size(chrom)

    # TODO: change the algorithm of 'anno_peak' to reduce I/O.
    for peak in peaks:
        anno_peak(anno_dir, peak, chr_size)

    # saved the result
    chr_result_dir = '%s/%s' % (result_dir, chrom)
    os.makedirs(chr_result_dir, exist_ok=True)

    result_filename = os.path.basename(bed_file_path).replace('.bed', '.dat')
    result_file_path = '%s/%s' % (chr_result_dir, result_filename)

    with open(result_file_path, 'wb') as result_file:
        pickle.dump(peaks, result_file)


def anno_peak(anno_dir, peak, chr_size):
    """
    Read data about the genic regions of the peak and save the data on the object for the peak

    :param anno_dir: a directory genome annotation data is saved.
    :param peak: a 'RBPPeak' object
    :param chr_size: a size of chromosome ID same with the 'chrom' of the peak
                     this parameter is necessary to know the end position of the last chromosomal fragment.
    """
    assert peak.__class__.__name__ == 'MutPeak'

    peak_start, peak_end = peak.get_position()
    assert peak_end <= chr_size

    if peak.strand == '+':
        peak_strand = 'top'
    else:  # '-'
        peak_strand = 'btm'

    bin_size = 1000000
    start_digit = int(peak_start / bin_size)  # Most significant digit
    end_digit = int(peak_end / bin_size)
    digit_diff = end_digit - start_digit

    # Split the peak with consideration of the bin sizes of the genome annotation files
    peak_frags = []  # A list of tuple (peak_frag_start, peak_frag_end)
    peak_frag_start = peak_start
    peak_frag_end = start_digit * bin_size

    for i in range(digit_diff):
        peak_frag_end += bin_size
        peak_frags.append((peak_frag_start, peak_frag_end))
        peak_frag_start = peak_frag_end

    peak_frag_end = peak_end
    peak_frags.append((peak_frag_start, peak_frag_end))

    # Read the array of genic regions and concatenates them
    chr_anno_dir = '%s/%s' % (anno_dir, peak.chrom)

    if not os.path.isdir(anno_dir):
        eprint('[ERROR] in %s: %s does not exist.' % (caller_file_and_line(), anno_dir))
        sys.exit()

    total_frag_len = 0
    anno_val_arr = []

    for peak_frag_start, peak_frag_end in peak_frags:
        # Chromosomal bins: 0-1000000, 1000000-2000000, ...
        chr_bin_start = int(peak_frag_start / bin_size) * bin_size
        chr_bin_end = chr_bin_start + bin_size

        if chr_bin_end > chr_size:
            chr_bin_end = chr_size

        if peak_frag_end < chr_bin_end:
            frag_len = peak_frag_end - peak_frag_start
        else:
            frag_len = chr_bin_end - peak_frag_start

        total_frag_len += frag_len

        frag_anno_val_arr = array.array('i', [])
        arr_item_size = frag_anno_val_arr.itemsize

        anno_file = open('%s/%d_%d_%s.dat' % (chr_anno_dir, chr_bin_start, chr_bin_end, peak_strand), 'rb')
        anno_file.seek((peak_frag_start - chr_bin_start) * arr_item_size)  # relative bin peak_start position
        frag_anno_val_arr.fromfile(anno_file, frag_len)
        anno_file.close()

        anno_val_arr += frag_anno_val_arr

    peak.gene_based_anno(list(anno_val_arr))
    assert total_frag_len == (peak_end - peak_start)


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
