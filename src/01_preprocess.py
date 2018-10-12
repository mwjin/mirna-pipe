#!/opt/ohpc/pub/bin/anaconda3/envs/minu/bin/python3
"""
Extract only human miRNAs from mirBase's miRNAs and make the list of miRNAs
"""
from settings import *
from lab.utils import eprint

import os
import re

# input path settings
mirbase_file_path = '/extdata6/Minwoo/data/miRNA/mature.fa'
targetscan_fam_info_path = '/extdata6/Minwoo/data/miRNA/miR_Family_Info.txt'

# result path settings
pre_proc_dir = '%s/data' % PROJECT_DIR
os.makedirs(pre_proc_dir, exist_ok=True)

human_mir_file_path = '%s/mature_human.fa' % pre_proc_dir
human_mir_list_path = '%s/mature_human_list.txt' % pre_proc_dir
cons_human_mir_list_path = '%s/cons_mature_human_list.txt' % pre_proc_dir  # conserved miRNAs

# Parse the mature.fa
eprint('[LOG] Extract human miRNAs from mirBase')
mirna_info = []  # element: (header, seq)

with open(mirbase_file_path, 'r') as mirbase_file:
    while True:
        header = mirbase_file.readline()

        if not header:  # EOF
            break

        seq = mirbase_file.readline()
        mirna_info.append((header.strip(), seq.strip()))

eprint('[LOG] The number of miRNAs from mirBase: %d' % len(mirna_info))

# Collect only human miRNAs with an annotation number less than 1000
human_mir_info_dict = {}  # key: miRNA name, value: (header, seq)
human_mir_acc_dict = {}  # key: miRNA name, value: mirBase accession ID

for mir_header, mir_seq in mirna_info:
    fields = mir_header.split()
    mir_name = fields[0][1:]
    mir_acc_id = fields[1]

    if mir_name.startswith('hsa'):  # human
        name_fields = mir_name.split('-')
        anno_id = name_fields[2]

        num_match = re.match('[0-9]+', anno_id)
        assert num_match is not None
        anno_num = int(num_match.group())

        if anno_num < 1000:
            human_mir_info_dict[mir_name] = (mir_header, mir_seq)
            human_mir_acc_dict[mir_name] = mir_acc_id
    else:
        continue

eprint('[LOG] The number of human miRNAs from mirBase: %d' % len(human_mir_info_dict.keys()))

# filter the miRNAs by their conservations
mir_acc_to_cons = {}  # key: accID, value: 2 (Highly conserved), 1 (Conserved), 0 (Poorly), -1 (Poorly, Misannotated)

with open(targetscan_fam_info_path, 'r') as mirna_fam_info_file:
    mirna_fam_info_file.readline()  # remove a header

    for line in mirna_fam_info_file.readlines():
        fields = line.strip().split('\t')

        if len(fields) == 6:  # no mirBase accession ID
            continue

        cons_val = int(fields[5])
        mir_acc_id = fields[6]
        mir_acc_to_cons[mir_acc_id] = cons_val

cons_human_mirnas = []

for mir_name in human_mir_info_dict:
    mir_acc_id = human_mir_acc_dict[mir_name]
    cons_val = mir_acc_to_cons.get(mir_acc_id)

    if cons_val is None:
        continue
    elif cons_val >= 2:
        cons_human_mirnas.append(mir_name)

eprint('[LOG] The number of conserved human miRNAs from mirBase: %d' % len(cons_human_mirnas))

# save the result
eprint('[LOG] Save the results')

with open(human_mir_file_path, 'w') as human_mir_file, open(human_mir_list_path, 'w') as human_mir_list_file:
    for mir_name in human_mir_info_dict:
        mir_header, mir_seq = human_mir_info_dict[mir_name]
        print(mir_name, file=human_mir_list_file)
        print(mir_header, file=human_mir_file)
        print(mir_seq, file=human_mir_file)

with open(cons_human_mir_list_path, 'w') as cons_human_mir_file:
    for mir_name in cons_human_mirnas:
        print(mir_name, file=cons_human_mir_file)
