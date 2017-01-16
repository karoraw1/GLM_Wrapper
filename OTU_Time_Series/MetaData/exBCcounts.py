"""
These show the first lines of the files from which we need to extract the barcode information: 
1.@HWUSI-EAS1563_0022:1:1:1171:1779#TNTNTATANNNNNN/1
2.@MISEQ578:1:1101:15698:1856#NNNNNNNN/1
3.@MISEQ:1:1101:11244:2795#TGGGACCT/1
4.@MISEQ:1:1101:12371:2841#AAACAAGA/1
5.@MISEQ:1:1101:11435:2210#CTCTCCCC/1
6.@MISEQ:1:1101:11165:2419#ATATCGCC/1
7.@MISEQ:1:1101:13137:1805#AGGTTTAT/1
8.@HWI-M01799:1:1101:16019:1402#CGAATATT/1
9.@FCD19KAACXX:6:1101:1703:2235#CCGACAAA/1
"""

import os, sys
from collections import Counter

base_dir = '/scratch/groups/sprehei1/data/Mystic_lake_time_series_dir/raw_data'
ind_files = ['100714_B/100702_MystikLake_amaterna_Alm_L1_1_sequence.txt', 
             '121114Alm_dir/121114Alm_D12-4491_1_sequence.fastq',
             '130719Alm_dir/130719Alm_D13-3437_1_sequence.fastq',
             '130823Alm_dir/130823Alm_D13-4258_1_sequence.fastq', 
             '131001Alm_dir/131001Alm_D13-4961_phiX_best_1.fastq', 
             '131011Alm_dir/131011Alm_D13-5267_phiX_best_1.fastq',
             '131114Alm_dir/131114Alm_D13-6069_1_sequence.fastq',
             '131126Alm_dir/131126Alm_D13-6337_1_sequence.fastq', 
             'BGI_092012_dir/newsplit.Sarah_ML.1']

path_list = [os.path.join(base_dir, i) for i in ind_files]

def collectBCcounts(a_path):
    bc_list = []
    print "Mapping {}".format(os.path.basename(a_path))
    with open(a_path, 'r') as f:
        for line in f:
            this_l = line.strip()
            if line[0] == "@":
                this_bc = line.split("#")[1].split("/")[0]
                bc_len = len(this_bc)
                if bc_len != 14 and bc_len !=8:
                    print this_bc, bc_len
                    sys.exit("abnormal length detected")
                bc_list.append(this_bc)
    bc_dict = Counter(bc_list)
    return bc_dict 

import cPickle as pickle

print os.path.basename(path_list[int(sys.argv[1])])

bc_dict = collectBCcounts(path_list[int(sys.argv[1])])

out_name = "barcode_counts" + str(sys.argv[1]) + ".p"

pickle.dump(bc_dict, open( out_name, "wb" ))

print out_name, "saved"
