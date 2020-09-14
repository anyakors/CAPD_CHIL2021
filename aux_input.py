from utils import *

import numpy as np
import pandas as pd
import os
import json
import argparse


parser = argparse.ArgumentParser()
parser.add_argument("-g", "--gencode_list", help="GENCODE comprehensive list with exon info",
                    type=str)
parser.add_argument("-rb", "--rbp_list", help="List of RBP transcript IDs",
                    type=str)
parser.add_argument("-rm", "--rmod_list", help="List of RNA-modif protein transcript IDs",
                    type=str)
parser.add_argument("-i", "--input", help="Count tables path",
                    type=str)
parser.add_argument("-t", "--tr_dict", help="Transcript length dict",
                    type=str)
parser.add_argument("-o", "--output", help="Savepath for aux inputs",
                    type=str)
args = parser.parse_args()


# 0 - tr. name, 1 - chr, 2 - strand, 3,4 - tr. start/end, 5,6 - exon starts/ends, 7 - gene name
transcript_file = np.genfromtxt(args.gencode_list, usecols=(1, 2, 3, 4, 5, 9, 10, 12), skip_header=1, dtype='str')

l1 = pd.read_csv(args.rbp_list, sep=',', header=0, index_col=0)
l2 = pd.read_csv(args.rmod_list, sep=',', header=0, index_col=0)
l = l1.append(l2)
print('Aux transcripts IDs:', len(l))

dirs = os.listdir(args.input)
init = 0

print('Importing count tables from the input path...')

for file in dirs:
    if 'csv' in file:
        if init == 1:
            counts = counts.append(pd.read_csv(os.path.join(args.input+file), sep=',', header=0, index_col=0, low_memory=False))
        else:
            counts = pd.read_csv(os.path.join(args.input+file), sep=',', header=0, index_col=0, low_memory=False)
            init = 1

# normalize to cpm
#F = counts.sum(axis = 1)/10**6
#counts = counts.divide(F, axis='index')

print('Normalizing counts to rpkm...')

to_drop = []
for tr in list(counts.columns):
    if tr[:15] not in list(l['transcript_ID']):
        to_drop.append(tr)

counts = counts.drop(columns=to_drop)

tr_dict = {}
with open(args.tr_dict) as json_file:
    json_list = list(json_file)

for json_str in json_list:
    result = json.loads(json_str)
    tr_dict.update(result)

nf = 0
for tr in counts.columns:
    if tr in tr_dict.keys():
        counts[tr] = counts[tr]*1000/tr_dict[tr]
    else:
        counts = counts.drop(columns=tr)
        nf += 1

print("Transcripts found and retained in the aux counts table:", len(counts.columns))

if not os.path.exists(args.output):
    os.makedirs(args.output)

counts.to_csv(os.path.join(args.output, 'aux_RBP_RNAmod.csv'))