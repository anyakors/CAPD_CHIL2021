from Bio import SeqIO
from utils import *

import numpy as np
import pandas as pd

import time
import os
import argparse


parser = argparse.ArgumentParser()
parser.add_argument("-g", "--gencode_list", help="GENCODE comprehensive list with exon info",
                    type=str)
parser.add_argument("-hg", "--hg38", help="Path to hg38.fa",
                    type=str)
parser.add_argument("-r", "--genome_region", help="Genome region to extract genes from; either 'all' or a specific chr: 'chr22'",
                    type=str)
parser.add_argument("-d", "--genes_to_keep", help="Genes to keep in the labels",
                    type=str)
#parser.add_argument("-t", "--test_ratio", help="Ratio (0 to 1) of data to test model on",
#                    type=float)
parser.add_argument("-gl", "--length_limit", help="Gene length limit: above this length -> testing batch",
                    type=int)
parser.add_argument("-c", "--context", help="Flanking ends lengths (context), on each side",
                    type=int)
parser.add_argument("-t", "--test", help="Whether it's a test chr or no",
                    type=int)
parser.add_argument("-i", "--input", help="Count tables path",
                    type=str)
parser.add_argument("-ot", "--out_train", help="Savepath for main inputs to train",
                    type=str)
parser.add_argument("-os", "--out_test", help="Savepath for main inputs to test",
                    type=str)
args = parser.parse_args()


if args.genome_region=='all':
    region = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7',
              'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13',
              'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19',
              'chr20', 'chr21', 'chr22', 'chrX', 'chrY']
else:
    region = args.genome_region.split(",")

print("====================================================================")
print("Main inputs generation; genome region defined:", region)
print("====================================================================")

transcript_file = np.genfromtxt(args.gencode_list, usecols=(1, 2, 3, 4, 5, 9, 10, 12), skip_header=1, dtype='str')

print("Processing", region)
gene_dict = make_gene_dict(transcript_file, region)

to_keep = pd.read_csv(args.genes_to_keep, sep=',', header=0, index_col=0)
for key in list(gene_dict.keys()):
    if key not in list(to_keep['gene_symbol']):
        gene_dict.pop(key)

print("Omitted genes with missing exon info (GENCODEv26), left:", len(gene_dict.keys()))

if args.test:
    gene_dict_test = gene_dict
    print("{} genes selected for testing".format(len(gene_dict_test)))
else:
    gene_dict_train, gene_dict_test = train_test_sep_limit(gene_dict, args.length_limit)
    print("{} genes selected for training, {} (longer than threshold) for testing".format(len(gene_dict_train), len(gene_dict_test)))


hg38 = {}

with open(args.hg38, mode='r') as handle:
    for record in SeqIO.parse(handle, 'fasta'):
        identifier = record.id
        description = record.description
        sequence = record.seq

        hg38[identifier] = sequence


# construct sample library from count file and save into sep jsonl files
dirs = os.listdir(args.input)

start_time = time.time()

for file in dirs:
    if 'csv' in file:
        counts = pd.read_csv(os.path.join(args.input, file), sep=',', header=0, index_col=0, low_memory=False)
        print("+++++++++++++++++++++++++++++++++++")
        print('Generating inputs from', file, '; shape:', np.shape(counts))
        print("+++++++++++++++++++++++++++++++++++")

        if not args.test:
            save_labels_jsonl(counts, gene_dict_train, hg38, args.out_train, args.context, region)
            save_labels_jsonl(counts, gene_dict_test, hg38, args.out_test, args.context, region,
                                      long_from_train=True)
        else:
            save_labels_jsonl(counts, gene_dict_test, hg38, args.out_test, args.context, region)

        print('Saved .jsonl dicts for training/testing with {} label; one file = one sample'.format(file[:3]))


time_to_human(time.time() - start_time)