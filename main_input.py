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
parser.add_argument("-t", "--test_ratio", help="Ratio (0 to 1) of data to test model on",
                    type=float)
parser.add_argument("-c", "--context", help="Flanking ends lengths (context), on each side",
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

print("Genome region defined:", region)

gene_dict = {}

# Import GENCODE v34 & construct gene dict
transcript_file = np.genfromtxt(args.gencode_list, usecols=(1, 2, 3, 4, 5, 9, 10, 12), skip_header=1, dtype='str')

start_time = time.time()

# each key is a gene
for row in transcript_file:
    if row[1] in region:
        if row[7] not in gene_dict.keys():
            gene_dict[row[7]] = {}
            gene_dict[row[7]]['chr'] = row[1]
            gene_dict[row[7]]['strand'] = row[2]
            gene_dict[row[7]][row[0]] = {}
            gene_dict[row[7]][row[0]]['starts'] = [int(x) for x in row[5].split(',')[:-1]]
            gene_dict[row[7]][row[0]]['ends'] = [int(x) for x in row[6].split(',')[:-1]]
            gene_dict[row[7]][row[0]]['length'] = mRNA_len(gene_dict[row[7]][row[0]]['starts'],
                                                                gene_dict[row[7]][row[0]]['ends'])
            gene_dict[row[7]]['global_start'] = int(row[3])
            gene_dict[row[7]]['global_end'] = int(row[4])
        else:
            gene_dict[row[7]][row[0]] = {}
            gene_dict[row[7]][row[0]]['starts'] = [int(x) for x in row[5].split(',')[:-1]]
            gene_dict[row[7]][row[0]]['ends'] = [int(x) for x in row[6].split(',')[:-1]]
            gene_dict[row[7]][row[0]]['length'] = mRNA_len(gene_dict[row[7]][row[0]]['starts'],
                                                                gene_dict[row[7]][row[0]]['ends'])
            if int(row[3]) < gene_dict[row[7]]['global_start']:
                gene_dict[row[7]]['global_start'] = int(row[3])
            if int(row[4]) > gene_dict[row[7]]['global_end']:
                gene_dict[row[7]]['global_end'] = int(row[4])


print("Took {} seconds to construct the gene dict".format(time.time() - start_time))
print("Number of genes in the region of interest:", len(gene_dict.keys()))

to_keep = pd.read_csv(args.genes_to_keep, sep=',', header=0, index_col=0)
for key in list(gene_dict.keys()):
    if key not in list(to_keep['gene_symbol']):
        gene_dict.pop(key)

print("Omitted genes with missing exon info (GENCODEv26), left:", len(gene_dict.keys()))

gene_dict_train, gene_dict_test = train_test_sep(gene_dict, args.test_ratio)

print("{} genes selected for training, {} for testing".format(len(gene_dict_train), len(gene_dict_test)))
# need to download and unpack the genome file into the ./data directory:
# wget --timestamping 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz' -O hg38.fa.gz
# gunzip hg38.fa.gz

hg38 = {}

with open(args.hg38, mode='r') as handle:
    for record in SeqIO.parse(handle, 'fasta'):
        identifier = record.id
        description = record.description
        sequence = record.seq

        hg38[identifier] = sequence


# construct sample library from count file and save into sep jsonl files
dirs = os.listdir(args.input)
print("Input directories:", dirs)

start_time = time.time()

for file in dirs:
    if 'csv' in file:
        counts = pd.read_csv(os.path.join(args.input, file), sep=',', header=0, index_col=0, low_memory=False)
        print(file, '; shape:', np.shape(counts))
        save_labels_jsonl(counts, gene_dict_train, hg38, args.out_train, args.context)
        print('Saved .jsonl dicts for training in {} with {} label; one file = one sample'.format(args.out_train,
                                                                                        file[:3]))
        save_labels_jsonl(counts, gene_dict_test, hg38, args.out_test, args.context)
        print('Saved .jsonl dicts for testing in {} with {} label; one file = one sample'.format(args.out_test,
                                                                                                  file[:3]))

time_to_human(time.time() - start_time)