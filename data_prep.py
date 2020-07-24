from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
import numpy as np
from utils import *

transcript_file = np.genfromtxt('./data/GENCODE_v34_hg38_comprehensive', usecols=(1, 2, 3, 4, 5, 9, 10), skip_header=1, dtype='str')
DE_tr = np.genfromtxt('./lists/limma_DE_AML_RAN', usecols=(1,), skip_header=1, dtype='str')

print(DE_tr)

path_to_file = './data/hg38.fa'
chr_names = ['chr'+str(x) for x in [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,'X','Y']]

chr_dict = {}

with open(path_to_file, mode='r') as handle:

    for record in SeqIO.parse(handle, 'fasta'):

        identifier = record.id
        description = record.description
        sequence = record.seq

        if identifier in chr_names:
            chr_dict[identifier] = sequence

print(chr_dict.keys())

transcripts = []

for t in DE_tr:
    found = 0
    # explicitly checking transcript_name
    for row in transcript_file:
        if row[0][:15]==t[:15]:
            found = 1
            s = chr_dict[row[1]][int(row[3]): int(row[4])].upper()
            print(len(s))
            # adding the transcripts of the sense strand: whole transcript + flanks + zero-padded, labels + zero-padded
            if row[2] == '+':
                # extract the transcript sequence with 1k flanks
                if 'N' not in s:
                    # padding labels here
                    es, ee = row[5].split(',')[:-1], row[6].split(',')[:-1]
                    # decrease the pad length from both sides because the context-1 and context+sequence+1 sites are
                    # donor and acceptor, respectively
                    s = make_mRNA(s, es, ee)
                    print(len(s))
                    transcripts.append(s)
            # adding the transcripts of the antisense strand
            if row[2] == '-':
                if 'N' not in s:
                    # padding labels here
                    es, ee = row[5].split(',')[:-1], row[6].split(',')[:-1]
                    # decrease the pad length from both sides because the context-1 and context+sequence+1 sites are
                    # donor and acceptor, respectively
                    s = make_mRNA(s, es, ee)
                    s = ''.join([complementary(x) for x in s])
                    print(len(s))
                    transcripts.append(s)
    if found==0:
        print(t, 'not found')


handle = open('./data/DE_AML_transcripts', "w")
records = []

for i in range(len(transcripts)):
    records.append(SeqRecord(Seq(transcripts[i], generic_dna), DE_tr[i]))

SeqIO.write(records, handle, "fasta")