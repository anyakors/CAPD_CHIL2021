from utils import *
import h5py
import argparse
import time
import random
import numpy as np


def find_indx(all_sources, keyword):
    ind = []
    if isinstance(keyword, list):
        for k in range(len(all_sources)):
            if all([x.lower() in all_sources[k].decode('utf-8').lower() for x in keyword]):
                ind.append(k)
    else:
        for k in range(len(all_sources)):
            if keyword.lower() in all_sources[k].decode('utf-8').lower():
                ind.append(k)
    print('Sources with {}: {}'.format(str(keyword), len(ind)))
    return ind

# link for file download:
# javascript:downloadFile('https://s3.amazonaws.com/mssm-seq-matrix/human_transcript_v8.h5','human_transcript.h5','8')

parser = argparse.ArgumentParser()
parser.add_argument("-h5", "--h5file", help="ARCHS4 human h5 table directory",
                    type=str)
parser.add_argument("-s", "--savedir", help="Path to save countfiles",
                    type=str)
parser.add_argument("-m", "--max_samples", help="Max samples per tissue",
                    type=int)
parser.add_argument("-n", "--min_samples", help="Min samples per tissue",
                    type=int)

# ! New args for shuffling samples- start
parser.add_argument("--seed", default=42, help="Random seed for reproducibility",
                    type=int)
parser.add_argument("--shuffle", action="store_true", help="Whether to shuffle samples")
# ! New args for shuffling samples- end

args = parser.parse_args()

# ! Set random seed for reproducibility- start
random.seed(args.seed)
np.random.seed(args.seed)
# ! Set random seed for reproducibility- end

f = h5py.File(args.h5file, 'r')
# list of all sample sources and transcripts
all_sources = list(f['meta']['Sample_source_name_ch1'][()])
all_transcripts = list(f['meta']['transcripts'][()])

print('No of sources (samples): {}'.format(len(all_sources)))
print('No of transcripts: {}'.format(len(all_transcripts)))

# ADP ADR BLD BRN BRS CLN HRT KDN LNG LMP OVR PRS SKM TST THR AML
tissues = ['adipose', 'adrenal', 'blood', 'brain', 'breast', 'colon',
           'heart', 'kidney', 'liver', 'lung', 'lymph', 'ovary',
           'prostate', ['skeletal', 'muscle'], 'testis', 'thyroid']

start = time.time()

indx = {}
for tissue in tissues:
    if not isinstance(tissue, list):
        indx[tissue] = find_indx(all_sources, tissue)
    else:
        # the keyword is only the first word
        indx[tissue[0]] = find_indx(all_sources, tissue)

# looking separately for AML samples as there are at least two ways to write it down
indx['AML'] = find_indx(all_sources, ['acute', 'myeloid', 'leukemia'])
indx['AML'].extend(find_indx(all_sources, 'AML'))

# ! Shuffle sample indices- start
if args.shuffle:
    for tissue in indx.keys():
        random.shuffle(indx[tissue])
# ! Shuffle sample indices- end

# create a table with all the sample expression arrays (AML):
# rows - samples, cols - transcripts, elements - No of reads, not normalised
column_names = ['samples']
column_names.extend([x.decode('utf-8') for x in all_transcripts])

data = []
min_samples = args.min_samples
max_samples = args.max_samples

if min_samples:
    for key in indx.keys():
        if len(indx[key])<min_samples:
            print('{} is not represented enough; excluded'.format(key))
            indx.pop(key)

for tissue in indx.keys():
    counts_to_csv(f, column_names, indx, tissue, args.savedir, max_samples)

time_to_human(time.time() - start)