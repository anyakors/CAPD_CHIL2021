import os
import pandas as pd
import json
import time
import numpy as np


def tissue_to_label(tissue):
    tissue_lbls = {'adipose': 'ADP',
                   'adrenal': 'ADR',
                   'blood': 'BLD',
                   'brain': 'BRN',
                   'breast': 'BRS',
                   'colon': 'CLN',
                   'heart': 'HRT',
                   'kidney': 'KDN',
                   'liver': 'LVR',
                   'lung': 'LNG',
                   'lymph': 'LMP',
                   'ovary': 'OVR',
                   'prostate': 'PRS',
                   'skeletal': 'SKM',
                   'testis': 'TST',
                   'thyroid': 'THR',
                   'AML': 'AML'}
    return tissue_lbls[tissue]


def counts_to_csv(f, column_names, indx, tissue, savepath, max_samples=None):

    if not os.path.exists(savepath):
        os.makedirs(savepath)

    data = []
    label = tissue_to_label(tissue)
    if max_samples:
        for i in range(min(len(indx[tissue]), max_samples)):
            data.append([label+str(i+1)] + list( f['data']['expression'][indx[tissue][i]]))
    else:
        for i in range(len(indx[tissue])):
            data.append([label+str(i+1)] + list( f['data']['expression'][indx[tissue][i]]))

    df = pd.DataFrame(columns=column_names, data=data)
    df = df.set_index('samples')
    df.to_csv(os.path.join(savepath, label+'.csv'))
    print('{} counts saved to .csv in {}'.format(tissue, savepath))

    return


def time_to_human(time):
    hrs = time//3600
    mins = (time - hrs*3600)//60
    secs = time - hrs*3600 - mins*60
    print('Overall time elapsed: {} hrs {} mins {} seconds'.format(int(hrs), int(mins), round(secs)))
    return


def mRNA_len(starts, ends):
    l = 0
    for x in zip(starts, ends):
        l += x[1] - x[0]
    return l


def construct_library(gene_dict, counts):
    library = {}
    start_time = time.time()

    for sample in counts.index:
        library[sample] = {}
        for gene in gene_dict.keys():
            for tr in gene_dict[gene].keys():

                if tr in counts.columns:
                    if gene not in library[sample].keys():
                        library[sample][gene] = {}
                        library[sample][gene] = {}
                        length = int(gene_dict[gene]['global_end']) - int(gene_dict[gene]['global_start']) + 1
                        library[sample][gene]['alabels'] = np.zeros(length)
                        library[sample][gene]['dlabels'] = np.zeros(length)
                        library[sample][gene]['norm_factor'] = 0
                        gs = int(gene_dict[gene]['global_start'])
                        ge = int(gene_dict[gene]['global_end'])
                        l = gene_dict[gene][tr]['length']
                        for s in gene_dict[gene][tr]['starts']:
                            # normalise to tpkm
                            library[sample][gene]['alabels'][int(s) - gs] += float(counts.at[sample, tr] / l)
                        for s in gene_dict[gene][tr]['ends']:
                            # normalise to tpkm
                            library[sample][gene]['dlabels'][int(s) - gs] += float(counts.at[sample, tr] / l)
                        # normalise to tpkm
                        library[sample][gene]['norm_factor'] += float(counts.at[sample, tr] / l)
                    else:
                        gs = int(gene_dict[gene]['global_start'])
                        ge = int(gene_dict[gene]['global_end'])
                        l = gene_dict[gene][tr]['length']
                        for s in gene_dict[gene][tr]['starts']:
                            # normalise to tpkm
                            library[sample][gene]['alabels'][int(s) - gs] += float(counts.at[sample, tr] / l)
                        for s in gene_dict[gene][tr]['ends']:
                            # normalise to tpkm
                            library[sample][gene]['dlabels'][int(s) - gs] += float(counts.at[sample, tr] / l)
                        # normalise to tpkm
                        library[sample][gene]['norm_factor'] += float(counts.at[sample, tr] / l)

    #print("All samples calc in {} seconds".format(time.time() - start_time))
    return library


def labels_squeeze(labels, c):
    ind = np.nonzero(labels)[0]
    if c == 0:
        c = 1
    labels = [round(labels[i] / c, 4) for i in ind]
    labels_ = {}
    for i, l in zip(ind, labels):
        labels_[str(i)] = l
    return labels_


def exons(gene):
    exons = ''
    gs = int(gene['global_start'])
    for key in gene.keys():
        if 'ENST' in key:
            starts = gene[key]['starts']
            ends = gene[key]['ends']
            for x in zip(starts, ends):
                exons += str(x[0] - gs) + ' ' + str(x[1] - gs) + ','
            if exons:
                exons = exons[:-1] + ';'
    return exons


def train_test_sep(gene_dict, test_ratio):

    testN = int(test_ratio * len(gene_dict.keys()))

    genes = []
    lengths = []
    for key in gene_dict.keys():
        genes.extend([key])
        l = gene_dict[key]['global_end'] - gene_dict[key]['global_start']
        lengths.extend([l])
    genes = np.array(genes)
    lengths = np.array(lengths)

    ind = np.argsort(lengths)
    genes_test = genes[ind][-testN:]

    gene_dict_test = {}
    for gene in genes_test:
        val = gene_dict.pop(gene)
        gene_dict_test.update({gene: val})

    return (gene_dict, gene_dict_test)


def save_labels_jsonl(counts, gene_dict, hg38, out_dir, context=1000):
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # normalise to cpm first ALRD NORM IN BATCH NORM R
    #F = counts.sum(axis=1) / 10 ** 6
    #counts = counts.divide(F, axis='index')

    start_time = time.time()
    library = construct_library(gene_dict, counts)

    for sample in library.keys():
        with open(os.path.join(out_dir, sample + '.jsonl'), 'a') as f1:
            # genes
            for gene in library[sample].keys():
                jsonl_dict = {}
                jsonl_dict['gene'] = gene

                c = library[sample][gene]['norm_factor']
                jsonl_dict['alabels'] = labels_squeeze(library[sample][gene]['alabels'], c)
                jsonl_dict['dlabels'] = labels_squeeze(library[sample][gene]['dlabels'], c)

                jsonl_dict['exons'] = exons(gene_dict[gene])

                json.dump(jsonl_dict, f1)
                f1.write('\n')

    if not os.path.exists(os.path.join(out_dir, 'gene_dict.jsonl')):
        with open(os.path.join(out_dir, 'gene_dict.jsonl'), 'a') as f2:
            for gene in library[sample].keys():
                gs = int(gene_dict[gene]['global_start'])
                ge = int(gene_dict[gene]['global_end'])

                seq = str(hg38[gene_dict[gene]['chr']][gs - context: ge + context + 1])

                gene_jsonl_dict = {}
                gene_jsonl_dict[gene] = seq
                json.dump(gene_jsonl_dict, f2)
                f2.write('\n')

    print("All samples saved in {} seconds".format(time.time() - start_time))
    return