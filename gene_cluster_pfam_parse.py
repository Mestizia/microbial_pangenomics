import pandas as pd
import numpy as np
import sys
import seaborn as sns
import matplotlib.pyplot as plt
from itertools import islice
import os

def parse_table(table_file):
    annot_dict = {}
    last_gc = ''
    with open(table_file, 'r') as f:
        for line in islice(f, 3, None):
            try:
                no_space_l = line.split()
                eval = float(no_space_l[4])
                if eval < 0.05:
                    #record only significant hits
                    gene_cluster = line.split('|')[1]
                    pfam = no_space_l[0]
                    if gene_cluster != last_gc:
                        annot_dict[gene_cluster] = [pfam]
                        last_gc = gene_cluster
                    annot_dict[gene_cluster].append(pfam)
                else:
                    continue
            except (ValueError, IndexError) as e:
                #line to discard, no pfam hit
                continue
    return annot_dict


def condense_annotation(annot_dict):
    gc_annotation = {}
    for gc in annot_dict.keys():
        pfam_list = annot_dict[gc]
        hit_count = pd.Series(pfam_list).value_counts()
        annot = get_consensus_annot(hit_count)
        #best_hit = hit_count.loc[hit_count.index[0]]
        gc_annotation[gc] = annot
        if annot == ['no consensus']:
            print (hit_count)
            print (sum(hit_count))
    return gc_annotation


def get_consensus_annot(hit_count):
    tot_annot = hit_count.sum()
    if hit_count.iloc[0] / tot_annot > 0.7:
        hits = [hit_count.index[0]]
        for x in range(1, len(hit_count)):
            if hit_count.iloc[x] / tot_annot > 0.5:
                hits.append(hit_count.index[x])
        return hits
    else:
        if hit_count.iloc[0] / tot_annot > 0.2:
            hits = [hit_count.index[0]]
            #if at least 30% of hits are the same, check if there are other similar hits
            for x in range(1, len(hit_count)):
                if hit_count.iloc[x] / tot_annot > 0.2:
                    hits.append(hit_count.index[x])
                    #if another hit is at least 30% of total hits, return also that. stop if not
                else:
                    return hits
    return ['no consensus']

if __name__ == "__main__":
    table_file = '/DATA_RAID2/vtracann/shared/db/isolates/anvio/filtered_panpangenome/gene_clusters_pfam_hits35.tsv'
    annot_dict = parse_table(table_file)
    table_file2 = '/DATA_RAID2/vtracann/shared/db/isolates/anvio/filtered_panpangenome/gene_clusters_pfam_hits35_2.tsv'
    annot_dict2 = parse_table(table_file2)

    gc_annotation = condense_annotation(annot_dict)
    gc_annotation2 = condense_annotation(annot_dict2)

    #print (annot_dict['gene_cluster:GC_00007329'])
    #count the occurrences of every element in the list
    #hit_count = pd.Series(annot_dict['gene_cluster:GC_00007329']).value_counts()
    #print (hit_count)
    #print (len(annot_dict['gene_cluster:GC_00007329']))
    #print (annot_dict['gene_cluster:GC_00027313'])
    #print (len(annot_dict['gene_cluster:GC_00027313']))
    #hit_count = pd.Series(annot_dict['gene_cluster:GC_00027313']).value_counts()
    # print (hit_count)

    #print (annot_dict['gene_cluster:GC_00007329'])
    print (gc_annotation['gene_cluster:GC_00007329'])
    print (pd.Series(annot_dict['gene_cluster:GC_00007329']).value_counts())
    print (len(annot_dict['gene_cluster:GC_00007329']))


    if not os.path.isfile('a'):

        out_txt = ''
        for annot in gc_annotation.keys():
            out_txt += annot + '\t' + '\t'.join(gc_annotation[annot]) + '\n'
        for annot in gc_annotation2.keys():
            out_txt += annot + '\t' + '\t'.join(gc_annotation2[annot]) + '\n'

        outfile = '/DATA_RAID2/vtracann/shared/db/isolates/anvio/filtered_panpangenome/gene_clusters_pfam_hits35_condensed.tsv'
        
        with open(outfile, 'w') as f:
            f.write(out_txt)
    
