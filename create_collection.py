#IMPORTS

import subprocess
import os
import sys
import urllib
import argparse
import pandas as pd
import numpy as np
import seaborn as sns
from collections import Counter
from Bio import Entrez
from scipy.spatial.distance import squareform
from scipy.spatial.distance import pdist
import matplotlib.pyplot as plt


#FUNCTIONS
def commandLineParser():
    #make an argument parser that loads metadata, ncbi, and miscellanous files
    parser = argparse.ArgumentParser(description="Create a collection of genomes from NCBI")
    local = os.getcwd()
    parser.add_argument('-m', '--metadata', help='metadata file from GTDB', type=str, required=True)
    parser.add_argument('-n', '--ncbi', help='ncbi file from NCBI', type=str, required=True)
    parser.add_argument('-l', '--leaf', help='leaf file from NCBI', type=str, required=True)
    parser.add_argument('-p', '--phyllo', help='phyllosphere file from NCBI', type=str, required=True)
    parser.add_argument('-d', '--distance', help='distance file from NCBI', type=str, required=True)
    parser.add_argument('-pnas', '--pnas', help='collection from pnas paper', type=str, required=True)
    parser.add_argument('-mash', '--mash', help='mash distance file', type=str, required=True')
    parser.add_argument('-o', '--output', help='output file name', type=str, required=True)
    return vars(parser.parse_args())

def loadVariablesFromArgParse(commandLineArgs):
    #metadata for this work is a bit messy because of the different sources and origins of the data
    metadata = pd.read_csv(commandLineArgs['metadata'], sep='\t', index_col=0, low_memory=False)
    ncbi = pd.read_csv(commandLineArgs['ncbi'], sep=',', index_col=5)
    pnas = pd.read_csv(commandLineArgs['pnas'], sep='\t', header=6, index_col=0)
    leaf = open(commandLineArgs['leaf']).read().split('\n')
    phyllo = open(commandLineArgs['phyllo']).read().split('\n')
    distance = open(commandLineArgs['distance']).read().split('\n')[1:]
    mash = pd.read_csv(commandLineArgs['mash'], sep='\t', header=None)
    #return all variables
    return metadata, ncbi, pnas, leaf, phyllo, distance, mash

def get_assembly_summary(id):
    """Get esummary for an entrez id"""
    esummary_handle = Entrez.esummary(db="assembly", id=id, report="full")
    esummary_record = Entrez.read(esummary_handle)
    return esummary_record

def get_assemblies_refseq(term, download=True, path='assemblies'):
    """Download genbank assemblies for a given search term.
    Args:
        term: search term, usually organism name
        download: whether to download the results
        path: folder to save to
    """

    #provide your own mail here
    handle = Entrez.esearch(db="assembly", term=term, retmax='200')
    record = Entrez.read(handle)
    ids = record['IdList']
    print (f'found {len(ids)} ids')
    links = []
    for id in ids:
        #get summary
        summary = get_assembly_summary(id)
        #get ftp link
        url = summary['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_RefSeq']
        if url == '':
            continue
        label = os.path.basename(url)
        #get the fasta link - change this to get other formats
        link = os.path.join(url,label+'_genomic.fna.gz')
        links.append(link)
        if download == True:
            #download link
            if not os.path.isfile('{}{}.fna.gz'.format(path, label)):
                urllib.request.urlretrieve(link, '{}{}.fna.gz'.format(path, label))
    return links

def get_assemblies_genbank(term, download=True, path='assemblies'):
    """Download genbank assemblies for a given search term.
    Args:
        term: search term, usually organism name
        download: whether to download the results
        path: folder to save to
    """

    #provide your own mail here
    handle = Entrez.esearch(db="assembly", term=term, retmax='200')
    record = Entrez.read(handle)
    ids = record['IdList']
    print (f'found {len(ids)} ids')
    links = []
    for id in ids:
        #get summary
        summary = get_assembly_summary(id)
        #get ftp link
        url = summary['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_GenBank']
        if url == '':
            continue
        label = os.path.basename(url)
        #get the fasta link - change this to get other formats
        link = os.path.join(url,label+'_genomic.fna.gz')
        links.append(link)
        if download == True:
            #download link
            if not os.path.isfile('{}{}.fna.gz'.format(path, label)):
                urllib.request.urlretrieve(link, '{}{}.fna.gz'.format(path, label))
    return links

def find_samID(samID):
    #search for the samID in assembly and biosample
    search_biosample = Entrez.esearch(db='assembly', term=samID)
    assembly_summary = dict(Entrez.read(search_biosample))
    if assembly_summary['Count'] != '0':
        print ('found {} matches in assembly for {}'.format(assembly_summary['Count'], samID))
        print (assembly_summary)
        return assembly_summary
    
    
    search_biosample = Entrez.esearch(db='biosample', term=samID)
    biosample_summary = dict(Entrez.read(search_biosample))
    if biosample_summary['Count'] != '0':
        print ('found {} matches in biosample for {}'.format(assembly_summary['Count'], samID))
        return biosample_summary

def downloadAssemblyFromSource(IDlist, source):
    #download from genbank or refseq (whichever is available)
    #Use the samID to download the assembly
    #create a command that downloads the assembly

    download_cmd = ''
    for ID in IDlist:
        #genbank
        cmd = 'efetch -db BioSample -id {} | elink -db BioSample -target assembly | efetch -format docsum | xtract -pattern DocumentSummary -if BioSampleAccn -equals "{}" -element FtpPath_GenBank | awk -F"/" '.format(samID, samID) + "'" + '{print $0"/"$NF"_genomic.fna.gz"}' + "' | xargs -n1 wget -nc -P /DATA_RAID2/vtracann/shared/db/ncbi/{}".format(source)
        download_cmd += cmd + '\n'
        #refseq
        cmd = 'efetch -db BioSample -id {} | elink -db BioSample -target assembly | efetch -format docsum | xtract -pattern DocumentSummary -if BioSampleAccn -equals "{}" -element FtpPath_RefSeq | awk -F"/" '.format(samID, samID) + "'" + '{print $0"/"$NF"_genomic.fna.gz"}' + "' | xargs -n1 wget -nc -P /DATA_RAID2/vtracann/shared/db/isolates/ncbi/{}".format(source)
        download_cmd += cmd + '\n'

    if not os.path.isfile('/DATA_RAID2/vtracann/shared/db/isolates/ncbi/{}_ncbi_wget.sh'.format(source)):
        outfile = open('/DATA_RAID2/vtracann/shared/db/isolates/ncbi/{}_ncbi_wget.sh'.format(source), 'w')
        outfile.write(tot_cmd)
        outfile.close()

def make_dataFrame_distance(distance):
    if not os.path.isfile('/DATA_RAID2/vtracann/shared/db/isolates/global/unzipped/pd_complete_collection_distance.tsv'):
        Index = []
        for I in distance:
            i = I.split('\t')[0]
            Index.append(i)

        df = pd.DataFrame(0, index=Index, columns=Index)
        for I in range(len(distance)):
            vals = distance[I].split('\t')[1:]
            for v in range(len(vals)):
                df.iloc[I, v] = vals[v]
                df.iloc[v, I] = vals[v]

        df.to_csv('/DATA_RAID2/vtracann/shared/db/isolates/global/unzipped/pd_complete_collection_distance.tsv',sep='\t')
    else:
        df = pd.read_csv('/DATA_RAID2/vtracann/shared/db/isolates/global/unzipped/pd_complete_collection_distance.tsv', sep='\t', index_col=0)
    return df

def make_dataFrame_mash(mash):
    if not os.path.isfile('/DATA_RAID2/vtracann/shared/db/isolates/global/unzipped/pd_mash_distance.tsv'):
        value_df = mash.loc[:, [0,1,4]]
        value_df.columns = ['id1', 'id2', 'kmers']
        multiI = value_df.set_index(['id1','id2'])

        nonzero_df = value_df[value_df.loc[:, 'kmers']!= '0/1000']

        I = list(set(nonzero_df.loc[:, 'id1'].values))
        DF = pd.DataFrame(0, index=I, columns=I)

        for i in multiI.index:
            i1 = i[0]
            i2 = i[1]
            DF.loc[i1, i2] = multiI.loc[(i1, i2)].values[0]

        print (DF)
        DF.to_csv('/DATA_RAID2/vtracann/shared/db/isolates/global/unzipped/pd_mash_distance.tsv',sep='\t')

    else:
        DF = pd.read_csv('/DATA_RAID2/vtracann/shared/db/isolates/global/unzipped/pd_mash_distance.tsv', sep='\t', index_col=0)

def demultiplex_mash(mash_distance_df):
    mash_distance_df = mash_distance_df.loc[mash_distance_df.index.dropna(), mash_distance_df.index.dropna()]

    grouped_entries = {}
    excluded = []
    included = []

    for i in mash_distance_df.index:
        if (sum(mash_distance_df.loc[i] == 0)*1) > 1:
            duplicates = mash_distance_df[mash_distance_df.loc[i] == 0]
            excluded.extend(list(duplicates.index))
            already_found = False
            for d in duplicates.index:
                if d in included:
                    already_found = True
            if not already_found:
                included.append(i)
        else:
            included.append(i)
            
    demux_df = mash_distance_df.loc[included,included]
    
    plt.figure(figsize=(100,100))
    g = sns.clustermap(demux_df, cmap='Greens_r', xticklabels=False, yticklabels=False)
    order = g.dendrogram_row.reordered_ind
    order_names = [demux_df.index[x] for x in order]
    plt.show()

    sns.clustermap(demux_df, cmap='Greens_r', xticklabels=False, yticklabels=False)
    plt.show()

    demux_df.to_csv('/DATA_RAID2/vtracann/shared/db/isolates/global/unzipped/pd_demux_collection_distance.tsv',sep='\t')

    return demux_df

if __name__ == "__main__":
    #get the metadata
    commandLineArgs = commandLineParser()
    metadata, ncbi, pnas, leaf, phyllo, distance, mash = loadVariablesFromArgParse(commandLineArgs)

    #get the number of isolates per source
    isolates_by_isolation_source = Counter(metadata.loc[:,'ncbi_isolation_source'].values)
    isolates_by_isolation_source = {k: v for k, v in sorted(isolates_by_isolation_source.items(), key=lambda item: item[1], reverse=False)}

    #Choose which sources from the metadata to include in the collection
    selected_sources = ['leaf', 'root', 'tomato', 'Plant', 'plant', 'Rice seed', 'legume-root nodule', 'Plant material', 'citrus infected tissue', 'root nodules', 'tomato rhizosphere', 'rhizosphere', 'Leaf Litter', 'sugarcane', 'rhizosphere soil', 'Walnut blight lesion', 'internal stem tissue', 'roots', 'root nodule', 'grass', 'potato rhizosphere', 'rice leaves', 'flower', 'stem', 'plant root', 'Root nodule', 'Rice', 'sugar beet rhizosphere', 'Citrus sp.', 'plant (cereals)']
    selected_sources = ['leaf', 'root', 'tomato', 'Plant', 'plant', 'Rice seed', 'Plant material', 'citrus infected tissue','sugarcane','Walnut blight lesion', 'internal stem tissue', 'roots', 'rice leaves', 'flower', 'stem', 'plant root', 'Rice', 'Citrus sp.', 'plant (cereals)']

    #turn the gtdb genome id into ncbi genome id
    index_GC = [x[3:] for x in metadata.index]
    metadata.index = index_GC

    pairwise_distance_df = make_dataFrame_distance(distance)
    mash_distance_df = make_dataFrame_mash(mash)

    demultiplex_mash(mash_distance_df)


