import pandas as pd

def parsePfam2GO(pfam2goDB):
    pfam2go_dict = {}
    for pfamLine in pfam2goDB.split('\n')[6:-1]:
        PFAM = pfamLine.split(' ')[0].split(':')[1]
        PFAM_description = pfamLine.split(' ')[1].split(' >')[0]
        GO_ID = pfamLine.split('; ')[-1]
        GO_description = pfamLine.split('> ')[1].split(':')[1].split(';')[0][:-1]
        if PFAM in pfam2go_dict.keys():
            pfam2go_dict[PFAM_description].append([PFAM, GO_ID, GO_description])
        else:
            pfam2go_dict[PFAM_description] = [[PFAM, GO_ID, GO_description]]
    return pfam2go_dict



if __name__ == "__main__":
    output_dir = '/DATA_RAID2/vtracann/shared/db/isolates/'
    pfam2go_raw = open('/DATA_RAID2/vtracann/shared/db/pfam2go'.format(output_dir)).read()

    enriched_pfam = pd.read_csv('{}pfam/consensus_enriched_df.csv'.format(output_dir), sep=',', index_col=0)
    depleted_pfam = pd.read_csv('{}pfam/consensus_depleted_df.csv'.format(output_dir), sep=',', index_col=0)

    pfam2go_dict = parsePfam2GO(pfam2go_raw)

    for pfam in enriched_pfam.index:
        if pfam in pfam2go_dict.keys():
            for go in pfam2go_dict[pfam]:
                print (go)
                enriched_pfam.loc[pfam, 'GO_ID'] = go[1]
                enriched_pfam.loc[pfam, 'GO_description'] = go[2]
        else:
            enriched_pfam.loc[pfam, 'GO_ID'] = 'Missing'
            enriched_pfam.loc[pfam, 'GO_description'] = 'Missing'

    
    for pfam in depleted_pfam.index:
        if pfam in pfam2go_dict.keys():
            for go in pfam2go_dict[pfam]:
                depleted_pfam.loc[pfam, 'GO_ID'] = go[1]
                depleted_pfam.loc[pfam, 'GO_description'] = go[2]
        else:
            depleted_pfam.loc[pfam, 'GO_ID'] = 'Missing'
            depleted_pfam.loc[pfam, 'GO_description'] = 'Missing'
    
    enriched_pfam['ratio'] = enriched_pfam['sphingo_incidence'] - enriched_pfam['non_sphingo_incidence']
    depleted_pfam['ratio'] = depleted_pfam['sphingo_incidence'] - depleted_pfam['non_sphingo_incidence']
    enriched_pfam.sort_values(by='ratio', ascending=False).to_csv('{}pfam/GO_consensus_enriched_df.csv'.format(output_dir), sep=',')
    depleted_pfam.sort_values(by='ratio', ascending=False).to_csv('{}pfam/GO_consensus_depleted_df.csv'.format(output_dir), sep=',')

    print (enriched_pfam.sort_values(by='ratio', ascending=False))


import pandas as pd
import os
import sys
from collections import Counter

eggnogg_annot = pd.read_csv('/DATA_RAID2/vtracann/shared/db/isolates/anvio/filtered_panpangenome/eggnogg/filtered_panpangenome.emapper.annotations', sep='\t', index_col=0, skiprows=[0,1,2,3])
eggnogg_annot = eggnogg_annot.iloc[:-3]
eggnogg_annot = eggnogg_annot.reset_index()

go_terms = eggnogg_annot['GOs']
