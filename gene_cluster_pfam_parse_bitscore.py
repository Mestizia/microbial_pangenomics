import pandas as pd
import numpy as np
import sys
import seaborn as sns
import matplotlib.pyplot as plt
from itertools import islice
import os

def parse_table(table_file):
    genome_pfam_dict = {}
    with open(table_file, 'r') as f:
        for line in islice(f, 3, None):
            try:
                no_space_l = line.split()
                bitscore = float(no_space_l[5])
                if bitscore >= 10:
                    #record only significant hits
                    pfam = no_space_l[0]
                    genomeID = no_space_l[2].split('|')[2].split(":")[1]
                    if genomeID in genome_pfam_dict.keys():
                        genome_pfam_dict[genomeID].append(pfam)
                    else:
                        genome_pfam_dict[genomeID] = [pfam]
            except (ValueError, IndexError) as e:
                #line to discard, no pfam hit
                continue
    return genome_pfam_dict


if __name__ == "__main__":
    table_file = '/DATA_RAID2/vtracann/shared/db/isolates/anvio/filtered_panpangenome/gene_clusters_pfam_hits35.tsv'
    annot_dict = parse_table(table_file)
    table_file2 = '/DATA_RAID2/vtracann/shared/db/isolates/anvio/filtered_panpangenome/gene_clusters_pfam_hits35_2.tsv'
    annot_dict2 = parse_table(table_file2)

    print (len(annot_dict.keys()))
    print (len(annot_dict2.keys()))


    pfams = [set(annot_dict[k]) for k in annot_dict.keys()]
    pfams2 = [set(annot_dict2[k]) for k in annot_dict2.keys()]

    pfams = set.union(*pfams)
    pfams2 = set.union(*pfams2)
    print (len(pfams), len(pfams2), len(pfams.union(pfams2)))
    pfam_list = list(pfams.union(pfams2))

    #turn the dictionary into a dataframe of binary presence/absence
    pfam_representation_genomes = pd.DataFrame(0, columns=pfam_list, index=annot_dict.keys())
    print (pfam_representation_genomes.shape)


    selected_pfams = '5-FTHF_cyc-lig,AAA_25,AAA_assoc_C,AcetylCoA_hyd_C,ASH,Big_3_5,CitMHS,CpsB_CapC,CpXC,CsbD,Cys_rich_CPXG,Cytidylate_kin2,CytoC_RC,DHquinase_I,Exo_endo_phos,FAD_binding_7,Glucodextran_N,GWxTD_dom,HipA_2,HTH_37,HupF_HypC,HycI,HypD,Ig_GlcNase,KdpA,KdpC,KdpD,Lactonase,MCRA,Met_gamma_lyase,Methyltransf_4,MgtE,MNHE,MrpF_PhaF,OprB,Paired_CXXCH_1,PDZ_2,PGI,PhaG_MnhG_YufB,Phenol_MetA_deg,Phosphoesterase,Polbeta,PQQ,Pro-kuma_activ,SelO,SNase,SOUL,TctA,TelA,TFR_dimer,TPP_enzyme_C,TrbI,UvdE,WXXGXW,YHS,zf-CDGSH'.split(',')

    for genome in annot_dict.keys():
        for pfam in selected_pfams:
            if pfam in annot_dict[genome]:
                pfam_representation_genomes.loc[genome, pfam] = 1
    

    print (pfam_representation_genomes.loc[:, selected_pfams])
    print (pfam_representation_genomes.loc[:, selected_pfams].sum(axis=1))
    print (pfam_representation_genomes.loc[:, selected_pfams].sum(axis=0))

    pfam_representation_genomes.loc[:, selected_pfams].to_csv('collection_ph_genes.csv', sep='\t')