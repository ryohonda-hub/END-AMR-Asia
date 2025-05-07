#==============================================================================
# ARGI_prof.py ver.1.0 / created by Ryo Honda, Last updated: 2025-05-07
#==============================================================================
# This python script creates ARG profile data by merging gene information from the CARD catalog (card.json) with read count data created by RGI by:
#	$ python3 ARGIprof.py card_json rgi_results dir_out
#
#  card_json = reference CARD catalog (card.json)
#  rgi_results = an RGI output file
#  dir_out = output directory, or use "." to specify the current directory
#
# The output file with suffix of ".ARGIprof.tsv" is created in the output directory. 
# The output file contains gene information and proportion of RPK.
#------------------------------------------------------------------------------
import os
import sys
import json
import pandas as pd

def main(file_cat, file_rgi, dir_out):
    #=== cutoff thresholds of the blast result ===
    th_mlen=100     # (bp) hits with matched length below this value are excluded.
    th_pident=90    # (%) hits with Best_Identities below this value are excluded.
    #=============================================
    print(f"ARGIprof.py - threshold: length>{th_mlen} bp, pident>{th_pident}%")
    
    # import arg catalog and rgi results
    df_rgi=pd.read_table(file_rgi,header=0)
    # append gene length in the dataframe
    df_rgi['qlen'] = df_rgi['Stop']-df_rgi['Start']+1
    df_rgi['Gene Length'] = df_rgi['Predicted_DNA'].str.len()
    
    # cutoff the rgi results with 'loose" and lower length and identity than the cutoff thresholds
    cond1 = df_rgi[df_rgi['Cut_Off'].isin(['Perfect', 'Strict'])].copy()
    cond2 = df_rgi[
        (df_rgi['Cut_Off'] == 'Loose') &
        (df_rgi['Best_Identities'] >= th_pident) &
        (df_rgi['qlen'] >= th_mlen)
    ].copy()
    df_rgi = pd.concat([cond1, cond2], ignore_index=True)
    
    # count rgi hits, get the average length of the aggregated genes
    df_count = (
        df_rgi.groupby(['Model_ID'],as_index=False)
            .agg(reads=('Contig', 'count'), slen=('Gene Length', 'mean'))
    )
    # calculate RPK and proportion of RPK
    df_count['RPK']=df_count['reads'] / df_count['slen'] * 1000
    df_count['prop_RPK']=df_count['RPK'] / df_count['RPK'].sum()
    df_prof=df_count.sort_values('RPK',ascending=False)
    
    # map ARO accession, gene names and metadata in the dataframe
    unique_rgi = df_rgi.drop_duplicates('Model_ID').set_index('Model_ID')
    df_prof['ARO Accession'] = df_prof['Model_ID'].map(unique_rgi['ARO'])
    df_prof['CARD Short Name'] = df_prof['Model_ID'].map(unique_rgi['Best_Hit_ARO'])
    df_prof['AMR Gene Family'] = df_prof['Model_ID'].map(unique_rgi['AMR Gene Family'])
    df_prof['Contig'] = df_prof['Model_ID'].map(unique_rgi['Contig'])
    df_prof['Drug Class'] = df_prof['Model_ID'].map(unique_rgi['Drug Class'])
    df_prof['Resistance Mechanism'] = df_prof['Model_ID'].map(unique_rgi['Resistance Mechanism'])
    # append gene symbols and slim the description of drug classes and mechanisms
    df_prof['gene symbol']=df_prof['CARD Short Name'].str.split(' ').str[0]
    df_prof['Drug Class']=df_prof['Drug Class'].str.replace(' antibiotic','', regex=False)
    df_prof['Resistance Mechanism']=df_prof['Resistance Mechanism'].str.replace('antibiotic ','', regex=False)
    # count the resistant drug classes as MAR
    df_prof['MAR']=df_prof['Drug Class'].str.count(';')+1
    
    # create the output file
    if not os.path.exists(dir_out):
        os.makedirs(dir_out)
    df_out=df_prof.reindex(columns=['Contig','ARO Accession','gene symbol','CARD Short Name', 'AMR Gene Family','Drug Class','MAR','Resistance Mechanism','slen','reads','RPK','prop_RPK'])
    file_out=os.path.splitext(os.path.basename(file_rgi))[0]
    file_out=os.path.join(dir_out,file_out.removesuffix('._ARGI')+".ARGIprof.tsv")
    df_out.to_csv(file_out,sep='\t',index=False)

if __name__=="__main__":
    args=sys.argv
    main(args[1],args[2],args[3])
