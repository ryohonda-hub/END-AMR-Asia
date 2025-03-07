#==============================================================================
# make_MGE_prof.py ver.2.1 / created by Ryo Honda, Last updated: 2025-03-07
#==============================================================================
# This python script creates ARG profile data by merging gene information from the CARD catalog with read count data using ARO as index by:
#	$ python3 make_MGE_prof2.py catalog_file blast_results dir_out
#
#  catalog_file = reference MGEDB catalog (MGE_tax_table_trim.txt)
#  blast_results = a blast output file
#  dir_out = output directory, or use "." to specify the current directory
#
# The blast output format should be in -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen".
# The output file with suffix of ".MGE_profile.tsv" is created in the output directory. 
# The output file contains gene information and proportion of RPK.
#------------------------------------------------------------------------------
import os
import re
import sys
import pandas as pd

def main(file_cat, file_blast, dir_out):
    # cutoff conditions of the blast result
    th_mlen=50     # (bp) length of the matched sequence below this value will be excluded.
    th_pident=90   # (%) the pident below this value will be excluded.
    # import the database catalog and blast results
    df_catalog=pd.read_table(file_cat,header=0)
    df_blast=pd.read_table(file_blast,header=0)
    # remove empty columns in MGEDB catalog and define column titles.
    df_catalog=df_catalog.dropna(how='all', axis=1) 
    df_catalog.columns=['sseqid','Function','gene name']
    
    # exclude the blast results with pident or length below the cutoff thresholds
    df_blast.columns=['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'qlen', 'slen']
    df_blast=df_blast.query('pident >= @th_pident & length >= @th_mlen')
    # count blast hits
    df_blast['reads']=1
    df_blast=df_blast.reindex(columns=['sseqid','slen','reads'])
    df_count=df_blast.groupby(['sseqid','slen'], as_index=False)['reads'].sum()
    # calculate RPK and proportion of RPK
    df_count['RPK']=df_count['reads'] / df_count['slen'] * 1000
    df_count['prop_RPK']=df_count['RPK'] / df_count['RPK'].sum()
    df_count=df_count.sort_values('RPK',ascending=False)
    
    # separate MGEDB ID and gene symbol from the sseqid in the read count file.
    df_count.insert(0,'MGEDB ID','')
    df_count.insert(1,'gene symbol','')
    for i, gene in df_count.iterrows():
        sseqid=gene.sseqid
        df_count.at[i, 'MGEDB ID']= re.split('_',sseqid)[0] # get MGEDB ID
        df_count.at[i, 'gene symbol'] = re.split('_',sseqid)[1] # get gene symbol 
    
    # lookup and merge gene information from CARD catalog
    df_prof_raw=df_count.merge(df_catalog, on='sseqid', how='left')
    # re-aggregate the count data.
    df_prof=df_prof_raw.reindex(columns=['Function', 'gene symbol','reads','RPK','prop_RPK'])
    df_prof=df_prof.fillna({'Function':'NA'})
    df_prof=df_prof.groupby(['gene symbol','Function'], as_index=False).sum()
    df_prof=df_prof.sort_values('RPK',ascending=False)
    
    # create the output file
    if not os.path.exists(dir_out):
        os.makedirs(dir_out)
    file_out_base=os.path.splitext(os.path.basename(file_blast))[0].rstrip('.blast')
    
    df_out=df_prof_raw.reindex(columns=['sseqid','MGEDB ID', 'gene symbol','Function', 'gene name','slen','reads','RPK','prop_RPK'])
    file_out=os.path.join(dir_out,file_out_base+".MGE_profile.tsv")
    df_out.to_csv(file_out,sep='\t',index=False)
    
    df_out=df_prof
    file_out=os.path.join(dir_out,file_out_base+".MGE_profile.summary.tsv")
    df_out.to_csv(file_out,sep='\t',index=False)
    
if __name__=="__main__":
    args=sys.argv
    main(args[1],args[2],args[3])
