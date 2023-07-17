#==============================================================================
# make_MGE_prof.py ver.2/ created by Ryo Honda, 2023-07-09
#==============================================================================
# This python script creates ARG profile data by merging gene information from the CARD catalog with read count data using ARO as index by:
#	$ python3 make_MGE_prof2.py arg1 arg2 arg3
#
#  arg1 = reference MGEDB catalog (MGE_tax_table_trim.txt)
#  arg2 = a blast output file
#  arg3 = output directory, or use "." to specify the current directory
#
# The blast output format should be in -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen".
# The output file with suffix of ".MGE_profile.tsv" is created in the output directory. 
# The output file contains gene information and proportion of RPK.
#------------------------------------------------------------------------------
import os
import re
import sys
import pandas as pd

args=sys.argv
df_catalog=pd.read_table(args[1],header=None)
df_blast=pd.read_table(args[2],header=None)
dir_out=args[3]

# remove empty columns in MGEDB catalog and define column titles.
df_catalog=df_catalog.dropna(how='all', axis=1) 
df_catalog.columns=['sseqid','Function','gene name']

# exclude the blast results with pident <90% or length <25bp
df_blast.columns=['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'qlen', 'slen']
df_blast=df_blast.query('pident >=90 & length >=25')

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
file_out_base=os.path.splitext(os.path.basename(args[2]))[0].rstrip('.blast')

df_out=df_prof_raw.reindex(columns=['sseqid','MGEDB ID', 'gene symbol','Function', 'gene name','slen','reads','RPK','prop_RPK'])
file_out=os.path.join(dir_out,file_out_base+".MGE_profile.tsv")
df_out.to_csv(file_out,sep='\t',index=False)

df_out=df_prof
file_out=os.path.join(dir_out,file_out_base+".MGE_profile.summary.tsv")
df_out.to_csv(file_out,sep='\t',index=False)