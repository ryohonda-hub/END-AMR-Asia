#==============================================================================
# make_ARG_prof.py ver.2 / created by Ryo Honda, 2023-07-10
#==============================================================================
# This python script creates ARG profile data by merging gene information from the CARD catalog with read count data using ARO as index by:
#	$ python3 make_ARG_prof2.py arg1 arg2 arg3
#
#  arg1 = reference CARD catalog (aro_index.tsv)
#  arg2 = a blast output file
#  arg3 = output directory, or use "." to specify the current directory
#
# The blast output format should be in -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen".
# The output file with suffix of ".ARG_profile.tsv" is created in the output directory. 
# The output file contains gene information and proportion of RPK.
#------------------------------------------------------------------------------
import os
import re
import sys
import pandas as pd

args=sys.argv
df_catalog=pd.read_table(args[1],header=0)
df_blast=pd.read_table(args[2],header=None)
dir_out=args[3]

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

# separate DNA accession and gene symbol from the sseqid in the read count file.
df_count.insert(0,'ARO Accession','')
df_count.insert(1,'gene symbol','')
for i, gene in df_count.iterrows():
    sseqid=gene.sseqid
    df_count.at[i, 'ARO Accession']= re.split('\|', sseqid)[4] # get ARO accession ID
    str = re.split('\|', sseqid)[5] # get gene symbol 
    # remove allele prefix, which is numbers and alphabets in lower cases followed by a hyphen '-' or a period '.' .
    df_count.at[i, 'gene symbol'] = re.sub(r'[-|.][0-9a-z]+$','',str)

# lookup and merge gene information from CARD catalog
df_prof=df_count.merge(df_catalog, on='ARO Accession', how='left')
# slim the description of drug classes and mechanisms
df_prof['Drug Class']=df_prof['Drug Class'].str.replace(' antibiotic','')
df_prof['Resistance Mechanism']=df_prof['Resistance Mechanism'].str.replace('antibiotic ','')
# count the resistant drug classes as MAR
df_prof['MAR']=df_prof['Drug Class'].str.count(';')+1

# create the output file
if not os.path.exists(dir_out):
    os.makedirs(dir_out)
df_out=df_prof.reindex(columns=['sseqid','ARO Accession', 'gene symbol','CARD Short Name', 'AMR Gene Family','Drug Class','MAR','Resistance Mechanism','slen','reads','RPK','prop_RPK'])
file_out=os.path.splitext(os.path.basename(args[2]))[0]
file_out=os.path.join(dir_out,file_out.rstrip('.blast')+".ARG_profile.tsv")
df_out.to_csv(file_out,sep='\t',index=False)


