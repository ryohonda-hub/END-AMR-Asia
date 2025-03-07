#==============================================================================
# make_ARG_profile.py / created by Ryo Honda (END-AMR-Asia), 2023-06-05
#==============================================================================
# This python script creates ARG profile data by merging gene information from the CARD catalog with read count data using ARO as index by:
#	$ python make_ARG_profile.py arg1 arg2 arg3
#
#  arg1 = reference CARD catalog (aro_index.tsv)
#  arg2 = a table with read count of ARGs and sseqid of genbank
#  arg3 = output directory, or use "." to specify the current directory
#
# The output file with suffix of ".ARG_profile.tsv" is created in the output directory. 
# The output file contains gene information and proportion of RPK.
#------------------------------------------------------------------------------
import os
import re
import sys
import pandas as pd

args=sys.argv
df_catalog=pd.read_table(args[1],header=0)
df_count=pd.read_table(args[2],header=0)

# separate DNA accession and gene symbol from the sseqid in the read count file.
df_count.insert(0,'ARO Accession','')
df_count.insert(1,'gene symbol','')
for i, gene in df_count.iterrows():
    sseqid=gene.sseqid
    df_count.at[i, 'ARO Accession']= re.split('\|', sseqid)[4] # get ARO accession ID
    str = re.split('\|', sseqid)[5] # get gene symbol 
    # remove allele prefix, which is numbers and alphabets in lower cases followed by a hyphen '-' or a period '.' .
    df_count.at[i, 'gene symbol'] = re.sub(r'[-|.][0-9a-z]+$','',str)

# calculate the RPK proportion of each ARG to all ARGs. 
sumRPK=df_count['RPK'].sum()
df_count['prop_RPK']=df_count['RPK']/sumRPK
# sort the data by RPK
df_count=df_count.sort_values('RPK',ascending=False)

# lookup and merge gene information from CARD catalog
df_prof=df_count.merge(df_catalog, on='ARO Accession', how='left')
# slim the description of drug classes and mechanisms
df_prof['Drug Class']=df_prof['Drug Class'].str.replace(' antibiotic','')
df_prof['Resistance Mechanism']=df_prof['Resistance Mechanism'].str.replace('antibiotic ','')
# count the resistant drug classes as MAR
df_prof['MAR']=df_prof['Drug Class'].str.count(';')+1

# create the output file
df_out=df_prof.reindex(columns=['sseqid','ARO Accession', 'gene symbol','CARD Short Name', 'AMR Gene Family','Drug Class','MAR','Resistance Mechanism','slen','reads','RPK','prop_RPK'])
file_out=os.path.splitext(os.path.basename(args[2]))[0]
file_out=os.path.join(args[3],file_out.rstrip('.count')+".ARG_profile.tsv")
df_out.to_csv(file_out,sep='\t',index=False)


