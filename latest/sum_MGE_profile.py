#========================================================================
# sum_MGE_profile.py / created by Ryo Honda, 2023-07-17
#========================================================================
# This python script creates a profile comparison table of multiple samples by:
#	$ python3 sum_MGE_profile.py dir_in dir_out
#
# [IMPORTANT] The input file must be ".tsv" files
#  dir_in = directory of sample data files. (all .tsv files in the directory will be merged.)
#  dir_out = directory to output the data file

#  If the dir_in contains "_sample_names.tsv", which has sequence names in the first column and  sample names in the second column, the columns of the output table are labelled with the corresponding sample names.
#------------------------------------------------------------------------------
import glob
import os
import sys
import pandas as pd

#======= Parameters for merging ========================
# A param name should be identical to a column name in the input tables.
# 'param' =value data to be extracted from the sample data. 
param='prop_RPK' 
# 'key' = list of the index parameters to merge
key=['sseqid','MGEDB ID','gene symbol', 'Function', 'slen']
# suffix of the sample data file. 
suffix='.MGE_profile.tsv'
# (Sample names in the output table are automatically named after the input file names by removing this suffix.)

# name of the file containing list of sample names corresponding to sequence names
# (list of allowable file names for compatibility)
file_sample_name=['_sample_names.tsv','sample_names.tsv','_sample_name.tsv','sample_name.tsv']
#======= Parameters for summation ========================
# 'cats' = list of ARG categories to be summed up.
cats=['gene symbol', 'Function']
#=========================================================

####### Creat a merged ARG profiles ###################################
# import arguments 引数処理
args=sys.argv
dir_in=args[1]
dir_out=args[2]

# get the list of sequence files in the input directory
files_in=sorted(glob.glob(os.path.join(dir_in,"*"+suffix)))
# create the dictionary of sample names corresponding to sequence names
dic_sample=pd.DataFrame() # create an empty dataframe
for f in file_sample_name:
    if os.path.isfile(os.path.join(dir_in,f)):
        dic_sample=pd.read_table(os.path.join(dir_in,f), header=None, index_col=0).squeeze(axis=1)
        dic_sample.to_dict()
        break
    else:
        pass

# merge data from input files
df_joined=pd.DataFrame()
key.append(param) # create the list of columns to output
for f in files_in: 
    df_sample=pd.read_table(f,header=0)
    sample_name=os.path.basename(f).replace(suffix,"")
    df_sample=df_sample.reindex(columns=key).rename(columns={param: sample_name})
    if df_joined.empty:
        df_joined=df_sample # for the first data file
    else:
        df_joined=pd.merge(df_joined, df_sample, on=key[:-1], how='outer')
    print("["+args[0]+"] "+f+" merged.")
    
# fill out NaN with zero
df_joined=df_joined.fillna(0)
# rename the columns as sample name
if not dic_sample.empty:
    df_joined=df_joined.rename(columns=dic_sample)

# create the output file
if not os.path.isdir(dir_out):
    os.makedirs(dir_out) # create the output directory if not existed.
file_out='_merged'+suffix[0:-4]+'.'+param+suffix[-4:]
df_joined.to_csv(os.path.join(dir_out,file_out),sep='\t',index=False)
print("["+args[0]+"] "+"Merged profile was created.")
######### Summation by args categories ###############################
# list of the sample columns to aggregate
col_sum=df_joined.columns[len(key)-len(df_joined.columns)-1:]
# sum up by category
for cat in cats:
    # sum up by the key 集計
    df_sum=df_joined.groupby(cat)[col_sum].sum() 

    # output the summary file.
    file_out="MGE."+cat.replace(' ','_')+"."+param+".csv"
    df_sum.to_csv(os.path.join(dir_out,file_out))
    print("["+args[0]+"] "+file_out+" was created.")

print("["+args[0]+"] "+"completed.")
