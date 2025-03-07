#========================================================================
# agg_ARG_profile.py / created by Ryo Honda (END-AMR-Asia), 2023-06-05
#========================================================================
# This python script creates a profile comparison table of multiple samples by:
#	$ python agg_ARG_profile.py dir_in dir_out
#
# [IMPORTANT] The input file must be ".tsv" files
#  dir_in = directory of sample data files. (all .tsv files in the directory will be merged.)
#  dir_out = directory to output the data file
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
key=['ARO Accession','gene symbol','CARD Short Name', 'Drug Class', 'MAR', 'Resistance Mechanism', 'slen']
# suffix of the sample data file. 
suffix='.ARG_profile.tsv'
# (Sample names in the output table are automatically named after the input file names by removing this suffix.)

#======= Parameters for summation ========================
# 'cats' = list of ARG categories to be summed up.
cats=['MAR', 'Resistance Mechanism','gene symbol','CARD Short Name']
#=========================================================

####### Creat a merged ARG profiles ###################################
# import arguments 引数処理
args=sys.argv
dir_in=args[1]
dir_out=args[2]

# get the file list in the input directory
files_in=sorted(glob.glob(os.path.join(dir_in,"*.tsv")))

# merge data from input files
df_joined=pd.DataFrame()
key.append(param) # create the list of columns to output
for f in files_in: 
    df_sample=pd.read_table(f,header=0)
    sample_name=os.path.basename(f).rstrip(suffix)
    df_sample=df_sample.reindex(columns=key).rename(columns={param: sample_name})
    if df_joined.empty:
        df_joined=df_sample # for the first data file
    else:
        df_joined=pd.merge(df_joined, df_sample, on=key[:-1], how='outer')
    print(f+" merged.")
    
# fill out NaN with zero
df_joined=df_joined.fillna(0)

# create the output file
file_out='merged.'+param+suffix
df_joined.to_csv(os.path.join(dir_out,file_out),sep='\t',index=False)
print("Merged profile was created.")

######### Summation by args categories ###############################

for cat in cats:
    # sum up by the key 集計
    df_sum=df_joined.groupby(cat).sum() 
    # delete unnecessary columns 不要な列を削除
    #df_sum=df_sum.drop(df_sum.columns[[0,1,2,3,4,5]],axis=1) 
    df_sum=df_sum.drop(df_sum.iloc[:,0:len(key)-2],axis=1) 

    # output the summary file.
    file_out=cat.replace(' ','_')+"."+param+".sum.csv"
    df_sum.to_csv(os.path.join(dir_out,file_out))
    print(file_out+" was created.")

print(args[0]+" completed.")
