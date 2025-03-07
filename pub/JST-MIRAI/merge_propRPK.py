#==============================================================================
# merge_propRPK.py, for ARG profiling / created by Ryo Honda, 2023-05-28
#==============================================================================
# This python script creates a comparison table of a parameter from data tables of multiple samples by:
#	$ python merge_propRPK.py file_out dir_in, file_in1, file_in2, ... 
#
#  dir_in = directory of sample data files (all input files should be in a single directory.)
#  file_out = output file: a comparison table of the samples
#  file_in = input files: data tables of the samples including the target parameter to compare
#
# [IMPORTANT]
# The following params should be set in the below script.
#  'key' = the index parameter to merge
#  'param' = the parameter of value data to be extracted from the sample data
#------------------------------------------------------------------------------
import os
import sys
import pandas as pd

#======= SET index and data parameters to merge ==========
# The param names should be identical to the target column names in the input tables.
key='CARD Short Name'
param='prop_RPK'
# suffix of the sample data file (to be removed from the column title, which is automatically named after the input data file.)
suffix='.ARG_profile.tsv'
#=========================================================

#########################################################
# import arguments 引数処理
args=sys.argv
file_out=args[1]
dir=args[2].rstrip('/')+'/'

# merge data from files
df_joined=pd.DataFrame()
for data in args[3:]: 
    df_sample=pd.read_table(dir+data,header=0)
    sample_name=os.path.basename(data).rstrip(suffix)
    df_sample=df_sample.reindex(columns=[key,param]).rename(columns={param: sample_name})
    if df_joined.empty:
        df_joined=df_sample # for the first data file
    else:
        df_joined=pd.merge(df_joined, df_sample, on=key, how='outer')
# fill out NaN with zero
df_joined=df_joined.fillna(0)

# create the output file
df_joined.to_csv(file_out,sep='\t',index=False)
