#========================================================================
# curate_16S_taxon.py ver.2 / created by Ryo Honda, Last updated: 2023-07-21
#========================================================================
# This python script creates a profile comparison table of multiple samples from mpa-style taxonomy read-count data by:
#	$ python3 crt16S.py dir_in dir_out
#
# [IMPORTANT] The input file must be ".tsv" files
#  dir_in = directory of sample data files. (all .tsv files in the directory will be merged.)
#  dir_out = directory to output the data file
#
#  If the dir_in contains "_sample_names.tsv", which has sequence names in the first column and  sample names in the second column, the columns of the output table are labelled with the corresponding sample names.
#------------------------------------------------------------------------------
import glob
import re
import os
import sys
import pandas as pd

def main(dir_in, dir_out):
    #======= Parameter settings ========================
    # A param name should be identical to a column name in the input tables.
    # 'param' =value data to be extracted from the sample data. 
    param='abundance' 
    # 'key' = list of the index parameters to merge
    key=['taxonomy']
    # suffix of the sample data file. 
    suffix='.mpa.tsv'
    # The taxonomy with lower abundance than the below threshold will be ignored.
    thres= 10 ** -5
    
    # name of the file containing list of sample names corresponding to sequence names
    # (list of allowable file names for compatibility)
    file_sample_name=['_sample_names.tsv','sample_names.tsv','_sample_name.tsv','sample_name.tsv']
    #======================================================
    # specify data type of the dataframe
    dtype_dict = {'taxonomy': 'category', param: 'float32'}

    # get the list of the profile files in the input directory
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
    
    ####### Creat a merged ARG profiles ###################################
    # merge data from input files
    df_joined=pd.DataFrame()
    key.append(param) # create the list of columns to output
    for f in files_in: 
        df_sample=pd.read_table(f,header=None, index_col=False, dtype=dtype_dict)
        df_sample.columns=key
        # Remove Eukaryota domain
        df_sample=df_sample[~df_sample[key[0]].str.contains('d__Eukaryota')]
        # calculate the relative abundance
        df_sample[param]=df_sample[param] / df_sample[param].sum()
        # exclude the taxonomy with lower abundance than the threshold
        df_sample=df_sample[df_sample[param] > thres]
        # rename the column with the sequence file name
        sample_name=os.path.basename(f).rstrip(suffix)
        df_sample=df_sample.reindex(columns=key).rename(columns={param: sample_name})
        # merge the sample data
        if df_joined.empty:
            df_joined=df_sample # for the first data file
        else:
            df_joined=pd.merge(df_joined, df_sample, on=key[:-1], how='outer', copy=False)
        print("["+args[0]+"] "+f+" merged.")
        
    # fill out NaN with zero
    df_joined=df_joined.fillna(0)
    # rename the columns as sample name
    if not dic_sample.empty:
        df_joined.rename(columns=dic_sample, inplace=True)
    
    # re-arrange taxonomy columns
    df_taxon=df_joined[key[0]].str.split(r'\|', expand=True, n=5)
    df_taxon=df_taxon.dropna(how='all', axis=1) # remove empty columns
    df_taxon=df_taxon.replace('^[dpcofg]__',r'',regex=True) # remove class initials 
    df_taxon=df_taxon.fillna("-") 
    df_taxon.columns=['domain','phylum','class','order','family','genus']
    df_joined=pd.concat([df_taxon,df_joined],axis=1)
    
    # create the output file
    if not os.path.isdir(dir_out):
        os.makedirs(dir_out) # create the output directory if not existed.
    file_out='_merged.16S.'+param+suffix[-4:]
    df_joined.to_csv(os.path.join(dir_out,file_out),sep='\t',index=False)
    print("["+args[0]+"] "+file_out+" was created.")
    
    # sort the row indices by sum of each row
    df_joined=df_joined.assign(sum=df_joined.sum(axis=1, numeric_only=True))
    df_joined.sort_values('sum',ascending=False, inplace=True)
    df_joined=df_joined.drop('sum', axis=1)
    # Summary table for pca analysis
    df_pca=df_joined.drop(df_taxon.columns, axis=1)
    df_pca=df_pca.set_index('taxonomy').T
    file_out='16S.'+param+'.mpa.csv'
    df_pca.to_csv(os.path.join(dir_out,file_out))
    print("["+args[0]+"] "+file_out+" was created.")
    
    ## summary table for each taxonomy level
    for i, tx in enumerate(df_taxon.columns[1:]):
        df_pca=df_joined[df_joined[tx] != "-"] # remove unclassified rows
        df_pca=df_pca.groupby(tx).sum(numeric_only=True)
        # re-calculate the abundance
        for s in df_pca.columns[1:]:
            df_pca[s]=df_pca[s]/df_pca[s].sum()
        # sort the row indices by sum of each row
        df_pca=df_pca.assign(sum=df_pca.sum(axis=1, numeric_only=True))
        df_pca.sort_values('sum',ascending=False, inplace=True)
        df_pca=df_pca.drop('sum', axis=1).T
        # output the summary table
        file_out='16S.'+param+'.'+str(i+1)+'_'+tx+'.csv'
        df_pca.to_csv(os.path.join(dir_out,file_out),index=True)
        print("["+args[0]+"] "+file_out+" was created.")
    print("["+args[0]+"] completed.")
    
if __name__=="__main__":
    args=sys.argv
    main(args[1],args[2])

