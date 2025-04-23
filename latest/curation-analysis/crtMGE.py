#========================================================================
# curate_MGE_profile.py ver.2 / created by Ryo Honda, Last updated: 2025-04-23
#========================================================================
# This python script creates a profile comparison table of multiple samples by:
#	$ python3 crtMGE.py dir_in dir_out
#
# [IMPORTANT] The input file must be ".tsv" files
#  dir_in = directory of sample data files. (all .tsv files in the directory will be merged.)
#  dir_out = directory to output the data file

#  If the dir_in contains "_sample_names.tsv", which has sequence names in the first column and sample names in the second column, the columns of the output table are labelled with the corresponding sample names.
#------------------------------------------------------------------------------
import glob
import os
import sys
import pandas as pd

def main(dir_in, dir_out):
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
    # specify data type of the dataframe
    dtype_dict = {
    'sseqid': 'category',
    'MGEDB ID': 'category',
    'gene symbol': 'category',
    'Function': 'category',
    'slen': 'int32',
    param: 'float32'
    }
    ####### Creat a merged ARG profiles ###################################    
    n_warn=0 # count warning
    # get the list of sequence files in the input directory
    files_in=sorted(glob.glob(os.path.join(dir_in,"*"+suffix)))
    # create the dictionary of sample names corresponding to sequence names
    dic_sample=pd.DataFrame() # create an empty dataframe
    for f in file_sample_name:
        if os.path.isfile(os.path.join(dir_in,f)):
            print(f"'{f}' is found. The results will be output with sample names. ")
            dic_sample=pd.read_table(os.path.join(dir_in,f), header=None, index_col=0).squeeze(axis=1)
            # warn if any duplicated sample names
            duplicate_names = dic_sample[dic_sample.duplicated(keep=False)]
            if not duplicate_names.empty:
                print("Warning: There are duplicated sample names. Check the list in 'duplicated_sample_names.csv'")
                duplicate_names.sort_values().to_csv(os.path.join(dir_out,'duplicated_sample_names.csv'), header=False)
                n_warn+=1
            # convert the list to the dict type.
            dic_sample=dic_sample.to_dict()
            break
        else:
            dic_sample={}
    
    # merge data from input files
    df_joined=pd.DataFrame()
    key.append(param) # create the list of columns to output
    for f in files_in: 
        df_sample=pd.read_table(f,header=0, dtype=dtype_dict)
        sample_name=os.path.basename(f).removesuffix(suffix)
        df_sample=df_sample.reindex(columns=key).rename(columns={param: sample_name})
        #df_sample = df_sample.drop_duplicates(subset=key[:-1])  # remove duplicated data
        if df_joined.empty:
            df_joined=df_sample # for the first data file
            n_file=0
        else:
            df_joined=pd.merge(df_joined, df_sample, on=key[:-1], how='outer')
            n_file+=1
        print(f"[{args[0]}] {f} merged.")
        
    # fill out NaN with zero (for the numeric columns only)
    df_joined[df_joined.select_dtypes(include=['number']).columns] = df_joined.select_dtypes(include=['number']).fillna(0)

    # sort the row indices by sum of each row
    df_joined=df_joined.assign(sum=df_joined.sum(axis=1, numeric_only=True))
    df_joined.sort_values('sum',ascending=False, inplace=True)
    df_joined=df_joined.drop('sum', axis=1)
    
    # rename the columns as sample name
    if dic_sample:
        df_joined_out=df_joined.rename(columns=dic_sample, inplace=False)
    else:
        df_joined_out=df_joined
    # create the output file
    if not os.path.isdir(dir_out):
        os.makedirs(dir_out) # create the output directory if not existed.
    file_out='_merged'+suffix[0:-4]+'.'+param+suffix[-4:]
    df_joined_out.to_csv(os.path.join(dir_out,file_out),sep='\t',index=False)
    print(f"[{args[0]}] Merged profile was created.")
    ######### Summation by args categories ###############################
    # list of the sample columns to aggregate
    col_sum=df_joined.columns[len(key)-len(df_joined.columns)-1:]
    # sum up by category
    for cat in cats:
        # sum up by the key 集計
        df_sum=df_joined.groupby(cat,observed=False)[col_sum].sum() 
        # sort the row indices by sum of each row
        df_sum=df_sum.assign(sum=df_sum.sum(axis=1, numeric_only=True))
        df_sum.sort_values('sum',ascending=False, inplace=True)
        df_sum=df_sum.drop('sum', axis=1)
        
        # rename the columns as sample name
        if dic_sample:
            df_sum.rename(columns=dic_sample, inplace=True)
        # arrange the table for output
        df_sum=df_sum.sort_index(axis=1,ascending=False).T
        # output the summary file.
        file_out="MGE."+cat.replace(' ','_')+"."+param+".csv"
        df_sum.to_csv(os.path.join(dir_out,file_out))
        print(f"[{args[0]}] {file_out} was created.")
    print(f"{args[0]} completed. {n_file} files were curated. There are {n_warn} warnings.")

if __name__=="__main__":
    args=sys.argv
    main(args[1],args[2])

