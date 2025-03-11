#========================================================================
# curate_drug_class.py / created by Ryo Honda, Last updated: 2025-03-11
#========================================================================
# This python script creates a summary table of sequence reads by:
#	$ python3 crtDrugClass.py summary_reads.csv dir_arg dir_out
#
# [IMPORTANT] 
#  summary_reads.csv = summary file containing RPK of total ARGs and 16S of each sample, created by summary_reads_16S_ARG.py or summary_reads_16S_ARG_MGE.py 
#  dir_arg = directory of ARG profiles created by make_ARG_profile.py
#  dir_out = directory for output files. (Specify as "." for the current diretory)
#
#  If the dir_in contains "_sample_names.tsv", which has sequence names in the first column and  sample names in the second column, the columns of the output table are labelled with the corresponding sample names.
#------------------------------------------------------------------------------
import glob
import os
import sys
import pandas as pd
import openpyxl

def main(file_reads, dir_in, dir_out):
    #======= Parameters settings ========================
    # set "True" if you want to output drug class profile of each sample
    output_each_profile=False
    # suffix of the sample data file. 
    suffix='.ARG_profile.tsv'
    # name of the file containing list of sample names corresponding to sequence names
    # (list of allowable file names for compatibility)
    file_sample_name=['_sample_names.tsv','sample_names.tsv','_sample_name.tsv','sample_name.tsv']
    
    ####### Main ###################################    
    # get the list of the profile files in the input directory
    files_in=sorted(glob.glob(os.path.join(dir_in,"*"+suffix)))
    n=len(files_in)
    print("["+args[0]+"] started. ("+str(n)+" profiles are found.)")
    
    # create the output directory if not existed.
    if not os.path.isdir(dir_out):
        os.makedirs(dir_out) 
    
    # create the dictionary of sample names corresponding to sequence names
    dic_sample=pd.DataFrame() # create an empty dataframe
    for f in file_sample_name:
        if os.path.isfile(os.path.join(dir_in,f)):
            dic_sample=pd.read_table(os.path.join(dir_in,f), header=None, index_col=0).squeeze(axis=1)
            dic_sample=dic_sample.to_dict()
            break
        else:
            print("Warning: _sample_names.csv is not found.")
            dic_sample={}
    #======================================================
    
    ## Summarize ARG profile data
    df_sum_all=pd.DataFrame() # for summation of all ARGs
    df_sum_single=pd.DataFrame() # for summation of single-drug ARGs
    df_sum_nonefflux=pd.DataFrame() # for summation of non-efflux ARGs
    df_sum_efflux=pd.DataFrame() # for summation of efflux ARGs
    
    for i, f in enumerate(files_in, 1): 
        # get sequence name
        sample_name=os.path.basename(f).replace(suffix,'')
        # get the ratio of total RPK to RPK-16S
        df_reads=pd.read_csv(file_reads, header=0, index_col=0)
        rrpk=df_reads.at[sample_name,'Total RPK of ARG'] / df_reads.at[sample_name,'Total RPK of 16S']
        
        ##======= append drug-class counts in ARG profile data =======
        f_arg=os.path.join(dir_in,sample_name+suffix)
        df_arg=pd.read_table(f_arg, header=0)
        # convert unit in ARG profiles from RPK/RPK-total to RPK/RPK-16S
        df_arg['RPK/RPK-16S']=df_arg['prop_RPK'] * rrpk
        # create the list of drug classes
        drugs=df_arg['Drug Class'].str.split(';')
        drugs=drugs.explode()
        drugs=pd.unique(drugs.dropna())
        
        # count RPK/RPK-16S of each drug class included in each ARG
        df_arg[drugs]=0
        for d in drugs:
            rows_d=df_arg['Drug Class'].str.contains(d)
            df_arg.loc[rows_d, d]=df_arg.loc[rows_d, 'RPK/RPK-16S']
        
        # output the new profile data to a tsv file
        if output_each_profile:
            file_out=os.path.join(dir_out, sample_name+'.ARG_profile+drug_class.tsv')
            df_arg.to_csv(file_out,sep='\t',index=False)
            print("["+args[0]+"] "+file_out+" was created. ("+str(i)+"/"+str(n)+")")
        
        ##======= create the aggregated summary file ============================
        # calculate the sum for each drug class, and merge in the summary file
        # summation of all ARGs
        sum_drug=df_arg[drugs].sum()
        sum_drug=sum_drug.rename(sample_name)
        df_sum_all=pd.merge(df_sum_all, sum_drug, how='outer', left_index=True, right_index=True)
        # summation of single-drug ARGs
        df_tmp=df_arg[df_arg['MAR'] == 1]
        sum_drug=df_tmp[drugs].sum()
        sum_drug=sum_drug.rename(sample_name)
        df_sum_single=pd.merge(df_sum_single, sum_drug, how='outer', left_index=True, right_index=True)
        # summation of non-efflux ARGs
        df_tmp=df_arg[df_arg['Resistance Mechanism'] != 'efflux']
        sum_drug=df_tmp[drugs].sum()
        sum_drug=sum_drug.rename(sample_name)
        df_sum_nonefflux=pd.merge(df_sum_nonefflux, sum_drug, how='outer', left_index=True, right_index=True)
        # summation of efflux ARGs
        df_tmp=df_arg[df_arg['Resistance Mechanism'] == 'efflux']
        sum_drug=df_tmp[drugs].sum()
        sum_drug=sum_drug.rename(sample_name)
        df_sum_efflux=pd.merge(df_sum_efflux, sum_drug, how='outer', left_index=True, right_index=True)
    
    # rename the sequence names in the columns to sample names
    if dic_sample:
        df_sum_all.rename(columns=dic_sample, inplace=True)
        df_sum_single.rename(columns=dic_sample, inplace=True)
        df_sum_nonefflux.rename(columns=dic_sample, inplace=True)
        df_sum_efflux.rename(columns=dic_sample, inplace=True)
    df_sum_all=df_sum_all.fillna(0).sort_index().T
    df_sum_single=df_sum_single.fillna(0).sort_index().T
    df_sum_nonefflux=df_sum_nonefflux.fillna(0).sort_index().T
    df_sum_efflux=df_sum_efflux.fillna(0).sort_index().T
    
    # output the aggregated summary to an Excel file
    if "openpyxl" in sys.modules:
        file_xlsx=os.path.join(dir_out, 'ARG.drug_class.per_16S.xlsx')
        with pd.ExcelWriter(file_xlsx) as writer:
            df_sum_all.to_excel(writer, sheet_name='all')
            df_sum_single.to_excel(writer, sheet_name='single_ARG')
            df_sum_nonefflux.to_excel(writer, sheet_name='non-efflux')
            df_sum_efflux.to_excel(writer, sheet_name='efflux')
        print("["+args[0]+"] "+file_xlsx+" was created.")
    else:
        # output the aggregated summary to csv files
        file_out=os.path.join(dir_out, 'ARG.drug_class.all.per_16S.csv')
        df_sum_all.to_csv(file_out)
        print("["+args[0]+"] "+file_out+" was created.")
        file_out=os.path.join(dir_out, 'ARG.drug_class.single.per_16S.csv')
        df_sum_single.to_csv(file_out)
        print("["+args[0]+"] "+file_out+" was created.")
        file_out=os.path.join(dir_out, 'ARG.drug_class.nonefflux.per_16S.csv')
        df_sum_nonefflux.to_csv(file_out)
        print("["+args[0]+"] "+file_out+" was created.")
        file_out=os.path.join(dir_out, 'ARG.drug_class.efflux.per_16S.csv')
        df_sum_efflux.to_csv(file_out)
        print("["+args[0]+"] "+file_out+" was created.")

if __name__=="__main__":
    args=sys.argv
    main(args[1],args[2],args[3])

