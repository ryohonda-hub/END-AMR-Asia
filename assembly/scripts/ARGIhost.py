#==============================================================================
# ARGIhost.py ver.2 / created by Ryo Honda, Last updated: 2025-10-22
#==============================================================================
# This python script creates a ARG-host table from ARG and taxonomy profiles of a contig by:
#	$ python3 ARGIhost.py kraken2_index rgi_results.txt taxonomy_results.tsv dir_out
#
#  kraken2_index = the Inspect file of Kraken2 index. (refer https://benlangmead.github.io/aws-indexes/k2)
#  rgi_results.txt = RGI output file for ARG
#  taxonomy_results.tsv = Kraken taxonomy table
#  dir_out = output directory, or use "." to specify the current directory
#
# The output file with suffix of ".ARGhost.tsv" is created in the output directory. 
# The output file contains gene information and proportion of RPK.
#------------------------------------------------------------------------------
import os
import sys
import pandas as pd

def main(idx_file, arg_file, taxon_file, dir_out):
    #====== load taxonomy file ======
    df_tax = pd.read_csv(taxon_file, sep='\t', header=None, names=['status', 'contig_id', 'taxid', 'kmer_count', 'classification_info'])
    # load taxonomy index of kraken2 and create the dictionary
    dict_taxon= {}
    with open(idx_file, 'r') as f:
        for i, line in enumerate(f):
            if i<7: continue    # ignore header lines
            cols = line.strip().split('\t')
            _, _, _, rank, taxid, name_txt = cols
            dict_taxon[int(taxid.strip())] = [rank, name_txt.strip()]
                
    ## create contig-taxonomy table
    df_tax = df_tax[df_tax['status'] == 'C'].copy() # Remove Unclassified (U)
    df_tax['contig_id'] = df_tax['contig_id'].str.strip()
    # convert taxid into taxonomy names
#    df_tax['Species'] = df_tax['taxid'].map(dict_taxon)
    df[['rank', 'taxonomy']] = pd.DataFrame(df['taxid'].map(dict_taxon).tolist(), index=df.index)
    df_tax['Species'] = df_tax['Species'].fillna('Unknown')
    
    # output the contig-taxonomy table
    if not os.path.exists(dir_out):
        os.makedirs(dir_out)
    df_out = df_tax[['contig_id', 'Species']]
    file_out=os.path.splitext(os.path.basename(taxon_file))[0]
    file_out=os.path.join(dir_out,file_out.removesuffix('.kraken2')+".ctg.taxon.tsv")
    df_out.to_csv(file_out, sep='\t', index=False)

    #====== load RGI results ======
    df_arg = pd.read_csv(arg_file, sep='\t',header=0)
    # extract gene symbol and contig ID
    df_arg.insert(0,'contig_id','')
    df_arg.insert(1,'gene symbol','')
    df_arg['gene symbol']=df_arg['Best_Hit_ARO'].str.split(' ').str[0]
    df_arg['contig_id']=df_arg['Contig'].str.replace(r'_[0-9 ]+$','',regex=True)
    
    #======= Merge ARG and taxonomy =====
    df_merged = pd.merge(df_arg, df_tax, on='contig_id', how='left')
    df_merged['Species'] = df_merged['Species'].fillna('Unknown')
    
    # output the ARG-taxonomy table (full-list)
    file_out=os.path.splitext(os.path.basename(arg_file))[0]
    file_out=os.path.join(dir_out,file_out.removesuffix('.ARGI')+".ctg.ARGhost.tsv")
    df_merged.to_csv(file_out, sep='\t', index=False)
    
    #===== ARG-host summary table =====
    df_arg_sum=df_arg.groupby('contig_id')['gene symbol'].apply(lambda x: ','.join(sorted(set(x)))).reset_index()
    df_arg_sum.columns = ['contig_id', 'ARGs']
    df_arg_sum['num_ARG']=df_arg_sum['ARGs'].str.count(',')+1
    
    #======= Merge ARG and taxonomy =====
    df_merged = pd.merge(df_arg_sum, df_tax, on='contig_id', how='left')
    df_merged['Species'] = df_merged['Species'].fillna('Unknown')

    # output the ARG-taxonomy table (summary)
    file_out=os.path.splitext(os.path.basename(arg_file))[0]
    file_out=os.path.join(dir_out,file_out.removesuffix('.ARG_profile')+".ctg.ARGhost.summary.tsv")
    df_merged.to_csv(file_out, sep='\t', index=False)

if __name__ == "__main__":
    args=sys.argv
    if len(args) != 5:
        print(f"Usage: python {sys.argv[0]} taxonomy_names.dump rgi_results taxonomy_results dir_out")
        sys.exit(1)
    main(args[1],args[2],args[3],args[4])

