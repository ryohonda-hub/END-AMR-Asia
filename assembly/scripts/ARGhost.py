#==============================================================================
# ARGhost.py ver.1.0 / created by Ryo Honda, Last updated: 2025-04-29
#==============================================================================
# This python script creates a ARG-host table from ARG and taxonomy profiles of a contig by:
#	$ python3 ARGIhost.py names.dump blast_results.txt taxonomy_results.tsv dir_out
#
#  names.dump = NCBI taxdump database 
#  blast_results.txt = blast output file for ARG
#  taxonomy_results.tsv = Kraken taxonomy table
#  dir_out = output directory, or use "." to specify the current directory
#
# The output file with suffix of ".ARG-host.tsv" is created in the output directory. 
# The output file contains gene information and proportion of RPK.
#------------------------------------------------------------------------------
import os
import re
import sys
import pandas as pd

def main(names_dmp_file, arg_file, taxon_file, dir_out):
    #=== cutoff thresholds of the blast result ===
    th_mlen=100     # (bp) hits with matched length below this value are excluded.
    th_pident=90    # (%) hits with pident below this value are excluded.
    #=============================================
    
    #====== load taxonomy file ======
    df_tax = pd.read_csv(taxon_file, sep='\t', header=None, names=['status', 'contig_id', 'taxid', 'kmer_count', 'classification_info'])
    # load taxonomy names.dump
    dict_taxon= {}
    with open(names_dmp_file, 'r') as f:
        for line in f:
            cols = line.strip().split('\t|\t')
            taxid, name_txt, _, name_class = cols
            if name_class.strip() == "scientific name":
                dict_taxon[int(taxid.strip())] = name_txt.strip()
                
    ## create contig-taxonomy table
    df_tax = df_tax[df_tax['status'] == 'C'].copy() # Remove Unclassified (U)
    df_tax['contig_id'] = df_tax['contig_id'].str.strip()
    # convert taxid into taxonomy names
    df_tax['Species'] = df_tax['taxid'].map(dict_taxon)
    df_tax['Species'] = df_tax['Species'].fillna('Unknown')
    
    # output the contig-taxonomy table
    if not os.path.exists(dir_out):
        os.makedirs(dir_out)
    df_out = df_tax[['contig_id', 'Species']]
    file_out=os.path.splitext(os.path.basename(taxon_file))[0]
    file_out=os.path.join(dir_out,file_out.removesuffix('.kraken2')+".ctg.taxon.tsv")
    df_out.to_csv(file_out, sep='\t', index=False)

    #====== load blast results ======
    df_arg = pd.read_csv(arg_file, sep='\t',header=None)
    # exclude the blast results with pident or length below the cutoff thresholds
    df_arg.columns=['contig_id', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'qlen', 'slen']
    df_arg=df_arg.query('pident >= @th_pident & length >= @th_mlen').copy()

    # extract ARO accession and gene symbol from the sseqid 
    df_arg.insert(1,'ARO Accession','')
    df_arg.insert(2,'gene symbol','')
    for i, gene in df_arg.iterrows():
        sseqid=gene.sseqid
        df_arg.at[i, 'ARO Accession']= re.split(r'\|', sseqid)[4] # get ARO accession ID
        str = re.split(r'\|', sseqid)[5] # get gene symbol 
        # remove allele prefix, which is numbers and alphabets in lower cases followed by a hyphen '-' or a period '.' .
        df_arg.at[i, 'gene symbol'] = re.sub(r'[-|.][0-9a-z]+$','',str)

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
    file_out=os.path.join(dir_out,file_out.removesuffix('.ARGI')+".ctg.ARGhost.summary.tsv")
    df_merged.to_csv(file_out, sep='\t', index=False)

if __name__ == "__main__":
    args=sys.argv
    if len(args) != 5:
        print(f"Usage: python {sys.argv[0]} taxonomy_names.dump rgi_results taxonomy_results dir_out")
        sys.exit(1)
    main(args[1],args[2],args[3],args[4])

