#========================================================================
# summary_reads_16S_ARG_MGE.py / created by Ryo Honda, 2023-07-15
#========================================================================
# This python script creates a summary table of sequence reads by:
#	$ python3 summary_reads_16S_ARG.py dir_qt dir_16s dir_arg dir_mge dir_out
#
# [IMPORTANT] 
#  dir_qt = directory of json reports of quality trimming by fastp
#  dir_16s = directory of 16S read counts data by kraken2 & bracken (which is classified by '(D)omain')
#  dir_arg = directory of ARG profiles created by make_ARG_profile.py
#  dir_mge = directory of MGE profiles created by make_MGE_profile.py
#  dir_out = directory to output the summary file. (Specify as "." for the current diretory)
#------------------------------------------------------------------------------
import glob
import json
import os
import sys
import pandas as pd

#======= Parameters settings ========================
## suffix of the report files
sfx_qt=".report.json"
sfx_16s="_0.domain.tsv"
sfx_arg=".ARG_profile.tsv"
sfx_mge=".MGE_profile.tsv"

## length of 16S rRNA gene of E.coli for calculation of RPK
## [Reference] Brosius, J., Palmer, M. L., Kennedy, P. J. & Noller, H. F. Complete nucleotide sequence of a 16S ribosomal RNA gene from Escherichia coli. Proc. Natl. Acad. Sci. 75, 4801–4805 (1978).
slen_16s=1541

####### Summarize sequence reads data ###################################
## import arguments 引数処理
args=sys.argv
dir_qt=args[1]; dir_16s=args[2]; dir_arg=args[3]; dir_mge=args[4]; dir_out=args[5]
## get the file list in the input directory
files_qt=sorted(glob.glob(os.path.join(dir_qt,"*"+sfx_qt)))
i=0; n=len(files_qt)
print("["+args[0]+"] started. ("+str(i)+"/"+str(n)+")")

## create an output dataframe
df_out=pd.DataFrame(columns=['Raw sequence reads','Quality sequence reads','Total 16S reads','Total ARG reads','Total MGE reads','Total RPK of 16S','Total RPK of ARG','Total ARG abundance (RPK/RPK-16S)','Total RPK of MGE','Total MGE abundance (RPK/RPK-16S)'])

## Summarize reads data
for f in files_qt: 
    # get sequence name
    name=os.path.basename(f).rstrip(sfx_qt) 
    # get total reads before and after quality trimming
    with open(f, 'r') as data:
        js_qt=json.load(data)
    reads_before_qt=js_qt['summary']['before_filtering']['total_reads']
    reads_after_qt=js_qt['summary']['after_filtering']['total_reads']
    # get total 16S reads and RPK
    f_16s=os.path.join(dir_16s,name+sfx_16s)
    if os.path.isfile(f_16s):
        df_16s=pd.read_table(f_16s, header=0,index_col=0)
        reads_16s=df_16s.loc['Bacteria','new_est_reads']
        rpk_16s=reads_16s / slen_16s * 1000
    else:
        reads_16s=None; rpk_16s=None
    # get total ARG reads and RPK
    f_arg=os.path.join(dir_arg,name+sfx_arg)
    if os.path.isfile(f_arg):
        df_arg=pd.read_table(f_arg, header=0)
        reads_arg=df_arg['reads'].sum()
        rpk_arg=df_arg['RPK'].sum()
        if rpk_16s is not None and rpk_16s != 0:
            abundance_arg=rpk_arg / rpk_16s
        else:
            abundance_arg=None
    else:
        reads_arg=None; reads_arg=None; abundance_arg=None
    # get total MGE reads and RPK
    f_mge=os.path.join(dir_mge,name+sfx_mge)
    if os.path.isfile(f_mge):
        df_mge=pd.read_table(f_mge, header=0)
        reads_mge=df_mge['reads'].sum()
        rpk_mge=df_mge['RPK'].sum()
        if rpk_16s is not None and rpk_16s != 0:
            abundance_mge=rpk_mge / rpk_16s
        else:
            abundance_mge=None
    else:
        reads_mge=None; reads_mge=None; abundance_mge=None
        
    # add in the output dataframe
    df_out.loc[name]=[reads_before_qt, reads_after_qt,reads_16s,reads_arg,reads_mge,rpk_16s,rpk_arg,abundance_arg,rpk_mge,abundance_mge]
    i+=1
    print("["+args[0]+"] "+name+" is done.("+str(i)+"/"+str(n)+")")

# output to a csv file
if not os.path.isdir(dir_out):
    os.makedirs(dir_out) # create the output directory if not existed.
file_out=os.path.join(dir_out,'_summary_reads_16S_ARG_MGE.csv')
df_out.to_csv(file_out,sep=',',index=True)
print("["+args[0]+"] "+file_out+" was created. ("+str(i)+"/"+str(n)+")")
