#========================================================================
# pca_hc_scale.py ver.2 / created by Ryo Honda, Last updated: 2025-03-07
#========================================================================
# This python script perform principal component analysis and hierarchic cluster analysis by:
#	$ python3 pca_hc_scale.py data.csv dir_out
#
# [IMPORTANT] 
#  data.csv : data file in csv format (comma-delimited). 
#   In the input csv, each sample data should be contained in each column.

# The scrpit requires: pandas, numpy, scikit-learn, scipy, matplotlib, openpyxl
#------------------------------------------------------------------------------
import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import openpyxl
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from scipy.cluster.hierarchy import dendrogram, linkage

def main(f_in, dir_out):
    scaling=1           # switch on/off data scaling in pca. True=1, False=0
    cal_proportion=0    # switch on/off if you want to convert the values into proportion. True=1, False=0

    # Load the CSV file / CSVファイルの読み込み
    data_in = pd.read_csv(f_in, header=0, index_col=0)
    
    ##------ preparation of data ------
    # Select numeric columns only
    data_in = data_in.select_dtypes(include=[np.number])
    # Get sample names / サンプル名を取得
    samples = data_in.index
    params=data_in.columns
    # convert the values into proportion / 値を割合に変換
    if cal_proportion:
        data_in = data_in.div(data_in.sum(axis=1),axis=0)
    # Scale the data / データのスケーリング
    if scaling:
        scaler = StandardScaler()
        data_num = scaler.fit_transform(data_in)
    else:
        data_num=data_in
    
    ## ------ Principal component analysis ------
    # Perform PCA
    pca = PCA()
    pca.fit(data_num)
    # get the results of pca / 主成分分析の結果を取得
    pca_scores = pca.transform(data_num)  # PC scores / PCスコア（各主成分の値）
    explained_variance_ratio = pca.explained_variance_ratio_  # explained variance of each PC / 各主成分の寄与度
    loadings = pca.components_.T * np.sqrt(pca.explained_variance_)  # PC loadings / 各成分の負荷量
    
    # --- Output principal components and explained variance ratio as CSV files
    # prepare dataframe of PC scores. / 主成分のスコア値をdataframeに用意
    pca_scores_df = pd.DataFrame(pca_scores, columns=['PC{}'.format(i+1) for i in range(len(data_num))], index=samples)
    # add explained variance / 各成分の寄与率を追加
    pca_scores_df.loc['Explained Variance Ratio'] = explained_variance_ratio
    # add cumulative explained variance / 各成分の累積寄与率を追加
    pca_scores_df.loc['Cumulative Variance Ratio'] = np.cumsum(explained_variance_ratio)
    # output PC scores in a csv file / 主成分のスコア値をcsvファイル書き出し
    sfx_scaling = ".scaling" if scaling else ".no_scaling"
    file_scores=os.path.join(dir_out,'pca_scores.'+os.path.splitext(os.path.basename(f_in))[0]+sfx_scaling+'.csv')
    pca_scores_df.to_csv(file_scores)
    print(f"[{args[0]}] {file_scores} was created.")
    
    # output PC loadings in a csv file / 各成分の負荷量をCSVファイルに出力
    loadings_df = pd.DataFrame(loadings, index=params, columns=['PC{}'.format(i+1) for i in range(len(data_num))])
    loadings_df=loadings_df.sort_values('PC1',ascending=False)
    file_loadings=os.path.join(dir_out, 'pca_loadings.'+os.path.splitext(os.path.basename(f_in))[0]+sfx_scaling+'.csv')
    loadings_df.to_csv(file_loadings)
    print(f"[{args[0]}] {file_loadings} was created.")
    
    # Output as an Excel file / Excelファイルに出力
    file_xlsx=os.path.join(dir_out, 'pca.'+os.path.splitext(os.path.basename(f_in))[0]+sfx_scaling+'.xlsx')
    with pd.ExcelWriter(file_xlsx) as writer:
        pca_scores_df.to_excel(writer, sheet_name='PC scores')
        loadings_df.to_excel(writer, sheet_name='loadings')
    
    ## ------ Hierarchic cluster analysis -----
    # Prepare the data / データの準備
    X = data_in.values
    # perform hierarchical clustering analysis / 階層クラスタリング解析を実行
    Z = linkage(X, method='ward')
    
    # Plot the dendrogram / 樹形図をプロット
    plt.figure(figsize=(10, 5))
    dendrogram(Z, labels=samples, color_threshold=0, above_threshold_color='k',orientation='left')
    plt.xlabel('Distance')
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['left'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.subplots_adjust(right=0.7)
    
    # Create a folder for saving the dendrogram / 出力用のディレクトリを作成
    os.makedirs(dir_out, exist_ok=True)
    # Save the dendrogram as a pdf file / 樹形図をpdfに保存
    output_path = os.path.join(dir_out, 'hc.dendrogram.'+os.path.splitext(os.path.basename(f_in))[0]+'.pdf')
    plt.savefig(output_path, bbox_inches='tight', pad_inches=0)
    plt.show()
    
if __name__=="__main__":
    args=sys.argv
    main(args[1],args[2])

