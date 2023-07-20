#========================================================================
# pca_hc_scaling.py / created by Ryo Honda, 2023-06-21
#========================================================================
# This python script perform principal component analysis and hierarchic cluster analysis by:
#	$ python3 pca_hc_scaling.py data.csv dir_out
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

scaling=1 # switch on/off data scaling in pca. True=1, False=0

# get arguments / 引数の読み込み
args=sys.argv
f_in=args[1]; dir_out=args[2]
# Load the CSV file / CSVファイルの読み込み
data_in = pd.read_csv(f_in, header=0, index_col=0)
# Get sample names / サンプル名を取得
samples = data_in.index
params=data_in.columns

## ------ Principal component analysis ------
# Select numeric columns only
#numeric_data = data_in.select_dtypes(include=[np.number])

# Scale the data / データのスケーリング
if scaling:
    scaler = StandardScaler()
    data_num = scaler.fit_transform(data_in)

# Perform PCA
pca = PCA()
pca.fit(data_num)
# 主成分分析の結果を取得
pca_scores = pca.transform(data_num)  # PCスコア（各主成分の値）
explained_variance_ratio = pca.explained_variance_ratio_  # 各主成分の寄与度
loadings = pca.components_.T * np.sqrt(pca.explained_variance_)  # 各成分の負荷量

# --- Output principal components and explained variance ratio as CSV files
# 主成分の値をCSVファイルとして出力
pca_scores_df = pd.DataFrame(pca_scores, columns=['PC{}'.format(i+1) for i in range(pca_scores.shape[0])], index=samples)
# 各成分の寄与率を追加
pca_scores_df.loc['Explained Variance Ratio'] = explained_variance_ratio
# 各成分の累積寄与率を追加
pca_scores_df.loc['Cumulative Variance Ratio'] = np.cumsum(explained_variance_ratio)
# csvファイル書き出し
sfx_scaling = ".scaling" if scaling else ".no_scaling"
file_scores=os.path.join(dir_out,'pca_scores.'+os.path.splitext(os.path.basename(f_in))[0]+sfx_scaling+'.csv')
pca_scores_df.to_csv(file_scores)
print("["+args[0]+"] "+file_scores+" was created.")

# 各成分の負荷量をCSVファイルに出力
loadings_df = pd.DataFrame(loadings, index=params, columns=['PC{}'.format(i+1) for i in range(loadings.shape[1])])
loadings_df=loadings_df.sort_values('PC1',ascending=False)
file_loadings=os.path.join(dir_out, 'pca_loadings.'+os.path.splitext(os.path.basename(f_in))[0]+sfx_scaling+'.csv')
loadings_df.to_csv(file_loadings)
print("["+args[0]+"] "+file_loadings+" was created.")

# Excelファイルに出力
file_xlsx=os.path.join(dir_out, 'pca.'+os.path.splitext(os.path.basename(f_in))[0]+sfx_scaling+'.xlsx')
with pd.ExcelWriter(file_xlsx) as writer:
    pca_scores_df.to_excel(writer, sheet_name='PC scores')
    loadings_df.to_excel(writer, sheet_name='loadings')

## ---PCスコアと各パラメータの寄与度をプロット
#fig, ax = plt.subplots(figsize=(10, 6))

## PCスコアをプロット
#ax.scatter(pca_scores[:, 0], pca_scores[:, 1])
#ax.set_xlabel('PC1')
#ax.set_ylabel('PC2')
#ax.set_title('PCA - PC Scores')

## 各成分の負荷量を矢印で示す
#for i, (x, y) in enumerate(zip(pca_scores[:, 0], pca_scores[:, 1])):
#    ax.text(pca_scores[i, 0] * 1.1, pca_scores[i,1] * 1.1, samples[i], color='black')
#    ax.arrow(0, 0, loadings[0, i] * 10, loadings[1, i] * 10, color='r', alpha=0.5)
#    ax.text(loadings[0, i] * 10.1, loadings[1, i] * 10.1, params[i], color='r')

## Create a folder for saving the plot
#os.makedirs(dir_out, exist_ok=True)
## Save the plot as an image file
#output_path = os.path.join(dir_out, 'pcaplot.'+os.path.splitext(os.path.basename(f_in))[0]+'.png')
#plt.savefig(output_path)

## ------ Hierarchic cluster analysis -----
# Prepare the data / データの準備
X = data_in.iloc[:, 1:].values

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

# Create a folder for saving the dendrogram
os.makedirs(dir_out, exist_ok=True)
# Save the dendrogram as an image file
output_path = os.path.join(dir_out, 'hc.dendrogram.'+os.path.splitext(os.path.basename(f_in))[0]+'.pdf')
plt.savefig(output_path, bbox_inches='tight', pad_inches=0)
plt.show()
