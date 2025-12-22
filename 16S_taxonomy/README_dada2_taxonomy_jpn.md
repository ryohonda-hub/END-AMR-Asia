# 16S rRNA V3–V4（341f/805r）ASV解析 NIGスパコン実行マニュアル
 
**対象:** SILVA v138.2 + DADA2による16S rRNA（V3–V4領域）ASV解析

---

##  目次

1. [解析の全体像](#1-解析の全体像)
2. [事前準備（必須）](#2-事前準備必須)
3. [プロジェクト構成](#3-プロジェクト構成)
4. [スクリプトの準備と設定](#4-スクリプトの準備と設定)
5. [ジョブの投入と実行](#5-ジョブの投入と実行)
6. [結果の確認と後処理の実行](#6-結果の確認と後処理の実行)
7. [パラメータ調整ガイド](#7-パラメータ調整ガイド)
8. [出力ファイル一覧](#8-出力ファイル一覧)
9. [SILVAデータベースのダウンロード_先生用](#9-SILVAデータベースのダウンロード_先生用)

---

## 1. 解析の全体像

### 1-1. DADA2解析パイプラインの流れ

```
FASTQ配置 → メイン解析(js1) → 結果確認 → 後処理(js2) → 完了
```

**メイン解析で実行される処理:**
1. FASTQ読み込み（R1/R2自動検出）
2. 品質プロット作成（PDF出力）
3. フィルタリング（低品質領域除去、プライマー除去）
4. エラー学習（learnErrors）
5. ASV推定（dada）
6. ペアエンドマージ（mergePairs）
7. キメラ除去（removeBimeraDenovo）
8. ASVテーブル作成（RDS保存）
9. タクソノミー付与（assignTaxonomy）
10. Read tracking（全工程のread数記録）

**後処理で実行される処理:**
- ASV配列のFASTA出力
- カウントテーブル・タクソノミーテーブル作成
- 相対存在量（%）計算
- 階層別集計（Kingdom～Genus）

### 1-2. 使用ツール

- **解析手法:** DADA2
- **データベース:** SILVA v138.2
- **実行環境:** NIGスパコン

---

## 2. 事前準備（必須）

### 2-1.  FASTQファイルの配置（最重要）

**手順:**

1. 自分のパソコンで、生物技研から返却された **「16S解析結果フォルダ」** を開く
2. その中の **`raw_fastq/`** フォルダを開く
3. 中にある **すべての `*.fastq.gz` ファイル** を、NIGスパコンの以下のディレクトリにコピー:

```bash
~/16S_SILVA/raw_fastq/
```

**注意事項:**
-  すべてのFASTQファイルを **同じ `raw_fastq/` フォルダ** に配置すること
-  サンプルごと、実験ごとにフォルダを分けない
-  追加データが来ても **同じフォルダに追加する**
-  ファイル名は生物技研のままでOK（`*_R1_001.fastq.gz`, `*_R2_001.fastq.gz`）

### 2-2. NIGスパコンへのログイン

```bash
# Gateway nodeにログイン
↓

# 計算ノードへログイン
ex：ssh a001
```

### 2-3. 必要なディレクトリの作成

```bash
# プロジェクトルートの作成
mkdir -p ~/16S_SILVA/raw_fastq
mkdir -p ~/16S_SILVA/scripts

# ログディレクトリの作成
mkdir -p ~/log
```

** 注意:** `~/log/` を作成しないとジョブ投入時にエラーになります。

---

## 3. プロジェクト構成

### 3-1. ディレクトリ構造

**構成例:**

```
~/16S_SILVA/                    ← プロジェクトルート
├── raw_fastq/                  ← FASTQファイル置き場
│   ├── Sample1_R1_001.fastq.gz
│   ├── Sample1_R2_001.fastq.gz
│   ├── Sample2_R1_001.fastq.gz
│   └── Sample2_R2_001.fastq.gz
├── scripts/                    ← スクリプト類
│   ├── s1_dada2_pipeline.r
│   ├── s2_asv_postprocess.r
│   ├── js1_dada2_pipeline.sh
│   └── js2_asv_postprocess.sh
└── output/                     ← 解析結果（自動作成）
    ├── filtered/               ← フィルタリング済みFASTQ
    ├── ASV_nochim.rds
    ├── taxonomy.rds
    └── ... (その他出力ファイル)
```

**SILVAデータベースの場所(例えば。先生のSILVAのdbの場所が決まったらパスを変更):**
```
/home/ryohonda/db/dada2/SILVA_138.2/silva_nr99_v138.2_toGenus_trainset.fa.gz
```
（本多先生のディレクトリ内）

### 3-2. プロジェクトルートとは

**プロジェクトルート** = 解析に必要な素材がまとまっている最上位フォルダのこと。

スクリプト内の以下の部分を、自分のプロジェクトルートに合わせて変更します:

```r
PROJECT_DIR <- "/home/username/16S_SILVA"
```

---

## 4. スクリプトの準備と設定

### 4-1. スクリプトファイルの配置

以下の4つのファイルを `~/16S_SILVA/scripts/` に配置:

1. **s1_dada2_pipeline.r** - メイン解析用Rスクリプト
2. **s2_asv_postprocess.r** - 後処理用Rスクリプト
3. **js1_dada2_pipeline .sh** - メイン解析ジョブ投入用シェルスクリプト
4. **js2_asv_postprocess.sh** - 後処理ジョブ投入用シェルスクリプト

### 4-2. Rスクリプトの編集

#### s1_dada2_pipeline.r の編集箇所

**1. プロジェクトディレクトリの指定:**

```r
# ★自分のユーザー名に変更
PROJECT_DIR <- "/home/username/16S_SILVA"
FASTQ_DIR   <- file.path(PROJECT_DIR, "raw_fastq")
OUT_DIR     <- file.path(PROJECT_DIR, "output")
```

**2. SILVAデータベースのパス指定:**
```r
# ★本多先生のパスに変更（または自分でダウンロードした場所）
silva_file <- "/home/ryohonda/db/dada2/SILVA_138.2/silva_nr99_v138.2_toGenus_trainset.fa.gz"   # ★ここを先生のdbのパスに変更
```

**3. フィルタリングパラメータ（必要に応じて調整）:**

```r
out <- filterAndTrim(
    fwd = fnFs, filt = filtFs,
    rev = fnRs, filt.rev = filtRs,
    trimLeft = c(17, 21),      # プライマー長（V3-V4: 341F=17bp, 805R=21bp）
    truncLen = c(280, 220),    # 品質プロットを見て調整
    maxN = 0, maxEE = c(2, 2), 
    truncQ = 2, compress = TRUE, multithread = TRUE
)
```

#### s2_asv_postprocess.r の編集箇所

**出力ディレクトリの指定:**

```r
# ★自分のユーザー名に変更
OUT_DIR <- "/home/username/16S_SILVA/output"
```

### 4-3. ジョブ投入スクリプトの作成

#### js1_dada2_pipeline.sh（メイン解析用）
それなりにリソースをつかうので，medium か rome 推奨。
8コアとメモリ32GB要求（サンプル数が100を超えるようなら 64GB要求。ただしキュー待ち時間は長くなる）
```bash
#!/bin/bash
#SBATCH --output=/home/username/log/%x_%j.out.log    # ★ユーザー名変更
#SBATCH --error=/home/username/log/%x_%j.err.log     # ★ユーザー名変更
#SBATCH -p medium　　 # 必要に応じて変更
#SBATCH -N 1-1
#SBATCH -n 8
#SBATCH --mem=32G         　   # 45サンプルなら十分(必要に応じて変更)

# スクリプト実行（★パス変更）
Rscript /home/username/16S_SILVA/scripts/s1_dada2_pipeline.r
```

#### js2_asv_postprocess.sh（後処理用）
軽いのでshortでOK。自分のPCでもOK。
```bash
#!/bin/bash

#SBATCH --output=/home/username/log/%x_%j.out.log    # ★ユーザー名変更
#SBATCH --error=/home/username/log/%x_%j.err.log     # ★ユーザー名変更
#SBATCH -p short

# スクリプト実行（★パス変更）
Rscript /home/username/16S_SILVA/scripts/s2_asv_postprocess.r
```

### 4-4. 実行権限の付与

```bash
cd ~/16S_SILVA/scripts
chmod u+x js1_dada2_pipeline.sh
chmod u+x js2_asv_postprocess.sh
```

---

## 5. ジョブの投入と実行

### 5-1. メイン解析の投入

```bash
cd ~/16S_SILVA/scripts
sbatch js1_dada2_pipeline.sh
```

**成功時の表示:**
```
Submitted batch job 123456
```

この `123456` がジョブIDです。

### 5-2. ジョブの状態確認

```bash
# 自分のジョブの状態を確認
squeue -u username
```

**出力例:**
```
JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
123456    medium dada2_me username  R       1:23      1 node01
```

**ステータスの意味:**
- `R` (Running): 実行中
- `PD` (Pending): 待機中（リソース空き待ち）
- 何も表示されない: 完了または終了

### 5-3. ログのリアルタイム確認

```bash
# 標準出力ログ（解析の進行状況）
tail -f ~/log/dada2_medium_123456.out.log

# エラーログ（エラーメッセージ）
tail -f ~/log/dada2_medium_123456.err.log
```

**終了方法:** `Ctrl + C`

### 5-4. ジョブのキャンセル（必要時）

```bash
scancel 123456
```

---

## 6. 結果の確認と後処理の実行

### 6-1.  メイン解析完了後の必須確認項目

**後処理（js2）を実行する前に、必ず以下を確認してください:**

####  確認1: 品質プロットの確認

-R1_quality_before.pdf 
-R2_quality_before.pdf

フィルタリング後（効果確認用）:

-R1_quality_after.pdf
-R2_quality_after.pdf

**PDFをダウンロードして確認:**

**確認ポイント:**
- R1とR2の品質スコア（緑の帯）がQ30以上を保っている領域を確認
- `truncLen` で指定した位置が適切か判断

####  確認2: Read trackingの確認

```bash
cat read_tracking.tsv
```

**確認ポイント:**
- 各ステップでreadが極端に減っていないか
- 特に `merged` と `nonchim` のread数が input よりだいぶ減っていないか
- **`final_perc_reads_retained`（最終残存率）を計算してあります**
  - この値は「最初のリード数に対して、最終的に何％残ったか」を示します
  - 計算式: `(nonchim / input) × 100`


**良い例:**
```
           input filtered denoisedF denoisedR merged nonchim final_perc_reads_retained
Sample1   150000   140000    135000    130000  90000   88000  58.7
Sample2   120000   110000    105000    100000  70000   68000  56.7
```
→ 残存率が50%以上あり、各ステップで極端な減少なし

**悪い例（要再調整）:**
```
           input filtered denoisedF denoisedR merged nonchim final_perc_reads_retained
Sample1   150000   140000    135000    130000   5000    4800   3.2
Sample2   120000   110000    105000    100000   3000    2900   2.4
```

```
実際の作業フロー:
1回目の解析:
  ① beforeプロット確認 → truncLen決定
  ② スクリプト実行（filterAndTrim含む）
  ③ afterプロット確認 → フィルタリング効果を評価
  ④ read_tracking.tsv確認 → 各ステップのread数推移を確認
     ↓ もし不満なら（merged readsが少ない等）
2回目の解析:
  ① truncLenを変更してs1スクリプト再実行
  ② afterプロット再確認
     ↓ OKなら
  次のステップへ進む


```

####  確認3: ログファイルの確認

```bash
# "解析完了"のメッセージを確認
grep "解析完了" ~/log/dada2_medium_123456.out.log

# エラーがないか確認
cat ~/log/dada2_medium_123456.err.log
```


### 6-2. 問題がある場合の対処

**問題: マージ率が低い（merged reads が少ない）**

→ `truncLen` を調整して **メイン解析をやり直す**

詳細は「[7. パラメータ調整ガイド](#7-パラメータ調整ガイド)」を参照

### 6-3. 後処理の投入

**すべての確認が完了し、問題がなければ:**

```bash
cd ~/16S_SILVA/scripts
sbatch js2_asv_postprocess.sh
```

**状態確認:**
```bash
squeue -u username
tail -f ~/log/dada2_post_*.out.log
```

**完了確認:**
```bash
# "ASV後処理完了"のメッセージを確認
grep "後処理完了" ~/log/dada2_post_*.out.log
```

---

## 7. パラメータ調整ガイド

### 7-1. truncLen（末端切り落とし長）の決め方

```r
truncLen = c(280, 220)  # R1=280bp, R2=220bp
```

**決定方法:**

1. **品質プロットを確認**
   - 緑の帯（Q30ライン）を下回る手前で切る

2. **オーバーラップを確保**
   - R1とR2の末端が重なる領域（overlap）が **50bpくらいあれば十分** 
   - 重なりすぎも重ならなすぎも良くない(低くても20くらい)
   - 計算式: `overlap = truncLen_R1 + truncLen_R2 - amplicon_length`
   - 例: `280 + 220 - 460 = 40bp`  OK
   - `R1＋R2=510bp`が目安。(460＋50)　　 
   - 例えば、R1:280のR2:230とか
   - R1の方が長くなりがち


**調整例:**

| パターン | R1 | R2 | Overlap | 判定 |
|---------|----|----|---------|------|
| 推奨 | 280 | 220 | 40bp |  Good |
| 短すぎ | 240 | 180 | -40bp |  マージ失敗 |
| 長すぎ | 300 | 260 | 100bp |  低品質領域を含む可能性 |

### 7-2. trimLeft（プライマー除去）の決め方

```r
trimLeft = c(17, 21)  # 341F=17bp, 805R=21bp
```

**決定方法:**

使用したプライマー配列の **ハイフン（-）以降の塩基数** を数える。

**生物技研のV3-V4プライマー（作業報告書より）:**

```
1st-341f: ACACTCTTTCCCTACACGACGCTCTTCCGATCTNNNNN-CCTACGGGNGGCWGCAG
                                                  └─────────┬─────────┘
                                                         17 bp

1st-805r: GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTNNNNN-GACTACHVGGGTATCTAATCC
                                                  └──────────┬──────────┘
                                                          21 bp
```

**解説:**
- ハイフン `-` より **前**: アダプター配列（DADA2実行前に除去済み）
- ハイフン `-` より **後**: 実際のプライマー配列（**DADA2で除去する必要がある**）

**計算:**
- 341F: `CCTACGGGNGGCWGCAG` → 17塩基
- 805R: `GACTACHVGGGTATCTAATCC` → 21塩基

→ `trimLeft = c(17, 21)`

**注意:** プライマーが異なる場合は、必ずこの値を変更すること。

### 7-3. その他のパラメータ

```r
maxN = 0           # N（不明塩基）を含むreadを除外
maxEE = c(2, 2)    # 期待エラー数の上限（R1, R2）
truncQ = 2         # 品質スコア2以下の位置で切断
```

**通常は変更不要です。**

---

## 8. 出力ファイル一覧

### 8-1. メイン解析の出力（js1実行後）

| ファイル名 | 説明 |
|-----------|------|
| `ASV_seqtab.rds` | キメラ除去前ASVテーブル |
| `ASV_nochim.rds` | キメラ除去後ASVテーブル |
| `taxonomy.rds` | ASV × Taxonomyテーブル |
| `filter_stats.tsv` | フィルタリング統計 |
| `read_tracking.tsv` | 全工程のread数推移 |
| `R1_quality_before.pdf` | R1品質プロット・フィルタリング前（★要確認） |
| `R2_quality_before.pdf` | R2品質プロット・フィルタリング前（★要確認） |
| `R1_quality_after.pdf` | R1品質プロット・フィルタリング後（★要確認） |
| `R2_quality_after.pdf` | R2品質プロット・フィルタリング後（★要確認） |
| `filtered/` | フィルタリング済みFASTQ（中間ファイル） |

### 8-2. 後処理の出力（js2実行後）

| ファイル名 | 説明 |
|-----------|------|
| `ASV_sequences.fasta` | 全ASVのDNA配列（FASTA形式） |
| `ASV_counts.tsv` | ASV × サンプルのカウントテーブル |
| `ASV_taxonomy.tsv` | ASV × 分類情報テーブル |
| `ASV_results.tsv` | カウント + タクソノミー統合テーブル |
| `ASV_results_relabund.tsv` | 相対存在量（%）テーブル |
| `ASV_Kingdom.tsv` | Kingdomレベル集計 |
| `ASV_Phylum.tsv` | Phylumレベル集計 |
| `ASV_Class.tsv` | Classレベル集計 |
| `ASV_Order.tsv` | Orderレベル集計 |
| `ASV_Family.tsv` | Familyレベル集計 |
| `ASV_Genus.tsv` | Genusレベル集計 |
| `dada2_post_workspace.RData` | Rワークスペース（全データ保存） |

---

## 9. SILVAデータベースのダウンロード_先生用

**ダウンロード元:**  
https://www.arb-silva.de/current-release/DADA2/1.36.0/SSU

**必要なファイル:**
```
silva_nr99_v138.2_toGenus_trainset.fa.gz
```
---

## 改訂履歴

- 2025/12/15: 初版作成（小林美穂）

---













## 自分用のメモ(小林)

# ペアエンドシーケンシングとオーバーラップの図解説明

## 1. ペアエンドシーケンシングの基本

### 1-1. イルミナシーケンサーがやっていること

```
【実際のDNA断片（V3-V4領域）】
5'==================================================3'  約460bp
   ↑                                              ↑
  341F                                          805R
  プライマー                                    プライマー


【イルミナの読み方】
       R1（Forward）→→→→→→→→→
5'==================================================3'
3'==================================================5'
       ←←←←←←←←←R2（Reverse）

・R1: 5'端から順方向に読む（280～300bp読める）
・R2: 3'端から逆方向に読む（220～280bp読める）
```

### 1-2. なぜ両方から読むのか？

**理由1:** イルミナは長い配列を一度に読めない（1回で300bp前後が限界）  
**理由2:** V3-V4領域は約460bpなので、両端から読んで「真ん中で合わせる」必要がある

---

## 2. オーバーラップ（重なり）とは

### 2-1. 視覚的な説明

```
【オーバーラップがある場合（正常）】

R1: 280bp読んだ
5'━━━━━━━━━━━━━━━━━━━━━━━━━━━→
                          ■■■■■■■  ← この部分が重なり（overlap）
                    ←━━━━━━━━━━━━━━━━━━━━━━3'
                    R2: 220bp読んだ

重なり領域 = 280 + 220 - 460 = 40bp


【オーバーラップが無い場合（マージ失敗）】

R1: 240bp読んだ
5'━━━━━━━━━━━━━━━━→
                      ?????? ← 読めていない領域（ギャップ）
                          ←━━━━━━━━━━━━━3'
                          R2: 180bp読んだ

重なり = 240 + 180 - 460 = -40bp （マイナス！）
→ つながらない！
```

### 2-2. オーバーラップの計算式

```
オーバーラップ（bp） = R1の長さ + R2の長さ - DNA断片の実際の長さ

例1: truncLen = c(280, 220), V3-V4 = 460bp
     → 280 + 220 - 460 = 40bp  OK

例2: truncLen = c(240, 180), V3-V4 = 460bp
     → 240 + 180 - 460 = -40bp  NG（ギャップができる）

例3: truncLen = c(300, 260), V3-V4 = 460bp
     → 300 + 260 - 460 = 100bp  重なりすぎ（低品質領域を含む可能性）
```
---

---

## 3. マージ（merge）とは

### 3-1. マージの意味

**マージ = R1とR2を1本の配列につなぐこと**

```
【マージ前】
R1: ATCGATCGATCG...（Forward read）
R2: CGTAGCTAGCTA...（Reverse read）

【マージ後】
1本の完全な配列: ATCGATCGATCG...CGTAGCTAGCTA
```

### 3-2. DADA2のマージプロセス

```
1. R1とR2の重なり部分を探す
2. 重なり部分の配列が一致するか確認
   ┌─────┐
   │ATCGC│ ← R1の末尾
   │ATCGC│ ← R2の先頭（相補鎖なので逆向き変換後）
   └─────┘
   一致！ → マージ成功

3. 一致しない場合 → マージ失敗（そのreadは捨てられる）
```

---

End of document.

```
