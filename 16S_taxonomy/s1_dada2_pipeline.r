############################################################
# 16S rRNA V3–V4（341f/805r）ASV 解析テンプレート（HPC版）
# R 4.5+ / DADA2 / SILVA v138.2
############################################################

library(dada2)
library(ggplot2)
library(dplyr)
library(readr)

# ------------------------------------------------------------------
# 1. ディレクトリ設定（★必ず変更）
# ------------------------------------------------------------------
## プロジェクトの基準ディレクトリ（★必ず自分の場所に変更）→ 全ての解析結果(output)がここにまとまる
PROJECT_DIR <- "/home/username/16S_SILVA"  
## FASTQ が置いてあるディレクトリ(必ず一か所にまとめる)
FASTQ_DIR   <- file.path(PROJECT_DIR, "raw_fastq")
# → 出力ディレクトリ
OUT_DIR     <- file.path(PROJECT_DIR, "output")  # ★自動作成

## Silva データベースのパス
silva_file <- "/home/ryohonda/db/dada2/SILVA_138.2/silva_nr99_v138.2_toGenus_trainset.fa.gz"   # ★ここを先生のdbのパスに変更

if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)
setwd(OUT_DIR)

# ------------------------------------------------------------------
# 2. FASTQ 読み込み
# ------------------------------------------------------------------
fnFs <- sort(list.files(FASTQ_DIR, pattern="_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(FASTQ_DIR, pattern="_R2_001.fastq.gz", full.names = TRUE))

sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
cat("検出サンプル数:", length(sample.names), "\n")  # ★サンプル数自動検出

# ------------------------------------------------------------------
# 3. フィルタリングフォルダ作成
# ------------------------------------------------------------------
filt_path <- file.path(OUT_DIR, "filtered")
if (!dir.exists(filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

# ------------------------------------------------------------------
# 4. 品質プロット（フィルタリング前）
# ------------------------------------------------------------------
n_samples <- length(fnFs)
plot_height <- max(10, n_samples * 0.8)  # サンプル数に応じて高さ調整（最低10）

p1 <- plotQualityProfile(fnFs[1:n_samples])
p2 <- plotQualityProfile(fnRs[1:n_samples])
ggsave("R1_quality_before.pdf", p1, width = 20, height = plot_height, dpi = 300)
ggsave("R2_quality_before.pdf", p2, width = 20, height = plot_height, dpi = 300)

# ------------------------------------------------------------------
# 5. フィルタリング
# ------------------------------------------------------------------
out <- filterAndTrim(
    fwd = fnFs, filt = filtFs,
    rev = fnRs, filt.rev = filtRs,
    trimLeft = c(17, 21),      # ★プライマー長(生物技研の作業報告書より)
    truncLen = c(280, 220),    # ★品質プロットで調整
    maxN = 0, maxEE = c(2, 2), 
    truncQ = 2, compress = TRUE, multithread = TRUE
)
write.table(out, "filter_stats.tsv", sep="\t", quote=FALSE, col.names=NA)

# ------------------------------------------------------------------
# 5.5. フィルタリング後の品質プロット
# ------------------------------------------------------------------
p3 <- plotQualityProfile(filtFs[1:n_samples])
p4 <- plotQualityProfile(filtRs[1:n_samples])
ggsave("R1_quality_after.pdf", p3, width = 20, height = plot_height, dpi = 300)
ggsave("R2_quality_after.pdf", p4, width = 20, height = plot_height, dpi = 300)


# ------------------------------------------------------------------
# 6. エラー率学習
# ------------------------------------------------------------------
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

# ------------------------------------------------------------------
# 7. ASV 推定
# ------------------------------------------------------------------
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

# ------------------------------------------------------------------
# 8. マージ
# ------------------------------------------------------------------
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

# ------------------------------------------------------------------
# 9. ASV テーブル
# ------------------------------------------------------------------
seqtab <- makeSequenceTable(mergers)
saveRDS(seqtab, "ASV_seqtab.rds")

# ------------------------------------------------------------------
# 10. キメラ除去
# ------------------------------------------------------------------
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE)
saveRDS(seqtab.nochim, "ASV_nochim.rds")

# ------------------------------------------------------------------
# 11. リードトラッキング
# ------------------------------------------------------------------
getN <- function(x) sum(getUniques(x))

if (length(sample.names) == 1) {
    track <- data.frame(
        input = out[1],
        filtered = out[2],
        denoisedF = getN(dadaFs[[1]]),
        denoisedR = getN(dadaRs[[1]]),
        merged = getN(mergers[[1]]),
        nonchim = sum(seqtab.nochim),
        final_perc_reads_retained = round((sum(seqtab.nochim)/out[1])*100, 1)  
    )
    rownames(track) <- sample.names
} else {
    track <- cbind(
        out,
        sapply(dadaFs, getN),
        sapply(dadaRs, getN),
        sapply(mergers, getN),
        rowSums(seqtab.nochim),
        round((rowSums(seqtab.nochim)/out[,1])*100, 1)  
    )
    colnames(track) <- c("input","filtered","denoisedF","denoisedR","merged","nonchim","final_perc_reads_retained")  
    rownames(track) <- sample.names
}
write.table(track, "read_tracking.tsv", sep="\t", quote=FALSE, col.names=NA)

# ------------------------------------------------------------------
# 12. SILVA タクソノミー付与
# ------------------------------------------------------------

taxa <- assignTaxonomy(seqtab.nochim, silva_file, multithread=TRUE)
saveRDS(taxa, "taxonomy.rds")

# ------------------------------------------------------------------
# 13〜19. ASV fasta・テーブル・相対存在量・集計
# （後処理スクリプト s2_asv_postprocess.R で実行）
# ------------------------------------------------------------------
cat("\n============ 解析完了 ============\n")
