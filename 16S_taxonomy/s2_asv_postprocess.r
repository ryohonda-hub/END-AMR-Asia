#!/usr/bin/env Rscript

library(dada2)
library(dplyr)
library(readr)

# ------------------------------------------------------------------
# ディレクトリとファイル
# ------------------------------------------------------------------
OUT_DIR <- "/home/username/16S_SILVA/output"
ASV_RDS <- file.path(OUT_DIR, "ASV_nochim.rds")
TAX_RDS <- file.path(OUT_DIR, "taxonomy.rds")

seqtab.nochim <- readRDS(ASV_RDS)
taxa <- readRDS(TAX_RDS)

# ------------------------------------------------------------------
# 13. ASV fasta 出力
# ------------------------------------------------------------------
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- paste0(">ASV_", seq_len(length(asv_seqs)))
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, file.path(OUT_DIR, "ASV_sequences.fasta"))

# ------------------------------------------------------------------
# 14. ASV カウントテーブル出力
# ------------------------------------------------------------------
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, file.path(OUT_DIR, "ASV_counts.tsv"), sep="\t", quote=FALSE, col.names=NA)

# ------------------------------------------------------------------
# 15. タクソノミーテーブル出力
# ------------------------------------------------------------------
asv_tax <- taxa
row.names(asv_tax) <- sub(">", "", asv_headers)
write.table(asv_tax, file.path(OUT_DIR, "ASV_taxonomy.tsv"), sep="\t", quote=FALSE, col.names=NA)

# ------------------------------------------------------------------
# 16. 結合テーブル作成
# ------------------------------------------------------------------
counts <- read.delim(file.path(OUT_DIR, "ASV_counts.tsv"), sep="\t", header=TRUE, check.names=FALSE)
tax    <- read.delim(file.path(OUT_DIR, "ASV_taxonomy.tsv"), sep="\t", header=TRUE, check.names=FALSE)
colnames(counts)[1] <- "ASV"
colnames(tax)[1] <- "ASV"

ASV_results <- tax %>% inner_join(counts, by="ASV")
write.table(ASV_results, file.path(OUT_DIR, "ASV_results.tsv"), sep="\t", row.names=FALSE, quote=FALSE)

# ------------------------------------------------------------------
# 17. 相対存在量計算 + 階層別集計
# ------------------------------------------------------------------
df <- read_tsv(file.path(OUT_DIR, "ASV_results.tsv"), show_col_types = FALSE)
tax_cols <- c("Kingdom","Phylum","Class","Order","Family","Genus","ASV")
sample_cols <- setdiff(colnames(df), tax_cols)

# 数値化
df[sample_cols] <- lapply(df[sample_cols], function(x){
    as_num <- suppressWarnings(as.numeric(x))
    as_num[is.na(as_num)] <- 0
    as_num
})

# 相対存在割合の計算
sample_sums <- colSums(df[sample_cols], na.rm=TRUE)
df_relabund <- df
df_relabund[sample_cols] <- sweep(df[sample_cols], 2, sample_sums, FUN="/")
write_tsv(df_relabund, file.path(OUT_DIR, "ASV_results_relabund.tsv"))

# 階層別集計関数
safe_sum_tax <- function(df, group_cols, sample_cols){
    missing_cols <- setdiff(group_cols, colnames(df))
    if(length(missing_cols)>0){ return(tibble()) }
    df %>% group_by(across(all_of(group_cols))) %>%
        summarise(across(all_of(sample_cols), ~sum(.x, na.rm=TRUE)), .groups="drop")
}

ASV_Kingdom <- safe_sum_tax(df_relabund, c("Kingdom"), sample_cols)
ASV_Phylum  <- safe_sum_tax(df_relabund, c("Kingdom","Phylum"), sample_cols)
ASV_Class   <- safe_sum_tax(df_relabund, c("Kingdom","Phylum","Class"), sample_cols)
ASV_Order   <- safe_sum_tax(df_relabund, c("Kingdom","Phylum","Class","Order"), sample_cols)
ASV_Family  <- safe_sum_tax(df_relabund, c("Kingdom","Phylum","Class","Order","Family"), sample_cols)
ASV_Genus   <- safe_sum_tax(df_relabund, c("Kingdom","Phylum","Class","Order","Family","Genus"), sample_cols)

write_tsv(ASV_Kingdom,  file.path(OUT_DIR,"ASV_Kingdom.tsv"))
write_tsv(ASV_Phylum,  file.path(OUT_DIR,"ASV_Phylum.tsv"))
write_tsv(ASV_Class,   file.path(OUT_DIR,"ASV_Class.tsv"))
write_tsv(ASV_Order,   file.path(OUT_DIR,"ASV_Order.tsv"))
write_tsv(ASV_Family,  file.path(OUT_DIR,"ASV_Family.tsv"))
write_tsv(ASV_Genus,   file.path(OUT_DIR,"ASV_Genus.tsv"))

# ワークスペース保存
save.image(file.path(OUT_DIR,"dada2_post_workspace.RData"))

cat("\n============ ASV後処理完了 ============\n")
cat("出力ファイル:\n")
cat("ASV_sequences.fasta, ASV_counts.tsv, ASV_taxonomy.tsv, ASV_results.tsv, ASV_results_relabund.tsv\n")  
cat("階層別集計: ASV_Kingdom~ASV_Genus.tsv\n")
