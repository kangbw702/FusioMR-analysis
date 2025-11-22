library(data.table); library(tidyverse)

# merge
# 1) Find all matching files in the working directory
files = Sys.glob("sumtbl_bryois/sumtbl_ad_0.1_0.001_50_0.1_0.9_0.6_50_1_*_*.txt")
cat(length(files), "files\n")

# 2) Read + annotate each file
dt_list = lapply(files, function(f) {
dt = fread(f)
# Parse cell_type and chr from the file name
parts = strsplit(basename(f), "_")[[1]]
# last 2 tokens are: [cell_type] and [chr].txt
cell_type = parts[length(parts) - 1]
chr = sub("\\.txt$", "", parts[length(parts)])
# Add columns
dt[, `:=`(cell_type = cell_type, chr = as.integer(chr))]
})

# 3) Merge (rbind) them all; fill=TRUE handles column mismatches gracefully
merged = rbindlist(dt_list, use.names = TRUE, fill = TRUE)

# 4) Check chromosomes are 1–22
if (any(!merged$chr %in% 1:22)) warning("Some rows have chr outside 1–22; check filenames.")

# 5) Save the merged file
fwrite(merged, "merged/sumtbl_ad_merged.csv")


# Summary
# cd /gpfs/data/linchen-lab/Bowei/fusiomr_sc_bk/r1/
library(data.table); library(tidyverse)

# function to calculate genomic inflation factor
calculate_lambda <- function(p_values) {
p_values = p_values[!is.na(p_values)]
chi_sq_values <- qchisq(1 - p_values, df=1)
lambda <- median(chi_sq_values) / qchisq(0.5, df=1)
return(lambda)
}

niv_cut1 = 2; niv_cut2 = 5; niv_cut3 = 8
disease_short = 'ad'
ctypes_short = c("end","ast","exc","inh","mic","opc","oli","per")
tmp1 = fread("merged/sumtbl_ad_merged.csv") %>% dplyr::filter(niv >= niv_cut3) %>% drop_na(gene)

tbl = NULL
for (i in 1:length(ctypes_short)) {
cat(i, ctypes_short[i], "\n")
df = tmp1 %>% dplyr::filter(cell_type == ctypes_short[i]) 
K = nrow(df)
lambdas = sapply(df %>% dplyr::select(contains('p_')), calculate_lambda)
tmp = c(ctypes_short[i], K, round(lambdas, 3), sum(df$p_ivw<(0.05/K), na.rm = T), sum(df$p_egger<(0.05/K), na.rm = T), sum(df$p_cml<(0.05/K), na.rm = T), sum(df$p_cml_dp<(0.05/K), na.rm = T), sum(df$p_fusio<(0.05/K), na.rm = T))
tbl = rbind(tbl, tmp)
}

tmp2 = tmp1 %>% dplyr::filter(niv < niv_cut2) 
tbl_smallk = NULL
for (i in 1:length(ctypes_short)) {
cat(i, ctypes_short[i], "\n")
df = tmp2 %>% dplyr::filter(cell_type == ctypes_short[i]) 
K = nrow(df)
lambdas = sapply(df %>% dplyr::select(contains('p_')), calculate_lambda)
tmp = c(ctypes_short[i], K, round(lambdas, 3), sum(df$p_ivw<(0.05/K), na.rm = T), sum(df$p_egger<(0.05/K), na.rm = T), sum(df$p_cml<(0.05/K), na.rm = T), sum(df$p_cml_dp<(0.05/K), na.rm = T), sum(df$p_fusio<(0.05/K), na.rm = T))
tbl_smallk = rbind(tbl_smallk, tmp)
}

colnames(tbl) = colnames(tbl_smallk) = c('ctype', 'ntest', 'lambda_ivw', 'lambda_egger', 'lambda_cml', 'lambda_cml_dp', 'lambda_fusio', 'n_sig_ivw', 'n_sig_egger', 'n_sig_cml', 'n_sig_cml_dp', 'n_sig_fusio')
tbl; tbl_smallk


# Fig2
## get sig genes list. p.adj.Bonf < 0.05 in at least one cell type

library(data.table); library(tidyverse); library(pheatmap); library(viridis); library(EnsDb.Hsapiens.v79); library(grid)
df = fread("sumtbl_ad_merged.csv") %>% dplyr::select(gene,niv,p_fusio,cell_type,chr) %>% dplyr::filter(niv>=5) %>% group_by(cell_type) %>% 
mutate(nn=n()) %>% ungroup() %>% mutate(sig = p_fusio<(0.05/nn), nlogp = -log(p_fusio,10)) 
genes = unique(df$gene[df$sig]) 
df = df %>% dplyr::filter(gene %in% genes)
df_wide = df %>% dplyr::select(gene, nlogp, cell_type) %>% pivot_wider(names_from = cell_type, values_from = nlogp) 
# rank genes
summary(df$nlogp) 
df_wide = df_wide %>% rowwise() %>% 
mutate(n5 = sum(c_across(end:per) > 5, na.rm=T), n4 = sum(c_across(end:per) > 4, na.rm=T), n3 = sum(c_across(end:per) > 4, na.rm=T)) %>% 
ungroup() %>% arrange(desc(n5), desc(n4), desc(n3)) %>% column_to_rownames(var = "gene") 
# rank ct
mat = t(as.matrix(df_wide %>% dplyr::select(ast:per)))
orders = order(apply(mat, 1, function(t) sum(t>4, na.rm=T)), decreasing = T)
ctypes_short = c("ast","end","exc","inh","mic","opc","oli","per")
mat = t(as.matrix(df_wide %>% dplyr::select(ctypes_short[orders])))
mat[mat>10] = 10

## build significance matrix ("*" for sig, "" otherwise) 
# wide matrix of sig (logical) with same structure as df_wide
sig_wide = df %>% dplyr::select(gene, cell_type, sig) %>% pivot_wider(names_from = cell_type, values_from = sig)
# align rows (genes) to df_wide row order
sig_wide = sig_wide %>% dplyr::filter(gene %in% rownames(df_wide)) %>% arrange(match(gene, rownames(df_wide))) %>% column_to_rownames("gene")

## heatmap 
tissues = c("Amygdala","BA24","Caudate_basal_ganglia","Cerebellar_Hemi","Cerebellum","Cortex","Hippocampus","Hypothalamus","BA9","Nucleus_accumbens_basal_ganglia","Putamen_basal_ganglia","Spinal_cord","Substantia_nigra")
df_all = as.data.frame(t(mat)) %>% rownames_to_column("gene") 
df_all_sig = sig_wide %>% rownames_to_column("gene")
for (i in 1:length(tissues)) {
print(tissues[i])
df_bk = fread(paste0("~/Dropbox/Lin/FusioMR/af_predict/scmr/merged/bk_s2_ad_0.001_50_0.1_", tissues[i], ".txt")) %>% dplyr::select(gene, niv, p_fusio) %>% dplyr::filter(niv>=5) %>% drop_na() 
nn = nrow(df_bk)
df_bk = df_bk %>% mutate(nlogp = -log(p_fusio,10), sig = p_fusio < (0.05/nn))
# Convert from ensembl.gene to gene.symbol
ensembl.genes = df_bk$gene
geneIDs1 = ensembldb::select(EnsDb.Hsapiens.v79, keys= ensembl.genes, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
df_bk = df_bk %>% left_join(geneIDs1, by = c("gene"="GENEID"))
df_all = df_all %>% left_join(df_bk[,c("SYMBOL","nlogp")], by=c("gene"="SYMBOL")) %>% dplyr::rename(!!tissues[i] := nlogp)
## join significance
df_all_sig = df_all_sig %>% left_join(df_bk[, c("SYMBOL", "sig")], by = c("gene" = "SYMBOL")) %>% dplyr::rename(!!tissues[i] := sig)
}

# build combined matrix
df_all = df_all %>% column_to_rownames("gene")
df_all_sig = df_all_sig %>% column_to_rownames("gene")
## matrices: rows = cell types or tissues, cols = genes
sc_mat = t(as.matrix(df_all[, ctypes_short]))
bk_mat = t(as.matrix(df_all[, tissues]))
## same for significance
sc_sig_mat = t(as.matrix(df_all_sig[, ctypes_short]))
bk_sig_mat = t(as.matrix(df_all_sig[, tissues]))
## combined
mat = rbind(sc_mat, bk_mat)
sig_mat = rbind(sc_sig_mat, bk_sig_mat)
## star matrix: * for TRUE, "" for FALSE/NA
star_mat = ifelse(is.na(sig_mat), "", ifelse(sig_mat, "*", ""))

#Row annotation: data source (“sc” vs “bk”)
annotation_row = data.frame(data_source = c(rep("sc", nrow(sc_mat)), rep("bulk", nrow(bk_mat))))
rownames(annotation_row) = rownames(mat)
ann_colors = list(data_source = c(sc = "#f781bf", bulk = "#80b1d3"))

ph = pheatmap(
mat,
cluster_rows   = FALSE,
cluster_cols   = FALSE,
color = colorRampPalette(c("white","lightblue","royalblue","blue"))(100),
main = "",
display_numbers = star_mat,
number_color   = "black",
fontsize_number = 12,
annotation_row = annotation_row,
annotation_colors = ann_colors,
gaps_row = nrow(sc_mat),         # visual gap between sc and bulk
border_color = NA,
legend_title = expression(-log[10](p))
)

pdf("fig2_heatmap.pdf", width = 10, height = 6) 
grid.newpage()
grid.draw(ph$gtable)
grid.text(
expression(-log[10](p)),
x = 0.8, y = 0.95,
gp = gpar(fontsize = 10)
)

dev.off()


