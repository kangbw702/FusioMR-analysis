# summarize
library(data.table); library(tidyverse)
# function to calculate the genomic inflation factor
calculate_lambda <- function(p_values) {
p_values = p_values[!is.na(p_values)]
chi_sq_values = qchisq(1 - p_values, df=1)
lambda = median(chi_sq_values) / qchisq(0.5, df=1)
return(lambda)
}
tissue_short = 
c("Amygdala","BA24","Caudate_basal_ganglia","Cerebellar_Hemi","Cerebellum","Cortex","BA9","Hippocampus","Hypothalamus","Nucleus_accumbens_basal_ganglia","Putamen_basal_ganglia","Spinal_cord","Substantia_nigra","Blood","Aorta","Coronary","Left_Ventricle","Atrial_Appendage")
res = NULL
for (i in 1:length(tissue_short)) {
#for (i in 5:5) {
aa = fread(paste0("res_r1/fusio_res_", tissue_short[i], "_s7.txt"), fill=T) %>% distinct(V3, .keep_all = TRUE)
colnames(aa) = 
c("tissue","chr","gene","K", "b_fusiom1", "se_fusiom1", "p_fusiom1", "b_fusiom2", "se_fusiom2", "p_fusiom2", "b_fusios", "se_fusios", "p_fusios")
K = nrow(aa)
res = rbind(res, c(i, tissue_short[i], K, 
round(sapply(aa %>% dplyr::select(contains('p_')), calculate_lambda), 3), 
sapply(aa %>% dplyr::select(contains('p_')), function(t) sum((t*K) < 0.05, na.rm = T))
))}
res

# significant genes in at least one tissue
genes_fusiom1 = genes_fusiom2 = genes_fusios = NULL
for (i in 1:length(tissue_short)) {
aa = fread(paste0("res_r1/fusio_res_", tissue_short[i], "_s7.txt"), fill=T) %>% distinct(V3, .keep_all = TRUE)
colnames(aa) = 
c("tissue","chr","gene","K", "b_fusiom1", "se_fusiom1", "p_fusiom1", "b_fusiom2", "se_fusiom2", "p_fusiom2", "b_fusios", "se_fusios", "p_fusios")
K = nrow(aa)
genes_fusiom1 = union(genes_fusiom1, aa$gene[aa$p_fusiom1 < (0.05/K)])
genes_fusiom2 = union(genes_fusiom2, aa$gene[aa$p_fusiom2 < (0.05/K)])
genes_fusios  = union(genes_fusios,  aa$gene[aa$p_fusios  < (0.05/K)])
}


# sig genes
library(data.table); library(tidyverse)
tissue_short = c("Amygdala","BA24","Caudate_basal_ganglia","Cerebellar_Hemi","Cerebellum","Cortex","BA9","Hippocampus","Hypothalamus","Nucleus_accumbens_basal_ganglia","Putamen_basal_ganglia","Spinal_cord","Substantia_nigra","Blood","Aorta","Coronary","Left_Ventricle","Atrial_Appendage")
#traits = c("stroke", "af")
tbls = tbl1m = tbl2m =NULL
# get list of the union of sig. genes (either traits)
sig_gene_lists = sig_gene_list1m = sig_gene_list2m = NULL
for (j in 1:length(tissue_short)) {
tissue = tissue_short[j]
# tbl = fread(paste0("merged/apa_", traits[1], "_", traits[2], "_", tissue_short[j], "_", pcut, "_", clump_kb, "_", clump_r2, ".txt")) %>% drop_na(gene) 
tbl = fread(paste0("apa_result/fusio_res_", tissue, ".txt"), fill=T) %>% drop_na(V3) %>% distinct(V3, .keep_all = TRUE)
colnames(tbl) = c("tissue","chr","gene","niv", "b_fusiom1", "se_fusiom1", "p_fusiom1", "b_fusiom2", "se_fusiom2", "p_fusiom2", "b_fusios", "se_fusios", "p_fusios")
K = sum(tbl$niv >=5)
tbls = tbl %>% dplyr::filter(p_fusios<(0.05/K), niv>=5) %>% dplyr::select(gene) 
tbl1m = tbl %>% dplyr::filter(p_fusiom1<(0.05/K), niv>=5) %>% dplyr::select(gene) 
tbl2m = tbl %>% dplyr::filter(p_fusiom2<(0.05/K), niv>=5) %>% dplyr::select(gene) 
sig_gene_lists = c(sig_gene_lists, tbls$gene)
sig_gene_list1m = c(sig_gene_list1m, tbl1m$gene)
sig_gene_list2m = c(sig_gene_list2m, tbl2m$gene)
}

nshare_cut = 1
sig_gene_lists = names(which(table(sig_gene_lists)>=nshare_cut))
sig_gene_list1m = names(which(table(sig_gene_list1m)>=nshare_cut))
sig_gene_list2m = names(which(table(sig_gene_list2m)>=nshare_cut))
# unique
length(sig_gene_lists)  
length(sig_gene_list1m)
length(sig_gene_list2m)
length(union(sig_gene_list1m, sig_gene_list2m)) 

# get pval of genes in the list (union)
tbl = NULL
for (j in 1:length(tissue_short)) {
print(paste0(j, ", ", tissue_short[j]))
file_path1 = paste0("~/Dropbox/Lin/FusioMR/af_predict/r1_analysis1/apa_result/fusio_res_", tissue_short[j], "_s7.txt")
tmp = fread(file_path1, fill=T) %>% drop_na(V3) %>% distinct(V3, .keep_all = TRUE)
colnames(tmp) = c("tissue","chr","gene","niv", "b_fusiom1", "se_fusiom1", "p_fusiom1", "b_fusiom2", "se_fusiom2", "p_fusiom2", "b_fusios", "se_fusios", "p_fusios")
K = sum(tmp$niv >=5)
tmp = tmp %>% dplyr::filter(gene %in% union(sig_gene_list1m, sig_gene_list2m), niv >= 5) %>% mutate(sig1 = ifelse(p_fusiom1<(0.05/K) & niv>=5, 1, 0), sig2 = ifelse(p_fusiom2<(0.05/K) & niv>=5, 1, 0), tissue = tissue_short[j], nlog10p1 = -log10(p_fusiom1), nlog10p2 = -log10(p_fusiom2)) %>% dplyr::rename(b1m=b_fusiom1, b2m=b_fusiom2, b2s=b_fusios, p1m=p_fusiom1, p2m=p_fusiom2, p2s=p_fusios) %>% dplyr::select(tissue, gene, b1m, b2m, b2s, p1m, p2m, p2s, nlog10p1, nlog10p2, sig1, sig2) %>% drop_na() 
tbl = rbind(tbl, tmp)
} # end loop j
fwrite(tbl, 'apa_sig_genes.txt')

            
# is a sig gene tissue shared or specific?
library(data.table); library(tidyverse)
tissue_short = c("Amygdala","BA24","Caudate_basal_ganglia","Cerebellar_Hemi","Cerebellum","Cortex","BA9","Hippocampus","Hypothalamus","Nucleus_accumbens_basal_ganglia","Putamen_basal_ganglia","Spinal_cord","Substantia_nigra","Blood","Aorta","Coronary","Left_Ventricle","Atrial_Appendage")
tbl2m = sig_gene_list2m = NULL
# generate list of sig genes per tissue
sig_genes = lapply(tissue_short, function(t) {
tbl = fread(paste0("apa_result/fusio_res_", t, ".txt"), fill=T) %>% drop_na(V3) %>% distinct(V3, .keep_all = TRUE)
colnames(tbl) = c("tissue","chr","gene","niv", "b_fusiom1", "se_fusiom1", "p_fusiom1", "b_fusiom2", "se_fusiom2", "p_fusiom2", "b_fusios", "se_fusios", "p_fusios")
K2 = sum(tbl$niv >=5)
tbl2m = tbl %>% dplyr::filter(p_fusiom2<(0.05/K2), niv>=5) %>% dplyr::select(gene) 
return(tbl2m$gene)
})
names(sig_genes) = tissue_short
lengths(sig_genes)
# All unique genes across tissues
all_genes = unique(unlist(sig_genes))
length(all_genes)
# Build presence/absence matrix
gene_mat = sapply(sig_genes, function(g) all_genes %in% g)
rownames(gene_mat) = all_genes
# Count number of tissues where each gene is significant
tissue_counts = rowSums(gene_mat)
# Identify tissue-specific and shared genes
shared_genes = all_genes[tissue_counts > 1]
length(shared_genes)
length(shared_genes)/length(all_genes)
tissue_specific_genes = all_genes[tissue_counts == 1]
length(tissue_specific_genes)
length(tissue_specific_genes)/length(all_genes)
# Count per tissue
df = data.frame(tissue = names(sig_genes), shared = NA_integer_, specific = NA_integer_)
for (t in names(sig_genes)) {
g_set = sig_genes[[t]]
df[df$tissue == t, "shared"] = sum(g_set %in% shared_genes)
df[df$tissue == t, "specific"] =sum(g_set %in% tissue_specific_genes)
}
df = df %>% mutate(spec_percent = specific/(shared + specific))
sum(df$specific)
df


# Figure S1: distribution of n.tissue and direction consistency
# for af only
library(tidyverse); library(data.table); library(ggplot2)
tissue_short = c("Amygdala","BA24","Caudate_basal_ganglia","Cerebellar_Hemi","Cerebellum","Cortex","BA9","Hippocampus","Hypothalamus","Nucleus_accumbens_basal_ganglia","Putamen_basal_ganglia","Spinal_cord","Substantia_nigra","Blood","Aorta","Coronary","Left_Ventricle","Atrial_Appendage")
tbls = tblm = NULL
# get list of sig. genes for AF, fusios and fusiom
sig_gene_lists = sig_gene_listm = NULL
for (j in 1:length(tissue_short)) {
print(tissue_short[j])
tbl = fread(paste0("apa_result/fusio_res_", tissue_short[j], ".txt"), fill=T) %>% drop_na(V3) %>% distinct(V3, .keep_all = TRUE)
colnames(tbl) = c("tissue","chr","gene","niv", "b_fusiom1", "se_fusiom1", "p_fusiom1", "b_fusiom2", "se_fusiom2", "p_fusiom2", "b_fusios", "se_fusios", "p_fusios")
K2 = sum(tbl$niv >=5)
tbls = tbl %>% dplyr::filter(p_fusios<(0.05/K2), niv>=5) 
tblm = tbl %>% dplyr::filter(p_fusiom2<(0.05/K2), niv>=5) 
sig_gene_lists = c(sig_gene_lists, tbls$gene)
sig_gene_listm = c(sig_gene_listm, tblm$gene)
}
nshare_cut = 1
sig_gene_lists = names(which(table(sig_gene_lists)>=nshare_cut)) # 94
sig_gene_listm = names(which(table(sig_gene_listm)>=nshare_cut)) # 108
length(intersect(sig_gene_lists, sig_gene_listm)) # 87

# get pval of genes in the list (union)
tbl = NULL
for (j in 1:length(tissue_short)) {
print(paste0(j, ", ", tissue_short[j]))
tmp = fread(paste0("apa_result/fusio_res_", tissue_short[j], "_s7.txt"), fill=T) %>% drop_na(V3) %>% distinct(V3, .keep_all = TRUE)
colnames(tmp) = c("tissue","chr","gene","niv", "b_fusiom1", "se_fusiom1", "p_fusiom1", "b_fusiom2", "se_fusiom2", "p_fusiom2", "b_fusios", "se_fusios", "p_fusios")
K2 = sum(tmp$niv >=5)
tmp = tmp %>% dplyr::filter(gene %in% union(sig_gene_lists, sig_gene_listm)) %>% mutate(sigs = ifelse(p_fusios<(0.05/K2) & niv>=5, 1, 0), sigm = ifelse(p_fusiom2<(0.05/K2) & niv>=5, 1, 0), tissue = tissue_short[j], nlog10ps = -log10(p_fusios), nlog10pm = -log10(p_fusiom2)) %>% dplyr::rename(bs=b_fusios, bm=b_fusiom2) %>% dplyr::select(tissue, gene, bs, bm, nlog10ps, nlog10pm, sigs, sigm) %>% drop_na() 
tbl = rbind(tbl, tmp)
} # end loop j
tbls = tbl %>% dplyr::filter(sigs==1 & gene %in% sig_gene_lists)
tblm = tbl %>% dplyr::filter(sigm==1 & gene %in% sig_gene_listm)
table(tbls$tissue)
table(tblm$tissue)

# across tissues distribution
tmp = tblm %>% group_by(gene) %>% summarise(nn=n()); table(tmp$nn)

# effect direction consistency 
# only genes in sig_gene_listm
df = tbl %>% dplyr::filter(gene %in% sig_gene_listm) %>% dplyr::filter(gene %in% sig_gene_listm) %>% group_by(gene) %>% mutate(mean_b = mean(bm)) %>% ungroup() %>% mutate(gene = reorder(gene, mean_b))
pp = ggplot(df, aes(x = gene, y = bm)) + geom_boxplot(outlier.size = 1, fill = "skyblue", color = "blue") + theme_bw() + theme(text = element_text(size = 14), axis.text.x = element_text(size = 8, angle = 45, hjust = 1)) + geom_hline(yintercept = 0) + labs(x = "Gene", y = "Effect size")

ggsave("effect_size_dotplot.pdf", plot = pp, width = 12, height = 7)


# compare with eGenes (GTEx)
library(data.table); library(tidyverse); library(EnsDb.Hsapiens.v79)
tissue_short = c("Amygdala","BA24","Caudate_basal_ganglia","Cerebellar_Hemi","Cerebellum","Cortex","BA9","Hippocampus","Hypothalamus","Nucleus_accumbens_basal_ganglia","Putamen_basal_ganglia","Spinal_cord","Substantia_nigra","Blood","Aorta","Coronary","Left_Ventricle","Atrial_Appendage")
apagenes_combined = egenes_combnied = NULL
for (j in 1:length(tissue_short)) {
tissue = tissue_short[j]
# print(tissue)
# get apa genes
df_apa = fread(paste0("apa_result/fusio_res_", tissue, "_s7.txt"), fill=T) %>% drop_na(V3) %>% distinct(V3, .keep_all = TRUE)
colnames(df_apa) = c("tissue","chr","gene","niv", "b_fusiom1", "se_fusiom1", "p_fusiom1", "b_fusiom2", "se_fusiom2", "p_fusiom2", "b_fusios", "se_fusios", "p_fusios")
df_apa = df_apa[df_apa$niv >= 5,]
K_apa = nrow(df_apa)
df_apa = df_apa %>% mutate(sig = p_fusiom2 < (0.05/K_apa))
apagenes = df_apa$gene[df_apa$sig]
apagenes_combined = c(apagenes_combined, apagenes)
}
length(unique(apagenes_combined))

# get egenes
for (j in 1:length(tissue_short)) {
tissue = tissue_short[j]
df_e = fread(paste0("merged/bkeqtl_stroke_af_", tissue, "_0.005_50_0.01.txt")) %>% dplyr::select(gene, niv, b_fusiom2:p_fusiom2) %>% drop_na(gene) 
# Convert from ensembl.gene to gene.symbol
ensembl.genes = df_e$gene
geneIDs1 = ensembldb::select(EnsDb.Hsapiens.v79, keys= ensembl.genes, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
df_e = df_e %>% left_join(geneIDs1, by = c("gene"="GENEID")) %>% drop_na() %>% dplyr::filter(niv >= 5) 
K_e = nrow(df_e)
# df_e = df_e %>% mutate(sig = p_fusiom2 < (0.05/K_e))
df_e = df_e %>% mutate(sig = p_fusiom2 < (0.05))
egenes = df_e$SYMBOL[df_e$sig]
egenes_combined = c(egenes_combined, egenes)
}

length(unique(egenes_combined))
length(intersect(unique(apagenes_combined), unique(egenes_combined)))


# Figure 3, heatmap
library(tidyverse); library(data.table); library(pheatmap); library(viridis); library(scales); library(colorspace);  library(grid)
tbl = fread("apa_sig_genes.txt")
# only keep AF sig genes in >=1 tissue
tbl = tbl %>% group_by(gene) %>% mutate(nsig1 = sum(sig1), nsig2 = sum(sig2), nisg_both = sum(sig1) + sum(sig2), strong1 = sum(p1m<0.05), strong2 = sum((p2m<0.05))) %>% dplyr::filter(nisg_both > 0, strong1  >1 | strong2 >1) 
df_wide = tbl %>% dplyr::select(gene,tissue,nlog10p1,nlog10p2) %>% pivot_wider(names_from = tissue, values_from = c(nlog10p1,nlog10p2)) 
# rank genes
summary(tbl$nlog10p1) 
df_wide = df_wide %>% rowwise() %>% mutate(n5_12 = sum(c_across(nlog10p1_Amygdala:nlog10p2_Atrial_Appendage) > 5, na.rm=T), n5_1 = sum(c_across(nlog10p1_Amygdala:nlog10p1_Atrial_Appendage) > 5, na.rm=T), n5_2 = sum(c_across(nlog10p2_Amygdala:nlog10p2_Atrial_Appendage) > 5, na.rm=T)) %>% ungroup() %>% arrange(desc(n5_12), desc(n5_1), desc(n5_2)) %>% column_to_rownames(var = "gene") 
# rank ct
mat = t(as.matrix(df_wide %>% dplyr::select(nlog10p1_Amygdala:nlog10p2_Atrial_Appendage)))

## clean + smooth scale BEFORE heatmap
# 1) replace Inf with NA first
mat[!is.finite(mat)] <- NA  # catches Inf, -Inf, NaN but keeps NA for now
# 2) compute quantiles on the vectorized matrix
v = as.numeric(mat)
mx = quantile(v, 0.98, na.rm = TRUE)
mn = quantile(v, 0.02, na.rm = TRUE)
# 3) cap extremes
mat[mat > mx] = mx
mat[mat < mn] = mn
# mat_raw is your original matrix with NA's
mat_raw = mat
# 4) now deal with NA: send them to the minimum (or 0 if you prefer)
mat[is.na(mat)] = mn   # or 0
dist_cols = dist(t(mat), method = "euclidean")
# 5) final sanity check: everything must be finite
stopifnot(all(is.finite(mat)))

# add disease group
annotation_row = data.frame(Disease = factor(rep(c("IS", "AF"), c(18,18))), Organ = factor(rep(c("Brain", "Blood", "Artery", "Heart", "Brain", "Blood", "Artery", "Heart"), c(13,1,2,2,13,1,2,2))))
# Add row names to match row names of the matrix
rownames(annotation_row) <- rownames(mat)
ann_colors <- list(Organ = lighten(c(Brain = "#efee00", Blood = "#ff00bb", Artery = "tomato", Heart = "#9900ff"), amount = 0.4), Disease = lighten(c(IS = "#006d2c", AF = "gray"), amount = 0.4))
names(ann_colors$Organ) = c("Brain", "Blood", "Artery", "Heart")
names(ann_colors$Disease) = c("IS", "AF")
custom_labels = sub("nlog10p1_", "", rownames(mat))
custom_labels = sub("nlog10p2_", "", custom_labels)
custom_labels = sub("_basal_ganglia", "_BG", custom_labels)

## build star matrix for significance (sig1 for IS, sig2 for AF)
# long table of unique gene–tissue significance
sig_long = tbl %>% distinct(gene, tissue, sig1, sig2)
# rows corresponding to IS (nlog10p1_) and AF (nlog10p2_)
sig_is = sig_long %>% transmute(row  = paste0("nlog10p1_", tissue), gene, star = ifelse(sig1, "*", ""))
sig_af = sig_long %>% transmute(row  = paste0("nlog10p2_", tissue), gene, star = ifelse(sig2, "*", ""))
# wide matrix with rows = rownames(mat) (nlog10p1_xxx / nlog10p2_xxx) and cols = genes
sig_mat_df = bind_rows(sig_is, sig_af) %>% pivot_wider(names_from = gene, values_from = star, values_fill = "") %>% column_to_rownames("row")
# align to heatmap matrix order
sig_mat = sig_mat_df[rownames(mat), colnames(mat)] %>% as.matrix()
sig_mat[is.na(sig_mat)] = ""   

# plot
n_IS <- sum(grepl("^nlog10p1_", rownames(mat)))
ph = pheatmap(mat_raw, cluster_rows = F, cluster_cols = T, 
clustering_distance_cols = dist_cols,
na_col = "gray95",    
clustering_method = "complete", 
#clustering_distance_cols = "euclidean", 
clustering_distance_rows = "manhattan", fontsize_col = 7, fontsize_row = 7, gaps_row = n_IS, color = colorRampPalette(c("white", "royalblue", "blue1"))(100), main = "", annotation_row = annotation_row, show_colnames = T, labels_row = custom_labels, treeheight_row = 0, treeheight_col = 0, border_color = NA, annotation_colors = ann_colors, display_numbers = sig_mat, number_color  = "black",fontsize_number = 10)

png("heatmap_apa.png", width = 1600, height = 800, res = 150)
grid.newpage()
grid.draw(ph$gtable)
grid.text(
expression(-log[10](p)),
x = 0.86, y = 0.94,
gp = gpar(fontsize = 10)
)
dev.off()

                  

            
            

