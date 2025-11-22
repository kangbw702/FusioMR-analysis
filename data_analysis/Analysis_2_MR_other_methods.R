# module load R/4.5.1

library(data.table)
library(tidyverse)
library(MendelianRandomization)
library(MR2)

args <- commandArgs(trailingOnly = T)
tissue_id = as.numeric(args[1])
chr = as.numeric(args[2])

print(tissue_id)
print(chr)

#===============================================================================================================================================

# j: index of the target gene 
run_mr_others <- function(tissue, chr, j, df_clump_list, niv_cut=5, out_dir) {
gg = gene_list[j]
df = df_clump_list[[gg]]
b_exp = df$beta; se_exp = df$se
b_out_1 = df$b_out_1; se_out_1 = df$se_out_1
b_out_2 = df$b_out_2; se_out_2 = df$se_out_2
K_post = nrow(df)
cat("Processing gene", j, gg, "out of", length(df_clump_list), K_post, "IVs to process\n")

if (K_post >= niv_cut) {
# run MR
# competing methods (single-outcome MR, outcome 2)
mr.obj = MendelianRandomization::mr_input(bx = b_exp, bxse = se_exp, by = b_out_2, byse = se_out_2)
# IVW fixed
IVW_f = MendelianRandomization::mr_ivw(mr.obj, model = 'fixed')
b_ivw = IVW_f$Estimate; se_ivw = IVW_f$StdError; p_ivw = IVW_f$Pvalue
# MR-Egger
Egger = try(MendelianRandomization::mr_egger(mr.obj))
if(class(Egger) != 'try-error' & !is.null(Egger)) {
b_egger = Egger$Estimate; se_egger = Egger$StdError.Est; p_egger = Egger$Causal.pval
} else {b_egger = se_egger = p_egger = NA}
# cML
cml = try(MendelianRandomization::mr_cML(mr.obj, MA = T, DP = F, num_pert = 200, n = 500))
if(class(cml) != 'try-error' & !is.null(cml)) {
b_cml = cml$Estimate; se_cml = cml$StdError; p_cml = cml$Pvalue
} else {b_cml = se_cml = p_cml = NA}
# cML-DP
cml_dp = try(MendelianRandomization::mr_cML(mr.obj, MA = T, DP = T, num_pert = 200, n = 500))
if(class(cml_dp) != 'try-error' & !is.null(cml_dp)) {
b_cml_dp = cml_dp$Estimate; se_cml_dp = cml_dp$StdError; p_cml_dp = cml_dp$Pvalue
} else {b_cml_dp = se_cml_dp = p_cml_dp = NA}

# MR2, multi-response
betaHat_Y = as.matrix(cbind(b_out_1, b_out_2))
betaHat_X = as.matrix(cbind(b_exp))
colnames(betaHat_Y) = colnames(betaHat_X) = NULL
MR2_output = MR2(betaHat_Y, betaHat_X, EVgamma = 0.5, niter = 7500, burnin = 2500, thin = 5, monitor = 1000)
# head(MR2_output$postMean$theta)
p_MR2_1 = MR2_output$samplerPar$pval[1]
p_MR2_2 = MR2_output$samplerPar$pval[2]
PostProc_output = PostProc(MR2_output, betaHat_Y, betaHat_X)
b_MR2_1 = PostProc_output$thetaPost[1]
b_MR2_2 = PostProc_output$thetaPost[2]
#bci_MR2_1 = PostProc_output$thetaPost_CI[1,1,]
#bci_MR2_2 = PostProc_output$thetaPost_CI[1,2,]

# write out
tmp = c(tissue, chr, gg, K_post, b_ivw, se_ivw, p_ivw, b_egger, se_egger, p_egger, b_cml, se_cml, p_cml, b_cml_dp, se_cml_dp, p_cml_dp, b_MR2_1, p_MR2_1, b_MR2_2, p_MR2_2)
write.table(t(tmp), paste0(out_dir, '/others_res_', tissue, '.txt'), quote=F, row.names=F, col.names=F, append=T)
} # end if
}

#===============================================================================================================================================

pcut = 0.001; clump_kb = 50; clump_r2 = 0.01
niv_cut = 5
out_dir = "res_r1"

traits = c("stroke", "af") 
tissue_short = c("Amygdala","BA24","Caudate_basal_ganglia","Cerebellar_Hemi","Cerebellum","Cortex","BA9","Hippocampus","Hypothalamus","Nucleus_accumbens_basal_ganglia","Putamen_basal_ganglia","Spinal_cord","Substantia_nigra","Blood","Aorta","Coronary","Left_Ventricle","Atrial_Appendage")

tissue = tissue_short[tissue_id]
df_clump_list = readRDS(paste0('df_clump/dfclump_', pcut, '_', clump_kb, '_', clump_r2, '_', tissue, '_', chr, '.rds')) 
gene_list = names(df_clump_list)
cat(tissue, "chr", chr, ",", length(gene_list), "genes are tested.\n")
start_g = 1; end_g = length(gene_list)
#start_g = 106; end_g = 106

for (j in start_g:end_g) {
cat(j, "out of", length(gene_list), "\n")
run_mr_others(tissue, chr, j, df_clump_list, niv_cut, out_dir)
}

cat("\nFinish Others for", tissue, "chr", chr, "genes", start_g, "-", end_g, "\n")




