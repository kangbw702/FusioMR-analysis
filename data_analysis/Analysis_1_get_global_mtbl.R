# calculate global mom R/4.2.1

library(data.table)
library(tidyverse)
library(pbapply)

source("functions/set_var_prior.R")

#==============================================================================================

# Build a global EB center from many genes' local MoMs (vector),
# using a trimmed median of positive values (robust, no recomputation inside priors).
build_global_mom <- function(local_mom_vec, q = 0.5, trim = 0.05, min_pos = 20) {
x = local_mom_vec[is.finite(local_mom_vec) & local_mom_vec > 0]
if (length(x) == 0) return(0)
if (length(x) < min_pos) return(stats::median(x, na.rm = TRUE))
lo = stats::quantile(x, trim, na.rm = TRUE)
hi = stats::quantile(x, 1 - trim, na.rm = TRUE)
x = x[x >= lo & x <= hi]
if (!length(x)) return(0)
stats::quantile(x, q, na.rm = TRUE)
}

#==============================================================================================

pcut = 0.001; clump_kb = 50; clump_r2 = 0.1
glob_q = 0.9
ctypes_short = c("end","ast","exc","inh","mic","opc","oli","per")
diseases_short = "ad"

for (i in 1:length(ctypes_short)) {
print(paste0("ctype-", i, ", ", ctypes_short[i]))
# calculate global moms
res = NULL
for (chr in 1:22) {
cat(chr, "\n")
df_clump_list = readRDS(paste0("df_clump/dfclump_", pcut, "_", clump_kb, "_", clump_r2, "_", ctypes_short[i], "_", chr, '.rds'))

tbl = do.call(rbind, pblapply(df_clump_list, function(df) {
ghat = df$beta; gse = df$se; Ghat = df$b_out; Gse = df$se_out
moms = set_variance_priors(ghat, gse, Ghat, Gse, beta0 = NULL, K = NULL, Kmin = 5, Kmax = 20, 
rho_ov = 0, c_gamma = 0.5, c_theta = 0.8, global_mean_gamma = NULL, global_mean_theta = NULL, 
hybrid = FALSE, kappa_hybrid = 5, z_thresh = NULL, trim = 0.1)
c(mom_gamma = moms$gamma$mom_local, mom_theta = moms$theta$mom_local, beta0 = moms$beta0, K = moms$K)
}))

tbl = as.data.frame(tbl)
tbl$gene_symbol = names(df_clump_list)
tbl$ctype = ctypes_short[i]
rownames(tbl) = NULL
res = bind_rows(res, tbl) 
} # end loop chr

# Globals per cell type using the same aux pipeline:
# trimmed median of positive locals
global_mean_gamma = build_global_mom(res$mom_gamma, q = glob_q, trim = 0.05, min_pos = 20)  
global_mean_theta = build_global_mom(res$mom_theta, q = glob_q, trim = 0.05, min_pos = 20)
res$global_m_gamma = global_mean_gamma
res$global_m_theta = global_mean_theta

fwrite(res, paste0("mtbls/mtbls_testclump_", diseases_short, "_", glob_q, "_", pcut, "_", clump_kb, "_", clump_r2, "_", ctypes_short[i], ".txt"))
} # end loop cell type

cat("step3, calculating local and global mom from aux clump df, finished!\n")




