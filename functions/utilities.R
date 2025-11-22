bern_prob <- function(rho, p, q) {
p00 = rho * sqrt(p*q*(1-p)*(1-q)) + (1-p)*(1-q)
p10 = 1 - q - p00
p01 = 1 - p - p00
p11 = p + q + p00 -  1
p00 = max(0, p00); p10 = max(0, p10); p01 = max(0, p01); p11 = max(0, p11)
return(c(p00, p01, p10, p11))
}

# handle flip

label_flip = function(niter, res_prop, eta_true) {
ids = (niter/2 + 1):niter
eta_post_mode = apply(res_prop$eta_tk[ids, ], 2, function(x) { uniqx = unique(x); uniqx[which.max(tabulate(match(x, uniqx)))] }) # estimated eta, posterior mode 
qq = mean(res_prop$q_tk[ids]) 
if (qq > 0.5) {
eta_post = 1 - eta_post_mode
b_mean = mean(res_prop$beta_tk[ids] + res_prop$alpha_tk[ids])
b_sd = sd(res_prop$beta_tk[ids] + res_prop$alpha_tk[ids])
ab_mean = mean(res_prop$beta_tk[ids])
ab_sd = sd(res_prop$beta_tk[ids])
bci = quantile(res_prop$beta_tk[ids] + res_prop$alpha_tk[ids], probs=c(0.025,0.975), na.rm=T)
}
if (qq <=0.5) {
eta_post = eta_post_mode
b_mean = mean(res_prop$beta_tk[ids])
b_sd = sd(res_prop$beta_tk[ids])
ab_mean = mean(res_prop$alpha_tk[ids] + res_prop$beta_tk[ids])
ab_sd = sd(res_prop$alpha_tk[ids] + res_prop$beta_tk[ids])
bci = quantile(res_prop$beta_tk[ids], probs=c(0.025,0.975), na.rm=T)
}
eta_match = mean(eta_post == eta_true)
eta_est = eta_post

return(list(qq = qq, eta_match = eta_match, eta_est = eta_est, b_mean = b_mean, b_sd = b_sd, bci = bci))

}

# handle flip - joint

label_flip_joint = function(niter, res, eta_true_1, eta_true_2) {
  ids = (niter/2 + 1):niter
  p00_post = res$pst_tk[ids, 1]
  p01_post = res$pst_tk[ids, 2]
  p10_post = res$pst_tk[ids, 3]
  p11_post = res$pst_tk[ids, 4]
  eta_1_post_mode = apply(res$eta_1_tk[ids, ], 2, function(x) { uniqx = unique(x); uniqx[which.max(tabulate(match(x, uniqx)))] })
  eta_2_post_mode = apply(res$eta_2_tk[ids, ], 2, function(x) { uniqx = unique(x); uniqx[which.max(tabulate(match(x, uniqx)))] })
  qq1 = mean(p10_post + p11_post) 
  qq2 = mean(p01_post + p11_post)
  if (qq1 > 0.5) {
  eta_1_post = 1 - eta_1_post_mode
  b1_mean = mean(res$beta_1_tk[ids] + res$alpha_1_tk[ids])
  b1_sd = sd(res$beta_1_tk[ids] + res$alpha_1_tk[ids])
  bci1 = quantile(res$beta_1_tk[ids] + res$alpha_1_tk[ids], probs=c(0.025,0.975), na.rm=T)
  }
  if (qq1 <=0.5) {
  eta_1_post = eta_1_post_mode
  b1_mean = mean(res$beta_1_tk[ids])
  b1_sd = sd(res$beta_1_tk[ids])
  bci1 = quantile(res$beta_1_tk[ids], probs=c(0.025,0.975), na.rm=T)
  }
  if (qq2 > 0.5) {
  eta_2_post = 1 - eta_2_post_mode
  b2_mean = mean(res$beta_2_tk[ids] + res$alpha_2_tk[ids])
  b2_sd = sd(res$beta_2_tk[ids] + res$alpha_2_tk[ids])
  bci2 = quantile(res$beta_2_tk[ids] + res$alpha_2_tk[ids], probs=c(0.025,0.975), na.rm=T)
  }
  if (qq2 <=0.5) {
  eta_2_post = eta_2_post_mode
  b2_mean = mean(res$beta_2_tk[ids])
  b2_sd = sd(res$beta_2_tk[ids])
  bci2 = quantile(res$beta_2_tk[ids], probs=c(0.025,0.975), na.rm=T)
  }
  
  eta_1_match = mean(eta_1_post == eta_true_1)
  eta_2_match = mean(eta_2_post == eta_true_2)
  eta_1_est = eta_1_post
  eta_1_est = eta_1_post
  
  return(list(qq1 = qq1, qq2 = qq2, eta_1_match = eta_1_match, eta_2_match = eta_2_match, 
  eta_1_est = eta_1_post, eta_2_est = eta_2_post, b1_mean = b1_mean, b1_sd = b1_sd, b2_mean = b2_mean, b2_sd = b2_sd, bci1 = bci1, bci2 = bci2))
  
}

clump <- function(dat, SNP_col = "eQTL_variant_id", pval_col = "rowmeta", clump_kb = 250, clump_r2 = 0.1, clump_p = 0.999, bfile = "/scratch/t.phs.yihaolu/GWAS_summary_statistics/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.MAF01", plink_bin = "plink", pop="EUR") {
df <- data.frame(rsid = dat[, ..SNP_col], pval = dat[,..pval_col])
colnames(df) = c("rsid", "pval")
out <- tryCatch({
ieugwasr::ld_clump(df, clump_kb=clump_kb, clump_r2=clump_r2, clump_p=clump_p, bfile=bfile, plink_bin = plink_bin, pop = pop)
}, silent = TRUE, error = function(x) return(NA)
)
if(length(out)==1) return(NA)
MRdat <- dat[which(unlist(dat[,..SNP_col]) %in% out$rsid),]
return(MRdat)
}

