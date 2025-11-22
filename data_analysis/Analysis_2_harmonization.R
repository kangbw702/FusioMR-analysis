library(data.table); library(tidyverse)
args <- commandArgs(trailingOnly = T)
chr = as.numeric(args[1])
print(chr)

tissues = c("Brain_Amygdala","Brain_Anterior_cingulate_cortex_BA24","Brain_Caudate_basal_ganglia","Brain_Cerebellar_Hemisphere","Brain_Cerebellum","Brain_Cortex","Brain_Frontal_Cortex_BA9","Brain_Hippocampus","Brain_Hypothalamus","Brain_Nucleus_accumbens_basal_ganglia","Brain_Putamen_basal_ganglia","Brain_Spinal_cord_cervical_c-1","Brain_Substantia_nigra","Whole_Blood","Artery_Aorta","Artery_Coronary","Artery_Tibial","Heart_Left_Ventricle","Heart_Atrial_Appendage") # 18
tissue_short = c("Amygdala","BA24","Caudate_basal_ganglia","Cerebellar_Hemi","Cerebellum","Cortex","BA9","Hippocampus","Hypothalamus","Nucleus_accumbens_basal_ganglia","Putamen_basal_ganglia","Spinal_cord","Substantia_nigra","Blood","Aorta","Coronary","Left_Ventricle","Atrial_Appendage")
traits = c("stroke", "af")

# read in out gwas
out1 = fread("gwas/stroke.tsv.gz") %>% dplyr::select(rsid,effect_allele:standard_error) %>% drop_na() %>% dplyr::rename(b_out_1=beta, se_out_1=standard_error) %>% dplyr::filter(se_out_1!=Inf, se_out_1!=0) %>% distinct(rsid, .keep_all=T) # 7,465,328 snps 
out2 = fread("gwas/af.Nielsen.tsv.gz") %>% dplyr::select(-hm_chrom, -hm_pos, -p_value) %>% drop_na() %>% dplyr::rename(rsid=hm_rsid, effect_allele=hm_effect_allele, other_allele=hm_other_allele, b_out_2=hm_beta, se_out_2=standard_error) %>% dplyr::filter(se_out_2!=Inf, se_out_2!=0) %>% distinct(rsid, .keep_all=T) # 32,555,902 snps
# harmonize two outcome gwas
out = out1 %>% left_join(out2, by = 'rsid') %>% mutate(align = ifelse(effect_allele.x==effect_allele.y & other_allele.x==other_allele.y, 1, ifelse(effect_allele.x==other_allele.y & other_allele.x==effect_allele.y, -1, NA))) %>% dplyr::filter(!(is.na(align) & !(is.na(b_out_2)))) %>% mutate(b_out_2 = b_out_2*align) %>% drop_na(b_out_2, se_out_2) %>% dplyr::select(-other_allele.y, -effect_allele.y, -align) %>% dplyr::rename(EFF=effect_allele.x, REF=other_allele.x) # 7,429,631 snps
for (j in 1:length(tissues)) {
print(paste0(j, ", ", tissues[j]))
apaqtl = fread(paste0("apaqtl/", tissues[j], ".cis_3aQTL.txt.gz")) %>% dplyr::filter(chromosome==paste0("chr",chr)) %>% separate(gene, into = c("gene_name", "gene_symbol", "chr", "direction"), sep = "\\|") %>% dplyr::select(rsid, gene_symbol, beta:se) %>% drop_na() %>% dplyr::filter(se!=Inf, se!=0, p!=1)
# harmonize apaqtl and gwas
tmp = apaqtl %>% left_join(out, by = 'rsid') %>% mutate(align = ifelse(alt==EFF & ref==REF, 1, ifelse(alt==REF & ref==EFF, -1, NA))) %>% dplyr::filter(!(is.na(align) & !(is.na(b_out_1)) & !(is.na(b_out_2)))) %>% mutate(b_out_1 = b_out_1*align, b_out_2 = b_out_2*align) %>% drop_na(b_out_1, b_out_2, se_out_1, se_out_2) %>% dplyr::select(-REF, -EFF, -align) 
# write out
fwrite(tmp, paste0('harmo/apaqtl_', traits[1], '_', traits[2], '_', tissue_short[j], '_', chr, '.txt.gz'))
#} # end loop






