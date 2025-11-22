library(data.table); library(tidyverse)
snps = fread("bryois/snp_pos.txt.gz", select = c("SNP", "effect_allele", "other_allele"))
ctypes = c("Endothelial.cells","Astrocytes","Excitatory.neurons","Inhibitory.neurons","Microglia","OPCs...COPs","Oligodendrocytes","Pericytes")
ctypes_short = c("end","ast","exc","inh","mic","opc","oli","per")
disease = c("dcfdx_ad.assoc.logistic")
disease_short = c("ad")

args <- commandArgs(trailingOnly = T)
chr = as.numeric(args[1])
print(chr)

out = fread(paste0("../ad_educ/ad_data/", disease, ".gz_all.txt")) %>% dplyr::select(-pos) %>% dplyr::rename(b_out=beta, se_out=se)

for (i in 1:length(ctypes)) {
print(paste0(i, ", ", ctypes[i]))
# colnames: Gene_id, SNP_id, Distance to TSS, Nominal p-value, Beta
eqtl = fread(paste0("../single_cell_eQTL/bryois/", ctypes[i], ".", chr, ".gz")) %>% 
inner_join(snps, by = c("V2"="SNP")) %>% drop_na() %>% separate(V1, into = c("gene_symbol", "gene_name"), sep = "_") %>% 
dplyr::rename(rsid=V2, p=V4, beta=V5) %>% mutate(se = abs(beta/qnorm(p/2, lower.tail = F))) %>% 
dplyr::select(-V3)
tmp = eqtl %>% left_join(out, by = c('rsid'='rsid')) %>% 
mutate(align = ifelse(effect_allele==EFF & other_allele==REF, 1, ifelse(effect_allele==REF & other_allele==EFF, -1, NA))) %>% 
dplyr::filter(!(is.na(align) & !(is.na(b_out)))) %>% mutate(b_out = b_out*align) %>% drop_na(b_out, se_out) %>% 
dplyr::select(-REF, -EFF, -align, -chr)
fwrite(tmp, paste0('harmo_out/eqtl_', disease_short, '_', ctypes_short[i], "_", chr, '.txt.gz'))
} # end loop cell types




