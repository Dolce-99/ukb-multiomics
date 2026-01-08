# use somascan to validate PRG3 to skin cancer MR result
library(dplyr)
library(data.table)
library(MendelianRandomization)
library(MRPRESSO)
library(janitor)
library(stringr)
library(parallel)
library(readxl)
source("mr_utility_functions.R")

# outnames for skin cancer: C3_OTHER_SKIN_EXALLC
# both are hg38 based

somascan = read_excel("41588_2021_978_MOESM4_ESM.xlsx", sheet = "ST02", skip=1)
prg3.vars = somascan[which(somascan[[3]] == "PRG3"), ]
indep.vars = prg3.vars %>% select(`pQTL_ID \r\n(global)`, variant, `pos\r\n(var.)`, `beta\r\n(adj.)`, `-Log10(P)\r\n(adj.)`) %>% 
                mutate(chr = as.integer(sub("chr","", sub("_.*", "", `pQTL_ID \r\n(global)`))), BETA = `beta\r\n(adj.)`, P = 10^(-`-Log10(P)\r\n(adj.)`), SE = abs(BETA / qnorm(P/2)))

out.name = "C3_OTHER_SKIN_EXALLC"
gwas = fread(paste0("/proj/yunligrp/users/djw/UKB_multiomics/protein_contextualize/GWASsumstats/finngen_R12_",out.name,".gz"))
final = indep.vars %>% left_join(gwas, by=c("chr" = "#chrom", "pos\r\n(var.)" = "pos")) %>%
        rename(out.beta = beta, out.se = sebeta) %>% filter(pval > 5e-8) 
    
mr_input_obj = mr_input(
    bx = final$BETA,
    bxse = final$SE,
    by = final$out.beta,
    byse = final$out.se,
    exposure = "PRG3",
    outcome = "Skin cancer"
)
ivw_result <- mr_ivw(mr_input_obj)
egger_result <- mr_egger(mr_input_obj)
median_result <- mr_median(mr_input_obj)

IVW.df = ivw_results_to_df(ivw_result)
egger.df = egger_results_to_df(egger_result)
median.df = median_results_to_df(median_result)

combined.df = bind_rows(IVW.df, egger.df, median.df)
print(combined.df)
