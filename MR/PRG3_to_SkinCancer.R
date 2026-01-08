# PRG3 (UKB) -> skin cancer
library(dplyr)
library(data.table)
library(MendelianRandomization)
library(MRPRESSO)
library(janitor)
library(stringr)
library(parallel)
source("mr_utility_functions.R")

process.mr = function(exp.name, out.name){
    indep.vars = NULL
    for (i in c(1:22)){
        path = paste0(exp.name,"_variants/",exp.name,"_chr",i,".significant.snps.indep.txt")
        if (file.exists(path)){
            tmp = fread(path)
            indep.vars = rbind(indep.vars, tmp)
        }
    }
    gwas = fread(paste0("finngen_R12_",out.name,".gz"))
    final = indep.vars %>% left_join(gwas, by=c("CHROM" = "#chrom", "GENPOS" = "pos", "ALLELE0" = "ref", "ALLELE1" = "alt")) %>%
            rename(out.beta = beta, out.se = sebeta) %>% filter(pval > 5e-8) # remove those significant in outcome

    mr_input_obj = mr_input(
        bx = final$BETA,
        bxse = final$SE,
        by = final$out.beta,
        byse = final$out.se,
        exposure = exp.name,
        outcome = out.name
    )
    ivw_result <- mr_ivw(mr_input_obj)
    egger_result <- mr_egger(mr_input_obj)
    median_result <- mr_median(mr_input_obj)

    IVW.df = ivw_results_to_df(ivw_result)
    egger.df = egger_results_to_df(egger_result)
    median.df = median_results_to_df(median_result)

    combined.df = bind_rows(IVW.df, egger.df, median.df)
    print(paste0("outcome = ",out.name))
    print(combined.df)
    return(combined.df)
}

prg3.df = process.mr("PRG3", "C3_OTHER_SKIN_EXALLC")

