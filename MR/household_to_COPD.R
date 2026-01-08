library(dplyr)
library(data.table)
library(MendelianRandomization)
library(MRPRESSO)
library(janitor)
library(stringr)
library(parallel)
source("mr_utility_functions.R")

# outnames for for COPD: J10_COPD hg38 based
out.name = "J10_COPD"
process.mr = function(out.name){
    indep.vars = NULL
    for (i in c(1:20)){
        path = paste0("household_indep_hg38_chr",i)
        if (file.exists(path)){
            tmp = fread(path)
            indep.vars = rbind(indep.vars, tmp)
        }
    }

    gwas = fread(paste0("finngen_R12_",out.name,".gz"))
    final = indep.vars %>% left_join(gwas, by = c("SNP" = "rsids")) %>%
            rename(out.beta = beta, out.se = sebeta) %>% filter(pval > 5e-8)
    
    mr_input_obj = mr_input(
        bx = final$Beta,
        bxse = final$Standard_Error_of_Beta,
        by = final$out.beta,
        byse = final$out.se
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
