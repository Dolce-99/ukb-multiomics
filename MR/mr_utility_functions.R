library(ggplot2)
library(stringr)

my_theme <- theme_bw() +
  # set transparency
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA),
    legend.position = "none"
  )


egger_results_to_df <- function(mr_egger_results, presso_corrected = F){
  if(presso_corrected){
    estimand_vector = c("PRESSO Corrected\nEgger",
                        "PRESSO Corrected\nEgger-Int.")
  } else{
    estimand_vector <- c("Egger", 
                         "Egger-Int.")
  }
  tibble(
    estimand = estimand_vector,
    est = c(mr_egger_results@Estimate, mr_egger_results@Intercept),
    est_ci_lower = c(mr_egger_results@CILower.Est, mr_egger_results@CILower.Int),
    est_ci_upper = c(mr_egger_results@CIUpper.Est,mr_egger_results@CIUpper.Int), 
    pval = c(mr_egger_results@Pvalue.Est, mr_egger_results@Pvalue.Int)
  ) %>% 
    return()
}

ivw_results_to_df <- function(ivw_results){
  tibble(
    estimand = c("IVW"),
    est = c(ivw_results@Estimate),
    est_ci_lower = c(ivw_results@CILower),
    est_ci_upper = c(ivw_results@CIUpper), 
    pval = c(ivw_results@Pvalue)
  )
}

median_results_to_df <- function(median_results){
  tibble(
    estimand = c("Weighted Median"),
    est = c(median_results@Estimate),
    est_ci_lower = c(median_results@CILower),
    est_ci_upper = c(median_results@CIUpper), 
    pval = c(median_results@Pvalue)
  )
}

presso_ivw_results_to_df <- function(presso_results){
  tibble(
    estimand = c("MR-PRESSO Outlier\nw/o Corrected IVW"), 
    est = presso_results$`Main MR results`[1,3], 
    est_ci_lower = presso_results$`Main MR results`[1,3] - 1.96*presso_results$`Main MR results`[1,4],
    est_ci_upper = presso_results$`Main MR results`[1,3] + 1.96*presso_results$`Main MR results`[1,4],
    pval = presso_results$`Main MR results`[1,6]
    
  )
}

combine_mr_results <- function(iwv_results, 
                               mr_egger_results, 
                               mr_median_results,
                               presso_results){
  bind_rows(
    ivw_results_to_df(iwv_results),
    egger_results_to_df(mr_egger_results),
    median_results_to_df(mr_median_results),
    presso_ivw_results_to_df(presso_results)
  ) %>% 
    mutate(nom_sig = if_else(pval < 0.05, "sig", "not"), 
           nom_sig = case_when(
             str_detect(estimand, "Int.") & nom_sig == "not" ~ "sig",
             str_detect(estimand, "Int.") & nom_sig == "sig" ~ "not",
             TRUE ~ nom_sig
             )
           )
}

combine_mr_results_nopresso <- function(iwv_results, 
                               mr_egger_results, 
                               mr_median_results){
  bind_rows(
    ivw_results_to_df(iwv_results),
    egger_results_to_df(mr_egger_results),
    median_results_to_df(mr_median_results),
  ) %>% 
    mutate(nom_sig = if_else(pval < 0.05, "sig", "not"), 
           nom_sig = case_when(
             str_detect(estimand, "Int.") & nom_sig == "not" ~ "sig",
             str_detect(estimand, "Int.") & nom_sig == "sig" ~ "not",
             TRUE ~ nom_sig
             )
           )
}

plot_mr_results <- function(iwv_results, 
                            mr_egger_results, 
                            presso_results, 
                            presso_egger_results){
  combine_mr_results(iwv_results, 
                     mr_egger_results, 
                     presso_results, 
                     presso_egger_results) %>% 
    ggplot(aes(x = est, y = estimand, alpha = nom_sig)) + 
    geom_point(aes(color = nom_sig)) + 
    geom_errorbar(aes(xmin = est_ci_lower,
                      xmax = est_ci_upper, 
                      color = nom_sig)) + 
    geom_vline(xintercept = 0, linetype = "dashed", 
               alpha= 0.5) + 
    labs(x = "Estimate and 95% CI", 
         y = "Estimand") +
    scale_color_manual(values=c("#000000", "#27A4F2")) +
    scale_alpha_manual(values=c(0.25, 1)) + my_theme
}
