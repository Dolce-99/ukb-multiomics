library(dplyr)
library(data.table)
library(ggplot2)
library(tidyr)

t_test = function(covarset){
    df = fread(paste0("all_results_",covarset,"_5rs_tables.txt"))
    df$model <- gsub(" ", "_", df$model)
    df$model <- gsub("-", "_", df$model)
    df = df %>% pivot_wider(names_from = model, values_from = value)
    res = NULL
    for (t in unique(df$Trait)){
        tmp = df %>% filter(Trait == t)
        diff1 = mean(tmp$Metabolomics_only) - mean(tmp$Baseline_covariates)
        pval1 = t.test(tmp$Metabolomics_only, tmp$Baseline_covariates, paired = TRUE)$p.value
        diff2 = mean(tmp$Proteomics_only) - mean(tmp$Baseline_covariates)
        pval2 = t.test(tmp$Proteomics_only, tmp$Baseline_covariates, paired = TRUE)$p.value
        diff3 = mean(tmp$Combined) - mean(tmp$Baseline_covariates)
        pval3 = t.test(tmp$Combined, tmp$Baseline_covariates, paired = TRUE)$p.value
        diff4 = mean(tmp$Proteomics_only) - mean(tmp$Metabolomics_only)
        pval4 = t.test(tmp$Proteomics_only, tmp$Metabolomics_only, paired = TRUE)$p.value
        diff5 = mean(tmp$Combined) - mean(tmp$Metabolomics_only)
        pval5 = t.test(tmp$Combined, tmp$Metabolomics_only, paired = TRUE)$p.value
        diff6 = mean(tmp$Combined) - mean(tmp$Proteomics_only)
        pval6 = t.test(tmp$Combined, tmp$Proteomics_only, paired = TRUE)$p.value
        res = rbind(res, data.frame(
            Trait = t,
            C.covar = mean(tmp$Baseline_covariates),
            C.meta  = mean(tmp$Metabolomics_only),
            C.prote = mean(tmp$Proteomics_only),
            C.mp    = mean(tmp$Combined),
            diff1 = diff1, pval1 = pval1,
            diff2 = diff2, pval2 = pval2,
            diff3 = diff3, pval3 = pval3,
            diff4 = diff4, pval4 = pval4,
            diff5 = diff5, pval5 = pval5,
            diff6 = diff6, pval6 = pval6
        ))
    }
    return(res)
}
res1 = t_test("age_sex")
res2 = t_test("ASCVD")
res3 = t_test("PANEL")

save_function = function(df, covarset){
    df = df %>% mutate(Predictor_set = covarset) %>%
        rename(`C.index_baseline` = C.covar,
        `C.index_metabolomic-sonly` = C.meta,
        `C.index_proteomics-only` = C.prote,
        `C.index_combined` = C.mp,
        `C.index_diff_metabolomics_vs_baseline` = diff1,
        `pval_metabolomics_vs_baseline` = pval1,
        `C.index_diff_proteomics_vs_baseline` = diff2,
        `pval_proteomics_vs_baseline` = pval2,
        `C.index_diff_combined_vs_baseline` = diff3,
        `pval_combined_vs_baseline` = pval3,
        `C.index_diff_proteomics_vs_metabolomics` = diff4,
        `pval_proteomics_vs_metabolomics` = pval4,
        `C.index_diff_combined_vs_metabolomics` = diff5,
        `pval_combined_vs_metabolomics` = pval5,
        `C.index_diff_combined_vs_proteomics` = diff6,
        `pval_combined_vs_proteomics` = pval6
        )
    return(df)
}
res1.out = save_function(res1, "Age+Sex")
res2.out = save_function(res2, "ASCVD")
res3.out = save_function(res3, "PANEL")
allres = rbind(res1.out, res2.out, res3.out)
trait_order <- c("MACE",'Heart_Failure','CHD','PAD', 'Atrial_Fibrillation', 'T2_Diabetes',
                 'Prostate_Cancer', 'Skin_Cancer', 'Breast_Cancer', 'Dementia', 
                 'Cataracts', 'Glaucoma', 'COPD', 'Asthma', 'Renal_Disease', 
                 'Liver_Disease', 'Fractures')
pred_order <- c("Age+Sex", "ASCVD", "PANEL")
allres = allres %>% select(Trait, Predictor_set, starts_with("C.index_"),starts_with("pval_")) %>%
        mutate(across(starts_with("C.index_"), ~round(.,4))) %>%
        mutate(across(starts_with("pval_"), ~signif(.,4))) %>%
        arrange(match(Trait, trait_order),match(Predictor_set, pred_order))
write.table(as.data.frame(allres), "04.new.perf/ttest_results.txt", sep = "\t", quote = F, row.names = F)


# draw figure that columns are 17 traits and rows are 3*3
# color is mean C-index increment
process_df = function(inputdf, covar){
    output_df = inputdf %>% select(Trait, C.meta,C.prote,C.mp, pval4, pval5, pval6) %>% 
            mutate(best = case_when(C.meta > C.prote & C.meta > C.mp ~ "meta",
                                    C.prote > C.meta & C.prote > C.mp ~ "prote",
                                    C.mp > C.meta & C.mp > C.prote ~ "mp"),
                    second_best = case_when(
                        (C.meta > C.prote & C.meta < C.mp) | (C.meta < C.prote & C.meta > C.mp) ~ "meta",
                        (C.prote > C.meta & C.prote < C.mp) | (C.prote < C.meta & C.prote > C.mp) ~ "prote",
                        (C.mp > C.meta & C.mp < C.prote) | (C.mp < C.meta & C.mp > C.prote) ~ "mp"
                    ),
                    p_best2 = case_when(
                        (best=="meta"  & second_best=="prote") |
                        (best=="prote" & second_best=="meta")   ~ pval4, 
                        
                        (best=="meta"  & second_best=="mp")    |
                        (best=="mp"    & second_best=="meta")   ~ pval5,
                        
                        (best=="prote" & second_best=="mp")    |
                        (best=="mp"    & second_best=="prote")  ~ pval6 
                        )) %>% 
        select(Trait, C.meta, C.prote, C.mp, best, second_best, p_best2) %>%
        pivot_longer(cols = starts_with("C."), names_to = "type", names_prefix= "C.", values_to = "C.index") %>%
        mutate(
            is_best   = (type == best),
            is_second = (type == second_best & p_best2 > 0.05),
            covarset = covar
        )
    return(output_df)
}
plotdf1 = process_df(res1, "Age+Sex")
plotdf2 = process_df(res2, "ASCVD")
plotdf3 = process_df(res3, "PANEL")
plot_df <- bind_rows(plotdf1, plotdf2, plotdf3)

plot_df$type = recode(plot_df$type,
    "meta" = "Metabolomics",
    "prote" = "Proteomics",
    "mp" = "Combined"
)
plot_df$type = factor(plot_df$type, levels = c("Metabolomics", "Proteomics", "Combined"))
plot_df$Trait = factor(plot_df$Trait, levels = c(
    "MACE",'Heart_Failure','CHD','PAD', 'Atrial_Fibrillation', 'T2_Diabetes',
    'Prostate_Cancer', 'Skin_Cancer', 'Breast_Cancer',
    'Dementia', 
    'Cataracts',  'Glaucoma',
    'COPD', 'Asthma', 
    'Renal_Disease', 'Liver_Disease', 'Fractures'
))
plot_df = plot_df %>% mutate(num.type = match(type, c("Metabolomics", "Proteomics", "Combined")),
                            num.covar = match(covarset, c("Age+Sex","ASCVD", "PANEL")),
                            y.order = num.type + (num.covar - 1) * (3+0.3)) 
sec_axis_data = plot_df %>%
            group_by(covarset) %>%
            summarise(y.center = mean(y.order), .groups = 'drop')
minC = min(plot_df$C.index, na.rm = TRUE)
maxC = max(plot_df$C.index, na.rm = TRUE)
meanC = mean(plot_df$C.index, na.rm = TRUE)
p = ggplot(plot_df, aes(x = Trait, y = y.order, fill = C.index)) +
    geom_tile(color = "lightgrey", linewidth = 1) +
    geom_text(data = filter(plot_df, is_best),label = "**",size = 6) +
    geom_text(data = filter(plot_df, is_second),label = "*",size  = 6) +
    scale_fill_gradient(low = "white",high= "red",limits = c(minC,maxC)) +
    scale_x_discrete(position = "top") +
    scale_y_reverse(
        breaks = plot_df$y.order, 
        labels = plot_df$type,
        expand = expansion(mult = 0, add = 0),
        sec.axis = sec_axis(
            transform = ~.,
            breaks = sec_axis_data$y.center,
            labels = sec_axis_data$covarset,
            name = NULL
        )
    ) +
    coord_fixed(ratio = 1) +
    theme_bw() +
    labs(x = NULL, y = NULL, fill = "C.index") +
    theme(
        axis.text.x = element_text(angle = 30, hjust = 0, size = 12),
        axis.text.y = element_text(size = 12),
        axis.text.y.right = element_text(size = 12),
        panel.grid  = element_blank()
    )
ggsave("suppfigure1.pdf", p, width = 10, height = 10, dpi = 300)
