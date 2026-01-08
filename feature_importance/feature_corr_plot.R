library(dplyr)
library(data.table)
library(tidyr)
library(ggplot2)
library(tibble)
library(patchwork)
library(pheatmap)
library(RColorBrewer)

for (t in c("MACE", "Heart_Failure", "CHD", "Atrial_Fibrillation", "PAD", "T2_Diabetes", "Prostate_Cancer", 
           "Skin_Cancer", "Breast_Cancer", "Dementia", "Cataracts", "Glaucoma", "COPD", 
           "Asthma", "Renal_Disease", "Liver_Disease", "Fractures")) {
    feat = fread(paste0("mp_", t, "_PANEL_weighted_abs.txt")) %>% 
            select(-V1) %>% mutate(!!sym(paste0("rank_",t)) := row_number(), !!sym(paste0("norm_",t)) := scale(Final_Score, center = FALSE, scale = TRUE))
    feat = feat %>% select(Feature, !!(sym(paste0("norm_",t)))) %>% rename(!!sym(t) := paste0("norm_", t))
    if (t == "MACE") {
        all_feat = feat
    } else {
        all_feat = all_feat %>% left_join(feat, by = "Feature")
    }
}
cormatrix = cor(all_feat[,2:ncol(all_feat)], method = "pearson",use = "pairwise.complete.obs")

# draw correlation plot
cor.plotdf = cormatrix %>% as.data.frame() %>% rownames_to_column("feat1") %>%
    pivot_longer(-feat1, names_to = "feat2", values_to = "correlation") %>%
    filter(!is.na(correlation)) %>% 
    mutate(
        feat1 = factor(feat1, levels = c("MACE", "Heart_Failure", "CHD", "Atrial_Fibrillation", "PAD","T2_Diabetes", 
                                          "Prostate_Cancer", "Skin_Cancer", "Breast_Cancer", "Dementia", 
                                          "Cataracts", "Glaucoma", "COPD", "Asthma", 
                                          "Renal_Disease", "Liver_Disease", "Fractures")),
        feat2 = factor(feat2, levels = c("MACE", "Heart_Failure", "CHD", "Atrial_Fibrillation", "PAD", "T2_Diabetes", 
                                          "Prostate_Cancer", "Skin_Cancer", "Breast_Cancer", "Dementia", 
                                          "Cataracts", "Glaucoma", "COPD", "Asthma", 
                                          "Renal_Disease", "Liver_Disease", "Fractures"))
    )
p = ggplot(cor.plotdf, aes(x = feat1, y = feat2, fill = correlation)) +
    geom_tile(color = "lightgrey", linewidth = 0.3) +
    geom_text(aes(label = round(correlation, 3)), size = 2) +
    scale_fill_gradient2(low = "blue", mid = "snow", high = "red", midpoint = 0, name = "Features\nImportance\nCorrelation", limits = c(-1, 1), breaks = c(-1, 0, 1)) +
    theme_bw() +
    labs(x = "", y = "") +
    scale_y_discrete(limits = rev) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          aspect.ratio = 1) 
    # guides(fill = guide_colorbar(title = "Protein\nImportance\nCorrelation"))
ggsave("all.feat.correlation.plot.pdf", p, width = 8, height = 8)
