library(here)
library(tidyverse)
library(janitor)
library(ape)
library(ggtree)
library(viridis)
library(brms)
library(patchwork)
library(abind)
here()

australian_data <- read.csv("Australia_Data/australian_data.csv")

australia_data_feed <- australian_data %>% select(1, 5:10, 12, 33:35, 48:52, 67:69, 104:116, 151:165, 240:244, 
                                                 267:270, 294:301, 311:325)
  
  
australia_data_feed_clean <- australia_data_feed %>% filter(!if_all(7:79, is.na)) %>%
  filter(Order %in% c("Odonata", "Plecoptera", "Trichoptera", "Ephemeroptera", "Hemiptera", 
                      "Diptera", "Coleoptera",  "Megaloptera", "Neuroptera", "Lepidoptera"))




arthropod_tree <- read.tree("arthropods_genus.nwk")
arthropod_tree_extract <- sub("-.*", "", arthropod_tree$tip.label)

target_taxa <- c(unique(australia_data_feed_clean$Genus))

matching_tips <- arthropod_tree$tip.label[arthropod_tree_extract %in% target_taxa]
subtree <- drop.tip(arthropod_tree, setdiff(arthropod_tree$tip.label, matching_tips))  
base_labels <- sub("-\\d+$", "", subtree$tip.label)

dup_tips <- subtree$tip.label[duplicated(base_labels)]
subtree_cleaned <- drop.tip(subtree, dup_tips)
subtree_cleaned$tip.label <- sub("-\\d+$", "", subtree_cleaned$tip.label)

subtree_labels <- c(subtree_cleaned$tip.label)
labels_df <- tibble(label = sub("[^A-Za-z].*", "", subtree_cleaned$tip.label))
labels_df <- tibble(label = subtree_cleaned$tip.label)

tip_data_all <- labels_df %>% 
  left_join(australia_data_feed_clean, by = c("label" = "Genus"))



tip_data_1 <- tip_data_all %>% select(1,5, 12:16)
tip_data_1_condensed <- tip_data_1 %>%
  group_by(label) %>%
  summarise(
    `predator_Marchant` = sum(`predator_Marchant`, na.rm = TRUE),
    `detritivore_Marchant` = sum(`detritivore_Marchant`, na.rm = TRUE),
    `grazer_Marchant` = sum(`grazer_Marchant`, na.rm = TRUE),
    `shredder_Marchant` = sum(`shredder_Marchant`, na.rm = TRUE), 
    `filterer_Marchant` = sum(`filterer_Marchant`, na.rm = TRUE)
  ) %>%
  pivot_longer(
    cols = -label,
    names_to = "trait",
    values_to = "count"
  ) %>%
  filter(count > 0) %>%
  mutate(trait = gsub("_Marchant", "", trait)) %>% 
  mutate(trait = recode(trait,
                        predator = "PR",
                        detritivore = "DE",
                        grazer = "GR",
                        shredder = "SH",
                        filterer = "FF"))


tip_data_2 <- tip_data_all %>% select(1,5, 48:52)
tip_data_2_condensed <- tip_data_2 %>%
  group_by(label) %>%
  summarise(
    `Trop1_botwe` = sum(`Trop1_botwe`, na.rm = TRUE),
    `Trop2_botwe` = sum(`Trop2_botwe`, na.rm = TRUE),
    `Trop3_botwe` = sum(`Trop3_botwe`, na.rm = TRUE),
    `Trop4_botwe` = sum(`Trop4_botwe`, na.rm = TRUE), 
    `Trop5_botwe` = sum(`Trop5_botwe`, na.rm = TRUE)
  ) %>%
  pivot_longer(
    cols = -label,
    names_to = "trait",
    values_to = "count"
  ) %>%
  filter(count > 0) %>%
  mutate(trait = recode(trait,
                        Trop1_botwe = "CG",
                        Trop2_botwe = "CF",
                        Trop3_botwe = "HB",
                        Trop4_botwe = "PR",
                        Trop5_botwe = "SH"))



tip_data_3 <- tip_data_all %>% select(1,5, 57:64)
tip_data_3_condensed <- tip_data_3 %>%
  group_by(label) %>%
  summarise(
    `C_Maxwell` = sum(`C_Maxwell`, na.rm = TRUE),
    `P_Maxwell` = sum(`P_Maxwell`, na.rm = TRUE),
    `SH_Maxwell` = sum(`SH_Maxwell`, na.rm = TRUE),
    `C_SH_Maxwell` = sum(`C_SH_Maxwell`, na.rm = TRUE), 
    `SC_Maxwell` = sum(`SC_Maxwell`, na.rm = TRUE),
    `PA_Maxwell` = sum(`PA_Maxwell`, na.rm = TRUE),
    `C_SC_Maxwell` = sum(`C_SC_Maxwell`, na.rm = TRUE),
    `F_Maxwell` = sum(`F_Maxwell`, na.rm = TRUE)
  ) %>%
  pivot_longer(
    cols = -label,
    names_to = "trait",
    values_to = "count"
  ) %>%
  filter(count > 0) %>%
  mutate(trait = recode(trait,
                        C_Maxwell = "CG",
                        P_Maxwell = "PR",
                        SH_Maxwell = "SH",
                        C_SH_Maxwell = "CSH",
                        SC_Maxwell = "SC", 
                        PA_Maxwell = "PA", 
                        C_SC_Maxwell = "CSC",
                        F_Maxwell = "FF")) 



tip_data_4 <- tip_data_all %>% select(1,5, 65, 68, 71, 74, 77)
tip_data_4_condensed <- tip_data_4 %>%
  group_by(label) %>%
  summarise(
    `Shredder_proportion_of_feeding_genus_Chessman2017` = sum(`Shredder_proportion_of_feeding_genus_Chessman2017`, na.rm = TRUE),
    `Scraper_proportion_of_feeding_genus_Chessman2017` = sum(`Scraper_proportion_of_feeding_genus_Chessman2017`, na.rm = TRUE),
    `Predator_proportion_of_feeding_genus_Chessman2017` = sum(`Predator_proportion_of_feeding_genus_Chessman2017`, na.rm = TRUE),
    `Gatherer_proportion_of_feeding_genus_Chessman2017` = sum(`Gatherer_proportion_of_feeding_genus_Chessman2017`, na.rm = TRUE), 
    `Filterer_proportion_of_feeding_genus_Chessman2017` = sum(`Filterer_proportion_of_feeding_genus_Chessman2017`, na.rm = TRUE),
  )

tip_4_mat <- as.matrix(tip_data_4_condensed[2:6])
tip_4_mat <- tip_4_mat / rowSums(tip_4_mat)
row.names(tip_4_mat) <- tip_data_4_condensed$label
colnames(tip_4_mat) <- c("SH", "SC", "PR", "CG", "FF")
tip_4_mat <- tip_4_mat[, c("PR", "CG", "SH", "SC", "FF")]
tip_4_mat <- na.omit(tip_4_mat)


tree <- ladderize(subtree_cleaned)
tip_labels <- tree$tip.label
tree_plot <- ggtree(subtree_cleaned)


target_taxa_1 <- c(unique(tip_data_1_condensed$label))
matching_tips_1 <- arthropod_tree$tip.label[arthropod_tree_extract %in% target_taxa_1]
subtree_1 <- drop.tip(arthropod_tree, setdiff(arthropod_tree$tip.label, matching_tips_1)) 
base_labels_1 <- sub("-\\d+$", "", subtree_1$tip.label)
dup_tips_1 <- subtree_1$tip.label[duplicated(base_labels_1)]
subtree_cleaned_1 <- drop.tip(subtree_1, dup_tips_1)
subtree_cleaned_1$tip.label <- sub("-\\d+$", "", subtree_cleaned_1$tip.label)


ffg_phylo_1 <- ggtree(subtree_cleaned_1, layout = "fan") %<+% tip_data_1_condensed +
  geom_tippoint(aes(fill = trait), shape = 21, size = 4) +
  scale_fill_viridis_d(
    name = NULL,
    option = "viridis") +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 20))


target_taxa_2 <- c(unique(tip_data_2_condensed$label))
matching_tips_2 <- arthropod_tree$tip.label[arthropod_tree_extract %in% target_taxa_2]
subtree_2 <- drop.tip(arthropod_tree, setdiff(arthropod_tree$tip.label, matching_tips_2)) 
base_labels_2 <- sub("-\\d+$", "", subtree_2$tip.label)
dup_tips_2 <- subtree_2$tip.label[duplicated(base_labels_2)]
subtree_cleaned_2 <- drop.tip(subtree_2, dup_tips_2)
subtree_cleaned_2$tip.label <- sub("-\\d+$", "", subtree_cleaned_2$tip.label)


ffg_phylo_2 <- ggtree(subtree_cleaned_2, layout = "fan") %<+% tip_data_2_condensed +
  geom_tippoint(aes(fill = trait), shape = 21, size = 4) +
  scale_fill_viridis_d(
    name = NULL,
    option = "viridis") +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 20))



target_taxa_3 <- c(unique(tip_data_3_condensed$label))
matching_tips_3 <- arthropod_tree$tip.label[arthropod_tree_extract %in% target_taxa_3]
subtree_3 <- drop.tip(arthropod_tree, setdiff(arthropod_tree$tip.label, matching_tips_3)) 
base_labels_3 <- sub("-\\d+$", "", subtree_3$tip.label)
dup_tips_3 <- subtree_3$tip.label[duplicated(base_labels_3)]
subtree_cleaned_3 <- drop.tip(subtree_3, dup_tips_3)
subtree_cleaned_3$tip.label <- sub("-\\d+$", "", subtree_cleaned_3$tip.label)


ffg_phylo_3 <- ggtree(subtree_cleaned_3, layout = "fan") %<+% tip_data_3_condensed +
  geom_tippoint(aes(fill = trait), shape = 21, size = 4) +
  scale_fill_viridis_d(
    name = NULL,
    option = "viridis") +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 20))


target_taxa_4 <- c(unique(rownames(tip_4_mat)))
matching_tips_4 <- unique(arthropod_tree$tip.label[arthropod_tree_extract %in% target_taxa_4])
subtree_4 <- drop.tip(arthropod_tree, setdiff(arthropod_tree$tip.label, matching_tips_4)) 
base_labels_4 <- sub("-\\d+$", "", subtree_4$tip.label)
dup_tips_4 <- subtree_4$tip.label[duplicated(base_labels_4)]
subtree_cleaned_4 <- drop.tip(subtree_4, dup_tips_4)
subtree_cleaned_4$tip.label <- sub("-\\d+$", "", subtree_cleaned_4$tip.label)

ggtree_4 <- ggtree(subtree_cleaned_4)

ffg_phylo_4 <- gheatmap(
  ggtree_4,                
  tip_4_mat,
  offset = 0.3,  
  width = 0.5,
  colnames = FALSE,
  font.size = 8
) +
  scale_fill_viridis(option = "viridis") +
  theme(
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.key.size = unit(0.5, "cm"))





corr_matrix_1 <- vcv.phylo(subtree_cleaned_1, corr = TRUE)
corr_matrix_2 <- vcv.phylo(subtree_cleaned_2, corr = TRUE)
corr_matrix_3 <- vcv.phylo(subtree_cleaned_3, corr = TRUE)
corr_matrix_4 <- vcv.phylo(subtree_cleaned_4, corr = TRUE)

tip_data_1_condensed$label <- factor(tip_data_1_condensed$label, levels = rownames(corr_matrix_1))
tip_data_2_condensed$label <- factor(tip_data_2_condensed$label, levels = rownames(corr_matrix_2))
tip_data_3_condensed$label <- factor(tip_data_3_condensed$label, levels = rownames(corr_matrix_3))
rownames(tip_4_mat_cut) <- factor(rownames(tip_4_mat_cut), levels = rownames(corr_matrix_4))


epsilon <- 1e-6

tip_1_mat_cut_df <- as.data.frame(tip_1_mat_cut)
tip_2_mat_cut_df <- as.data.frame(tip_2_mat_cut)
tip_3_mat_cut_df <- as.data.frame(tip_3_mat_cut)


tip_4_mat_cut <- tip_4_mat[!apply(is.nan(tip_4_mat), 1, any), , drop = FALSE]
tip_4_mat_cut_df <- as.data.frame(tip_4_mat_cut) %>%
  rownames_to_column(var = "label") %>%
  mutate(across(all_of(2:6), ~ .x + epsilon)) %>%
  rowwise() %>%
  mutate(
    total = sum(c_across(2:6)),                  
    across(2:6, ~ .x / total)                    
  ) %>%
  ungroup() %>%
  select(-total) 


tip_4_mat_cut_df$diet <- with(tip_4_mat_cut_df,
                              cbind(PR, CG, SH, SC, FF))



genera_1 <- rownames(tip_1_mat)
keep_1 <- intersect(genera_1, rownames(corr_matrix))
corr_matrix_1 <- corr_matrix[keep_1, keep_1, drop = FALSE]
corr_matrix_1 <- corr_matrix_1 + diag(1e-5, nrow(corr_matrix_1))

genera_2 <- rownames(tip_2_mat)
keep_2 <- intersect(genera_2, rownames(corr_matrix))
corr_matrix_2 <- corr_matrix[keep_2, keep_2, drop = FALSE]
corr_matrix_2 <- corr_matrix_2 + diag(1e-5, nrow(corr_matrix_2))

genera_3 <- rownames(tip_3_mat)
keep_3 <- intersect(genera_3, rownames(corr_matrix))
corr_matrix_3 <- corr_matrix[keep_3, keep_3, drop = FALSE]
corr_matrix_3 <- corr_matrix_3 + diag(1e-5, nrow(corr_matrix_3))

genera_4 <- rownames(tip_4_mat)
keep_4 <- intersect(genera_4, rownames(corr_matrix))
corr_matrix_4 <- corr_matrix[keep_4, keep_4, drop = FALSE]
corr_matrix_4 <- corr_matrix_4 + diag(1e-5, nrow(corr_matrix_4))



form_1 <- bf(
  trait ~ 1 + (1 | gr(label, cov = corr_matrix_1)))

form_2 <- bf(
  trait ~ 1 + (1 | gr(label, cov = corr_matrix_2)))

form_3 <- bf(
  trait ~ 1 + (1 | gr(label, cov = corr_matrix_3)))

form_4 <- bf(
  diet ~ 1 + (1 | gr(label, cov = corr_matrix_4)))


prior_1 <- get_prior(
  form_1,
  data = tip_data_1_condensed,
  family = categorical (link = "logit"), 
  data2 = list(corr_matrix_1 = corr_matrix_1)
)

prior_2 <- get_prior(
  form_2,
  data = tip_data_2_condensed,
  family = categorical (link = "logit"), 
  data2 = list(corr_matrix_2 = corr_matrix_2)
)

prior_3 <- get_prior(
  form_3,
  data = tip_data_3_condensed,
  family = categorical (link = "logit"), 
  data2 = list(corr_matrix_3 = corr_matrix_3)
)

prior_4 <- get_prior(
  form_4,
  data = tip_4_mat_cut_df,
  family = dirichlet(),
  data2 = list(corr_matrix_4 = corr_matrix_4)
)


fit_1 <- brm(form_1,
           data = tip_data_1_condensed,
           family = categorical (link = "logit"), 
           data2 = list(corr_matrix_1 = corr_matrix_1),
           prior = prior_1,
           cores = 4, chains = 4,  iter = 4000,  
           control = list(
             adapt_delta = 0.90,
             max_treedepth = 10)) 


fit_2 <- brm(form_2,
             data = tip_data_2_condensed,
             family = categorical (link = "logit"), 
             data2 = list(corr_matrix_2 = corr_matrix_2),
             prior = prior_2,
             cores = 4, chains = 4,  iter = 4000,  
             control = list(
               adapt_delta = 0.80,
               max_treedepth = 10))


fit_3 <- brm(form_3,
             data = tip_data_3_condensed,
             family = categorical (link = "logit"), 
             data2 = list(corr_matrix_3 = corr_matrix_3),
             prior = prior_3,
             cores = 4, chains = 4,  iter = 4000,  
             control = list(
               adapt_delta = 0.80,
               max_treedepth = 10))


fit_4 <- brm(form_4,
             data = tip_4_mat_cut_df,
             family = dirichlet(),
             data2 = list(corr_matrix_4 = corr_matrix_4),
             prior = prior_4,
             cores = 4, chains = 4,  iter = 4000,  
             control = list(
               adapt_delta = 0.90,
               max_treedepth = 10))
 

summary(fit_1)
summary(fit_2)
summary(fit_3)
summary(fit_4)

posterior_1 <- as_draws_df(fit_1)
phylo_cols_1 <- grepl("^sd_label__mu.*_Intercept", colnames(posterior_1))
sd_phylo_mat_1 <- posterior_1[, phylo_cols_1]
var_phylo_1 <- sd_phylo_mat_1^2
var_phylo_total_1 <- rowSums(var_phylo_1)
sigma_resid_1 <- (pi^2) / 3
sigma_resid_total_1 <- ncol(var_phylo_1) * sigma_resid_1
VPC_post_total_1 <- var_phylo_total_1 / (var_phylo_total_1 + sigma_resid_total_1)
VPC_summar_1y <- tibble(
  mean_VPC = mean(VPC_post_total_1),
  lower_CI = quantile(VPC_post_total_1, 0.025),
  upper_CI = quantile(VPC_post_total_1, 0.975)
)
print(VPC_summary_1)


posterior_2 <- as_draws_df(fit_2)
phylo_cols_2 <- grepl("^sd_label__mu.*_Intercept", colnames(posterior_2))
sd_phylo_mat_2 <- posterior_2[, phylo_cols_2]
var_phylo_2 <- sd_phylo_mat_2^2
var_phylo_total_2 <- rowSums(var_phylo_2)
sigma_resid_2 <- (pi^2) / 3
sigma_resid_total_2 <- ncol(var_phylo_2) * sigma_resid_2
VPC_post_total_2 <- var_phylo_total_2 / (var_phylo_total_2 + sigma_resid_total_2)
VPC_summary_2 <- tibble(
  mean_VPC = mean(VPC_post_total_2),
  lower_CI = quantile(VPC_post_total_2, 0.025),
  upper_CI = quantile(VPC_post_total_2, 0.975)
)
print(VPC_summary_2)



posterior_3 <- as_draws_df(fit_3)
phylo_cols_3 <- grepl("^sd_label__mu.*_Intercept", colnames(posterior_3))
sd_phylo_mat_3 <- posterior_3[, phylo_cols_3]
var_phylo_3 <- sd_phylo_mat_3^2
var_phylo_total_3 <- rowSums(var_phylo_3)
sigma_resid_3 <- (pi^2) / 3
sigma_resid_total_3 <- ncol(var_phylo_3) * sigma_resid_3
VPC_post_total_3 <- var_phylo_total_3 / (var_phylo_total_3 + sigma_resid_total_3)
VPC_summary_3 <- tibble(
  mean_VPC = mean(VPC_post_total_3),
  lower_CI = quantile(VPC_post_total_3, 0.025),
  upper_CI = quantile(VPC_post_total_3, 0.975)
)
print(VPC_summary_3)


posterior_4 <- as_draws_df(fit_4)
phylo_cols_4 <- grepl("^sd_label__", colnames(posterior_4))
sd_phylo_mat_4 <- posterior_4[, phylo_cols_4]
var_phylo_mat_4 <- sd_phylo_mat_4^2
phylo_var_4 <- rowSums(var_phylo_mat_4)
phi_samples_4 <- posterior_4$phi
mu_mat_4 <- tip_4_mat_cut_df %>% select(2:6) %>% as.matrix()
mu_4 <- colMeans(mu_mat_4)
resid_var_4 <- sapply(phi_samples_4, function(phi_s) {
  sum(mu_4 * (1 - mu_4) / (1 + phi_s))
})
VPC_post_4 <- phylo_var_4 / (phylo_var_4 + resid_var_4)
VPC_summary_4 <- tibble(
  mean_VPC = mean(VPC_post_4),
  lower_CI = quantile(VPC_post_4, 0.025),
  upper_CI = quantile(VPC_post_4, 0.975)
)
print(VPC_summary_4)




categorical_vpc <- function(fit) {
  
  # --- STEP 1: posterior predictions ---
  y_with_phy <- posterior_epred(fit)
  y_without_phy <- posterior_epred(fit, re_formula = NA)
  
  # dims: draws × observations × categories
  n_cat <- dim(y_with_phy)[3]
  cat_names <- dimnames(y_with_phy)[[3]]
  
  # --- STEP 2: isolate phylogenetic effects ---
  phylo_effects <- sapply(1:n_cat, function(k) {
    rowMeans(y_with_phy[,,k]) - rowMeans(y_without_phy[,,k])
  })
  colnames(phylo_effects) <- cat_names
  
  # --- STEP 3: variance partition ---
  vpc_draws_matrix <- t(apply(phylo_effects, 1, function(x) {
    
    var_trait <- x^2
    total_var <- sum(var_trait)
    
    if (total_var == 0) return(rep(0, length(x)))
    
    pmin(pmax(var_trait / total_var * 100, 0), 100)
  }))
  
  # --- STEP 4: summarise ---
  vpc_summary <- data.frame(
    trait = colnames(vpc_draws_matrix),
    mean = apply(vpc_draws_matrix, 2, mean),
    sd   = apply(vpc_draws_matrix, 2, sd),
    lci  = apply(vpc_draws_matrix, 2, quantile, probs = 0.025),
    uci  = apply(vpc_draws_matrix, 2, quantile, probs = 0.975)
  )
  
  return(vpc_summary)
}
dirichlet_vpc <- function(fit) {
  
  # --- STEP 1: posterior predictions ---
  y_with_phy <- posterior_epred(fit)
  y_without_phy <- posterior_epred(fit, re_formula = NA)
  
  # dims: draws × observations × components
  n_comp  <- dim(y_with_phy)[3]
  comp_names <- dimnames(y_with_phy)[[3]]
  
  # --- STEP 2: isolate phylogenetic effects ---
  phylo_effects <- sapply(1:n_comp, function(k) {
    rowMeans(y_with_phy[,,k]) - rowMeans(y_without_phy[,,k])
  })
  colnames(phylo_effects) <- comp_names
  
  # --- STEP 3: variance partition ---
  vpc_draws_matrix <- t(apply(phylo_effects, 1, function(x) {
    
    var_trait <- x^2
    total_var <- sum(var_trait)
    
    if (total_var == 0) return(rep(0, length(x)))
    
    pmin(pmax(var_trait / total_var * 100, 0), 100)
  }))
  
  # --- STEP 4: summarise ---
  vpc_summary <- data.frame(
    trait = colnames(vpc_draws_matrix),
    mean = apply(vpc_draws_matrix, 2, mean),
    sd   = apply(vpc_draws_matrix, 2, sd),
    lci  = apply(vpc_draws_matrix, 2, quantile, probs = 0.025),
    uci  = apply(vpc_draws_matrix, 2, quantile, probs = 0.975)
  )
  
  return(vpc_summary)
}

vpc_multivariate_1 <- categorical_vpc(fit_1)
vpc_multivariate_2 <- categorical_vpc(fit_2)
vpc_multivariate_3 <- categorical_vpc(fit_3)
vpc_multivariate_4 <- dirichlet_vpc(fit_4)


png(width = 500, height = 500)
ggplot(vpc_multivariate_1, aes(x = trait, y = mean)) +
  geom_point(aes(color = trait == "PR"), size = 10) +      
  geom_errorbar(
    aes(ymin = lci, ymax = uci, color = trait == "PR"),
    width = 0.4,
    linewidth = 2   # thicker error bars
  ) + 
  scale_color_manual(
    values = c("FALSE" = "steelblue", "TRUE" = "red"),
    guide = "none"
  ) +                       
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 25, face = "bold"),  
    axis.text.y = element_text(size = 30, face = "bold"),                          
    legend.position = "bottom",
    legend.text = element_text(size = 30, face = "bold"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()                          
  ) +
  ylim(0, 100)
dev.off()

vpc_2 <- ggplot(vpc_multivariate_2, aes(x = trait, y = mean)) +
  geom_point(aes(color = trait == "PR"), size = 10) +      
  geom_errorbar(
    aes(ymin = lci, ymax = uci, color = trait == "PR"),
    width = 0.4,
    linewidth = 2   # thicker error bars
  ) + 
  scale_color_manual(
    values = c("FALSE" = "steelblue", "TRUE" = "red"),
    guide = "none"
  ) +                       
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 30, face = "bold"),  
    axis.text.y = element_text(size = 30, face = "bold"),                          
    legend.position = "bottom",
    legend.text = element_text(size = 30, face = "bold"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()                          
  ) +
  ylim(0, 100)


vpc_3 <- ggplot(vpc_multivariate_3, aes(x = trait, y = mean)) +
  geom_point(aes(color = trait == "PR"), size = 10) +      
  geom_errorbar(
    aes(ymin = lci, ymax = uci, color = trait == "PR"),
    width = 0.4,
    linewidth = 2   # thicker error bars
  ) + 
  scale_color_manual(
    values = c("FALSE" = "steelblue", "TRUE" = "red"),
    guide = "none"
  ) +                       
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 30, face = "bold"),  
    axis.text.y = element_text(size = 30, face = "bold"),                          
    legend.position = "bottom",
    legend.text = element_text(size = 30, face = "bold"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()                          
  ) +
  ylim(0, 100)

vpc_4 <- ggplot(vpc_multivariate_4, aes(x = trait, y = mean)) +
  geom_point(aes(color = trait == "PR"), size = 5) +      
  geom_errorbar(
    aes(ymin = lci, ymax = uci, color = trait == "PR"),
    width = 0.4,
    linewidth = 2   
  ) + 
  scale_color_manual(
    values = c("FALSE" = "steelblue", "TRUE" = "red"),
    guide = "none"
  ) +                       
  theme_minimal() +
  scale_x_discrete(labels =  
                     c("muCG" = "CG",
                       "muFF" = "FF",
                       "muSC" = "SC",
                       "muSH" = "SH",
                       "ref_trait" = "PR")) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 15, face = "bold"),  
    axis.text.y = element_text(size = 15, face = "bold"),                          
    legend.position = "bottom",
    legend.text = element_text(size = 30, face = "bold"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()                          
  ) +
  ylim(0, 100)


vpc_1 + vpc_2 + vpc_3 + vpc_4
