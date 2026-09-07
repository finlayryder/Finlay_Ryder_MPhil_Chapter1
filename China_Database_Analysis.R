library(here)
library(tidyverse)
library(janitor)
library(ape)
library(ggtree)
library(viridis)
library(brms)
here()

china_data_all <- read.csv("China_Database/China_Data.csv") 

china_feed <- china_data_all %>% select(1:3, 16) %>%
  filter(!is.na(Feed1), Feed1 != "#N/A")



arthropod_tree <- read.tree("TimeTree_Phylogenies/arthropods_genus.nwk")
arthropod_tree_extract <- sub("-.*", "", arthropod_tree$tip.label)

target_taxa <- c(unique(china_feed$Genus))


matching_tips <- arthropod_tree$tip.label[arthropod_tree_extract %in% target_taxa]
subtree <- drop.tip(arthropod_tree, setdiff(arthropod_tree$tip.label, matching_tips))  
base_labels <- sub("-\\d+$", "", subtree$tip.label)

dup_tips <- subtree$tip.label[duplicated(base_labels)]
subtree_cleaned <- drop.tip(subtree, dup_tips)
subtree_cleaned$tip.label <- sub("-\\d+$", "", subtree_cleaned$tip.label)

subtree_labels <- c(subtree_cleaned$tip.label)
labels_df <- tibble(label = sub("[^A-Za-z].*", "", subtree_cleaned$tip.label))
labels_df <- tibble(label = subtree_cleaned$tip.label)


tip_data <- labels_df %>% 
  left_join(china_feed, by = c("label" = "Genus"))

ffg_long <- tip_data %>%
  mutate(across(everything(), ~ as.numeric(trimws(.)))) %>%
  tibble::rownames_to_column("Taxon") %>%
  pivot_longer(-Taxon, values_to = "Trait") %>%
  filter(!is.na(Trait))

ffg_long <- ffg_long %>%
  group_by(Taxon) %>%
  mutate(Weight = 1 / n())

ffg_wide <- ffg_long %>%
  select(-name) %>%         
  distinct() %>%             
  pivot_wider(
    names_from = Trait,
    values_from = Weight
  )


ffg_df <- as.data.frame(ffg_wide[,2:6])
rownames(ffg_df) <- tip_data$label

ffg_df <- ffg_df %>%
  dplyr::rename(
    PR = `1`,
    CG = `2`,
    CF = `3`,
    HB = `4`,
    SH = `5`,
  )

ffg_matrix <- as.matrix(ffg_df)


tree_plot <- ggtree(subtree_cleaned)


ffg_phylo <- ggtree(subtree_cleaned, layout = "fan") %<+% tip_data +
  geom_tippoint(aes(fill = Feed1), shape = 21, size = 5) +
  scale_fill_viridis_d(
    name = NULL,
    option = "viridis",
    labels = c(
      `1` = "PR",
      `2` = "CG",
      `3` = "CF",
      `4` = "HB",
      `5` = "SH"
    )
  ) +
  theme(legend.position = "bottom") + 
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 20)    # makes the legend labels bigger
  )


corr_matrix <- vcv.phylo(subtree_cleaned, corr = TRUE)
tip_data$label <- factor(tip_data$label, levels = rownames(corr_matrix))

#ffg_df_zeros <- ffg_df %>% mutate(across(everything(), ~replace_na(.x, 0)))

#dat_cols <- c("SH","PR","HB","CG","CF")
#epsilon <- 1e-6
##ffg_df_zeros[, dat_cols] <- ffg_df_zeros[, dat_cols] + epsilon
##ffg_df_zeros[, dat_cols] <- ffg_df_zeros[, dat_cols] / rowSums(ffg_df_zeros[, dat_cols])


##tip_data$diet <- with(ffg_df_zeros,
                      #cbind(PR, SH, HB, CG, CF))

##colSums(tip_data[, c("diet")], na.rm = TRUE)


tip_data$Feed1 <- as.character(tip_data$Feed1)  # ensure it's character first
tip_data$Feed1[tip_data$Feed1 == "1"] <- "PR"
tip_data$Feed1[tip_data$Feed1 == "2"] <- "CG"
tip_data$Feed1[tip_data$Feed1 == "3"] <- "CF"
tip_data$Feed1[tip_data$Feed1 == "4"] <- "HB"
tip_data$Feed1[tip_data$Feed1 == "5"] <- "SH"


prior <- get_prior(bf(Feed1 ~ 1 + (1|gr(label, cov = corr_matrix))), 
                   data = tip_data, 
                   family = categorical (link = "logit"), 
                   data2 = list(corr_matrix = corr_matrix))

fit <- brm(bf(Feed1 ~ 1+(1|gr(label, cov = corr_matrix))),
           data = tip_data, 
           family = categorical(link = "logit"), 
           data2 = list(corr_matrix = corr_matrix), 
           cores = 4, chains = 4, iter = 4000,
           control = list(adapt_delta = 0.95,
                          max_treedepth = 10),
           prior = prior)


summary(fit)

posterior <- as_draws_df(fit)
phylo_cols <- grepl("^sd_label__mu.*_Intercept", colnames(posterior))
sd_phylo_mat <- posterior[, phylo_cols]
var_phylo <- sd_phylo_mat^2
var_phylo_total <- rowSums(var_phylo)
sigma_resid <- (pi^2) / 3
sigma_resid_total <- ncol(var_phylo) * sigma_resid
VPC_post_total <- var_phylo_total / (var_phylo_total + sigma_resid_total)
VPC_summary <- tibble(
  mean_VPC = mean(VPC_post_total),
  lower_CI = quantile(VPC_post_total, 0.025),
  upper_CI = quantile(VPC_post_total, 0.975)
)
print(VPC_summary)



categorical_vpc_absolute <- function(fit, ref_name = NULL) {
  
  # --- STEP 0: predictions ---
  y_with_phy <- posterior_epred(fit)          
  y_without_phy <- posterior_epred(fit, re_formula = NA)
  
  n_draws <- dim(y_with_phy)[1]
  n_obs   <- dim(y_with_phy)[2]
  n_cat   <- dim(y_with_phy)[3]
  
  # --- category names ---
  cat_names <- dimnames(y_with_phy)[[3]]
  modelled_names <- cat_names[-n_cat]
  ref_default <- cat_names[n_cat]
  if (is.null(ref_name)) ref_name <- ref_default
  
  # --- STEP 1: phylogenetic effects ---
  phylo_effects <- sapply(1:(n_cat - 1), function(k) {
    y_with_phy[, , k] - y_without_phy[, , k]
  })
  phylo_effects <- array(phylo_effects, dim = c(n_draws, n_obs, n_cat - 1))
  dimnames(phylo_effects)[[3]] <- modelled_names
  
  # --- STEP 2: reference category ---
  phylo_ref <- -apply(phylo_effects, c(1,2), sum)
  phylo_effects_full <- abind::abind(phylo_effects, phylo_ref, along = 3)
  dimnames(phylo_effects_full)[[3]][n_cat] <- ref_name
  
  # --- STEP 3: absolute VPC per category ---
  vpc_draws_matrix <- matrix(NA, nrow = n_draws, ncol = n_cat)
  colnames(vpc_draws_matrix) <- dimnames(phylo_effects_full)[[3]]
  
  for (i in 1:n_draws) {
    phylo_i <- phylo_effects_full[i,,]           # n_obs x n_cat
    y_i     <- y_with_phy[i,,]                    # n_obs x n_cat
    
    total_var <- sum(apply(y_i, 2, var))         # sum of per-category variances across observations
    phylo_var <- apply(phylo_i, 2, var)          # variance of phylo effect per category across observations
    
    vpc_draws_matrix[i, ] <- if (total_var == 0) rep(0, n_cat) else phylo_var / total_var * 100
  }
  
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
vpc_multivariate <- categorical_vpc_absolute(fit)
print(vpc_multivariate)

ggplot(vpc_multivariate, aes(x = trait, y = mean)) +
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
  scale_x_discrete(
    labels = c(
      "muSH" = "SH",
      "muHB" = "HB",
      "muCG" = "CG",
      "muCF" = "CF",
      "PR" = "PR"
    )
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















