library(here)
library(tidyverse)
library(janitor)
library(ape)
library(ggtree)
library(brms)

here()

all_sa_data <- read.csv("SA_Database.csv")

sa_data_cleaned <- all_sa_data %>%
  slice(-(1:11)) %>%
  select(1:6, 20) %>%
  clean_names() %>%
  mutate(
    across(everything(), ~ na_if(.x, "")),
    ffg = functional_feeding_group_ffg_schael_2006 %>%
      str_replace("\\s*\\(.*$", "") %>%  
      str_squish()) %>%
  select(-(7)) %>%
  filter(!is.na(ffg) & !is.na(genus)) %>%
  filter(order %in% c(
    "Coleoptera", "Diptera", "Odonata", "Plecoptera",
    "Trichoptera", "Hemiptera", "Megaloptera", "Ephemeroptera"
  )) %>%
  filter(!ffg %in% c("grazer", "Varied"))


arthropod_tree <- read.tree("arthropods_genus.nwk")
arthropod_tree_extract <- sub("-.*", "", arthropod_tree$tip.label)

target_taxa <- c(unique(sa_data_cleaned$genus))

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
  left_join(sa_data_cleaned, by = c("label" = "genus"))


ffg_phylo <- ggtree(subtree_cleaned, layout = "fan") %<+% tip_data +
  geom_tippoint(aes(fill = ffg), shape = 21, size = 4) +
  scale_fill_viridis_d(
    name = NULL,  # removes the legend title
    option = "viridis",
    labels = c(
      `Predator 1` = "PR1",
      `Predator 2` = "PR2",
      `Grazer 1` = "GR1",
      `Scraper 1` = "SC1",
      `Filter feeder` = "FF",
      `Scraper 2` = "SCR2",
      `Deposit feeder 2` = "DF2",
      `Shredder` = "SH",
      `Deposit feeder 1` = "DF1",
      `Grazer 2` = "GR2"
    )
  ) +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 20)  # just the text size
  )


corr_matrix <- vcv.phylo(subtree_cleaned, corr = TRUE)
corr_matrix_pd <- corr_matrix + diag(1e-6, nrow(corr_matrix))
tip_data$label <- factor(tip_data$label, levels = rownames(corr_matrix_pd))

tip_data %>% count(ffg)

tip_data$ffg <- factor(tip_data$ffg,
                       levels = c("Predator 2", "Deposit feeder 1", "Grazer 1", "Filter feeder", "Predator 1",
                                  "Shredder", "Scraper 1", "Scraper 2", "Grazer 2", "Deposit feeder 2"))

prior <- get_prior(bf(ffg ~ 1 + (1|gr(label, cov = corr_matrix_pd))), 
                   data = tip_data, 
                   family = categorical (link = "logit"), 
                   data2 = list(corr_matrix_pd = corr_matrix_pd))

fit <- brm(bf(ffg ~ 1+(1|gr(label, cov = corr_matrix_pd))),
               data = tip_data, 
               family = categorical(link = "logit"), 
               data2 = list(corr_matrix_pd = corr_matrix_pd), 
               cores = 4, chains = 4, iter = 4000,
               control = list(adapt_delta = 0.9,
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




categorical_vpc <- function(fit, ref_name = NULL) {
  
  ##Posterior predictions
  y_with_phy <- posterior_epred(fit)
  y_without_phy <- posterior_epred(fit, re_formula = NA)
  
  ##Dimensions: draws × observations × categories
  n_draws <- dim(y_with_phy)[1]
  n_cat   <- dim(y_with_phy)[3]
  cat_names <- dimnames(y_with_phy)[[3]]

  modelled_names <- cat_names[-n_cat]
  ref_default <- cat_names[n_cat]
  
  if (is.null(ref_name)) {
    ref_name <- ref_default
  }
  
  ##Isolating phylogenetic effects only 
  phylo_effects <- sapply(1:(n_cat - 1), function(k) {
    rowMeans(y_with_phy[,,k]) - rowMeans(y_without_phy[,,k])
  })
  
  colnames(phylo_effects) <- modelled_names
  
  phylo_ref <- -rowSums(phylo_effects)
  
  phylo_effects_full <- cbind(phylo_effects, ref = phylo_ref)
  colnames(phylo_effects_full)[ncol(phylo_effects_full)] <- ref_name
  
  ##Partitioning variance by that explained by phylo
  vpc_draws_matrix <- t(apply(phylo_effects_full, 1, function(x) {
    
    var_trait <- x^2
    total_var <- sum(x^2)
    
    if (total_var == 0) return(rep(0, length(x)))
    
    pmin(pmax(var_trait / total_var * 100, 0), 100)
  }))
  
  ##Summarise outputs
  vpc_summary <- data.frame(
    trait = colnames(vpc_draws_matrix),
    mean = apply(vpc_draws_matrix, 2, mean),
    sd   = apply(vpc_draws_matrix, 2, sd),
    lci  = apply(vpc_draws_matrix, 2, quantile, probs = 0.025),
    uci  = apply(vpc_draws_matrix, 2, quantile, probs = 0.975)
  )
  
  return(vpc_summary)
}


vpc_multivariate <- categorical_vpc(fit)
print(vpc_multivariate)



ggplot(vpc_multivariate, aes(x = trait, y = mean)) +
  geom_point(aes(color = trait == "Predator 2"), size = 10) +      
  geom_errorbar(
    aes(ymin = lci, ymax = uci, color = trait == "Predator 2"),
    width = 0.4,
    linewidth = 2   # thicker error bars
  ) + 
  scale_color_manual(
    values = c("FALSE" = "steelblue", "TRUE" = "red"),
    guide = "none"
  ) +
  scale_x_discrete(
    labels = c(
      "Deposit feeder 1" = "DF1",
      "Grazer 1" = "GR1",
      "Filter feeder" = "FF",
      "Predator 1" = "PR1",
      "Shredder" = "SH",
      "Scraper 1" = "SC1",
      "Scraper 2" = "SC2",
      "Grazer 2" = "GR",
      "Deposit feeder 2" = "DF2",
      "Predator 2" = "PR2"
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

