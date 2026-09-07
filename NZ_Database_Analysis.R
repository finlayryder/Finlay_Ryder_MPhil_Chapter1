library(here)
library(tidyverse)
library(janitor)
library(ape)
library(ggtree)
library(viridis)
library(brms)


all_NZ_data <- read.csv("NZ_data_raw.csv")

NZ_ffg_data <- all_NZ_data %>%
  select(1:4, 48:53) %>%
  row_to_names(row_number = 3) %>%
  filter(Order %in% c("Coleoptera", "Diptera", "Ephemeroptera", 
                      "Hemiptera", "Megaloptera", "Lepidoptera", 
                      "Neuroptera", "Trichoptera", "Plecoptera", 
                      "Trichoptera")) %>%
  filter(Family != `Genus or Higher`)

NZ_ffg_data <- NZ_ffg_data %>%
  filter(!`Genus or Higher` %in% c("Stictocladius", "Paucispinigera"))

NZ_ffg_data <- NZ_ffg_data %>%
  mutate(across(
    c(shredders, scrapers, `deposit-feeders`, `filter-feeders`, `algal piercer`, predator),
    ~ as.numeric(.)
  ))

NZ_ffg_data <- NZ_ffg_data %>% rename(
  SH = shredders, 
  SC = scrapers, 
  DF = `deposit-feeders`, 
  FF = `filter-feeders`, 
  PR = predator, 
  AP = `algal piercer`, 
  Genus = `Genus or Higher`) %>%
  filter(if_any(everything(), ~ !is.na(.)))




arthropod_tree <- read.tree("arthropods_genus.nwk")
arthropod_tree_extract <- sub("-.*", "", arthropod_tree$tip.label)
target_taxa <- unique(NZ_ffg_data$Genus)

matching_tips <- arthropod_tree$tip.label[arthropod_tree_extract %in% target_taxa]
subtree <- drop.tip(arthropod_tree, setdiff(arthropod_tree$tip.label, matching_tips)) 
base_labels <- sub("-\\d+$", "", subtree$tip.label)
dup_tips <- subtree$tip.label[duplicated(base_labels)]
subtree_cleaned <- drop.tip(subtree, dup_tips)
subtree_cleaned$tip.label <- sub("-\\d+$", "", subtree_cleaned$tip.label)

subtree_labels <- c(subtree_cleaned$tip.label)
labels_df <- tibble(label = sub("[^A-Za-z].*", "", subtree_cleaned$tip.label))
labels_df <- tibble(label = subtree_cleaned$tip.label)


tip_data <- NZ_ffg_data %>% filter(NZ_ffg_data$Genus %in% labels_df$label)

tip_data_condensed <- tip_data %>%
  group_by(Genus) %>%
  summarise(
    SH = sum(SH, na.rm = TRUE),
    SC = sum(SC, na.rm = TRUE),
    DF = sum(DF, na.rm = TRUE),
    FF = sum(FF, na.rm = TRUE), 
    AP = sum(`AP`, na.rm = TRUE),
    PR = sum(PR, na.rm = TRUE)
  )

colSums(tip_data_condensed[sapply(tip_data_condensed, is.numeric)])


ffg_matrix_norm <- as.matrix(tip_data_condensed[,2:7])
ffg_matrix_norm <- ffg_matrix_norm / rowSums(ffg_matrix_norm)
rownames(ffg_matrix_norm) <- tip_data_condensed$Genus




tree_plot <- ggtree(subtree_cleaned)


png("NZ_PhyloHeat.png", width = 1400, height = 2000)
ggtree::gheatmap(
  tree_plot,
  ffg_matrix_norm,
  offset = 100,    
  width = 0.5,
  colnames_angle = 0, 
  colnames_position = "top",
  colnames_offset_y = 0.5,
  font.size = 3) + 
  geom_tiplab(align = TRUE, size = 3, offset = 0.01) +
  scale_fill_viridis(option = "viridis") +
  theme(    legend.text = element_text(size = 12), 
            legend.title = element_text(size = 14),  
            legend.key.size = unit(1.5, "cm"))
dev.off()

corr_matrix <- vcv.phylo(subtree_cleaned, corr = TRUE)
corr_matrix_pd <- corr_matrix + diag(1e-5, nrow(corr_matrix))


epsilon <- 1e-4
cols <- c("SH","SC","DF","FF","AP","PR")
tip_data_condensed[cols] <- tip_data_condensed[cols] + epsilon
tip_data_condensed[cols] <- tip_data_condensed[cols] /
  rowSums(tip_data_condensed[cols])

tip_data_condensed$diet <- with(tip_data_condensed,
                                cbind(SC, SH, DF, FF, AP, PR))
tip_data_condensed <- tip_data_condensed %>%
  select(-(2:7))

form <- bf(
  diet ~ 1 + (1 | gr(Genus, cov = corr_matrix_pd))
)

prior <- get_prior(
  form,
  data = tip_data_condensed,
  family = dirichlet(),
  data2 = list(corr_matrix_pd = corr_matrix_pd)
)

fit <- brm(form,
           data = tip_data_condensed,
           family = dirichlet(),
           data2 = list(corr_matrix_pd = corr_matrix_pd),
           prior = prior,
           cores = 4, chains = 4,  iter = 4000,  
           control = list(
             adapt_delta = 0.9,
             max_treedepth = 10)) 

summary(fit)


posterior <- as_draws_df(fit)
phylo_cols <- grepl("^sd_Genus__", colnames(posterior))
sd_phylo_mat <- posterior[, phylo_cols]
var_phylo_mat <- sd_phylo_mat^2
phylo_var <- rowSums(var_phylo_mat)
phi_samples <- posterior$phi

mu_mat <- tip_data_condensed %>% select(2) %>% as.matrix()
mu <- colMeans(mu_mat)

resid_var <- sapply(phi_samples, function(phi_s) {
  sum(mu * (1 - mu) / (1 + phi_s))
})

VPC_post <- phylo_var / (phylo_var + resid_var)

VPC_summary <- tibble(
  mean_VPC = mean(VPC_post),
  lower_CI = quantile(VPC_post, 0.025),
  upper_CI = quantile(VPC_post, 0.975)
)

print(VPC_summary)





dirichlet_vpc <- function(fit) {
  
  ##Posterior predictions
  y_with_phy <- posterior_epred(fit)
  y_without_phy <- posterior_epred(fit, re_formula = NA)
  
  ##Dimensions: draws × observations × categories
  n_comp  <- dim(y_with_phy)[3]
  comp_names <- dimnames(y_with_phy)[[3]]
  
  ##Isolating phylogenetic effects only 
  phylo_effects <- sapply(1:n_comp, function(k) {
    rowMeans(y_with_phy[,,k]) - rowMeans(y_without_phy[,,k])
  })
  colnames(phylo_effects) <- comp_names
  
  ##Partitioning variance by that explained by phylo
  vpc_draws_matrix <- t(apply(phylo_effects, 1, function(x) {
    
    var_trait <- x^2
    total_var <- sum(var_trait)
    
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
vpc_multivariate <- dirichlet_vpc(fit)
print(vpc_multivariate)

ggplot(vpc_multivariate, aes(x = trait, y = mean)) +
  geom_point(aes(color = trait == "SC"), size = 5) +      
  geom_errorbar(aes(ymin = lci, ymax = uci, color = trait == "SC"), width = 0.4, 
                linewidth = 2) + 
  scale_color_manual(values = c("FALSE" = "steelblue", "TRUE" = "red"), guide = "none") +
  scale_x_discrete(labels =  
                     c("muSH" = "SH",
                       "muDF" = "DF",
                       "muFF" = "FF",
                       "muAP" = "AP",
                       "muPR" = "PR",
                       "ref_trait" = "SC")) +                       
  ylab("VPC (%)") +
  xlab("Trait") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 15, face = "bold"),  
    axis.text.y = element_text(size = 15, face = "bold"),                          
    legend.position = "bottom",
    legend.text = element_text(size = 15, face = "bold"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()                          
  ) +
  ylim(c(0,100))



