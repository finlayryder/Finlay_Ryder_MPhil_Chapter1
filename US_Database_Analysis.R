library(tidyverse)
library(here)
library(ape)
library(ggtree)
library(viridis)
library(brms)
library(ggplot2)

US_data_TSN <- read.csv("US_TSN_Data.csv")
US_trait_data <- read.csv("US_Trait_Data.csv")


US_data_TSN <- US_data_TSN %>% filter(Taxonomic_resolution == "Genus") %>%
  select(3,4) %>% distinct() %>%
  rename(Genus = Submitted_name_trait)

feed_trait_data <- US_trait_data %>%
  filter(Trait_group == "Feed_prim_abbrev") %>%
  left_join(US_data_TSN, by = "Genus", relationship = "many-to-many") %>%
  drop_na()

tsns <- feed_trait_data$Submitted_TSN
hierarchy_full <- function(Genus, db = "itis") {
  uids <- get_tsn(Genus, messages = FALSE)
  classifications <- taxize::classification(tsns, db = db)
  taxonomy_df <- do.call(rbind, lapply(names(classifications), function(Genus) {
    cl <- classifications[[Genus]]
    if (!is.null(cl)) {
      cl$Genus <- Genus
      return(cl)
    } else {
      return(NULL)
    }
  }))
  
  return(taxonomy_df)
}

taxonomy <- taxize::classification(as.character(tsns), db = "itis")

out_dir <- here("Data/Cleaned_RDS")
tax_file <- "CONUS_Full_Taxonomy_Genus.rds" 
file_path <- file.path(out_dir, tax_file)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(file_path)) {
  saveRDS(taxonomy, file_path)
  message("Full Taxonomy Saved")
} else {
  message("File already exists")
}


extract_lineage <- function(hierarchy_df) {
  desired_ranks <- c("kingdom", "phylum", "class", "order", "family", "genus")
  lineage <- setNames(rep(NA_character_, length(desired_ranks)), desired_ranks)
  present_ranks <- hierarchy_df$rank
  present_names <- hierarchy_df$name
  for (i in seq_along(present_ranks)) {
    rank <- present_ranks[i]
    if (rank %in% desired_ranks) {
      lineage[rank] <- present_names[i]
    }
  }
  return(as.data.frame(t(lineage), stringsAsFactors = FALSE))
}

taxa_recall <- readRDS(here("CONUS_Full_Taxonomy_Genus.rds"))
lineage_df <- bind_rows(lapply(taxa_recall, extract_lineage))
lineage_df_clean <- lineage_df %>%
  filter(if_all(everything(), ~ !is.na(.)))

taxonomy_trait_matched <- feed_trait_data %>%
  left_join(lineage_df_clean, by  = c("Genus" = "genus")) %>%
  filter(!is.na(kingdom), !is.na(phylum), !is.na(class), !is.na(order), !is.na(family)) %>%
  distinct()



Taxon_name <- c("Diptera", "Ephemeroptera", "Plecoptera", "Trichoptera", 
                "Odonata", "Coleoptera", "Megaloptera", "Neuroptera", 
                "Lepidoptera", "Hemiptera")

Taxon_Filter <- taxonomy_trait_matched %>%
  filter(order %in% Taxon_name)

arthropod_tree <- read.tree(here("arthropods_genus.nwk"))
arthropod_tree_extract <- sub("-.*", "", arthropod_tree$tip.label)

target_taxa <- taxonomy_trait_matched %>% 
  mutate(Genus = sub("[^A-Za-z].*", "", Genus)) %>% 
  pull(Genus) %>% 
  unique()


target_taxa <- taxonomy_trait_matched %>% 
  mutate(Genus = gsub(" ", "_", Genus)) %>%
  pull(Genus) %>% 
  unique()

matching_tips <- arthropod_tree$tip.label[arthropod_tree_extract %in% target_taxa]
subtree <- drop.tip(arthropod_tree, setdiff(arthropod_tree$tip.label, matching_tips)) 
base_labels <- sub("-\\d+$", "", subtree$tip.label)
dup_tips <- subtree$tip.label[duplicated(base_labels)]
subtree_cleaned <- drop.tip(subtree, dup_tips)
subtree_cleaned$tip.label <- sub("-\\d+$", "", subtree_cleaned$tip.label)

subtree_labels <- c(subtree_cleaned$tip.label)

labels_df <- tibble(label = sub("[^A-Za-z].*", "", subtree_cleaned$tip.label))
labels_df <- tibble(label = subtree_cleaned$tip.label)

taxonomy_trait_matched_clean <- taxonomy_trait_matched %>%
  mutate(Genus = sub("[^A-Za-z].*", "", Genus))



tip_data <- labels_df %>% 
  left_join(taxonomy_trait_matched_clean, by = c("label" = "Genus"), 
            relationship = "many-to-many") %>% 
  filter(!is.na(Trait))

tips_to_drop <- subtree_cleaned$tip.label[is.na(tip_data$Trait)]
final_tree <- drop.tip(subtree_cleaned, tips_to_drop)


tip_data <- tip_data %>%
  filter(!is.na(Trait)) %>%
  filter(label %in% final_tree$tip.label)

tip_traits_wide <- tip_data %>%
  select(label, Trait, Trait_affinity) %>%
  group_by(label, Trait) %>%
  summarise(Trait_affinity = mean(Trait_affinity), .groups = "drop") %>%
  pivot_wider(
    names_from = Trait,
    values_from = Trait_affinity,
    values_fill = 0
  )

ffg_matrix_norm <- as.matrix(tip_traits_wide[,2:7])
rownames(ffg_matrix_norm) <- tip_traits_wide$label
tree_plot <- ggtree(subtree_cleaned)

png("US_PhyloHeat.png", width = 1400, height = 3000)
ggtree::gheatmap(
  tree_plot,
  ffg_matrix_norm,
  offset = 200,    
  width = 0.5,
  colnames_angle = 0, 
  colnames_position = "top",
  colnames_offset_y = 0.5,
  font.size = 3) +
  geom_tiplab(align = TRUE, size = 3, offset = 0.01) +
  scale_fill_viridis(option = "viridis") +
  theme(    legend.text = element_text(size = 12), 
            legend.title = element_text(size = 14),  
            legend.key.size = unit(0.2, "cm"))
dev.off()




corr_matrix <- vcv.phylo(subtree_cleaned, corr = TRUE)
corr_matrix_pd <- corr_matrix + diag(1e-5, nrow(corr_matrix))
tip_traits_wide$label <- factor(tip_traits_wide$label, levels = rownames(corr_matrix_pd))

cols <- c("PR","CG","CF","HB","SH","PA")
epsilon <- 1e-4
tip_traits_wide[cols] <- tip_traits_wide[cols] + epsilon
tip_traits_wide[cols] <- tip_traits_wide[cols] /
  rowSums(tip_traits_wide[cols])

tip_traits_wide$diet <- as.matrix(tip_traits_wide[cols])

##tip_traits_wide[cols] <- tip_traits_wide[cols] / 100

tip_traits_wide$diet <- with(tip_traits_wide,
                                cbind(PR, CG, CF, HB, SH, PA))

tip_traits_wide <- tip_traits_wide %>%
  select(-(2:7))



form <- bf(
  diet ~ 1 + (1 | gr(label, cov = corr_matrix_pd))
)

prior <- get_prior(
  form,
  data = tip_traits_wide,
  family = dirichlet(),
  data2 = list(corr_matrix_pd = corr_matrix_pd)
)


fit <- brm(form,
           data = tip_traits_wide,
           family = dirichlet(),
           data2 = list(corr_matrix_pd = corr_matrix_pd),
           prior = prior,
           cores = 4, chains = 4,  iter = 4000,  
           control = list(
             adapt_delta = 0.8,
             max_treedepth = 10)) 

summary(fit)

posterior <- as_draws_df(fit)

phylo_cols <- grepl("^sd_label__", colnames(posterior))
sd_phylo_mat <- posterior[, phylo_cols]
var_phylo_mat <- sd_phylo_mat^2
phylo_var <- rowSums(var_phylo_mat)
phi_samples <- posterior$phi

mu_mat <- tip_traits_wide %>% select(2) %>% as.matrix()
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
  geom_point(aes(color = trait == "PR"), size = 5) +      
  geom_errorbar(aes(ymin = lci, ymax = uci, color = trait == "PR"), width = 0.4, 
                linewidth = 2) + 
  scale_color_manual(values = c("FALSE" = "steelblue", "TRUE" = "red"), guide = "none") +
  scale_x_discrete(labels =  
                     c("muCG" = "CG",
                   "muCF" = "CF",
                   "muHB" = "HB",
                   "muSH" = "SH",
                   "muPA" = "PA",
                   "ref_trait" = "PR")) +                       
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


