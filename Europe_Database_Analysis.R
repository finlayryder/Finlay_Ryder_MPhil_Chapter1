
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

##BiocManager::install("ggtree")


library("here")
library("tidyverse")
library("ggExtra")
library("ggtree")
library("pheatmap")
library("viridis")
library("vegan")
library("ape")
library("brms")


european_data_all <- read.csv("Europe_Data_Clean.csv")

european_data_all <- european_data_all %>%
  mutate(Genus = word(Taxon, 1))



##Tree Construction

arthropod_tree <- read.tree("arthropods_genus.nwk")
arthropod_tree_extract <- sub("-.*", "", arthropod_tree$tip.label)
target_taxa <- unique(european_data_all$Genus)

matching_tips <- arthropod_tree$tip.label[arthropod_tree_extract %in% target_taxa]
subtree <- drop.tip(arthropod_tree, setdiff(arthropod_tree$tip.label, matching_tips)) 
base_labels <- sub("-\\d+$", "", subtree$tip.label)
dup_tips <- subtree$tip.label[duplicated(base_labels)]
subtree_cleaned <- drop.tip(subtree, dup_tips)
subtree_cleaned$tip.label <- sub("-\\d+$", "", subtree_cleaned$tip.label)

subtree_labels <- c(subtree_cleaned$tip.label)
labels_df <- tibble(label = sub("[^A-Za-z].*", "", subtree_cleaned$tip.label))
labels_df <- tibble(label = subtree_cleaned$tip.label)


tip_data <- european_data_all %>% filter(european_data_all$Genus %in% labels_df$label)

tip_data_condensed <- tip_data %>%
  group_by(Genus) %>%
  summarise(
    gra = sum(gra, na.rm = TRUE),
    min = sum(min, na.rm = TRUE),
    xyl = sum(xyl, na.rm = TRUE),
    shr = sum(shr, na.rm = TRUE), 
    gat = sum(gat, na.rm = TRUE),
    aff = sum(aff, na.rm = TRUE),
    pff = sum(pff, na.rm = TRUE),
    pre = sum(pre, na.rm = TRUE), 
    par = sum(par, na.rm = TRUE),
    oth = sum(oth, na.rm = TRUE)
  )

ffg_cols <- c(
  gra = "#1b9e77",   
  min = "#d95f02",   
  xyl = "#7570b3",  
  shr = "#e7298a",  
  gat = "#66a61e",  
  aff = "#e6ab02",   
  pff = "#a6761d",   
  pre = "#e41a1c",  
  par = "#377eb8",   
  oth = "#999999")


ffg_matrix_norm <- as.matrix(tip_data_condensed[,2:11])
ffg_matrix_norm <- ffg_matrix_norm / rowSums(ffg_matrix_norm)
rownames(ffg_matrix_norm) <- tip_data_condensed$Genus

tree_plot <- ggtree(subtree_cleaned)

png("All_FreshEcol_PhyloHeat.png", width = 1400, height = 2000)
ggtree::gheatmap(
tree_plot,
ffg_matrix_norm,
offset = 0.3,    
width = 0.5,
colnames_angle = 0, 
colnames_position = "top",
colnames_offset_y = 0.5,
font.size = 3) +
  ggtitle("All FFG Trait Affinity - Species Level (European)") +
  geom_tiplab(align = TRUE, size = 2, offset = 0.01) +
  scale_fill_viridis(option = "viridis") +
  theme(    legend.text = element_text(size = 12), 
            legend.title = element_text(size = 14),  
            legend.key.size = unit(1.5, "cm"))
dev.off()





arthropod_tree <- read.tree(here("arthropods_genus.nwk"))
arthropod_tree_extract <- sub("-.*", "", arthropod_tree$tip.label)
target_taxa <- ffg_indices %>% 
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
tip_data <- labels_df %>% 
  left_join(ffg_indices, by = c("label" = "Genus"))

ffg_indices_filtered <- ffg_indices %>% filter(ffg_indices$Genus %in% matching_tips)


ffg_indices_matrix <- as.matrix(ffg_indices_filtered)
ffg_indices_matrix <- data.frame(FFG_index = ffg_indices_matrix)
rownames(ffg_indices_matrix) <- subtree_cleaned$tip.label
ffg_indices_matrix[] <- lapply(ffg_indices_matrix, as.numeric)


tree_plot <- ggtree(subtree_cleaned)

png("FFG_Index_FreshEcol_PhyloHeat.png", width = 1400, height = 2000)
gheatmap(
  tree_plot,
  ffg_indices_matrix,
  offset = 150,    
  width = 0.1,
  colnames_angle = 0, 
  colnames_position = "top",
  colnames_offset_y = 0.5,
  font.size = 5) +
  ggtitle("FFG Affinity Index (European)") +
  scale_fill_viridis(option = "viridis", limits = c(0,1))
dev.off()





## PGLMM

corr_matrix <- vcv.phylo(subtree_cleaned, corr = TRUE)
corr_matrix_pd <- corr_matrix + diag(1e-5, nrow(corr_matrix))
tip_data_condensed$Genus <- factor(tip_data_condensed$Genus, levels = rownames(corr_matrix_pd))

tip_data_condensed %>%
  summarise(across(c(gra,min,xyl,shr,gat,aff,pff,pre,par,oth),
                   sum, na.rm = TRUE))


cols <- c("pre","min","xyl","shr","gat","aff","pff","gra","par","oth")
epsilon <- 1e-4
tip_data_condensed[cols] <- tip_data_condensed[cols] + epsilon
tip_data_condensed[cols] <- tip_data_condensed[cols] /
  rowSums(tip_data_condensed[cols])

tip_data$diet <- as.matrix(tip_data[cols])

tip_data_condensed[cols] <- tip_data_condensed[cols] / 100

tip_data_condensed$diet <- with(tip_data_condensed,
                      cbind(pre, min, xyl, shr, gat, aff, pff, gra, par, oth))

tip_data_condensed <- tip_data_condensed %>%
  select(-(2:11))



form <- bf(
  diet ~ 1 + (1 | gr(Genus, cov = corr_matrix_pd))
)

prior <- get_prior(
  form,
  data = tip_data_condensed,
  family = dirichlet(),
  data2 = list(corr_matrix_pd = corr_matrix_pd)
)


prior_set <- c(
  set_prior("normal(1, 0.5)", class = "Intercept", dpar = "muaff"),
  set_prior("normal(1, 0.5)", class = "Intercept", dpar = "mugat"),
  set_prior("normal(1, 0.5)", class = "Intercept", dpar = "mugra"),
  set_prior("normal(1, 0.5)", class = "Intercept", dpar = "muoth"),
  set_prior("normal(1, 0.5)", class = "Intercept", dpar = "mupff"),
  
  set_prior("normal(-1.5, 0.5)", class = "Intercept", dpar = "mumin"),
  set_prior("normal(-1.5, 0.5)", class = "Intercept", dpar = "mushr"),
  set_prior("normal(-1.5, 0.5)", class = "Intercept", dpar = "mupar"),
  set_prior("normal(-1.5, 0.5)", class = "Intercept", dpar = "muxyl"), 
  
  set_prior("gamma(5, 1)", class = "phi"),

  set_prior("exponential(2)", class = "sd", group = "Genus", dpar = "muaff"),
  set_prior("exponential(2)", class = "sd", group = "Genus", dpar = "mugat"),
  set_prior("exponential(2)", class = "sd", group = "Genus", dpar = "mugra"),
  set_prior("exponential(2)", class = "sd", group = "Genus", dpar = "muoth"),
  set_prior("exponential(2)", class = "sd", group = "Genus", dpar = "mupff"),
  set_prior("exponential(2)", class = "sd", group = "Genus", dpar = "mumin"),
  set_prior("exponential(2)", class = "sd", group = "Genus", dpar = "mushr"),
  set_prior("exponential(2)", class = "sd", group = "Genus", dpar = "mupar"),
  set_prior("exponential(2)", class = "sd", group = "Genus", dpar = "muxyl")
)

fit <- brm(form,
           data = tip_data_condensed,
           family = dirichlet(),
           data2 = list(corr_matrix_pd = corr_matrix_pd),
           prior = prior,
           cores = 4, chains = 4,  iter = 4000,  
           control = list(
             adapt_delta = 0.8,
             max_treedepth = 10)) 



summary(fit)





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
  geom_point(aes(color = trait == "pre"), size = 5) +      
  geom_errorbar(aes(ymin = lci, ymax = uci, color = trait == "pre"), width = 0.4, 
                linewidth = 2) + 
  scale_color_manual(values = c("FALSE" = "steelblue", "TRUE" = "red"), guide = "none") +
  scale_x_discrete(labels =  
                     c("min" = "MI",
                       "xyl" = "XY",
                       "shr" = "SH",
                       "gat" = "CG",
                       "aff" = "AFF",
                       "pff" = "PFF", 
                       "gra" = "GR",
                       "oth" = "OT",
                       "par" = "PA",
                       "pre" = "PR")) +                       
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 30, face = "bold"),  
    axis.text.y = element_text(size = 30, face = "bold"),                          
    legend.position = "bottom",
    legend.text = element_text(size = 30, face = "bold"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()                          
  ) +
  ylim(c(0,100))






