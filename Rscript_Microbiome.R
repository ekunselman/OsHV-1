# ANALYSIS OF OSHV-1 EXPERIMENT 16S SEQUENCING DATA

setwd("/Users/emilykunselman/d_OsHV1_Exp2_16S/")
library(qiime2R)
physeq<-qza_to_phyloseq(features="exp_analysis/table-oshv1-dada2-filt.qza", 
                        tree="sepp-rooted-tree.qza", 
                        taxonomy="taxonomy-gg2.qza", 
                        metadata= "Oshv1_exp2_16S_sample_metadata.tsv")

# Alpha Diversity -----------------

library(phyloseq)
library(RColorBrewer)
# plot shannon diversity 
plot_richness(physeq, color = "condition_ABCD", x = "condition_ABCD", measures = c("Shannon", "Simpson"))+
  geom_boxplot(aes(fill = condition_ABCD))+
  scale_colour_brewer(palette = "PuOr", direction = -1)+
  scale_fill_brewer(palette = "PuOr", direction = -1)+
  theme_bw(base_size = 15)
# plot richness by eveness values of different samples: https://journals.asm.org/doi/10.1128/msphere.01019-20
 #First, estimate richness and evenness
alpha1<-estimate_richness(physeq, measures = c("Shannon", "Simpson", "Observed"))
library(microbiome)
pielou<-evenness(physeq, index = "pielou", zeroes = TRUE, detection = 0)
 #Then, bind data frames together and bind them with metadata
alpha <- cbind(alpha1, pielou)
alpha.metadata <- cbind(alpha, sample_data(physeq))
 #Plot the richness on x and eveness on y and color dots by condition
library(ggplot2)
ggplot(data = alpha.metadata, aes(x = Observed, y = pielou, color = condition_ABCD))+
  geom_point(size = 5)+
  xlab("Richness (Observed ASVs)")+
  ylab("Evenness (Pielou's)")+
  scale_colour_brewer(palette = "PuOr", direction = -1)+
  scale_fill_brewer(palette = "PuOr", direction = -1)+
  theme_bw(base_size = 15)
 #statistical tests between richness and evenness
 #test for normality and even variance
hist(alpha.metadata$Observed)
shapiro.test(alpha.metadata$Observed) #normally distributed
bartlett.test(Observed ~ condition, data = alpha.metadata) #even variance
hist(alpha.metadata$pielou)
shapiro.test(alpha.metadata$pielou) # not normal
hist(alpha.metadata$Shannon)
shapiro.test(alpha.metadata$Shannon) # normally distributed
bartlett.test(Shannon ~ condition, data = alpha.metadata) #even variance

library(dunn.test)
kruskal.test(Observed~condition, data = alpha.metadata)
dunn.test(alpha.metadata$Observed, alpha.metadata$condition, method = "Holm")
kruskal.test(pielou~condition, data = alpha.metadata)
dunn.test(alpha.metadata$pielou, alpha.metadata$condition, method = "Holm")
kruskal.test(Shannon~condition, data = alpha.metadata)
dunn.test(alpha.metadata$Shannon, alpha.metadata$condition, method = "Holm")
kruskal.test(Simpson~condition, data = alpha.metadata)
dunn.test(alpha.metadata$Simpson, alpha.metadata$condition, method = "Holm")


# Beta Diversity ---------

count_to_prop <- function(x) {return( x / sum(x) )}
physeq_proportions<- transform_sample_counts(physeq, count_to_prop)
#check to confirm sample sums are 1
sample_sums(physeq_proportions)[1:5] 

wunifracD <- distance(physeq, method = "wunifrac")
wunifracO <- ordinate(physeq, method="PCoA", distance=wunifracD)
plot_ordination(physeq, wunifracO, color="condition_ABCD", title="Weighted Unifrac PCoA")+
  theme_bw(base_size = 15)+
  geom_point(size=5)+
  scale_colour_brewer(palette = "PuOr", direction = -1)

unifracD <- distance(physeq, method = "unifrac")
unifracO <- ordinate(physeq, method="PCoA", distance=unifracD)
plot_ordination(physeq, unifracO, color="condition_ABCD", title="Unweighted Unifrac PCoA")+
  theme_bw(base_size = 15)+
  geom_point(size=5)+
  scale_colour_brewer(palette = "PuOr", direction = -1)+
  stat_ellipse(type = "norm")

unifracD <- distance(physeq_proportions, method = "unifrac")
unifracO <- ordinate(physeq_proportions, method="PCoA", distance=unifracD)
plot_ordination(physeq_proportions, unifracO, color="condition_ABCD", title="Unweighted Unifrac PCoA")+
  theme_bw(base_size = 15)+
  geom_point(size=5)+
  scale_colour_brewer(palette = "PuOr", direction = -1)+
  stat_ellipse(type = "norm")

brayD <- distance(physeq, method = "bray")
brayO <- ordinate(physeq, method="PCoA", distance=brayD)
plot_ordination(physeq, brayO, color="condition", title="Unweighted Unifrac PCoA")


# statistical tests
library(vegan)
library(devtools)
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)

adonis2(wunifracD ~ sample_data(physeq)$condition)
pairwise.adonis(wunifracD, sample_data(physeq)$condition, p.adjust.m = "holm")

adonis2(unifracD ~ sample_data(physeq)$condition)
pairwise.adonis(unifracD, sample_data(physeq)$condition, p.adjust.m = "holm")
adonis2(unifracD ~ sample_data(physeq_proportions)$condition)
pairwise.adonis(unifracD, sample_data(physeq_proportions)$condition, p.adjust.m = "holm")

adonis2(brayD ~ sample_data(physeq)$condition)
pairwise.adonis(brayD, sample_data(physeq)$condition, p.adjust.m = "holm")

# Taxa Bar Plots ---------

top12 <- names(sort(taxa_sums(physeq), decreasing=TRUE))[1:12]
ps.top12 <- transform_sample_counts(physeq, function(OTU) OTU/sum(OTU))
ps.top12 <- prune_taxa(top12, ps.top12)
plot_bar(ps.top12, x="sample_identity", fill="Genus") + facet_wrap(~condition_ABCD, scales="free_x")
# 12 - 20, Family or Genus

# MIA - TSE -------

#devtools::install_github("microbiome/mia")

library(mia)
tse <- loadFromQIIME2(
  featureTableFile = "exp_analysis/table-oshv1-dada2-filt.qza",
  taxonomyTableFile = "taxonomy-gg2.qza",
  sampleMetaFile = "Oshv1_exp2_16S_sample_metadata_mia.txt",
  featureNamesAsRefSeq = FALSE,
  phyTreeFile = "sepp-rooted-tree.qza"
)

# Hierarchical clustering ----------

# distance matrix is required for hclust
wunifracD <- distance(physeq, method = "wunifrac") #terrible for clustering
unifracD <- distance(physeq, method = "unifrac")
brayD <- distance(physeq, method = "bray")

# Unweighted Unifrac distances
hclust_avg <- hclust(unifracD, method = 'complete')
plot(hclust_avg)

# cut dendogram to create desired number of clusters (4)
cut_avg <- cutree(hclust_avg, k = 5)
plot(hclust_avg)
rect.hclust(hclust_avg , k = 5)

# color lines
avg_dend_obj <- as.dendrogram(hclust_avg)
#install.packages('dendextend', dependencies = TRUE)
library(dendextend)
avg_col_dend <- color_branches(avg_dend_obj, k = 5)
plot(avg_col_dend)

# Bray Curtis Distances
hclust_avg <- hclust(brayD, method = 'complete')
plot(hclust_avg)

# cut dendogram to create desired number of clusters (4)
cut_avg <- cutree(hclust_avg, k = 4)
plot(hclust_avg)
rect.hclust(hclust_avg , k = 4)

# color lines
avg_dend_obj <- as.dendrogram(hclust_avg)
avg_col_dend <- color_branches(avg_dend_obj, k = 4)
plot(avg_col_dend)

# ANCOMBC- Differential Abundance Analysis -------------

# ANCOM-BC corrects bias from sampling fraction and sequencing efficiency
# accounts for sampling fraction by introducing a sample-specific offset term in linear regression, which is estimated from observed data
# linear regression is conducted in log scale (to deal with compositionality/ sparse matrices)
  # structural zeros are determined to be zero across that whole group and removed from the analysis for that group, but automatically declared diff. abundant and the othe two groups abundances of that taxa are compared
  # These sampling zeros are imputed by using a small pseudo-count value (e.g., 1) before analyzing the data.
# multiple pairwise comparisons are conducted and multi direction false discovery rate is controlled
# includes dunnets test which is more powerful
# https://bioconductor.org/packages/release/bioc/vignettes/ANCOMBC/inst/doc/ANCOMBC2.html

# clump alive + sick because not differentially abundant?

library(mia)
tse <- loadFromQIIME2(
  featureTableFile = "exp_analysis/table-oshv1-dada2-filt.qza",
  taxonomyTableFile = "taxonomy-gg2.qza",
  sampleMetaFile = "Oshv1_exp2_16S_sample_metadata_mia.txt",
  featureNamesAsRefSeq = FALSE,
  phyTreeFile = "sepp-rooted-tree.qza"
)
library(ANCOMBC)
library(knitr)
library(tidyverse)
library(ggplot2)

tse$condition = factor(tse$condition, levels = c("before", "alive", "sick", "dead"))
asv <- t(assay(tse))
meta_data <- data.frame(colData(tse))

out <- ancombc2(
  data = tse,
  assay_name = "counts",
  tax_level="Genus",
  fix_formula = "condition", 
  p_adj_method = "holm", 
  pseudo_sens = TRUE,
  s0_perc = 0.05,
  prv_cut = 0.1,
  lib_cut = 0, 
  group = "condition", 
  struc_zero = TRUE, 
  neg_lb = TRUE,
  iter_control = list(tol = 1e-5, max_iter = 1000, verbose = FALSE),
  em_control = list(tol = 1e-5, max_iter = 1000), # use max_iter >= 100 on real data 
  mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
  alpha = 0.05, 
  global = TRUE,
  pairwise = TRUE,
  dunnet = TRUE
)

out_no_agglom <- ancombc2(
  data = tse,
  assay_name = "counts",
  fix_formula = "condition", 
  p_adj_method = "holm", 
  pseudo_sens = TRUE,
  s0_perc = 0.05,
  prv_cut = 0.1,
  lib_cut = 0, 
  group = "condition", 
  struc_zero = TRUE, 
  neg_lb = TRUE,
  iter_control = list(tol = 1e-5, max_iter = 1000, verbose = FALSE),
  em_control = list(tol = 1e-5, max_iter = 1000), # use max_iter >= 100 on real data 
  mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
  alpha = 0.05, 
  global = TRUE,
  pairwise = TRUE,
  dunnet = TRUE
)

res <- out$res #get results
res <- out_no_agglom$res #get second test results

# df_fig_1 for filtering by q value
# filter out non significantly different taxa (q > 0.05)
res_TRUE<-filter(res, diff_conditionalive == TRUE | diff_conditionsick == TRUE | diff_conditiondead == TRUE)

# mutate is used to add new variables in a data frame
# create columns 'lfc1, lfc2, lfc3' which correspond to the rounded values of lfc for each condition if significantly different
# or 0 if not significantly different
df_fig_1 = mutate(res_TRUE, lfc1 = ifelse(diff_conditionalive == TRUE, 
                                   round(lfc_conditionalive, 2), 0),
                  lfc2 = ifelse(diff_conditionsick == TRUE, 
                                round(lfc_conditionsick, 2), 0),
                  lfc3 = ifelse(diff_conditiondead == TRUE, 
                                round(lfc_conditiondead, 2), 0),)
# convert from wide to long format
# name columns "group" and "value"
df_fig_1 = pivot_longer(df_fig_1, cols = lfc1:lfc3, 
                        names_to = "group", values_to = "value")
#re-arrange rows by taxon
df_fig = arrange(df_fig_1, taxon)

# recode lfc1,2,3 to what comparison they actually show
df_fig$group = recode(df_fig$group, 
                        'lfc1' = "Alive - Control",
                        'lfc2' = "Sick - Control",
                        'lfc3' = "Dead - Control")
df_fig$group = factor(df_fig$group, 
                          levels = c("Alive - Control",
                                     "Sick - Control",
                                     "Dead - Control"))
# set axis limits for gradient scale
lo = floor(min(df_fig$value))
up = ceiling(max(df_fig$value))
mid = (lo + up)/2

# plot heatmap
fig = df_fig %>%
  ggplot(aes(x = group, y = taxon, fill = value)) + 
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       na.value = "white", midpoint = mid, limit = c(lo, up),
                       name = NULL) +
  geom_text(aes(group, taxon, label = value), size = 4) +
  scale_color_identity(guide = FALSE) +
  labs(x = NULL, y = NULL, title = "Log fold changes from unexposed to OsHV-1-exposed oysters") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
fig


# Random Forest ------------
library(qiime2R)
physeq<-qza_to_phyloseq(features="exp_analysis/table-oshv1-dada2-filt.qza", 
                        tree="sepp-rooted-tree.qza", 
                        taxonomy="taxonomy-gg2.qza", 
                        metadata= "Oshv1_exp2_16S_sample_metadata.tsv")

#https://rpubs.com/michberr/randomforestmicrobe
#install.packages("randomForest")

library(ggplot2)
library(vegan)
library(dplyr)
library(magrittr)
library(scales)
library(grid)
library(reshape2)
library(phyloseq)
library(knitr)
library(randomForest)


# keep only samples that we want to test in this model

physeqA<-subset_samples(physeq, condition == "alive" | condition == "dead")

# prune (remove) any taxa not present in current sample set
physeq2 <- prune_taxa(taxa_sums(physeqA) > 0, physeqA)

# prune rare taxa (found in less than 10% of samples)
taxa_to_keep <- genefilter_sample(physeq2, filterfun_sample(function(x) x > 1), A=0.1*nsamples(physeq2))
physeq2_filt = prune_taxa(taxa_to_keep, physeq2)

#count_to_prop <- function(x) {return( x / sum(x) )}
#data_proportions<- transform_sample_counts(physeq2, count_to_prop)
#taxa_to_keep <- taxa_names(data_proportions)[taxa_sums(data_proportions)>0.01]
#physeq_pruned <- prune_taxa(taxa_to_keep, physeq2)

ntaxa(physeq2)
ntaxa(physeq2_filt)
sample_sums(physeq2)
sample_sums(physeq2_filt)
taxa_sums(physeq2)
taxa_sums(physeq2_filt)


# Make a dataframe of training data with OTUs as column and samples as rows
predictors <- t(otu_table(physeq2_filt))
dim(predictors)

# Make one column for our outcome/response variable 
response <- as.factor(sample_data(physeq2_filt)$condition)

# Combine them into 1 data frame
rf.data <- data.frame(response, predictors)

# run random forest package
# uses 2/3 of data for training and 1/3 for validation/testing
# the 1/3 is referred to as "out of bag"

set.seed(2)
erie.classify <- randomForest(response~., data = rf.data, ntree = 100)
print(erie.classify)

# out of bag error = percentage of time the model is wrong?
Call:
  randomForest(formula = response ~ ., data = rf.data, ntree = 100) 
Type of random forest: classification
Number of trees: 100
No. of variables tried at each split: 44

OOB estimate of  error rate: 10%
Confusion matrix:
  alive dead class.error
alive     8    2         0.2
dead      0   10         0.0

# now lets see which OTUs are most important for predicting condition
# GINI coefficient is a measure of node purity (decrease in GINI due to condition)
# Make a data frame with predictor names and their importance
imp <- importance(erie.classify)
imp <- data.frame(predictors = rownames(imp), imp)

# Order the predictor levels by importance
imp.sort <- arrange(imp, desc(MeanDecreaseGini))
imp.sort$predictors <- factor(imp.sort$predictors, levels = imp.sort$predictors)

# Select the top 20 predictors
imp.20 <- imp.sort[1:20, ]
write_csv(imp.20, "top20_RF_dead_alive_filt.csv")

# ggplot
ggplot(imp.20, aes(x = predictors, y = MeanDecreaseGini)) +
  geom_bar(stat = "identity", fill = "indianred") +
  coord_flip() +
  ggtitle("Most important OTUs for classifying oyster samples\n into their alive or dead")

# What are those OTUs?
otunames <- imp.20$predictors
otunames
# remove X from names
# Check the levels before modification
print(levels(otunames))
# Remove "X" from the level names
levels(otunames) <- sub("X", "", levels(otunames))
# Check the levels after modification
print(levels(otunames))

r <- rownames(tax_table(physeq_pruned)) %in% otunames
kable(tax_table(physeq_pruned)[r, ])

write.csv(RF, "RF_dead_alive.csv", row.names = FALSE)

# Repeat for before vs after (but removing dead/decaying oysters)_________________

# filter out samples that we don't want to test in this model

physeqB<-subset_samples(physeq,condition != "dead")

# prune (remove) any taxa not present in current sample set
physeq2 <- prune_taxa(taxa_sums(physeqB) > 0, physeqB)

# remove OTUs not found in more than 10% of samples
taxa_to_keep <- genefilter_sample(physeq2, filterfun_sample(function(x) x > 1), A=0.1*nsamples(physeq2))
physeq2_filt = prune_taxa(taxa_to_keep, physeq2)

# Make a dataframe of training data with OTUs as column and samples as rows
predictors <- t(otu_table(physeq2_filt))
dim(predictors)

# Make one column for our outcome/response variable 
response <- as.factor(sample_data(physeq2_filt)$OsHV1_exp)

# Combine them into 1 data frame
rf.data <- data.frame(response, predictors)

# run random forest package
# uses 2/3 of data for training and 1/3 for validation/testing
# the 1/3 is referred to as "out of bag"

set.seed(2)
erie.classify <- randomForest(response~., data = rf.data, ntree = 100)
print(erie.classify)

# out of bag error = percentage of time the model is wrong?
Call:
  randomForest(formula = response ~ ., data = rf.data, ntree = 100) 
Type of random forest: classification
Number of trees: 100
No. of variables tried at each split: 42

OOB estimate of  error rate: 3.45%
Confusion matrix:
  After Before class.error
After     20      0   0.0000000
Before     1      8   0.1111111

# now lets see which OTUs are most important for predicting condition
# GINI coefficient is a measure of node purity (decrease in GINI due to condition)
# Make a data frame with predictor names and their importance
imp <- importance(erie.classify)
imp <- data.frame(predictors = rownames(imp), imp)

# Order the predictor levels by importance
imp.sort <- arrange(imp, desc(MeanDecreaseGini))
imp.sort$predictors <- factor(imp.sort$predictors, levels = imp.sort$predictors)

# Select the top 20 predictors
imp.20 <- imp.sort[1:20, ]
write_csv(imp.20, "top20_RF_before_after_filt.csv")

# ggplot
ggplot(imp.20, aes(x = predictors, y = MeanDecreaseGini)) +
  geom_bar(stat = "identity", fill = "indianred") +
  coord_flip() +
  ggtitle("Most important OTUs for classifying oyster samples\n into before or after OsHV-1 exposure")

# What are those OTUs?
otunames <- imp.20$predictors
otunames
# remove X from names
# Check the levels before modification
print(levels(otunames))
# Remove "X" from the level names
levels(otunames) <- sub("X", "", levels(otunames))
# Check the levels after modification
print(levels(otunames))

r <- rownames(tax_table(physeq_pruned)) %in% otunames
kable(tax_table(physeq_pruned)[r, ])

# From ChatGPT
In the context of Random Forest, the "Mean Decrease in Gini" is a measure used to assess the importance of a variable (feature) in making predictions. Gini impurity is a criterion used for splitting nodes in decision trees. The Mean Decrease in Gini for a variable provides insights into how much the variable contributes to the overall reduction in impurity across all decision trees in the Random Forest.

a breakdown of the concept:

Gini Impurity:
Gini impurity is a measure of how often a randomly chosen element would be incorrectly classified. In decision trees, it is used as a criterion for choosing the best split at each node. A lower Gini impurity indicates a better split.

Mean Decrease in Gini:
In a Random Forest, multiple decision trees are trained on different subsets of the data. For each tree, the Gini impurity is measured for each split, and the decrease in Gini impurity caused by a variable is recorded. The "Mean Decrease in Gini" is the average of these decreases across all trees.

A higher Mean Decrease in Gini for a variable indicates that the variable is more important in making accurate predictions.
A lower Mean Decrease in Gini suggests that the variable is less influential.
Interpretation:

Higher values: Variables with higher Mean Decrease in Gini are more important for making predictions. They contribute more to the reduction of impurity, and their presence in a split provides more discriminatory power.
Lower values: Variables with lower Mean Decrease in Gini are less important. Their contribution to impurity reduction is relatively smaller.
Feature Importance Ranking:
The Mean Decrease in Gini can be used to rank the variables in terms of importance. Variables with higher values are often considered more influential in the model.

In summary, the Mean Decrease in Gini helps to identify the most informative variables in a Random Forest model. It is a useful metric for feature selection and understanding the impact of different features on the overall model performance.


# Prepare Tables for manuscript

# merge taxonomic classification with predictor OTU IDs
BeforeAfter<- read.csv("top20_RF_before_after_filt.csv", row.names = 1)
DeadAlive<- read.csv("top20_RF_dead_alive_filt.csv" ,row.names = 1)
T<-as.data.frame(tax_table(physeq2_filt))

BA<-merge(BeforeAfter, T, by = 0, all.x = FALSE)
DA<-merge(DeadAlive, T, by = 0, all.x = FALSE)

write_csv(BA, "top20_RF_before_after_filt.csv")
write_csv(DA, "top20_RF_dead_alive_filt.csv")

# SPARCC -----------

table <- read.csv("exp_analysis/exported-feature-table/table.from_biom.csv", check.names = FALSE)

#change first column name to OTU.ID
colnames(table)[1]  <- "OTU.ID"
#insert csv file of your taxonomy file 
#I literally had to go into the file and add ;p__ and ;c__ to every taxa that only had d__ or d and p
#then I also removed taxa strings like "; o__; f__; g__; s__"
taxa <- read.csv("exported-gg-taxonomy/sparcc-taxonomy.csv")
#again change first column name to OTU.ID to match the other file
colnames(taxa)[1]  <- "OTU.ID"
#merge the two files
data <- merge(taxa, table, by.x="OTU.ID")
#get rid of OTU.ID column
data <- data[-c(1)]
#combine rows with duplicate taxonomy
data <- data %>%
  group_by(Taxon) %>%
  summarise_all((sum))
# Run this function - this selects the species part of the taxonomy string (s__ and then whatever the name is)

select_label_sp <- function(taxa_vector)
{
  #   Initializing variables genus_string, family_string, order_string &
  # class_string to FALSE.
  species_string <- FALSE
  genus_string <- FALSE
  family_string <- FALSE
  order_string <- FALSE
  class_string <- FALSE
  # Initializing taxa_label to empty space.
  taxa_label <- ""
  # Search for the species string pattern "s__".
  if(TRUE %in% grepl("s__", taxa_vector))
  {
    # Obtaininig the label.
    taxa_label <- taxa_vector[grepl("s__", taxa_vector) == TRUE]
    # Setting species_string to TRUE because the species is present in the string.
    species_string <- TRUE
    # Return the taxa label.
    return(taxa_label)
  }
  #search for genus string pattern "g__"
  if((species_string == FALSE) & (TRUE %in% grepl("g__", taxa_vector)))
  {
    # Obtaininig the label.
    taxa_label <- taxa_vector[grepl("g__", taxa_vector) == TRUE]
    # Setting genus_string to TRUE because the genus is present in the string.
    genus_string <- TRUE
    # Return the taxa label.
    return(taxa_label)
  }
  # Search for the family string pattern "f__".
  if((genus_string == FALSE) & (TRUE %in% grepl("f__", taxa_vector)))
  {
    # Obtaininig the label.
    taxa_label <- taxa_vector[grepl("f__", taxa_vector) == TRUE]
    # Setting family_string to TRUE because the family is present in the string.
    family_string <- TRUE
    # Return the taxa label.
    return(taxa_label)
  }
  # Search for the order string pattern "o__".
  if((family_string == FALSE) & (TRUE %in% grepl("o__", taxa_vector)))
  {
    # Obtaininig the label.
    taxa_label <- taxa_vector[grepl("o__", taxa_vector) == TRUE]
    # Setting order_string to TRUE because the order is present in the string.
    order_string <- TRUE
    # Return the taxa label.
    return(taxa_label)
  }
  # Search for the class string pattern "c__".
  if((order_string == FALSE) & (TRUE %in% grepl("c__", taxa_vector)))
  {
    # Obtaininig the label.
    taxa_label <- taxa_vector[grepl("c__", taxa_vector) == TRUE]
    # Setting class_string to TRUE because the class is present in the string.
    class_string <- TRUE
    # Return the taxa label.
    return(taxa_label)
  }
  # search for phylum string pattern "p__"
  if((species_string == FALSE) & (TRUE %in% grepl("p__", taxa_vector)))
  {
    # Obtaininig the label.
    taxa_label <- taxa_vector[grepl("p__", taxa_vector) == TRUE]
    # Setting genus_string to TRUE because the genus is present in the string.
    genus_string <- TRUE
    # Return the taxa label.
    return(taxa_label)
  }
  #search for kingdom string pattern "d_"
  if((species_string == FALSE) & (TRUE %in% grepl("d_", taxa_vector)))
  {
    # Obtaininig the label.
    taxa_label <- taxa_vector[grepl("d_", taxa_vector) == TRUE]
    # Setting genus_string to TRUE because the genus is present in the string.
    genus_string <- TRUE
    # Return the taxa label.
    return(taxa_label)
  }
}

# Defining an empty vector to add collected taxa labels.
labels <- vector()

# Defining an empty data.frame to convert labels to a data.frame. 
taxa_labels <- data.frame()

# Traverse the taxonomy information.
for(Taxon in data$Taxon)
{
  # Split taxonomy using ; as the split character.
  taxa_list <- strsplit(Taxon, ';')
  # Convert the taxa_list into a vector (taxa_vector).
  taxa_vector <- unlist(taxa_list, recursive = TRUE, use.names = TRUE)
  
  # Obtain the available label according to the information provided.
  #   The search order is as follows:
  #     1) genus string pattern "g__".
  #     2) family string pattern "f__".
  #     3) order string pattern "o__".
  #     4) class string pattern "c__".
  taxa_label <- select_label_sp(taxa_vector)
  # Add the taxa_label to the vector labels that contains all the labels
  # for the data.
  labels <- c(labels, taxa_label)
  # Covert the vector labels into the dataframe taxa_labels.
  taxa_labels <- as.data.frame(labels)
  # Set the name of the column of the dataframe taxa_labels.
  colnames(taxa_labels) <- c("taxa_labels")
}

# Add the column taxalabels to the biom_file.
data$taxa_labels <- taxa_labels$taxa_labels

#make a new column called OTU_ID that gets rid of the s__
data$OTU_ID <- substring(data$taxa_labels, 4)

#relocate all character columns to front
data <- data %>% relocate(Taxon, taxa_labels, OTU_ID)

#get rid of first two columns
data <- data[-c(1)]
data <- data[-c(1)]

#combine rows with duplicate species
data <- data %>%
  group_by(OTU_ID) %>%
  summarise_all((sum))

#save csv
write.csv(data, "species_otu_table.csv", row.names = FALSE)

# download cytoscape

library(tidyverse)
#install.packages("SpiecEasi")
library(SpiecEasi)
library(phyloseq)
#install.packages("RCy3") - not available
#library(RCy3) 
library(igraph)

species <- read.csv("exp_analysis/sparcc/species_otu_table.csv", check.names = FALSE)
cytoscapePing()
