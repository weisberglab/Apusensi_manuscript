library(data.table)
library(tidyverse)
library(combinat)
library(hrbrthemes)
library(forcats)
library(ggforce)
library(ggimage)
library(ggtree)
library(ape)
library(phytools)
library(ggplot2)
library(ggpubr)
library(perm)
library(vegan)
library(doParallel)
library(foreach)

#with the double entered Af100-12_9 removed
PIRATE.gene_families <- as.data.frame(fread("data/PIRATE.gene_families.ordered.tsv")) 

PIRATE.gene_families <-  PIRATE.gene_families[,!names(PIRATE.gene_families) %in% c("Agrobacterium_fabrum_C58", "Agrobacterium_tumefaciens_MH_0_5_111223_17_1677754708", 'Agrobacterium_pusense_MGBC108980', 'Agrobacterium_pusense_SCN18_30_10_14_R3_B_60_7')]
dim(PIRATE.gene_families)

# Need to update the "number_genomes" column
PIRATE.gene_families$number_genomes <- 0
PIRATE.gene_families$number_genomes <- rowSums( !(PIRATE.gene_families[ , 23:length(colnames(PIRATE.gene_families))] == "") )
PIRATE.gene_families <- PIRATE.gene_families[ which(PIRATE.gene_families$number_genomes > 0), ]


# Distribution of genes across the strains
gene_total <- colSums( !(PIRATE.gene_families[ , 23:length(colnames(PIRATE.gene_families))] == ""))
sort(gene_total)
summary(as.vector(gene_total))
sort(as.vector(gene_total))


sd(gene_total)


##basic stats:
#how many gene families are there?
n_gene_fams<- nrow(PIRATE.gene_families)
n_gene_fams
# 24618

#how many strains?

#names of cols to exclude
cols_to_exclude<- c("allele_name",
                    "gene_family",
                    "consensus_gene_name",
                    "consensus_product",
                    "threshold",
                    "alleles_at_maximum_threshold",
                    "number_genomes",
                    "average_dose",
                    "min_dose",
                    "max_dose",
                    "genomes_containing_fissions",
                    "genomes_containing_duplications",
                    "number_fission_loci",
                    "number_duplicated_loci",
                    "no_loci",
                    "products",
                    "gene_names",
                    "min_length(bp)",
                    "max_length(bp)",
                    "average_length(bp)", 
                    "cluster", 
                    "cluster_order"
                    )

strains_only <- PIRATE.gene_families[,!names(PIRATE.gene_families) %in% cols_to_exclude]
strain_names <- names(strains_only)


#length(strain_names)
ngenomes<- length(unique(strain_names)) 
ngenomes
#84


#print to cross ref for tree building
#write.table(strain_names, "strain_names_27Sep2021.txt", sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

#how many of the gene families are in every genome (of 261)
n_gene_fams_core_all<- sum(PIRATE.gene_families$number_genomes == ngenomes)
n_gene_fams_core_all
# 1250

#that's what percent out of the total?
(n_gene_fams_core_all*100)/n_gene_fams
# 1.71895

#present in 95% of genomes (n >=248)
cutoff<- round(.95* ngenomes)
n_gene_fams_core_w95per<- sum(PIRATE.gene_families$number_genomes >= cutoff)
n_gene_fams_core_w95per
# 3785

#that's what percent out of the total?
(n_gene_fams_core_w95per*100)/n_gene_fams
# 15.37493%

#how many of the gene families are singletons (accessory)
n_gene_fams_singletons<- sum(PIRATE.gene_families$number_genomes == 1)
n_gene_fams_singletons
# 11026

#that's what percent out of the total?
(n_gene_fams_singletons*100)/n_gene_fams
# 45.08391

#how many accessory?)
n_accessory <- n_gene_fams - (n_gene_fams_singletons + n_gene_fams_core_w95per)
n_accessory
# 9807

#that's what percent out of the total?
(n_accessory*100)/n_gene_fams
# 41.18998%

#get average per genome 
n_accessory / ngenomes
# 121.5862

##graph the distribution of gene presence in a gene family (distribution of core to accessory genes)
#plot
gene_fam_totals <- as.data.frame(PIRATE.gene_families$number_genomes)
colnames(gene_fam_totals) <- 'count'

#set groups
gene_fam_totals$group = 0                        
for (i in 1:nrow(gene_fam_totals)){
  if (gene_fam_totals$count[i] == 1) {
    gene_fam_totals$group[i] = "Singleton"
  } else if (gene_fam_totals$count[i] >= .95*ngenomes) {
    gene_fam_totals$group[i] = "Core"
  } else {
    gene_fam_totals$group[i] = "Accessory"
  }
}



#plot
#pallet: 
#singelton = "#316A6E",
#accessory = "#BA9141",
#core = "#6E572C")

p1 <- gene_fam_totals %>%
  ggplot(aes(x = count)) +
  geom_bar(aes(fill = group), 
           position = "identity") +
  ggtitle("Pangenome Distribution by gene ortholog group") +
  ylab("n gene families") + xlab("n genomes in family") +
  theme_classic() +
  theme(plot.title = element_text(size=15), legend.title = element_blank(), legend.position = c(0.85, 0.85)) +
  scale_fill_manual(name = "", values = c("blueviolet", "coral",  "darkcyan")) +
  theme_minimal()
  
p1


#ggsave("dist_by_gene_fam.pdf",p, width=6, height=4, units="in")


###plot as doughnut chart 
#make df of totals
fam_dist_df <- data.frame(
  category=c("Singleton", "Accessory", "Core"),
  count=c(n_gene_fams_singletons, n_gene_fams -(n_gene_fams_singletons + n_gene_fams_core_w95per), n_gene_fams_core_w95per)
)

# Compute percentages
fam_dist_df$fraction <- fam_dist_df$count / sum(fam_dist_df$count)
# Compute the cumulative percentages (top of each rectangle)
fam_dist_df$ymax <- cumsum(fam_dist_df$fraction)
# Compute the bottom of each rectangle
fam_dist_df$ymin <- c(0, head(fam_dist_df$ymax, n=-1))
# Compute label position
fam_dist_df$labelPosition <- (fam_dist_df$ymax + fam_dist_df$ymin) / 2
# Compute a good label
fam_dist_df$label <- paste0(fam_dist_df$category, "\n", fam_dist_df$count)
#plot
p2 <- ggplot(fam_dist_df, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category)) +
  geom_rect() +
  geom_text( x=1.5, aes(y=labelPosition, label=label, color=category), size=4) + # x here controls label position (inner / outer)
  scale_fill_manual(values=c("blueviolet", "coral",  "darkcyan"))+
  scale_color_manual(values=c("blueviolet", "coral",  "darkcyan"))+
  coord_polar(theta="y") +
  xlim(c(-1, 4)) +
  theme_void() +
  theme(legend.position = "none")

p2

##################################
# Make a gene accumulation curve #
##################################

binary_df <- strains_only %>% 
  mutate_all(~ as.numeric(nzchar(.)))

#
# Another Try
#

count_all_ones <- function(df) {
  num_cols <- ncol(df)
  count_cases <- df %>%
    filter(rowSums(.) == num_cols) %>%
    nrow()
  return(count_cases)
}

find_common_genes <- function(df, cols) {
  common_genes <- count_all_ones(df[cols])
  return(common_genes)
}

# Initialize dataframe to store results
results <- data.frame(
  NumGenomes = integer(),
  CommonGenes = integer(),
  TotalGenes = integer()
)

# Initialize an empty list to store results
result_list <- list()

#setup parallel backend to use many processors
cores=detectCores()
cl <- makeCluster(cores[1]-2) #not to overload your computer

registerDoParallel(cl)


# Generate combinations and count common genes
# Initialize an empty list to store results
for(num_genomes in 2:ncol(binary_df)) {
  rand_combs <- replicate(1000, sample(names(binary_df), num_genomes), simplify = FALSE)
  print(num_genomes)
  # if (length(combs) > random_sample_size) {
  #   random_combs <- sample(combs, random_sample_size, replace = FALSE)
  # } else {
  #   random_combs <- combs
  local_results <- foreach(i = 1:length(rand_combs), .combine = rbind, .packages = 'dplyr') %dopar% {
    common_genes_count <- find_common_genes(binary_df, rand_combs[[i]])
    tota_genes_count <- sum(rowSums(binary_df[rand_combs[[i]]]) > 0)
    data.frame(NumGenomes = num_genomes, CommonGenes = common_genes_count, TotalGenes = tota_genes_count)
  }
  result_list[[num_genomes - 1]] <- local_results  # Store results in list
}


# Write a new 

# Stop the parallel backend
stopCluster(cl)

# Combine all results into one data frame
results <- do.call(rbind, result_list)

write.csv(results, file='results/accumulation_curve_data3.csv', quote = FALSE, row.names = FALSE)


# Fit power law

# Log-transform the data
results <- read.csv('results/accumulation_curve_data3.csv', header=TRUE)

head(results)

results <- results %>% 
  mutate(log_NumGenomes = log(NumGenomes), log_CommonGenes = log(CommonGenes), log_TotalGenes = log(TotalGenes))


head(results)

# Fit a linear model to the log-transformed data
lm_model_common <- lm(log_CommonGenes ~ log_NumGenomes, data = results)
lm_model_all <- lm(log_TotalGenes ~ log_NumGenomes, data = results)

# Get the coefficients
coefficients_common <- coef(lm_model_common)
intercept_common <- coefficients_common[1]
slope_common <- coefficients_common[2]

coefficients_all <- coef(lm_model_all)
intercept_all <- coefficients_all[1]
slope_all <- coefficients_all[2]

# Create a data frame for the fitted power law curve
fitted_curve <- results %>%
  mutate(FittedCommonGenes = exp(intercept_common) * NumGenomes^slope_common)

# Core genome size based on accumulation curve
exp(intercept_common) * max(results$NumGenomes)^slope_common

fitted_curve <- fitted_curve %>%
  mutate(FittedAllGenes = exp(intercept_all) * NumGenomes^slope_all)

# Plot the results
#p3 <- ggplot(results, aes(x = NumGenomes, y = CommonGenes)) +
#  geom_point(color='grey') +
#  geom_line(data = fitted_curve, aes(x = NumGenomes, y = FittedCommonGenes), color = "red") +
#  scale_color_manual(name = "",
#                     values = c("Common Genes" = "grey", "Fitted Curve" = "red"))
#  labs(title = "Core gene accumulation curve",
#       x = "Number of Genomes Compared",
#       y = "Number of Common Genes") +
#  theme_minimal()

p3 <- ggplot(results, aes(x = NumGenomes, y = CommonGenes)) +
  geom_point(aes(color = "Common Genes"), size = 3) +
  geom_line(data = fitted_curve, aes(x = NumGenomes, y = FittedCommonGenes, color = "Fitted Curve"), size = 1) +
  scale_color_manual(name = "",
                     values = c("Common Genes" = "grey", "Fitted Curve" = "red")) +
  labs(title = "Core Gene Accumulation Curve",
       x = "Number of Genomes Compared",
       y = "Number of Common Genes") +
  theme_minimal() +
  theme(legend.position = c(0.95, 0.95), # Position inside top right corner
        legend.justification = c("right", "top"))

p4 <- ggplot(results, aes(x = NumGenomes, y = TotalGenes)) +
  geom_point(aes(color = "Total Genes"), size = 3) +
  geom_line(data = fitted_curve, aes(x = NumGenomes, y = FittedAllGenes, color = "Fitted Curve"), size = 1) +
  scale_color_manual(name = "",
                     values = c("Total Genes" = "grey", "Fitted Curve" = "red")) +
  labs(title = "Pan Genome Accumulation Curve",
       x = "Number of Genomes Compared",
       y = "Number of Total Genes") +
  theme_minimal() +
  theme(legend.position = c(0.95, 0.05), # Position inside bottom right corner
        legend.justification = c("right", "bottom"))


ggarrange(
  p1,                # First row with line plot
  # Second row with box and dot plots
  ggarrange(p3, p4, ncol = 2, labels = c("B", "C")), 
  nrow = 2, 
  labels = "A"       # Label of the line plot
) 


##which strain has the highest/lowest number of singletons and accessory gene fams? 
gene_fam_by_strain<-as.data.frame(PIRATE.gene_families[,23:ncol(PIRATE.gene_families)])
ncol(gene_fam_by_strain)

##make binary (if gene = 1, if not = 0)
#fill in zeros
gene_fam_by_strain_zeros<- sapply(gene_fam_by_strain, gsub, pattern = "^\\s*$" , replacement = 0 )
#fill in ones
gene_fam_by_strain_ones<- as.data.frame(replace(gene_fam_by_strain_zeros, gene_fam_by_strain_zeros!="0", 1))
#change to numeric
gene_fam_by_strain_ones_num <- mutate_all(gene_fam_by_strain_ones, function(x) as.numeric(as.character(x)))

#subset to remove core genes from accessory and singletons 
all_accessory_1<- gene_fam_by_strain_ones_num[rowSums(gene_fam_by_strain_ones_num) > 1,]
in_95_percent<- .95 * ngenomes
all_accessory<- all_accessory_1[rowSums(all_accessory_1) < in_95_percent,]

#subset to get only singletons 
singletons_only<- gene_fam_by_strain_ones_num[rowSums(gene_fam_by_strain_ones_num) == 1,]

#get average
ave_accessory<- colSums(all_accessory)
mean(ave_accessory) #1558.598

#get average
ave_singletons<- colSums(singletons_only)
mean(ave_singletons) #133.0805




# Which strain has the largest accessory genome (genes not in the core?) 
accessory_by_strain<- as.data.frame(colSums(all_accessory))
colnames(accessory_by_strain) <- "totals"
accessory_by_strain$strain<- row.names(accessory_by_strain)
#get max
accessory_by_strain[which.max(accessory_by_strain$totals),]
#get min
accessory_by_strain[which.min(accessory_by_strain$totals),]
#sort 
accessory_by_strain<- accessory_by_strain[order(accessory_by_strain$totals),]
#mean
mean(accessory_by_strain$totals)

#fix names so that they match 
name_map<-read.delim("results/PIRATE_isolate_clade_map.tsv", header = TRUE, sep = "\t", fill = TRUE, strip.white = TRUE)
  row.names(accessory_by_strain) <- name_map$strain[match(row.names(accessory_by_strain), name_map$strain)]
accessory_by_strain$pop_name <- row.names(accessory_by_strain)
#remove "DMC2" for graphing 
accessory_by_strain$pop_name<- sapply(accessory_by_strain$pop_name, gsub, pattern = "DMC2_", replacement = "")


#graph accessory genome size by strain
p <- accessory_by_strain %>%
  mutate(name = fct_reorder(pop_name, totals)) %>%
  ggplot( aes(x=name, y=totals))+
  geom_bar(stat="identity", fill="blueviolet", alpha=5, width=1, position = position_dodge(width=0.6)) +
  xlab("") + ylab("n accessory gene families") +
  ggtitle("Accessory genome size by strain") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
        #axis.text.x = element_text(size = 1.5, angle=90, hjust=1, vjust=0.5, margin=margin(-3,0,0,0)), legend.position = "none")
        axis.text.x = element_text(size = 6, angle=75, hjust=1, vjust=0.9, margin=margin(0,0,0,0)), legend.position = "none")
p


#  theme(text=element_text(size=9), 
#        axis.text.x = element_text(size = 2, angle=90, hjust=1), legend.position = "none")+
#  facet_zoom(ylim = c(min(accessory_by_strain$totals), max(accessory_by_strain$totals)), zoom.data = ifelse(a <= 6000,  FALSE))

#ggsave("accessory.pdf",p, width=6.9, height=3, units="in")



#for singletons 
#Which strain has the largest accessory genome (genes not in the core?) 
singletons_only_by_strain <- as.data.frame(colSums(singletons_only))
colnames(singletons_only_by_strain) <- "totals"
singletons_only_by_strain$strain<- row.names(singletons_only_by_strain)
#get max
singletons_only_by_strain[which.max(singletons_only_by_strain$totals),]
#get min
singletons_only_by_strain[which.min(singletons_only_by_strain$totals),]
#how many are zero?
no_singeltons<-data.frame(singletons_only_by_strain[singletons_only_by_strain$totals == 0,])
nrow(no_singeltons)
#how many are 1?
one_singeltons<-data.frame(singletons_only_by_strain[singletons_only_by_strain$totals == 1,])
nrow(one_singeltons)
#get average 
mean(singletons_only_by_strain$totals)
#135.4773

#View(singletons_only_by_strain)
#sort 
singletons_only_by_strain<- singletons_only_by_strain[order(singletons_only_by_strain$totals),]
#fix names
row.names(singletons_only_by_strain) <- name_map$strain[match(row.names(singletons_only_by_strain), name_map$strain)]
singletons_only_by_strain$pop_name<- row.names(singletons_only_by_strain)
#remove "DMC2" for graphing 
#singletons_only_by_strain$pop_name<- sapply(singletons_only_by_strain$pop_name, gsub, pattern = "DMC2_", replacement = "")


#graph singleton genome size by strain
p <- singletons_only_by_strain %>%
  mutate(name = fct_reorder(pop_name, totals)) %>%
  ggplot( aes(x=name, y=totals)) +
  geom_bar(stat="identity", fill="darkcyan", alpha=2, width=1) +
  xlab("") + ylab("n singleton gene families") +
  ggtitle("Singleton genome size by strain") +
  theme(text=element_text(size=9), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
        #theme(text=element_text(size=9), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
        axis.text.x = element_text(size = 6, angle=75, hjust=1), legend.position = "none")
p
