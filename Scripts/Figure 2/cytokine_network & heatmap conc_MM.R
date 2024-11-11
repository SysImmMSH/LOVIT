###### LOVIT data#######

## 1. Cytokine network
## 2. Heatmap of cytokine concentrations by subtype


#read df if not in environmment already
cytokine_df_LOVIT <- read.csv("CSV files/patient_df_3clusters.csv")

##need to make an adjacency matrix
head(cytokine_df_LOVIT)

cytokine_df_LOVIT <- cytokine_df_LOVIT[, -c(2, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 
                                            45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59,
                                            60, 61, 62, 63, 64)]

#add cluster list to cytokine dataframe



head(cytokine_df_LOVIT)
colnames(cytokine_df_LOVIT)


#split dataframe by cluster
library(dplyr)
split_df <- cytokine_df_LOVIT %>%
  group_by(groups) %>%
  group_split()


cluster1_df_L <- split_df[[1]]
cluster2_df_L <- split_df[[2]]
cluster3_df_L <- split_df[[3]]

head(cluster1_df_L)

#remove sample list and cluster list from dataframes
cluster1_df_L <- cluster1_df_L[-c(1, 2, 3, 30)]
cluster2_df_L <- cluster2_df_L[-c(1, 2, 3, 30)]
cluster3_df_L <- cluster3_df_L[-c(1, 2, 3, 30)]


#cytokine network plot 
library("qgraph")
library("Hmisc")

colnames(cluster1_df_L)
colnames(cluster2_df_L)
colnames(cluster3_df_L)


#calculate the correlation matrix and p-values
cluster1_df_L_res <- rcorr(as.matrix(cluster1_df_L))
cluster2_df_L_res <- rcorr(as.matrix(cluster2_df_L))
cluster3_df_L_res <- rcorr(as.matrix(cluster3_df_L))


#extract the correlation matrix and p-values
cormat_1_r_L <- cluster1_df_L_res$r
pvals_1_L <- cluster1_df_L_res$P

cormat_2_r_L <- cluster2_df_L_res$r
pvals_2_L <- cluster2_df_L_res$P

cormat_3_r_L <- cluster3_df_L_res$r
pvals_3_L <- cluster3_df_L_res$P


#set the p-value threshold 
pval_threshold <- 0.001

#######------------------------ P value filtering only -------------------##########
#filter the correlation matrix based on the p-value threshold
cormat_1_filtered_L <- cormat_1_r_L
cormat_1_filtered_L[pvals_1_L > pval_threshold] <- 0

cormat_2_filtered_L <- cormat_2_r_L
cormat_2_filtered_L[pvals_2_L > pval_threshold] <- 0

cormat_3_filtered_L <- cormat_3_r_L
cormat_3_filtered_L[pvals_3_L > pval_threshold] <- 0

#use cormat on each dataframs to make adjacency matrix
cormat_1_L <- cor(cluster1_df_L, method = "spearman")
cormat_2_L <- cor(cluster2_df_L, method = "spearman")
cormat_3_L <- cor(cluster3_df_L, method = "spearman")


#############################################################################################
#############################################################################################
##### To make the node weighted by concentration I need to create a vector that contains the concentration levels for each node. 
#i need to aggregate the concentrations of each cytokine by median

head(cluster1_df_L)

library(data.table)
library(tidyr)
#make the data long
long<- cluster1_df_L %>%
  pivot_longer(
    cols = `IFNy`:`IL.33`,
    names_to = "cytokines",
    values_to = "concentrations"
  )

#find the median of each cytokine
aggregated_data <- long%>%
  group_by(cytokines) %>%
  summarise(concentrations = median(concentrations, na.rm = TRUE))

#extract the concentraion levels
concentration_levels <- aggregated_data$concentrations

#now i need to scale the nodes by the concentration levels
node_sizes <- log2(concentration_levels)


#############################################################################################
###### Cluster 2 
#make the data long
long2<- cluster2_df_L %>%
  pivot_longer(
    cols = `IFNy`:`IL.33`,
    names_to = "cytokines",
    values_to = "concentrations"
  )

#find the median of each cytokine
aggregated_data2 <- long2%>%
  group_by(cytokines) %>%
  summarise(concentrations = median(concentrations, na.rm = TRUE))

#extract the concentraion levels
concentration_levels2 <- aggregated_data2$concentrations

#now i need to scale the nodes by the concentration levels
node_sizes2 <- log2(concentration_levels2)

#############################################################################################
###### Cluster 3 
#make the data long
long3<- cluster3_df_L %>%
  pivot_longer(
    cols = `IFNy`:`IL.33`,
    names_to = "cytokines",
    values_to = "concentrations"
  )

#find the median of each cytokine
aggregated_data3 <- long3%>%
  group_by(cytokines) %>%
  summarise(concentrations = median(concentrations, na.rm = TRUE))

#extract the concentraion levels
concentration_levels3 <- aggregated_data3$concentrations

#now i need to scale the nodes by the concentration levels
node_sizes3 <- log2(concentration_levels3)

#############################################################################################
#Rename cytokine columsn that have greek letters
#cluster 1
colnames(cormat_1_L)
rownames(cormat_1_L)
str(cormat_1_L)

colnames(cormat_1_L)[c(1, 6, 14:26)] <- c("IFN𝛄", "TNFα", "IL-6", "IL-5", "IL-1α", "IL-18",
                                                 "IL-4", "IFNα2", "IL-10", "IL-1RA", "sTNF-RI", "IL-12p40",
                                                 "IFNβ", "IL-1β", "IL-33")

rownames(cormat_1_L)[c(1, 6, 14:26)] <- c("IFN𝛄", "TNFα", "IL-6", "IL-5", "IL-1α", "IL-18",
                                          "IL-4", "IFNα2", "IL-10", "IL-1RA", "sTNF-RI", "IL-12p40",
                                          "IFNβ", "IL-1β", "IL-33")
#cluster 2
#Rename cytokine columsn that have greek letters
colnames(cormat_2_L)

colnames(cormat_2_L)[c(1, 6, 14:26)] <- c("IFN𝛄", "TNFα", "IL-6", "IL-5", "IL-1α", "IL-18",
                                          "IL-4", "IFNα2", "IL-10", "IL-1RA", "sTNF-RI", "IL-12p40",
                                          "IFNβ", "IL-1β", "IL-33")

rownames(cormat_2_L)[c(1, 6, 14:26)] <- c("IFN𝛄", "TNFα", "IL-6", "IL-5", "IL-1α", "IL-18",
                                          "IL-4", "IFNα2", "IL-10", "IL-1RA", "sTNF-RI", "IL-12p40",
                                          "IFNβ", "IL-1β", "IL-33")

#cluster 3
#Rename cytokine columsn that have greek letters
colnames(cormat_3_L)

colnames(cormat_3_L)[c(1, 6, 14:26)] <- c("IFN𝛄", "TNFα", "IL-6", "IL-5", "IL-1α", "IL-18",
                                          "IL-4", "IFNα2", "IL-10", "IL-1RA", "sTNF-RI", "IL-12p40",
                                          "IFNβ", "IL-1β", "IL-33")

rownames(cormat_3_L)[c(1, 6, 14:26)] <- c("IFN𝛄", "TNFα", "IL-6", "IL-5", "IL-1α", "IL-18",
                                          "IL-4", "IFNα2", "IL-10", "IL-1RA", "sTNF-RI", "IL-12p40",
                                          "IFNβ", "IL-1β", "IL-33")
#############################################################################################
# below only plots those that are significant
## THIS IS USED FOR LOVIT - 16.10.24

#cluster 1
png("Figures/cytokine_network_clus1.png", width = 275, height = 225, units='mm', res = 300) 
qgraph(cormat_1_L,graph = "cor", shape="circle", posCol="red", negCol="blue", layout="spring", repulsion=0.6, vsize=node_sizes, labels = colnames(cormat_1_L),
       label.scale.equal = TRUE, label.cex = 3, threshold = "bonferroni",alpha = 0.0001, sampleSize = nrow(cluster1_df_L),
       minimum = "sig", edge.width = 1)
dev.off()

 #cluster2
png("Figures/cytokine_network_clus2.png", width = 275, height = 225, units='mm', res = 300) 
qgraph(cormat_2_L,graph = "cor", shape="circle", posCol="red", negCol="navy", layout="spring", repulsion=0.65, vsize=node_sizes2, labels = colnames(cormat_2_L),
       label.scale.equal = TRUE, label.cex = 35, threshold = "bonferroni",alpha = 0.0001, sampleSize = nrow(cluster2_df_L),
       minimum = "sig", edge.width = 1)
dev.off()

#cluster 3
png("Figures/cytokine_network_clus3.png", width = 275, height = 225, units='mm', res = 300) 
qgraph(cormat_3_L,graph = "cor", shape="circle", posCol="red", negCol="navy", layout="spring", repulsion=0.85, vsize=node_sizes3, labels = colnames(cormat_3_L),
       label.scale.equal = TRUE, label.cex = 2, threshold = "bonferroni",alpha = 0.0001, sampleSize = nrow(cluster3_df_L),
       minimum = "sig", edge.width = 1)
dev.off()


########################################################################################################
#######################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################

#Make a heatmap of the different concentrations of cytokines by subtype

#Need to combine concentration lists - the aggregated data
library(dplyr)

Concentrations_bySubtype <- bind_cols(aggregated_data, aggregated_data2, aggregated_data3)

Concentrations_bySubtype <- Concentrations_bySubtype %>%
  select(-cytokines...3, -cytokines...5) %>%
  rename(Cytokine = cytokines...1,
         `Subtype 1` = concentrations...2,
         `Subtype 2` = concentrations...4,
         `Subtype 3` = concentrations...6)

head(Concentrations_bySubtype_hm)

head(Concentrations_bySubtype)

#####change the cytokine names to geek letters, 

# Replace "Alice" with "Alicia"
Concentrations_bySubtype$Cytokine <- replace(Concentrations_bySubtype$Cytokine, Concentrations_bySubtype$Cytokine == "IFNy", "IFN𝛄")
Concentrations_bySubtype$Cytokine <- replace(Concentrations_bySubtype$Cytokine, Concentrations_bySubtype$Cytokine == "TNFa", "TNFα")
Concentrations_bySubtype$Cytokine <- replace(Concentrations_bySubtype$Cytokine, Concentrations_bySubtype$Cytokine == "IL.6", "IL-6")
Concentrations_bySubtype$Cytokine <- replace(Concentrations_bySubtype$Cytokine, Concentrations_bySubtype$Cytokine == "IL.5", "IL-5")
Concentrations_bySubtype$Cytokine <- replace(Concentrations_bySubtype$Cytokine, Concentrations_bySubtype$Cytokine == "IL.1a", "IL-1α")
Concentrations_bySubtype$Cytokine <- replace(Concentrations_bySubtype$Cytokine, Concentrations_bySubtype$Cytokine == "IL.18", "IL-18")
Concentrations_bySubtype$Cytokine <- replace(Concentrations_bySubtype$Cytokine, Concentrations_bySubtype$Cytokine == "IL.4", "IL-4")
Concentrations_bySubtype$Cytokine <- replace(Concentrations_bySubtype$Cytokine, Concentrations_bySubtype$Cytokine == "IFNa2", "IFNα2")
Concentrations_bySubtype$Cytokine <- replace(Concentrations_bySubtype$Cytokine, Concentrations_bySubtype$Cytokine == "IL.10", "IL-10")
Concentrations_bySubtype$Cytokine <- replace(Concentrations_bySubtype$Cytokine, Concentrations_bySubtype$Cytokine == "IL.1RA", "IL-1RA")
Concentrations_bySubtype$Cytokine <- replace(Concentrations_bySubtype$Cytokine, Concentrations_bySubtype$Cytokine == "STNF.RI", "sTNF-RI")
Concentrations_bySubtype$Cytokine <- replace(Concentrations_bySubtype$Cytokine, Concentrations_bySubtype$Cytokine == "IL.12p40", "IL-12p40")
Concentrations_bySubtype$Cytokine <- replace(Concentrations_bySubtype$Cytokine, Concentrations_bySubtype$Cytokine == "IFNb", "IFNβ")
Concentrations_bySubtype$Cytokine <- replace(Concentrations_bySubtype$Cytokine, Concentrations_bySubtype$Cytokine == "IL.1b", "IL-1β")
Concentrations_bySubtype$Cytokine <- replace(Concentrations_bySubtype$Cytokine, Concentrations_bySubtype$Cytokine == "IL.33", "IL-33")


print(Concentrations_bySubtype)

#make heatmap
library(tibble)
library(tidyverse)
library(ggplot2)
library(MEM)
library(heatmaply)

Concentrations_bySubtype_hm <- Concentrations_bySubtype %>% 
  column_to_rownames(var = "Cytokine")

##################################################
##### standardising independently ###############
Concentrations_bySubtype_hm <- Concentrations_bySubtype_hm %>% select(1:3)%>% as.matrix()

#Concentrations_bySubtype <- asNumericMatrix(Concentrations_bySubtype)
Concentrations_bySubtype_hm <- Concentrations_bySubtype_hm %>%(log2)
#test
conc_hm <- heatmap(Concentrations_bySubtype_hm, 
                   show_heatmap_legend = TRUE) 
#colorRampPalette(c("blue", "white", "red" ))(1000 ) ,)

Concentrations_bySubtype_hm[1:3] <- scale(Concentrations_bySubtype_hm[1:3])

############################################################
##### standardising against a baseline column ##########v

# Extract the reference column (column 3)
reference_column <- Concentrations_bySubtype_hm[[3]]
# Compute the mean and standard deviation of the reference column
ref_mean <- mean(reference_column, na.rm = TRUE)
ref_sd <- sd(reference_column, na.rm = TRUE)

# Normalize columns 1:3 using the reference column's mean and standard deviation
Concentrations_bySubtype_hm[1:3] <- apply(Concentrations_bySubtype_hm[1:3], 2, function(x) (x - ref_mean) / ref_sd)

############################################################

str(Concentrations_bySubtype_hm)

#sort out colour scheme
col_order <- Concentrations_bySubtype%>%
  select(Cytokine, `Subtype 1`, `Subtype 2`, `Subtype 3`)


library(pheatmap)
library(RColorBrewer)


#conc_long <- col_order%>% 
#  pivot_longer(cols = -Cytokine, names_to = "Subtype", values_to = "Value")


#hm <- ggplot(conc_long, aes(x = Subtype, y= Cytokine, fill = Value))+
#  geom_tile()+
#  scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
#  theme_minimal()+
#  labs(title = "Concentration heatmap", x = "Subtype", y = "Cytokine", fill= "Value")


hm <- pheatmap(Concentrations_bySubtype_hm,
               scale = "row",#name = "Unknown",
               Colv = col_order, 
               angle_col = 0,
               fontsize_row = 18,
               fontsize_col = 18)

png("FIGURES/cytokine_concentration_heatmap.png",width=6,height=10,units="in",res=1200)
hm
dev.off()

#legend = TRUE)
#col = col_fun,
#cluster_columns = TRUE,
#cluster_rows = H.fit,
#show_heatmap_legend = TRUE,
#clustering_distance_rows = "euclidean",
# clustering_method_rows = "ward",
#bottom_annotation = columnAnnotation(mark=anno),
# width = ncol(df4)*unit(6, "mm"), 
# height = nrow(df4)*unit(4, "mm"),
#column_title_gp = gpar(fontsize = 3, fontface = "bold"),
#border = TRUE,
#row_km = 5,
#column_km = 3,
#show_column_names = TRUE,
#heatmap_legend_param = list(title = "Scale"),
#column_names_gp = grid::gpar(fontsize = 14),
#row_names_gp = grid::gpar(fontsize = 1),
#right_annotation = ha,
#row_dend_width = unit(40, "mm"),
#column_dend_side = c("top", "bottom"),
#column_dend_height = unit(10, "mm"))#+


