### LOVIT

###PCA plot 
###Bioplot
## PCA centroid 


################################################################################
################################################################################

# exploring the use of the longitudinal analysis significant gene set in PCA analysis. 
# use df that was created already. 

# 1) PCA
#runPCA
########################### 
#####################

#read df if not in environmment already
cytokine_df_LOVIT <- read.csv("CSV files/patient_df_log10_3clusters.csv")

##need to make an adjacency matrix
head(cytokine_df_LOVIT)

cytokine_df_LOVIT <- cytokine_df_LOVIT[, -c(1, 28:62)]

head (cytokine_df_LOVIT)
colnames(cytokine_df_LOVIT)
str(cytokine_df_LOVIT)

cytokine_df_log10<- cytokine_df_LOVIT %>% mutate_at(vars(2:27), log10) %>% 
  select(2:27, groups) %>% as.data.frame() %>% na.omit()

cytokine_df_log10<- cytokine_df_LOVIT %>%
  select(1:26, groups) %>% as.data.frame() %>% na.omit()

pca <- prcomp(cytokine_df_log10[1:26], scale. = F)
pca1 <- as.data.frame(pca$x)
pca1

cytokine_df_LOVIT_PCA <- cbind(cytokine_df_LOVIT, pca1[1:4])
cytokine_df_LOVIT_PCA$groups

colnames(cytokine_df_LOVIT_PCA)

################################################################################
################################################################################

## PLOT PCA

#################### see what it looks like with PCA  ##################################
####### run PCA #### 
library(NbClust)
library(factoextra)


###### PCA plot #####
str(cytokine_pca1)
#cytokine_pca$sample_id<- as.numeric(cytokine_pca$sample_id)
cytokine_pca1 <- cytokine_pca1[-c(1)]
#colnames(cytokine_pca) <- make.unique(names(cytokine_pca))

# subtype colours
#"1" = "#1BD76E", "2" = "#D76E1B", "3" = "purple4"

PCA_plot <- ggplot(cytokine_df_LOVIT_PCA, aes(x = PC1, y = PC2)) +
  geom_point(aes(fill = as.factor(groups)),pch =21, alpha = 1, size = 3, colour = "black")+
  theme_bw()+ 
  guides(fill = guide_legend(override.aes = list(shape=21)))+
  #scale_fill_manual(values = c("1" = "#1BD76E", "2" = "#D76E1B", "3" = "purple4"))+
  #scale_fill_manual(values = c("A" = "orchid1", "B" = "gold1"))+
  #scale_fill_manual(values = c("Yes" = "firebrick2", "No" = "dodgerblue"))+
  #scale_fill_gradient(low="blue", high="red")+
  #scale_fill_viridis_c()+
  scale_fill_manual(values = c("1" = "#1BD76E", 
                               "2" = "#D76E1B", 
                               "3" = "purple4"),
                    name = "Subtype")+
  scale_shape_manual(values = c(21, 22))+
  theme(aspect.ratio = 1)
#guides(fill=guide_legend(title="Cluster"))
PCA_plot

ggsave(plot= PCA_plot, "FIGURES/Figure_1/PCA.png", width = 7, height = 7, units = "in", dpi = 1200)


#######################
######################
## biplot

library(factoextra)

biplot <- fviz_pca_biplot(pca, repel = TRUE,
                          select.var = list(contrib = 10), #select number of variables - top 10 contributing factors
                          col.var = "contrib", # Variables color
                          gradient.cols = c("grey", "blue", "red"), #gradient colour scheme
                          col.ind = "grey",  # Individuals color
                          alpha.ind = 0.3,
                          label = "var",
                          title = "",
                          pointsize = 3,labelsize = 5)+
  theme_bw()+
  theme(aspect.ratio = 1)+
  #ylab("PC2")+
  #xlab("PC1")+
  theme(text = element_text(size=16),
        legend.text=element_text(size=16))

biplot

ggsave(plot= biplot, "Figures/Figure_1/biplot.png", width = 7, height = 7, units = "in", dpi = 1200)
