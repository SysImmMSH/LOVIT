## LOVIT 

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
cytokine_df_LOVIT <- read.csv("CSV files/patient_df_3clusters.csv")

##need to make an adjacency matrix
head(cytokine_df_LOVIT)

cytokine_df_LOVIT <- cytokine_df_LOVIT[, -c(2, 3,4,31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 
                                            45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59,
                                            60, 61, 62, 63)]

head (cytokine_df_LOVIT)
colnames(cytokine_df_LOVIT)
str(cytokine_df_LOVIT)

cytokine_df_log10<- cytokine_df_LOVIT %>% mutate_at(vars(2:27), log10) %>% 
  select(2:27, groups) %>% as.data.frame() %>% na.omit()


pca <- prcomp(cytokine_df_log10[2:27], scale. = F)
pca1 <- as.data.frame(pca$x)
pca1

cytokine_df_LOVIT_PCA <- cbind(cytokine_df_LOVIT, pca1[1:4])
cytokine_df_LOVIT_PCA$groups

##########################################################################################


# calculate centroids of Cohort_time in the PCA

colnames(cytokine_df_LOVIT_PCA)
PCA_centroids <- aggregate(cytokine_df_LOVIT_PCA[,30:31], list(Type = cytokine_df_LOVIT_PCA$groups), median)

PCA_centroids$Type <- factor(PCA_centroids$Type, levels = c("1", "2", "3"))



########################################################################################################
PCA_centroids_plot <- ggplot(PCA_centroids, aes(x = PC1, y = PC2))+
  geom_density_2d(data = cytokine_df_LOVIT_PCA, aes(x = PC1, y = PC2, color = as.factor(groups)), 
                  colour = "grey", alpha = 0.5)+
  geom_point(aes(fill = as.factor(Type)), alpha = 0.8, size = 5, stroke = 0.7)+
  #define fill colours
  scale_fill_manual(values = c(" 1" = "#efe350ff",
                               " 2" = "#f68f46ff",
                               " 3" = "#a65c85ff"))+
  #scale_shape_manual(values = c("Sepsis" = 21, "Cardiac" = 24))+
  #theme settings
  theme_bw()+ 
  guides(colour =guide_legend(override.aes=list(shape=21)))+ #legend for fill
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(legend.position = "bottom")+ #sets the legend at the bottom
  #ylim(c(-20,25))+
  #xlim(c(-40,45))+
  theme(aspect.ratio = 1)+ #keeps the aspect ratio square
  ylab("PC2")+
  xlab("PC1")+ theme(text = element_text(size = 12))+    
  # Add paths with arrows for each Type
  geom_path(data = PCA_centroids %>% filter(Type == 1), aes(group = Type), 
            color = "green", 
            arrow = arrow(type = "closed", length = unit(0.15, "inches"), ends = "last")) +
  
  geom_path(data = PCA_centroids %>% filter(Type == 2), aes(group = Type), 
            color = "black", 
            arrow = arrow(type = "closed", length = unit(0.15, "inches"), ends = "last")) +
  
  geom_path(data = PCA_centroids %>% filter(Type == 3), aes(group = Type), 
            color = "black", 
            arrow = arrow(type = "closed", length = unit(0.15, "inches"), ends = "last"))

  
PCA_centroids_plot


summary(cytokine_df_LOVIT_PCA$PC1)
summary(cytokine_df_LOVIT_PCA$PC2)

# Check for NA or Inf values in PC1 and PC2
sum(is.na(cytokine_df_LOVIT_PCA$PC1))  # Count of NAs in PC1
sum(is.na(cytokine_df_LOVIT_PCA$PC2))  # Count of NAs in PC2

sum(is.infinite(cytokine_df_LOVIT_PCA$PC1))  # Count of Inf values in PC1
sum(is.infinite(cytokine_df_LOVIT_PCA$PC2))  # Count of Inf values in PC2

  
########################################################################################################
PCA_centroids_plot <- ggplot(PCA_centroids, aes(x = PC1, y = PC2))+
  geom_point(data = cytokine_df_LOVIT_PCA, aes(x = PC1, y = PC2, color = as.factor(groups)), 
                  colour = "pink", alpha = 0.5)+ #replaces density with scatter plot
  geom_point(aes(fill = as.factor(Type)), colour = "black", alpha = 0.8, size = 5, stroke = 0.7)+
  #define fill colours
  scale_fill_manual(values = c(" 1" = "#efe350ff",
                               " 2" = "#f68f46ff",
                               " 3" = "#a65c85ff"))+
  #scale_shape_manual(values = c("Sepsis" = 21, "Cardiac" = 24))+
  #theme settings
  theme_bw()+ 
  guides(fill=guide_legend(override.aes=list(shape=21)))+ #legend for fill
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(legend.position = "bottom")+ #sets the legend at the bottom
  #ylim(c(-20,25))+
  #xlim(c(-40,45))+
  theme(aspect.ratio = 1)+ #keeps the aspect ratio square
  ylab("PC2")+
  xlab("PC1")+ theme(text = element_text(size = 12))+    
  # Add paths with arrows for each Type
  geom_path(data = PCA_centroids %>% filter(Type == 1), aes(group = Type), 
            color = "black", 
            arrow = arrow(type = "closed", length = unit(0.15, "inches"), ends = "last")) +
  
  geom_path(data = PCA_centroids %>% filter(Type == 2), aes(group = Type), 
            color = "black", 
            arrow = arrow(type = "closed", length = unit(0.15, "inches"), ends = "last")) +
  
  geom_path(data = PCA_centroids %>% filter(Type == 3), aes(group = Type), 
            color = "black", 
            arrow = arrow(type = "closed", length = unit(0.15, "inches"), ends = "last"))


PCA_centroids_plot


########################################################################################################

PCA_centroids_plot <- ggplot(PCA_centroids, aes(x = PC1, y = PC2)) +
  # Add 2D density plot for background data
  geom_density_2d(data = cytokine_df_LOVIT_PCA, aes(x = PC1, y = PC2, color = as.factor(groups)), 
                  colour = "grey", alpha = 0.5) +
  
  # Scatter plot for centroids, coloring by 'Type' (use color, not fill)
  geom_point(aes(color = as.factor(Type)), alpha = 0.8, size = 5, stroke = 0.7) +
  
  # Define custom colors for each 'Type'
  scale_color_manual(values = c("1" = "#289e65", 
                                "2" = "#d97416", 
                                "3" = "#4e08a3")) +
  
  # Theme and legend adjustments
  theme_bw() +  
  guides(color = guide_legend(override.aes = list(shape = 21))) +  # Adjust legend for color
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(legend.position = "bottom") +  # Legend position
  theme(aspect.ratio = 1) +  # Keep aspect ratio square
  ylab("PC2") + 
  xlab("PC1") + 
  theme(text = element_text(size = 12)) +
  
  # Add paths with arrows for each Type
  geom_path(data = PCA_centroids %>% filter(Type == 1), aes(group = Type), 
            color = "black", 
            arrow = arrow(type = "closed", length = unit(0.15, "inches"), ends = "last")) +
  
  geom_path(data = PCA_centroids %>% filter(Type == 2), aes(group = Type), 
            color = "black", 
            arrow = arrow(type = "closed", length = unit(0.15, "inches"), ends = "last")) +
  
  geom_path(data = PCA_centroids %>% filter(Type == 3), aes(group = Type), 
            color = "black", 
            arrow = arrow(type = "closed", length = unit(0.15, "inches"), ends = "last"))

# Display the plot
PCA_centroids_plot

ggsave(plot = PCA_centroids_plot, "FIGURES/PCA_centroids_plot.png", dpi=300, height = 5, width = 5, units = "in" )

