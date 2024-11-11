##LOVIT
#Figure 3
### make overall boxplots of each cytokine split by treatment

#load packages
library(ggplot2)
library(dplyr)
library(tidyverse)
library(tidyr)
###############################################################
###############################################################
# load dataframe - cytokine_df_LOVIT

# need to show proportions of each of each cluster
# Stacked + percent

#####################################
######### stacked proportions #######
#####################################

##IF NEEDED LOAD FILE 
Full_df <- read_csv("CSV files/full_df_MM/Complete_df.csv")
##### relabel cluster columns 

colnames(Full_df)
#counts_percofsample_spread <- counts_percofsample_spread[ -c(1) ]

#select the columns of interest
selected_cyto <- Full_df %>% dplyr::select(c(1,14,19,26,27,21, 32, 37,66))

colnames(selected_cyto)
head(selected_cyto)


colnames(selected_cyto)[c(3, 4, 5, 6, 7)] <- c("IL-6", "IL-10", "IL-1RA", "IL-1⍺", "IL-33")

# Reshape data to long format
df_long <- selected_cyto %>%
  pivot_longer(cols = c("CCL20", "IL-6", "IL-10", "IL-1RA", "IL-1⍺", "IL-33"),
               names_to = "Cytokine",
               values_to = "Value")

# Change names of timepoints and treatments
df_long <- df_long %>%
  mutate(
    timepoint = recode(timepoint, "S1" = "Day 1", "S3" = "Day 7"),
    Treatment = recode(Treatment, "A" = "Vitamin C", "B" = "Placebo"))

# Calculate log2 fold change relative to the first timepoint within each Treatment and Cytokine
df_log2_fc <- df_long %>%
  group_by(Treatment, Cytokine) %>%
  mutate(FoldChange = Value / first(Value[timepoint == min(timepoint)]),
         Log2FoldChange = log2(FoldChange)) %>%
  ungroup()


#################################################################################
# Create a new Treatment_Timepoint factor that combines Treatment and Timepoint
df_log2_fc <- df_log2_fc %>%
  mutate(Treatment_Timepoint = interaction(Treatment, timepoint, sep = " ")) %>%
  # Order levels to keep Placebo and Vitamin C together but within each timepoint
  mutate(Treatment_Timepoint = factor(Treatment_Timepoint, 
                                      levels = c("Placebo Day 1", "Placebo Day 7", 
                                                 "Vitamin C Day 1", "Vitamin C Day 7")))

# Set the custom order for cytokines
custom_cytokine_order <- c("IL-6", "IL-10", "IL-1RA", "CCL20", "IL-1⍺", "IL-33")


# Modify the Cytokine factor in your data frame to use the custom order
df_log2_fc <- df_log2_fc %>%
  mutate(Cytokine = factor(Cytokine, levels = custom_cytokine_order))


#################################################################################
# Plot with the new Treatment_Timepoint variable
overall_boxplot <- ggplot(df_log2_fc, aes(x = Treatment_Timepoint, y = Log2FoldChange, fill = Treatment))+
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA, alpha = 0.6) +  # Set alpha for boxplot fill
  geom_jitter(position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.8), 
              size = 1.5, alpha = 0.4) +  # Add individual points with more jitter
  facet_grid(. ~ Cytokine, scales = "free_y") +  # Facet cytokines as columns with labels at the top
  labs(x = "Timepoint", y = "Log2 Fold Change", fill = "Treatment") +
  scale_fill_manual(values = c("Placebo" = "gold1", "Vitamin C" = "orchid1")) +  # Set custom colors
  scale_x_discrete(labels = c("Placebo Day 1" = "Day 1", "Placebo Day 7" = "Day 7", 
                              "Vitamin C Day 1" = "Day 1", "Vitamin C Day 7" = "Day 7")) +  # Custom x-axis labels (only Day 1, Day 7)
  theme_classic() +
  theme(
    strip.text = element_text(size = 15, face = "bold"),  # Format cytokine labels at the top
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1), # Rotate x-axis labels for readability
    axis.title.x = element_text(size = 15, face = "bold"),
    axis.title.y = element_text(size = 15, face = "bold"),
    legend.text = element_text(size = 12),  # Increase legend text size
    legend.title = element_text(size = 14)  # Increase legend title size
  )

overall_boxplot


png("Figures/Figure_3/overall_boxplots.png", width=15,height=7,units="in",res=1200)
overall_boxplot
dev.off()


##################################################################################################################################################################
#################################################################################
