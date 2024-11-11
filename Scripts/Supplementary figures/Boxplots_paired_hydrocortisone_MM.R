### make stacked boxplots of cytokines grouped by cluster

#load packages
library(ggplot2)
library(dplyr)
library(tidyverse)
###############################################################
###############################################################
# load dataframe - cytokine_df_LOVIT

#####################################
######### stacked proportions #######
#####################################

#########################

#### systemic corticosteroids vs none

### compare the changes in biomarkers over time by concomitant hydrocortisone therapy 

##IF NEEDED LOAD FILE 
Full_df <- read_csv("CSV files/full_df_MM/Complete_df.csv")
##### relabel cluster columns 

colnames(Full_df)
#select the columns of interest
steroids <- Full_df %>% dplyr::select(c(1,5:17,19:32,36:37, 50, 66))

head(steroids)
colnames(steroids)

#remove VEGF from df
steroids <- steroids[, -c(21)] 

##filter the treatment
##########filter the steroids ######
steroids_Y <- steroids %>% dplyr::filter(Systemic_corticoster=="Yes")
steroids_N <- steroids %>% dplyr::filter(Systemic_corticoster=="No")

head(steroids_Y)
head(steroids_N)

colnames(steroids_Y)
colnames(steroids_N)


#############################################################
## NO Corticosteroid

# Reshape data to long format
df_long_N <- steroids_N %>%
  pivot_longer(cols = c(2:27),
               names_to = "Cytokine",
               values_to = "Value")

# Change names of timepoints and treatments
df_long_N <- df_long_N %>%
  mutate(
    timepoint = recode(timepoint, "S1" = "Day 1", "S3" = "Day 7"),
    Treatment = recode(Treatment, "A" = "Vitamin C", "B" = "Placebo"))

# Calculate log2 fold change relative to the first timepoint within each Treatment and Cytokine
df_log2_fc <- df_long_N %>%
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


#################################################################################
# Plot with the new Treatment_Timepoint variable

overall_boxplot <- ggplot(df_log2_fc, aes(x = Treatment_Timepoint, y = Log2FoldChange)) +
  geom_line(aes(group = patient_id), color = "black", alpha = 0.5, size = 0.2) +
  geom_boxplot(aes(fill = Treatment), position = position_dodge(width = 0.8), outlier.shape = NA, alpha = 0.6) +  # Boxplot with transparency
  geom_jitter(aes(color = factor(timepoint)),  # Color points by Timepoint factor
              position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.8),
              shape = 21, size = 1.5, stroke = 0.5, alpha = 0.9) +  # Use shape 21 with black border and thin stroke
  facet_wrap(~ Cytokine, scales = "free_y") +  # Split cytokines into individual panels
  labs(x = "Timepoint", y = "Log2 Fold Change", fill = "Treatment", color = "Timepoint") +  # Add color legend label
  scale_fill_manual(values = c("Placebo" = "gold1", "Vitamin C" = "orchid1")) +  # Custom colors for boxplot fill
  scale_color_manual(values = c("Day 1" = "black", "Day 7" = "grey")) +  # Custom colors for points by Timepoint
  scale_x_discrete(labels = c("Placebo Day 1" = "Day 1", "Placebo Day 7" = "Day 7",
                              "Vitamin C Day 1" = "Day 1", "Vitamin C Day 7" = "Day 7")) +  # Custom x-axis labels
  expand_limits(y = range(df_log2_fc$Log2FoldChange, na.rm = TRUE)) +  # Set y-axis to fit all points
  theme_classic() +
  theme(
    strip.text = element_text(size = 8, face = "bold"),  # Format cytokine panel labels
    axis.text.x = element_blank(),  # Remove x-axis tick labels
    axis.title.x = element_text(size = 15, face = "bold"),
    axis.title.y = element_text(size = 15, face = "bold"),
    legend.text = element_text(size = 12),  # Increase legend text size
    legend.title = element_text(size = 14)  # Increase legend title size
  )

overall_boxplot



png("Figures/Supplementary/Boxplots_no_corticosteroids.png", width=15,height=7,units="in",res=1200)
overall_boxplot
dev.off()
###########################################################################

#############################################################
## Corticosteroid

colnames(steroids_Y)

# Reshape data to long format
df_long_Y <- steroids_Y %>%
  pivot_longer(cols = c(2:27),
               names_to = "Cytokine",
               values_to = "Value")

head(df_long_Y)

# Change names of timepoints and treatments
df_long_Y <- df_long_Y %>%
  mutate(
    timepoint = recode(timepoint, "S1" = "Day 1", "S3" = "Day 7"),
    Treatment = recode(Treatment, "A" = "Vitamin C", "B" = "Placebo"))


# Calculate log2 fold change relative to the first timepoint within each Treatment and Cytokine
df_log2_fc <- df_long_Y %>%
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


#################################################################################
# Plot with the new Treatment_Timepoint variable

overall_boxplot <- ggplot(df_log2_fc, aes(x = Treatment_Timepoint, y = Log2FoldChange)) +
  geom_line(aes(group = patient_id), color = "black", alpha = 0.5, size = 0.2) +
  geom_boxplot(aes(fill = Treatment), position = position_dodge(width = 0.8), outlier.shape = NA, alpha = 0.6) +  # Boxplot with transparency
  geom_jitter(aes(color = factor(timepoint)),  # Color points by Timepoint factor
              position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.8),
              shape = 21, size = 1.5, stroke = 0.5, alpha = 0.9) +  # Use shape 21 with black border and thin stroke
  facet_wrap(~ Cytokine, scales = "free_y") +  # Split cytokines into individual panels
  labs(x = "Timepoint", y = "Log2 Fold Change", fill = "Treatment", color = "Timepoint") +  # Add color legend label
  scale_fill_manual(values = c("Placebo" = "gold1", "Vitamin C" = "orchid1")) +  # Custom colors for boxplot fill
  scale_color_manual(values = c("Day 1" = "black", "Day 7" = "grey")) +  # Custom colors for points by Timepoint
  scale_x_discrete(labels = c("Placebo Day 1" = "Day 1", "Placebo Day 7" = "Day 7",
                              "Vitamin C Day 1" = "Day 1", "Vitamin C Day 7" = "Day 7")) +  # Custom x-axis labels
  expand_limits(y = range(df_log2_fc$Log2FoldChange, na.rm = TRUE)) +  # Set y-axis to fit all points
  theme_classic() +
  theme(
    strip.text = element_text(size = 8, face = "bold"),  # Format cytokine panel labels
    axis.text.x = element_blank(),  # Remove x-axis tick labels
    axis.title.x = element_text(size = 15, face = "bold"),
    axis.title.y = element_text(size = 15, face = "bold"),
    legend.text = element_text(size = 12),  # Increase legend text size
    legend.title = element_text(size = 14)  # Increase legend title size
  )

overall_boxplot



png("Figures/Supplementary/Boxplots_corticosteroids.png", width=15,height=7,units="in",res=1200)
overall_boxplot
dev.off()
###########################################################################

