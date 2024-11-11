##LOVIT
#Figure 3

## 1. Volcano plots - Overall VitC vs Placebo
## 2. Volcano plots - Systemic corticosteroids vs none a) Placebo
#                                                      b) Vitamin C

##load files, 

all_cytokines <- read.csv("CSV files/OneDrive_1_04-11-2024/all_cytokines.csv")

clinical_data_S1 <- read_excel("CSV files/LOVIT_MSH_EXPORT_DATA.xlsx")

All_timepoints <- read.csv("CSV files/OneDrive_1_04-11-2024/LOVIT_finalmerge.csv")

samples <- read.csv("CSV files/patient_df_3clusters.csv")

names(samples)
head(samples)

library(readxl)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(ggrepel)

#merge cytokines with the clinical data
Full_df <- merge(all_cytokines, All_timepoints, by.x = "sample_ID", by.y = "sample_ID") 

colnames(Full_df)

#select the columns of interest
Day7 <- Full_df %>% dplyr::select(c(1,5:17,19:32,36:37, 66))

colnames(Day7)

#remove VEGF from df
Day7 <- Day7[, -c(21)] 

names(Day7)
head(Day7)

##filter the treatment
#Day7 <- Day7 %>% dplyr::filter(timepoint=="S3")
Day7_Placebo <- Day7 %>% dplyr::filter(Treatment=="B")
Day7_VitC <- Day7 %>% dplyr::filter(Treatment=="A")

colnames(Day7_Placebo)
colnames(Day7_VitC)

##log2 the cytokines
df_log_Plac <- Day7_Placebo %>% mutate_at(vars(2:27), log2) %>% 
  select(2:27, sample_ID, timepoint) %>% as.data.frame() %>% 
  na.omit()
colnames(df_log_Plac)

df_log_VitC <- Day7_VitC %>% mutate_at(vars(2:27), log2) %>% 
  select(2:27, sample_ID, timepoint) %>% as.data.frame() %>% 
  na.omit()
colnames(df_log_VitC)


names(df_log_Plac)
analysis_Plac <- df_log_Plac %>% select("sample_ID", 1:26) 
colnames(analysis_Plac)

names(df_log_VitC)
analysis_VitC <- df_log_VitC %>% select("sample_ID", 1:26) 
colnames(analysis_VitC)


################################################################################
#differential cytokine plot

#make matrix of values samples columns and frequencies as rows
count_comparison_tmydf = setNames(data.frame(t(analysis_Plac[,-1])), analysis_Plac[,1]) %>% as.matrix() 
str(count_comparison_tmydf)

dim(count_comparison_tmydf)


library(limma)

f <- factor(df_log_Plac$timepoint, levels=c("S1", "S3"))
design <- model.matrix(~0+f)
dim(design)
colnames(design) <- c("S1","S3")
#To make all pair-wise comparisons between the three groups one could proceed
fit <- lmFit(count_comparison_tmydf, design)
contrast.matrix <- makeContrasts(S3-S1,
                                 levels=design)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
topTable(fit2)
?eBayes

df_stats_Plac <- topTable(fit2, coef=1, adjust="BH", number = "inf") %>% rownames_to_column("Cytokine")

write_csv(df_stats_Plac, "CSV files/Output_files/Limma/Stats/Placebo_D7vsD1.csv")
# add a column of NAs
df_stats_Plac$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
df_stats_Plac$diffexpressed[df_stats_Plac$logFC > 0.5 & df_stats_Plac$adj.P.Val < 0.01] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
df_stats_Plac$diffexpressed[df_stats_Plac$logFC < -0.5 & df_stats_Plac$adj.P.Val < 0.01] <- "DOWN"

# Create a new column "delabel" to de, that will contain the name of genes differentially expressed (NA in case they are not)
df_stats_Plac$delabel <- NA
df_stats_Plac$delabel[df_stats_Plac$diffexpressed != "NO"] <- df_stats_Plac$Cytokine[df_stats_Plac$diffexpressed != "NO"]

head(df_stats_Plac)

df_stats_Plac<- df_stats_Plac%>%
  mutate(delabel = recode(delabel, 
                          "IFNy" = "IFN𝛄", 
                          "TNFa" = "TNF⍺", 
                          "IL.6" = "IL-6", 
                          "IL.5" = "IL-5", 
                          "IL.1a" = "IL-1⍺", 
                          "IL.18" = "IL-18",
                          "IL.4" = "IL-4", 
                          "IFNa2" = "IFN⍺2", 
                          "IL.10" = "IL-10",
                          "IL.1RA" = "IL-1RA", 
                          "STNF.RI" = "sTNF-RI", 
                          "IL.12p40" = "IL-12p40",
                          "IFNb" = "IFNβ", 
                          "IL.1b" = "IL-1β", 
                          "IL.33" = "IL-33"))


mycolors <- c("#709eeb", "#b02428", "grey")
names(mycolors) <- c("DOWN", "UP", "NO")


Time_D7vsD1 <- ggplot(data=df_stats_Plac, aes(x=logFC, y=-log10(adj.P.Val), label=delabel)) + 
  geom_point(aes(fill = diffexpressed), pch = 21, colour = "black", size = 5)+
  theme_classic()+
  geom_vline(xintercept=c(-0.5, 0.5), col="black", linetype = "dashed") +
  geom_hline(yintercept=-log10(0.01), col="black", linetype = "dashed")+
  scale_fill_manual(values=mycolors)+
  theme(axis.text.x = element_text(size=18,  colour = "black"),
        axis.text.y = element_text(size=18, colour = "black"))+
  theme(axis.title.x = element_text(size=19, face= "bold", colour = "black"))+
  theme(axis.title.y = element_text(size=19, face= "bold", colour = "black"))+
  theme(legend.text=element_text(size=14))+
  theme(legend.title=element_text(size=15))+
  ylim(0,30)+
  xlim(-5.5,5.5)+
  geom_text_repel(size = 5)+
  ggtitle("Placebo - D7 vs D1")

Time_D7vsD1

ggsave(plot = Time_D7vsD1, "Figures/Figure_3/limmaplot_Placebo_D7vsD1.png", dpi = 1200,  height = 7,  width = 8)


#########################
################################################################################
#differential cytokine plot - Vitamin C

#make matrix of values samples columns and frequencies as rows
count_comparison_tmydf = setNames(data.frame(t(analysis_VitC[,-1])), analysis_VitC[,1]) %>% as.matrix() 
str(count_comparison_tmydf)

dim(count_comparison_tmydf)


library(limma)

f <- factor(df_log_VitC$timepoint, levels=c("S1", "S3"))
design <- model.matrix(~0+f)
dim(design)
colnames(design) <- c("S1","S3")
#To make all pair-wise comparisons between the three groups one could proceed
fit <- lmFit(count_comparison_tmydf, design)
contrast.matrix <- makeContrasts(S3-S1,
                                 levels=design)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
topTable(fit2)
?eBayes

df_stats_VitC <- topTable(fit2, coef=1, adjust="BH", number = "inf") %>% rownames_to_column("Cytokine")

write_csv(df_stats_VitC, "CSV files/Output_files/Limma/Stats/VitaminC_D7vsD1.csv")
# add a column of NAs
df_stats_VitC$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
df_stats_VitC$diffexpressed[df_stats_VitC$logFC > 0.5 & df_stats_VitC$adj.P.Val < 0.01] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
df_stats_VitC$diffexpressed[df_stats_VitC$logFC < -0.5 & df_stats_VitC$adj.P.Val < 0.01] <- "DOWN"

# Create a new column "delabel" to de, that will contain the name of genes differentially expressed (NA in case they are not)
df_stats_VitC$delabel <- NA
df_stats_VitC$delabel[df_stats_VitC$diffexpressed != "NO"] <- df_stats_VitC$Cytokine[df_stats_VitC$diffexpressed != "NO"]

head(df_stats_VitC)

df_stats_VitC<- df_stats_VitC%>%
  mutate(delabel = recode(delabel, 
                          "IFNy" = "IFN𝛄", 
                          "TNFa" = "TNF⍺", 
                          "IL.6" = "IL-6", 
                          "IL.5" = "IL-5", 
                          "IL.1a" = "IL-1⍺", 
                          "IL.18" = "IL-18",
                          "IL.4" = "IL-4", 
                          "IFNa2" = "IFN⍺2", 
                          "IL.10" = "IL-10",
                          "IL.1RA" = "IL-1RA", 
                          "STNF.RI" = "sTNF-RI", 
                          "IL.12p40" = "IL-12p40",
                          "IFNb" = "IFNβ", 
                          "IL.1b" = "IL-1β", 
                          "IL.33" = "IL-33"))


mycolors <- c("#709eeb", "#b02428", "grey")
names(mycolors) <- c("DOWN", "UP", "NO")


Time_D7vsD1 <- ggplot(data=df_stats_VitC, aes(x=logFC, y=-log10(adj.P.Val), label=delabel)) + 
  geom_point(aes(fill = diffexpressed), pch = 21, colour = "black", size = 5)+
  theme_classic()+
  geom_vline(xintercept=c(-0.5, 0.5), col="black", linetype = "dashed") +
  geom_hline(yintercept=-log10(0.01), col="black", linetype = "dashed")+
  scale_fill_manual(values=mycolors)+
  theme(axis.text.x = element_text(size=18,  colour = "black"),
        axis.text.y = element_text(size=18, colour = "black"))+
  theme(axis.title.x = element_text(size=19, face= "bold", colour = "black"))+
  theme(axis.title.y = element_text(size=19, face= "bold", colour = "black"))+
  theme(legend.text=element_text(size=14))+
  theme(legend.title=element_text(size=15))+
  ylim(0,32)+
  xlim(-5.5,5.5)+
  geom_text_repel(size = 5)+
  ggtitle("VitaminC - D7 vs D1")

Time_D7vsD1

ggsave(plot = Time_D7vsD1, "Figures/Figure_3/limmaplot_VitC_D7vsD1.png", dpi = 1200,  height = 7,  width = 8)


#########################

#### systemic corticosteroids vs none

### compare the changes in biomarkers over time by concomitant hydrocortisone therapy 

colnames(Full_df)
#select the columns of interest
steroids <- Full_df %>% dplyr::select(c(1,5:17,19:32,36:37, 50, 66))

head(steroids)
#remove VEGF from df
steroids <- steroids[, -c(21)] 

##filter the treatment
steroids_Placebo <- steroids %>% dplyr::filter(Treatment=="B")
steroids_VitC <- steroids %>% dplyr::filter(Treatment=="A")

colnames(steroids_Placebo)
colnames(steroids_VitC)

##########filter the steroids ######
#Placebo
steroids_Plac_Y <- steroids_Placebo %>% dplyr::filter(Systemic_corticoster=="Yes")
steroids_Plac_N <- steroids_Placebo %>% dplyr::filter(Systemic_corticoster=="No")

head(steroids_Plac_Y)
head(steroids_Plac_N)

#VitaminC
steroids_VitC_Y <- steroids_VitC %>% dplyr::filter(Systemic_corticoster=="Yes")
steroids_VitC_N <- steroids_VitC %>% dplyr::filter(Systemic_corticoster=="No")

head(steroids_VitC_Y)
head(steroids_VitC_N)

##log2 the cytokines

df_log_Plac_Y <- steroids_Plac_Y %>% mutate_at(vars(2:27), log2) %>% 
  select(2:27, sample_ID, timepoint, Systemic_corticoster) %>% as.data.frame() %>% 
  na.omit()
colnames(df_log_Plac_Y)

df_log_Plac_N <- steroids_Plac_N %>% mutate_at(vars(2:27), log2) %>% 
  select(2:27, sample_ID, timepoint, Systemic_corticoster) %>% as.data.frame() %>% 
  na.omit()
colnames(df_log_Plac_N)

df_log_VitC_Y <- steroids_VitC_Y %>% mutate_at(vars(2:27), log2) %>% 
  select(2:27, sample_ID, timepoint, Systemic_corticoster) %>% as.data.frame() %>% 
  na.omit()
colnames(df_log_VitC_Y)

df_log_VitC_N <- steroids_VitC_N %>% mutate_at(vars(2:27), log2) %>% 
  select(2:27, sample_ID, timepoint, Systemic_corticoster) %>% as.data.frame() %>% 
  na.omit()
colnames(df_log_VitC_N)

# split_df<- df_log %>%
#   group_by(groups) %>%
#   group_split()
# 
# df_log_1 <- split_df[[1]]

names(df_log_Plac_Y)
analysis_Plac_Y <- df_log_Plac_Y %>% select("sample_ID", 1:26) 
colnames(analysis_Plac_Y)

names(df_log_Plac_N)
analysis_Plac_N <- df_log_Plac_N %>% select("sample_ID", 1:26) 
colnames(analysis_Plac_N)

names(df_log_VitC_Y)
analysis_VitC_Y <- df_log_VitC_Y %>% select("sample_ID", 1:26) 
colnames(analysis_VitC_Y)

names(df_log_VitC_N)
analysis_VitC_N <- df_log_VitC_N %>% select("sample_ID", 1:26) 
colnames(analysis_VitC_N)

################################################################################
################################################################################
################################################################################
#differential cytokine plot

#Placebo with corticosteroids
#make matrix of values samples columns and frequencies as rows
count_comparison_tmydf = setNames(data.frame(t(analysis_Plac_Y[,-1])), analysis_Plac_Y[,1]) %>% as.matrix() 
str(count_comparison_tmydf)

dim(count_comparison_tmydf)


library(limma)

f <- factor(df_log_Plac_Y$timepoint, levels=c("S1", "S3"))
design <- model.matrix(~0+f)
dim(design)
colnames(design) <- c("S1","S3")
#To make all pair-wise comparisons between the three groups one could proceed
fit <- lmFit(count_comparison_tmydf, design)
contrast.matrix <- makeContrasts(S3-S1,
                                 levels=design)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
topTable(fit2)
#?eBayes

df_stats_Plac_Y <- topTable(fit2, coef=1, adjust="BH", number = "inf") %>% rownames_to_column("Cytokine")

write_csv(df_stats_Plac_Y, "CSV files/Output_files/Limma/Stats/Placebo_D7vsD1_Cortico.csv")
# add a column of NAs
df_stats_Plac_Y$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
df_stats_Plac_Y$diffexpressed[df_stats_Plac_Y$logFC > 0.5 & df_stats_Plac_Y$adj.P.Val < 0.01] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
df_stats_Plac_Y$diffexpressed[df_stats_Plac_Y$logFC < -0.5 & df_stats_Plac_Y$adj.P.Val < 0.01] <- "DOWN"

# Create a new column "delabel" to de, that will contain the name of genes differentially expressed (NA in case they are not)
df_stats_Plac_Y$delabel <- NA
df_stats_Plac_Y$delabel[df_stats_Plac_Y$diffexpressed != "NO"] <- df_stats_Plac_Y$Cytokine[df_stats_Plac_Y$diffexpressed != "NO"]

head(df_stats_Plac_Y)

df_stats_Plac_Y<- df_stats_Plac_Y%>%
  mutate(delabel = recode(delabel, 
                          "IFNy" = "IFN𝛄", 
                          "TNFa" = "TNF⍺", 
                          "IL.6" = "IL-6", 
                          "IL.5" = "IL-5", 
                          "IL.1a" = "IL-1⍺", 
                          "IL.18" = "IL-18",
                          "IL.4" = "IL-4", 
                          "IFNa2" = "IFN⍺2", 
                          "IL.10" = "IL-10",
                          "IL.1RA" = "IL-1RA", 
                          "STNF.RI" = "sTNF-RI", 
                          "IL.12p40" = "IL-12p40",
                          "IFNb" = "IFNβ", 
                          "IL.1b" = "IL-1β", 
                          "IL.33" = "IL-33"))


mycolors <- c("#709eeb", "#b02428", "grey")
names(mycolors) <- c("DOWN", "UP", "NO")


Steroids_PlacY <- ggplot(data=df_stats_Plac_Y, aes(x=logFC, y=-log10(adj.P.Val), label=delabel)) + 
  geom_point(aes(fill = diffexpressed), pch = 21, colour = "black", size = 5)+
  theme_classic()+
  geom_vline(xintercept=c(-0.5, 0.5), col="black", linetype = "dashed") +
  geom_hline(yintercept=-log10(0.01), col="black", linetype = "dashed")+
  scale_fill_manual(values=mycolors)+
  theme(axis.text.x = element_text(size=18,  colour = "black"),
        axis.text.y = element_text(size=18, colour = "black"))+
  theme(axis.title.x = element_text(size=19, face= "bold", colour = "black"))+
  theme(axis.title.y = element_text(size=19, face= "bold", colour = "black"))+
  theme(legend.text=element_text(size=14))+
  theme(legend.title=element_text(size=15))+
  ylim(0,25)+
  xlim(-5.5,5.5)+
  geom_text_repel(size = 5)+
  ggtitle("Placebo - corticosteroids")

Steroids_PlacY

ggsave(plot = Steroids_PlacY, "Figures/Figure_3/limmaplot_Placebo_D7vsD1_cortico.png", dpi = 1200,  height = 7,  width = 8)

################################################################################
#differential cytokine plot

#Placebo without corticosteroids
#make matrix of values samples columns and frequencies as rows
count_comparison_tmydf = setNames(data.frame(t(analysis_Plac_N[,-1])), analysis_Plac_N[,1]) %>% as.matrix() 
str(count_comparison_tmydf)

dim(count_comparison_tmydf)


library(limma)

f <- factor(df_log_Plac_N$timepoint, levels=c("S1", "S3"))
design <- model.matrix(~0+f)
dim(design)
colnames(design) <- c("S1","S3")
#To make all pair-wise comparisons between the three groups one could proceed
fit <- lmFit(count_comparison_tmydf, design)
contrast.matrix <- makeContrasts(S3-S1,
                                 levels=design)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
topTable(fit2)
?eBayes

df_stats_Plac_N <- topTable(fit2, coef=1, adjust="BH", number = "inf") %>% rownames_to_column("Cytokine")

write_csv(df_stats_Plac_N, "CSV files/Output_files/Limma/Stats/Placebo_D7vsD1_NoCortico.csv")
# add a column of NAs
df_stats_Plac_N$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
df_stats_Plac_N$diffexpressed[df_stats_Plac_N$logFC > 1 & df_stats_Plac_N$adj.P.Val < 0.01] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
df_stats_Plac_N$diffexpressed[df_stats_Plac_N$logFC < -1 & df_stats_Plac_N$adj.P.Val < 0.01] <- "DOWN"

# Create a new column "delabel" to de, that will contain the name of genes differentially expressed (NA in case they are not)
df_stats_Plac_N$delabel <- NA
df_stats_Plac_N$delabel[df_stats_Plac_N$diffexpressed != "NO"] <- df_stats_Plac_N$Cytokine[df_stats_Plac_N$diffexpressed != "NO"]

head(df_stats_Plac_N)

df_stats_Plac_N<- df_stats_Plac_N%>%
  mutate(delabel = recode(delabel, 
                          "IFNy" = "IFN𝛄", 
                          "TNFa" = "TNF⍺", 
                          "IL.6" = "IL-6", 
                          "IL.5" = "IL-5", 
                          "IL.1a" = "IL-1⍺", 
                          "IL.18" = "IL-18",
                          "IL.4" = "IL-4", 
                          "IFNa2" = "IFN⍺2", 
                          "IL.10" = "IL-10",
                          "IL.1RA" = "IL-1RA", 
                          "STNF.RI" = "sTNF-RI", 
                          "IL.12p40" = "IL-12p40",
                          "IFNb" = "IFNβ", 
                          "IL.1b" = "IL-1β", 
                          "IL.33" = "IL-33"))


mycolors <- c("#709eeb", "#b02428", "grey")
names(mycolors) <- c("DOWN", "UP", "NO")


Steroids_PlacN <- ggplot(data=df_stats_Plac_N, aes(x=logFC, y=-log10(adj.P.Val), label=delabel)) + 
  geom_point(aes(fill = diffexpressed), pch = 21, colour = "black", size = 5)+
  theme_classic()+
  geom_vline(xintercept=c(-0.5, 0.5), col="black", linetype = "dashed") +
  geom_hline(yintercept=-log10(0.01), col="black", linetype = "dashed")+
  scale_fill_manual(values=mycolors)+
  theme(axis.text.x = element_text(size=18,  colour = "black"),
        axis.text.y = element_text(size=18, colour = "black"))+
  theme(axis.title.x = element_text(size=19, face= "bold", colour = "black"))+
  theme(axis.title.y = element_text(size=19, face= "bold", colour = "black"))+
  theme(legend.text=element_text(size=14))+
  theme(legend.title=element_text(size=15))+
  ylim(0,25)+
  xlim(-5.5,5.5)+
  geom_text_repel(size = 5)+
  ggtitle("Placebo - No corticosteroids")

Steroids_PlacN

ggsave(plot = Steroids_PlacN, "Figures/Figure_3/limmaplot_Placebo_D7vsD1_Nocortico.png", dpi = 1200,  height = 7,  width = 8)

################################################################################
################################################################################
#differential cytokine plot

#Vitamin C with corticosteroids
#make matrix of values samples columns and frequencies as rows
count_comparison_tmydf = setNames(data.frame(t(analysis_VitC_Y[,-1])), analysis_VitC_Y[,1]) %>% as.matrix() 
str(count_comparison_tmydf)

dim(count_comparison_tmydf)


library(limma)

f <- factor(df_log_VitC_Y$timepoint, levels=c("S1", "S3"))
design <- model.matrix(~0+f)
dim(design)
colnames(design) <- c("S1","S3")
#To make all pair-wise comparisons between the three groups one could proceed
fit <- lmFit(count_comparison_tmydf, design)
contrast.matrix <- makeContrasts(S3-S1,
                                 levels=design)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
topTable(fit2)
?eBayes

df_stats_VitC_Y <- topTable(fit2, coef=1, adjust="BH", number = "inf") %>% rownames_to_column("Cytokine")

write_csv(df_stats_VitC_Y, "CSV files/Output_files/Limma/Stats/VitC_D7vsD1_Cortico.csv")
# add a column of NAs
df_stats_VitC_Y$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
df_stats_VitC_Y$diffexpressed[df_stats_VitC_Y$logFC > 0.5 & df_stats_VitC_Y$adj.P.Val < 0.01] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
df_stats_VitC_Y$diffexpressed[df_stats_VitC_Y$logFC < -0.5 & df_stats_VitC_Y$adj.P.Val < 0.01] <- "DOWN"

# Create a new column "delabel" to de, that will contain the name of genes differentially expressed (NA in case they are not)
df_stats_VitC_Y$delabel <- NA
df_stats_VitC_Y$delabel[df_stats_VitC_Y$diffexpressed != "NO"] <- df_stats_VitC_Y$Cytokine[df_stats_VitC_Y$diffexpressed != "NO"]

head(df_stats_VitC_Y)

df_stats_VitC_Y<- df_stats_VitC_Y%>%
  mutate(delabel = recode(delabel, 
                          "IFNy" = "IFN𝛄", 
                          "TNFa" = "TNF⍺", 
                          "IL.6" = "IL-6", 
                          "IL.5" = "IL-5", 
                          "IL.1a" = "IL-1⍺", 
                          "IL.18" = "IL-18",
                          "IL.4" = "IL-4", 
                          "IFNa2" = "IFN⍺2", 
                          "IL.10" = "IL-10",
                          "IL.1RA" = "IL-1RA", 
                          "STNF.RI" = "sTNF-RI", 
                          "IL.12p40" = "IL-12p40",
                          "IFNb" = "IFNβ", 
                          "IL.1b" = "IL-1β", 
                          "IL.33" = "IL-33"))


mycolors <- c("#709eeb", "#b02428", "grey")
names(mycolors) <- c("DOWN", "UP", "NO")


Steroids_VitCY <- ggplot(data=df_stats_VitC_Y, aes(x=logFC, y=-log10(adj.P.Val), label=delabel)) + 
  geom_point(aes(fill = diffexpressed), pch = 21, colour = "black", size = 5)+
  theme_classic()+
  geom_vline(xintercept=c(-0.5, 0.5), col="black", linetype = "dashed") +
  geom_hline(yintercept=-log10(0.01), col="black", linetype = "dashed")+
  scale_fill_manual(values=mycolors)+
  theme(axis.text.x = element_text(size=18,  colour = "black"),
        axis.text.y = element_text(size=18, colour = "black"))+
  theme(axis.title.x = element_text(size=19, face= "bold", colour = "black"))+
  theme(axis.title.y = element_text(size=19, face= "bold", colour = "black"))+
  theme(legend.text=element_text(size=14))+
  theme(legend.title=element_text(size=15))+
  ylim(0,25)+
  xlim(-5.5,5.5)+
  geom_text_repel(size = 5)+
  ggtitle("VitaminC - corticosteroids")

Steroids_VitCY

ggsave(plot = Steroids_VitCY, "Figures/Figure_3/limmaplot_VitC_D7vsD1_cortico.png", dpi = 1200,  height = 7,  width = 8)

################################################################################
#differential cytokine plot

#VitaminC without corticosteroids
#make matrix of values samples columns and frequencies as rows
count_comparison_tmydf = setNames(data.frame(t(analysis_VitC_N[,-1])), analysis_VitC_N[,1]) %>% as.matrix() 
str(count_comparison_tmydf)

dim(count_comparison_tmydf)

library(limma)

f <- factor(df_log_VitC_N$timepoint, levels=c("S1", "S3"))
design <- model.matrix(~0+f)
dim(design)
colnames(design) <- c("S1","S3")
#To make all pair-wise comparisons between the three groups one could proceed
fit <- lmFit(count_comparison_tmydf, design)
contrast.matrix <- makeContrasts(S3-S1,
                                 levels=design)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
topTable(fit2)
?eBayes

df_stats_VitC_N <- topTable(fit2, coef=1, adjust="BH", number = "inf") %>% rownames_to_column("Cytokine")

write_csv(df_stats_VitC_N, "CSV files/Output_files/Limma/Stats/VitC_D7vsD1_NoCortico.csv")
# add a column of NAs
df_stats_VitC_N$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
df_stats_VitC_N$diffexpressed[df_stats_VitC_N$logFC > 1 & df_stats_VitC_N$adj.P.Val < 0.01] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
df_stats_VitC_N$diffexpressed[df_stats_VitC_N$logFC < -1 & df_stats_VitC_N$adj.P.Val < 0.01] <- "DOWN"

# Create a new column "delabel" to de, that will contain the name of genes differentially expressed (NA in case they are not)
df_stats_VitC_N$delabel <- NA
df_stats_VitC_N$delabel[df_stats_VitC_N$diffexpressed != "NO"] <- df_stats_VitC_N$Cytokine[df_stats_VitC_N$diffexpressed != "NO"]

head(df_stats_VitC_N)

df_stats_VitC_N<- df_stats_VitC_N%>%
  mutate(delabel = recode(delabel, 
                          "IFNy" = "IFN𝛄", 
                          "TNFa" = "TNF⍺", 
                          "IL.6" = "IL-6", 
                          "IL.5" = "IL-5", 
                          "IL.1a" = "IL-1⍺", 
                          "IL.18" = "IL-18",
                          "IL.4" = "IL-4", 
                          "IFNa2" = "IFN⍺2", 
                          "IL.10" = "IL-10",
                          "IL.1RA" = "IL-1RA", 
                          "STNF.RI" = "sTNF-RI", 
                          "IL.12p40" = "IL-12p40",
                          "IFNb" = "IFNβ", 
                          "IL.1b" = "IL-1β", 
                          "IL.33" = "IL-33"))


mycolors <- c("#709eeb", "#b02428", "grey")
names(mycolors) <- c("DOWN", "UP", "NO")


Steroids_VitCN <- ggplot(data=df_stats_VitC_N, aes(x=logFC, y=-log10(adj.P.Val), label=delabel)) + 
  geom_point(aes(fill = diffexpressed), pch = 21, colour = "black", size = 5)+
  theme_classic()+
  geom_vline(xintercept=c(-0.5, 0.5), col="black", linetype = "dashed") +
  geom_hline(yintercept=-log10(0.01), col="black", linetype = "dashed")+
  scale_fill_manual(values=mycolors)+
  theme(axis.text.x = element_text(size=18,  colour = "black"),
        axis.text.y = element_text(size=18, colour = "black"))+
  theme(axis.title.x = element_text(size=19, face= "bold", colour = "black"))+
  theme(axis.title.y = element_text(size=19, face= "bold", colour = "black"))+
  theme(legend.text=element_text(size=14))+
  theme(legend.title=element_text(size=15))+
  ylim(0,25)+
  xlim(-5.5,5.5)+
  geom_text_repel(size = 5)+
  ggtitle("VitmanC - No corticosteroids")

Steroids_VitCN

 ggsave(plot = Steroids_VitCN, "Figures/Figure_3/limmaplot_VitC_D7vsD1_Nocortico.png", dpi = 1200,  height = 7,  width = 8)
