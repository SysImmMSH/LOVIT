## LOVIT 
##figure 2
#Volcano plots - primary outcome

samples <- read.csv("CSV files/patient_df_3clusters.csv")

names(samples)
head(samples)
library(readxl)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(ggrepel)

pod <- read_excel("CSV files/LOVIT_MSH_EXPORT_DATA.xlsx") %>% dplyr::select("PatientID", "Primary")
names(pod)
head(pod)

samples <- merge(samples, pod, by.x = "patient_id", by.y = "PatientID")
head(samples)
names(samples)

df_log <- samples %>% mutate_at(vars(6:31), log2) %>% 
  select(6:31, sample_ID, groups, Primary.y) %>% as.data.frame() %>% 
  na.omit()
colnames(df_log)


df_log_1 <- df_log %>% dplyr::filter(groups == 1)
df_log_2 <- df_log %>% dplyr::filter(groups == 2)
df_log_3 <- df_log %>% dplyr::filter(groups == 3)

# split_df<- df_log %>%
#   group_by(groups) %>%
#   group_split()
# 
# df_log_1 <- split_df[[1]]

names(df_log_1)
analysis1 <- df_log_1 %>% select("sample_ID", 1:26) 
colnames(analysis1)

names(df_log_2)
analysis2 <- df_log_2 %>% select("sample_ID", 1:26) 
colnames(analysis2)

names(df_log_3)
analysis3 <- df_log_3 %>% select("sample_ID", 1:26) 
colnames(analysis3)

################################################################################
#differential cytokine plot

#make matrix of values samples columns and frequencies as rows
count_comparison_tmydf = setNames(data.frame(t(analysis1[,-1])), analysis1[,1]) %>% as.matrix() 
str(count_comparison_tmydf)

dim(count_comparison_tmydf)
dim(design)

library(limma)

f <- factor(df_log_1$Primary.y, levels=c("Yes", "No"))
design <- model.matrix(~0+f)
colnames(design) <- c("Yes","No")
#To make all pair-wise comparisons between the three groups one could proceed
fit <- lmFit(count_comparison_tmydf, design)
contrast.matrix <- makeContrasts(Yes-No,
                                 levels=design)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
topTable(fit2)
?eBayes
df_stats1 <- topTable(fit2, coef=1, adjust="BH", number = "inf") %>% rownames_to_column("Cytokine")
#write_csv(df_stats1, "Stats/CC_T0vSC_T0.csv")
# add a column of NAs
df_stats1$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
df_stats1$diffexpressed[df_stats1$logFC > 0.5 & df_stats1$adj.P.Val < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
df_stats1$diffexpressed[df_stats1$logFC < -0.5 & df_stats1$adj.P.Val < 0.05] <- "DOWN"

# Create a new column "delabel" to de, that will contain the name of genes differentially expressed (NA in case they are not)
df_stats1$delabel <- NA
df_stats1$delabel[df_stats1$diffexpressed != "NO"] <- df_stats1$Cytokine[df_stats1$diffexpressed != "NO"]

head(df_stats1)

df_stats1<- df_stats1%>%
  mutate(delabel = recode(delabel, 
                          "IL.18" = "IL-18", 
                          "STNF.RI" = "sTNF-RI", 
                          "IL.10" = "IL-10"))

mycolors <- c("#709eeb", "#b02428", "grey")
names(mycolors) <- c("DOWN", "UP", "NO")


Treat_YvsN <- ggplot(data=df_stats1, aes(x=logFC, y=-log10(adj.P.Val), label=delabel)) + 
  geom_point(aes(fill = diffexpressed), pch = 21, colour = "black", size = 5)+
  theme_classic()+
  geom_vline(xintercept=c(-0.5, 0.5), col="black", linetype = "dashed") +
  geom_hline(yintercept=-log10(0.05), col="black", linetype = "dashed")+
  scale_fill_manual(values=mycolors)+
  theme(axis.text.x = element_text(size=18,  colour = "black"),
        axis.text.y = element_text(size=18, colour = "black"))+
  theme(axis.title.x = element_text(size=19, face= "bold", colour = "black"))+
  theme(axis.title.y = element_text(size=19, face= "bold", colour = "black"))+
  theme(legend.text=element_text(size=14))+
  theme(legend.title=element_text(size=15))+
  ylim(0,6)+
  xlim(-2.2,2.2)+
  geom_text_repel(size = 8)+
  ggtitle("Subtype 1 - Y vs N")

Treat_YvsN

ggsave(plot = Treat_YvsN, "Figures/Figure_2/limmaplot_C1_primaryoutcome.png", dpi = 1200,  height = 7,  width = 8)


#########################
##### subtype 2

#make matrix of values samples columns and frequencies as rows
count_comparison_tmydf = setNames(data.frame(t(analysis2[,-1])), analysis2[,1]) %>% as.matrix() 
str(count_comparison_tmydf)

dim(count_comparison_tmydf)
dim(design)

library(limma)

f <- factor(df_log_2$Primary.y, levels=c("Yes", "No"))
design <- model.matrix(~0+f)
colnames(design) <- c("Yes","No")
#To make all pair-wise comparisons between the three groups one could proceed
fit <- lmFit(count_comparison_tmydf, design)
contrast.matrix <- makeContrasts(Yes-No,
                                 levels=design)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
topTable(fit2)
?eBayes
df_stats2 <- topTable(fit2, coef=1, adjust="BH", number = "inf") %>% rownames_to_column("Cytokine")
#write_csv(df_stats1, "Stats/CC_T0vSC_T0.csv")
# add a column of NAs
df_stats2$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
df_stats2$diffexpressed[df_stats2$logFC > 0.5 & df_stats2$adj.P.Val < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
df_stats2$diffexpressed[df_stats2$logFC < -0.5 & df_stats2$adj.P.Val < 0.05] <- "DOWN"

# Create a new column "delabel" to de, that will contain the name of genes differentially expressed (NA in case they are not)
df_stats2$delabel <- NA
df_stats2$delabel[df_stats2$diffexpressed != "NO"] <- df_stats2$Cytokine[df_stats2$diffexpressed != "NO"]

Treat_YvsN <- ggplot(data=df_stats2, aes(x=logFC, y=-log10(adj.P.Val), label=delabel)) + 
  geom_point(aes(fill = diffexpressed), pch = 21, colour = "black", size = 5)+
  theme_classic()+
  geom_vline(xintercept=c(-0.5, 0.5), col="black", linetype = "dashed") +
  geom_hline(yintercept=-log10(0.05), col="black", linetype = "dashed")+
  scale_fill_manual(values=mycolors)+
  theme(axis.text.x = element_text(size=18,  colour = "black"),
        axis.text.y = element_text(size=18, colour = "black"))+
  theme(axis.title.x = element_text(size=19, face= "bold", colour = "black"))+
  theme(axis.title.y = element_text(size=19, face= "bold", colour = "black"))+
  theme(legend.text=element_text(size=14))+
  theme(legend.title=element_text(size=15))+
  ylim(0,6)+
  xlim(-2.2,2.2)+
  geom_text_repel(size = 8)+
  ggtitle("Subtype 2 - Y vs N")

Treat_YvsN

ggsave(plot = Treat_YvsN, "Figures/Figure_2/limmaplot_C2_primaryoutcome.png", dpi = 1200,  height = 7,  width = 8)


#########################
##### subtype 3

#make matrix of values samples columns and frequencies as rows
count_comparison_tmydf = setNames(data.frame(t(analysis3[,-1])), analysis3[,1]) %>% as.matrix() 
str(count_comparison_tmydf)

dim(count_comparison_tmydf)
dim(design)

library(limma)

f <- factor(df_log_3$Primary.y, levels=c("Yes", "No"))
design <- model.matrix(~0+f)
colnames(design) <- c("Yes","No")
#To make all pair-wise comparisons between the three groups one could proceed
fit <- lmFit(count_comparison_tmydf, design)
contrast.matrix <- makeContrasts(Yes-No,
                                 levels=design)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
topTable(fit2)
?eBayes
df_stats3 <- topTable(fit2, coef=1, adjust="BH", number = "inf") %>% rownames_to_column("Cytokine")
#write_csv(df_stats1, "Stats/CC_T0vSC_T0.csv")
# add a column of NAs
df_stats3$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
df_stats3$diffexpressed[df_stats3$logFC > 0.5 & df_stats3$adj.P.Val < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
df_stats3$diffexpressed[df_stats3$logFC < -0.5 & df_stats3$adj.P.Val < 0.05] <- "DOWN"

# Create a new column "delabel" to de, that will contain the name of genes differentially expressed (NA in case they are not)
df_stats3$delabel <- NA
df_stats3$delabel[df_stats3$diffexpressed != "NO"] <- df_stats3$Cytokine[df_stats3$diffexpressed != "NO"]


df_stats3<- df_stats3%>%
  mutate(delabel = recode(delabel, "STNF.RI" = "sTNF-RI"))

Treat_YvsN <- ggplot(data=df_stats3, aes(x=logFC, y=-log10(adj.P.Val), label=delabel)) + 
  geom_point(aes(fill = diffexpressed), pch = 21, colour = "black", size = 5)+
  theme_classic()+
  geom_vline(xintercept=c(-0.5, 0.5), col="black", linetype = "dashed") +
  geom_hline(yintercept=-log10(0.05), col="black", linetype = "dashed")+
  scale_fill_manual(values=mycolors)+
  theme(axis.text.x = element_text(size=18,  colour = "black"),
        axis.text.y = element_text(size=18, colour = "black"))+
  theme(axis.title.x = element_text(size=19, face= "bold", colour = "black"))+
  theme(axis.title.y = element_text(size=19, face= "bold", colour = "black"))+
  theme(legend.text=element_text(size=14))+
  theme(legend.title=element_text(size=15))+
  ylim(0,6)+
  xlim(-2.2,2.2)+
  geom_text_repel(size = 8)+
  ggtitle("Subtype 3 - Y vs N")

Treat_YvsN

ggsave(plot = Treat_YvsN, "Figures/Figure_2/limmaplot_C3_primaryoutcome.png", dpi = 1200,  height = 7,  width = 8)

