### make stacked boxplots of cytokines grouped by cluster

#load packages
library(ggplot2)
library(dplyr)
library(tidyverse)
###############################################################
###############################################################
# load dataframe - cytokine_df_LOVIT

# need to show proportions of each of each cluster
# Stacked + percent

#####################################
######### stacked proportions #######
#####################################

##IF NEEDED LOAD FILE 
cytokine_df_LOVIT <- read_csv("CSV files/patient_df_3clusters.csv")


colnames(cytokine_df_LOVIT)

cytokine_df_LOVIT <- cytokine_df_LOVIT[, -c(2, 3, 4, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 
                                            45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59,
                                            60, 61, 62, 63, 64)]

colnames(cytokine_df_LOVIT)

colnames(cytokine_df_LOVIT)[c(2, 7, 15:27)] <- c("IFN𝛄", "TNF⍺", "IL-6", "IL-5", "IL-1⍺", "IL-18",
                                                 "IL-4", "IFN⍺2", "IL-10", "IL-1RA", "sTNF-RI", "IL-12p40",
                                                 "IFNβ", "IL-1β", "IL-33")

colnames(cytokine_df_LOVIT)

######## ######## ######## ######## ######## ######## ######## ######## ######## 
## This bit reorganizes the dataframe
#if this doesn't work, break it down step by step with the code below 
#fold_change <- cytokine_df_LOVIT %>% dplyr::select(groups, 2:27) %>%
#  gather("Cytokine", "value", 2:27) %>%
#  group_by(cytokine_df_LOVIT[,28])%>%
#  na.omit #%>% #remove this next line for cell stuff because it logs the data and divides by the mean
#mutate_at(c("value"),funs(./median(.[Timepoint == "D1"]))) %>% as.data.frame()
 #also remove this line as we don’t want to log proportion data

test <- cytokine_df_LOVIT%>% group_by(cytokine_df_LOVIT[,28])

######## ######## ######## ######## ######## ######## 
######## the break down of the code above ########

# Select columns 2 to 24 and the 'groups' column
selected_data <- cytokine_df_LOVIT %>%
  dplyr::select(groups, 2:27)

# Check the structure of the selected data
str(selected_data)

# Reshape the data to long format
long_data <- selected_data %>%
  gather("Cytokine", "value", -groups)

str(long_data)

# Group by the 'groups' column and remove rows with NA values
fold_change<- long_data %>%
  group_by(groups)%>%
  na.omit()

########### end of breakdown ###################

fold_change$l2fc<-log(fold_change$value,10)

fold_change


## change the column name of 'groups' to 'subtype'
colnames(fold_change) [1] <- "Subtype" 
colnames(fold_change) [2] <- "Protein" 

#convert subtype into a factor 
fold_change$Subtype <- as.factor(fold_change$Subtype)
levels(fold_change$Subtype)  # Check the levels

library(ggplot2)
library(wesanderson)
library(RColorBrewer)

# plot Log2 fold change or (raw) value
foldchange <- ggplot(fold_change,aes(x=l2fc, y=Protein))+
  geom_boxplot(aes(fill = Subtype),outlier.shape = NA,alpha=1, position = position_dodge(width =1))+
  scale_fill_manual(values = c("1" = "#1BD76E", 
                               "2" = "#D76E1B", 
                               "3" = "purple4"))+
  #geom_point(aes(fill = Cluster))+
  #theme(axis.text.x = element_text(angle = 45, hjust = 0.2,vjust=0.4,size=4))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 0, hjust = 1,vjust=1,size=9))+
  #rotate_x_text(angle = 45, align = 1, valign = 0.25)+
  geom_hline(yintercept=0,linetype="dashed")+
  #geom_point(aes(fill= Subtype), shape= 21, position = position_jitterdodge(dodge.width = 1), size= 0.5)+
  #facet_grid(Cell_Type, scales = "", space = "free_y", ncol = 3)+
  #coord_flip()+
  theme(legend.position = "bottom")+
  theme(axis.text.x = element_text(size=15,  colour = "black"),
        axis.text.y = element_text(size=15, colour = "black"))+
  theme(axis.title.x = element_text(size=18, face= "bold", colour = "black"))+
  theme(axis.title.y = element_text(size=18, face= "bold", colour = "black"))+
  ylab("Protein")+
  xlab("Log10 Foldchange")+
  theme(legend.text=element_text(size=14))+
  theme(legend.title=element_text(size=15))+
  #coord_flip()
  xlim(0,5)

foldchange

png("Figures/Figure_1/stacked_boxplots.png", width=12,height=20,units="in",res=1200)
foldchange
dev.off()
########


