#### Heatmap cytokines with clustering



#¢##########################################################################################
#¢##########################################################################################

#define what makes the clusters different
install.packages("MEM")

# install MEM
devtools::install_github("cytolab/mem")


library(MEM)
library(heatmaply)
library(dplyr)

library(ggrepel)
library(tidyverse)


cytokine_df_LOVIT_log10 <- read.csv("CSV files/patient_df_log10_3clusters.csv")

########################################################################################################################
########################################################################################################################
########################################################################################################################
####.     LOVIT

# draw heatmap
colnames(cytokine_df_LOVIT_log10)

install.packages("ComplexHeatmap")
install.packages("ddply")
# draw heatmap
library(tidyverse)
library(RColorBrewer)
library(heatmaply)
library(colorRamp2)
library(ComplexHeatmap)
library(dplyr)
library(tibble)


cytokine_df_LOVIT_log10[2:27] <- scale(cytokine_df_LOVIT_log10[2:27])

names(cytokine_df_LOVIT_log10)[names(cytokine_df_LOVIT_log10) == "groups"] <- "Subtype"
names(cytokine_df_LOVIT_log10)[names(cytokine_df_LOVIT_log10) == "Primary"] <- "Primary trial outcome"
colnames(cytokine_df_LOVIT_log10)

colnames(cytokine_df_LOVIT_log10)[c(2, 7, 15:27)] <- c("IFN𝛄", "TNF⍺", "IL-6", "IL-5", "IL-1⍺", "IL-18",
                                                       "IL-4", "IFN⍺2", "IL-10", "IL-1RA", "sTNF-RI", "IL-12p40",
                                                       "IFNβ", "IL-1β", "IL-33")

#replace Treatment letters with actual treatment
cytokine_df_LOVIT_log10$Treatment <- gsub("A", "Vitamin C", cytokine_df_LOVIT_log10$Treatment)
cytokine_df_LOVIT_log10$Treatment <- gsub("B", "Placebo", cytokine_df_LOVIT_log10$Treatment)

cytokine_df_LOVIT_log10_11 <- cytokine_df_LOVIT_log10 %>% column_to_rownames(var = "sample_ID")


cytokine_df_LOVIT_log10_11 <- cytokine_df_LOVIT_log10_11 %>% select(1:26)%>% as.matrix()



#sort out colour scheme
col_fun = colorRamp2(c(-1.5, 0, 1.5), c("blue", "white", "red"))
col_fun(seq(-1.5, 1.5))

str(cytokine_df_LOVIT_log10)
class(cytokine_df_LOVIT_log10)


#annotation for row labels 
row_anno <- anno_mark(at = c(1:457), labels = rownames(cytokine_df_LOVIT_log10), which = "row",
                      link_width = unit(1, "mm"),
                      link_height = unit(1.5, "mm"),
                      padding = unit(0.55, "mm"),
                      labels_rot =  0,
                      extend = unit(0, "mm"),
                      labels_gp = gpar(fontsize = 2))
#draw label annotations for column 

colnames(cytokine_df_LOVIT_log10)
anno = anno_mark(at = c(2:27), labels = colnames(cytokine_df_LOVIT_log10), which = "column", side = "bottom",
                 link_width = unit(0.05, "mm"), 
                 link_height = unit(0.05, "mm"),
                 labels_rot = 90,
                 padding = unit(0, "mm"),
                 extend = unit(0, "mm"),
                 labels_gp = gpar(fontsize = 14))

#make annotation rows (patients)
colnames(cytokine_df_LOVIT_log10)

#anno_df <- cytokine_df_LOVIT_log10() %>%
#  data.frame() %>%
#  select("Subtype") %>%
#  mutate_if(is.factor,as.character)

anno_df <- cytokine_df_LOVIT_log10 %>% dplyr::select("Subtype", "Treatment",#"X28.day_mortality", 
                                        "Primary trial outcome"#,"covid"
) %>% mutate_if(is.character,as.factor)

######################################################################################################
library(ComplexHeatmap)
library(circlize)  # for color functions

# Assuming `m1mat` is the data matrix and `cytokine_df_clus2` includes a `Subtypes` column
# Convert `Subtypes` to a factor if it isn’t already, to ensure proper grouping
cytokine_df_LOVIT_log10$Subtype <- as.factor(cytokine_df_LOVIT_log10$Subtype)

# Define annotations in a DataFrame, `anno_df`, and ensure `Subtype` is a factor
anno_df$Subtype <- as.factor(anno_df$Subtype)  # Ensure Subtype is a factor for split
subtype_colors <- c("1" = "#1BD76E", "2" = "#D76E1B", "3" = "purple4")

# Create row annotations for `Subtype`, `Primary trial outcome`, and `Treatment`
ha <- HeatmapAnnotation(
  df = anno_df,
  col = list(
    Subtype = subtype_colors,
    "Primary trial outcome" = c("Yes" = "firebrick2", "No" = "dodgerblue"),
    Treatment = c("Vitamin C" = "orchid1", "Placebo" = "gold1")
  ),
  which = "row",               # Row annotations
  gp = gpar(col = "black", lwd = 0.01)
)


#table(cytokine_df_LOVIT_log10$)
 
colnames(anno_df)
head(anno_df)

##subtype colours
#"1" = "#1BD76E", "2" = "#D76E1B", "3" = "purple4"

ha = HeatmapAnnotation(df = anno_df, 
                       col = list("Subtype" = c("1" = "#1BD76E", "2" = "#D76E1B", "3" = "purple4"),
                                  #timepoint = c("S1" = "#2596be","S3" ="#be4d25"),
                                  #X28.day_mortality = c("Yes" = "firebrick2", "No" = "dodgerblue"),
                                  "Primary trial outcome" = c("Yes" = "firebrick2", "No" = "dodgerblue"),
                                  "Treatment" = c("Vitamin C" = "orchid1", "Placebo" = "gold1")
                                  #covid = c("Yes" = "darkslateblue","No" = "khaki4")
                       ),
                       which = "row",
                       gp = gpar(col = "black", lwd = 0.01))
ha
draw(ha, 1:457)

########################################################################

str(cytokine_df_LOVIT_log10_11)

m1 <- as.data.frame(cytokine_df_LOVIT_log10_11)
str(m1)
head(m1)
m1 [1:26]<- lapply(m1[1:26], function (x) as.numeric(as.character(x)))
#m1 <- m1[-c(24)]

m1mat <- as.matrix(m1)

# Perform hierarchical clustering on the rows of m1mat
H.fit <- hclust(dist(m1mat), method = "ward.D2")  # using Ward's method as an example

test <- heatmap(m1mat)

#################################


library(heatmaply)
######################################
#### draw heatmap and include clustering from above
hm <- Heatmap(m1mat, name = "Expression",
              col = col_fun,
              cluster_columns = TRUE,
              cluster_rows = TRUE, #enables clustering
              show_heatmap_legend = TRUE,
              #clustering_distance_rows = "euclidean",
              # clustering_method_rows = "ward",
              #bottom_annotation = columnAnnotation(mark=anno),
              # width = ncol(df4)*unit(6, "mm"), 
              # height = nrow(df4)*unit(4, "mm"),
              column_title_gp = gpar(fontsize = 3, fontface = "bold"),
              border = TRUE,
              #row_km = 5,
              #column_km = 3,
              show_column_names = TRUE,
              heatmap_legend_param = list(title = "Scale"),
              column_names_gp = grid::gpar(fontsize = 14),
              row_names_gp = grid::gpar(fontsize = 1),
              row_dend_width = unit(40, "mm"),
              column_dend_side = c("top", "bottom"),
              column_dend_height = unit(10, "mm"),
              right_annotation = ha,
              split = anno_df$Subtype)#+

#rowAnnotation(mark=row_anno)
hm



png("FIGURES/Figure_1/heatmap.png",width=12,height=12,units="in",res=1200)
hm
dev.off()
