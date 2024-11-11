# LOVIT
#Table of concentrations
install.packages("dplyr")
library(tidyverse)
library(dplyr)
library(readxl)
library(writexl)
library(officer)

all_cytokines <- read_csv("CSV files/full_df_MM/Complete_df.csv")
LOD <- read_excel ("CSV files/cyotkine LOD.xlsx")

colnames(all_cytokines)
head(LOD)
str(all_cytokines)

#select the columns of interest
selected_cyto <- all_cytokines %>% dplyr::select(c(5:17, 19:32, 66))

colnames(selected_cyto)

# Convert the wide format to long format
cytokine_data_long <- selected_cyto %>%
  pivot_longer(
    cols = -Treatment,             # All columns except Treatment
    names_to = "Cytokine",          # Name for the new cytokine column
    values_to = "Value"             # Name for the values column
  )

# View the long format data
print(cytokine_data_long)


# Function to calculate overall Median (IQR)
calculate_median_iqr <- function(column) {
  median_value <- median(column, na.rm = TRUE)
  q1 <- quantile(column, 0.25, na.rm = TRUE)
  q3 <- quantile(column, 0.75, na.rm = TRUE)
  paste0(round(median_value, 2), " (", round(q1, 2), " - ", round(q3, 2), ")")
}


# Calculate overall Median (IQR) for each cytokine
overall_summary <- cytokine_data_long %>%
  group_by(Cytokine) %>%
  summarize(Median_IQR = calculate_median_iqr(Value), .groups = "drop")

# Calculate Median (IQR) by Treatment for each cytokine
treatment_summary <- cytokine_data_long %>%
  group_by(Cytokine, Treatment) %>%
  summarize(Treatment_Median_IQR = calculate_median_iqr(Value), .groups = "drop") %>%
  pivot_wider(names_from = Treatment, values_from = Treatment_Median_IQR, names_prefix = "Treatment_")

# Check the results
print(overall_summary)
print(treatment_summary)

# make sure all the cytokines are labelled consistently in each df
overall_summary<- overall_summary%>%
  mutate(Cytokine = recode(Cytokine, 
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

treatment_summary<- treatment_summary%>%
  mutate(Cytokine = recode(Cytokine, 
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


LOD<- LOD%>%
  mutate(Cytokine = recode(Cytokine, 
                           "IFNy" = "IFN𝛄", 
                           "TNFa" = "TNF⍺", 
                           "IL.6" = "IL-6", 
                           "IL.5" = "IL-5", 
                           "IL-1a" = "IL-1⍺", 
                           "IL.18" = "IL-18",
                           "IL.4" = "IL-4", 
                           "IFNa2" = "IFN⍺2", 
                           "IL.10" = "IL-10",
                           "IL.1RA" = "IL-1RA", 
                           "STNF-RI" = "sTNF-RI", 
                           "IL.12p40" = "IL-12p40",
                           "IFNb" = "IFNβ", 
                           "IL-1b" = "IL-1β", 
                           "IL.33" = "IL-33"))

# Join all tables to create the final summary table
final_table <- overall_summary %>%
  # Step 1: Join with LLOD and ULOD data to create LLOD_ULOD column
  left_join(LOD, by = "Cytokine") %>%
  mutate(LLOD_ULOD = paste0(round(LLOD, 4), " - ", ULOD)) %>%
  # Step 2: Join with treatment summary data to add Treatment_A and Treatment_B
  left_join(treatment_summary, by = "Cytokine") %>%
  # Select and reorder columns for final table
  select(Cytokine, LLOD_ULOD, Median_IQR, Treatment_B, Treatment_A)

# Print the final table
print(final_table)


# Export to CSV
write.csv(final_table, "CSV files/Output_files/cytokine_concentration_table.csv", row.names = FALSE)

# Export to Excel
write_xlsx(final_table, "CSV files/Output_files/cytokine_concentration_table.xlsx")

# Create Word document
doc <- read_docx()
doc <- doc %>%
  body_add_par("Cytokine Concentrations", style = "heading 1") %>%
  body_add_table(final_table, style = "table_template")
print(doc, target = "cytokine_concentration_table.docx")


####################################