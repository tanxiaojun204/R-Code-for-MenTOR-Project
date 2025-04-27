########################################
# SECTION 1: SETUP WORKING DIRECTORY 
########################################

setwd("/Users/xiaojun/OneDrive - Imperial College London/Desktop")
library(tidyverse)
library(dplyr)
library(readxl)


########################################
# SECTION 2: IMPORT QUANTIFICATION LIMIT DATA
########################################

### Quantification limit data file
LLOQ = read.csv("LLOQ.csv") |> select(!X) 


########################################
# SECTION 3: IMPORT AND INITIAL CLEANING OF MenTOR 36-PLEX DATA
########################################

### MenTOR data file
path <- "/Users/xiaojun/OneDrive - Imperial College London/Desktop/MenTOR 36-plex results 16062023_vs sent 30 June.xlsx"

### Import and cleaning
RawData = read_excel("MenTOR 36-plex results 16062023_vs sent 30 June.xlsx", sheet = "Data")|> 
  select(PlateName = `Plate Name`, Assay, Sample, CalcConcMean = `Calc. Conc. Mean`, CV) |> # please include Well to check layout
  ### Column correction, uniform variable names
  ### Replace space with _ in Sample and Plate_Name columns
  mutate(
    Sample = gsub("\\s+", "_", Sample),
    Plate_Name = gsub("\\s+", "_", PlateName)
  ) |>
  ### Add panel and plate number according to the excel file
  mutate(
    Panel = case_when(
      Plate_Name == "Plate_2BM3RA6337" ~ "Chemokine",
      Plate_Name == "Plate_2BM3RAS369" ~ "Chemokine",
      Plate_Name == "Plate_2BM0JAR039" ~ "TH17",
      Plate_Name == "Plate_2BM0JAV589" ~ "TH17",
      Plate_Name == "Plate_2BLHIAUA71" ~ "Cytokine",
      Plate_Name == "Plate_2BLHIAN565" ~ "Cytokine",
      Plate_Name == "Plate_2BM1MAM546" ~ "Pro-Inflammatory",
      Plate_Name == "Plate_2BM1MAS549" ~ "Pro-Inflammatory"
    ),
    # Notice the chemokine plates have wrong label (i.e. Chemokine plate 1 should be plate 2)
    # the chemokine panels layouts were flipped horizontally (please refer to excel file)
    Plate = case_when(
      Plate_Name == "Plate_2BM3RA6337" ~ 2,
      Plate_Name == "Plate_2BM3RAS369" ~ 1,
      Plate_Name == "Plate_2BM0JAR039" ~ 1,
      Plate_Name == "Plate_2BM0JAV589" ~ 2,
      Plate_Name == "Plate_2BLHIAUA71" ~ 1,
      Plate_Name == "Plate_2BLHIAN565" ~ 2,
      Plate_Name == "Plate_2BM1MAM546" ~ 1,
      Plate_Name == "Plate_2BM1MAS549" ~ 2
    ),
    ### Correct analyte names
    Assay = case_when(
      Assay == "IL-12" ~ "IL-12/IL-23p40",
      Assay == "IL-17" ~ "IL-17A",
      Assay == "IL17A" ~ "IL-17A",
      Assay == "VEGF" ~ "VEGF-A",
      str_detect(Assay, "^IL\\d{2}$") ~ str_replace(Assay, "^IL(\\d{2})$", "IL-\\1"),
      TRUE ~ Assay
    )
  ) |> 
  ### Replace string NaN to NA, and change Calc_Conc_Mean to numerical
  mutate(
    Calc_Conc_Mean = ifelse(CalcConcMean == "NaN", NA, CalcConcMean),
    CV = ifelse(CV == "NaN", NA, CV),
    Calc_Conc_Mean = as.numeric(Calc_Conc_Mean),
    CV = as.numeric(CV)
  )

### Add quantification limits and intraasay CV threshold
RawData <- RawData |> 
  left_join(LLOQ, by = c("Panel", "Assay")) |> # NEED TO ADD LLOQ AND ULOQ OF IL-17A IN TH17 PANEL !!!
  mutate(
    Q_Range = case_when(
      is.na(Calc_Conc_Mean) ~ "<LLOD",
      Calc_Conc_Mean < LLOQ ~ "<LLOQ",
      Calc_Conc_Mean > ULOQ ~ ">ULOQ",
      Calc_Conc_Mean >= LLOQ & Calc_Conc_Mean <= ULOQ ~ "Within"
    ),
    CV_Lim = ifelse((Assay == "Eotaxin-3" | Assay == "MCP-4"), 26, 21),
  ) |>
  ### Make them single replicate 
  filter(row_number() %% 2 == 1) |> relocate(Panel, Plate) |> select(!Plate_Name) |>
  ### Remove IL-17A on TH17 panel
  filter(!(Panel == "TH17" & Assay == "IL-17A"))

# Check assay names
unique(RawData$Assay) |> view() #35


########################################
# SECTION 4: GENERATE IMPUTATION DATA
########################################

### Generate Imputation data
# lowest detectable reading per plate per analyte
Imp_Data = RawData |> group_by(Panel, Plate, Assay) |>
  filter(!is.na(Calc_Conc_Mean)) |>
  summarise(LLOD_Imp = min(Calc_Conc_Mean)) |> 
  left_join(LLOQ) |>
  mutate(
    # LLOQ imp as mean of LLOD estimate and LLOQ
    LLOQ_Imp = (LLOD_Imp + LLOQ)/2,
    LLOD_Imp = if_else(LLOD_Imp > LLOQ_Imp, LLOQ_Imp, LLOD_Imp)
  ) |>
  select(Panel, Plate, Assay, LLOD_Imp, LLOQ_Imp)


########################################
# SECTION 5: MenTOR DATA (FILTERING, HANDLING MISSING CV, IMPUTATION)
########################################

MenTOR_Data = RawData |> filter(grepl("MEN", Sample)) |> 
  filter(!(Q_Range == "Within" & CV > CV_Lim)) # remove unreliable measurements
MenTOR_Data |> filter(is.na(CV)) |> nrow()
MenTOR_Data |> group_by(Sample, Assay, Panel) |> filter(n()>0) |> view() #MEN2013 Chemokine Plate 2 


### Handle NA values in CV column
MenTOR_Data = MenTOR_Data |> 
  mutate(CV = replace_na(CV, Inf)) |> # Replace NA with Inf for consistent filtering
  group_by(Sample, Assay, Panel) |> 
  filter(
    if (sum(is.na(Calc_Conc_Mean)) == 1) { 
      !is.na(Calc_Conc_Mean) # Keep non-NA value
    } else if (sum(Q_Range == "Within") == 1) {
      Q_Range == "Within" # Prefer values within quantification range
    } else {
      CV == min(CV, na.rm = TRUE) # Select based on minimum CV
    }
  ) |> 
  ungroup() |> 
  mutate(CV = ifelse(CV == Inf, NA, CV)) # Restore NA for CV after filtering

### Imputation 
MenTOR_Data = MenTOR_Data  |> 
  left_join(Imp_Data) |>
  mutate(
    Calc_Conc_Mean = case_when(
      Q_Range == "<LLOD" ~ LLOD_Imp, # Use the lowest detectable reading
      Q_Range == "<LLOQ" ~ 0.5 * LLOQ, # Fixed-point imputation for <LLOQ samples
      Q_Range == ">ULOQ" ~ ULOQ, # Use ULOQ for >ULOQ samples
      TRUE ~ Calc_Conc_Mean # Keep original value for "Within" range
    )
  )

MenTOR_Data |> filter(Q_Range == "<LLOQ") |> summarise(Count = n())
MenTOR_Data |> filter(Q_Range == "<LLOQ") |> view()
MenTOR_Data |> filter(Q_Range == "<LLOQ") |> group_by(Panel, Assay) |> summarise(Count = n())


########################################
# SECTION 6: COMPARATOR DATA
########################################

Comparator_Data = RawData |> filter(grepl("SF", Sample) | grepl("RA", Sample) | grepl("OMB", Sample)) |>
  mutate(
    Matrix = case_when(
      str_detect(Sample, "SF") ~ "SF",  
      TRUE ~ str_extract(Sample, "[^_]+$")
    ),
    Sample = case_when(
      str_detect(Sample, "RA") ~ "RA_POOL",
      str_detect(Sample, "OMB") ~ "OMB106",
      Sample == "SF10_L" ~ "10L",
      str_detect(Sample, "SF") ~ str_replace(Sample, "SF", ""),
    )
  )

Comparator_Data |> group_by(Sample, Assay, Panel) |>
  summarise(N = n()) |> filter(N>1) |> view()


########################################
# SECTION 7: SUMMARIZE DATA BY PANEL AND ASSAY
########################################

# Summarize data by Panel and Assay
Protein_Summary <- MenTOR_Data |> 
  group_by(Panel, Assay) |> 
  summarise(
    Total_Samples = n(),
    Passed_QC = sum(CV <= CV_Lim, na.rm = TRUE), # Count samples passing CV threshold
    `% Passed Sample QC` = round((Passed_QC / Total_Samples) * 100, 1),
    `% Within` = round((sum(Q_Range == "Within", na.rm = TRUE) / Total_Samples) * 100, 1),
    `% <LLOQ` = round((sum(Q_Range == "<LLOQ", na.rm = TRUE) / Total_Samples) * 100, 1),
    `% >ULOQ` = round((sum(Q_Range == ">ULOQ", na.rm = TRUE) / Total_Samples) * 100, 1)
  ) |> 
  arrange(Panel, Assay) # Arrange for clarity

# Reorder "Within" to descend
Protein_Summary <- Protein_Summary |> 
  group_by(Panel) |>  # Group by Panel
  arrange(desc(`% Within`), .by_group = TRUE) # Arrange within each Panel by descending % Within

Protein_Summary |> view()

# Update Protein_Summary to exclude <LLOD and recalculate metrics
Protein_Summary <- MenTOR_Data |> 
  filter(Q_Range != "<LLOD") |>  # Remove <LLOD values
  group_by(Panel, Assay) |> 
  summarise(
    Total_Samples = n(),
    Passed_QC = sum(CV <= CV_Lim, na.rm = TRUE), # Count samples passing CV threshold
    `% Passed Sample QC` = round((Passed_QC / Total_Samples) * 100, 1),
    `% Within` = round((sum(Q_Range == "Within", na.rm = TRUE) / Total_Samples) * 100, 1),
    `% <LLOQ` = round((sum(Q_Range == "<LLOQ", na.rm = TRUE) / Total_Samples) * 100, 1),
    `% >ULOQ` = round((sum(Q_Range == ">ULOQ", na.rm = TRUE) / Total_Samples) * 100, 1),
    .groups = "drop"  # Drop the grouping after summarisation
  ) |> 
  group_by(Panel) |> 
  arrange(Panel, desc(`% Within`), .by_group = TRUE) # Group by Panel and sort by % Within descending

# Step 1: Calculate Intra-assay CV
Intra_Assay_CV <- RawData |> 
  filter(grepl("MEN", Sample)) |>  
  filter(Q_Range == "Within") |>   # only include values within range
  group_by(Panel, Plate, Assay) |> 
  summarise(
    Intra_CV_Plate = ifelse(n() > 1, 
                            sd(Calc_Conc_Mean, na.rm = TRUE) / mean(Calc_Conc_Mean, na.rm = TRUE) * 100, 
                            NA), 
    .groups = "drop"
  ) |>
  group_by(Panel, Assay) |>
  summarise(Intra_CV = mean(Intra_CV_Plate, na.rm = TRUE), .groups = "drop")


# Debugging: Check Intra-assay CV output
print("Intra_Assay_CV:")
Intra_Assay_CV |> print()

# Step 2: Calculate Inter-assay CV
Inter_Assay_CV <- RawData |> 
  filter(grepl("MEN", Sample)) |> 
  filter(Q_Range == "Within") |> 
  group_by(Panel, Assay) |> 
  summarise(
    Inter_CV = ifelse(n() > 1, 
                      sd(Calc_Conc_Mean, na.rm = TRUE) / mean(Calc_Conc_Mean, na.rm = TRUE) * 100, 
                      NA), 
    .groups = "drop"
  )


# Check Inter-assay CV output
print("Inter_Assay_CV:")
Inter_Assay_CV |> print()

# Step 3: Merge CV Metrics into Protein Summary
Protein_Summary <- Protein_Summary |> 
  left_join(Intra_Assay_CV, by = c("Panel", "Assay")) |> 
  left_join(Inter_Assay_CV, by = c("Panel", "Assay"))

# Check final Protein_Summary
print("Updated Protein_Summary:")
Protein_Summary |> print()

# Filter out IL-8 from Chemokine Panel
Protein_Summary <- Protein_Summary |> 
  filter(!(str_trim(Assay) == "IL-8" & str_trim(Panel) == "Chemokine"))

# Filter out the TH17 panel
Protein_Summary <- Protein_Summary |> 
  filter(Panel != "TH17")

# View the updated Protein Summary 
Protein_Summary |> view()

#Calculate unique sample
unique_samples <- RawData |> 
  filter(grepl("MEN", Sample)) |> # Ensure you're filtering correctly for participant samples
  summarise(Unique_Participants = n_distinct(Sample)) |> # Count unique values
  pull(Unique_Participants) # Extract the number from the summarised result
print(unique_samples)

# Summarize minimum, median, and maximum Total_Samples from Protein Summary table
protein_summary_summary <- Protein_Summary |> 
  summarise(
    Min_Samples = min(Total_Samples, na.rm = TRUE),
    Median_Samples = median(Total_Samples, na.rm = TRUE),
    Max_Samples = max(Total_Samples, na.rm = TRUE)
  )

# Display the summarized table
print(protein_summary_summary)

# Summarize median, max, and min of % Within by Panel
Within_Summary <- Protein_Summary |> 
  group_by(Panel) |> 
  summarise(
    Median_Within = median(`% Within`, na.rm = TRUE),
    Max_Within = max(`% Within`, na.rm = TRUE),
    Min_Within = min(`% Within`, na.rm = TRUE)
  )

# View the summary table
Within_Summary |> view()



########################################
# SECTION 8: FINAL SAMPLE SUMMARY
########################################

# Calculate total samples before filtering
Total_Samples_Before <- RawData |> 
  filter(grepl("MEN", Sample)) |> 
  summarise(Total_Before_Filtering = n()) |> 
  pull(Total_Before_Filtering)

# Calculate total samples after filtering using Protein Summary
Total_Samples_After <- Protein_Summary |> 
  summarise(Total_After_Filtering = sum(Total_Samples, na.rm = TRUE)) |> 
  pull(Total_After_Filtering)

# Calculate quantifiable samples before and after filtering
Quantifiable_Samples <- Protein_Summary |> 
  summarise(
    Quantifiable_Before_Filter = sum(Total_Samples, na.rm = TRUE),
    Quantifiable_After_Filter = sum(Passed_QC, na.rm = TRUE) # Using % Passed Sample QC
  )

# Combine these metrics into one table
Sample_Summary <- tibble(
  Total_Before_Filtering = Total_Samples_Before,
  Total_After_Filtering = Total_Samples_After,
  Quantifiable_Before_Filter = Quantifiable_Samples$Quantifiable_Before_Filter,
  Quantifiable_After_Filter = Quantifiable_Samples$Quantifiable_After_Filter
)

# View the summary table
Sample_Summary |> view()


########################################
# SECTION 9: GENERATE HTML TABLE
########################################

# 1) Install (once) if you haven't yet:
# install.packages("kableExtra")

# 2) Load the libraries
library(knitr)
library(kableExtra)

# 3) Make a nice HTML table from Protein_Summary:
Protein_Summary %>%
  # If you want to reorder columns, do it here
  select(Panel, Assay, Total_Samples, `% Passed Sample QC`,
         `% Within`, `% <LLOQ`, `% >ULOQ`, Intra_CV, Inter_CV) %>%
  # Create a basic kable output
  kable(
    format = "html",                # Produce an HTML table
    caption = "Protein Summary Table", 
    digits = 2                      # Number of decimals for numeric columns
  ) %>%
  # Add some styling
  kable_styling(
    bootstrap_options = c("striped", "hover", "condensed", "responsive"),
    full_width = FALSE
  )


########################################
# SECTION 10: POST QC TABLE AND WIDE FORMAT
########################################

Post_QC_Table <- Protein_Summary |>
  filter(`% Within` > 75) |>
  select(Panel, Assay, `% Within`) |>
  left_join(MenTOR_Data, by = c("Panel", "Assay")) |>
  filter(CV <= CV_Lim) |>
  select(Sample, Panel, Assay, Calc_Conc_Mean, `% Within`)

# View the resulting Post QC table
Post_QC_Table |> view()

# Reshape Post QC Table and remove duplicates
Post_QC_Wide <- Post_QC_Table %>%
  select(Sample, Assay, Calc_Conc_Mean) %>% # Keep only relevant columns
  distinct(Sample, Assay, .keep_all = TRUE) %>% # Remove duplicate Sample-Assay pairs
  pivot_wider(
    names_from = Assay,          # Use Assay as column headers
    values_from = Calc_Conc_Mean # Use Calc_Conc_Mean for values
  ) %>%
  arrange(Sample)                # Order rows by Sample

# View the resulting table
Post_QC_Wide |> view()

