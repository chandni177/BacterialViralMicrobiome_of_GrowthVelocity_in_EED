library("phyloseq")
library("tidyverse")

physeqPriorAllEEPhage <- readRDS("../data/physeqObjects/physeqPriorAllEEPhage.RDS") # 1004 x 83
phageSampleData <- as(sample_data(physeqPriorAllEEPhage), "data.frame")
phageSampleData <- phageSampleData %>% 
  mutate("unique_id" = str_extract(sample_name, 
                                   pattern = "(-|_)[:digit:]+(-|_)LNS_?[:digit:]+(-|_)"),
         "unique_id" = str_remove_all(unique_id, pattern = "(_|-)"))
phageSampleDataFinal <- select(phageSampleData,
                               unique_id, sample_name, Growth_Status_Prior,
                               age_months, sex, WHZ, HAZ, breastfed_months,
                               intervention)

bacterialSampleData <- read.delim("../documents/analyzedBacterialMetadata.txt")
bacterialSampleDataFinal <- bacterialSampleData %>% 
  mutate("unique_id" = str_extract(sample_name, 
                                   pattern = "(-|_)[:digit:]+(-|_)LNS_?[:digit:]+(-|_)"),
         "unique_id" = str_remove_all(unique_id, pattern = "(_|-)"))

combinedData <- merge(phageSampleDataFinal,
                      bacterialSampleDataFinal,
                      by = "unique_id",
                      all = TRUE)

write.table(combinedData, "../documents/combinedMetadata.txt",
           quote = FALSE, row.names = FALSE, sep = "\t")

#---- Generate Statistics -----#

# n
#   Growth_Status_Prior

# Median, Q1, Q3: 
#   age_months
#   WAZ
#   HAZ
#   breastfed_months

# Percentage
#   gender (Male)
#   intervention

analyzedData <- read.delim("../documents/all_analyzed_sample_data.txt")

levels(analyzedData$Growth_Status_Prior)
analyzedData$Growth_Status_Prior <- factor(analyzedData$Growth_Status_Prior,
                                           levels = c("Poor", "Moderate", "Adequate"))

table(analyzedData$Growth_Status_Prior)

table <- analyzedData %>%
  dplyr::group_by(Growth_Status_Prior) %>% 
  dplyr::summarize("median_age" = median(age_months),
                   "Q1_age" = quantile(age_months, 0.25),
                   "Q3_age" = quantile(age_months, 0.75),
                   "median_WAZ" = median(WHZ),
                   "Q1_WHZ" = quantile(WHZ, 0.25),
                   "Q3_WHZ" = quantile(WHZ, 0.75),
                   "median_HAZ" = median(HAZ),
                   "Q1_HAZ" = quantile(HAZ, 0.25),
                   "Q3_HAZ" = quantile(HAZ, 0.75),
                   "median_breastfed" = median(breastfed_months, na.rm = TRUE),
                   "Q1_breastfed" = quantile(breastfed_months, 0.25, na.rm = TRUE),
                   "Q3_breastfed" = quantile(breastfed_months, 0.75, na.rm = TRUE))

# gender
analyzedData %>% 
  dplyr::group_by(Growth_Status_Prior, gender) %>% 
  dplyr::summarize(n())

# intervention
analyzedData %>% 
  dplyr::group_by(Growth_Status_Prior, intervention) %>% 
  dplyr::summarise(n())

# delta HAZ
analyzedSamplesWithDeltas <- read.delim("../documents/analyzedSamplesWithDeltas.txt",
                                        sep = "\t", header = TRUE)
analyzedSamplesWithDeltas %>% 
  dplyr::group_by(Growth_Status_Prior) %>% 
  dplyr::summarize("median_delta" = median(delta_HAZ),
                   "Q1_delta" = quantile(delta_HAZ, 0.25),
                   "Q2_delta" = quantile(delta_HAZ, 0.75))
