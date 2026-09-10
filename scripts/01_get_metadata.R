# Load libraries
library(recount3)
library(dplyr)
library(stringr)

# Get full sample metadata for GSE104704 (SRP119561) via recount3
human_projects <- available_projects()
proj_info <- subset(human_projects, project == "SRP119561")
rse <- create_rse(proj_info)

sample_info <- as.data.frame(colData(rse)[, c("external_id", "sra.sample_attributes")])

# Extract condition (Young/Old/AD) from sample attributes
sample_info$condition <- str_extract(sample_info$sra.sample_attributes, 
                                     "(?<=study group;;)[^|]+")
sample_info$condition <- ifelse(sample_info$condition == "Aged, diseased", "AD", sample_info$condition)

table(sample_info$condition)

# Save clean metadata for use in downstream scripts
write.csv(sample_info[, c("external_id", "condition")], 
          "data/sample_metadata_full.csv", row.names = FALSE)
