#! /usr/bin/env Rscript

# Title: Calculating bacterial richness from Bracken Report
# Name: Jamal Sheriff
# Date: 13 May 2025 

# DESCRIPTION: R Script will return the number of species within the Bacteria Domain.
# All species with less than 10 fragment alignments to reference genome has been 
# removed by default. 
#
# This function requires two arguments. args[1] = input file, args[2] = output file name with csv extension 
 
# Load packages 
library.dynam("dplyr", package = c("dplyr"), lib.loc = "/home/ac.jsheriff/dev_container/modules/app_service/lib/")
library.dynam("tibble", package = c("tibble"), lib.loc = "/home/ac.jsheriff/dev_container/modules/app_service/lib/")


library("dplyr")
library("tibble")

# (.packages())

# set command line arguments ----
args <- commandArgs(trailingOnly = TRUE)
filename <- args[1]
outputfile <- args[2]

# print(filename)
#stop the script if no command line argument
if(length(args)==0){
  print("Please select bracken_report.txt input file.")
  stop("Requires command line argument.")
}
# Calculating species richness

calc_richness_bracken <- function(file) {
  df = read.delim(file, header = FALSE)
  
  df %>% # keep only rows above where the second domain is introduced -- this part needs revisiting
    dplyr::filter(! dplyr::cumany(grepl("Eukaryota", V6, fixed = TRUE))) ->
    df_filtered
  
  data.frame() %>% # create DF with a row for each sample 
    matrix(ncol = 0, nrow = 1) -> richness_df
  
  row.names(richness_df) = strsplit(file, "\\_bracken")[[1]][1] # set rownames after beginning of filename 
  
  for(rows in df_filtered$V4) { # for rows in V4 column, count rows containing  "S" values
    if(nrow(df_filtered[df_filtered$V4 == "S",]) > 0){ #If S is present, count S for species richness (For WGS Samples)
    Species_Richness = nrow(df_filtered[df_filtered$V4 == "S",])
    }
    else { # If S is not present, count G instead (for 16S Samples)
      Species_Richness = nrow(df_filtered[df_filtered$V4 == "G",])
    }
    }
  richness_df <- as.data.frame(cbind(richness_df, Species_Richness)) # add richness values to output DF
  richness_df <- tibble::rownames_to_column(richness_df, "Sample")
  return(richness_df)
} 

richness_bracken <- 
calc_richness_bracken(filename)

richness_bracken$Sample <- as.character(richness_bracken$Sample)
richness_bracken$Species_Richness <- as.numeric(richness_bracken$Species_Richness)

# str(richness_bracken)

write.csv(richness_bracken, file = outputfile, row.names = FALSE) 
