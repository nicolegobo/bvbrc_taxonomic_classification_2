#!/usr/bin/env Rscript

library(argparse)
library(DT)
library(pavian)
library(htmltools)

parser <- ArgumentParser(description= 'This script creates an multi sample comparison table HTML from a Kracken2 report')

parser$add_argument('--input', '-i', help= 'Path one or multiple Kraken2 Report')
parser$add_argument('--output_dir', '-o', help= 'Path where you would like all output files (multisample_comarison.HTML, summary.HTML, summary.CSV)')
xargs<- parser$parse_args()

### Get sample names from input files to name columns of report ###
input_list <- unlist(strsplit(xargs$input, " "))

# Initialize an empty list
sample_name_list <- list()
input_file_list <- list()
for (file in input_list) {
  if (substring(file, nchar(file)-2, nchar(file)) == "txt") {
   input_file_list <- c(input_file_list, file)
  # get file name
  pattern <- "^(.*/)?([^/]*)_k2_report\\.txt$"
  sample_name <- sub(pattern, "\\2", file)
  # append to list 
  sample_name_list <- c(sample_name_list, sample_name)
  } else {
   next
  }
}

# Convert the list to a character vector
input_file_vector <- unlist(input_file_list)
reports <- pavian::read_reports(input_file_vector)

### Prep strings for output files
out_file_path <- xargs$output
multisample_comparison_html <- paste0(out_file_path, 'multisample_comparison.html')
summary_table_html <- paste0(out_file_path, 'summary_table.html')
summary_table_csv <- paste0(out_file_path, 'summary_table.CSV')


### this is where the analysis begins
### from pavain https://github.com/fbreitwieser/pavian
merged_reports <- pavian::merge_reports2(reports, col_names = sample_name_list)
summary_report <- pavian::summarize_reports(reports)

tax_data <- merged_reports[[1]]
clade_reads <- merged_reports[[2]]
taxon_reads <- merged_reports[[3]]

colSums(clade_reads,na.rm = T)
colSums(taxon_reads,na.rm = T)

sel_rows <- pavian::filter_taxa(tax_data,
                                rm_clades = c("Chordata", "artificial sequences", "unclassified"),
                                taxRank = "S")
summary(sel_rows)
filtered_clade_reads <- pavian::filter_cladeReads(clade_reads, tax_data, c("Chordata", "artificial sequences", "unclassified"))

tax_data1 <- tax_data[sel_rows,]
filtered_clade_reads1 <- filtered_clade_reads[sel_rows, ]
taxon_reads1 <- taxon_reads[sel_rows, ]

head(cbind(tax_data1[,1:3],clade_reads[sel_rows, ])[order(-apply(filtered_clade_reads1,1,max, na.rm=T)),])
normalized_clade_reads <- normalize(filtered_clade_reads1)
normalized_taxon_reads <- normalize(taxon_reads[sel_rows,], sums = colSums(filtered_clade_reads1,na.rm = T))
head(cbind(tax_data1[,1:3],max=apply(cbind(normalized_clade_reads),1,max, na.rm=T), normalized_clade_reads)[order(-apply(cbind(normalized_clade_reads),1,max, na.rm=T)),])

reads_zscore <- robust_zscore(100*cbind(normalized_clade_reads,normalized_taxon_reads), 0.001)
clade_reads_zscore <- reads_zscore[,1:length(sample_name_list)]
# ^^ could add other functions here 
reads_zscore_df <- cbind(tax_data1[,1:3],max=apply(clade_reads_zscore,1,max, na.rm=T), clade_reads_zscore)[order(-apply(clade_reads_zscore,1,max, na.rm=T)),]

### set up table coloring/stype
# subsample the non-heatmap columns
std_cols <- reads_zscore_df[,c("name", "taxRank", "taxID", "max")]
# subsection sample columns
sample_cols <- reads_zscore_df[, -c(1:4)]

#Set up breaks and colors
brks <- quantile(sample_cols, probs = seq(.05, .95, .05), na.rm = TRUE, names = FALSE)
clrs <- round(seq(180, 60, length.out = length(brks) + 1), 0) %>%
  {paste0("rgb(255",",", .,",", .,")")}

### apply table coloring/style to the data 
multi_sample_table <- DT::datatable(reads_zscore_df,
                filter = "bottom", selection = "single", escape = FALSE,
                caption = "This table shows the taxa at the species level. The color is weighted according to z-score.",
                options = list(
                  autoWidth = TRUE, columnDefs = list(list(width = '100px',targets = 1:5)))) %>%
                formatStyle(names(sample_cols),background = styleInterval(brks,clrs))

### make output files ###
# export multisample_compairision as .HTML
htmlwidgets::saveWidget(multi_sample_table, multisample_comparison_html )
# export report summary as .HTML
summary_report_dt <- DT::datatable(summary_report)
htmlwidgets::saveWidget(summary_report_dt, file = summary_table_html)

### export report summary as .CSV ###
write.csv(summary_report, summary_table_csv, row.names=TRUE)