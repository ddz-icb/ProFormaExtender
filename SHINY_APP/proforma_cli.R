#!/usr/bin/env Rscript

library(optparse)

source("core_calculation.R")

option_list <- list(
  make_option(c("-i", "--input"), type="character", help="Input Excel file (.xlsx)"),
  make_option(c("-o", "--output"), type="character", default=NULL,
              help="Optional output filename (.tsv)")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$input)) {
  stop("Error: Please provide --input FILE.xlsx")
}

df <- run_proforma_pipeline(opt$input)

outfile <- if(!is.null(opt$output)) {
  opt$output
} else {
  paste0(format(Sys.Date(), "%Y%m%d"), "_proforma-ext.output.tsv")
}

write.table(df, outfile, sep="\t", row.names=FALSE, quote=FALSE)

cat("DONE. Output written to:", outfile, "\n")
