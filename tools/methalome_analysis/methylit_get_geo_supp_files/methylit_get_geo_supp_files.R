#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("MethylIT"))
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("RCurl"))

option_list <- list(
    make_option(c("--geo"), action="store", dest="geo", help="Comma-separated list of GEO accession numbers"),
    make_option(c("--pattern"), action="store", dest="pattern", default=NULL, help="pattern for GEO accession files"),
)

parser <- OptionParser(usage="%prog [options] file", option_list=option_list);
args <- parse_args(parser, positional_arguments=TRUE);
opt <- args$options;

# Convert GEO accession numbers into a character vector.
geo_str <- as.character(opt$geo);
geo_list <- strsplit(geo_str, ",")[[1]];
geo <- c(unlist(geo_list, use.names=FALSE));


df <- getGEOSuppFiles(GEO=c("GSM1279513"), makeDirectory=FALSE, baseDir='output_data_dir', pattern=opt$pattern, verbose=TRUE);
cat("\ndf:\n");
df
cat("\ntypeof(df):\n");
typeof(df)

