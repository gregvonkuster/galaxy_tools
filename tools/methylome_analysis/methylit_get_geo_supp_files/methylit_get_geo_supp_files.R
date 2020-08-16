#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("fs"))
suppressPackageStartupMessages(library("MethylIT"))
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("RCurl"))
suppressPackageStartupMessages(library("R.utils"))

option_list <- list(
    make_option(c("--geo"), action="store", dest="geo", help="Comma-separated list of GEO accession numbers"),
    make_option(c("--output"), action="store", dest="output", defaul=NULL, help="Single output file"),
    make_option(c("--output_dir"), action="store", dest="output_dir", help="Output directory for downloading files"),
    make_option(c("--pattern"), action="store", dest="pattern", default=NULL, help="pattern for GEO accession files"),
    make_option(c("--script_dir"), action="store", dest="script_dir", help="R script source directory")
)

parser <- OptionParser(usage="%prog [options] file", option_list=option_list);
args <- parse_args(parser, positional_arguments=TRUE);
opt <- args$options;

# Import the shared utility functions.
utils_path <- paste(opt$script_dir, "utils.R", sep="/");
source(utils_path);

geo <- string_to_character_vector(opt$geo);

# Download the file(s) to the output directory.
df <- getGEOSuppFiles(GEO=geo, makeDirectory=FALSE, baseDir=opt$output_dir, pattern=opt$pattern, verbose=TRUE);
# If we have a single file, move it to the output.
if (!is.null(opt$output)) {
    for (f in list.files(opt$output_dir, full.names=TRUE)) {
        destname <- gsub("[.]gz$", "", f)
        # We have to gunzip the file since Galaxy doesn't
        # seem to support the tabular.gz format on outputs.
        gunzip(f, destname);
        file_move(destname, opt$output);
    }
}

