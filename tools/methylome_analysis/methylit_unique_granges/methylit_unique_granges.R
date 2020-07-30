#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("Biobase"))
suppressPackageStartupMessages(library("BiocGenerics"))
suppressPackageStartupMessages(library("BiocManager"))
suppressPackageStartupMessages(library("BiocParallel"))
suppressPackageStartupMessages(library("GenomeInfoDb"))
suppressPackageStartupMessages(library("GenomicRanges"))
suppressPackageStartupMessages(library("IRanges"))
suppressPackageStartupMessages(library("MethylIT"))
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("S4Vectors"))

option_list <- list(
    make_option(c("--chromosomes"), action="store", dest="chromosomes", default=NULL, help="Comma-separated list of the chromosomes to use from the GRange objects"),
    make_option(c("--columns"), action="store", dest="columns", default=NULL, help="Comma-separated list of numbers of each column to use from the GRange objects"),
    make_option(c("--ignore_strand"), action="store", dest="ignore_strand", help="Ignore strand"),
    make_option(c("--input_dir"), action="store", dest="input_dir", help="Directory of GRanges files"),
    make_option(c("--maxgap"), action="store", dest="maxgap", type="integer", default=NULL, help="Maximum gap for finding interval overlaps between GRange objects"),
    make_option(c("--minoverlap"), action="store", dest="minoverlap", type="integer", default=NULL, help="Minimum gap for finding interval overlaps between GRange objects"),
    make_option(c("--missing"), action="store", dest="missing", type="integer", default=NULL, help="Value to use for missing values in GRange objects"),
    make_option(c("--ncols"), action="store", dest="ncols", type="integer", help="Number of columns to use from the metadata of each GRange"),
    make_option(c("--num_cores"), action="store", dest="num_cores", type="integer", help="The number of cores to use"),
    make_option(c("--output"), action="store", dest="output", help="Output file"),
    make_option(c("--select"), action="store", dest="select", help="Select value"),
    make_option(c("--overlap_type"), action="store", dest="overlap_type", help="Overlap type")
)

parser <- OptionParser(usage="%prog [options] file", option_list=option_list);
args <- parse_args(parser, positional_arguments=TRUE);
opt <- args$options;

if (is.null(opt$chromosomes)) {
    chromosomes = NULL;
} else {
    # Convert chromosomes into a character vector.
    chromosomes_str <- as.character(opt$chromosomes);
    chromosomes_list <- strsplit(chromosomes_str, ",")[[1]];
    chromosomes <- c(unlist(chromosomes_list, use.names=FALSE));
}

# Convert columns to an integer or integer vector if not NULL.
if (!is.null(opt$columns)) {
    columns_list <- strsplit(opt$columns, ",")[[1]];
    num_columns <- length(columns_list)[[1]];
    if (num_columns == 1) {
        columns <- columns_list[[1]];
    } else {
        columns <- integer(num_columns);
        for (i in 1:num_columns) {
            column[[i]] <- columns_list[[i]];
        }
    }
} else {
    columns <- NULL;
}

# Convert ignore_strand to boolean.
if (opt$ignore_strand == 'no') {
    ignore_strand <- FALSE;
} else {
    ignore_strand <- TRUE;
}

if (is.null(opt$maxgap)) {
    max_gap <- -1;
} else {
    max_gap <- opt$max_gap;
}

if (is.null(opt$minoverlap)) {
    minoverlap <- 1;
} else {
    minoverlap <- opt$minoverlap;
}

if (is.null(opt$missing)) {
    missing <- 0;
} else {
    missing <- opt$missing;
}

# Read the list of input GRange files.
input_files <- list.files(path=opt$input_dir, full.names=TRUE);
num_input_files <- length(input_files);

grange_list <- list();
for (i in 1:num_input_files) {
    input_file <- normalizePath(input_files[[i]]);
    grange_list[[i]] <- readRDS(input_file);
}

# Get the unique GRange object from the list of GRange objects.
grange <- uniqueGRanges(grange_list,
                        ncols=opt$ncols,
                        columns=columns,
                        chromosomes=chromosomes,
                        maxgap=max_gap,
                        minoverlap=minoverlap,
                        missing=missing,
                        type=opt$overlap_type,
                        select=opt$select,
                        ignore.strand=ignore_strand,
                        num.cores=opt$num_cores,
                        verbose=TRUE);

saveRDS(grange, file=opt$output, compress=TRUE);

