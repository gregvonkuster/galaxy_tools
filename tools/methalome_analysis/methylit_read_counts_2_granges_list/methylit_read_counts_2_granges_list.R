#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("Biobase"))
suppressPackageStartupMessages(library("BiocManager"))
suppressPackageStartupMessages(library("GenomeInfoDb"))
suppressPackageStartupMessages(library("GenomicRanges"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("MethylIT"))
suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option(c("--chromosomes"), action="store", dest="chromosomes", default=NULL, help="Chromosome names to include"),
    make_option(c("--chromosome_names"), action="store", dest="chromosome_names", default=NULL, help="Optionally change chromosome names to those in this list"),
    make_option(c("--context"), action="store", dest="context", type="integer", default=NULL, help="Integer number of the context column in the inputs"),
    make_option(c("--coverage"), action="store", dest="coverage", type="integer", default=NULL, help="Integer number of the coverage column in the inputs"),
    make_option(c("--end"), action="store", dest="end", type="integer", default=NULL, help="Integer number of the end column in the inputs"),
    make_option(c("--fraction"), action="store", dest="fraction", type="integer", default=NULL, help="Integer number of the fraction column in the inputs"),
    make_option(c("--input_data_dir"), action="store", dest="input_data_dir", help="Directory containing the input files"),
    make_option(c("--mC"), action="store", dest="mC", type="integer", default=NULL, help="Integer number of the mC column in the inputs"),
    make_option(c("--output_log"), action="store", dest="output_log", help="Process log file"),
    make_option(c("--pattern"), action="store", dest="pattern", default=NULL, help="Chromosome name pattern"),
    make_option(c("--percent"), action="store", dest="percent", type="integer", default=NULL, help="Integer number of the percent column in the inputs"),
    make_option(c("--sample_id"), action="store", dest="sample_id", default=NULL, help="Names of the samples corresponding to each file"),
    make_option(c("--seqnames"), action="store", dest="seqnames", type="integer", default=NULL, help="Integer number of the seqnames column in the inputs"),
    make_option(c("--start"), action="store", dest="start", type="integer", default=NULL, help="Integer number of the start column in the inputs"),
    make_option(c("--strand"), action="store", dest="strand", type="integer", default=NULL, help="Integer number of the strand column in the inputs"),
    make_option(c("--uC"), action="store", dest="uC", type="integer", default=NULL, help="Integer number of the uC column in the inputs")
)

parser <- OptionParser(usage="%prog [options] file", option_list=option_list);
args <- parse_args(parser, positional_arguments=TRUE);
opt <- args$options;

get_columns <- function(seqnames, start, end, strand, fraction, percent, mC, uC, coverage, context) {
    columns <- integer();
    if (!is.null(seqnames)) {
        columns['seqnames'] <- seqnames;
    }
    if (!is.null(start)) {
        columns['start'] <- start;
    }
    if (!is.null(end)) {
        columns['end'] <- end;
    }
    if (!is.null(strand)) {
        columns['strand'] <- strand;
    }
    if (!is.null(fraction)) {
        columns['fraction'] <- fraction;
    }
    if (!is.null(percent)) {
        columns['percent'] <- percent;
    }
    if (!is.null(mC)) {
        columns['mC'] <- mC;
    }
    if (!is.null(uC)) {
        columns['uC'] <- uC;
    }
    if (!is.null(coverage)) {
        columns['coverage'] <- coverage;
    }
    if (!is.null(context)) {
        columns['context'] <- context;
    }
    return (columns)
}

# Get the list of input data files.
input_data_files <- list.files(path=opt$input_data_dir, full.names=TRUE);
#f = getGEOSuppFiles(GEO=c("GSM1279513"), makeDirectory=FALSE, baseDir='output_data_dir', pattern=NULL, verbose=TRUE);

if (is.null(opt$chromosomes)) {
    chromosomes = NULL;
} else {
    # Convert chromosomes into a character vector.
    chromosomes_str <- as.character(opt$chromosomes);
    chromosomes_list <- strsplit(chromosomes_str, ",")[[1]];
    chromosomes <- c(unlist(chromosomes_list, use.names=FALSE));
}

if (is.null(opt$chromosome_names)) {
    chromosome_names = NULL;
} else {
    # Convert chromosome_names into a character vector.
    chromosome_names_str <- as.character(opt$chromosome_names);
    chromosome_names_list <- strsplit(chromosome_names_str, ",")[[1]];
    chromosome_names <- c(unlist(chromosome_names_list, use.names=FALSE));
}

columns <- get_columns(opt$seqnames, opt$start, opt$end, opt$strand, opt$fraction, opt$percent, opt$mC, opt$uC, opt$coverage, opt$context);

if (is.null(opt$sample_id)) {
    sample_id = NULL;
} else {
    # Convert sample_id into a character vector.
    sample_id_str <- as.character(opt$sample_id);
    sample_id_list <- strsplit(sample_id_str, ",")[[1]];
    sample_id <- c(unlist(sample_id_list, use.names=FALSE));
}

# Create the GRangesList data frame.
meth_list <- readCounts2GRangesList(filenames=input_data_files,
                                    sample.id=sample_id,
                                    pattern=opt$pattern,
                                    remove=FALSE, 
                                    columns=columns,
                                    chromosome.names=chromosome_names,
                                    chromosomes=chromosomes,
                                    verbose=TRUE);
num_dfs <- length(meth_list)[[1]];

for (i in 1:num_dfs) {
    df <- meth_list[i];
    # Save the data frame.
    file_name <- paste(sample_id[i], ".csv", sep="");
    file_path = paste("output_data_dir", file_name, sep="/");
    write.csv(df, file=file_path, row.names=F);
}

