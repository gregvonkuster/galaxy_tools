#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("Biobase"))
suppressPackageStartupMessages(library("BiocManager"))
suppressPackageStartupMessages(library("GenomeInfoDb"))
suppressPackageStartupMessages(library("GenomicRanges"))
suppressPackageStartupMessages(library("MethylIT"))
suppressPackageStartupMessages(library("Repitools"))
suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option(c("--chromosomes"), action="store", dest="chromosomes", default=NULL, help="Chromosome names to include"),
    make_option(c("--chromosome_names"), action="store", dest="chromosome_names", default=NULL, help="Optionally change chromosome names to those in this list"),
    make_option(c("--context"), action="store", dest="context", type="integer", default=NULL, help="Integer number of the context column in the inputs"),
    make_option(c("--coverage"), action="store", dest="coverage", type="integer", default=NULL, help="Integer number of the coverage column in the inputs"),
    make_option(c("--end"), action="store", dest="end", type="integer", default=NULL, help="Integer number of the end column in the inputs"),
    make_option(c("--fraction"), action="store", dest="fraction", type="integer", default=NULL, help="Integer number of the fraction column in the inputs"),
    make_option(c("--input"), action="store", dest="input", default=NULL, help="Input file"),
    make_option(c("--mC"), action="store", dest="mC", type="integer", default=NULL, help="Integer number of the mC column in the inputs"),
    make_option(c("--output"), action="store", dest="output", default=NULL, help="Output file"),
    make_option(c("--output_log"), action="store", dest="output_log", help="Process log file"),
    make_option(c("--pattern"), action="store", dest="pattern", default=NULL, help="Chromosome name pattern"),
    make_option(c("--percent"), action="store", dest="percent", type="integer", default=NULL, help="Integer number of the percent column in the inputs"),
    make_option(c("--sample_id"), action="store", dest="sample_id", default=NULL, help="Names of the samples corresponding to each file"),
    make_option(c("--seqnames"), action="store", dest="seqnames", type="integer", default=NULL, help="Integer number of the seqnames column in the inputs"),
    make_option(c("--start"), action="store", dest="start", type="integer", default=NULL, help="Integer number of the start column in the inputs"),
    make_option(c("--strand"), action="store", dest="strand", type="integer", default=NULL, help="Integer number of the strand column in the inputs"),
    make_option(c("--uC"), action="store", dest="uC", type="integer", default=NULL, help="Integer number of the uC column in the inputs"),
    make_option(c("--utils"), action="store", dest="utils", help="utils.R script")
)

parser <- OptionParser(usage="%prog [options] file", option_list=option_list);
args <- parse_args(parser, positional_arguments=TRUE);
opt <- args$options;

# Import the shared utility functions.
source(opt$utils);

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

chromosomes <- val_to_null_or_cvector(opt$chromosomes);
chromosome_names <- val_to_null_or_cvector(opt$chromosome_names);
columns <- get_columns(opt$seqnames, opt$start, opt$end, opt$strand, opt$fraction, opt$percent, opt$mC, opt$uC, opt$coverage, opt$context);
sample_id <- val_to_null_or_cvector(opt$sample_id);

# Create the GRanges list.
grange_list <- readCounts2GRangesList(filenames=opt$input,
                                      sample.id=sample_id,
                                      pattern=opt$pattern,
                                      remove=FALSE, 
                                      columns=columns,
                                      chromosome.names=chromosome_names,
                                      chromosomes=chromosomes,
                                      verbose=TRUE);

############
# Debugging.
cat("input: ", opt$input, "\n");
cat("chromosomes: ", toString(chromosomes), "\n");
cat("chromosome_names: ", toString(chromosome_names), "\n");
cat("columns: ", toString(columns), "\n");
cat("sample_id: ", toString(sample_id), "\n");
cat("grange_list: \n");
num_granges <- length(grange_list)[[1]];
for (i in 1:num_granges) {
    grange <- grange_list[[i]];
    show(grange);
    cat("\n\n");
}
############

saveRDS(grange_list, file=opt$output, compress=TRUE);

