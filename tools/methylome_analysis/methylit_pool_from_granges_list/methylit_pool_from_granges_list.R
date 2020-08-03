#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("Biobase"))
suppressPackageStartupMessages(library("BiocManager"))
suppressPackageStartupMessages(library("GenomicRanges"))
suppressPackageStartupMessages(library("matrixStats"))
suppressPackageStartupMessages(library("MethylIT"))
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("Repitools"))
suppressPackageStartupMessages(library("S4Vectors"))

option_list <- list(
    make_option(c("--column_number"), action="store", dest="column_number", type="integer", default=1, help="LR column containing probability values"),
    make_option(c("--jstat"), action="store", dest="jstat", default=NULL, help="If stat is jackmean, then any of the stat possible values"),
    make_option(c("--input_data_dir"), action="store", dest="input_data_dir", help="Directory of GRanges files"),
    make_option(c("--num_cores"), action="store", dest="num_cores", help="The number of cores to use"),
    make_option(c("--output"), action="store", dest="output", help="Output Granges converted to data frame file"),
    make_option(c("--prob"), action="store", dest="prob", help="Whether the variable for pooling is between 0 and 1"),
    make_option(c("--stat"), action="store", dest="stat", help="Statistic used to estimate the methylation pool")
)

parser <- OptionParser(usage="%prog [options] file", option_list=option_list);
args <- parse_args(parser, positional_arguments=TRUE);
opt <- args$options;

# Get the list of input data files.
input_data_files <- list.files(path=opt$input_data_dir, full.names=TRUE);
num_input_files <- length(input_data_files)[[1]];

# Load the list of intputs into data frames
grange_list <- list();
for (i in 1:num_input_files) {
    input_data_file <- normalizePath(input_data_files[[i]]);
    grange_list[[i]] <- readRDS(file=input_data_file);
}

# Convert prob to boolean.
if (opt$prob == 'yes') {
    prob <- TRUE;
} else {
    prob <- FALSE;
}

############
# Debugging.
cat("input_data_files: ", toString(input_data_files), "\n");
cat("stat: ", opt$stat, "\n");
cat("prob: ", prob, "\n");
cat("column: ", opt$column_number, "\n");
cat("jstat: ", opt$jstat, "\n");
cat("LR: \n");
for (i in 1:num_input_files) {
    grange <- grange_list[[i]];
    show(grange);
    cat("\n\n");
}
############

# Create a unique GRange from the GRanges list.
unique_grange <- poolFromGRlist(LR=grange_list,
                                stat=opt$stat,
                                num.cores=opt$num_cores,
                                tasks=0L,
                                prob=prob,
                                column=opt$column_number,
                                jstat=opt$jstat,
                                verbose=TRUE);
saveRDS(unique_grange, file=opt$output, compress=TRUE);

