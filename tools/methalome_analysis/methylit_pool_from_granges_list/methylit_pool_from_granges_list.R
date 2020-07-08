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
num_input_files <- length(input_data_files);

granges_list <- list();
for (i in 1:num_input_files) {
    input_data_file <- input_data_files[i];
    df = read.csv(file=input_data_file, header=T, strip.white=TRUE, stringsAsFactors=FALSE, sep=",");
    granges_list[[i]] <- df;
}

# Convert prob to boolean.
if (opt$prob == 'yes') {
    prob = TRUE;
} else {
    prob = FALSE;
}

# Create a unique GRange object from the GRanges list.
unique_grange <- poolFromGRlist(LR=granges_list,
                                stat=opt$stat,
                                num.cores=opt$num_cores,
                                tasks=0L,
                                prob=prob,
                                column=opt$column_number,
                                jstat=opt$jstat,
                                verbose=TRUE);
# Convert the GRange object to a data frame for saving
# to a file.
df <- annoGR2DF(unique_grange);
write.csv(df, file=opt$output, row.names=F);

