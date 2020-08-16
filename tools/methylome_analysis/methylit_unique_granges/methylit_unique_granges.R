#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("Biobase"))
suppressPackageStartupMessages(library("BiocGenerics"))
suppressPackageStartupMessages(library("BiocManager"))
suppressPackageStartupMessages(library("BiocParallel"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("GenomeInfoDb"))
suppressPackageStartupMessages(library("GenomicRanges"))
suppressPackageStartupMessages(library("IRanges"))
suppressPackageStartupMessages(library("MethylIT"))
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("reshape2"))
suppressPackageStartupMessages(library("S4Vectors"))

option_list <- list(
    make_option(c("--chromosomes"), action="store", dest="chromosomes", default=NULL, help="Comma-separated list of the chromosomes to use from the GRange objects"),
    make_option(c("--columns"), action="store", dest="columns", default=NULL, help="Comma-separated list of numbers of each column to use from the GRange objects"),
    make_option(c("--ignore_strand"), action="store", dest="ignore_strand", help="Ignore strand"),
    make_option(c("--input"), action="store", dest="input", help="Input information divergence estimator file"),
    make_option(c("--maxgap"), action="store", dest="maxgap", type="integer", default=NULL, help="Maximum gap for finding interval overlaps between GRange objects"),
    make_option(c("--minoverlap"), action="store", dest="minoverlap", type="integer", default=NULL, help="Minimum gap for finding interval overlaps between GRange objects"),
    make_option(c("--missing"), action="store", dest="missing", type="integer", default=NULL, help="Value to use for missing values in GRange objects"),
    make_option(c("--ncols"), action="store", dest="ncols", type="integer", help="Number of columns to use from the metadata of each GRange"),
    make_option(c("--num_cores"), action="store", dest="num_cores", type="integer", help="The number of cores to use"),
    make_option(c("--output_data_frame"), action="store", dest="output_data_frame", help="Output csv file"),
    make_option(c("--output_grange"), action="store", dest="output_grange", help="Output GRange file"),
    make_option(c("--select"), action="store", dest="select", help="Select value"),
    make_option(c("--overlap_type"), action="store", dest="overlap_type", help="Overlap type"),
    make_option(c("--script_dir"), action="store", dest="script_dir", help="R script source directory")
)

parser <- OptionParser(usage="%prog [options] file", option_list=option_list);
args <- parse_args(parser, positional_arguments=TRUE);
opt <- args$options;

# Import the shared utility functions.
utils_path <- paste(opt$script_dir, "utils.R", sep="/");
source(utils_path);

if (is.null(opt$chromosomes)) {
    chromosomes = NULL;
} else {
    chromosomes <- string_to_character_vector(opt$chromosomes);
}

if (is.null(opt$columns)) {
    columns <- NULL;
} else {
    columns <- string_to_integer_or_vector(opt$columns);
}

ignore_strand <- string_to_boolean(opt$ignore_strand, default=FALSE);

if (is.null(opt$maxgap)) {
    maxgap <- -1;
} else {
    maxgap <- opt$maxgap;
}

if (is.null(opt$minoverlap)) {
    minoverlap <- 1;
} else {
    minoverlap <- opt$minoverlap;
}

if (is.null(opt$missing) || opt$missing == '0') {
    missing <- 0;
} else {
    # The value of missing will be NA here.
    missing <- NA;
}

if (opt$ncols == 0) {
    ncols <- NULL;
} else {
    ncols <- opt$ncols;
}

inf_div <- readRDS(opt$input);
inf_div_names <- names(inf_div);
############
# Debugging.
cat("\ninf_div_names: ", toString(inf_div_names), "\n");
############
num_inf_divs <- length(inf_div_names)[[1]];
hdiv_list <- list();
for (i in 1:num_inf_divs) {
    grange <- inf_div[[i]];
    hdiv_list[[i]] <- grange[, "hdiv"];
}

############
# Debugging.
cat("num_inf_divs: ", num_inf_divs, "\n");
cat("ncols: ", ncols, "\n");
cat("columns: ", toString(columns), "\n");
cat("chromosomes: ", toString(chromosomes), "\n");
cat("maxgap: ", maxgap, "\n");
cat("minoverlap: ", minoverlap, "\n");
cat("missing: ", missing, "\n");
cat("overlap_type: ", opt$overlap_type, "\n");
cat("opt$select: ", opt$select, "\n");
cat("ignore_strand: ", ignore_strand, "\n");
cat("granges:\n");
for (i in 1:num_inf_divs) {
    grange <- inf_div[[i]];
    show(grange);
    cat("\n\n");
}
cat("inf_divs: \n");
for (i in 1:num_inf_divs) {
    hdiv <- hdiv_list[[i]];
    show(hdiv);
    cat("\n\n");
}
############

# Get the unique GRange object from the list of GRange objects.
unique_grange <- uniqueGRanges(hdiv_list,
                               ncols=ncols,
                               columns=columns,
                               chromosomes=chromosomes,
                               maxgap=maxgap,
                               minoverlap=minoverlap,
                               missing=missing,
                               type=opt$overlap_type,
                               select=opt$select,
                               ignore.strand=ignore_strand,
                               num.cores=opt$num_cores,
                               verbose=TRUE);
colnames(mcols(unique_grange)) <- c(unlist(inf_div_names, use.names=FALSE));

############
# Debugging.
cat("\nunique_grange: \n");
show(unique_grange);
cat("\n\n");
############

# Copy the GRange object to a data frame for persisting.
df <- data.frame(matrix(ncol=num_inf_divs, nrow=nrow(mcols(unique_grange))));
colnames(df) <- c(unlist(inf_div_names, use.names=FALSE));
for (i in 1:num_inf_divs) {
    column_name <- inf_div_names[[i]];
    column_index <- grep(column_name, colnames(df))[[1]];
    ############
    # Debugging.
    cat("\ncolumn_index: ", column_index, "\n");
    cat("column_name: ", column_name, "\n\n");
    ############
    df[,column_index] <- mcols(unique_grange)[[column_name]];
}

############
# Debugging.
cat("\nInitial data frame after GRange conversion: \n");
cat("dim(df):\n");
dim(df);
cat("str(df):\n");
str(df);
cat("summary(df):\n");
summary(df);
cat("colnames(df):\n");
colnames(df);
cat("head(df):\n");
head(df);
cat("\n\n");
############

df <- melt(df);
colnames(df) <- c("Category", "HellingerDivergence")
df <- df[df$HellingerDivergence > 0, ]

############
# Debugging.
cat("\nData frame after melting: \n");
cat("dim(df):\n");
dim(df);
cat("str(df):\n");
str(df);
cat("summary(df):\n");
summary(df);
cat("colnames(df):\n");
colnames(df);
cat("head(df):\n");
head(df);
cat("\n\n");
############

# Save the unique GRange object.
saveRDS(grange, file=opt$output_grange, compress=TRUE);

# Save the data frame.
write.table(df, file=opt$output_data_frame, sep='\t', quote=FALSE, row.names=FALSE);

