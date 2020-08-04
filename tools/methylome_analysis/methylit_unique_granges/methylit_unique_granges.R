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
    make_option(c("--input_dir"), action="store", dest="input_dir", help="Directory of GRanges files"),
    make_option(c("--maxgap"), action="store", dest="maxgap", type="integer", default=NULL, help="Maximum gap for finding interval overlaps between GRange objects"),
    make_option(c("--minoverlap"), action="store", dest="minoverlap", type="integer", default=NULL, help="Minimum gap for finding interval overlaps between GRange objects"),
    make_option(c("--missing"), action="store", dest="missing", type="integer", default=NULL, help="Value to use for missing values in GRange objects"),
    make_option(c("--ncols"), action="store", dest="ncols", type="integer", help="Number of columns to use from the metadata of each GRange"),
    make_option(c("--num_cores"), action="store", dest="num_cores", type="integer", help="The number of cores to use"),
    make_option(c("--output_data_frame"), action="store", dest="output_data_frame", help="Output csv file"),
    make_option(c("--output_grange"), action="store", dest="output_grange", help="Output GRange file"),
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

# Read the list of input GRange files.
input_files <- list.files(path=opt$input_dir, full.names=TRUE);
num_input_files <- length(input_files);

column_list <- list();
grange_list <- list();
hdiv_list <- list();
for (i in 1:num_input_files) {
    input_file <- input_files[[i]];
    column <- sub(pattern="(.*)\\..*$", replacement="\\1", basename(input_file))
    column_list[[i]] <- column;
    grange <- readRDS(input_file);
    grange_list[[i]] <- grange;
    hdiv_list[[i]] <- grange[, "hdiv"];
}

############
# Debugging.
cat("\ninput_files: ", toString(input_files), "\n");
cat("\nopt$ncols: ", ncols, "\n");
cat("\ncolumns: ", toString(columns), "\n");
cat("\nchromosomes: ", toString(chromosomes), "\n");
cat("\nmaxgap: ", maxgap, "\n");
cat("\nminoverlap: ", minoverlap, "\n");
cat("\nmissing: ", missing, "\n");
cat("\noverlap_type: ", opt$overlap_type, "\n");
cat("\nopt$select: ", opt$select, "\n");
cat("\nignore_strand: ", ignore_strand, "\n");
cat("\n\n");
cat("grange_list: \n");
for (i in 1:num_input_files) {
    grange <- grange_list[[i]];
    show(grange);
    cat("\n\n");
}
cat("hdiv_list: \n");
for (i in 1:num_input_files) {
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
colnames(mcols(unique_grange)) <- c(unlist(column_list, use.names=FALSE));

############
# Debugging.
cat("\nunique_grange: \n");
show(unique_grange);
cat("\n\n");
############

# Copy the GRange object to a data frame
# for persisting.
df <- data.frame(matrix(ncol=num_input_files, nrow=nrow(mcols(unique_grange))));
colnames(df) <- c(unlist(column_list, use.names=FALSE));
for (i in 1:num_input_files) {
    column_name <- column_list[[i]];
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

