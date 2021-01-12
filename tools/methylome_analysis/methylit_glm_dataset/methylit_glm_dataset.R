#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("Biobase"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("GenomicRanges"))
suppressPackageStartupMessages(library("MethylIT"))
suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option(c("--input_grange"), action="store", dest="input_grange", help="GRanges file"),
    make_option(c("--output_glm"), action="store", dest="output_glm", help="Output glmDataset file")
)

parser <- OptionParser(usage="%prog [options] file", option_list=option_list);
args <- parse_args(parser, positional_arguments=TRUE);
opt <- args$options;

# Here gr is a GRanges object with the
# count matrix in the meta columns.
gr <- readRDS(opt$input_grange);
gr_names <- names(mcols(gr));

############
# Debugging.
cat("\nInput gr:\n\n");
gr
cat("\n\n");
cat("\ngr_names: ", toString(gr_names), "\n\n");
cat("\n\n");
############

# Create the colData data frame, which must have 1 column
# named "condition", which must be a factor with 2 levels.
colData <- data.frame(condition=factor(c(gr_names), levels=c(1, 2)),
                      gr_names,
                      row.names=2);

############
# Debugging.
cat("\ncolData: ",toString(colData), "\n");
cat("\n\n");
############

# Create a RangedGlmDataset.
glm_dataset <- glmDataSet(GR=gr, colData=colData);

############
# Debugging.
cat("\nglm_dataset: \n");
glm_dataset
cat("\n\n");
############

# Save the glm_dataset.
saveRDS(glm_dataset, file=opt$output_glm, compress=TRUE);
