#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("Biobase"))
suppressPackageStartupMessages(library("BiocManager"))
suppressPackageStartupMessages(library("GenomicRanges"))
suppressPackageStartupMessages(library("MethylIT"))
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("xtable"))

option_list <- list(
    make_option(c("--absolute"), action="store", dest="absolute", help="Flag to transform total variation (TV) to |TV|"),
    make_option(c("--div_col"), action="store", dest="div_col", type="integer", default=NULL, help="GRange column number containing the divergence variable for which the estimation of the cutpoint will be performed"),
    make_option(c("--input_cutpoint"), action="store", dest="input_cutpoint", help="File containing an estimated cutpoint"),
    make_option(c("--input_pdmp"), action="store", dest="input_pdmp", help="File containing GRange objects that include selected cytosine sites and specified divergence probabilities"),
    make_option(c("--output_crc"), action="store", dest="output_crc", help="Output cytosince read counts file"),
    make_option(c("--output_pdmpdmp"), action="store", dest="output_pdmpdmp", help="Output pdmp file"),
    make_option(c("--pval_col"), action="store", dest="pval_col", type="integer", default=NULL, help="Index of the GRanges column containing the p-values"),
    make_option(c("--script_dir"), action="store", dest="script_dir", help="R script source directory"),
    make_option(c("--tv_col"), action="store", dest="tv_col", type="integer", default=NULL, help="Index of the GRanges column containing the total variation for filtering cytosine positions"),
    make_option(c("--tv_cut"), action="store", dest="tv_cut", type="double", default=0.25, help="Total variation cutoff")
)

parser <- OptionParser(usage="%prog [options] file", option_list=option_list);
args <- parse_args(parser, positional_arguments=TRUE);
opt <- args$options;

# Import the shared utility functions.
utils_path <- paste(opt$script_dir, "utils.R", sep="/");
source(utils_path);

absolute <- string_to_boolean(opt$absolute, default=FALSE);
l <- readRDS(opt$input_cutpoint);
cutpoint <- l$cutpoint;
LR <- readRDS(opt$input_pdmp);

############
# Debugging.
cat("\nabsolute: ", absolute, "\n");
cat("\ncutpoint:\n");
cutpoint
cat("\nopt$div_col: ", opt$div_col, "\n");
cat("\nopt$pval_col: ", opt$pval_col, "\n");
cat("\nopt$tv_col: ", opt$tv_col, "\n");
cat("\nopt$tv_cut: ", opt$tv_cut, "\n");
cat("\n\n");
############

pDMP <- selectDIMP(LR=LR,
                   div.col=opt$div_col,
                   pval.col=opt$pval_col,
                   absolute=absolute,
                   cutpoint=cutpoint,
                   tv.col=opt$tv_col,
                   tv.cut=opt$tv_cut);

############
# Debugging.
cat("pDMP:\n");
pDMP
cat("\n\n");
############

# Save the potential_methylation_signal.
saveRDS(pDMP, file=opt$output_pdmpdmp, compress=TRUE);

# Output statistics.
mrs_df <- get_methylated_read_statistics(pDMP);
output_statistics(opt$output_crc, mrs_df=mrs_df);

