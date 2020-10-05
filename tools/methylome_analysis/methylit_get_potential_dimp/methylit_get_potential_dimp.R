#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("Biobase"))
suppressPackageStartupMessages(library("BiocManager"))
suppressPackageStartupMessages(library("GenomicRanges"))
suppressPackageStartupMessages(library("MethylIT"))
suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option(c("--absolute"), action="store", dest="absolute", help="If TRUE, then TV is transformed into |TV|"),
    make_option(c("--alpha"), action="store", dest="alpha", type="double", help="Value to select cytosine sites with information divergence for which the probabilities hold"),
    make_option(c("--dist_name"), action="store", dest="dist_name", default=NULL, help="A comma-separated string naming the models to fit"),
    make_option(c("--div_column"), action="store", dest="div_column", type="integer", help="Index of the GRanges column where the information divergence is given"),
    make_option(c("--gof_report"), action="store", dest="gof_report", default=NULL, help="File containing a Goodness of Fit (GOF) report"),
    make_option(c("--hdiv_col"), action="store", dest="hdiv_col", type="integer", default=NULL, help="Index of the GRanges column containing the Hellinger distance for filtering cytosince positions"),
    make_option(c("--hdiv_cut"), action="store", dest="hdiv_cut", type="double", default=NULL, help="Hellinger distance for filtering cytosince positions cutoff"),
    make_option(c("--input"), action="store", dest="input", help="File containing an information divergence estimator or a GRange object that includes Fisher's columns"),
    make_option(c("--min_coverage"), action="store", dest="min_coverage", type="integer", default=NULL, help="Minimum coverage for cysosine sites"),
    make_option(c("--output_potdimp"), action="store", dest="output_potdimp", help="Output potdimp file"),
    make_option(c("--padjustmethod"), action="store", dest="padjustmethod", default=NULL, help="Method for adjusting the p-values"),
    make_option(c("--pval_column"), action="store", dest="pval_column", type="integer", default=NULL, help="Index of the GRanges column containing the p-values"),
    make_option(c("--script_dir"), action="store", dest="script_dir", help="R script source directory"),
    make_option(c("--tv_col"), action="store", dest="tv_col", type="integer", default=NULL, help="Index of the GRanges column containing the total variation for filtering cytosine positions"),
    make_option(c("--tv_cut"), action="store", dest="tv_cut", type="double", default=NULL, help="Total variation cutoff")
)

parser <- OptionParser(usage="%prog [options] file", option_list=option_list);
args <- parse_args(parser, positional_arguments=TRUE);
opt <- args$options;

cat("args:\n", toString(args), "\n\n");

# Import the shared utility functions.
utils_path <- paste(opt$script_dir, "utils.R", sep="/");
source(utils_path);

absolute <- string_to_boolean(opt$absolute, default=FALSE);

if (!is.null(opt$gof_report)) {
    gof_report <- readRDS(opt$gof_report);
    dist_name <- gof_report$bestModel;
    nlms <- gof_report$nlms;
} else {
    nlms <- NULL;
    if (is.null(opt$dist_name)) {
        dist_name <- NULL;
    } else {
        dist_name <- string_to_character_vector(opt$dist_name);
    }
}

LR <- readRDS(opt$input);

############
# Debugging.
cat("\nabsolute: ", absolute, "\n");
cat("\nopt$alpha: ", opt$alpha, "\n");
cat("\ndist_name: ", dist_name, "\n");
cat("\nopt$div_column: ", opt$div_column, "\n");
cat("\nopt$gof_report: ", opt$gof_report, "\n");
cat("\nopt$hdiv_col: ", opt$hdiv_col, "\n");
cat("\nopt$hdiv_cut: ", opt$hdiv_cut, "\n");
cat("\nnlms:\n");
nlms
cat("\n\n");
cat("\nopt$input: ", opt$input, "\n");
cat("\nopt$min_coverage: ", opt$min_coverage, "\n");
cat("\nopt$padjustmethod: ", opt$padjustmethod, "\n");
cat("\nopt$pval_column: ", opt$pval_column, "\n");
cat("\nopt$tv_col: ", opt$tv_col, "\n");
cat("\nopt$tv_cut: ", opt$tv_cut, "\n");
cat("\n\n");
############

potential_methylation_signal <- getPotentialDIMP(LR,
                                                 nlms=nlms,
                                                 div.col=opt$div_column,
                                                 dist.name=dist_name,
                                                 absolute=absolute,
                                                 alpha=opt$alpha,
                                                 pval.col=opt$pval_column,
                                                 tv.col=opt$tv_col,
                                                 tv.cut=opt$tv_cut,
                                                 min.coverage=opt$min_coverage,
                                                 hdiv.col=opt$hdiv_col,
                                                 hdiv.cut=opt$hdiv_cut,
                                                 pAdjustMethod=opt$padjustmethod);

############
# Debugging.
cat("\npotential_methylation_signal:\n");
potential_methylation_signal
cat("\n\n");
############

# Save the potential_methylation_signal.
saveRDS(potential_methylation_signal, file=opt$output_potdimp, compress=TRUE);

