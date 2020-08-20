#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("Biobase"))
suppressPackageStartupMessages(library("BiocManager"))
suppressPackageStartupMessages(library("BiocParallel"))
suppressPackageStartupMessages(library("caret"))
suppressPackageStartupMessages(library("GenomicRanges"))
suppressPackageStartupMessages(library("MASS"))
suppressPackageStartupMessages(library("MethylIT"))
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("xtable"))

option_list <- list(
    make_option(c("--center"), action="store", dest="center", default=NULL, help="Flag to shift variables to be zero-centered"),
    make_option(c("--classifier"), action="store", dest="classifier", help="A comma-separated string of classification models"),
    make_option(c("--column"), action="store", dest="column", help="Column names in the input for the control samples"),
    make_option(c("--control_names"), action="store", dest="control_names", help="Column names in the input for the control samples"),
    make_option(c("--input"), action="store", dest="input", help="File containing a pDMP object"),
    make_option(c("--interactions"), action="store", dest="interactions", default=NULL, help="Interactions to consider in the logistics regression model"),
    make_option(c("--n_pc"), action="store", dest="n_pc", type="integer", help="Number of principal components to use if the classifier is not 'logistic'"),
    make_option(c("--num_boot"), action="store", dest="num_boot", type="integer", help="Number of bootstrap validations to perform in the evaluation of the logistic regression 'group versus divergence' at DMPs"),
    make_option(c("--num_cores"), action="store", dest="num_cores", type="integer", help="Number of processors to use"),
    make_option(c("--output_html"), action="store", dest="output_html", help="Output html file"),
    make_option(c("--output_type"), action="store", dest="output_type", help="Output type"),
    make_option(c("--prop"), action="store", dest="prop", type="double", help="Proportion to split the dataset used in the logistic regression into two subsets, training and testing"),
    make_option(c("--pval_column"), action="store", dest="pval_column", type="integer", default=NULL, help="Index of the GRanges column containing the p-values"),
    make_option(c("--scale"), action="store", dest="scale", help="Flag to scale variables to have unit variance before the analysis takes place"),
    make_option(c("--script_dir"), action="store", dest="script_dir", help="R script source directory"),
    make_option(c("--seed"), action="store", dest="seed", type="integer", help="Seed value for random number generation"),
    make_option(c("--treatment_names"), action="store", dest="treatment_names", help="Column names in the input for the treatment samples")
)

parser <- OptionParser(usage="%prog [options] file", option_list=option_list);
args <- parse_args(parser, positional_arguments=TRUE);
opt <- args$options;

# Import the shared utility functions.
utils_path <- paste(opt$script_dir, "utils.R", sep="/");
source(utils_path);

center <- string_to_boolean(opt$center, default=FALSE);
classifier <- string_to_character_vector(opt$classifier);
# Convert opt$column to a named logical vector.
names <- string_to_character_vector(opt$column);
num_names <- length(names);
column <- logical(length=num_names);
for (i in 1:num_names) {
    column[[i]] <- TRUE;
}
names(column) <- names;
control_names <- string_to_character_vector(opt$control_names);
LR <- readRDS(opt$input);
scale <- string_to_boolean(opt$scale, default=FALSE);
treatment_names <- string_to_character_vector(opt$treatment_names);

############
# Debugging.
cat("\ncenter: ", center, "\n");
cat("\nclassifier: ", toString(classifier), "\n");
cat("\ncolumn: ", toString(column), "\n");
cat("\ncontrol_names: ", toString(control_names), "\n");
cat("\nopt$input: ", opt$input, "\n");
cat("\nopt$interactions: ", opt$interactions, "\n");
cat("\nopt$n_pc: ", opt$n_pc, "\n");
cat("\nopt$num_boot: ", opt$num_boot, "\n");
cat("\nopt$output_type: ", opt$output_type, "\n");
cat("\nopt$prop: ", opt$prop, "\n");
cat("\nopt$pval_column: ", opt$pval_column, "\n");
cat("\nscale: ", scale, "\n");
cat("\nseed: ", opt$seed, "\n");
cat("\ntreatment_names: ", toString(treatment_names), "\n");
cat("\n\n");
############

performance <- evaluateDIMPclass(LR,
                                 control.names=control_names,
                                 treatment.names=treatment_names,
                                 column=column,
                                 classifier=classifier,
                                 pval.col=opt$pval_column,
                                 n.pc=opt$n_pc,
                                 center=center,
                                 scale=scale,
                                 interactions=opt$interactions,
                                 output=opt$output_type,
                                 prop=opt$prop,
                                 num.boot=opt$num_boot,
                                 num.cores=opt$num_cores,
                                 seed=opt$seed,
                                 verbose=TRUE);

############
# Debugging.
cat("performance:\n");
performance
cat("\ntypeof(performance): ", typeof(performance), "\n");
cat("\n\n");
############

# Output the performance statistics.
output_statistics(opt$output_html, perf_df=performance);

