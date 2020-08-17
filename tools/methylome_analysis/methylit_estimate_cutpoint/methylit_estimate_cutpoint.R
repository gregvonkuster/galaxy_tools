#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("Biobase"))
suppressPackageStartupMessages(library("BiocManager"))
suppressPackageStartupMessages(library("caret"))
suppressPackageStartupMessages(library("GenomicRanges"))
suppressPackageStartupMessages(library("MethylIT"))
suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option(c("--clas_perf"), action="store", dest="clas_perf", help="Flag to evaluate the classification performance for the estimated cutpoint using a model classifier"),
    make_option(c("--classifier1"), action="store", dest="classifier1", help="A comma-separated string of classification models"),
    make_option(c("--classifier2"), action="store", dest="classifier2", default=NULL, help="A comma-separated string of classification models"),
    make_option(c("--column"), action="store", dest="column", help="Column names in the input for the control samples"),
    make_option(c("--control_names"), action="store", dest="control_names", help="Column names in the input for the control samples"),
    make_option(c("--cut_values_by"), action="store", dest="cut_values_by", type="integer", default=NULL, help="Increment of the sequence"),
    make_option(c("--cut_values_from"), action="store", dest="cut_values_from", type="integer", default=NULL, help="Starting point of the sequence"),
    make_option(c("--cut_values_to"), action="store", dest="cut_values_to", type="integer", default=NULL, help="Ending point of the sequence"),
    make_option(c("--cutp_data"), action="store", dest="cutp_data", help="Flag to output a data frame for further analysis or estimation of the optimal cutpoint based  only on the selected divergence"),
    make_option(c("--div_col"), action="store", dest="div_col", type="integer", default=NULL, help="GRange column number containing the divergence variable for which the estimation of the cutpoint will be performed"),
    make_option(c("--input"), action="store", dest="input", help="File containing GRange objects that include selected cytosine sites and specified divergence probabilities"),
    make_option(c("--interactions"), action="store", dest="interactions", default=NULL, help="Interactions to consider in the logistics regression model"),
    make_option(c("--n_pc"), action="store", dest="n_pc", type="integer", help="Number of principal components to use if the classifier is not 'logistic'"),
    make_option(c("--num_cores"), action="store", dest="num_cores", type="integer", help="Number of processors to use"),
    make_option(c("--output_rdata"), action="store", dest="output_rdata", help="Output rdata file"),
    make_option(c("--post_cut"), action="store", dest="post_cut", type="double", default=0.5, help="Posterior probability to decide whether a DMP belongs to treatment group"),
    make_option(c("--prop"), action="store", dest="prop", type="double", help="Proportion to split the dataset used in the logistic regression into two subsets, training and testing"),
    make_option(c("--script_dir"), action="store", dest="script_dir", help="R script source directory"),
    make_option(c("--simple"), action="store", dest="simple", default=NULL, help="Flag to use the Youden Index to estimate the cutpoint"),
    make_option(c("--stat"), action="store", dest="stat", type="integer", default=1, help="Set the number indicating the statistic to be used in the testing testing"),
    make_option(c("--treatment_names"), action="store", dest="treatment_names", help="Column names in the input for the treatment samples"),
    make_option(c("--tv_col"), action="store", dest="tv_col", type="integer", default=NULL, help="Index of the GRanges column containing the total variation for filtering cytosine positions"),
    make_option(c("--tv_cut"), action="store", dest="tv_cut", type="double", default=0.25, help="Total variation cutoff")
)

parser <- OptionParser(usage="%prog [options] file", option_list=option_list);
args <- parse_args(parser, positional_arguments=TRUE);
opt <- args$options;

# Import the shared utility functions.
utils_path <- paste(opt$script_dir, "utils.R", sep="/");
source(utils_path);

clas_perf <- string_to_boolean(opt$clas_perf, default=FALSE);
classifier1 <- string_to_character_vector(opt$classifier1);
if (is.null(opt$classifier2)) {
    classifier2 <- NULL;
} else {
    classifier2 <- string_to_character_vector(opt$classifier2);
}

# Convert opt$column to a named logical vector.
names <- string_to_character_vector(opt$column);
num_names <- length(names);
column <- logical(length=num_names);
for (i in 1:num_names) {
    column[[i]] <- TRUE;
}
names(column) <- names;

control_names <- string_to_character_vector(opt$control_names);
if (is.null(opt$cut_values_by) && is.null(opt$cut_values_from) && is.null(opt$cut_values_to)) {
    cut_values <- NULL;
} else {
    cut_values <- seq(opt$cut_values_from, opt$cut_values_to, opt$cut_values_by);
}
cutp_data <- string_to_boolean(opt$cutp_data, default=FALSE);
LR <- readRDS(opt$input);
simple <- string_to_boolean(opt$simple, default=TRUE);
treatment_names <- string_to_character_vector(opt$treatment_names);


############
# Debugging.
cat("\nclas_perf: ", clas_perf, "\n");
cat("\nclassifier1: ", toString(classifier1), "\n");
cat("\nclassifier2: ", toString(classifier2), "\n");
cat("\ncolumn: ", toString(column), "\n");
cat("\ncontrol_names: ", toString(control_names), "\n");
cat("\nopt$cut_values_by: ", opt$cut_values_by, "\n");
cat("\nopt$cut_values_from: ", opt$cut_values_from, "\n");
cat("\nopt$cut_values_to: ", opt$cut_values_to, "\n");
cat("\ncutp_data: ", cutp_data, "\n");
cat("\nopt$div_col: ", opt$div_col, "\n");
cat("\nopt$input: ", opt$input, "\n");
cat("\nopt$interactions: ", opt$interactions, "\n");
cat("\nopt$n_pc: ", opt$n_pc, "\n");
cat("\nopt$post_cut: ", opt$post_cut, "\n");
cat("\nopt$prop: ", opt$prop, "\n");
cat("\nsimple: ", toString(simple), "\n");
cat("\nopt$stat: ", opt$stat, "\n");
cat("\ntreatment_names: ", toString(treatment_names), "\n");
cat("\nopt$tv_col: ", opt$tv_col, "\n");
cat("\nopt$tv_cut: ", opt$tv_cut, "\n");
cat("\n\n");
############

l <- estimateCutPoint(LR,
                      control.names=control_names,
                      treatment.names=treatment_names,
                      simple=simple,
                      column=column,
                      classifier1=classifier1,
                      classifier2=classifier2,
                      tv.cut=opt$tv_cut,
                      tv.col=opt$tv_col,
                      div.col=opt$div_col,
                      clas.perf=clas_perf,
                      post.cut=opt$post_cut,
                      prop=opt$prop,
                      n.pc=opt$n_pc,
                      interactions=opt$interactions,
                      cut.values=cut_values,
                      stat=opt$stat,
                      cutp_data=cutp_data,
                      num.cores=opt$num_cores);

############
# Debugging.
cat("l$cutpoint:\n");
l$cutpoint
cat("\n\n");
cat("l$testSetPerformance:\n");
l$testSetPerformance
cat("\n\n");
cat("l$testSetModel.FDR:\n");
l$testSetModel.FDR
cat("\n\n");
cat("l$model:\n");
l$model
cat("\n\n");
cat("l$modelConfMatrix:\n");
l$modelConfMatrix
cat("\n\n");
cat("l$initModel:\n");
l$initModel
cat("\n\n");
cat("l$postProbCut:\n");
l$postProbCut
cat("\n\n");
cat("l$postCut:\n");
l$postCut
cat("\n\n");
cat("l$classifier:\n");
l$classifier
cat("\n\n");
cat("l$classifier:\n");
l$classifier
cat("\n\n");
cat("l$statistic:\n");
l$statistic
cat("\n\n");
cat("l$optStatVal:\n");
l$optStatVal
cat("\n\n");
cat("l$cutpData:\n");
l$ocutpData
cat("\n\n");
############

# Save the potential_methylation_signal.
saveRDS(l, file=opt$output_rdata, compress=TRUE);

