#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("Biobase"))
suppressPackageStartupMessages(library("BiocGenerics"))
suppressPackageStartupMessages(library("BiocParallel"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("genefilter"))
suppressPackageStartupMessages(library("IRanges"))
suppressPackageStartupMessages(library("matrixStats"))
suppressPackageStartupMessages(library("MethylIT"))
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("S4Vectors"))

option_list <- list(
    make_option(c("--count_filter"), action="store", dest="count_filter", default=FALSE, help="Filter counts according to the minimum count per region per each sample"),
    make_option(c("--count_per_bp"), action="store", dest="count_per_bp", type="integer", default=NULL, help="The count per bp for each group must be greater than or equal to this value"),
    make_option(c("--filter_log_2_fc"), action="store", dest="filter_log_2_fc", default=FALSE, help="Filter results using the minimum absolute value of log2FoldChanges observed"),
    make_option(c("--input_glm"), action="store", dest="input_glm", help="GRanges file"),
    make_option(c("--maxgrpcv1"), action="store", dest="maxgrpcv1", type="integer", default=NULL, help="Low threshold of variance"),
    make_option(c("--maxgrpcv2"), action="store", dest="maxgrpcv2", type="integer", default=NULL, help="High threshold of variance"),
    make_option(c("--min_count_per_indiv"), action="store", dest="min_count_per_indiv", type="integer", help="Threshold for average group count for each gene or region per individual"),
    make_option(c("--min_log_2_fc"), action="store", dest="min_log_2_fc", type="double", help="Minimum logarithm base 2 of fold changes"),
    make_option(c("--mv_rate"), action="store", dest="mv_rate", type="double", help="Minimum Mean/Variance rate"),
    make_option(c("--num_cores"), action="store", dest="num_cores", type="integer", help="Number of processors to use"),
    make_option(c("--output"), action="store", dest="output", help="Output file"),
    make_option(c("--p_adjust_method"), action="store", dest="p_adjust_method", help="Method used to adjust the results"),
    make_option(c("--pval_cutoff"), action="store", dest="pval_cutoff", type="double", help="Cutoff used when a p-value adjustment is performed"),
    make_option(c("--save_all"), action="store", dest="save_all", default=FALSE, help="Include all temporal results in the output"),
    make_option(c("--scaling"), action="store", dest="scaling", type="integer", help="Scaling factor to estimate the signal density"),
    make_option(c("--test"), action="store", dest="test", help="Test")
)

parser <- OptionParser(usage="%prog [options] file", option_list=option_list);
args <- parse_args(parser, positional_arguments=TRUE);
opt <- args$options;

# Here glm is a glmDataset.
glm <- readRDS(opt$input_glm);

############
# Debugging.
cat("\nInput glm:\n\n");
glm
cat("\n\n");
############

if (!is.null(opt$maxgrpcv1)) {
    if (!is.null(opt$maxgrpcv2)) {
        max_grp_cv <- c(opt$maxgrpcv1, opt$maxgrpcv2);
    } else {
        max_grp_cv <- c(opt$maxgrpcv1, opt$maxgrpcv1);
    }
} else {
    max_grp_cv <- NULL;
}

# Create a data frame or GRanges object (if the DS contain
# the GRanges information for each gene) with the test results
# and original count matrix, plus control and treatment signal
# densities and their variation.
#                      CountPerBP=opt$count_per_bp,
#                      minCountPerIndiv=opt$min_count_per_indiv,
#                      tasks=0L,
#                      num.cores=opt$num_cores,
ret_val <- countTest2(glm,
                      countFilter=opt$count_filter,
                      maxGrpCV=max_grp_cv,
                      FilterLog2FC=opt$filter_log_2_fc,
                      pAdjustMethod=opt$p_adjust_method,
                      pvalCutOff=opt$pval_cutoff,
                      MVrate=opt$mv_rate,
                      Minlog2FC=opt$min_log_2_fc,
                      test=opt$test,
                      scaling=opt$scaling,
                      saveAll=opt$save_all);

############
# Debugging.
cat("\nret_val: \n\n");
ret_val
cat("\n\n");
############

# Save the glm_dataset.
saveRDS(ret_val, file=opt$output, compress=TRUE);

