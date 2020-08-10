#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("Biobase"))
suppressPackageStartupMessages(library("BiocManager"))
suppressPackageStartupMessages(library("GenomicRanges"))
suppressPackageStartupMessages(library("MethylIT"))
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("xtable"))

option_list <- list(
    make_option(c("--absolute"), action="store", dest="absolute", help="If TRUE, then TV is transformed into |TV|"),
    make_option(c("--column"), action="store", dest="column", type="integer", help="number denoting the index of the GRanges column where the information divergence is given"),
    make_option(c("--confl_model"), action="store", dest="confl_model", help="If TRUE, then the best model based on highest R.Cross.val is returned"),
    make_option(c("--input"), action="store", dest="input", help="Input information divergence estimator file"),
    make_option(c("--model"), action="store", dest="model", help="A comma-separated string naming the models to fit"),
    make_option(c("--npoints"), action="store", dest="npoints", type="integer", help="Number of points to be used in the fit"),
    make_option(c("--num_cores"), action="store", dest="num_cores", type="integer", help="The number of cores to use"),
    make_option(c("--output_gof"), action="store", dest="output_gof", help="Output rdata file"),
    make_option(c("--output_txt"), action="store", dest="output_txt", help="Output txt file"),
    make_option(c("--output_spec"), action="store", dest="output_spec", help="Output specifications")
)

parser <- OptionParser(usage="%prog [options] file", option_list=option_list);
args <- parse_args(parser, positional_arguments=TRUE);
opt <- args$options;

# Convert absolute to boolean.
if (opt$absolute == 'no') {
    absolute <- FALSE;
} else {
    absolute <- TRUE;
}

# Convert confl_model to boolean.
if (opt$confl_model == 'no') {
    confl_model <- FALSE;
} else {
    confl_model <- TRUE;
}

# Convert opt$model into a character vector.
model_str <- as.character(opt$model);
model_list <- strsplit(model_str, ",")[[1]];
model <- c(unlist(model_list, use.names=FALSE));

# If npoints is 0, comnvert to NULL.
if (opt$npoints == 0) {
    npoints <- NULL;
} else {
    npoints <- opt$npoints;
}

inf_div <- readRDS(opt$input);

############
# Debugging.
cat("\nopt$input: ", opt$input, "\n");
cat("\nabsolute: ", absolute, "\n");
cat("\ncolumn: ", opt$column, "\n");
cat("\nconfl_model: ", confl_model, "\n");
cat("\nmodel: ", toString(model), "\n");
cat("\nnpoints: ", npoints, "\n");
cat("\nopt$output_spec: ", opt$output_spec, "\n");
cat("\n\n");
cat("inf_div:\n");
show(inf_div);
cat("\n\n");

gof_report <- gofReport(HD=inf_div,
                        model=model,
                        column=opt$column,
                        absolute=absolute,
                        output=opt$output_spec,
                        confl_model=confl_model,
                        npoints=npoints,
                        num.cores=opt$num_cores,
                        verbose=TRUE);

############
# Debugging.
cat("\ngof_report:\n");
gof_report
cat("\n\n");
############

# Save the gof report.
saveRDS(gof_report, file=opt$output_gof, compress=TRUE);

# Output the gof report HTML file.
sink(opt$output_txt);
gof_report
cat("\n\n");
sink();

