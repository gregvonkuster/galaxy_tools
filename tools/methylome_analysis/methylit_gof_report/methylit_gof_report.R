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
    make_option(c("--output_spec"), action="store", dest="output_spec", help="Output specifications"),
    make_option(c("--script_dir"), action="store", dest="script_dir", help="R script source directory")
)

parser <- OptionParser(usage="%prog [options] file", option_list=option_list);
args <- parse_args(parser, positional_arguments=TRUE);
opt <- args$options;

# Import the shared utility functions.
utils_path <- paste(opt$script_dir, "utils.R", sep="/");
source(utils_path);

absolute <- string_to_boolean(opt$absolute, default=FALSE);
confl_model <- string_to_boolean(opt$confl_model, default=FALSE);
model <- string_to_character_vector(opt$model);
npoints <- zero_to_null(opt$npoints);
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
if (opt$output_spec == 'all') {
    cat("## Statistics ## \n");
    gof_report$stats
    cat("\n\n");
}
cat("## Best Model ##\n");
gof_report$bestModel
cat("\n\n## NLMS ##\n");
gof_report$nlms
cat("\n\n");
sink();

