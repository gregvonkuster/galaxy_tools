#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("adegenet"))
suppressPackageStartupMessages(library("ape"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("knitr"))
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("poppr"))
suppressPackageStartupMessages(library("RColorBrewer"))
suppressPackageStartupMessages(library("vcfR"))
suppressPackageStartupMessages(library("vegan"))

option_list <- list(
    make_option(c("--input_vcf"), action="store", dest="input_vcf", help="VCF input file"),
    make_option(c("--input_pop_info"), action="store", dest="input_pop_info", help="Population information input file"),
    make_option(c("--output_missing_data"), action="store", dest="output_missing_data", help="Missing data outputfile"),
    make_option(c("--output_mlg_id"), action="store", dest="output_mlg_id", help="Mlg Id data outputfile")
)

parser <- OptionParser(usage="%prog [options] file", option_list=option_list);
args <- parse_args(parser, positional_arguments=TRUE);
opt <- args$options;

get_file_path = function(file_name) {
    file_path = paste("output_plots_dir", file_name, sep="/");
    return(file_path);
}

# Read in VCF input file.
vcf <- read.vcfR(opt$input_vcf);

#Missing GT in samples submitted
gt <- extract.gt(vcf, element="GT", as.numeric=FALSE);
missing_gt <- apply(gt, MARGIN=2, function(x){ sum(is.na(x))});
missing_gt <- (missing_gt / nrow(vcf)) * 100;
missing_gt_data_frame <- data.frame(missing_gt);
write.table(missing_gt_data_frame, file=opt$output_missing_data, quote=FALSE);

# Convert VCF file into formats compatiable with the Poppr package.
genind <- vcfR2genind(vcf);
# Add population information to the genind object.
poptab <- read.table(opt$input_pop_info, check.names=FALSE, header=T, na.strings=c("", "NA"));
genind@pop <- as.factor(poptab$region);
# Convert genind to genclone object
gclo <- as.genclone(genind);
# Calculate the bitwise distance between individuals,
# the following is similar to Provesti's distance.
xdis <- bitwise.dist(gclo, missing_match=FALSE);

# Multilocus genotypes (threshold of 1.6%).
mlg.filter(gclo, distance=xdis) <- 0.016;
# Start PDF device driver.
dev.new(width=10, height=7);
file_path = get_file_path("mlg_table.pdf");
pdf(file=file_path, width=10, height=7);
m <- mlg.table(gclo, background=TRUE, color=TRUE);
dev.off();

# Create table of MLGs.
id <- mlg.id(gclo);
df <- data.frame(matrix((id), byrow=T));
write.table(df, file=opt$output_mlg_id);

# Rarifaction curve.
H.obj <- mlg.table(gclo, plot=TRUE);
# Start PDF device driver.
dev.new(width=10, height=7);
file_path = get_file_path("geno_rarifaction_curve.pdf");
pdf(file=file_path, width=10, height=7);
rarecurve(H.obj, ylab="Number of expected MLGs", sample=min(rowSums(H.obj)), border=NA, fill=NA, font=2, cex=1, col="blue");
dev.off()

# Create a phylogeny of samples based on distance matrices.
cols <- palette(brewer.pal(n=12, name='Set3'));
set.seed(999);
# Start PDF device driver.
dev.new(width=10, height=7);
file_path = get_file_path("nj_phylogeny.pdf");
pdf(file=file_path, width=10, height=7);
# Organize branches by clade.
tree <- gclo %>% aboot(dist=provesti.dist, sample=10, tree="nj", cutoff=50, quiet=TRUE) %>% ladderize();
plot.phylo(tree, tip.color=cols[obj2$pop],label.offset=0.0125, cex=0.7, font=2, lwd=4);
# Add a scale bar showing 5% difference..
add.scale.bar(length=0.05, cex=0.65);
nodelabels(tree$node.label, cex=.5, adj=c(1.5, -0.1), frame="n", font=3, xpd=TRUE);
dev.off()

