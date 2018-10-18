#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("vcfR"))
suppressPackageStartupMessages(library("poppr"))
suppressPackageStartupMessages(library("adegenet"))
suppressPackageStartupMessages(library("ape"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("knitr"))

option_list <- list(
    make_option(c("--input_vcf"), action="store", dest="input_vcf", help="VCF input file")
    make_option(c("--input_pop_info"), action="store", dest="input_pop_info", help="Population information input file")
)

parser <- OptionParser(usage="%prog [options] file", option_list=option_list);
args <- parse_args(parser, positional_arguments=TRUE);
opt <- args$options;

#extract Provesti's distance from the distance matrix
provesti_distance <- function(distance, selection) {
  eval(parse(text=paste("as.matrix(distance)[", selection, "]")));
}

#Read in VCF input file.
vcf <- read.vcfR(opts$input_vcf);

# Convert VCF file into formats compatiable with the Poppr package.
gind <- vcfR2genind(vcf);
# Add population information to the genind object.
poptab <- read.table(opt$input_pop_info, check.names=FALSE, header=T, na.strings = c("", "NA"));
gind@pop <- as.factor(poptab$region);
# Convert genind to genclone object
gclo <- as.genclone(gind);
# Calculate the bitwise distance between individuals,
# the following is similar to Provesti's distance.
xdis <- bitwise.dist(gclo);
# All alleles must match to make a unique multilocus
# genotype (“original” naive approach). This is the
# default behavior of poppr.
mll(gclo) <- "original";
# The method above does not take the genetic distance
# into account, but we can use this matrix to collapse
# MLGs that are under a specified distance threshold.
# To determine the distance threshold, we will generate
# a neighbor-joining tree for all samples.

# Create a phylogeny of samples based on distance matrices
# colors.
cols <- c("skyblue2","#C38D9E", '#E8A87C',"darkcyan","#e27d60");
set.seed(999);

theTree <- gclo %>%
  aboot(dist=provesti.dist, sample=50, tree="nj", cutoff=50, quiet=TRUE) %>%
   # Organize branches by clade.
  ladderize();
plot.phylo(theTree, tip.color=cols[gclo$pop], label.offset=0.0125, cex=0.7, font=2, lwd=4);
 # Add a scale bar showing 5% difference.
add.scale.bar(length=0.05, cex=0.65);
nodelabels(theTree$node.label, cex=.5, adj=c(1.5, -0.1), frame="n", font=3, xpd=TRUE);

# Use of mlg.filter() will create a dissimiliarity distance
# matrix from the data and then filter based off of that
# matrix. Here we will use the bitwise distance matrix
# calculated above.

# Multilocus genotypes (threshold of 1%).
mlg.filter(gclo, distance= xdis) <- 0.01;
m <- mlg.table(gclo, background=TRUE, color=TRUE);

# Different clustering methods for tie breakers used by
# mlg.filter, default is farthest neighbor.
gclo_filtered <- filter_stats(gclo, distance=xdis, plot=TRUE);

# Create table of MLGs.
id <- mlg.id(gclo);
df <- data.frame(matrix((id), byrow=T));

# We can use the genotype_curve() function to create a
# genotype accumulation curve to determine the minimum
# number of loci to identify unique MLGs.
gac <- genotype_curve(gind, sample=5, quiet=TRUE);

p <- last_plot();
p + geom_smooth() + xlim(0, 100) + theme_bw();

# From the collapsed MLGs, we can calculate genotypic
# richness, diversity and eveness.
kable(poppr(gclo));
kable(diversity_ci(gclo, n=100L, raw=FALSE ));

# Now we can correct the original data for clones using
# clonecorrect. This step will reduce the dataset to
# only have 1 representative genotype per multilocus
# lineages (MLL).
gclo_cor <- clonecorrect(gclo, strata=NA);

# Lastly, we can use a discriminant analysis of principal
# components to cluster genetically related individuals.
# This multivariate statistical approach partions the
# sample into a between-group and within- group component,
# in an effort to maximize discrimination between groups.
# Data is first transformed using a principal components
# analysis (PCA) and subsequently clusters are identified
# using discriminant analysis (DA).More information can be
# found here.
dapc.coral <- dapc(gclo_cor, var.contrib=TRUE, scale=FALSE, n.pca=62, n.da=nPop(gclo_cor)-1);
scatter(dapc.coral, cell=0, pch=18:23, cstar=0, lwd=2, lty=2, legend=TRUE, cleg=0.75, clabel=TRUE, col=cols);

