#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option(c("--input"), action="store", dest="input", help="Unique GRange file"),
    make_option(c("--output"), action="store", dest="output", help="Output file")
)

parser <- OptionParser(usage="%prog [options] file", option_list=option_list);
args <- parse_args(parser, positional_arguments=TRUE);
opt <- args$options;

# Here df should be a data frame with 2 columns,
# Category and HelingerDivergence
df <- read.table(opt$input, check.names=FALSE, header=T, na.strings=c("", "NA"), stringsAsFactors=FALSE, sep="\t", quote="");

############
# Debugging.
cat("\nInitial data frame after GRange conversion:\n\n");
str(df);
cat("\nsummary(df):\n");
summary(df);
cat("head(df):\n");
head(df);
cat("\n\n");
############

# Start PDF device driver.
dev.new(width=20, height=20);

pdf(file=opt$output, width=20, height=20);
ggplot(df, aes(x=HellingerDivergence, fill=factor(Category))) +
geom_histogram(alpha=0.5, bins=50, position="identity", na.rm=TRUE, size=0.7) +
ylab("Counts") +
xlim(0,30) +
theme(axis.title.x=element_text(face="bold", size=20),
      axis.text.x=element_text(face="bold", size=20, color="black", hjust=0.5, vjust=0.75),
      axis.text.y = element_text(face="bold", size=20, color="black"),
      axis.title.y = element_text(face="bold", size=20,color="black"),
      legend.text = element_text(size=20, face="bold"),
      legend.title = element_text(size=20, face="bold"),
      plot.title = element_text(size=30, face="bold")) +
ylab("Counts") +
ggtitle("Histogram for Potential methylation signal - GGamma3P");
dev.off();

