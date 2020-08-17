#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option(c("--input"), action="store", dest="input", help="Unique GRange file"),
    make_option(c("--title"), action="store", dest="title", help="Plot title"),
    make_option(c("--output"), action="store", dest="output", help="Output file")
)

parser <- OptionParser(usage="%prog [options] file", option_list=option_list);
args <- parse_args(parser, positional_arguments=TRUE);
opt <- args$options;

fun_length <- function(x){
    return(data.frame(y=median(x)+1, label=paste0("n=", length(x))))
}

# Here df should be a data frame with 2 columns,
# Category and HelingerDivergence
df <- read.table(opt$input, check.names=FALSE, header=T, na.strings=c("", "NA"), stringsAsFactors=FALSE, sep="\t", quote="");

############
# Debugging.
cat("\nInput data frame:\n\n");
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
ggplot(df, aes(x=Category, y=HellingerDivergence , fill=factor(Category))) +
geom_boxplot(na.rm=TRUE) +
ylim(-0.1, 20) +
stat_summary(fun.data=fun_length, geom="text", na.rm=TRUE, position=position_dodge(width=0.9), vjust=1, size=6, fontface="bold") +
theme(axis.title.x=element_text(face="bold", size=20),
      axis.text.x=element_text(face="bold", size=20, color="black", hjust=0.5, vjust=0.75),
      axis.text.y=element_text(face="bold", size=20, color="black"),
      axis.title.y=element_text(face="bold", size=20,color="black"),
      plot.title = element_text(size=30, face="bold"),
      legend.position="none") +
ggtitle(opt$title)
dev.off();

