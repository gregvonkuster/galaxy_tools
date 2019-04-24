#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("adegenet"))
suppressPackageStartupMessages(library("ape"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("dbplyr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("knitr"))
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("poppr"))
suppressPackageStartupMessages(library("RColorBrewer"))
suppressPackageStartupMessages(library("rworldmap"))
suppressPackageStartupMessages(library("RPostgres"))
suppressPackageStartupMessages(library("SNPRelate"))
suppressPackageStartupMessages(library("tidyr"))
suppressPackageStartupMessages(library("vcfR"))
suppressPackageStartupMessages(library("vegan"))
suppressPackageStartupMessages(library("yarrr"))
theme_set(theme_bw())

option_list <- list(
    make_option(c("--database_connection_string"), action="store", dest="database_connection_string", help="Corals (stag) database connection string"),
    make_option(c("--input_affy_metadata"), action="store", dest="input_affy_metadata", help="Affymetrix 96 well plate input file"),
    make_option(c("--input_pop_info"), action="store", dest="input_pop_info", help="Population information input file"),
    make_option(c("--input_vcf"), action="store", dest="input_vcf", help="VCF input file"),
    make_option(c("--output_stag_db_report"), action="store", dest="output_stag_db_report", help="stag db report output file")
)

parser <- OptionParser(usage="%prog [options] file", option_list=option_list);
args <- parse_args(parser, positional_arguments=TRUE);
opt <- args$options;

get_file_path = function(file_name) {
    file_path = paste("output_plots_dir", file_name, sep="/");
    return(file_path);
}

get_database_connection <- function(db_conn_string) {
    # Instantiate database connection.
    # The connection string has this format:
    # postgresql://user:password@host/dbname
    conn_items <- strsplit(db_conn_string, "://")[[1]];
    string_needed <- conn_items[2];
    items_needed <- strsplit(string_needed, "@")[[1]];
    user_pass_string <- items_needed[1];
    host_dbname_string <- items_needed[2];
    user_pass_items <- strsplit(user_pass_string, ":")[[1]];
    host_dbname_items <- strsplit(host_dbname_string, "/")[[1]];
    user <- user_pass_items[1];
    pass <- user_pass_items[2];
    host <- host_dbname_items[1];
    dbname <- host_dbname_items[2];
    # FIXME: is there a way to not hard-code the port?
    conn <- DBI::dbConnect(RPostgres::Postgres(), host=host, port="5432", dbname=dbname, user=user, password=pass);
    return (conn);
}

time_elapsed <- function(start_time) {
    cat("Elapsed time: ", Sys.time() - start_time, "\n\n");
}

time_start <- function(msg) {
    start_time <- Sys.time();
    cat("\n", msg, "...\n");
    return(start_time);
}

# Read in VCF input file.
start_time <- time_start("Reading VCF input");
vcf <- read.vcfR(opt$input_vcf);
time_elapsed(start_time);

# Convert VCF file into a genind for the Poppr package.
start_time <- time_start("Converting VCF to a genind object");
gind <- vcfR2genind(vcf);
time_elapsed(start_time);

# Add population information to the genind object.
poptab <- read.table(opt$input_pop_info, check.names=FALSE, header=F, na.strings=c("", "NA"), stringsAsFactors=FALSE, sep="\t");
colnames(poptab) <- c("row_id", "affy_id", "user_specimen_id", "region");
gind@pop <- as.factor(poptab$region);
strata(gind)<-data.frame(pop(gind));

# Convert genind object to a genclone object.
start_time <- time_start("Converting genind to genclone object");
obj2 <- as.genclone(gind);
time_elapsed(start_time);

# Calculate the bitwise distance between individuals.
start_time <- time_start("Calculating the bitwise distance between individuals");
xdis <- bitwise.dist(obj2);
time_elapsed(start_time);

# Multilocus genotypes (threshold of 3.2%).
# threshold doubled because of how the data
# is formatted in genind compared to genlight
mlg.filter(obj2, distance=xdis) <- 0.032;
m <- mlg.table(obj2, background=TRUE, color=TRUE);

# Create table of MLGs.
id <- mlg.id(obj2);

# Read user's Affymetrix 96 well plate tabular file.
pinfo <- read.table(opt$input_affy_metadata, header=FALSE, stringsAsFactors=FALSE, sep="\t", na.strings = c("", "NA"));
colnames(pinfo) <- c("user_specimen_id", "field_call", "bcoral_genet_id", "bsym_genet_id", "reef",
                     "region", "latitude", "longitude", "geographic_origin", "sample_location",
                     "latitude_outplant", "longitude_outplant", "depth", "disease_resist",
                     "bleach_resist", "mortality","tle", "spawning", "collector_last_name",
                     "collector_first_name", "organization", "collection_date", "email", "seq_facility",
                     "array_version", "public", "public_after_date", "sperm_motility", "healing_time",
                     "dna_extraction_method", "dna_concentration", "registry_id");
pinfo$user_specimen_id <- as.character(pinfo$user_specimen_id);
pinfo2 <- as.character(pinfo$user_specimen_id);
pi <- data.table(pinfo2, pinfo$field_call);
setnames(pi, c("pinfo2"), c("user_specimen_id"));
setnames(pi, c("V2"), c("field_call"));

# Connect to database.
conn <- get_database_connection(opt$database_connection_string);

# Import the sample table.
sample_table <- tbl(conn, "sample");

# Import the genotype table.
genotype_table <- tbl(conn, "genotype");

# Select columns from the sample table and the
# genotype table joined by genotype_id.
sample_table_columns <- sample_table %>% select(user_specimen_id, affy_id, genotype_id);
smlg <- sample_table_columns %>%
    left_join(genotype_table %>%
              select("id", "coral_mlg_clonal_id", "symbio_mlg_clonal_id"),
              by=c("genotype_id"="id"));

# Convert to dataframe.
sm <- data.frame(smlg);
# Name the columns.
colnames(sm) <- c("user_specimen_id", "affy_id", "genotype_id", "coral_mlg_clonal_id", "symbio_mlg_clonal_id");

# Missing GT in samples submitted.
gt <- extract.gt(vcf, element="GT", as.numeric=FALSE);
myMiss <- apply(gt, MARGIN=2, function(x){ sum(is.na(x))});
myMiss <- (myMiss / nrow(vcf)) * 100;
miss <- data.frame(myMiss);

# Convert missing data into data table.
mi <-setDT(miss, keep.rownames=TRUE)[];
setnames(mi, c("rn"), c("affy_id"));
setnames(mi, c("myMiss"), c("percent_missing_data_coral"));
# Round data to two digits.
mi$percent_missing_data_coral <- round(mi$percent_missing_data_coral, digits=2);

# Heterozygous alleles.
start_time <- time_start("Discovering heterozygous alleles");
hets <- apply(gt, MARGIN=2, function(x) {sum(lengths(regmatches(x, gregexpr("0/1", x))))});
hets <- (hets / nrow(vcf)) * 100;
ht <- data.frame(hets);
time_elapsed(start_time);

# Convert heterozygosity data into data table.
ht <- setDT(ht, keep.rownames=TRUE)[];
setnames(ht, c("rn"), c("affy_id"));
setnames(ht, c("hets"), c("percent_heterozygous_coral"));
# Round data to two digits.
ht$percent_heterozygous_coral <- round(ht$percent_heterozygous_coral, digits=2);

# Reference alleles.
start_time <- time_start("Discovering reference alleles");
refA <- apply(gt, MARGIN=2, function(x) {sum(lengths(regmatches(x, gregexpr("0/0", x))))});
refA <- (refA / nrow(vcf)) * 100;
rA <- data.frame(refA);
time_elapsed(start_time);

# Convert refA data into data.table.
rA  <- setDT(rA, keep.rownames=TRUE)[];
setnames(rA, c("rn"), c("affy_id"));
setnames(rA, c("refA"), c("percent_reference_coral"));
# Round data to two digits.
rA$percent_reference <- round(rA$percent_reference, digits=2);

# Alternative alleles
start_time <- time_start("Discovering alternative alleles");
altB <- apply(gt, MARGIN=2, function(x) {sum(lengths(regmatches(x, gregexpr("1/1", x))))});
altB <- (altB / nrow(vcf)) * 100;
aB <- data.frame(altB);
time_elapsed(start_time);

# Convert altB data into data table.
aB <- setDT(aB, keep.rownames=TRUE)[];
setnames(aB, c("rn"), c("affy_id"));
setnames(aB, c("altB"), c("percent_alternative_coral"));
# Round data to two digits.
aB$percent_alternative <- round(aB$percent_alternative, digits=2);

# Convert mlg id to data.table format.
dt <- data.table(id, keep.rownames=TRUE);
setnames(dt, c("id"), c("affy_id"));

# Transform.
df3 <- dt %>%
    group_by(row_number()) %>%
    dplyr::rename(group="row_number()") %>%
    unnest (affy_id) %>%
    # Join with mlg table.
    left_join(sm %>%
              select("affy_id","coral_mlg_clonal_id"),
              by="affy_id");

# If found in database, group members on previous mlg id.
uniques <- unique(df3[c("group", "coral_mlg_clonal_id")]);
uniques <- uniques[!is.na(uniques$coral_mlg_clonal_id),];
na.mlg <- which(is.na(df3$coral_mlg_clonal_id));
na.group <- df3$group[na.mlg];
df3$coral_mlg_clonal_id[na.mlg] <- uniques$coral_mlg_clonal_id[match(na.group, uniques$group)];

# Determine if the sample mlg matched previous genotyped sample.
df4 <- df3 %>%
    group_by(group) %>%
    mutate(DB_match = ifelse(is.na(coral_mlg_clonal_id),"no_match", "match"));

# Create new mlg id for samples that did not match those in the database.
none <- unique(df4[c("group", "coral_mlg_clonal_id")]);
none <- none[is.na(none$coral_mlg_clonal_id),];
na.mlg2 <- which(is.na(df4$coral_mlg_clonal_id));
n.g <- df4$group[na.mlg2];
ct <- length(unique(n.g));

# List of new group ids, the sequence starts at the number of
# ids present in df4$coral_mlg_clonal_ids plus 1.
# TODO: Not sure if # the df4 file contains all ids.  If it
# doesn't then look below to change the seq() function.
n.g_ids <- sprintf("HG%04d", seq((sum(!is.na(unique(df4["coral_mlg_clonal_id"]))) + 1), by=1, length=ct));
# Pair group with new ids.
rat <- cbind(unique(n.g), n.g_ids);
# Assign the new id iteratively for all that have NA.
for (i in 1:length(na.mlg2)) {
    df4$coral_mlg_clonal_id[na.mlg2[i]] <- n.g_ids[match(df4$group[na.mlg2[i]], unique(n.g))];
}

# Subset poptab for all samples.
subpop <- poptab[c(2, 3)];

# Merge data frames for final table.
start_time <- time_start("Merging data frames");
report_user <- pi %>%
    left_join(subpop %>%
        select("affy_id", "user_specimen_id"),
        by="user_specimen_id") %>%
    left_join(df4 %>%
        select("affy_id", "coral_mlg_clonal_id", "DB_match"),
        by="affy_id") %>%
    left_join(mi %>%
        select("affy_id", "percent_missing_data_coral"),
        by="affy_id") %>%
    left_join(ht %>%
        select("affy_id", "percent_heterozygous_coral"),
        by="affy_id") %>%
    left_join(rA %>%
        select("affy_id", "percent_reference_coral"),
        by="affy_id") %>%
    left_join(aB %>%
        select("affy_id", "percent_alternative_coral"),
        by="affy_id") %>%
    mutate(DB_match = ifelse(is.na(DB_match), "failed", DB_match))%>%
    mutate(coral_mlg_clonal_id = ifelse(is.na(coral_mlg_clonal_id), "failed", coral_mlg_clonal_id)) %>%
    mutate(genetic_coral_species_call = ifelse(percent_alternative_coral >= 40 & percent_alternative_coral <= 44.5, "A.palmata","other")) %>%
    mutate(genetic_coral_species_call = ifelse(percent_alternative_coral >= 45.5 & percent_alternative_coral <= 50, "A.cervicornis", genetic_coral_species_call)) %>%
    mutate(genetic_coral_species_call = ifelse(percent_heterozygous_coral > 40, "A.prolifera", genetic_coral_species_call)) %>%
    ungroup() %>%
    select(-group);
time_elapsed(start_time);

start_time <- time_start("Writing csv output");
write.csv(report_user, file=opt$output_stag_db_report, quote=FALSE);
time_elapsed(start_time);

# Database sample table.
sample_db <- pinfo %>%
  left_join(report_user %>%
      select("user_specimen_id","affy_id", "percent_missing_data_coral", "percent_heterozygous_coral", "percent_reference_coral", "percent_alternative_coral"),
      by='user_specimen_id');

# Representative clone for genotype.table.
start_time <- time_start("Creating representative clone for genotype table");
cc <- clonecorrect(obj2, strata = ~pop.gind.);
id_rep <- mlg.id(cc);
dt_cc <- data.table(id_rep, keep.rownames=TRUE);
setnames(dt_cc, c("id_rep"), c("affy_id"));
time_elapsed(start_time);

# Transform mlg data.table.
start_time <- time_start("Selecting from various database tables");
df_cc <- dt_cc %>%
    group_by(row_number()) %>%
    rename(group='row_number()') %>%
    unnest(affy_id) %>%
    left_join(report_user %>%
              select("coral_mlg_clonal_id", "user_specimen_id", "affy_id"),
              by='affy_id') %>%
    mutate(coral_mlg_rep_sample_id=ifelse(is.na(coral_mlg_clonal_id), "", affy_id)) %>%
    ungroup() %>%
    select(-group);

# Database genotype table.
geno_db <- df4 %>%
    left_join(df_cc %>%
    select("affy_id", "coral_mlg_rep_sample_id", "user_specimen_id"),
    by='affy_id') %>%
    ungroup() %>%
    select(-group);

# Database taxonomy table.
tax_db <- report_user %>%
    select(genetic_coral_species_call, affy_id) %>%
    mutate(genus_name = ifelse(genetic_coral_species_call == genetic_coral_species_call[grep("^A.*", genetic_coral_species_call)], "Acropora", "other")) %>%
    mutate(species_name = ifelse(genetic_coral_species_call == "A.palmata", "palmata", "other")) %>%
    mutate(species_name = ifelse(genetic_coral_species_call == "A.cervicornis", "cervicornis", species_name)) %>%
    mutate(species_name = ifelse(genetic_coral_species_call == "A.prolifera", "prolifera", species_name));
time_elapsed(start_time);

# Table of alleles for the new samples subset to new plate data.
# Create vector indicating number of individuals desired from
# affy_id column from report_user data table.
i <- ifelse(is.na(report_user[1]), "", report_user[[1]]);
i <- i[!apply(i== "", 1, all),];
sub96 <- obj2[i, mlg.reset=FALSE, drop=FALSE];

# Convert to data frame.
at_96 <- genind2df(sub96, sep="");
at_96 <- at_96 %>%
select(-pop);

# Allele string for Allele.table in database.
uat_96 <- unite(at_96, alleles, 1:19696, sep=" ", remove=TRUE);
uat_96 <- setDT(uat_96, keep.rownames=TRUE)[];
setnames(uat_96, c("rn"), c("user_specimen_id"));

# Create a phylogeny tree of samples based on distance matrices.
cols <- piratepal("basel");
set.seed(999);

# Start PDF device driver.
start_time <- time_start("Creating nj_phylogeny_tree.pdf");
dev.new(width=10, height=7);
file_path = get_file_path("nj_phylogeny_tree.pdf");
pdf(file=file_path, width=10, height=7);
# Organize branches by clade.
theTree <- sub96 %>%
    aboot(dist=provesti.dist, sample=100, tree="nj", cutoff=50, quiet=TRUE) %>%
    ladderize();
theTree$tip.label <- report_user$user_specimen_id[match(theTree$tip.label, report_user$affy_id)];
plot.phylo(theTree, tip.color=cols[sub96$pop], label.offset=0.0125, cex=0.3, font=2, lwd=4, align.tip.label=F, no.margin=T);
# Add a scale bar showing 5% difference.
add.scale.bar(0, 0.95, length=0.05, cex=0.65, lwd=3);
nodelabels(theTree$node.label, cex=.5, adj=c(1.5, -0.1), frame="n", font=3, xpd=TRUE);
legend("topright", legend=c(levels(sub96$pop)), text.col=cols, xpd=T, cex=0.8);
#write.tree(theTree, file=opt$output_nj_tree);
dev.off()
time_elapsed(start_time);

# Identity-by-state analysis.
# Subset VCF to the user samples.
start_time <- time_start("Subsetting vcf to the user samples");
l <- length(i);
n <- ncol(vcf@gt);
s <- n - l;
svcf <- vcf[, s:n];
write.vcf(svcf, "subset.vcf.gz");
vcf.fn <- "subset.vcf.gz";
snpgdsVCF2GDS(vcf.fn, "test3.gds", method="biallelic.only");

genofile <- snpgdsOpen(filename="test3.gds", readonly=FALSE);
hd <- read.gdsn(index.gdsn(genofile, "sample.id"));
hd <- data.frame(hd);
hd <- setDT(hd, keep.rownames = FALSE)[];
setnames(hd, c("hd"), c("user_specimen_id"));

subpop2 <- poptab[c(2,4)];
poptab_sub <- hd %>%
    left_join(subpop2 %>%
        select("affy_id","region"),
        by='affy_id')%>%
        drop_na();

samp.annot <- data.frame(pop.group=c(poptab_sub$region));
add.gdsn(genofile, "sample.annot", samp.annot);

pop_code <- read.gdsn(index.gdsn(genofile, path="sample.annot/pop.group"));
pop.group <- as.factor(read.gdsn(index.gdsn(genofile, "sample.annot/pop.group")));
time_elapsed(start_time);

# Distance matrix calculation.
start_time <- time_start("Calculating distance matrix");
ibs <- snpgdsIBS(genofile, num.thread=2, autosome.only=FALSE);
time_elapsed(start_time);

# Cluster analysis on the genome-wide IBS pairwise distance matrix.
start_time <- time_start("Clustering the genome-wide IBS pairwise distance matrix");
set.seed(100);
par(cex=0.6, cex.lab=1, cex.axis=1.5,cex.main=2);
ibs.hc <- snpgdsHCluster(snpgdsIBS(genofile, autosome.only=FALSE));
time_elapsed(start_time);

# Default clustering.
# Start PDF device driver.
start_time <- time_start("Creating ibs_default.pdf");
dev.new(width=10, height=7);
file_path = get_file_path("ibs_default.pdf");
pdf(file=file_path, width=10, height=7);
rv <- snpgdsCutTree(ibs.hc, col.list=cols, pch.list=15);
snpgdsDrawTree(rv, main="Color by Cluster", leaflab="perpendicular", y.label=0.2);
legend("topleft", legend=levels(rv$samp.group), xpd=T, col=cols[1:nlevels(rv$samp.group)], pch=15, ncol=4, cex=1.2);
dev.off()
time_elapsed(start_time);

# Color cluster by region.
# Start PDF device driver.
start_time <- time_start("Creating ibs_region.pdf");
dev.new(width=10, height=7);
file_path = get_file_path("ibs_region.pdf");
pdf(file=file_path, width=10, height=7);
race <- as.factor(pop_code);
rv2 <- snpgdsCutTree(ibs.hc, samp.group=race,col.list=cols, pch.list=15);
snpgdsDrawTree(rv2, main="Color by Region", leaflab="perpendicular", y.label=0.2);
legend("topleft", legend=levels(race), xpd=T, col=cols[1:nlevels(race)], pch=15, ncol=4, cex=1.2);
dev.off()
time_elapsed(start_time);

# close GDS file.
snpgdsClose(genofile);

# Sample MLG on a map.
start_time <- time_start("Creating mlg_map.pdf");
world_map_spatial_data_frame <- getMap(resolution="low");
world_map_data_frame <- as.data.frame(world_map_spatial_data_frame);
pinfo$mlg <- report_user$coral_mlg_clonal_id;
n <- nrow(pinfo);

mxlat <- max(pinfo$latitude, na.rm=TRUE);
mnlat <- min(pinfo$latitude, na.rm=TRUE);
mxlong <- max(pinfo$longitude, na.rm=TRUE);
mnlong <- min(pinfo$longitude, na.rm=TRUE);

# TODO: figoure out a way to replace the sf calls below.
#p5 <- ggplot(data=world_map_data_frame) + geom_sf() + coord_sf(xlim=c(mnlong-3, mxlong+3), ylim=c(mnlat-3, mxlat+3), expand=FALSE);
p5 <- ggplot(data=world_map_data_frame, ylim=c(mnlat-3, mxlat+3), expand=FALSE);

colourCount = length(unique(pinfo$mlg));
getPalette = colorRampPalette(piratepal("basel"));
# Start PDF device driver.
dev.new(width=10, height=7);
file_path = get_file_path("mlg_map.pdf");
pdf(file=file_path, width=10, height=7);
p6 <- p5 +
      geom_point(data=pinfo, aes(x=longitude, y=latitude, group=mlg, color=mlg), alpha=.7, size=3) +
      scale_color_manual(values=getPalette(colourCount)) +
      theme(legend.position="bottom") +
      guides(color=guide_legend(nrow=8,byrow=F));
dev.off()
time_elapsed(start_time);

# Missing data barplot.
start_time <- time_start("Creating missing_data.pdf");
poptab$miss <- report_user$percent_missing_data_coral[match(miss$affy_id, report_user$affy_id)];
test2 <- which(!is.na(poptab$miss));
miss96 <- poptab$miss[test2];
name96 <- poptab$user_specimen_id[test2];
# Start PDF device driver.
dev.new(width=10, height=7);
file_path = get_file_path("missing_data.pdf");
pdf(file=file_path, width=10, height=7);
par(mar = c(8, 4, 4, 2));
x <- barplot(miss96, las=2, col=cols, ylim=c(0, 3), cex.axis=0.8, space=0.8, ylab="Missingness (%)", xaxt="n");
text(cex=0.6, x=x-0.25, y=-.05, name96, xpd=TRUE, srt=60, adj=1);
dev.off()
time_elapsed(start_time);

# Generate a pie chart for each sample with a genotype.
# Store the numerical and user_specimen_id values from
# report_user for the charts (user_specimen_id names
# will be used to label each chart).
start_time <- time_start("Creating percent_breakdown.pdf");
dt1 <- data.table(report_user);
dt1 <- report_user[c(-2, -3, -4)];
dt1 <- na.omit(dt1);
# Translate to N (i.e., number of samples with a
# genotype) columns and 5 rows.
tdt1 <- t(dt1);
# Make another data table and transpose it the same as dt1 to
# get numerics. These will feed into the creation of N vectors.
dt2 <- data.table(report_user);
dt2 <- report_user[c(-1, -2, -3, -4)];
# Translate to N columns and 5 rows.
tdt2 <- t(dt2);
tdt1_matrix <- as.matrix(tdt1[-1,]);
# The number of columns is the number of samples with genotypes.
nc <- ncol(tdt1_matrix);
mode(tdt1_matrix) <- "numeric";
spy <- rowMeans(tdt1_matrix);
dev.new(width=10, height=7);
file_path = get_file_path("percent_breakdown.pdf");
pdf(file=file_path, width=10, height=7);
# Average pie of all samples.
labels <- paste(c("missing data", "mixed", "reference", "alternative"), " (", round(spy, 1), "%)", sep="");
col <- c("GREY", "#006DDB", "#24FF24", "#920000");
main <- "Average breakdown of SNP assignments across all samples";
pie(spy, labels=labels, radius=0.60, col=col, main=main, cex.main=.75);
par(mfrow=c(3, 2));
col <- c("GREY", "#006DDB", "#24FF24", "#920000");
for (i in 1:nc) {
    tmp_labels <- paste(labels, " (", round(tdt1_matrix[,i], 1), "%)", sep="");
    main <- paste("Breakdown of SNP assignments for", tdt1[1, i]);
    pie(tdt1_matrix[,i], labels=tmp_labels, radius=0.90, col=col, main=main, cex.main=.85, cex=0.75);
}
dev.off()
time_elapsed(start_time);
