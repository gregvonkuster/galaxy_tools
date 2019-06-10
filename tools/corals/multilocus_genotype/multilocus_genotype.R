#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("adegenet"))
suppressPackageStartupMessages(library("ape"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("dbplyr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("knitr"))
suppressPackageStartupMessages(library("maps"))
suppressPackageStartupMessages(library("mapproj"))
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("poppr"))
suppressPackageStartupMessages(library("RColorBrewer"))
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
    make_option(c("--output_nj_phylogeny_tree"), action="store", dest="output_nj_phylogeny_tree", default=NULL, help="Flag to plot neighbor-joining phylogeny tree"),
    make_option(c("--output_stag_db_report"), action="store", dest="output_stag_db_report", help="Flag to output stag db report file")
)

parser <- OptionParser(usage="%prog [options] file", option_list=option_list);
args <- parse_args(parser, positional_arguments=TRUE);
opt <- args$options;

get_file_path = function(dir, file_name) {
    file_path = paste(dir, file_name, sep="/");
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
    cat("Elapsed time: ", proc.time() - start_time, "\n\n");
}

time_start <- function(msg) {
    start_time <- proc.time();
    cat(msg, "...\n");
    return(start_time);
}

write_data_frame <- function(dir, file_name, data_frame) {
    file_path <- get_file_path(dir, file_name);
    write.table(data_frame, file=file_path, quote=FALSE, row.names=FALSE, sep="\t");
}

# Prepare for processing.
output_data_dir = "output_data_dir";
output_plots_dir = "output_plots_dir";
# Read in VCF input file.
start_time <- time_start("Reading VCF input");
vcf <- read.vcfR(opt$input_vcf);
time_elapsed(start_time);

# Convert VCF file into a genind for the Poppr package.
start_time <- time_start("Converting VCF data to a genind object");
genind_obj <- vcfR2genind(vcf);
time_elapsed(start_time);

# Add population information to the genind object.
population_info_data_table <- read.table(opt$input_pop_info, check.names=FALSE, header=F, na.strings=c("", "NA"), stringsAsFactors=FALSE, sep="\t");
colnames(population_info_data_table) <- c("row_id", "affy_id", "user_specimen_id", "region");
#write_data_frame(output_data_dir, "population_info_data_table", population_info_data_table);
genind_obj@pop <- as.factor(population_info_data_table$region);
strata(genind_obj) <- data.frame(pop(genind_obj));

# Convert genind object to a genclone object.
start_time <- time_start("Converting the genind object to a genclone object");
genind_clone <- as.genclone(genind_obj);
time_elapsed(start_time);

# Calculate the bitwise distance between individuals.
start_time <- time_start("Calculating the bitwise distance between individuals");
bitwise_distance <- bitwise.dist(genind_clone);
time_elapsed(start_time);

# Multilocus genotypes (threshold of 3.2%).
mlg.filter(genind_clone, distance=bitwise_distance) <- 0.032;
m <- mlg.table(genind_clone, background=TRUE, color=TRUE);

# Create list of MLGs.
mlg_ids <- mlg.id(genind_clone);

# Read user's Affymetrix 96 well plate tabular file.
affy_metadata_data_frame <- read.table(opt$input_affy_metadata, header=FALSE, stringsAsFactors=FALSE, sep="\t", na.strings=c("", "NA"));
colnames(affy_metadata_data_frame) <- c("user_specimen_id", "field_call", "bcoral_genet_id", "bsym_genet_id", "reef",
                                        "region", "latitude", "longitude", "geographic_origin", "colony_location",
                                        "depth", "disease_resist", "bleach_resist", "mortality","tle",
                                        "spawning", "collector_last_name", "collector_first_name", "organization", "collection_date",
                                        "email", "seq_facility", "array_version", "public", "public_after_date",
                                        "sperm_motility", "healing_time", "dna_extraction_method", "dna_concentration", "registry_id",
                                        "result_folder_name", "plate_barcode");
#write_data_frame(output_data_dir, "affy_metadata_data_frame", affy_metadata_data_frame);
affy_metadata_data_frame$user_specimen_id <- as.character(affy_metadata_data_frame$user_specimen_id);
user_specimen_ids <- as.character(affy_metadata_data_frame$user_specimen_id);
# The specimen_id_field_call_data_table looks like this:
# user_specimen_ids V2
# test_002          prolifera
# test_003          prolifera
specimen_id_field_call_data_table <- data.table(user_specimen_ids, affy_metadata_data_frame$field_call);
# Rename the user_specimen_ids column.
setnames(specimen_id_field_call_data_table, c("user_specimen_ids"), c("user_specimen_id"));
# Rename the V2 column.
setnames(specimen_id_field_call_data_table, c("V2"), c("field_call"));

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
        select("id", "coral_mlg_clonal_id", "coral_mlg_rep_sample_id",
                "genetic_coral_species_call", "bcoral_genet_id"),
        by=c("genotype_id"="id"));
# Name the columns.
smlg_data_frame <- as.data.frame(smlg);
colnames(smlg_data_frame) <- c("user_specimen_id", "affy_id", "genotype_id", "coral_mlg_clonal_id",
                               "coral_mlg_rep_sample_id",
                               "genetic_coral_species_call", "bcoral_genet_id");
#write_data_frame(output_data_dir, "smlg_data_frame", smlg_data_frame);

# Missing GT in samples submitted.
start_time <- time_start("Discovering missing GT in samples");
gt <- extract.gt(vcf, element="GT", as.numeric=FALSE);
missing_gt <- apply(gt, MARGIN=2, function(x){ sum(is.na(x))});
missing_gt <- (missing_gt / nrow(vcf)) * 100;
missing_gt_data_frame <- data.frame(missing_gt);
# The specimen_id_field_call_data_table looks like this:
# rn                                 missing_gt
# a100000-4368120-060520-256_I07.CEL 0.06092608
# a100000-4368120-060520-256_K07.CEL 0.05077173
missing_gt_data_table <-setDT(missing_gt_data_frame, keep.rownames=TRUE)[];
# Rename the rn column.
setnames(missing_gt_data_table, c("rn"), c("affy_id"));
# Rename the missing_gt column.
setnames(missing_gt_data_table, c("missing_gt"), c("percent_missing_data_coral"));
# Round data to two digits.
missing_gt_data_table$percent_missing_data_coral <- round(missing_gt_data_table$percent_missing_data_coral, digits=2);
#write_data_frame(output_data_dir, "missing_gt_data_table.tabular", missing_gt_data_table);
time_elapsed(start_time);

# Heterozygous alleles.
start_time <- time_start("Discovering heterozygous alleles");
heterozygous_alleles <- apply(gt, MARGIN=2, function(x) {sum(lengths(regmatches(x, gregexpr("0/1", x))))});
heterozygous_alleles <- (heterozygous_alleles / nrow(vcf)) * 100;
heterozygous_alleles_data_frame <- data.frame(heterozygous_alleles);
# The heterozygous_alleles_data_table looks like this:
# rn                                 heterozygous_alleles
# a100000-4368120-060520-256_I07.CEL 73.94903
# a100000-4368120-060520-256_K07.CEL 74.40089
heterozygous_alleles_data_table <- setDT(heterozygous_alleles_data_frame, keep.rownames=TRUE)[];
# Rename the rn column.
setnames(heterozygous_alleles_data_table, c("rn"), c("affy_id"));
# Rename the heterozygous_alleles column.
setnames(heterozygous_alleles_data_table, c("heterozygous_alleles"), c("percent_heterozygous_coral"));
# Round data to two digits.
heterozygous_alleles_data_table$percent_heterozygous_coral <- round(heterozygous_alleles_data_table$percent_heterozygous_coral, digits=2);
time_elapsed(start_time);

# Reference alleles.
start_time <- time_start("Discovering reference alleles");
reference_alleles <- apply(gt, MARGIN=2, function(x) {sum(lengths(regmatches(x, gregexpr("0/0", x))))});
reference_alleles <- (reference_alleles / nrow(vcf)) * 100;
reference_alleles_data_frame <- data.frame(reference_alleles);
# The reference_alleles_data_table looks like this:
# rn                                 reference_alleles
# a100000-4368120-060520-256_I07.CEL 11.60642
# a100000-4368120-060520-256_K07.CEL 11.45918
reference_alleles_data_table  <- setDT(reference_alleles_data_frame, keep.rownames=TRUE)[];
# Rename the rn column.
setnames(reference_alleles_data_table, c("rn"), c("affy_id"));
# Rename the reference_alleles column.
setnames(reference_alleles_data_table, c("reference_alleles"), c("percent_reference_coral"));
# Round data to two digits.
reference_alleles_data_table$percent_reference <- round(reference_alleles_data_table$percent_reference, digits=2);
time_elapsed(start_time);

# Alternative alleles
start_time <- time_start("Discovering alternative alleles");
alternative_alleles <- apply(gt, MARGIN=2, function(x) {sum(lengths(regmatches(x, gregexpr("1/1", x))))});
alternative_alleles <- (alternative_alleles / nrow(vcf)) * 100;
alternative_alleles_data_frame <- data.frame(alternative_alleles);
# The alternative_alleles_data_table looks like this:
# rn                                 alternative_alleles
# a100000-4368120-060520-256_I07.CEL 14.38363
# a100000-4368120-060520-256_K07.CEL 14.08916
alternative_alleles_data_table <- setDT(alternative_alleles_data_frame, keep.rownames=TRUE)[];
# Rename the rn column.
setnames(alternative_alleles_data_table, c("rn"), c("affy_id"));
# Rename the alternative_alleles column.
setnames(alternative_alleles_data_table, c("alternative_alleles"), c("percent_alternative_coral"));
# Round data to two digits.
alternative_alleles_data_table$percent_alternative <- round(alternative_alleles_data_table$percent_alternative, digits=2);
time_elapsed(start_time);

# The mlg_ids_data_table looks like this:
# mlg_ids
# a550962-4368120-060520-500_M23.CEL
# a550962-4368120-060520-256_A19.CEL
mlg_ids_data_table <- data.table(mlg_ids, keep.rownames=TRUE);
# Rename the mlg_ids column.
setnames(mlg_ids_data_table, c("mlg_ids"), c("affy_id"));

# sample_mlg_tibble looks like this:
# A tibble: 262 x 3
# Groups:   group [?]
# group affy_id                            coral_mlg_clonal_id
# <int> <chr>                              <chr>
# 1     a550962-4368120-060520-500_M23.CEL NA
# 2     a550962-4368120-060520-256_A19.CEL HG0006
sample_mlg_tibble <- mlg_ids_data_table %>%
    group_by(row_number()) %>%
    dplyr::rename(group="row_number()") %>%
    unnest (affy_id) %>%
    # Join with mlg table.
    left_join(smlg_data_frame %>%
              select("affy_id","coral_mlg_clonal_id", "coral_mlg_rep_sample_id",
                      "genetic_coral_species_call", "bcoral_genet_id"),
              by="affy_id");

# If found in database, group members on previous mlg id.
uniques <- unique(sample_mlg_tibble[c("group", "coral_mlg_clonal_id")]);
uniques <- uniques[!is.na(uniques$coral_mlg_clonal_id),];
na.mlg <- which(is.na(sample_mlg_tibble$coral_mlg_clonal_id));
na.group <- sample_mlg_tibble$group[na.mlg];
sample_mlg_tibble$coral_mlg_clonal_id[na.mlg] <- uniques$coral_mlg_clonal_id[match(na.group, uniques$group)];
#write_data_frame(output_data_dir, "sample_mlg_tibble.tabular", sample_mlg_tibble);

# Find out if the sample mlg matched a previous genotyped sample.
# sample_mlg_match_tibble looks like this:
# A tibble: 262 x 4
# Groups:   group [230]
# group affy_id                            coral_mlg_clonal_id DB_match
# <int> <chr>                              <chr>               <chr>
# 1     a550962-4368120-060520-500_M23.CEL NA                  no_match
# 2     a550962-4368120-060520-256_A19.CEL HG0006              match
sample_mlg_match_tibble <- sample_mlg_tibble %>%
    group_by(group) %>%
    mutate(DB_match = ifelse(is.na(coral_mlg_clonal_id), "no_match", "match"));

# Create new mlg id for samples with no matches in the database.
none <- unique(sample_mlg_match_tibble[c("group", "coral_mlg_clonal_id")]);
none <- none[is.na(none$coral_mlg_clonal_id),];
na.mlg2 <- which(is.na(sample_mlg_match_tibble$coral_mlg_clonal_id));
n.g <- sample_mlg_match_tibble$group[na.mlg2];
ct <- length(unique(n.g));

# List of new group ids, the sequence starts at the number of
# ids present in sample_mlg_match_tibble$coral_mlg_clonal_ids
# plus 1.
# FIXME: Not sure if # the sample_mlg_match_tibble file
# contains all ids.  If it doesn't then look below to change
# the seq() function.
n.g_ids <- sprintf("HG%04d", seq((sum(!is.na(unique(sample_mlg_match_tibble["coral_mlg_clonal_id"]))) + 1), by=1, length=ct));

# Assign the new id iteratively for all that have NA.
for (i in 1:length(na.mlg2)) {
    sample_mlg_match_tibble$coral_mlg_clonal_id[na.mlg2[i]] <- n.g_ids[match(sample_mlg_match_tibble$group[na.mlg2[i]], unique(n.g))];
}
#write_data_frame(output_data_dir, "sample_mlg_match_tibble.tabular", sample_mlg_match_tibble);

# Subset population_info_data_table for all samples.
# affy_id_user_specimen_id_vector looks like this:
# affy_id                            user_specimen_id
# a100000-4368120-060520-256_I07.CEL 13704
# a100000-4368120-060520-256_K07.CEL 13706
affy_id_user_specimen_id_vector <- population_info_data_table[c(2, 3)];

# Merge data frames for final table.
start_time <- time_start("Merging data frames");
stag_db_report <- specimen_id_field_call_data_table %>%
    left_join(affy_id_user_specimen_id_vector %>%
        select("affy_id", "user_specimen_id"),
        by="user_specimen_id") %>%
    left_join(sample_mlg_match_tibble %>%
        select("affy_id", "coral_mlg_clonal_id", "DB_match"),
        by="affy_id") %>%
    left_join(missing_gt_data_table %>%
        select("affy_id", "percent_missing_data_coral"),
        by="affy_id") %>%
    left_join(heterozygous_alleles_data_table %>%
        select("affy_id", "percent_heterozygous_coral"),
        by="affy_id") %>%
    left_join(reference_alleles_data_table %>%
        select("affy_id", "percent_reference_coral"),
        by="affy_id") %>%
    left_join(alternative_alleles_data_table %>%
        select("affy_id", "percent_alternative_coral"),
        by="affy_id") %>%
    mutate(DB_match = ifelse(is.na(DB_match), "failed", DB_match))%>%
    mutate(coral_mlg_clonal_id = ifelse(is.na(coral_mlg_clonal_id), "failed", coral_mlg_clonal_id)) %>%
    mutate(genetic_coral_species_call = ifelse(percent_alternative_coral >= 40 & percent_alternative_coral <= 44.5, "A.palmata","other")) %>%
    mutate(genetic_coral_species_call = ifelse(percent_alternative_coral >= 45.5 & percent_alternative_coral <= 50, "A.cervicornis", genetic_coral_species_call)) %>%
    mutate(genetic_coral_species_call = ifelse(percent_heterozygous_coral > 40, "A.prolifera", genetic_coral_species_call)) %>%
    ungroup() %>%
    select(-group);
<<<<<<< HEAD

write.csv(report_user, file=opt$output_stag_db_report, quote=FALSE);

# Generate database tables from the genotyping results.
# Parse the information for the Sample.table.
sample_db <- affy_metadata_data_frame %>%
  left_join(
    stag_db_report %>%
      select("user_specimen_id","affy_id",
             "percent_missing_data_coral","percent_heterozygous_coral","percent_reference_coral",
             "percent_alternative_coral"),
    by='user_specimen_id');

# Select the representative clone for the genotype.table.
no_dup_genotypes_genind<-clonecorrect(genind_clone, strata = ~pop.genind_obj.);
id_rep<-mlg.id(cc);
id_data_table <-data.table(id_rep,keep.rownames = TRUE);
setnames(id_data_table, c("id_rep"), c("affy_id"));

# Combine with previously genotyped samples in the database.
geno_data_frame <- id_data_table %>%
  group_by(row_number()) %>%
  dplyr::rename(group='row_number()') %>%
  unnest(affy_id) %>%
  left_join(smlg_data_frame %>%
            select("affy_id","coral_mlg_rep_sample_id","coral_mlg_clonal_id", "user_specimen_id"),
            by='affy_id');

# Confirm that the representative mlg is the same between runs.
uniques2 <- unique(geno_data_frame[c("group", "coral_mlg_rep_sample_id")]);
uniques2 <- uniques2[!is.na(uniques2$coral_mlg_rep_sample_id),];
na.mlg3 <- which(is.na(df5$coral_mlg_rep_sample_id));
na.group2 <- df5$group[na.mlg3];
geno_data_frame$coral_mlg_rep_sample_id[na.mlg3] <- uniques2$coral_mlg_rep_sample_id[match(na.group2, uniques2$group)];

# Transform the representative mlg column with new genotyped samples.
representative_mlg_tibble <- geno_data_frame %>%
  mutate(coral_mlg_rep_sample_id=ifelse(is.na(coral_mlg_rep_sample_id),affy_id,coral_mlg_rep_sample_id)) %>%
  ungroup() %>%
  select(-group);

# Parse the information needed to populate the genotype.table.
geno_db <- stag_db_report %>%
  select("affy_id","coral_mlg_clonal_id", "user_specimen_id", "DB_match") %>%
  left_join(representative_mlg_tibble %>%
    select("affy_id","coral_mlg_rep_sample_id"),
    by='affy_id');

# Table of alleles for the new samples
# First subset to only the new plate data.
# Then, create vector indicating number of individuals desired made from affy_id collumn from report_user data table.
i<-ifelse(is.na(stag_db_report[2]),"",stag_db_report[[2]]);
i<-i[!apply(i == "", 1, all),];

#subset the genclone object to the user data
sample_alleles_vector<-genind_clone[i, mlg.reset = FALSE, drop = FALSE];

# Convert to the subset genclone to a data frame.
sample_alleles_df<-genind2df(sample_alleles_vector, sep="");
sample_alleles_df<- sample_alleles_df %>%
  select(-pop);

# Allele string for Allele.table in database.
concat_sample_alleles<-unite(sample_alleles_df, alleles, 1:19696, sep = " ", remove = TRUE);
concat_sample_alleles<-setDT(concat_sample_alleles, keep.rownames = TRUE)[];
setnames(concat_sample_alleles, c("rn"), c("affy_id"));
# write.csv(uat_96,file=paste("Seed_genotype_alleles.csv",sep = ""),quote=FALSE,row.names=FALSE);

# Create a phylogeny of samples based on distance matrices.
cols <- piratepal("basel");
set.seed(999);
# Start PDF device driver.
dev.new(width=10, height=7);
file_path = get_file_path("nj_phylogeny.pdf");
pdf(file=file_path, width=10, height=7);
# Organize branches by clade.
theTree <- sub96 %>%
    aboot(dist=provesti.dist, sample=100, tree="nj", cutoff=50, quiet=TRUE) %>%
    ladderize();
theTree$tip.label <- report_user$user_specimen_id[match(theTree$tip.label, report_user$affy_id)];
plot.phylo(theTree, tip.color=cols[sub96$pop], label.offset=0.0125, cex=0.3, font=2, lwd=4, align.tip.label=F, no.margin=T);
# Add a scale bar showing 5% difference..
add.scale.bar(0, 0.95, length=0.05, cex=0.65, lwd=3);
nodelabels(theTree$node.label, cex=.5, adj=c(1.5, -0.1), frame="n", font=3, xpd=TRUE);
legend("topright", legend=c(levels(sub96$pop)), text.col=cols, xpd=T, cex=0.8);
dev.off()

write.tree(theTree, file =opt$nj_tree, quote=FALSE);

# Perform identity-by-state analysis.
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("SNPRelate", version = "3.8")

# Subset VCF to the new user samples.
l<-length(i);
n<-ncol(vcf@gt);
s<-n-l;
svcf<-vcf[,s:n];
=======
time_elapsed(start_time);

start_time <- time_start("Writing csv output");
write.csv(stag_db_report, file=opt$output_stag_db_report, quote=FALSE);
time_elapsed(start_time);

# Representative clone for genotype table.
start_time <- time_start("Creating representative clone for genotype table");
no_dup_genotypes_genind <- clonecorrect(genind_clone, strata = ~pop.genind_obj.);
id_rep <- mlg.id(no_dup_genotypes_genind);
id_data_table <- data.table(id_rep, keep.rownames=TRUE);
# Rename the id_rep column.
setnames(id_data_table, c("id_rep"), c("affy_id"));
time_elapsed(start_time);

# Table of alleles for the new samples subset to new plate data.
# Create vector indicating number of individuals desired from
# affy_id column of stag_db_report data table.
i <- ifelse(is.na(stag_db_report[1]), "", stag_db_report[[1]]);
i <- i[!apply(i== "", 1, all),];
sample_alleles_vector <- genind_clone[i, mlg.reset=FALSE, drop=FALSE];

# Subset VCF to the user samples.
start_time <- time_start("Subsetting vcf to the user samples");
l <- length(i);
n <- ncol(vcf@gt);
s <- n - l;
svcf <- vcf[, s:n];
>>>>>>> c00d25c955e8dbd820cdb705ebd97e1aa9d463d6
write.vcf(svcf, "subset.vcf.gz");
vcf.fn <- "subset.vcf.gz";
snpgdsVCF2GDS(vcf.fn, "test3.gds", method="biallelic.only");
genofile <- snpgdsOpen(filename="test3.gds", readonly=FALSE);
gds_array <- read.gdsn(index.gdsn(genofile, "sample.id"));
# gds_array looks like this:
# [1] "a550962-4368120-060520-500_A03.CEL" "a550962-4368120-060520-500_A05.CEL"
# [3] "a550962-4368120-060520-500_A09.CEL" "a550962-4368120-060520-500_A11.CEL"
gds_data_frame <- data.frame(gds_array);
# gds_data_frame looks like this:
# gds_array
# a550962-4368120-060520-500_A03.CEL
# a550962-4368120-060520-500_A05.CEL
gds_data_table <- setDT(gds_data_frame, keep.rownames=FALSE)[];
# Rename the gds_array column.
setnames(gds_data_table, c("gds_array"), c("affy_id"));
# affy_id_region_list looks like this:
# affy_id                            region
# a100000-4368120-060520-256_I07.CEL USVI
# a100000-4368120-060520-256_K07.CEL USVI
affy_id_region_list <- population_info_data_table[c(2,4)];
gds_data_table_join <- gds_data_table %>%
    left_join(affy_id_region_list %>%
        select("affy_id", "region"),
        by='affy_id')%>%
        drop_na();
samp.annot <- data.frame(pop.group=c(gds_data_table_join$region));
add.gdsn(genofile, "sample.annot", samp.annot);
# population_code looks like this:
# [1] 18.361733   18.361733   18.361733   18.361733   18.361733   18.361733
# [7] 25.11844009 25.11844009 25.11844009 25.11844009 25.11844009 25.11844009
population_code <- read.gdsn(index.gdsn(genofile, path="sample.annot/pop.group"));
pop.group <- as.factor(read.gdsn(index.gdsn(genofile, "sample.annot/pop.group")));
# pop.group looks like this:
# [1] 18.361733   18.361733   18.361733   18.361733   18.361733   18.361733
# [7] 25.11844009 25.11844009 25.11844009 25.11844009 25.11844009 25.11844009
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

# cols looks like this:
#       blue1         red       green        pink      orange       blue2
# "#0C5BB0FF" "#EE0011FF" "#15983DFF" "#EC579AFF" "#FA6B09FF" "#149BEDFF"
#      green2      yellow   turquoise        poop
# "#A1C720FF" "#FEC10BFF" "#16A08CFF" "#9A703EFF"
cols <- piratepal("basel");
set.seed(999);

# Generate plots.
# Default clustering.
start_time <- time_start("Creating ibs_default.pdf");
# Start PDF device driver.
dev.new(width=40, height=20);
file_path = get_file_path(output_plots_dir, "ibs_default.pdf");
pdf(file=file_path, width=40, height=20);
rv <- snpgdsCutTree(ibs.hc, col.list=cols, pch.list=15);
snpgdsDrawTree(rv, main="Color by Cluster", leaflab="perpendicular", y.label=0.2);
legend("topleft", legend=levels(rv$samp.group), xpd=T, col=cols[1:nlevels(rv$samp.group)], pch=15, ncol=4, cex=1.2);
dev.off()
time_elapsed(start_time);

# Color cluster by region.
start_time <- time_start("Creating ibs_region.pdf");
# Start PDF device driver.
dev.new(width=40, height=20);
file_path = get_file_path(output_plots_dir, "ibs_region.pdf");
pdf(file=file_path, width=40, height=20);
race <- as.factor(population_code);
rv2 <- snpgdsCutTree(ibs.hc, samp.group=race,col.list=cols, pch.list=15);
snpgdsDrawTree(rv2, main="Color by Region", leaflab="perpendicular", y.label=0.2);
legend("topleft", legend=levels(race), xpd=T, col=cols[1:nlevels(race)], pch=15, ncol=4, cex=1.2);
dev.off()
time_elapsed(start_time);

# Missing data barplot.
start_time <- time_start("Creating missing_data.pdf");
population_info_data_table$miss <- stag_db_report$percent_missing_data_coral[match(missing_gt_data_frame$affy_id, stag_db_report$affy_id)];
test2 <- which(!is.na(population_info_data_table$miss));
miss96 <- population_info_data_table$miss[test2];
name96 <- population_info_data_table$user_specimen_id[test2];
# Start PDF device driver.
dev.new(width=20, height=10);
file_path = get_file_path(output_plots_dir, "missing_data.pdf");
pdf(file=file_path, width=20, height=10);
par(mar = c(8, 4, 4, 2));
x <- barplot(miss96, las=2, col=cols, ylim=c(0, 3), cex.axis=0.8, space=0.8, ylab="Missingness (%)", xaxt="n");
text(cex=0.6, x=x-0.25, y=-.05, name96, xpd=TRUE, srt=60, adj=1);
dev.off()
time_elapsed(start_time);

# Sample MLG on a map.
start_time <- time_start("Creating mlg_map.pdf");
# Get the lattitude and longitude boundaries for rendering
# the map.  Tese boundaries will restrict the map to focus
# (i.e., zoom) on the region of the world map from which
# the samples were taken.
max_latitude <- max(affy_metadata_data_frame$latitude, na.rm=TRUE);
min_latitude <- min(affy_metadata_data_frame$latitude, na.rm=TRUE);
latitude_range_vector <- c(min_latitude-3, max_latitude+3);
max_longitude <- max(affy_metadata_data_frame$longitude, na.rm=TRUE);
min_longitude <- min(affy_metadata_data_frame$longitude, na.rm=TRUE);
longitude_range_vector <- c(min_longitude-3, max_longitude+3);
# Get the palette colors for rendering plots.
colors <- length(unique(stag_db_report$coral_mlg_clonal_id));
# Get a color palette.
palette <- colorRampPalette(piratepal("basel"));
# Start PDF device driver.
dev.new(width=20, height=20);
file_path = get_file_path(output_plots_dir, "mlg_map.pdf");
pdf(file=file_path, width=20, height=20);
world_data = map_data("world");
# Add the coral_mlg_clonal_id column from the stag_db_report
# data fram to the affy_metadata_data_frame.
affy_metadata_data_frame$mlg <- stag_db_report$coral_mlg_clonal_id;
# Get the number of colors needed from the palette for plotting
# the sample locations on the world map.
num_colors = length(unique(affy_metadata_data_frame$mlg));
# Get a color palette.
palette = colorRampPalette(piratepal("basel"));
ggplot() +
geom_map(data=world_data, map=world_data, aes(x=long, y=lat, group=group, map_id=region), fill="white", colour="#7f7f7f") +
coord_map(xlim=longitude_range_vector, ylim=latitude_range_vector) +
geom_point(data=affy_metadata_data_frame, aes(x=longitude, y=latitude, group=mlg, colour=mlg), alpha=.7, size=3) +
scale_color_manual(values=palette(num_colors)) +
theme(legend.position="bottom") +
guides(color=guide_legend(nrow=8, byrow=F));
dev.off()
time_elapsed(start_time);

if (!is.null(opt$output_nj_phylogeny_tree)) {
    # Create a phylogeny tree of samples based on distance matrices.
    # Start PDF device driver.
    start_time <- time_start("Creating nj_phylogeny_tree.pdf");
    # Table of alleles for the new samples subset to new plate data.
    # Create vector indicating number of individuals desired from
    # affy_id column of stag_db_report data table.
    i <- ifelse(is.na(stag_db_report[1]), "", stag_db_report[[1]]);
    i <- i[!apply(i== "", 1, all),];
    sample_alleles_vector <- genind_clone[i, mlg.reset=FALSE, drop=FALSE];
    dev.new(width=40, height=80);
    file_path = get_file_path(output_plots_dir, "nj_phylogeny_tree.pdf");
    pdf(file=file_path, width=40, height=80);
    # Organize branches by clade.
    nj_phylogeny_tree <- sample_alleles_vector %>%
        aboot(dist=provesti.dist, sample=100, tree="nj", cutoff=50, quiet=TRUE) %>%
        ladderize();
    nj_phylogeny_tree$tip.label <- stag_db_report$user_specimen_id[match(nj_phylogeny_tree$tip.label, stag_db_report$affy_id)];
    plot.phylo(nj_phylogeny_tree, tip.color=cols[sample_alleles_vector$pop], label.offset=0.0125, cex=0.3, font=2, lwd=4, align.tip.label=F, no.margin=T);
    # Add a scale bar showing 5% difference.
    add.scale.bar(0, 0.95, length=0.05, cex=0.65, lwd=3);
    nodelabels(nj_phylogeny_tree$node.label, cex=.5, adj=c(1.5, -0.1), frame="n", font=3, xpd=TRUE);
    legend("topright", legend=c(levels(sample_alleles_vector$pop)), text.col=cols, xpd=T, cex=0.8);
    dev.off()
    time_elapsed(start_time);
}

# Generate a pie chart for each sample with a genotype.
# Store the numerical and user_specimen_id values from
# stag_db_report for the charts (user_specimen_id names
# will be used to label each chart).
start_time <- time_start("Creating percent_breakdown.pdf");
stag_db_report_data_table <- stag_db_report[c(-2, -3, -4)];
# Remove NA and NaN values.
stag_db_report_data_table <- na.omit(stag_db_report_data_table);
# Translate to N (i.e., number of samples with a genotype)
# columns and 5 rows.
translated_stag_db_report_data_table <- t(stag_db_report_data_table);
translated_stag_db_report_matrix <- as.matrix(translated_stag_db_report_data_table[-1,]);
# Set the storage mode of the matrix to numeric.  In some
# cases this could result in the following:
# Warning message:
# In mde(x) : NAs introduced by coercion
mode(translated_stag_db_report_matrix) <- "numeric";
# Remove NA and NaN values that may have been introduced
# by coercion.
translated_stag_db_report_matrix <- na.omit(translated_stag_db_report_matrix);
tsdbrm_row_means <- rowMeans(translated_stag_db_report_matrix, na.rm=TRUE);
dev.new(width=10, height=7);
file_path = get_file_path(output_plots_dir, "percent_breakdown.pdf");
pdf(file=file_path, width=10, height=7);
# Average pie of all samples.
labels <- paste(c("missing data", "mixed", "reference", "alternative"), " (", round(tsdbrm_row_means, 1), "%)", sep="");
col <- c("GREY", "#006DDB", "#24FF24", "#920000");
main <- "Average breakdown of SNP assignments across all samples";
pie(tsdbrm_row_means, labels=labels, radius=0.60, col=col, main=main, cex.main=.75);
par(mfrow=c(3, 2));
col <- c("GREY", "#006DDB", "#24FF24", "#920000");
# Generate a pie chart for each sample with genotypes.
for (i in 1:ncol(translated_stag_db_report_matrix)) {
    tmp_labels <- paste(labels, " (", round(translated_stag_db_report_matrix[,i], 1), "%)", sep="");
    main <- paste("Breakdown of SNP assignments for", translated_stag_db_report_data_table[1, i]);
    pie(translated_stag_db_report_matrix[,i], labels=tmp_labels, radius=0.90, col=col, main=main, cex.main=.85, cex=0.75);
}
dev.off()
time_elapsed(start_time);

# close GDS file.
snpgdsClose(genofile);

# Prepare to output data frames for input to a downstream
# tool that will use them to update the stag database.
start_time <- time_start("Building data frames for insertion into database tables");
# sample_prep_data_frame looks like this (split across comment lines):
# user_specimen_id  field_call bcoral_genet_id bsym_genet_id reef
# test_002          prolifera  NA              NA            JohnsonsReef
# region  latitude longitude geographic_origin colony_location
# Bahamas 18.36173 -64.77430 Reef              NA
# depth disease_resist bleach_resist
# 5     NA             N
# mortality tle spawning collector_last_name collector_first_name organization
# NA        NA  False    Kitchen             Sheila               Penn State
# collection_date email       seq_facility array_version public
# 2018-11-08      k89@psu.edu Affymetrix   1             True
# public_after_date sperm_motility healing_time dna_extraction_method
# NA               -9             -9            NA
# dna_concentration registry_id mlg    affy_id
# NA                NA          HG0227 a550962-4368120-060520-500_A03.CEL
# percent_missing_data_coral percent_heterozygous_coral
# 1.06                       19.10
# percent_reference_coral percent_alternative_coral
# 40.10459                39.73396
sample_prep_data_frame <- affy_metadata_data_frame %>%
    left_join(stag_db_report %>%
        select("user_specimen_id",
               "affy_id",
               "percent_missing_data_coral",
               "percent_heterozygous_coral",
               "percent_reference_coral",
               "percent_alternative_coral"),
        by='user_specimen_id');
# Get the number of rows for all data frames.
num_rows <- nrow(sample_prep_data_frame);
# Set the column names so that we can extract only those columns
# needed for the sample table.
colnames(sample_prep_data_frame) <- c("user_specimen_id", "field_call", "bcoral_genet_id", "bsym_genet_id", "reef",
                                      "region", "latitude", "longitude", "geographic_origin", "colony_location",
                                      "depth", "disease_resist", "bleach_resist", "mortality", "tle",
                                      "spawning", "collector_last_name", "collector_first_name", "organization",
                                      "collection_date", "email", "seq_facility", "array_version", "public",
                                      "public_after_date", "sperm_motility", "healing_time", "dna_extraction_method",
                                      "dna_concentration", "registry_id", "result_folder_name", "plate_barcode",
                                       "mlg", "affy_id", "percent_missing_data_coral",
                                      "percent_heterozygous_coral", "percent_reference_coral", "percent_alternative_coral");
#write_data_frame(output_data_dir, "sample_prep_data_frame.tabular", sample_prep_data_frame);

# representative_mlg_tibble looks like this:
# # A tibble: 230 x 4
# affy_id                            coral_mlg_clonal_id user_specimen_id coral_mlg_rep_sample_id
# <chr>                              <chr>               <chr>            <chr>
# a550962-4368120-060520-500_A05.CEL HG0135              test_084         a550962-4368120-060520…
representative_mlg_tibble <- id_data_table %>%
    group_by(row_number()) %>%
    rename(group='row_number()') %>%
    unnest(affy_id) %>%
    left_join(stag_db_report %>%
              select("coral_mlg_clonal_id", "user_specimen_id", "affy_id"),
              by='affy_id') %>%
    mutate(coral_mlg_rep_sample_id=ifelse(is.na(coral_mlg_clonal_id), "", affy_id)) %>%
    ungroup() %>%
    select(-group);
#write_data_frame(output_data_dir, "representative_mlg_tibble.tabular", representative_mlg_tibble);

# FIXME: We have no data for populating the allele table.
# FIXME: We have no data for populating the experiment table.
# FIXME: We have no data for populating the fragment table.

# Output the data frame for updating the colony table.
colony_table_data_frame <- data.frame(matrix(ncol=3, nrow=num_rows));
colnames(colony_table_data_frame) <- c("latitude", "longitude", "depth");
for (i in 1:num_rows) {
    colony_table_data_frame$latitude[i] <- sample_prep_data_frame$latitude[i];
    colony_table_data_frame$longitude[i] <- sample_prep_data_frame$longitude[i];
    colony_table_data_frame$depth[i] <- sample_prep_data_frame$depth[i];
}
write_data_frame(output_data_dir, "colony.tabular", colony_table_data_frame);

# Output the data frame for populating the genotype table.
# FIXME: This genotyope data frame contains 263 lines (including
# the header), while all other data frames contain 97 lines
# (including the header).  We need to figure out how these genotype
# records map to the sample table.
# genotype_table_tibble looks like this:
# # A tibble: 262 x 5
# affy_id                            coral_mlg_clonal_id DB_match coral_mlg_rep_samp… user_specimen_id
# <chr>                              <chr>               <chr>    <chr>               <chr>
# a550962-4368120-060520-500_A05.CEL HG0135              no_match a550962-4368120-06… test_084
# affy_id                            coral_mlg_clonal_id.x coral_mlg_rep_sample_id.x symbio_mlg_clonal_id
# a550962-4368120-060520-500_M23.CEL HG0135                NA                        NA
# symbio_mlg_rep_sample_id genetic_coral_species_call bcoral_genet_id bsym_genet_id DB_match
# NA                       NA                         NA              NA            no_match
# coral_mlg_clonal_id.y coral_mlg_rep_sample_id.y
# HG0135                a550962-4368120-060520-500_M23.CEL
genotype_table_tibble <- stag_db_report %>%
    select("affy_id","coral_mlg_clonal_id", "user_specimen_id", "DB_match") %>%
    left_join(representative_mlg_tibble %>%
    select("coral_mlg_clonal_id", "coral_mlg_rep_sample_id", "affy_id"),
    by='affy_id');

write_data_frame(output_data_dir, "genotype.tabular", genotype_table_tibble);

# Output the file needed for populating the person table.
person_table_data_frame <- data.frame(matrix(ncol=4, nrow=num_rows));
colnames(person_table_data_frame) <- c("last_name", "first_name", "organization", "email");
# person_table_data_frame looks like this:
# last_name first_name organization email
# Kitchen Sheila Penn State s89@psu.edu
for (i in 1:num_rows) {
    person_table_data_frame$last_name[i] <- sample_prep_data_frame$collector_last_name[i];
    person_table_data_frame$first_name[i] <- sample_prep_data_frame$collector_first_name[i];
    person_table_data_frame$organization[i] <- sample_prep_data_frame$organization[i];
    person_table_data_frame$email[i] <- sample_prep_data_frame$email[i];
}
write_data_frame(output_data_dir, "person.tabular", person_table_data_frame);

# Output the file needed for populating the phenotype table.
phenotype_table_data_frame <- data.frame(matrix(ncol=4, nrow=num_rows));
colnames(phenotype_table_data_frame) <- c("last_name", "first_name", "organization", "email");
# phenotype_table_data_frame looks like this:
# disease_resist bleach_resist mortality tle spawning sperm_motility healing_time
# NA            NA             NA        NA  False   NA              NA
for (i in 1:num_rows) {
    phenotype_table_data_frame$disease_resist[i] <- sample_prep_data_frame$disease_resist[i];
    phenotype_table_data_frame$bleach_resist[i] <- sample_prep_data_frame$bleach_resist[i];
    phenotype_table_data_frame$mortality[i] <- sample_prep_data_frame$mortality[i];
    phenotype_table_data_frame$tle[i] <- sample_prep_data_frame$tle[i];
    phenotype_table_data_frame$spawning[i] <- sample_prep_data_frame$spawning[i];
    phenotype_table_data_frame$sperm_motility[i] <- sample_prep_data_frame$sperm_motility[i];
    phenotype_table_data_frame$healing_time[i] <- sample_prep_data_frame$healing_time[i];
}
write_data_frame(output_data_dir, "phenotype.tabular", phenotype_table_data_frame);

# FIXME: We have no data for populating the probe_annotation table.

# Output the file needed for populating the reef table.
reef_table_data_frame <- data.frame(matrix(ncol=5, nrow=num_rows));
colnames(reef_table_data_frame) <- c("name", "region", "latitude", "longitude", "geographic_origin");
# reef_table_data_frame looks like this:
# name         region  latitude  longitude geographic_origin
# JohnsonsReef Bahamas 18.361733 -64.7743  Reef
for (i in 1:num_rows) {
    reef_table_data_frame$name[i] <- sample_prep_data_frame$reef[i];
    reef_table_data_frame$region[i] <- sample_prep_data_frame$region[i];
    reef_table_data_frame$latitude[i] <- sample_prep_data_frame$latitude[i];
    reef_table_data_frame$longitude[i] <- sample_prep_data_frame$longitude[i];
    reef_table_data_frame$geographic_origin[i] <- sample_prep_data_frame$geographic_origin[i];
}
write_data_frame(output_data_dir, "reef.tabular", reef_table_data_frame);

# Output the file needed for populating the sample table.
sample_table_data_frame <- data.frame(matrix(ncol=19, nrow=num_rows));
colnames(sample_table_data_frame) <- c("affy_id", "colony_location", "collection_date", "user_specimen_id",
                                       "registry_id", "depth", "dna_extraction_method", "dna_concentration",
                                       "public", "public_after_date", "percent_missing_data_coral",
                                       "percent_missing_data_sym", "percent_reference_coral",
                                       "percent_reference_sym", "percent_alternative_coral",
                                       "percent_alternative_sym", "percent_heterozygous_coral",
                                       "percent_heterozygous_sym", "field_call");
for (i in 1:num_rows) {
    sample_table_data_frame$affy_id[i] <- sample_prep_data_frame$affy_id[i];
    sample_table_data_frame$colony_location[i] <- sample_prep_data_frame$colony_location[i];
    sample_table_data_frame$collection_date[i] <- sample_prep_data_frame$collection_date[i];
    sample_table_data_frame$user_specimen_id[i] <- sample_prep_data_frame$user_specimen_id[i];
    sample_table_data_frame$registry_id[i] <- sample_prep_data_frame$registry_id[i];
    sample_table_data_frame$depth[i] <- sample_prep_data_frame$depth[i];
    sample_table_data_frame$dna_extraction_method[i] <- sample_prep_data_frame$dna_extraction_method[i];
    sample_table_data_frame$dna_concentration[i] <- sample_prep_data_frame$dna_concentration[i];
    sample_table_data_frame$public[i] <- sample_prep_data_frame$public[i];
    sample_table_data_frame$public_after_date[i] <- sample_prep_data_frame$public_after_date[i];
    sample_table_data_frame$percent_missing_data_coral[i] <- sample_prep_data_frame$percent_missing_data_coral[i];
    # FIXME: we have no data for percent_missing_data_sym
    sample_table_data_frame$percent_missing_data_sym[i] <- "0.0";
    sample_table_data_frame$percent_reference_coral[i] <- sample_prep_data_frame$percent_reference_coral[i];
    # FIXME: we have no data for percent_reference_sym
    sample_table_data_frame$percent_reference_sym[i] <- "0.0";
    sample_table_data_frame$percent_alternative_coral[i] <- sample_prep_data_frame$percent_alternative_coral[i];
    # FIXME: we have no data for percent_alternative_sym
    sample_table_data_frame$percent_alternative_sym[i] <- "0.0";
    sample_table_data_frame$percent_heterozygous_coral[i] <- sample_prep_data_frame$percent_heterozygous_coral[i];
    # FIXME: we have no data for percent_heterozygous_sym.
    sample_table_data_frame$percent_heterozygous_sym[i] <- "0.0";
    sample_table_data_frame$field_call[i] <- sample_prep_data_frame$field_call[i];
}
write_data_frame(output_data_dir, "sample.tabular", sample_table_data_frame);

# Output the file needed for populating the taxonomy table.
# taxonomy_table_data_frame looks like this:
# genetic_coral_species_call affy_id                            genus_name species_name
# A.palmata                  a550962-4368120-060520-500_A05.CEL Acropora   palmata
taxonomy_table_data_frame <- stag_db_report %>%
    select(genetic_coral_species_call, affy_id) %>%
    mutate(genus_name = ifelse(genetic_coral_species_call == genetic_coral_species_call[grep("^A.*", genetic_coral_species_call)], "Acropora", "other")) %>%
    mutate(species_name = ifelse(genetic_coral_species_call == "A.palmata", "palmata", "other")) %>%
    mutate(species_name = ifelse(genetic_coral_species_call == "A.cervicornis", "cervicornis", species_name)) %>%
    mutate(species_name = ifelse(genetic_coral_species_call == "A.prolifera", "prolifera", species_name));
colnames(taxonomy_table_data_frame) <- c("genetic_coral_species_call", "affy_id", "genus_name", "species_name");
write_data_frame(output_data_dir, "taxonomy.tabular", taxonomy_table_data_frame);
time_elapsed(start_time);
