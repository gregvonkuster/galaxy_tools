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

DEFAULT_MISSING_NUMERIC_VALUE <- -9.000000;

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
    conn <- DBI::dbConnect(RPostgres::Postgres(), host=host, port="5432", dbname=dbname, user=user, password=pass);
    return (conn);
}

log_data_frame <- function(name, df) {
    cat("\n", name, ":\n");
    show(df);
    cat("\n\n");
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
    cat("\nWriting file: ", file_name, "\n");
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
cat("\ngenind_obj:\n");
genind_obj
cat("\n\n");
time_elapsed(start_time);

# Add population information to the genind object.
population_info_data_table <- read.table(opt$input_pop_info, check.names=FALSE, header=F, na.strings=c("", "NA"), stringsAsFactors=FALSE, sep="\t", quote="");
colnames(population_info_data_table) <- c("row_id", "affy_id", "user_specimen_id", "region");
cat("\npopulation_info_data_table:\n");
population_info_data_table
cat("\n\n");
#write_data_frame(output_data_dir, "population_info_data_table", population_info_data_table);
genind_obj@pop <- as.factor(population_info_data_table$region);
strata(genind_obj) <- data.frame(pop(genind_obj));

# Convert genind object to a genclone object.
start_time <- time_start("Converting the genind object to a genclone object");
genind_clone <- as.genclone(genind_obj);
cat("\ngenind_clone:\n");
genind_clone
cat("\n\n");
time_elapsed(start_time);
# Remove genind object from memory.
rm(genind_obj);

# Calculate the bitwise distance between individuals.
start_time <- time_start("Calculating the bitwise distance between individuals");
bitwise_distance <- bitwise.dist(genind_clone);
time_elapsed(start_time);

# Multilocus genotypes (threshold of 3.2%).
cat("\nFiltering multilocus genotypes with threshold of 3.2%...\n\n");
mlg.filter(genind_clone, distance=bitwise_distance) <- 0.032;

# Create list of MLGs.
cat("\nCreating list of mlg_ids...\n\n");
mlg_ids <- mlg.id(genind_clone);
cat("\nCreated list of mlg_ids...\n\n");

# Read user's Affymetrix 96 well plate tabular file.
cat("\nCreating affy_metadata_data_frame...\n\n");
affy_metadata_data_frame <- read.table(opt$input_affy_metadata, header=FALSE, stringsAsFactors=FALSE, sep="\t", na.strings=c("", "NA"), quote="");
cat("\nCreated affy_metadata_data_frame...\n\n");
colnames(affy_metadata_data_frame) <- c("user_specimen_id", "field_call", "bcoral_genet_id", "bsym_genet_id", "reef",
                                        "region", "latitude", "longitude", "geographic_origin", "colony_location",
                                        "depth", "disease_resist", "bleach_resist", "mortality","tle",
                                        "spawning", "collector_last_name", "collector_first_name", "organization", "collection_date",
                                        "email", "seq_facility", "array_version", "public", "public_after_date",
                                        "sperm_motility", "healing_time", "dna_extraction_method", "dna_concentration", "registry_id",
                                        "result_folder_name", "plate_barcode");
affy_metadata_data_frame$user_specimen_id <- as.character(affy_metadata_data_frame$user_specimen_id);
log_data_frame("affy_metadata_data_frame", affy_metadata_data_frame);
user_specimen_ids <- as.character(affy_metadata_data_frame$user_specimen_id);
cat("\nuser_specimen_ids:\n", toString(user_specimen_ids), "\n\n");
# The specimen_id_field_call_data_table looks like this:
# user_specimen_ids V2
# 1090              prolifera
# 1091              prolifera
cat("\nCreating specimen_id_field_call_data_table...\n");
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
# Import the probe_annotation table.
probe_annotation_table <- tbl(conn, "probe_annotation");
# Select columns from the sample table and the
# genotype table joined by genotype_id.
sample_table_columns <- sample_table %>% select(user_specimen_id, affy_id, bcoral_genet_id, genotype_id);
smlg <- sample_table_columns %>%
    left_join(genotype_table %>%
        select("id", "coral_mlg_clonal_id", "coral_mlg_rep_sample_id", "genetic_coral_species_call"),
        by=c("genotype_id"="id"));
# Name the columns.
smlg_data_frame <- as.data.frame(smlg);
colnames(smlg_data_frame) <- c("user_specimen_id", "affy_id", "bcoral_genet_id", "genotype_id",
		               "coral_mlg_clonal_id", "coral_mlg_rep_sample_id", "genetic_coral_species_call");
log_data_frame("smlg_data_frame", smlg_data_frame);
# Missing GT in samples submitted.
start_time <- time_start("Discovering missing GT in samples");
gt <- extract.gt(vcf, element="GT", as.numeric=FALSE);
missing_gt <- apply(gt, MARGIN=2, function(x){ sum(is.na(x))});
missing_gt <- (missing_gt / nrow(vcf)) * 100;
missing_gt_data_frame <- data.frame(missing_gt);
log_data_frame("missing_gt_data_frame", missing_gt_data_frame);
# The specimen_id_field_call_data_table looks like this:
# rn                                 missing_gt
# a100000-4368120-060520-256_I07.CEL 0.06092608
# a100000-4368120-060520-256_K07.CEL 0.05077173
cat("\nCreating missing_gt_data_table...\n");
missing_gt_data_table <- setDT(missing_gt_data_frame, keep.rownames=TRUE)[];
# Rename the rn column.
setnames(missing_gt_data_table, c("rn"), c("affy_id"));
# Rename the missing_gt column.
setnames(missing_gt_data_table, c("missing_gt"), c("percent_missing_data_coral"));
# Round data to two digits.
missing_gt_data_table$percent_missing_data_coral <- round(missing_gt_data_table$percent_missing_data_coral, digits=2);
time_elapsed(start_time);

# Subset genotypes for the fixed SNPs by probe id.
# Select columns from the probe_annotation table.
probe_annotation_table_columns <- probe_annotation_table %>% select(probe_set_id, custid, fixed_status, acerv_allele);
# Convert to data frame.
fixed_snp_data_frame <- as.data.frame(probe_annotation_table_columns);
# Name the columns.
colnames(fixed_snp_data_frame) <- c("probe_set_id", "custid", "fixed_status", "acerv_allele");
# Filter unwanted rows.
fixed_snp_data_frame <- subset(fixed_snp_data_frame, fixed_snp_data_frame$fixed_status=="keep");
log_data_frame("fixed_snp_data_frame", fixed_snp_data_frame);
gt_fixed <- gt[rownames(gt) %in% fixed_snp_data_frame$probe_set_id, ];
log_data_frame("gt_fixed", gt_fixed);

# Missing GT in fixed SNPs.
missing_gt_fixed <- apply(gt_fixed, MARGIN=2, function(x){ sum(is.na(x))});
missing_gt_fixed <- (missing_gt_fixed / nrow(gt_fixed)) * 100;
missing_gt_fixed_data_frame <- data.frame(missing_gt_fixed);
log_data_frame("missing_gt_fixed_data_frame", missing_gt_fixed_data_frame);
cat("\nCreating missing_gt_fixed_data_table...\n");
missing_gt_fixed_data_table <- setDT(missing_gt_fixed_data_frame, keep.rownames=TRUE)[];
# Rename the rn column.
setnames(missing_gt_fixed_data_table, c("rn"), c("affy_id"));
# Rename the missing_gt column.
setnames(missing_gt_fixed_data_table, c("missing_gt_fixed"), c("percent_missing_data_fixed"));
# Round data to two digits.
missing_gt_fixed_data_table$percent_missing_data_fixed <- round(missing_gt_fixed_data_table$percent_missing_data_fixed, digits=2);

# Heterozygous alleles for fixed SNPs.
start_time <- time_start("Discovering heterozygous alleles");
heterozygous_alleles <- apply(gt_fixed, MARGIN=2, function(x) {sum(lengths(regmatches(x, gregexpr("0/1", x))))});
heterozygous_alleles <- (heterozygous_alleles / nrow(gt_fixed)) * 100;
heterozygous_alleles_data_frame <- data.frame(heterozygous_alleles);
log_data_frame("heterozygous_alleles_data_frame", heterozygous_alleles_data_frame);
# The heterozygous_alleles_data_table looks like this:
# rn                                 heterozygous_alleles
# a100000-4368120-060520-256_I07.CEL 73.94903
# a100000-4368120-060520-256_K07.CEL 74.40089
cat("\nCreating heterozygous_alleles_data_table...\n");
heterozygous_alleles_data_table <- setDT(heterozygous_alleles_data_frame, keep.rownames=TRUE)[];
# Rename the rn column.
setnames(heterozygous_alleles_data_table, c("rn"), c("affy_id"));
# Rename the heterozygous_alleles column.
setnames(heterozygous_alleles_data_table, c("heterozygous_alleles"), c("percent_heterozygous_coral"));
# Round data to two digits.
heterozygous_alleles_data_table$percent_heterozygous_coral <- round(heterozygous_alleles_data_table$percent_heterozygous_coral, digits=2);
time_elapsed(start_time);

# Create list of Acerv reference and alternative probes.
rAC <- subset(fixed_snp_data_frame, fixed_snp_data_frame$acerv_allele=="reference");
aAC <- subset(fixed_snp_data_frame, fixed_snp_data_frame$acerv_allele=="alternative");

# Subset probes for the reference and alternative SNPs in Acerv.
ref_ac <- gt_fixed[rownames(gt_fixed) %in% rAC$probe,];
alt_ac <- gt_fixed[rownames(gt_fixed) %in% aAC$probe,];

# Reference alleles for each species.
reference_alleles_ac <- apply(ref_ac, MARGIN=2, function(x) {sum(lengths(regmatches(x, gregexpr("0/0", x))))});
reference_alleles_ap <- apply(alt_ac, MARGIN=2, function(x) {sum(lengths(regmatches(x, gregexpr("0/0", x))))});

# Alternative alleles for each species.
alternative_alleles_ac <- apply(alt_ac, MARGIN=2, function(x) {sum(lengths(regmatches(x, gregexpr("1/1", x))))});
alternative_alleles_ap <- apply(ref_ac, MARGIN=2, function(x) {sum(lengths(regmatches(x, gregexpr("1/1", x))))});

# Apalm alleles.
start_time <- time_start("Discovering reference alleles");
ap_sum <- rowSums(cbind(reference_alleles_ap,alternative_alleles_ap));
ap_alleles <- (ap_sum / nrow(gt_fixed)) * 100;
ap_alleles_data_frame <- data.frame(ap_alleles);
log_data_frame("ap_alleles_data_frame", ap_alleles_data_frame);
# The reference_alleles_data_table looks like this:
# rn                                 reference_alleles
# a100000-4368120-060520-256_I07.CEL 11.60642
# a100000-4368120-060520-256_K07.CEL 11.45918
cat("\nCreating ap_alleles_data_table...\n");
ap_alleles_data_table <- setDT(ap_alleles_data_frame, keep.rownames=TRUE)[];
# Rename the rn column.
setnames(ap_alleles_data_table, c("rn"), c("affy_id"));
# Rename the reference_alleles column.
setnames(ap_alleles_data_table, c("ap_alleles"), c("percent_apalm_coral"));
# Round data to two digits.
ap_alleles_data_table$percent_apalm_coral <- round(ap_alleles_data_table$percent_apalm_coral, digits=2);
time_elapsed(start_time);

# Acerv alleles.
start_time <- time_start("Discovering alternative alleles");
ac_sum <- rowSums(cbind(reference_alleles_ac,alternative_alleles_ac));
ac_alleles <- (ac_sum  / nrow(gt_fixed)) * 100;
ac_alleles_data_frame <- data.frame(ac_alleles);
log_data_frame("ac_alleles_data_frame", ac_alleles_data_frame);
# The alternative_alleles_data_table looks like this:
# rn                                 alternative_alleles
# a100000-4368120-060520-256_I07.CEL 14.38363
# a100000-4368120-060520-256_K07.CEL 14.08916
cat("\nCreating ac_alleles_data_table...\n");
ac_alleles_data_table <- setDT(ac_alleles_data_frame, keep.rownames=TRUE)[];
# Rename the rn column.
setnames(ac_alleles_data_table, c("rn"), c("affy_id"));
# Rename the alternative_alleles column.
setnames(ac_alleles_data_table, c("ac_alleles"), c("percent_acerv_coral"));
# Round data to two digits.
ac_alleles_data_table$percent_acerv_coral <- round(ac_alleles_data_table$percent_acerv_coral, digits=2);
time_elapsed(start_time);

# The mlg_ids_data_table looks like this:
# mlg_ids
# a550962-4368120-060520-500_M23.CEL
# a550962-4368120-060520-256_A19.CEL
cat("\nCreating mlg_ids_data_table...\n");
mlg_ids_data_table <- data.table(mlg_ids, keep.rownames=TRUE);
# Rename the mlg_ids column.
setnames(mlg_ids_data_table, c("mlg_ids"), c("affy_id"));

# sample_mlg_tibble looks like this:
# A tibble: 262 x 3
# Groups:   group [?]
# group affy_id          coral_mlg_clonal_id coral_mlg_rep_sample_id
# <int> <chr>            <chr>               <chr>
# 1     a550962-4368.CEL NA                  13905
sample_mlg_tibble <- mlg_ids_data_table %>%
    group_by(row_number()) %>%
    dplyr::rename(group="row_number()") %>%
    unnest (affy_id) %>%
    # Join with mlg table.
    left_join(smlg_data_frame %>%
              select("affy_id","coral_mlg_clonal_id", "coral_mlg_rep_sample_id"),
              by="affy_id");

# If found in database, group members on previous mlg id.
uniques <- unique(sample_mlg_tibble[c("group", "coral_mlg_clonal_id")]);
uniques <- uniques[!is.na(uniques$coral_mlg_clonal_id),];
na.mlg <- which(is.na(sample_mlg_tibble$coral_mlg_clonal_id));
na.group <- sample_mlg_tibble$group[na.mlg];
sample_mlg_tibble$coral_mlg_clonal_id[na.mlg] <- uniques$coral_mlg_clonal_id[match(na.group, uniques$group)];

# Find out if the sample mlg matched a previous genotyped sample.
# sample_mlg_match_tibble looks like this:
# A tibble: 262 x 4
# Groups:   group [230]
# group affy_id         coral_mlg_clonal_id db_match
# <int> <chr>           <chr>               <chr>
# 1     a550962-436.CEL NA                  no_match
sample_mlg_match_tibble <- sample_mlg_tibble %>%
    group_by(group) %>%
    mutate(db_match = ifelse(is.na(coral_mlg_clonal_id), "no_match", "match"));

# Create new mlg id for samples with no matches in the database.
none <- unique(sample_mlg_match_tibble[c("group", "coral_mlg_clonal_id")]);
none <- none[is.na(none$coral_mlg_clonal_id),];
na.mlg2 <- which(is.na(sample_mlg_match_tibble$coral_mlg_clonal_id));
n.g <- sample_mlg_match_tibble$group[na.mlg2];
ct <- length(unique(n.g));

# List of new group ids, the sequence starts at the number of
# ids present in sample_mlg_match_tibble$coral_mlg_clonal_ids
# plus 1.
n.g_ids <- sprintf("HG%04d", seq((sum(!is.na(unique(sample_mlg_match_tibble["coral_mlg_clonal_id"]))) + 1), by=1, length=ct));

# Assign the new id iteratively for all that have NA.
for (i in 1:length(na.mlg2)) {
    sample_mlg_match_tibble$coral_mlg_clonal_id[na.mlg2[i]] <- n.g_ids[match(sample_mlg_match_tibble$group[na.mlg2[i]], unique(n.g))];
}

# Subset population_info_data_table for all samples.
# affy_id_user_specimen_id_vector looks like this:
# affy_id         user_specimen_id
# a100000-432.CEL 13704
affy_id_user_specimen_id_vector <- population_info_data_table[c(2, 3)];

# Merge data frames for final table.
start_time <- time_start("Merging data frames");
stag_db_report <- specimen_id_field_call_data_table %>%
    left_join(affy_id_user_specimen_id_vector %>%
        select("affy_id", "user_specimen_id"),
        by="user_specimen_id") %>%
    mutate(db_record = ifelse(affy_id %in% smlg_data_frame$affy_id, "genotyped", "new")) %>%
    filter(db_record=="new") %>%
    left_join(sample_mlg_match_tibble %>%
        select("affy_id", "coral_mlg_clonal_id", "db_match"),
        by="affy_id") %>%
    left_join(missing_gt_data_table %>%
        select("affy_id", "percent_missing_data_coral"),
        by="affy_id") %>%
    left_join(missing_gt_fixed_data_table %>%
        select("affy_id", "percent_missing_data_fixed"),
        by="affy_id") %>%
    left_join(heterozygous_alleles_data_table %>%
        select("affy_id", "percent_heterozygous_coral"),
        by="affy_id") %>%
    left_join(ac_alleles_data_table %>%
        select("affy_id", "percent_acerv_coral"),
        by="affy_id") %>%
    left_join(ap_alleles_data_table %>%
        select("affy_id", "percent_apalm_coral"),
        by="affy_id") %>%
    mutate(db_match = ifelse(is.na(db_match), "failed", db_match))%>%
    mutate(coral_mlg_clonal_id = ifelse(is.na(coral_mlg_clonal_id), "failed", coral_mlg_clonal_id)) %>%
    mutate(genetic_coral_species_call = ifelse(percent_apalm_coral >= 85 & percent_acerv_coral <= 10, "A.palmata", "other")) %>%
    mutate(genetic_coral_species_call = ifelse(percent_acerv_coral >= 85 & percent_apalm_coral <= 10, "A.cervicornis", genetic_coral_species_call)) %>%
    mutate(genetic_coral_species_call = ifelse(percent_heterozygous_coral > 40, "A.prolifera", genetic_coral_species_call)) %>%
    ungroup() %>%
    select(-group, -db_record);
time_elapsed(start_time);

start_time <- time_start("Writing csv output");
write.csv(stag_db_report, file=opt$output_stag_db_report, quote=FALSE);
time_elapsed(start_time);

# Representative clone for genotype table.
start_time <- time_start("Creating representative clone for genotype table");
no_dup_genotypes_genind <- clonecorrect(genind_clone, strata=~pop.genind_obj.);
id_rep <- mlg.id(no_dup_genotypes_genind);
cat("\nCreating id_data_table...\n");
id_data_table <- data.table(id_rep, keep.rownames=TRUE);
# Rename the id_rep column.
setnames(id_data_table, c("id_rep"), c("affy_id"));
time_elapsed(start_time);
# Remove clonecorrect genind from memory.
rm(no_dup_genotypes_genind);

# Table of alleles for the new samples subset to new plate data.
# Create vector indicating number of individuals desired from
# affy_id column of stag_db_report data table.
i <- ifelse(is.na(stag_db_report[3]), "", stag_db_report[[3]]);
i <- i[!apply(i== "", 1, all), ];

# Subset VCF to the user samples.
start_time <- time_start("Subsetting vcf to the user samples");
affy_list <- append(stag_db_report$affy_id,"FORMAT");
svcf <- vcf[,colnames(vcf@gt) %in% affy_list];
write.vcf(svcf, "subset.vcf.gz");

# Remove original and subset VCFs written to file from R memory.
rm(svcf);
rm(vcf);

# Load in subset VCF.
vcf.fn <- "subset.vcf.gz";
snpgdsVCF2GDS(vcf.fn, "test3.gds", method="biallelic.only");
genofile <- snpgdsOpen(filename="test3.gds", readonly=FALSE);
gds_array <- read.gdsn(index.gdsn(genofile, "sample.id"));
# gds_array looks like this:
# [1] "a550962-4368120-060520-500_A03.CEL" "a550962-4368120-060520-500_A05.CEL"
gds_data_frame <- data.frame(gds_array);
log_data_frame("gds_data_frame", gds_data_frame);
# gds_data_frame looks like this:
# gds_array
# a550962-4368120-060520-500_A03.CEL
# a550962-4368120-060520-500_A05.CEL
cat("\nCreating gds_data_table...\n");
gds_data_table <- setDT(gds_data_frame, keep.rownames=FALSE)[];
# Rename the gds_array column.
setnames(gds_data_table, c("gds_array"), c("affy_id"));
# affy_id_region_list looks like this:
# affy_id                            region
# a100000-4368120-060520-256_I07.CEL USVI
# a100000-4368120-060520-256_K07.CEL USVI
affy_id_region_list <- population_info_data_table[c(2,3,4)];
gds_data_table_join <- gds_data_table %>%
    left_join(affy_id_region_list %>%
        select("affy_id", "user_specimen_id", "region"),
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

# Distance matrix calculation and sample labels change to user specimen ids.
start_time <- time_start("Calculating distance matrix");
ibs <- snpgdsIBS(genofile, num.thread=2, autosome.only=FALSE);
ibs$sample.id <-gds_data_table_join$user_specimen_id;
time_elapsed(start_time);

# Cluster analysis on the genome-wide IBS pairwise distance matrix.
start_time <- time_start("Clustering the genome-wide IBS pairwise distance matrix");
set.seed(100);
par(cex=0.6, cex.lab=1, cex.axis=1.5, cex.main=2);
ibs.hc <- snpgdsHCluster(ibs);
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
snpgdsDrawTree(rv, main="Color by Cluster", leaflab="perpendicular", yaxis.kinship=FALSE);
abline(h=0.032, lty=2);
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
rv2 <- snpgdsCutTree(ibs.hc, samp.group=race, col.list=cols, pch.list=15);
snpgdsDrawTree(rv2, main="Color by Region", leaflab="perpendicular", yaxis.kinship=FALSE);
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
text(cex=0.8, x=x-0.25, y=-.05, name96, xpd=TRUE, srt=60, adj=1);
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
coord_quickmap(xlim=longitude_range_vector, ylim=latitude_range_vector) +
geom_point(data=affy_metadata_data_frame, aes(x=longitude, y=latitude, group=mlg, colour=mlg), alpha=.7, size=3) +
scale_color_manual(values=palette(num_colors)) +
theme(legend.position="bottom") +
guides(color=guide_legend(nrow=8, byrow=F));

# Sample MLG on a map for each region.
for (i in unique(affy_metadata_data_frame$region)) {
    m <- i;
    num_colors_2 = length(unique(affy_metadata_data_frame$mlg[which(affy_metadata_data_frame$region == m)]));
    max_latitude_region <- max(affy_metadata_data_frame$latitude[which(affy_metadata_data_frame$region == m)], na.rm=TRUE);
    min_latitude_region <- min(affy_metadata_data_frame$latitude[which(affy_metadata_data_frame$region == m)], na.rm=TRUE);
    latitude_range_vector_region <- c(min_latitude_region-0.5, max_latitude_region+0.5);
    max_longitude_region <- max(affy_metadata_data_frame$longitude[which(affy_metadata_data_frame$region == m)], na.rm=TRUE);
    min_longitude_region <- min(affy_metadata_data_frame$longitude[which(affy_metadata_data_frame$region == m)], na.rm=TRUE);
    longitude_range_vector_region <- c(min_longitude_region-0.5, max_longitude_region+0.5);
    print(ggplot() +
        geom_map(data=world_data, map=world_data, aes(x=long, y=lat, group=group, map_id=region),
                 fill="grey", colour="#7f7f7f") +
        coord_quickmap(xlim=longitude_range_vector_region, ylim=latitude_range_vector_region, clip = "on") +
        geom_point(data=affy_metadata_data_frame[which(affy_metadata_data_frame$region == m),],
                   aes(x=longitude, y=latitude, group=mlg, colour=mlg), alpha=.5, size=3) +
        scale_color_manual(values=palette(num_colors_2)) +
        theme(legend.position="bottom") + labs(title=paste("MLG assignments for", m)) +
        guides(color=guide_legend(nrow=8, byrow=F)));
}
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
        aboot(dist=provesti.dist, sample=100, tree="nj", cutoff=50, quiet=TRUE, showtree = FALSE) %>%
        ladderize();
    nj_phylogeny_tree$tip.label <- stag_db_report$user_specimen_id[match(nj_phylogeny_tree$tip.label, stag_db_report$affy_id)];
    plot.phylo(nj_phylogeny_tree, tip.color=cols[sample_alleles_vector$pop], label.offset=0.0025, cex=0.6, font=2, lwd=4, align.tip.label=F, no.margin=T);
    # Add a scale bar showing 5% difference.
    add.scale.bar(0, 0.95, length=0.05, cex=0.65, lwd=2);
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
stag_db_report_data_table <- stag_db_report[c(-2, -3, -4, -5)];
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
translated_stag_db_report_matrix<-translated_stag_db_report_matrix[-c(1),];
#tsdbrm_row_means <- rowMeans(translated_stag_db_report_matrix, na.rm=TRUE);
dev.new(width=10, height=7);
file_path = get_file_path(output_plots_dir, "percent_breakdown.pdf");
pdf(file=file_path, width=10, height=7);
par(mfrow=c(3, 2));
col <- c("#A6A6A6","#FFA626","#EB0ACF", "#80FF00" );
# Generate a pie chart for each sample with genotypes.
for (i in 1:ncol(translated_stag_db_report_matrix)) {
    tmp_labels <- paste(c("no call", "heterozygous", "A. cervicornis", "A. palmata"), " (", round(translated_stag_db_report_matrix[,i], 1), "%)", sep="");
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
# dna_concentration registry_id result_folder_name       plate_barcode mlg
# NA                NA          PRO100175_PSU175_SAX_b02 P9SR10074     HG0227
# affy_id         percent_missing_data_coral percent_heterozygous_coral
# a550962-436.CEL 1.06                       19.10
# percent_acerv_coral percent_apalm_coral
# 40.10459                39.73396
sample_prep_data_frame <- affy_metadata_data_frame %>%
    left_join(stag_db_report %>%
        select("user_specimen_id", "affy_id", "percent_missing_data_coral", "percent_heterozygous_coral",
               "percent_acerv_coral", "percent_apalm_coral"),
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
                                      "mlg", "affy_id", "percent_missing_data_coral", "percent_heterozygous_coral",
                                      "percent_acerv_coral", "percent_apalm_coral");

# Output the data frame for updating the alleles table.
# Subset to only the new plate data.
i <- ifelse(is.na(stag_db_report[3]), "", stag_db_report[[3]]);
# Create a vector indicating the number of individuals desired
# from the affy_id collumn in the report_user data table.
i <- i[!apply(i=="", 1, all),];
# Subset the genclone object to the user data.
allele_vector <- genind_clone[i, mlg.reset=FALSE, drop=FALSE];
# Convert the subset genclone to a data frame.
allele_data_frame <- genind2df(allele_vector, sep="");
allele_data_frame <- allele_data_frame %>%
select(-pop);
# Allele string for Allele.table in database.
allele_table_data_frame <- unite(allele_data_frame, alleles, 1:19696, sep=" ", remove=TRUE);
allele_table_data_frame <- setDT(allele_table_data_frame, keep.rownames=TRUE)[];
setnames(allele_table_data_frame, c("rn"), c("affy_id"));
# write.csv(concat_sample_alleles,file=paste("Seed_genotype_alleles.csv",sep = ""),quote=FALSE,row.names=FALSE);
write_data_frame(output_data_dir, "allele.tabular", allele_table_data_frame);

# Output the data frame for updating the experiment table.
experiment_table_data_frame <- data.frame(matrix(ncol=4, nrow=num_rows));
colnames(experiment_table_data_frame) <- c("seq_facility", "array_version", "result_folder_name", "plate_barcode");
for (i in 1:num_rows) {
    experiment_table_data_frame$seq_facility[i] <- sample_prep_data_frame$seq_facility[i];
    experiment_table_data_frame$array_version[i] <- sample_prep_data_frame$array_version[i];
    experiment_table_data_frame$result_folder_name[i] <- sample_prep_data_frame$result_folder_name[i];
    experiment_table_data_frame$plate_barcode[i] <- sample_prep_data_frame$plate_barcode[i];
}
write_data_frame(output_data_dir, "experiment.tabular", experiment_table_data_frame);

# Output the data frame for updating the colony table.
# The geographic_origin value is used for deciding into which table
# to insert the latitude and longitude values.  If the geographic_origin
# is "reef", the values will be inserted into the reef table, and if it is
# "colony", the values will be inserted into the colony table.  We insert
# these values in both data frames so that the downstream tool that parses
# them can determine the appropriate table.
colony_table_data_frame <- data.frame(matrix(ncol=4, nrow=num_rows));
colnames(colony_table_data_frame) <- c("latitude", "longitude", "depth", "geographic_origin");
for (i in 1:num_rows) {
    colony_table_data_frame$latitude[i] <- sample_prep_data_frame$latitude[i];
    colony_table_data_frame$longitude[i] <- sample_prep_data_frame$longitude[i];
    colony_table_data_frame$depth[i] <- sample_prep_data_frame$depth[i];
    colony_table_data_frame$geographic_origin[i] <- sample_prep_data_frame$geographic_origin[i];
}
write_data_frame(output_data_dir, "colony.tabular", colony_table_data_frame);

# Output the data frame for populating the genotype table.
# Combine with previously genotyped samples.
# prep_genotype_tibble looks like this:
# A tibble: 220 x 7
# Groups:   group [?]
# group affy_id coral_mlg_clona… user_specimen_id db_match
# <int> <chr>   <chr>            <chr>            <chr>
# 1     a10000… 13905            HG0048           match
# genetic_coral_species_call coral_mlg_rep_sample_id
# <chr>                      <chr>
# A.palmata                  1104
prep_genotype_tibble <- id_data_table %>%
    group_by(row_number()) %>%
    dplyr::rename(group='row_number()') %>%
    unnest(affy_id) %>%
    left_join(sample_mlg_match_tibble %>%
        select("affy_id", "coral_mlg_rep_sample_id", "coral_mlg_clonal_id",
               "genetic_coral_species_call", "bcoral_genet_id", "db_match"),
        by='affy_id') %>%
    right_join(sample_mlg_match_tibble %>%
        select("coral_mlg_rep_sample_id", "coral_mlg_clonal_id"),
        by='coral_mlg_clonal_id') %>%
        mutate(coral_mlg_rep_sample_id=ifelse(is.na(coral_mlg_rep_sample_id.x),coral_mlg_rep_sample_id.y,coral_mlg_rep_sample_id.x)) %>%
    ungroup() %>%
    dplyr::select(-coral_mlg_rep_sample_id.x,-coral_mlg_rep_sample_id.y, -group.x,-group.y) %>%
    distinct();

# Confirm that the representative mlg is the same between runs.
uniques2 <- unique(prep_genotype_tibble[c("group", "coral_mlg_rep_sample_id")]);
uniques2 <- uniques2[!is.na(uniques2$coral_mlg_rep_sample_id),];
na.mlg3 <- which(is.na(prep_genotype_tibble$coral_mlg_rep_sample_id));
na.group2 <- prep_genotype_tibble$group[na.mlg3];
prep_genotype_tibble$coral_mlg_rep_sample_id[na.mlg3] <- uniques2$coral_mlg_rep_sample_id[match(na.group2, uniques2$group)];
# Transform the representative mlg column with new genotyped samples.
# representative_mlg_tibble looks like this:
# A tibble: 220 x 5
# affy_id   coral_mlg_rep_sa… coral_mlg_clona… user_specimen_id
# <chr>     <chr>             <chr>            <chr>
# a100000-… 13905             HG0048           13905
# genetic_coral_species_call bcoral_genet_id
# <chr>                      <chr>
# A.palmata                  C1651
representative_mlg_tibble <- prep_genotype_tibble %>%
    mutate(coral_mlg_rep_sample_id=ifelse(is.na(coral_mlg_rep_sample_id) & (db_match =="no_match"), affy_id, coral_mlg_rep_sample_id)) %>%
    ungroup() %>%
    select(-group)%>%
    distinct();
# prep_genotype_table_tibble looks like this:
# affy_id       coral_mlg_clonal_id user_specimen_id db_match
# a550962...CEL HG0120              1090             match
# genetic_coral_species_call coral_mlg_rep_sample_id
# A.palmata                  1104
prep_genotype_table_tibble <- stag_db_report %>%
    select("affy_id", "coral_mlg_clonal_id", "user_specimen_id", "db_match", "genetic_coral_species_call") %>%
    left_join(representative_mlg_tibble %>%
        select("affy_id", "coral_mlg_rep_sample_id"),
        by='affy_id');
# genotype_table_tibble looks like this:
# affy_id         coral_mlg_clonal_id user_specimen_id db_match
# a550962-436.CEL HG0120              1090             match
# genetic_coral_species_call coral_mlg_rep_sample_id bcoral_genet_id
# A.palmata                  1104                    <NA>
genotype_table_tibble <- prep_genotype_table_tibble %>%
    left_join(affy_metadata_data_frame %>%
        select("user_specimen_id", "bcoral_genet_id"),
        by='user_specimen_id') %>%
    drop_na(coral_mlg_rep_sample_id);
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
phenotype_table_data_frame <- data.frame(matrix(ncol=7, nrow=num_rows));
colnames(phenotype_table_data_frame) <- c("disease_resist", "bleach_resist", "mortality", "tle",
                                          "spawning", "sperm_motility", "healing_time");
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

# Output the file needed for populating the reef table.
reef_table_data_frame <- data.frame(matrix(ncol=5, nrow=num_rows));
colnames(reef_table_data_frame) <- c("name", "region", "latitude", "longitude", "geographic_origin");
# The geographic_origin value is used for deciding into which table
# to insert the latitude and longitude values.  If the geographic_origin
# is "reef", the values will be inserted into the reef table, and if it is
# "colony", the values will be inserted into the colony table.  We insert
# these values in both data frames so that the downstream tool that parses
# them can determine the appropriate table.
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
sample_table_data_frame <- data.frame(matrix(ncol=20, nrow=num_rows));
colnames(sample_table_data_frame) <- c("affy_id", "colony_location", "collection_date", "user_specimen_id",
                                       "registry_id", "depth", "dna_extraction_method", "dna_concentration",
                                       "public", "public_after_date", "percent_missing_data_coral",
                                       "percent_missing_data_sym", "percent_acerv_coral",
                                       "percent_reference_sym", "percent_apalm_coral",
                                       "percent_alternative_sym", "percent_heterozygous_coral",
                                       "percent_heterozygous_sym", "field_call", "bcoral_genet_id");
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
    sample_table_data_frame$percent_missing_data_sym[i] <- DEFAULT_MISSING_NUMERIC_VALUE;
    sample_table_data_frame$percent_acerv_coral[i] <- sample_prep_data_frame$percent_acerv_coral[i];
    sample_table_data_frame$percent_reference_sym[i] <- DEFAULT_MISSING_NUMERIC_VALUE;
    sample_table_data_frame$percent_apalm_coral[i] <- sample_prep_data_frame$percent_apalm_coral[i];
    sample_table_data_frame$percent_alternative_sym[i] <- DEFAULT_MISSING_NUMERIC_VALUE;
    sample_table_data_frame$percent_heterozygous_coral[i] <- sample_prep_data_frame$percent_heterozygous_coral[i];
    sample_table_data_frame$percent_heterozygous_sym[i] <- DEFAULT_MISSING_NUMERIC_VALUE;
    sample_table_data_frame$field_call[i] <- sample_prep_data_frame$field_call[i];
	sample_table_data_frame$bcoral_genet_id[i] <- sample_prep_data_frame$bcoral_genet_id[i];
}
write_data_frame(output_data_dir, "sample.tabular", sample_table_data_frame);

# Output the file needed for populating the taxonomy table.
# taxonomy_table_data_frame looks like this:
# genetic_coral_species_call affy_id                            genus_name species_name
# A.palmata                  a550962-4368120-060520-500_A05.CEL Acropora   palmata
taxonomy_table_data_frame <- stag_db_report %>%
    select(genetic_coral_species_call, affy_id) %>%
    mutate(genus_name = ifelse(grepl("^A.*", genetic_coral_species_call), "Acropora",ifelse(!is.na(genetic_coral_species_call), "other", NA))) %>%
    mutate(species_name = ifelse(genetic_coral_species_call == "A.palmata", "palmata", "other")) %>%
    mutate(species_name = ifelse(genetic_coral_species_call == "A.cervicornis", "cervicornis", species_name)) %>%
    mutate(species_name = ifelse(genetic_coral_species_call == "A.prolifera", "prolifera", species_name));
colnames(taxonomy_table_data_frame) <- c("genetic_coral_species_call", "affy_id", "genus_name", "species_name");
write_data_frame(output_data_dir, "taxonomy.tabular", taxonomy_table_data_frame);
time_elapsed(start_time);
