#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("optparse"))

option_list <- list(
        make_option(c("--build"), action="store", dest="build", help="Genome build"),
        make_option(c("--chrom_len_file"), action="store", dest="chrom_len_file", help="Chromosome length file"),
        make_option(c("--email"),  action="store", dest="email", help="User email address"),
        make_option(c("--galaxy_url"),  action="store", dest="galaxy_url", help="Galaxy instance base URL"),
        make_option(c("--hub_long_label"), action="store", dest="hub_long_label", help="Hub long label"),
        make_option(c("--hub_name"),  action="store", dest="hub_name", default=NULL, help="Hub name without spaces"),
        make_option(c("--hub_short_label"), action="store", dest="hub_short_label", help="Hub short label"),
        make_option(c("--input_dir_para"), action="store", dest="input_dir_para", help="Directory containing .para outputs from IDEAS"),
        make_option(c("--input_dir_state"), action="store", dest="input_dir_state", help="Directory containing .state outputs from IDEAS"),
        make_option(c("--output_trackhub"),  action="store", dest="output_trackhub", help="Output hub file"),
        make_option(c("--output_trackhub_files_path"),  action="store", dest="output_trackhub_files_path", help="Output hub extra files path"),
        make_option(c("--output_trackhub_id"),  action="store", dest="output_trackhub_id", help="Encoded output_trackhub dataset id"),
        make_option(c("--script_dir"), action="store", dest="script_dir", help="R script source directory"),
        make_option(c("--state_colors"), action="store", dest="state_colors", default=NULL, help="List of state colors"),
        make_option(c("--state_indexes"), action="store", dest="state_indexes", default=NULL, help="List of state indexes"),
        make_option(c("--state_names"), action="store", dest="state_names", default=NULL, help="List of state names")
)

parser <- OptionParser(usage="%prog [options] file", option_list=option_list)
args <- parse_args(parser, positional_arguments=TRUE)
opt <- args$options

create_primary_html = function(output_trackhub, tracks_dir, build) {
    track_files <- list.files(path=tracks_dir);
    s <- paste('<html><head></head><body>', sep="\n");
    s <- paste(s, '<h3>Contents of directory ~/myHub/', build, ' required by UCSC TrackHub</h3>\n', sep="");
    s <- paste(s, '<ul>\n', sep="")
    for (i in 1:length(track_files)) {
        s <- paste(s, '<li><a href="', 'myHub/', build, "/", track_files[i], '">', track_files[i], '</a></li>\n', sep="");
    }
    s <- paste(s, '</ul>\n</body>\n</html>', sep="");
    cat(s, file=output_trackhub);
}

create_cells = function(input_dir_state, chrom_len_file, base_track_file_name, state_indexes, state_names, state_colors) {
    # Create everything needed, including the bigbed file,
    # to render the tracks within the UCSC track hub.
    state_files = list.files(path=input_dir_state, full.names=TRUE);
    cell_type_names = get_cell_type_names(state_files[1]);
    genome_size = read.table(chrom_len_file);
    g = NULL;
    for(i in state_files) {
        tg = as.matrix(fread(i));
        t = NULL;
        for(j in 1:dim(genome_size)[1]) {
            t = c(t, which(tg[,2]==as.character(genome_size[j, 1]) & as.numeric(tg[,4]) > as.numeric(genome_size[j, 2])));
        }
        if (length(t) > 0) {
            tg = tg[-t,];
        }
        t = which(is.na(match(tg[,2], genome_size[,1]))==T);
        if (length(t)>0) {
            tg = tg[-t,];
        }
        g = rbind(g, tg);
    }
    uchr = sort(unique(as.character(g[,2])));
    g1 = NULL;
    for(i in uchr) {
        t = which(g[,2]==i);
        g1 = rbind(g1, g[t[order(as.integer(g[t, 3]))],]);
    }
    g = NULL;
    chr = as.character(g1[,2]);
    posst = as.numeric(g1[,3]);
    posed = as.numeric(g1[,4]);
    state = as.matrix(g1[,5:(dim(g1)[2]-1)]);
    L = dim(g1)[1];
    # Here n will be the same length as the
    # list of cell_type_names defined above.
    n = dim(state)[2];
    cells = as.character(colnames(g1)[5:(dim(g1)[2]-1)]);
    options(scipen=999);
    tt = which(chr[2:L]!=chr[2:L-1]);
    tt = c(tt, which(posst[2:L]!=posed[2:L-1]));
    tt = sort(unique(tt));
    for(i in 1:n) {
        state_name = cell_type_names[i];
        tstate = state[,i];
        t = c(tt, which(tstate[2:L]!=tstate[2:L-1]));
        t = sort(unique(t));
        t0 = c(0, t) + 1;
        t = c(t, L);
        np = cbind(chr[t], posst[t0], posed[t], tstate[t]);
        #cat("JJJJJJJ as.numeric(np[,4]: ", as.numeric(np[,4]), "\n");
        track_file_name_bed_unsorted <- get_track_file_name(base_track_file_name, i, "bed_unsorted");
        track_file_name_bed <- get_track_file_name(base_track_file_name, i, "bed");
        track_file_name_bigbed <- get_track_file_name(base_track_file_name, i, "bigbed");
        x = cbind(np[, 1:3], state_name, 1000, ".", np[,2:3]);
        write.table(as.matrix(x), track_file_name_bed_unsorted, quote=F, row.names=F, col.names=F);
        cmd = paste("LC_COLLATE=C sort -k1,1 -k2,2n < ", track_file_name_bed_unsorted, " > ", track_file_name_bed);
        system(cmd);
        cmd = paste("bedToBigBed ", track_file_name_bed, chrom_len_file, " ", track_file_name_bigbed);
        system(cmd);
        system(paste("rm ", track_file_name_bed_unsorted));
        system(paste("rm ", track_file_name_bed));
    }
    return(cells);
}

create_track_db = function(galaxy_url, encoded_dataset_id, input_dir_para, input_dir_state, build, chrom_len_file, tracks_dir,
        hub_name, hub_short_label, hub_long_label, specified_state_indexes, specified_state_names, specified_state_colors) {
    # Create a trackDb.txt file that includes each state.
    #cat("\n\n\n\n\nXXXXXXXXXXXXXXXX Entered create_track_db\n");
    #cat("specified_state_indexes: ", specified_state_indexes, "\n");
    #cat("specified_state_names: ", specified_state_names, "\n");
    #cat("specified_state_colors: ", specified_state_colors, "\n");
    state_names = c();
    state_colors = list();
    # Generate the default state names and colors.
    # IDEAS state indexes are zero based.
    index = 0;
    para_files = list.files(path=input_dir_para, full.names=TRUE);
    #cat("para_files: ", toString(para_files), "\n\n\n");
    for (para_file in para_files) {
        #cat("para_file: ", para_file, "\n");
        # Append the default state name.
        state_names = c(state_names, toString(index));
        #cat("\n", index,index,index,index, "state_names: ", state_names, "\n");
        # Append the default color codes vector.
        data_frame = read.table(para_file, comment="!", header=T);
        # Here specified_state_colors will look something like this:
        # 248,239,0 254,153,95 255,255,255 245,245,0 #F8EF00 #FE995F #FFFFFF #F5F500
        state_color_codes_vector = create_heatmap(data_frame);
        #cat("state_color_codes_vector: ", state_color_codes_vector, "\n");
        #cat("typeof(state_color_codes_vector): ", typeof(state_color_codes_vector), "\n");
        #cat("length(state_color_codes_vector): ", length(state_color_codes_vector), "\n");
        state_colors[[index + 1]] = state_color_codes_vector;
        #cat("\n", index,index,index,index, "state_colors: ", state_colors, "\n");
        index = index + 1;
    }
    #cat("\n\nBBBBBBBB default state_names: ", state_names, "\n");
    #cat("BBBBBBBB default typeof(state_names): ", typeof(state_names), "\n");
    #cat("BBBBBBBB default length(state_names): ", length(state_names), "\n");
    #cat("BBBBBBBB default state_colors: ", toString(state_colors), "\n");
    #cat("BBBBBBBB default length(state_colors): ", length(state_colors), "\n\n\n");
    if (!is.null(specified_state_indexes)) {
        # Replace default name with specified name and default
        # color with specified color for each selected state.
        # Split specified_state_indexes into a list of strings.
        index_str = as.character(specified_state_indexes);
        specified_state_index_items = strsplit(index_str, ",")[[1]];
        # Split specified_state_names into a list of strings.
        name_str = as.character(specified_state_names);
        specified_state_name_items = strsplit(name_str, ",")[[1]];
        # Split specified_state_colors into a list of strings.
        color_str = as.character(specified_state_colors);
        specified_state_color_items = strsplit(color_str, ",")[[1]];
        # Replace default names and colors.
        loop_index = 1;
        for (specified_state_index_item in specified_state_index_items) {
            #cat("specified_state_index_item: ", specified_state_index_item, "\n");
            specified_index = as.integer(specified_state_index_item);
            #cat("specified_index: ", specified_index, "\n");
            # Replace default name with specified name.
            specified_name = specified_state_name_items[loop_index];
            # Handle the special string "use state index".
            if (identical(specified_name, "use state index")) {
                specified_name = specified_state_index_item;
            }
            #cat("specified_name: ", specified_name, "\n");
            #cat("state_names: ", state_names, "\n");
            #state_names[specified_index] = specified_name;
            state_names = replace(state_names, state_names==specified_index, specified_name);
            #cat("state_names: ", state_names, "\n");
            # Replace default color with specified color.
            #cat("specified_state_color_items[specified_item_index]: ", specified_state_color_items[specified_item_index], "\n");
            # FixMe: replacing default state colors is not yet working.
            #state_colors[specified_index] = specified_state_color_items[loop_index];
            #cat("state_colors: ", toString(state_colors), "\n");
            loop_index = loop_index + 1;
        }
    }
    #cat("\n\nCCCCCC reset state_names: ", state_names, "\n");
    #cat("CCCCCC reset state_colors: ", toString(state_colors), "\n");
    #cat("CCCCCC reset length(state_colors): ", length(state_colors), "\n\n\n");
    base_track_file_name <- paste(tracks_dir, hub_name, sep="");
    cells = create_cells(input_dir_state, chrom_len_file, base_track_file_name, state_indexes, state_names, state_colors);
    cell_info = cbind(cells, cells, cells, "#000000");
    cell_info = array(cell_info, dim=c(length(cells), 4));
    cell_info = as.matrix(cell_info);
    track_db = NULL;
    for (i in 1:length(cells)) {
        ii = which(cells[i] == cell_info[,1]);
        if(length(ii) == 0) {
            next;
        }
        ii = ii[1];
        # Get the color for the current state.
        color_vector = state_colors[i];
        #cat("color_vector: ", toString(color_vector), "\n");
        # We are handling hte hub label color setting here.
        # FixMe: Color selection is currently very confusing.  I think we may
        # need to add another color selector for the hub label.  The current
        # approach is to apply the specified color for the state to the hub label
        # if the list of specified_state_index_items contains the value of i - 1
        # since state indexes are zero based.
        #cat("specified_state_index_items: ", specified_state_index_items, "\n");
        #cat("length(specified_state_index_items): ", length(specified_state_index_items), "\n");
        state_index = i - 1;
        #cat("state_index: ", state_index, "\n");
        if (state_index %in% specified_state_index_items) {
            specified_state_color = specified_state_color_items[[i]];
            #cat("specified_state_color: ", specified_state_color, "\n");
            color = paste(c(col2rgb(specified_state_color)), collapse=",");
            #cat("color 1: ", color, "\n");
        } else {
            color = paste(c(col2rgb("#FFFFFF")), collapse=",")
            #cat("color 2: ", color, "\n");
        }
        # Get the bigDataUrl.
        big_data_url <- get_big_data_url(galaxy_url, encoded_dataset_id, tracks_dir, i, build);
        # trackDb.txt track hub entry.
        track_db = c(track_db, paste("track ", hub_name, "_track_", i, sep=""));
        track_db = c(track_db, "type bigBed");
        track_db = c(track_db, paste("bigDataUrl", big_data_url, sep=" "));
        track_db = c(track_db, paste("shortLabel", cell_info[ii, 2], sep=" "));
        track_db = c(track_db, paste("longLabel", paste(hub_name, cell_info[ii, 3], sep=" ")));
        track_db = c(track_db, paste("priority", i));
        track_db = c(track_db, "itemRgb on");
        track_db = c(track_db, "maxItems 100000");
        track_db = c(track_db, paste("color", color, sep=" "));
        track_db = c(track_db, "visibility dense");
        track_db = c(track_db, "");
    }
    return(track_db);
}

get_big_data_url = function(galaxy_url, encoded_dataset_id, tracks_dir, index, build) {
    track_files <- list.files(path=tracks_dir, pattern="\\.bigbed");
    s <- paste(galaxy_url, 'datasets/', encoded_dataset_id, '/display/myHub/', build, '/', track_files[index], sep="");
    return(s)
}

get_cell_type_names = function(state_file) {
    # The first line of a state file is a comment
    # that looks something like this:
    # ID CHR POSst POSed E001 E002 PosClass
    # The cell type names are the elemets whose
    # 1-based indexes start at 5 and end with -1,
    # which in the above case are E001 and E002.
    fh = file(state_file,"r");
    line = readLines(fh, n=1);
    close(fh);
    # Split line into a list of strings.
    items = strsplit(line, "\\s+")[[1]];
    # Extract the cell type names into a list.
    last_cell_type_name_index = length(items) -1;
    cell_type_names = c();
    for (i in 5:last_cell_type_name_index) {
        cell_type_names = c(cell_type_names, items[i]);
    }
    return(cell_type_names);
}

get_track_file_name = function(base_track_file_name, index, ext) {
    track_file_name <- paste(base_track_file_name, index, ext, sep=".");
    return(track_file_name);
}

# Create the hub.txt output.
hub_name_line <- paste("hub ", opt$hub_name, sep="");
hub_short_label_line <- paste("shortLabel ", opt$hub_short_label, sep="");
hub_long_label_line <- paste("longLabel ", opt$hub_long_label, sep="");
genomes_txt_line <- paste("genomesFile genomes.txt", sep="");
email_line <- paste("email ", opt$email, sep="");
contents <- paste(hub_name_line, hub_short_label_line, hub_long_label_line, genomes_txt_line, email_line, sep="\n");
hub_dir <- paste(opt$output_trackhub_files_path, "/", "myHub", "/", sep="");
dir.create(hub_dir, showWarnings=FALSE);
hub_file_path <- paste(hub_dir, "hub.txt", sep="");
write.table(contents, file=hub_file_path, quote=F, row.names=F, col.names=F);

# Create the genomes.txt output.
genome_line <- paste("genome ", opt$build, sep="");
track_db_line <- paste("trackDb ", opt$build, "/", "trackDb.txt", sep="");
contents <- paste(genome_line, track_db_line, sep="\n");
genomes_file_path <- paste(hub_dir, "genomes.txt", sep="");
write.table(contents, file=genomes_file_path, quote=F, row.names=F, col.names=F);

# Create the tracks.
heatmap_path <- paste(opt$script_dir, "create_heatmap.R", sep="/");
source(heatmap_path);
tracks_dir <- paste(hub_dir, opt$build, "/", sep="");
dir.create(tracks_dir, showWarnings=FALSE);
track_db <- create_track_db(opt$galaxy_url, opt$output_trackhub_id, opt$input_dir_para, opt$input_dir_state, opt$build,
        opt$chrom_len_file, tracks_dir, opt$hub_name, opt$hub_short_label, opt$hub_long_label, opt$state_indexes,
        opt$state_names, opt$state_colors);

# Create the trackDb.txt output.
track_db_file_path <- paste(tracks_dir, "trackDb.txt", sep="");
write.table(track_db, file=track_db_file_path, quote=F, row.names=F, col.names=F);

# Create the primary HTML dataset.
create_primary_html(opt$output_trackhub, tracks_dir, opt$build);
