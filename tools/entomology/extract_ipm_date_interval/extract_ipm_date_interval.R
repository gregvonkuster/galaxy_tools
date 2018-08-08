#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("hash"))
suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option(c("--input_data_dir"), action="store", dest="input_data_dir", help="Directory containing .csv outputs from insect_phenology_model"),
    make_option(c("--end_date"), action="store", dest="end_date", help="End date for date interval"),
    make_option(c("--start_date"), action="store", dest="start_date", help="Start date for date interval"),
    make_option(c("--script_dir"), action="store", dest="script_dir", help="R script source directory"),
    make_option(c("--tool_parameters"), action="store", dest="tool_parameters", help="Users defined parameters for executing the insect_phenology_model inputs")
)

parser <- OptionParser(usage="%prog [options] file", option_list=option_list);
args <- parse_args(parser, positional_arguments=TRUE);
opt <- args$options;

get_new_temperature_data_frame = function(input_data_file) {
    # Read a csv file to produce a data frame
    # consisting of the data which was produced
    # by the insect_phenology_model tool.
    temperature_data_frame = read.csv(file=input_data_file, header=T, strip.white=TRUE, stringsAsFactors=FALSE, sep=",");
    return(temperature_data_frame);
}

parse_tool_parameters = function(tool_parameters) {
    # Parse the tool parameters that were used to produce
    # the input datasets found in input_data_dir.  These
    # datasets were produced by the insect_phenology_model
    # tool.
    raw_params = sub("^__SeP__", "", tool_parameters);
    raw_param_items = strsplit(raw_params, "__SeP__")[[1]];
    keys = raw_param_items[c(T, F)];
    values = raw_param_items[c(F, T)];
    num_keys_and_vals = length(keys);
    for (i in 1:num_keys_and_vals) {
        values[i] = restore_text(values[[i]]);
    }
    for (i in 1:num_keys_and_vals) {
        key = keys[i];
        if (endsWith(key, "cond")) {
            value = values[i];
            # Galaxy passes some input job parameters as json-like strings
            # for complex objects like conditionals, so we should see if
            # we can re-implement this using r-jsonlite if possible.  An
            # exception is currently thrown when we do this:
            # params_hash = fromJSON(opt$tool_parameters);
            # Error: lexical error: invalid char in json text.
            #                   __SeP__adult_mortality__SeP____
            # (right here) ------^
            # Here is an example complex object parameter value, in
            # this case the parameter name is plot_nymph_life_stage_cond.
            # {"life_stages_nymph": ["Total"], "__current_case__": 0, "plot_nymph_life_stage": "yes"}
            # This code is somewhat brittle, so a better approach is
            # warranted if possible.
            if (key == "merge_ytd_temperature_data_cond") {
                val = grep("yes", value);
                if (length(val)>0) {
                    # Get the location.
                    items = strsplit(value, "\"location\": ")[[1]];
                    location_str = items[2];
                    val = grep("\",", location_str);
                    if (length(val)>0) {
                        items = strsplit(location_str, "\",")[[1]];
                        location = items[1];
                    } else {
                        location = items[1];
                    }
                    if (location == "\"") {
                        location = "";
                    }
                    keys[i] = "location";
                    values[i] = location;
                }
            } else if (key =="plot_nymph_life_stage_cond") {
                val = grep("yes", value);
                if (length(val)==0) {
                    keys[i] = "plot_nymph_life_stage";
                    values[i] = "no";
                } else {
                    # Get the value for "life_stages_nymph".
                    items = strsplit(value, "\"life_stages_nymph\": ")[[1]];
                    life_stages_nymph_str = items[2];
                    if (grep("],", life_stages_nymph_str)[[1]] > 0) {
                        items = strsplit(life_stages_nymph_str, "],")[[1]];
                        life_stages_nymph_str = items[1];
                        #life_stages_nymph_str = sub("^\\[", "", life_stages_nymph_str);
                        num_curent_keys = length(keys);
                        keys[num_curent_keys+1] = "life_stages_nymph";
                        values[num_curent_keys+1] = life_stages_nymph_str;
                    }
                    keys[i] = "plot_nymph_life_stage";
                    values[i] = "yes";
                }
            } else if (key =="plot_adult_life_stage_cond") {
                val = grep("yes", value);
                # The value of val is an integer if the pattern is not found.
                if (length(val)==0) {
                    keys[i] = "plot_adult_life_stage";
                    values[i] = "no";
                } else {
                    # Get the value for "life_stages_adult".
                    items = strsplit(value, "\"life_stages_adult\": ")[[1]];
                    life_stages_adult_str = items[2];
                    if (grep("],", life_stages_adult_str)[[1]] > 0) {
                        items = strsplit(life_stages_adult_str, "],")[[1]];
                        life_stages_adult_str = items[1];
                        #life_stages_adult_str = sub("^\\[", "", life_stages_adult_str);
                        num_curent_keys = length(keys);
                        keys[num_curent_keys+1] = "life_stages_adult";
                        values[num_curent_keys+1] = life_stages_adult_str;
                    }
                    keys[i] = "plot_adult_life_stage";
                    values[i] = "yes";
                }
            }
        }
    }
    # Strip all double qu0tes from values.
    for (i in 1:length(values)) {
        value = values[i];
        value = gsub("\"", "", value);
        values[i] = value;
    }
    return(hash(keys, values));
}

prepare_plot = function(life_stage, file_path, maxval, ticks, date_labels, chart_type, plot_std_error, insect, location,
    latitude, start_date, end_date, total_days_vector, replications, group, group_std_error, group2, group2_std_error,
    group3, group3_std_error, sub_life_stage=NULL) {
    # Start PDF device driver.
    dev.new(width=20, height=30);
    pdf(file=file_path, width=20, height=30, bg="white");
    par(mar=c(5, 6, 4, 4), mfrow=c(3, 1));
    render_chart(ticks, date_labels, chart_type, plot_std_error, insect, location, latitude, start_date, end_date,
        total_days_vector, maxval, replications, life_stage, group=group, group_std_error=group_std_error, group2=group2,
        group2_std_error=group2_std_error, group3=group3, group3_std_error=group3_std_error, sub_life_stage=sub_life_stage);
    # Turn off device driver to flush output.
    dev.off();
}

restore_text = function(text) {
    # Un-escape characters that are escaped by the
    # Galaxy tool parameter handlers.
    if (is.null(text) || length(text) == 0) {
        return(text);
    }
    chars = list(">", "<", "'", '"', "[", "]", "{", "}", "@", "\n", "\r", "\t", "#");
    mapped_chars = list("__gt__", "__lt__", "__sq__", "__dq__", "__ob__", "__cb__",
                        "__oc__", "__cc__", "__at__", "__cn__", "__cr__", "__tc__", "__pd__");
    for (i in 1:length(mapped_chars)) {
        char = chars[[i]];
        mapped_char = mapped_chars[[i]];
        text = gsub(mapped_char, char, text);
    }
    return(text);
}

# Import the shared utility functions.
utils_path <- paste(opt$script_dir, "utils.R", sep="/");
source(utils_path);

params_hash = parse_tool_parameters(opt$tool_parameters);

# Determine the data we need to generate for plotting.
if (params_hash$plot_generations_separately == "yes") {
    plot_generations_separately = TRUE;
} else {
    plot_generations_separately = FALSE;
}
if (params_hash$plot_std_error == "yes") {
    plot_std_error = TRUE;
} else {
    plot_std_error = FALSE;
}
process_eggs = FALSE;
process_nymphs = FALSE;
process_young_nymphs = FALSE;
process_old_nymphs = FALSE;
process_total_nymphs = FALSE;
process_adults = FALSE;
process_previttelogenic_adults = FALSE;
process_vittelogenic_adults = FALSE;
process_diapausing_adults = FALSE;
process_total_adults = FALSE;
if (params_hash$plot_egg_life_stage == "yes") {
    process_eggs = TRUE;
}
if (params_hash$plot_nymph_life_stage == "yes") {
    process_nymphs = TRUE;
    # Get the selected life stages.
    value = params_hash$life_stages_nymph;
    val = grep("Young", value);
    if (length(val)>0) {
        process_young_nymphs = TRUE;
    }
    val = grep("Old", value);
    if (length(val)>0) {
        process_old_nymphs = TRUE;
    }
    val = grep("Total", value);
    if (length(val)>0) {
        process_total_nymphs = TRUE;
    }
}
if (params_hash$plot_adult_life_stage == "yes") {
    process_adults = TRUE;
    # Get the selected life stages.
    value = params_hash$life_stages_adult;
    val = grep("Pre-vittelogenic", value);
    if (length(val)>0) {
        process_previttelogenic_adults = TRUE;
    }
    val = grep("Vittelogenic", value);
    if (length(val)>0) {
        process_vittelogenic_adults = TRUE;
    }
    val = grep("Diapausing", value);
    if (length(val)>0) {
        process_diapausing_adults = TRUE;
    }
    val = grep("Total", value);
    if (length(val)>0) {
        process_total_adults = TRUE;
    }
}

if (params_hash$plot_egg_life_stage == "yes" & params_hash$plot_nymph_life_stage == "yes" & params_hash$plot_adult_life_stage == "yes") {
    process_total = TRUE;
} else {
    process_total = FALSE;
}


# FIXME: currently custom date fields are free text, but
# Galaxy should soon include support for a date selector
# at which point this tool should be enhanced to use it.
# Validate start_date.
start_date = format(opt$start_date);
end_date = format(opt$end_date);

# Calaculate the number of days in the date interval.
start_date = validate_date(start_date);
# Validate end_date.
end_date = validate_date(end_date);
if (start_date >= end_date) {
    stop_err("The start date must be between 1 and 50 days before the end date when setting date intervals for plots.");
}
# Calculate the number of days in the date interval.
num_days = difftime(end_date, start_date, units=c("days"));
# Add 1 to the number of days to make the dates inclusive.  For
# example, if the user enters a date range of 2018-01-01 to
# 2018-01-31, they likely expect the end date to be included.
num_days = num_days + 1;
if (num_days > 50) {
    # We need to restrict date intervals since
    # plots render tick marks for each day.
    stop_err("Date intervals for plotting cannot exceed 50 days.");
}
# Display the total number of days in the Galaxy history item blurb.
cat("Number of days in date interval: ", num_days, "\n");

# Create the csv data files consisting of the date interval.
input_data_files = list.files(path=opt$input_data_dir, full.names=TRUE);
for (input_data_file in input_data_files) {
    file_name = basename(input_data_file);
    temperature_data_frame = get_new_temperature_data_frame(input_data_file);
    start_date_row = which(temperature_data_frame$DATE==start_date);
    end_date_row = which(temperature_data_frame$DATE==end_date);
    # Extract the date interval.
    temperature_data_frame = temperature_data_frame[start_date_row:end_date_row,];
    # Save the date interval data into an output file
    # named the same as the input.
    file_path = paste("output_data_dir", file_name, sep="/");
    write.csv(temperature_data_frame, file=file_path, row.names=F);
}

# Extract the vectors needed for plots from the input data files
# produced by the insect_phenology_model tool.
total_days_vector = NULL;
ticks_and_labels = NULL;
latitude = NULL;
input_data_files = list.files(path="output_data_dir", full.names=TRUE);
for (input_data_file in input_data_files) {
    file_name = basename(input_data_file);
    temperature_data_frame = get_new_temperature_data_frame(input_data_file);
    # Initialize the total_days_vector for later plotting.
    if (is.null(total_days_vector)) {
        total_days_vector = c(1:dim(temperature_data_frame)[1]);
    }
    if (is.null(ticks_and_labels)) {
        # Get the ticks date labels for later plotting
        ticks_and_labels = get_x_axis_ticks_and_labels(temperature_data_frame, date_interval=TRUE);
        ticks = c(unlist(ticks_and_labels[1]));
        date_labels = c(unlist(ticks_and_labels[2]));
    }
    if (is.null(latitude)) {
        # Get the latitude for later plotting.
        latitude = temperature_data_frame$LATITUDE[1];
    }

    if (file_name == "04_combined_generations.csv") {
        if (process_eggs) {
            eggs = temperature_data_frame$EGG;
            if (plot_std_error) {
                eggs.std_error = temperature_data_frame$EGGSE;
            }
        }
        if (process_young_nymphs) {
            young_nymphs = temperature_data_frame$YOUNGNYMPH;
            if (plot_std_error) {
                young_nymphs.std_error = temperature_data_frame$YOUNGNYMPHSE;
            }
        }
        if (process_old_nymphs) {
            old_nymphs = temperature_data_frame$OLDNYMPH;
            if (plot_std_error) {
                old_nymphs.std_error = temperature_data_frame$OLDNYMPHSE;
            }
        }
        if (process_total_nymphs) {
            total_nymphs = temperature_data_frame$TOTALNYMPH;
            if (plot_std_error) {
                total_nymphs.std_error = temperature_data_frame$TOTALNYMPHSE;
            }
        }
        if (process_previttelogenic_adults) {
            previttelogenic_adults = temperature_data_frame$PRE.VITADULT;
            if (plot_std_error) {
                previttelogenic_adults.std_error = temperature_data_frame$PRE.VITADULTSE;
            }
        }
        if (process_vittelogenic_adults) {
            vittelogenic_adults = temperature_data_frame$VITADULT;
            if (plot_std_error) {
                vittelogenic_adults.std_error = temperature_data_frame$VITADULTSE;
            }
        }
        if (process_diapausing_adults) {
            diapausing_adults = temperature_data_frame$DIAPAUSINGADULT;
            if (plot_std_error) {
                diapausing_adults.std_error = temperature_data_frame$DIAPAUSINGADULTSE;
            }
        }
        if (process_total_adults) {
            total_adults = temperature_data_frame$TOTALADULT;
            if (plot_std_error) {
                total_adults.std_error = temperature_data_frame$TOTALADULTSE;
            }
        }
    } else if (file_name == "01_generation_P.csv") {
        if (process_eggs) {
            P_eggs = temperature_data_frame$EGG.P;
            if (plot_std_error) {
                P_eggs.std_error = temperature_data_frame$EGG.P.SE;
            }
        }
        if (process_young_nymphs) {
            P_young_nymphs = temperature_data_frame$YOUNGNYMPH.P;
            if (plot_std_error) {
                P_young_nymphs.std_error = temperature_data_frame$YOUNGNYMPH.P.SE;
            }
        }
        if (process_old_nymphs) {
            P_old_nymphs = temperature_data_frame$OLDNYMPH.P;
            if (plot_std_error) {
                P_old_nymphs.std_error = temperature_data_frame$OLDNYMPH.P.SE;
            }
        }
        if (process_total_nymphs) {
            P_total_nymphs = temperature_data_frame$TOTALNYMPH.P;
            if (plot_std_error) {
                P_total_nymphs.std_error = temperature_data_frame$TOTALNYMPH.P.SE;
            }
        }
        if (process_previttelogenic_adults) {
            P_previttelogenic_adults = temperature_data_frame$PRE.VITADULT.P;
            if (plot_std_error) {
                P_previttelogenic_adults.std_error = temperature_data_frame$PRE.VITADULT.P.SE;
            }
        }
        if (process_vittelogenic_adults) {
            P_vittelogenic_adults = temperature_data_frame$VITADULT.P;
            if (plot_std_error) {
                P_vittelogenic_adults.std_error = temperature_data_frame$VITADULT.P.SE;
            }
        }
        if (process_diapausing_adults) {
            P_diapausing_adults = temperature_data_frame$DIAPAUSINGADULT.P;
            if (plot_std_error) {
                P_diapausing_adults.std_error = temperature_data_frame$DIAPAUSINGADULT.P.SE;
            }
        }
        if (process_total_adults) {
            P_total_adults = temperature_data_frame$TOTALADULT.P;
            if (plot_std_error) {
                P_total_adults.std_error = temperature_data_frame$TOTALADULT.P.SE;
            }
        }
        if (process_total) {
            P_all_total = temperature_data_frame$ALL.TOTAL.P;
            if (plot_std_error) {
                P_all_total.std_error = temperature_data_frame$ALL.TOTAL.P.SE;
            }
        }
    } else if (file_name == "02_generation_F1.csv") {
        if (process_eggs) {
            F1_eggs = temperature_data_frame$EGG.F1;
            if (plot_std_error) {
                F1_eggs.std_error = temperature_data_frame$EGG.F1.SE;
            }
        }
        if (process_young_nymphs) {
            F1_young_nymphs = temperature_data_frame$YOUNGNYMPH.F1;
            if (plot_std_error) {
                F1_young_nymphs.std_error = temperature_data_frame$YOUNGNYMPH.F1.SE;
            }
        }
        if (process_old_nymphs) {
            F1_old_nymphs = temperature_data_frame$OLDNYMPH.F1;
            if (plot_std_error) {
                F1_old_nymphs.std_error = temperature_data_frame$OLDNYMPH.F1.SE;
            }
        }
        if (process_total_nymphs) {
            F1_total_nymphs = temperature_data_frame$TOTALNYMPH.F1;
            if (plot_std_error) {
                F1_total_nymphs.std_error = temperature_data_frame$TOTALNYMPH.F1.SE;
            }
        }
        if (process_previttelogenic_adults) {
            F1_previttelogenic_adults = temperature_data_frame$PRE.VITADULT.F1;
            if (plot_std_error) {
                F1_previttelogenic_adults.std_error = temperature_data_frame$PRE.VITADULT.F1.SE;
            }
        }
        if (process_vittelogenic_adults) {
            F1_vittelogenic_adults = temperature_data_frame$VITADULT.F1;
            if (plot_std_error) {
                F1_vittelogenic_adults.std_error = temperature_data_frame$VITADULT.F1.SE;
            }
        }
        if (process_diapausing_adults) {
            F1_diapausing_adults = temperature_data_frame$DIAPAUSINGADULT.F1;
            if (plot_std_error) {
                F1_diapausing_adults.std_error = temperature_data_frame$DIAPAUSINGADULT.F1.SE;
            }
        }
        if (process_total_adults) {
            F1_total_adults = temperature_data_frame$TOTALADULT.F1;
            if (plot_std_error) {
                F1_total_adults.std_error = temperature_data_frame$TOTALADULT.F1.SE;
            }
        }
        if (process_total) {
            F1_all_total = temperature_data_frame$ALL.TOTAL.F1;
            if (plot_std_error) {
                F1_all_total.std_error = temperature_data_frame$ALL.TOTAL.F1.SE;
            }
        }
    } else if (file_name == "03_generation_F2.csv") {
        if (process_eggs) {
            F2_eggs = temperature_data_frame$EGG.F2;
            if (plot_std_error) {
                F2_eggs.std_error = temperature_data_frame$EGG.F2.SE;
            }
        }
        if (process_young_nymphs) {
            F2_young_nymphs = temperature_data_frame$YOUNGNYMPH.F2;
            if (plot_std_error) {
                F2_young_nymphs.std_error = temperature_data_frame$YOUNGNYMPH.F2.SE;
            }
        }
        if (process_old_nymphs) {
            F2_old_nymphs = temperature_data_frame$OLDNYMPH.F2;
            if (plot_std_error) {
                F2_old_nymphs.std_error = temperature_data_frame$OLDNYMPH.F2.SE;
            }
        }
        if (process_total_nymphs) {
            F2_total_nymphs = temperature_data_frame$TOTALNYMPH.F2;
            if (plot_std_error) {
                F2_total_nymphs.std_error = temperature_data_frame$TOTALNYMPH.F2.SE;
            }
        }
        if (process_previttelogenic_adults) {
            F2_previttelogenic_adults = temperature_data_frame$PRE.VITADULT.F2;
            if (plot_std_error) {
                F2_previttelogenic_adults.std_error = temperature_data_frame$PRE.VITADULT.F2.SE;
            }
        }
        if (process_vittelogenic_adults) {
            F2_vittelogenic_adults = temperature_data_frame$VITADULT.F2;
            if (plot_std_error) {
                F2_vittelogenic_adults.std_error = temperature_data_frame$VITADULT.F2.SE;
            }
        }
        if (process_diapausing_adults) {
            F2_diapausing_adults = temperature_data_frame$DIAPAUSINGADULT.F2;
            if (plot_std_error) {
                F2_diapausing_adults.std_error = temperature_data_frame$DIAPAUSINGADULT.F2.SE;
            }
        }
        if (process_total_adults) {
            F2_total_adults = temperature_data_frame$TOTALADULT.F2;
            if (plot_std_error) {
                F2_total_adults.std_error = temperature_data_frame$TOTALADULT.F2.SE;
            }
        }
        if (process_total) {
            F2_all_total = temperature_data_frame$ALL.TOTAL.F2;
            if (plot_std_error) {
                F2_all_total.std_error = temperature_data_frame$ALL.TOTAL.F2.SE;
            }
        }
    }
}

# Create the pdf plot files based on the date interval.
if (plot_generations_separately) {
    chart_type = "pop_size_by_generation";
    if (process_eggs) {
        # Total population size by generation.
        life_stage = "Egg";
        file_path = get_file_path(life_stage, "egg_pop_by_generation.pdf")
        maxval = max(P_eggs+F1_eggs+F2_eggs) + 100;
        prepare_plot(life_stage, file_path, maxval, ticks, date_labels, chart_type, params_hash$plot_std_error,
            params_hash$insect, params_hash$location, latitude, start_date, end_date, total_days_vector,
            params_hash$replications, group=P_eggs, group_std_error=P_eggs.std_error, group2=F1_eggs,
            group2_std_error=F1_eggs.std_error, group3=F2_eggs, group3_std_error=F2_eggs.std_error);
    }
    if (process_nymphs) {
        life_stage = "Nymph";
        if (process_young_nymphs) {
            # Young nymph population size by generation.
            sub_life_stage = "Young";
            file_path = get_file_path(life_stage, "nymph_pop_by_generation.pdf", sub_life_stage=sub_life_stage)
            maxval = max(P_young_nymphs+F1_young_nymphs+F2_young_nymphs) + 100;
            prepare_plot(life_stage, file_path, maxval, ticks, date_labels, chart_type, params_hash$plot_std_error,
                params_hash$insect, params_hash$location, latitude, start_date, end_date, total_days_vector,
                params_hash$replications, group=P_young_nymphs, group_std_error=P_young_nymphs.std_error,
                group2=F1_young_nymphs, group2_std_error=F1_young_nymphs.std_error, group3=F2_young_nymphs,
                group3_std_error=F2_young_nymphs.std_error, sub_life_stage=sub_life_stage);
        }
        if (process_old_nymphs) {
            # Old nymph population size by generation.
            sub_life_stage = "Old";
            file_path = get_file_path(life_stage, "nymph_pop_by_generation.pdf", sub_life_stage=sub_life_stage)
            maxval = max(P_old_nymphs+F1_old_nymphs+F2_old_nymphs) + 100;
            prepare_plot(life_stage, file_path, maxval, ticks, date_labels, chart_type, params_hash$plot_std_error,
                params_hash$insect, params_hash$location, latitude, start_date, end_date, total_days_vector,
                params_hash$replications, group=P_old_nymphs, group_std_error=P_old_nymphs.std_error,
                group2=F1_old_nymphs, group2_std_error=F1_old_nymphs.std_error, group3=F2_old_nymphs,
                group3_std_error=F2_old_nymphs.std_error, sub_life_stage=sub_life_stage);
        }
        if (process_total_nymphs) {
            # Total nymph population size by generation.
            sub_life_stage = "Total";
            file_path = get_file_path(life_stage, "nymph_pop_by_generation.pdf", sub_life_stage=sub_life_stage)
            maxval = max(P_total_nymphs+F1_total_nymphs+F2_total_nymphs) + 100;
            prepare_plot(life_stage, file_path, maxval, ticks, date_labels, chart_type, params_hash$plot_std_error,
                params_hash$insect, params_hash$location, latitude, start_date, end_date, total_days_vector,
                params_hash$replications, group=P_total_nymphs, group_std_error=P_total_nymphs.std_error,
                group2=F1_total_nymphs, group2_std_error=F1_total_nymphs.std_error, group3=F2_total_nymphs,
                group3_std_error=F2_total_nymphs.std_error, sub_life_stage=sub_life_stage);
        }
    }
    if (process_adults) {
        life_stage = "Adult";
        if (process_previttelogenic_adults) {
            # Pre-vittelogenic adult population size by generation.
            sub_life_stage = "Pre-vittelogenic";
            file_path = get_file_path(life_stage, "adult_pop_by_generation.pdf", sub_life_stage=sub_life_stage)
            maxval = max(P_previttelogenic_adults+F1_previttelogenic_adults+F2_previttelogenic_adults) + 100;
            prepare_plot(life_stage, file_path, maxval, ticks, date_labels, chart_type, params_hash$plot_std_error,
                params_hash$insect, params_hash$location, latitude, start_date, end_date, total_days_vector,
                params_hash$replications, group=P_previttelogenic_adults,
                group_std_error=P_previttelogenic_adults.std_error, group2=F1_previttelogenic_adults,
                group2_std_error=F1_previttelogenic_adults.std_error, group3=F2_previttelogenic_adults,
                group3_std_error=F2_previttelogenic_adults.std_error, sub_life_stage=sub_life_stage);
        }
        if (process_vittelogenic_adults) {
            # Vittelogenic adult population size by generation.
            sub_life_stage = "Vittelogenic";
            file_path = get_file_path(life_stage, "adult_pop_by_generation.pdf", sub_life_stage=sub_life_stage)
            maxval = max(P_vittelogenic_adults+F1_vittelogenic_adults+F2_vittelogenic_adults) + 100;
            prepare_plot(life_stage, file_path, maxval, ticks, date_labels, chart_type, params_hash$plot_std_error,
                params_hash$insect, params_hash$location, latitude, start_date, end_date, total_days_vector,
                params_hash$replications, group=P_vittelogenic_adults,
                group_std_error=P_vittelogenic_adults.std_error, group2=F1_vittelogenic_adults,
                group2_std_error=F1_vittelogenic_adults.std_error, group3=F2_vittelogenic_adults,
                group3_std_error=F2_vittelogenic_adults.std_error, sub_life_stage=sub_life_stage);
        }
        if (process_diapausing_adults) {
            # Diapausing adult population size by generation.
            sub_life_stage = "Diapausing";
            file_path = get_file_path(life_stage, "adult_pop_by_generation.pdf", sub_life_stage=sub_life_stage)
            maxval = max(P_diapausing_adults+F1_diapausing_adults+F2_diapausing_adults) + 100;
            prepare_plot(life_stage, file_path, maxval, ticks, date_labels, chart_type, params_hash$plot_std_error,
                params_hash$insect, params_hash$location, latitude, start_date, end_date, total_days_vector,
                params_hash$replications, group=P_diapausing_adults, group_std_error=P_diapausing_adults.std_error,
                group2=F1_diapausing_adults, group2_std_error=F1_diapausing_adults.std_error, group3=F2_diapausing_adults,
                group3_std_error=F2_diapausing_adults.std_error, sub_life_stage=sub_life_stage);
        }
        if (process_total_adults) {
            # Total adult population size by generation.
            sub_life_stage = "Total";
            file_path = get_file_path(life_stage, "adult_pop_by_generation.pdf", sub_life_stage=sub_life_stage)
            maxval = max(P_total_adults+F1_total_adults+F2_total_adults) + 100;
            prepare_plot(life_stage, file_path, maxval, ticks, date_labels, chart_type, params_hash$plot_std_error,
                params_hash$insect, params_hash$location, latitude, start_date, end_date, total_days_vector,
                params_hash$replications, group=P_total_adults, group_std_error=P_total_adults.std_error,
                group2=F1_total_adults, group2_std_error=F1_total_adults.std_error, group3=F2_total_adults,
                group3_std_error=F2_total_adults.std_error, sub_life_stage=sub_life_stage);
        }
    }
    if (process_total) {
        life_stage = "Total";
        # Total population size for egg, nymph and adult by generation.
        file_path = get_file_path(life_stage, "total_pop_by_generation.pdf")
        maxval = max(total_adults+eggs+total_nymphs) + 100;
        prepare_plot(life_stage, file_path, maxval, ticks, date_labels, chart_type, params_hash$plot_std_error,
            params_hash$insect, params_hash$location, latitude, start_date, end_date, total_days_vector,
            params_hash$replications, group=P_all_total, group_std_error=P_all_total.std_error, group2=F1_all_total,
            group2_std_error=F1_all_total.std_error, group3=F2_all_total, group3_std_error=F2_all_total.std_error);
    }
} else {
    chart_type = "pop_size_by_life_stage";
    if (process_eggs) {
        # Egg population size.
        life_stage = "Egg";
        file_path = get_file_path(life_stage, "egg_pop.pdf")
        maxval = max(eggs+eggs.std_error) + 100;
        prepare_plot(life_stage, file_path, maxval, ticks, date_labels, chart_type, params_hash$plot_std_error,
            params_hash$insect, params_hash$location, latitude, start_date, end_date, total_days_vector,
            params_hash$replications, group=eggs, group_std_error=eggs.std_error);
    }
    if (process_nymphs) {
        life_stage = "Nymph";
        if (process_young_nymphs) {
            # Young nymph population size.
            sub_life_stage = "Young";
            file_path = get_file_path(life_stage, "nymph_pop.pdf", sub_life_stage=sub_life_stage)
            maxval = max(young_nymphs+young_nymphs.std_error) + 100;
            prepare_plot(life_stage, file_path, maxval, ticks, date_labels, chart_type, params_hash$plot_std_error,
                params_hash$insect, params_hash$location, latitude, start_date, end_date, total_days_vector,
                params_hash$replications, group=young_nymphs, group_std_error=young_nymphs.std_error,
                sub_life_stage=sub_life_stage);
        }
        if (process_old_nymphs) {
            # Old nymph population size.
            sub_life_stage = "Old";
            file_path = get_file_path(life_stage, "nymph_pop.pdf", sub_life_stage=sub_life_stage)
            maxval = max(old_nymphs+old_nymphs.std_error) + 100;
            prepare_plot(life_stage, file_path, maxval, ticks, date_labels, chart_type, params_hash$plot_std_error,
                params_hash$insect, params_hash$location, latitude, start_date, end_date, total_days_vector,
                params_hash$replications, group=old_nymphs, group_std_error=old_nymphs.std_error,
                sub_life_stage=sub_life_stage);
        }
        if (process_total_nymphs) {
            # Total nymph population size.
            sub_life_stage = "Total";
            file_path = get_file_path(life_stage, "nymph_pop.pdf", sub_life_stage=sub_life_stage)
            maxval = max(total_nymphs+total_nymphs.std_error) + 100;
            prepare_plot(life_stage, file_path, maxval, ticks, date_labels, chart_type, params_hash$plot_std_error,
                params_hash$insect, params_hash$location, latitude, start_date, end_date, total_days_vector,
                params_hash$replications, group=total_nymphs, group_std_error=total_nymphs.std_error,
                sub_life_stage=sub_life_stage);
        }
    }
    if (process_adults) {
        life_stage = "Adult";
        if (process_previttelogenic_adults) {
            # Pre-vittelogenic adult population size.
            sub_life_stage = "Pre-vittelogenic";
            file_path = get_file_path(life_stage, "adult_pop.pdf", sub_life_stage=sub_life_stage)
            maxval = max(previttelogenic_adults+previttelogenic_adults.std_error) + 100;
            prepare_plot(life_stage, file_path, maxval, ticks, date_labels, chart_type, params_hash$plot_std_error,
                params_hash$insect, params_hash$location, latitude, start_date, end_date, total_days_vector,
                params_hash$replications, group=previttelogenic_adults,
                group_std_error=previttelogenic_adults.std_error, sub_life_stage=sub_life_stage);
        }
        if (process_vittelogenic_adults) {
            # Vittelogenic adult population size.
            sub_life_stage = "Vittelogenic";
            file_path = get_file_path(life_stage, "adult_pop.pdf", sub_life_stage=sub_life_stage)
            maxval = max(vittelogenic_adults+vittelogenic_adults.std_error) + 100;
            prepare_plot(life_stage, file_path, maxval, ticks, date_labels, chart_type, params_hash$plot_std_error,
                params_hash$insect, params_hash$location, latitude, start_date, end_date, total_days_vector,
                params_hash$replications, group=vittelogenic_adults,
                group_std_error=vittelogenic_adults.std_error, sub_life_stage=sub_life_stage);
        }
        if (process_diapausing_adults) {
            # Diapausing adult population size.
            sub_life_stage = "Diapausing";
            file_path = get_file_path(life_stage, "adult_pop.pdf", sub_life_stage=sub_life_stage)
            maxval = max(diapausing_adults+diapausing_adults.std_error) + 100;
            prepare_plot(life_stage, file_path, maxval, ticks, date_labels, chart_type, params_hash$plot_std_error,
                params_hash$insect, params_hash$location, latitude, start_date, end_date, total_days_vector,
                params_hash$replications, group=diapausing_adults, group_std_error=diapausing_adults.std_error,
                sub_life_stage=sub_life_stage);
        }
        if (process_total_adults) {
            # Total adult population size.
            sub_life_stage = "Total";
            file_path = get_file_path(life_stage, "adult_pop.pdf", sub_life_stage=sub_life_stage)
            maxval = max(total_adults+total_adults.std_error) + 100;
            prepare_plot(life_stage, file_path, maxval, ticks, date_labels, chart_type, params_hash$plot_std_error,
                params_hash$insect, params_hash$location, latitude, start_date, end_date, total_days_vector,
                params_hash$replications, group=total_adults, group_std_error=total_adults.std_error,
                sub_life_stage=sub_life_stage);
        }
    }
    if (process_total) {
        # Total population size.
        life_stage = "Total";
        file_path = get_file_path(life_stage, "total_pop.pdf")
        maxval = max(eggs+eggs.std_error, total_nymphs+total_nymphs.std_error, total_adults+total_adults.std_error) + 100;
        prepare_plot(life_stage, file_path, maxval, ticks, date_labels, chart_type, params_hash$plot_std_error,
            params_hash$insect, params_hash$location, latitude, start_date, end_date, total_days_vector,
            params_hash$replications, group=total_adults, group_std_error=total_adults.std_error,
            group2=total_nymphs, group2_std_error=total_nymphs.std_error, group3=eggs, group3_std_error=eggs.std_error);
    }
}

