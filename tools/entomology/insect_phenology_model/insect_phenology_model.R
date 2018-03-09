#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option(c("--adult_mortality"), action="store", dest="adult_mortality", type="integer", help="Adjustment rate for adult mortality"),
    make_option(c("--adult_accumulation"), action="store", dest="adult_accumulation", type="integer", help="Adjustment of degree-days accumulation (old nymph->adult)"),
    make_option(c("--egg_mortality"), action="store", dest="egg_mortality", type="integer", help="Adjustment rate for egg mortality"),
    make_option(c("--input"), action="store", dest="input", help="Temperature data for selected location"),
    make_option(c("--insect"), action="store", dest="insect", help="Insect name"),
    make_option(c("--insects_per_replication"), action="store", dest="insects_per_replication", type="integer", help="Number of insects with which to start each replication"),
    make_option(c("--life_stages"), action="store", dest="life_stages", help="Selected life stages for plotting"),
    make_option(c("--life_stages_adult"), action="store", dest="life_stages_adult", default=NULL, help="Adult life stages for plotting"),
    make_option(c("--life_stages_nymph"), action="store", dest="life_stages_nymph", default=NULL, help="Nymph life stages for plotting"),
    make_option(c("--location"), action="store", dest="location", help="Selected location"),
    make_option(c("--min_clutch_size"), action="store", dest="min_clutch_size", type="integer", help="Adjustment of minimum clutch size"),
    make_option(c("--max_clutch_size"), action="store", dest="max_clutch_size", type="integer", help="Adjustment of maximum clutch size"),
    make_option(c("--num_days"), action="store", dest="num_days", type="integer", help="Total number of days in the temperature dataset"),
    make_option(c("--nymph_mortality"), action="store", dest="nymph_mortality", type="integer", help="Adjustment rate for nymph mortality"),
    make_option(c("--old_nymph_accumulation"), action="store", dest="old_nymph_accumulation", type="integer", help="Adjustment of degree-days accumulation (young nymph->old nymph)"),
    make_option(c("--output"), action="store", dest="output", help="Dataset containing analyzed data"),
    make_option(c("--oviposition"), action="store", dest="oviposition", type="integer", help="Adjustment for oviposition rate"),
    make_option(c("--photoperiod"), action="store", dest="photoperiod", type="double", help="Critical photoperiod for diapause induction/termination"),
    make_option(c("--plot_generations_separately"), action="store", dest="plot_generations_separately", help="Plot Plot P, F1 and F2 as separate lines or pool across them"),
    make_option(c("--plot_std_error"), action="store", dest="plot_std_error", help="Plot Standard error"),
    make_option(c("--replications"), action="store", dest="replications", type="integer", help="Number of replications"),
    make_option(c("--young_nymph_accumulation"), action="store", dest="young_nymph_accumulation", type="integer", help="Adjustment of degree-days accumulation (egg->young nymph)")
)

parser <- OptionParser(usage="%prog [options] file", option_list=option_list);
args <- parse_args(parser, positional_arguments=TRUE);
opt <- args$options;

add_daylight_length = function(temperature_data_frame, num_rows) {
    # Return a vector of daylight length (photoperido profile) for
    # the number of days specified in the input temperature data
    # (from Forsythe 1995).
    p = 0.8333;
    latitude = temperature_data_frame$LATITUDE[1];
    daylight_length_vector = NULL;
    for (i in 1:num_rows) {
        # Get the day of the year from the current row
        # of the temperature data for computation.
        doy = temperature_data_frame$DOY[i];
        theta = 0.2163108 + 2 * atan(0.9671396 * tan(0.00860 * (doy - 186)));
        phi = asin(0.39795 * cos(theta));
        # Compute the length of daylight for the day of the year.
        darkness_length = 24 / pi * acos((sin(p * pi / 180) + sin(latitude * pi / 180) * sin(phi)) / (cos(latitude * pi / 180) * cos(phi)));
        daylight_length_vector[i] = 24 - darkness_length;
    }
    # Append daylight_length_vector as a new column to temperature_data_frame.
    temperature_data_frame = append_vector(temperature_data_frame, daylight_length_vector, "DAYLEN");
    return(temperature_data_frame);
}

append_vector = function(data_frame, vec, new_column_name) {
    num_columns = dim(data_frame)[2];
    current_column_names = colnames(data_frame);
    # Append vector vec as a new column to data_frame.
    data_frame[,num_columns+1] = vec;
    # Reset the column names with the additional column for later access.
    colnames(data_frame) = append(current_column_names, new_column_name);
    return(data_frame);
}

get_date_labels = function(temperature_data_frame, num_rows) {
    # Keep track of the years to see if spanning years.
    month_labels = list();
    current_month_label = NULL;
    for (i in 1:num_rows) {
        # Get the year and month from the date which
        # has the format YYYY-MM-DD.
        date = format(temperature_data_frame$DATE[i]);
        items = strsplit(date, "-")[[1]];
        month = items[2];
        month_label = month.abb[as.integer(month)];
        if (!identical(current_month_label, month_label)) {
            month_labels[length(month_labels)+1] = month_label;
            current_month_label = month_label;
        }
    }
    return(c(unlist(month_labels)));
}

get_file_path = function(life_stage, base_name, life_stage_nymph=NULL, life_stage_adult=NULL) {
    if (!is.null(life_stage_nymph)) {
        lsi = get_life_stage_index(life_stage, life_stage_nymph=life_stage_nymph);
        file_name = paste(lsi, tolower(life_stage_nymph), base_name, sep="_");
    } else if (!is.null(life_stage_adult)) {
        lsi = get_life_stage_index(life_stage, life_stage_adult=life_stage_adult);
        file_name = paste(lsi, tolower(life_stage_adult), base_name, sep="_");
    } else {
        lsi = get_life_stage_index(life_stage);
        file_name = paste(lsi, base_name, sep="_");
    }
    file_path = paste("output_dir", file_name, sep="/");
    return(file_path);
}

get_life_stage_index = function(life_stage, life_stage_nymph=NULL, life_stage_adult=NULL) {
    # Name collection elements so that they
    # are displayed in logical order.
    if (life_stage=="Egg") {
        lsi = "01";
    } else if (life_stage=="Nymph") {
        if (life_stage_nymph=="Young") {
            lsi = "02";
        } else if (life_stage_nymph=="Old") {
            lsi = "03";
        } else if (life_stage_nymph=="Total") {
            lsi="04";
        }
    } else if (life_stage=="Adult") {
        if (life_stage_adult=="Pre-vittelogenic") {
            lsi = "05";
        } else if (life_stage_adult=="Vittelogenic") {
            lsi = "06";
        } else if (life_stage_adult=="Diapausing") {
            lsi = "07";
        } else if (life_stage_adult=="Total") {
            lsi = "08";
        }
    } else if (life_stage=="Total") {
        lsi = "09";
    }
    return(lsi);
}

get_mean_and_std_error = function(p_replications, f1_replications, f2_replications) {
    # P mean.
    p_m = apply(p_replications, 1, mean);
    # P standard error.
    p_se = apply(p_replications, 1, sd) / sqrt(opt$replications);
    # F1 mean.
    f1_m = apply(f1_replications, 1, mean);
    # F1 standard error.
    f1_se = apply(f1_replications, 1, sd) / sqrt(opt$replications);
    # F2 mean.
    f2_m = apply(f2_replications, 1, mean);
    # F2 standard error.
    f2_se = apply(f2_replications, 1, sd) / sqrt(opt$replications);
    return(list(p_m, p_se, f1_m, f1_se, f2_m, f2_se))
}

get_temperature_at_hour = function(latitude, temperature_data_frame, row, num_days) {
    # Base development threshold for Brown Marmorated Stink Bug
    # insect phenology model.
    threshold = 14.17;
    # Minimum temperature for current row.
    curr_min_temp = temperature_data_frame$TMIN[row];
    # Maximum temperature for current row.
    curr_max_temp = temperature_data_frame$TMAX[row];
    # Mean temperature for current row.
    curr_mean_temp = 0.5 * (curr_min_temp + curr_max_temp);
    # Initialize degree day accumulation
    averages = 0;
    if (curr_max_temp < threshold) {
        averages = 0;
    }
    else {
        # Initialize hourly temperature.
        T = NULL;
        # Initialize degree hour vector.
        dh = NULL;
        # Daylight length for current row.
        y = temperature_data_frame$DAYLEN[row];
        # Darkness length.
        z = 24 - y;
        # Lag coefficient.
        a = 1.86;
        # Darkness coefficient.
        b = 2.20;
        # Sunrise time.
        risetime = 12 - y / 2;
        # Sunset time.
        settime = 12 + y / 2;
        ts = (curr_max_temp - curr_min_temp) * sin(pi * (settime - 5) / (y + 2 * a)) + curr_min_temp;
        for (i in 1:24) {
            if (i > risetime && i < settime) {
                # Number of hours after Tmin until sunset.
                m = i - 5;
                T[i] = (curr_max_temp - curr_min_temp) * sin(pi * m / (y + 2 * a)) + curr_min_temp;
                if (T[i] < 8.4) {
                    dh[i] = 0;
                }
                else {
                    dh[i] = T[i] - 8.4;
                }
            }
            else if (i > settime) {
                n = i - settime;
                T[i] = curr_min_temp + (ts - curr_min_temp) * exp( - b * n / z);
                if (T[i] < 8.4) {
                    dh[i] = 0;
                }
                else {
                    dh[i] = T[i] - 8.4;
                }
            }
            else {
                n = i + 24 - settime;
                T[i] = curr_min_temp + (ts - curr_min_temp) * exp( - b * n / z);
                if (T[i] < 8.4) {
                    dh[i] = 0;
                }
                else {
                    dh[i] = T[i] - 8.4;
                }
            }
        }
        averages = sum(dh) / 24;
    }
    return(c(curr_mean_temp, averages))
}

mortality.adult = function(temperature) {
    if (temperature < 12.7) {
        mortality.probability = 0.002;
    }
    else {
        mortality.probability = temperature * 0.0005 + 0.02;
    }
    return(mortality.probability)
}

mortality.egg = function(temperature) {
    if (temperature < 12.7) {
        mortality.probability = 0.8;
    }
    else {
        mortality.probability = 0.8 - temperature / 40.0;
        if (mortality.probability < 0) {
            mortality.probability = 0.01;
        }
    }
    return(mortality.probability)
}

mortality.nymph = function(temperature) {
    if (temperature < 12.7) {
        mortality.probability = 0.03;
    }
    else {
        mortality.probability = temperature * 0.0008 + 0.03;
    }
    return(mortality.probability);
}

parse_input_data = function(input_file, num_rows) {
    # Read in the input temperature datafile into a data frame.
    temperature_data_frame = read.csv(file=input_file, header=T, strip.white=TRUE, sep=",");
    num_columns = dim(temperature_data_frame)[2];
    if (num_columns == 6) {
        # The input data has the following 6 columns:
        # LATITUDE, LONGITUDE, DATE, DOY, TMIN, TMAX
        # Set the column names for access when adding daylight length..
        colnames(temperature_data_frame) = c("LATITUDE","LONGITUDE", "DATE", "DOY", "TMIN", "TMAX");
        current_column_names = colnames(temperature_data_frame);
        # Add a column containing the daylight length for each day.
        temperature_data_frame = add_daylight_length(temperature_data_frame, num_rows);
    }
    return(temperature_data_frame);
}

render_chart = function(date_labels, chart_type, plot_std_error, insect, location, latitude, start_date, end_date, days, maxval,
    replications, life_stage, group, group_std_error, group2=NULL, group2_std_error=NULL, group3=NULL, group3_std_error=NULL,
    life_stages_adult=NULL, life_stages_nymph=NULL) {
    if (chart_type=="pop_size_by_life_stage") {
        if (life_stage=="Total") {
            title = paste(insect, ": Reps", replications, ":", life_stage, "Pop :", location, ": Lat", latitude, ":", start_date, "-", end_date, sep=" ");
            legend_text = c("Egg", "Nymph", "Adult");
            columns = c(4, 2, 1);
            plot(days, group, main=title, type="l", ylim=c(0, maxval), axes=F, lwd=2, xlab="", ylab="", cex=3, cex.lab=3, cex.axis=3, cex.main=3);
            legend("topleft", legend_text, lty=c(1, 1, 1), col=columns, cex=3);
            lines(days, group2, lwd=2, lty=1, col=2);
            lines(days, group3, lwd=2, lty=1, col=4);
            axis(1, at=c(1:length(date_labels)) * 30 - 15, cex.axis=3, labels=date_labels);
            axis(2, cex.axis=3);
            if (plot_std_error=="yes") {
                # Standard error for group.
                lines(days, group+group_std_error, lty=2);
                lines(days, group-group_std_error, lty=2);
                # Standard error for group2.
                lines(days, group2+group2_std_error, col=2, lty=2);
                lines(days, group2-group2_std_error, col=2, lty=2);
                # Standard error for group3.
                lines(days, group3+group3_std_error, col=4, lty=2);
                lines(days, group3-group3_std_error, col=4, lty=2);
            }
        } else {
            if (life_stage=="Egg") {
                title = paste(insect,  ": Reps", replications, ":", life_stage, "Pop :", location, ": Lat", latitude, ":", start_date, "-", end_date, sep=" ");
                legend_text = c(life_stage);
                columns = c(4);
            } else if (life_stage=="Nymph") {
                stage = paste(life_stages_nymph, "Nymph Pop :", sep=" ");
                title = paste(insect, ": Reps", replications, ":", stage, location, ": Lat", latitude, ":", start_date, "-", end_date, sep=" ");
                legend_text = c(paste(life_stages_nymph, life_stage, sep=" "));
                columns = c(2);
            } else if (life_stage=="Adult") {
                stage = paste(life_stages_adult, "Adult Pop", sep=" ");
                title = paste(insect, ": Reps", replications, ":", stage, location, ": Lat", latitude, ":", start_date, "-", end_date, sep=" ");
                legend_text = c(paste(life_stages_adult, life_stage, sep=" "));
                columns = c(1);
            }
            plot(days, group, main=title, type="l", ylim=c(0, maxval), axes=F, lwd=2, xlab="", ylab="", cex=3, cex.lab=3, cex.axis=3, cex.main=3);
            legend("topleft", legend_text, lty=c(1), col="black", cex=3);
            axis(1, at=c(1:length(date_labels)) * 30 - 15, cex.axis=3, labels=date_labels);
            axis(2, cex.axis=3);
            if (plot_std_error=="yes") {
                # Standard error for group.
                lines(days, group+group_std_error, lty=2);
                lines(days, group-group_std_error, lty=2);
            }
        }
    } else if (chart_type=="pop_size_by_generation") {
        if (life_stage=="Total") {
            title_str = ": Total Pop by Gen :";
        } else if (life_stage=="Egg") {
            title_str = ": Egg Pop by Gen :";
        } else if (life_stage=="Nymph") {
            title_str = paste(":", life_stages_nymph, "Nymph Pop by Gen", ":", sep=" ");
        } else if (life_stage=="Adult") {
            title_str = paste(":", life_stages_adult, "Adult Pop by Gen", ":", sep=" ");
        }
        title = paste(insect, ": Reps", replications, title_str, location, ": Lat", latitude, ":", start_date, "-", end_date, sep=" ");
        legend_text = c("P", "F1", "F2");
        columns = c(1, 2, 4);
        plot(days, group, main=title, type="l", ylim=c(0, maxval), axes=F, lwd=2, xlab="", ylab="", cex=3, cex.lab=3, cex.axis=3, cex.main=3);
        legend("topleft", legend_text, lty=c(1, 1, 1), col=columns, cex=3);
        lines(days, group2, lwd=2, lty=1, col=2);
        lines(days, group3, lwd=2, lty=1, col=4);
        axis(1, at=c(1:length(date_labels)) * 30 - 15, cex.axis=3, labels=date_labels);
        axis(2, cex.axis=3);
        if (plot_std_error=="yes") {
            # Standard error for group.
            lines(days, group+group_std_error, lty=2);
            lines(days, group-group_std_error, lty=2);
            # Standard error for group2.
            lines(days, group2+group2_std_error, col=2, lty=2);
            lines(days, group2-group2_std_error, col=2, lty=2);
            # Standard error for group3.
            lines(days, group3+group3_std_error, col=4, lty=2);
            lines(days, group3-group3_std_error, col=4, lty=2);
        }
    }
}

# Determine if we're plotting generations separately.
if (opt$plot_generations_separately=="yes") {
    plot_generations_separately = TRUE;
} else {
    plot_generations_separately = FALSE;
}
# Read the temperature data into a data frame.
temperature_data_frame = parse_input_data(opt$input, opt$num_days);
# Get the date labels for plots.
date_labels = get_date_labels(temperature_data_frame, opt$num_days);
# All latitude values are the same, so get the value for plots from the first row.
latitude = temperature_data_frame$LATITUDE[1];
# Determine the specified life stages for processing.
# Split life_stages into a list of strings for plots.
life_stages_str = as.character(opt$life_stages);
life_stages = strsplit(life_stages_str, ",")[[1]];
# Determine the data we need to generate for plotting.
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
for (life_stage in life_stages) {
    if (life_stage=="Total") {
        process_eggs = TRUE;
        process_nymphs = TRUE;
        process_adults = TRUE;
    } else if (life_stage=="Egg") {
        process_eggs = TRUE;
    } else if (life_stage=="Nymph") {
        process_nymphs = TRUE;
    } else if (life_stage=="Adult") {
        process_adults = TRUE;
    }
}
if (process_nymphs) {
    # Split life_stages_nymph into a list of strings for plots.
    life_stages_nymph_str = as.character(opt$life_stages_nymph);
    life_stages_nymph = strsplit(life_stages_nymph_str, ",")[[1]];
    for (life_stage_nymph in life_stages_nymph) {
        if (life_stage_nymph=="Young") {
            process_young_nymphs = TRUE;
        } else if (life_stage_nymph=="Old") {
            process_old_nymphs = TRUE;
        } else if (life_stage_nymph=="Total") {
            process_total_nymphs = TRUE;
        }
    }
}
if (process_adults) {
    # Split life_stages_adult into a list of strings for plots.
    life_stages_adult_str = as.character(opt$life_stages_adult);
    life_stages_adult = strsplit(life_stages_adult_str, ",")[[1]];
    for (life_stage_adult in life_stages_adult) {
        if (life_stage_adult=="Pre-vittelogenic") {
            process_previttelogenic_adults = TRUE;
        } else if (life_stage_adult=="Vittelogenic") {
            process_vittelogenic_adults = TRUE;
        } else if (life_stage_adult=="Diapausing") {
            process_diapausing_adults = TRUE;
        } else if (life_stage_adult=="Total") {
            process_total_adults = TRUE;
        }
    }
}
# Initialize matrices.
if (process_eggs) {
    Eggs.replications = matrix(rep(0, opt$num_days*opt$replications), ncol=opt$replications);
}
if (process_young_nymphs | process_total_nymphs) {
    YoungNymphs.replications = matrix(rep(0, opt$num_days*opt$replications), ncol=opt$replications);
}
if (process_old_nymphs | process_total_nymphs) {
    OldNymphs.replications = matrix(rep(0, opt$num_days*opt$replications), ncol=opt$replications);
}
if (process_previttelogenic_adults | process_total_adults) {
    Previttelogenic.replications = matrix(rep(0, opt$num_days*opt$replications), ncol=opt$replications);
}
if (process_vittelogenic_adults | process_total_adults) {
    Vittelogenic.replications = matrix(rep(0, opt$num_days*opt$replications), ncol=opt$replications);
}
if (process_diapausing_adults | process_total_adults) {
    Diapausing.replications = matrix(rep(0, opt$num_days*opt$replications), ncol=opt$replications);
}
newborn.replications = matrix(rep(0, opt$num_days*opt$replications), ncol=opt$replications);
adult.replications = matrix(rep(0, opt$num_days*opt$replications), ncol=opt$replications);
death.replications = matrix(rep(0, opt$num_days*opt$replications), ncol=opt$replications);
if (plot_generations_separately) {
    # P is Parental, or overwintered adults.
    P.replications = matrix(rep(0, opt$num_days*opt$replications), ncol=opt$replications);
    # F1 is the first field-produced generation.
    F1.replications = matrix(rep(0, opt$num_days*opt$replications), ncol=opt$replications);
    # F2 is the second field-produced generation.
    F2.replications = matrix(rep(0, opt$num_days*opt$replications), ncol=opt$replications);
    if (process_eggs) {
        P_eggs.replications = matrix(rep(0, opt$num_days*opt$replications), ncol=opt$replications);
        F1_eggs.replications = matrix(rep(0, opt$num_days*opt$replications), ncol=opt$replications);
        F2_eggs.replications = matrix(rep(0, opt$num_days*opt$replications), ncol=opt$replications);
    }
    if (process_young_nymphs) {
        P_young_nymphs.replications = matrix(rep(0, opt$num_days*opt$replications), ncol=opt$replications);
        F1_young_nymphs.replications = matrix(rep(0, opt$num_days*opt$replications), ncol=opt$replications);
        F2_young_nymphs.replications = matrix(rep(0, opt$num_days*opt$replications), ncol=opt$replications);
    }
    if (process_old_nymphs) {
        P_old_nymphs.replications = matrix(rep(0, opt$num_days*opt$replications), ncol=opt$replications);
        F1_old_nymphs.replications = matrix(rep(0, opt$num_days*opt$replications), ncol=opt$replications);
        F2_old_nymphs.replications = matrix(rep(0, opt$num_days*opt$replications), ncol=opt$replications);
    }
    if (process_total_nymphs) {
        P_total_nymphs.replications = matrix(rep(0, opt$num_days*opt$replications), ncol=opt$replications);
        F1_total_nymphs.replications = matrix(rep(0, opt$num_days*opt$replications), ncol=opt$replications);
        F2_total_nymphs.replications = matrix(rep(0, opt$num_days*opt$replications), ncol=opt$replications);
    }
    if (process_previttelogenic_adults) {
        P_previttelogenic_adults.replications = matrix(rep(0, opt$num_days*opt$replications), ncol=opt$replications);
        F1_previttelogenic_adults.replications = matrix(rep(0, opt$num_days*opt$replications), ncol=opt$replications);
        F2_previttelogenic_adults.replications = matrix(rep(0, opt$num_days*opt$replications), ncol=opt$replications);
    }
    if (process_vittelogenic_adults) {
        P_vittelogenic_adults.replications = matrix(rep(0, opt$num_days*opt$replications), ncol=opt$replications);
        F1_vittelogenic_adults.replications = matrix(rep(0, opt$num_days*opt$replications), ncol=opt$replications);
        F2_vittelogenic_adults.replications = matrix(rep(0, opt$num_days*opt$replications), ncol=opt$replications);
    }
    if (process_diapausing_adults) {
        P_diapausing_adults.replications = matrix(rep(0, opt$num_days*opt$replications), ncol=opt$replications);
        F1_diapausing_adults.replications = matrix(rep(0, opt$num_days*opt$replications), ncol=opt$replications);
        F2_diapausing_adults.replications = matrix(rep(0, opt$num_days*opt$replications), ncol=opt$replications);
    }
    if (process_total_adults) {
        P_total_adults.replications = matrix(rep(0, opt$num_days*opt$replications), ncol=opt$replications);
        F1_total_adults.replications = matrix(rep(0, opt$num_days*opt$replications), ncol=opt$replications);
        F2_total_adults.replications = matrix(rep(0, opt$num_days*opt$replications), ncol=opt$replications);
    }
}
# Total population.
population.replications = matrix(rep(0, opt$num_days*opt$replications), ncol=opt$replications);

# Process replications.
for (current_replication in 1:opt$replications) {
    # Start with the user-defined number of insects per replication.
    num_insects = opt$insects_per_replication;
    # Generation, Stage, degree-days, T, Diapause.
    vector.ini = c(0, 3, 0, 0, 0);
    # Replicate to create a matrix where the columns are
    # Generation, Stage, degree-days, T, Diapause and the
    # rows are the initial number of insects per replication.
    vector.matrix = rep(vector.ini, num_insects);
    # Complete transposed matrix for the population, so now
    # the rows are Generation, Stage, degree-days, T, Diapause
    vector.matrix = base::t(matrix(vector.matrix, nrow=5));
    # Time series of population size.
    if (process_eggs) {
        Eggs = rep(0, opt$num_days);
    }
    if (process_young_nymphs | process_total_nymphs) {
        YoungNymphs = rep(0, opt$num_days);
    }
    if (process_old_nymphs | process_total_nymphs) {
        OldNymphs = rep(0, opt$num_days);
    }
    if (process_previttelogenic_adults | process_total_adults) {
        Previttelogenic = rep(0, opt$num_days);
    }
    if (process_vittelogenic_adults | process_total_adults) {
        Vittelogenic = rep(0, opt$num_days);
    }
    if (process_diapausing_adults | process_total_adults) {
        Diapausing = rep(0, opt$num_days);
    }
    N.newborn = rep(0, opt$num_days);
    N.adult = rep(0, opt$num_days);
    N.death = rep(0, opt$num_days);
    overwintering_adult.population = rep(0, opt$num_days);
    first_generation.population = rep(0, opt$num_days);
    second_generation.population = rep(0, opt$num_days);
    if (plot_generations_separately) {
        # P is Parental, or overwintered adults.
        # F1 is the first field-produced generation.
        # F2 is the second field-produced generation.
        if (process_eggs) {
            P.egg = rep(0, opt$num_days);
            F1.egg = rep(0, opt$num_days);
            F2.egg = rep(0, opt$num_days);
        }
        if (process_young_nymphs) {
            P.young_nymph = rep(0, opt$num_days);
            F1.young_nymph = rep(0, opt$num_days);
            F2.young_nymph = rep(0, opt$num_days);
        }
        if (process_old_nymphs) {
            P.old_nymph = rep(0, opt$num_days);
            F1.old_nymph = rep(0, opt$num_days);
            F2.old_nymph = rep(0, opt$num_days);
        }
        if (process_total_nymphs) {
            P.total_nymph = rep(0, opt$num_days);
            F1.total_nymph = rep(0, opt$num_days);
            F2.total_nymph = rep(0, opt$num_days);
        }
        if (process_previttelogenic_adults) {
            P.previttelogenic_adult = rep(0, opt$num_days);
            F1.previttelogenic_adult = rep(0, opt$num_days);
            F2.previttelogenic_adult = rep(0, opt$num_days);
        }
        if (process_vittelogenic_adults) {
            P.vittelogenic_adult = rep(0, opt$num_days);
            F1.vittelogenic_adult = rep(0, opt$num_days);
            F2.vittelogenic_adult = rep(0, opt$num_days);
        }
        if (process_diapausing_adults) {
            P.diapausing_adult = rep(0, opt$num_days);
            F1.diapausing_adult = rep(0, opt$num_days);
            F2.diapausing_adult = rep(0, opt$num_days);
        }
        if (process_total_adults) {
            P.total_adult = rep(0, opt$num_days);
            F1.total_adult = rep(0, opt$num_days);
            F2.total_adult = rep(0, opt$num_days);
        }
    }
    total.population = NULL;
    averages.day = rep(0, opt$num_days);
    # All the days included in the input temperature dataset.
    for (row in 1:opt$num_days) {
        # Get the integer day of the year for the current row.
        doy = temperature_data_frame$DOY[row];
        # Photoperiod in the day.
        photoperiod = temperature_data_frame$DAYLEN[row];
        temp.profile = get_temperature_at_hour(latitude, temperature_data_frame, row, opt$num_days);
        mean.temp = temp.profile[1];
        averages.temp = temp.profile[2];
        averages.day[row] = averages.temp;
        # Trash bin for death.
        death.vector = NULL;
        # Newborn.
        birth.vector = NULL;
        # All individuals.
        for (i in 1:num_insects) {
            # Individual record.
            vector.individual = vector.matrix[i,];
            # Adjustment for late season mortality rate (still alive?).
            if (latitude < 40.0) {
                post.mortality = 1;
                day.kill = 300;
            }
            else {
                post.mortality = 2;
                day.kill = 250;
            }
            if (vector.individual[2] == 0) {
                # Egg.
                death.probability = opt$egg_mortality * mortality.egg(mean.temp);
            }
            else if (vector.individual[2] == 1 | vector.individual[2] == 2) {
                # Nymph.
                death.probability = opt$nymph_mortality * mortality.nymph(mean.temp);
            }
            else if (vector.individual[2] == 3 | vector.individual[2] == 4 | vector.individual[2] == 5) {
                # Adult.
                if (doy < day.kill) {
                    death.probability = opt$adult_mortality * mortality.adult(mean.temp);
                }
                else {
                    # Increase adult mortality after fall equinox.
                    death.probability = opt$adult_mortality * post.mortality * mortality.adult(mean.temp);
                }
            }
            # Dependent on temperature and life stage?
            u.d = runif(1);
            if (u.d < death.probability) {
                death.vector = c(death.vector, i);
            }
            else {
                # End of diapause.
                if (vector.individual[1] == 0 && vector.individual[2] == 3) {
                    # Overwintering adult (pre-vittelogenic).
                    if (photoperiod > opt$photoperiod && vector.individual[3] > 68 && doy < 180) {
                        # Add 68C to become fully reproductively matured.
                        # Transfer to vittelogenic.
                        vector.individual = c(0, 4, 0, 0, 0);
                        vector.matrix[i,] = vector.individual;
                    }
                    else {
                        # Add average temperature for current day.
                        vector.individual[3] = vector.individual[3] + averages.temp;
                        # Add 1 day in current stage.
                        vector.individual[4] = vector.individual[4] + 1;
                        vector.matrix[i,] = vector.individual;
                    }
                }
                if (vector.individual[1] != 0 && vector.individual[2] == 3) {
                    # Not overwintering adult (pre-vittelogenic).
                    current.gen = vector.individual[1];
                    if (vector.individual[3] > 68) {
                        # Add 68C to become fully reproductively matured.
                        # Transfer to vittelogenic.
                        vector.individual = c(current.gen, 4, 0, 0, 0);
                        vector.matrix[i,] = vector.individual;
                    }
                    else {
                        # Add average temperature for current day.
                        vector.individual[3] = vector.individual[3] + averages.temp;
                        # Add 1 day in current stage.
                        vector.individual[4] = vector.individual[4] + 1;
                        vector.matrix[i,] = vector.individual;
                    }
                }
                # Oviposition -- where population dynamics comes from.
                if (vector.individual[2] == 4 && vector.individual[1] == 0 && mean.temp > 10) {
                    # Vittelogenic stage, overwintering generation.
                    if (vector.individual[4] == 0) {
                        # Just turned in vittelogenic stage.
                        num_insects.birth = round(runif(1, 2 + opt$min_clutch_size, 8 + opt$max_clutch_size));
                    }
                    else {
                        # Daily probability of birth.
                        p.birth = opt$oviposition * 0.01;
                        u1 = runif(1);
                        if (u1 < p.birth) {
                            num_insects.birth = round(runif(1, 2, 8));
                        }
                    }
                    # Add average temperature for current day.
                    vector.individual[3] = vector.individual[3] + averages.temp;
                    # Add 1 day in current stage.
                    vector.individual[4] = vector.individual[4] + 1;
                    vector.matrix[i,] = vector.individual;
                    if (num_insects.birth > 0) {
                        # Add new birth -- might be in different generations.
                        new.gen = vector.individual[1] + 1;
                        # Egg profile.
                        new.individual = c(new.gen, 0, 0, 0, 0);
                        new.vector = rep(new.individual, num_insects.birth);
                        # Update batch of egg profile.
                        new.vector = t(matrix(new.vector, nrow=5));
                        # Group with total eggs laid in that day.
                        birth.vector = rbind(birth.vector, new.vector);
                    }
                }
                # Oviposition -- for generation 1.
                if (vector.individual[2] == 4 && vector.individual[1] == 1 && mean.temp > 12.5 && doy < 222) {
                    # Vittelogenic stage, 1st generation
                    if (vector.individual[4] == 0) {
                        # Just turned in vittelogenic stage.
                        num_insects.birth = round(runif(1, 2+opt$min_clutch_size, 8+opt$max_clutch_size));
                    }
                    else {
                        # Daily probability of birth.
                        p.birth = opt$oviposition * 0.01;
                        u1 = runif(1);
                        if (u1 < p.birth) {
                            num_insects.birth = round(runif(1, 2, 8));
                        }
                    }
                    # Add average temperature for current day.
                    vector.individual[3] = vector.individual[3] + averages.temp;
                    # Add 1 day in current stage.
                    vector.individual[4] = vector.individual[4] + 1;
                    vector.matrix[i,] = vector.individual;
                    if (num_insects.birth > 0) {
                        # Add new birth -- might be in different generations.
                        new.gen = vector.individual[1] + 1;
                        # Egg profile.
                        new.individual = c(new.gen, 0, 0, 0, 0);
                        new.vector = rep(new.individual, num_insects.birth);
                        # Update batch of egg profile.
                        new.vector = t(matrix(new.vector, nrow=5));
                        # Group with total eggs laid in that day.
                        birth.vector = rbind(birth.vector, new.vector);
                    }
                }
                # Egg to young nymph.
                if (vector.individual[2] == 0) {
                    # Add average temperature for current day.
                    vector.individual[3] = vector.individual[3] + averages.temp;
                    if (vector.individual[3] >= (68+opt$young_nymph_accumulation)) {
                        # From egg to young nymph, degree-days requirement met.
                        current.gen = vector.individual[1];
                        # Transfer to young nymph stage.
                        vector.individual = c(current.gen, 1, 0, 0, 0);
                    }
                    else {
                        # Add 1 day in current stage.
                        vector.individual[4] = vector.individual[4] + 1;
                    }
                    vector.matrix[i,] = vector.individual;
                }
                # Young nymph to old nymph.
                if (vector.individual[2] == 1) {
                    # Add average temperature for current day.
                    vector.individual[3] = vector.individual[3] + averages.temp;
                    if (vector.individual[3] >= (250+opt$old_nymph_accumulation)) {
                        # From young to old nymph, degree_days requirement met.
                        current.gen = vector.individual[1];
                        # Transfer to old nym stage.
                        vector.individual = c(current.gen, 2, 0, 0, 0);
                        if (photoperiod < opt$photoperiod && doy > 180) {
                            vector.individual[5] = 1;
                        } # Prepare for diapausing.
                    }
                    else {
                        # Add 1 day in current stage.
                        vector.individual[4] = vector.individual[4] + 1;
                    }
                    vector.matrix[i,] = vector.individual;
                }
                # Old nymph to adult: pre-vittelogenic or diapausing?
                if (vector.individual[2] == 2) {
                    # Add average temperature for current day.
                    vector.individual[3] = vector.individual[3] + averages.temp;
                    if (vector.individual[3] >= (200+opt$adult_accumulation)) {
                        # From old to adult, degree_days requirement met.
                        current.gen = vector.individual[1];
                        if (vector.individual[5] == 0) {
                            # Previttelogenic.
                            vector.individual = c(current.gen, 3, 0, 0, 0);
                        }
                        else {
                            # Diapausing.
                            vector.individual = c(current.gen, 5, 0, 0, 1);
                        }
                    }
                    else {
                        # Add 1 day in current stage.
                        vector.individual[4] = vector.individual[4] + 1;
                    }
                    vector.matrix[i,] = vector.individual;
                }
                # Growing of diapausing adult (unimportant, but still necessary).
                if (vector.individual[2] == 5) {
                    vector.individual[3] = vector.individual[3] + averages.temp;
                    vector.individual[4] = vector.individual[4] + 1;
                    vector.matrix[i,] = vector.individual;
                }
            } # Else if it is still alive.
        } # End of the individual bug loop.

        # Number of deaths.
        num_insects.death = length(death.vector);
        if (num_insects.death > 0) {
            # Remove record of dead.
            vector.matrix = vector.matrix[-death.vector,];
        }
        # Number of births.
        num_insects.newborn = length(birth.vector[,1]);
        vector.matrix = rbind(vector.matrix, birth.vector);
        # Update population size for the next day.
        num_insects = num_insects - num_insects.death + num_insects.newborn;

        # Aggregate results by day.  Due to multiple transpose calls
        # on vector.matrix above, the columns of vector.matrix
        # are now Generation, Stage, degree-days, T, Diapause,
        if (process_eggs) {
            # For egg population size, column 2 (Stage), must be 0.
            Eggs[row] = sum(vector.matrix[,2]==0);
        }
        if (process_young_nymphs | process_total_nymphs) {
            # For young nymph population size, column 2 (Stage) must be 1.
            YoungNymphs[row] = sum(vector.matrix[,2]==1);
        }
        if (process_old_nymphs | process_total_nymphs) {
            # For old nymph population size, column 2 (Stage) must be 2.
            OldNymphs[row] = sum(vector.matrix[,2]==2);
        }
        if (process_previttelogenic_adults | process_total_adults) {
            # For pre-vittelogenic population size, column 2 (Stage) must be 3.
            Previttelogenic[row] = sum(vector.matrix[,2]==3);
        }
        if (process_vittelogenic_adults | process_total_adults) {
            # For vittelogenic population size, column 2 (Stage) must be 4.
            Vittelogenic[row] = sum(vector.matrix[,2]==4);
        }
        if (process_diapausing_adults | process_total_adults) {
            # For diapausing population size, column 2 (Stage) must be 5.
            Diapausing[row] = sum(vector.matrix[,2]==5);
        }

        # Newborn population size.
        N.newborn[row] = num_insects.newborn;
        # Adult population size.
        N.adult[row] = sum(vector.matrix[,2]==3) + sum(vector.matrix[,2]==4) + sum(vector.matrix[,2]==5);
        # Dead population size.
        N.death[row] = num_insects.death;

        total.population = c(total.population, num_insects);

        # For overwintering adult (P) population
        # size, column 1 (Generation) must be 0.
        overwintering_adult.population[row] = sum(vector.matrix[,1]==0);
        # For first field generation (F1) population
        # size, column 1 (Generation) must be 1.
        first_generation.population[row] = sum(vector.matrix[,1]==1);
        # For second field generation (F2) population
        # size, column 1 (Generation) must be 2.
        second_generation.population[row] = sum(vector.matrix[,1]==2);

        if (plot_generations_separately) {
            if (process_eggs) {
                # For egg life stage of generation P population size,
                # column 1 (generation) is 0 and column 2 (Stage) is 0.
                P.egg[row] = sum(vector.matrix[,1]==0 & vector.matrix[,2]==0);
                # For egg life stage of generation F1 population size,
                # column 1 (generation) is 1 and column 2 (Stage) is 0.
                F1.egg[row] = sum(vector.matrix[,1]==1 & vector.matrix[,2]==0);
                # For egg life stage of generation F2 population size,
                # column 1 (generation) is 2 and column 2 (Stage) is 0.
                F2.egg[row] = sum(vector.matrix[,1]==2 & vector.matrix[,2]==0);
            }
            if (process_young_nymphs) {
                # For young nymph life stage of generation P population
                # size, the following combination is required:
                # - column 1 (Generation) is 0 and column 2 (Stage) is 1 (Young nymph)
                P.young_nymph[row] = sum(vector.matrix[,1]==0 & vector.matrix[,2]==1);
                # For young nymph life stage of generation F1 population
                # size, the following combination is required:
                # - column 1 (Generation) is 1 and column 2 (Stage) is 1 (Young nymph)
                F1.young_nymph[row] = sum(vector.matrix[,1]==1 & vector.matrix[,2]==1);
                # For young nymph life stage of generation F2 population
                # size, the following combination is required:
                # - column 1 (Generation) is 2 and column 2 (Stage) is 1 (Young nymph)
                F2.young_nymph[row] = sum(vector.matrix[,1]==2 & vector.matrix[,2]==1);
            }
            if (process_old_nymphs) {
                # For old nymph life stage of generation P population
                # size, the following combination is required:
                # - column 1 (Generation) is 0 and column 2 (Stage) is 2 (Old nymph)
                P.old_nymph[row] = sum(vector.matrix[,1]==0 & vector.matrix[,2]==2);
                # For old nymph life stage of generation F1 population
                # size, the following combination is required:
                # - column 1 (Generation) is 1 and column 2 (Stage) is 2 (Old nymph)
                F1.old_nymph[row] = sum(vector.matrix[,1]==1 & vector.matrix[,2]==2);
                # For old nymph life stage of generation F2 population
                # size, the following combination is required:
                # - column 1 (Generation) is 2 and column 2 (Stage) is 2 (Old nymph)
                F2.old_nymph[row] = sum(vector.matrix[,1]==2 & vector.matrix[,2]==2);
            }
            if (process_total_nymphs) {
                # For total nymph life stage of generation P population
                # size, one of the following combinations is required:
                # - column 1 (Generation) is 0 and column 2 (Stage) is 1 (Young nymph)
                # - column 1 (Generation) is 0 and column 2 (Stage) is 2 (Old nymph)
                P.total_nymph[row] = sum((vector.matrix[,1]==0 & vector.matrix[,2]==1) | (vector.matrix[,1]==0 & vector.matrix[,2]==2));
                # For total nymph life stage of generation F1 population
                # size, one of the following combinations is required:
                # - column 1 (Generation) is 1 and column 2 (Stage) is 1 (Young nymph)
                # - column 1 (Generation) is 1 and column 2 (Stage) is 2 (Old nymph)
                F1.total_nymph[row] = sum((vector.matrix[,1]==1 & vector.matrix[,2]==1) | (vector.matrix[,1]==1 & vector.matrix[,2]==2));
                # For total nymph life stage of generation F2 population
                # size, one of the following combinations is required:
                # - column 1 (Generation) is 2 and column 2 (Stage) is 1 (Young nymph)
                # - column 1 (Generation) is 2 and column 2 (Stage) is 2 (Old nymph)
                F2.total_nymph[row] = sum((vector.matrix[,1]==2 & vector.matrix[,2]==1) | (vector.matrix[,1]==2 & vector.matrix[,2]==2));
            }
            if (process_previttelogenic_adults) {
                # For previttelogenic adult life stage of generation P population
                # size, the following combination is required:
                # - column 1 (Generation) is 0 and column 2 (Stage) is 3 (Pre-vittelogenic)
                P.previttelogenic_adult[row] = sum(vector.matrix[,1]==0 & vector.matrix[,2]==3);
                # For previttelogenic adult life stage of generation F1 population
                # size, the following combination is required:
                # - column 1 (Generation) is 1 and column 2 (Stage) is 3 (Pre-vittelogenic)
                F1.previttelogenic_adult[row] = sum(vector.matrix[,1]==1 & vector.matrix[,2]==3);
                # For previttelogenic adult life stage of generation F2 population
                # size, the following combination is required:
                # - column 1 (Generation) is 2 and column 2 (Stage) is 3 (Pre-vittelogenic)
                F2.previttelogenic_adult[row] = sum(vector.matrix[,1]==2 & vector.matrix[,2]==3);
            }
            if (process_vittelogenic_adults) {
                # For vittelogenic adult life stage of generation P population
                # size, the following combination is required:
                # - column 1 (Generation) is 0 and column 2 (Stage) is 4 (Vittelogenic)
                P.vittelogenic_adult[row] = sum(vector.matrix[,1]==0 & vector.matrix[,2]==4);
                # For vittelogenic adult life stage of generation F1 population
                # size, the following combination is required:
                # - column 1 (Generation) is 1 and column 2 (Stage) is 4 (Vittelogenic)
                F1.vittelogenic_adult[row] = sum(vector.matrix[,1]==1 & vector.matrix[,2]==4);
                # For vittelogenic adult life stage of generation F2 population
                # size, the following combination is required:
                # - column 1 (Generation) is 2 and column 2 (Stage) is 4 (Vittelogenic)
                F2.vittelogenic_adult[row] = sum(vector.matrix[,1]==2 & vector.matrix[,2]==4);
            }
            if (process_diapausing_adults) {
                # For diapausing adult life stage of generation P population
                # size, the following combination is required:
                # - column 1 (Generation) is 0 and column 2 (Stage) is 5 (Diapausing)
                P.diapausing_adult[row] = sum(vector.matrix[,1]==0 & vector.matrix[,2]==5);
                # For diapausing adult life stage of generation F1 population
                # size, the following combination is required:
                # - column 1 (Generation) is 1 and column 2 (Stage) is 5 (Diapausing)
                F1.diapausing_adult[row] = sum(vector.matrix[,1]==1 & vector.matrix[,2]==5);
                # For diapausing adult life stage of generation F2 population
                # size, the following combination is required:
                # - column 1 (Generation) is 2 and column 2 (Stage) is 5 (Diapausing)
                F2.diapausing_adult[row] = sum(vector.matrix[,1]==2 & vector.matrix[,2]==5);
            }
            if (process_total_adults) {
                # For total adult life stage of generation P population
                # size, one of the following combinations is required:
                # - column 1 (Generation) is 0 and column 2 (Stage) is 3 (Pre-vittelogenic)
                # - column 1 (Generation) is 0 and column 2 (Stage) is 4 (Vittelogenic)
                # - column 1 (Generation) is 0 and column 2 (Stage) is 5 (Diapausing)
                P.total_adult[row] = sum((vector.matrix[,1]==0 & vector.matrix[,2]==3) | (vector.matrix[,1]==0 & vector.matrix[,2]==4) | (vector.matrix[,1]==0 & vector.matrix[,2]==5));
                # For total adult life stage of generation F1 population
                # size, one of the following combinations is required:
                # - column 1 (Generation) is 1 and column 2 (Stage) is 3 (Pre-vittelogenic)
                # - column 1 (Generation) is 1 and column 2 (Stage) is 4 (Vittelogenic)
                # - column 1 (Generation) is 1 and column 2 (Stage) is 5 (Diapausing)
                F1.total_adult[row] = sum((vector.matrix[,1]==1 & vector.matrix[,2]==3) | (vector.matrix[,1]==1 & vector.matrix[,2]==4) | (vector.matrix[,1]==1 & vector.matrix[,2]==5));
                # For total adult life stage of generation F2 population
                # size, one of the following combinations is required:
                # - column 1 (Generation) is 2 and column 2 (Stage) is 3 (Pre-vittelogenic)
                # - column 1 (Generation) is 2 and column 2 (Stage) is 4 (Vittelogenic)
                # - column 1 (Generation) is 2 and column 2 (Stage) is 5 (Diapausing)
                F2.total_adult[row] = sum((vector.matrix[,1]==2 & vector.matrix[,2]==3) | (vector.matrix[,1]==2 & vector.matrix[,2]==4) | (vector.matrix[,1]==2 & vector.matrix[,2]==5));
            }
        }
    }   # End of days specified in the input temperature data.

    averages.cum = cumsum(averages.day);

    # Define the output values.
    if (process_eggs) {
        Eggs.replications[,current_replication] = Eggs;
    }
    if (process_young_nymphs | process_total_nymphs) {
        YoungNymphs.replications[,current_replication] = YoungNymphs;
    }
    if (process_old_nymphs | process_total_nymphs) {
        OldNymphs.replications[,current_replication] = OldNymphs;
    }
    if (process_previttelogenic_adults | process_total_adults) {
        Previttelogenic.replications[,current_replication] = Previttelogenic;
    }
    if (process_vittelogenic_adults | process_total_adults) {
        Vittelogenic.replications[,current_replication] = Vittelogenic;
    }
    if (process_diapausing_adults | process_total_adults) {
        Diapausing.replications[,current_replication] = Diapausing;
    }
    newborn.replications[,current_replication] = N.newborn;
    adult.replications[,current_replication] = N.adult;
    death.replications[,current_replication] = N.death;
    if (plot_generations_separately) {
        # P is Parental, or overwintered adults.
        P.replications[,current_replication] = overwintering_adult.population;
        # F1 is the first field-produced generation.
        F1.replications[,current_replication] = first_generation.population;
        # F2 is the second field-produced generation.
        F2.replications[,current_replication] = second_generation.population;
        if (process_eggs) {
            P_eggs.replications[,current_replication] = P.egg;
            F1_eggs.replications[,current_replication] = F1.egg;
            F2_eggs.replications[,current_replication] = F2.egg;
        }
        if (process_young_nymphs) {
            P_young_nymphs.replications[,current_replication] = P.young_nymph;
            F1_young_nymphs.replications[,current_replication] = F1.young_nymph;
            F2_young_nymphs.replications[,current_replication] = F2.young_nymph;
        }
        if (process_old_nymphs) {
            P_old_nymphs.replications[,current_replication] = P.old_nymph;
            F1_old_nymphs.replications[,current_replication] = F1.old_nymph;
            F2_old_nymphs.replications[,current_replication] = F2.old_nymph;
        }
        if (process_total_nymphs) {
            P_total_nymphs.replications[,current_replication] = P.total_nymph;
            F1_total_nymphs.replications[,current_replication] = F1.total_nymph;
            F2_total_nymphs.replications[,current_replication] = F2.total_nymph;
        }
        if (process_previttelogenic_adults) {
            P_previttelogenic_adults.replications[,current_replication] = P.previttelogenic_adult;
            F1_previttelogenic_adults.replications[,current_replication] = F1.previttelogenic_adult;
            F2_previttelogenic_adults.replications[,current_replication] = F2.previttelogenic_adult;
        }
        if (process_vittelogenic_adults) {
            P_vittelogenic_adults.replications[,current_replication] = P.vittelogenic_adult;
            F1_vittelogenic_adults.replications[,current_replication] = F1.vittelogenic_adult;
            F2_vittelogenic_adults.replications[,current_replication] = F2.vittelogenic_adult;
        }
        if (process_diapausing_adults) {
            P_diapausing_adults.replications[,current_replication] = P.diapausing_adult;
            F1_diapausing_adults.replications[,current_replication] = F1.diapausing_adult;
            F2_diapausing_adults.replications[,current_replication] = F2.diapausing_adult;
        }
        if (process_total_adults) {
            P_total_adults.replications[,current_replication] = P.total_adult;
            F1_total_adults.replications[,current_replication] = F1.total_adult;
            F2_total_adults.replications[,current_replication] = F2.total_adult;
        }
    }
    population.replications[,current_replication] = total.population;
    # End processing replications.
}

if (process_eggs) {
    # Mean value for eggs.
    eggs = apply(Eggs.replications, 1, mean);
    temperature_data_frame = append_vector(temperature_data_frame, eggs, "EGG");
    # Standard error for eggs.
    eggs.std_error = apply(Eggs.replications, 1, sd) / sqrt(opt$replications);
    temperature_data_frame = append_vector(temperature_data_frame, eggs.std_error, "EGGSE");
}
if (process_nymphs) {
    # Calculate nymph populations for selected life stage.
    for (life_stage_nymph in life_stages_nymph) {
        if (life_stage_nymph=="Total") {
            # Mean value for all nymphs.
            total_nymphs = apply((YoungNymphs.replications+OldNymphs.replications), 1, mean);
            temperature_data_frame = append_vector(temperature_data_frame, total_nymphs, "TOTALNYMPH");
            # Standard error for all nymphs.
            total_nymphs.std_error = apply((YoungNymphs.replications+OldNymphs.replications) / sqrt(opt$replications), 1, sd);
            temperature_data_frame = append_vector(temperature_data_frame, total_nymphs.std_error, "TOTALNYMPHSE");
        } else if (life_stage_nymph=="Young") {
            # Mean value for young nymphs.
            young_nymphs = apply(YoungNymphs.replications, 1, mean);
            temperature_data_frame = append_vector(temperature_data_frame, young_nymphs, "YOUNGNYMPH");
            # Standard error for young nymphs.
            young_nymphs.std_error = apply(YoungNymphs.replications / sqrt(opt$replications), 1, sd);
            temperature_data_frame = append_vector(temperature_data_frame, young_nymphs.std_error, "YOUNGNYMPHSE");
        } else if (life_stage_nymph=="Old") {
            # Mean value for old nymphs.
            old_nymphs = apply(OldNymphs.replications, 1, mean);
            temperature_data_frame = append_vector(temperature_data_frame, old_nymphs, "OLDNYMPH");
            # Standard error for old nymphs.
            old_nymphs.std_error = apply(OldNymphs.replications / sqrt(opt$replications), 1, sd);
            temperature_data_frame = append_vector(temperature_data_frame, old_nymphs.std_error, "OLDNYMPHSE");
        }
    }
}
if (process_adults) {
    # Calculate adult populations for selected life stage.
    for (life_stage_adult in life_stages_adult) {
        if (life_stage_adult=="Total") {
            # Mean value for all adults.
            total_adults = apply((Previttelogenic.replications+Vittelogenic.replications+Diapausing.replications), 1, mean);
            temperature_data_frame = append_vector(temperature_data_frame, total_adults, "TOTALADULT");
            # Standard error for all adults.
            total_adults.std_error = apply((Previttelogenic.replications+Vittelogenic.replications+Diapausing.replications), 1, sd) / sqrt(opt$replications);
            temperature_data_frame = append_vector(temperature_data_frame, total_adults.std_error, "TOTALADULTSE");
        } else if (life_stage_adult == "Pre-vittelogenic") {
            # Mean value for previttelogenic adults.
            previttelogenic_adults = apply(Previttelogenic.replications, 1, mean);
            temperature_data_frame = append_vector(temperature_data_frame, previttelogenic_adults, "PRE-VITADULT");
            # Standard error for previttelogenic adults.
            previttelogenic_adults.std_error = apply(Previttelogenic.replications, 1, sd) / sqrt(opt$replications);
            temperature_data_frame = append_vector(temperature_data_frame, previttelogenic_adults.std_error, "PRE-VITADULTSE");
        } else if (life_stage_adult == "Vittelogenic") {
            # Mean value for vittelogenic adults.
            vittelogenic_adults = apply(Vittelogenic.replications, 1, mean);
            temperature_data_frame = append_vector(temperature_data_frame, vittelogenic_adults, "VITADULT");
            # Standard error for vittelogenic adults.
            vittelogenic_adults.std_error = apply(Vittelogenic.replications, 1, sd) / sqrt(opt$replications);
            temperature_data_frame = append_vector(temperature_data_frame, vittelogenic_adults.std_error, "VITADULTSE");
        } else if (life_stage_adult == "Diapausing") {
            # Mean value for vittelogenic adults.
            diapausing_adults = apply(Diapausing.replications, 1, mean);
            temperature_data_frame = append_vector(temperature_data_frame, diapausing_adults, "DIAPAUSINGADULT");
            # Standard error for vittelogenic adults.
            diapausing_adults.std_error = apply(Diapausing.replications, 1, sd) / sqrt(opt$replications);
            temperature_data_frame = append_vector(temperature_data_frame, diapausing_adults.std_error, "DIAPAUSINGADULTSE");
        }
    }
}

if (plot_generations_separately) {
    m_se = get_mean_and_std_error(P.replications, F1.replications, F2.replications);
    P = m_se[[1]];
    P.std_error = m_se[[2]];
    F1 = m_se[[3]];
    F1.std_error = m_se[[4]];
    F2 = m_se[[5]];
    F2.std_error = m_se[[6]];
    if (process_eggs) {
        m_se = get_mean_and_std_error(P_eggs.replications, F1_eggs.replications, F2_eggs.replications);
        P_eggs = m_se[[1]];
        P_eggs.std_error = m_se[[2]];
        F1_eggs = m_se[[3]];
        F1_eggs.std_error = m_se[[4]];
        F2_eggs = m_se[[5]];
        F2_eggs.std_error = m_se[[6]];
    }
    if (process_young_nymphs) {
        m_se = get_mean_and_std_error(P_young_nymphs.replications, F1_young_nymphs.replications, F2_young_nymphs.replications);
        P_young_nymphs = m_se[[1]];
        P_young_nymphs.std_error = m_se[[2]];
        F1_young_nymphs = m_se[[3]];
        F1_young_nymphs.std_error = m_se[[4]];
        F2_young_nymphs = m_se[[5]];
        F2_young_nymphs.std_error = m_se[[6]];
    }
    if (process_old_nymphs) {
        m_se = get_mean_and_std_error(P_old_nymphs.replications, F1_old_nymphs.replications, F2_old_nymphs.replications);
        P_old_nymphs = m_se[[1]];
        P_old_nymphs.std_error = m_se[[2]];
        F1_old_nymphs = m_se[[3]];
        F1_old_nymphs.std_error = m_se[[4]];
        F2_old_nymphs = m_se[[5]];
        F2_old_nymphs.std_error = m_se[[6]];
    }
    if (process_total_nymphs) {
        m_se = get_mean_and_std_error(P_total_nymphs.replications, F1_total_nymphs.replications, F2_total_nymphs.replications);
        P_total_nymphs = m_se[[1]];
        P_total_nymphs.std_error = m_se[[2]];
        F1_total_nymphs = m_se[[3]];
        F1_total_nymphs.std_error = m_se[[4]];
        F2_total_nymphs = m_se[[5]];
        F2_total_nymphs.std_error = m_se[[6]];
    }
    if (process_previttelogenic_adults) {
        m_se = get_mean_and_std_error(P_previttelogenic_adults.replications, F1_previttelogenic_adults.replications, F2_previttelogenic_adults.replications);
        P_previttelogenic_adults = m_se[[1]];
        P_previttelogenic_adults.std_error = m_se[[2]];
        F1_previttelogenic_adults = m_se[[3]];
        F1_previttelogenic_adults.std_error = m_se[[4]];
        F2_previttelogenic_adults = m_se[[5]];
        F2_previttelogenic_adults.std_error = m_se[[6]];
    }
    if (process_vittelogenic_adults) {
        m_se = get_mean_and_std_error(P_vittelogenic_adults.replications, F1_vittelogenic_adults.replications, F2_vittelogenic_adults.replications);
        P_vittelogenic_adults = m_se[[1]];
        P_vittelogenic_adults.std_error = m_se[[2]];
        F1_vittelogenic_adults = m_se[[3]];
        F1_vittelogenic_adults.std_error = m_se[[4]];
        F2_vittelogenic_adults = m_se[[5]];
        F2_vittelogenic_adults.std_error = m_se[[6]];
    }
    if (process_diapausing_adults) {
        m_se = get_mean_and_std_error(P_diapausing_adults.replications, F1_diapausing_adults.replications, F2_diapausing_adults.replications);
        P_diapausing_adults = m_se[[1]];
        P_diapausing_adults.std_error = m_se[[2]];
        F1_diapausing_adults = m_se[[3]];
        F1_diapausing_adults.std_error = m_se[[4]];
        F2_diapausing_adults = m_se[[5]];
        F2_diapausing_adults.std_error = m_se[[6]];
    }
    if (process_total_adults) {
        m_se = get_mean_and_std_error(P_total_adults.replications, F1_total_adults.replications, F2_total_adults.replications);
        P_total_adults = m_se[[1]];
        P_total_adults.std_error = m_se[[2]];
        F1_total_adults = m_se[[3]];
        F1_total_adults.std_error = m_se[[4]];
        F2_total_adults = m_se[[5]];
        F2_total_adults.std_error = m_se[[6]];
    }
}

# Save the analyzed data.
write.csv(temperature_data_frame, file=opt$output, row.names=F);
# Display the total number of days in the Galaxy history item blurb.
cat("Number of days: ", opt$num_days, "\n");
# Information needed for plots plots.
days = c(1:opt$num_days);
start_date = temperature_data_frame$DATE[1];
end_date = temperature_data_frame$DATE[opt$num_days];

if (plot_generations_separately) {
    for (life_stage in life_stages) {
        if (life_stage == "Egg") {
            # Start PDF device driver.
            dev.new(width=20, height=30);
            file_path = get_file_path(life_stage, "egg_pop_by_generation.pdf")
            pdf(file=file_path, width=20, height=30, bg="white");
            par(mar=c(5, 6, 4, 4), mfrow=c(3, 1));
            # Egg population size by generation.
            maxval = max(P_eggs+F1_eggs+F2_eggs) + 100;
            render_chart(date_labels, "pop_size_by_generation", opt$plot_std_error, opt$insect, opt$location, latitude, start_date, end_date, days, maxval,
                opt$replications, life_stage, group=P_eggs, group_std_error=P_eggs.std_error, group2=F1_eggs, group2_std_error=F1_eggs.std_error, group3=F2_eggs,
                group3_std_error=F2_eggs.std_error);
            # Turn off device driver to flush output.
            dev.off();
        } else if (life_stage == "Nymph") {
            for (life_stage_nymph in life_stages_nymph) {
                # Start PDF device driver.
                dev.new(width=20, height=30);
                file_path = get_file_path(life_stage, "nymph_pop_by_generation.pdf", life_stage_nymph=life_stage_nymph)
                pdf(file=file_path, width=20, height=30, bg="white");
                par(mar=c(5, 6, 4, 4), mfrow=c(3, 1));
                if (life_stage_nymph=="Young") {
                    # Young nymph population size by generation.
                    maxval = max(P_young_nymphs+F1_young_nymphs+F2_young_nymphs) + 100;
                    group = P_young_nymphs;
                    group_std_error = P_young_nymphs.std_error;
                    group2 = F1_young_nymphs;
                    group2_std_error = F1_young_nymphs.std_error;
                    group3 = F2_young_nymphs;
                    group3_std_error = F2_young_nymphs.std_error;
                } else if (life_stage_nymph=="Old") {
                    # Total nymph population size by generation.
                    maxval = max(P_old_nymphs+F1_old_nymphs+F2_old_nymphs) + 100;
                    group = P_old_nymphs;
                    group_std_error = P_old_nymphs.std_error;
                    group2 = F1_old_nymphs;
                    group2_std_error = F1_old_nymphs.std_error;
                    group3 = F2_old_nymphs;
                    group3_std_error = F2_old_nymphs.std_error;
                } else if (life_stage_nymph=="Total") {
                    # Total nymph population size by generation.
                    maxval = max(P_total_nymphs+F1_total_nymphs+F2_total_nymphs) + 100;
                    group = P_total_nymphs;
                    group_std_error = P_total_nymphs.std_error;
                    group2 = F1_total_nymphs;
                    group2_std_error = F1_total_nymphs.std_error;
                    group3 = F2_total_nymphs;
                    group3_std_error = F2_total_nymphs.std_error;
                }
                render_chart(date_labels, "pop_size_by_generation", opt$plot_std_error, opt$insect, opt$location, latitude, start_date, end_date, days, maxval,
                    opt$replications, life_stage, group=group, group_std_error=group_std_error, group2=group2, group2_std_error=group2_std_error,
                    group3=group3, group3_std_error=group3_std_error, life_stages_nymph=life_stage_nymph);
                # Turn off device driver to flush output.
                dev.off();
            }
        } else if (life_stage == "Adult") {
            for (life_stage_adult in life_stages_adult) {
                # Start PDF device driver.
                dev.new(width=20, height=30);
                file_path = get_file_path(life_stage, "adult_pop_by_generation.pdf", life_stage_adult=life_stage_adult)
                pdf(file=file_path, width=20, height=30, bg="white");
                par(mar=c(5, 6, 4, 4), mfrow=c(3, 1));
                if (life_stage_adult=="Pre-vittelogenic") {
                    # Pre-vittelogenic adult population size by generation.
                    maxval = max(P_previttelogenic_adults+F1_previttelogenic_adults+F2_previttelogenic_adults) + 100;
                    group = P_previttelogenic_adults;
                    group_std_error = P_previttelogenic_adults.std_error;
                    group2 = F1_previttelogenic_adults;
                    group2_std_error = F1_previttelogenic_adults.std_error;
                    group3 = F2_previttelogenic_adults;
                    group3_std_error = F2_previttelogenic_adults.std_error;
                } else if (life_stage_adult=="Vittelogenic") {
                    # Vittelogenic adult population size by generation.
                    maxval = max(P_vittelogenic_adults+F1_vittelogenic_adults+F2_vittelogenic_adults) + 100;
                    group = P_vittelogenic_adults;
                    group_std_error = P_vittelogenic_adults.std_error;
                    group2 = F1_vittelogenic_adults;
                    group2_std_error = F1_vittelogenic_adults.std_error;
                    group3 = F2_vittelogenic_adults;
                    group3_std_error = F2_vittelogenic_adults.std_error;
                } else if (life_stage_adult=="Diapausing") {
                    # Diapausing adult population size by generation.
                    maxval = max(P_diapausing_adults+F1_diapausing_adults+F2_diapausing_adults) + 100;
                    group = P_diapausing_adults;
                    group_std_error = P_diapausing_adults.std_error;
                    group2 = F1_diapausing_adults;
                    group2_std_error = F1_diapausing_adults.std_error;
                    group3 = F2_diapausing_adults;
                    group3_std_error = F2_diapausing_adults.std_error;
                } else if (life_stage_adult=="Total") {
                    # Total adult population size by generation.
                    maxval = max(P_total_adults+F1_total_adults+F2_total_adults) + 100;
                    group = P_total_adults;
                    group_std_error = P_total_adults.std_error;
                    group2 = F1_total_adults;
                    group2_std_error = F1_total_adults.std_error;
                    group3 = F2_total_adults;
                    group3_std_error = F2_total_adults.std_error;
                }
                render_chart(date_labels, "pop_size_by_generation", opt$plot_std_error, opt$insect, opt$location, latitude, start_date, end_date, days, maxval,
                    opt$replications, life_stage, group=group, group_std_error=group_std_error, group2=group2, group2_std_error=group2_std_error,
                    group3=group3, group3_std_error=group3_std_error, life_stages_adult=life_stage_adult);
                # Turn off device driver to flush output.
                dev.off();
            }
        } else if (life_stage == "Total") {
            # Start PDF device driver.
            # Name collection elements so that they
            # are displayed in logical order.
            dev.new(width=20, height=30);
            file_path = get_file_path(life_stage, "total_pop_by_generation.pdf")
            pdf(file=file_path, width=20, height=30, bg="white");
            par(mar=c(5, 6, 4, 4), mfrow=c(3, 1));
            # Total population size by generation.
            maxval = max(P+F1+F2) + 100;
            render_chart(date_labels, "pop_size_by_generation", opt$plot_std_error, opt$insect, opt$location, latitude, start_date, end_date, days, maxval,
                opt$replications, life_stage, group=P, group_std_error=P.std_error, group2=F1, group2_std_error=F1.std_error, group3=F2, group3_std_error=F2.std_error);
            # Turn off device driver to flush output.
            dev.off();
        }
    }
} else {
    for (life_stage in life_stages) {
        if (life_stage == "Egg") {
            # Start PDF device driver.
            dev.new(width=20, height=30);
            file_path = get_file_path(life_stage, "egg_pop.pdf")
            pdf(file=file_path, width=20, height=30, bg="white");
            par(mar=c(5, 6, 4, 4), mfrow=c(3, 1));
            # Egg population size.
            maxval = max(eggs+eggs.std_error) + 100;
            render_chart(date_labels, "pop_size_by_life_stage", opt$plot_std_error, opt$insect, opt$location, latitude, start_date, end_date, days, maxval,
                opt$replications, life_stage, group=eggs, group_std_error=eggs.std_error);
            # Turn off device driver to flush output.
            dev.off();
        } else if (life_stage == "Nymph") {
            for (life_stage_nymph in life_stages_nymph) {
                # Start PDF device driver.
                dev.new(width=20, height=30);
                file_path = get_file_path(life_stage, "nymph_pop.pdf", life_stage_nymph=life_stage_nymph)
                pdf(file=file_path, width=20, height=30, bg="white");
                par(mar=c(5, 6, 4, 4), mfrow=c(3, 1));
                if (life_stage_nymph=="Total") {
                    # Total nymph population size.
                    group = total_nymphs;
                    group_std_error = total_nymphs.std_error;
                } else if (life_stage_nymph=="Young") {
                    # Young nymph population size.
                    group = young_nymphs;
                    group_std_error = young_nymphs.std_error;
                } else if (life_stage_nymph=="Old") {
                    # Old nymph population size.
                    group = old_nymphs;
                    group_std_error = old_nymphs.std_error;
                }
                maxval = max(group+group_std_error) + 100;
                render_chart(date_labels, "pop_size_by_life_stage", opt$plot_std_error, opt$insect, opt$location, latitude, start_date, end_date, days, maxval,
                    opt$replications, life_stage, group=group, group_std_error=group_std_error, life_stages_nymph=life_stage_nymph);
                # Turn off device driver to flush output.
                dev.off();
            }
        } else if (life_stage == "Adult") {
            for (life_stage_adult in life_stages_adult) {
                # Start PDF device driver.
                dev.new(width=20, height=30);
                file_path = get_file_path(life_stage, "adult_pop.pdf", life_stage_adult=life_stage_adult)
                pdf(file=file_path, width=20, height=30, bg="white");
                par(mar=c(5, 6, 4, 4), mfrow=c(3, 1));
                if (life_stage_adult=="Total") {
                    # Total adult population size.
                    group = total_adults;
                    group_std_error = total_adults.std_error
                } else if (life_stage_adult=="Pre-vittelogenic") {
                    # Pre-vittelogenic adult population size.
                    group = previttelogenic_adults;
                    group_std_error = previttelogenic_adults.std_error
                } else if (life_stage_adult=="Vittelogenic") {
                    # Vittelogenic adult population size.
                    group = vittelogenic_adults;
                    group_std_error = vittelogenic_adults.std_error
                } else if (life_stage_adult=="Diapausing") {
                    # Diapausing adult population size.
                    group = diapausing_adults;
                    group_std_error = diapausing_adults.std_error
                }
                maxval = max(group+group_std_error) + 100;
                render_chart(date_labels, "pop_size_by_life_stage", opt$plot_std_error, opt$insect, opt$location, latitude, start_date, end_date, days, maxval,
                    opt$replications, life_stage, group=group, group_std_error=group_std_error, life_stages_adult=life_stage_adult);
                # Turn off device driver to flush output.
                dev.off();
            }
        } else if (life_stage == "Total") {
            # Start PDF device driver.
            dev.new(width=20, height=30);
            file_path = get_file_path(life_stage, "total_pop.pdf")
            pdf(file=file_path, width=20, height=30, bg="white");
            par(mar=c(5, 6, 4, 4), mfrow=c(3, 1));
            # Total population size.
            maxval = max(eggs+eggs.std_error, total_nymphs+total_nymphs.std_error, total_adults+total_adults.std_error) + 100;
            render_chart(date_labels, "pop_size_by_life_stage", opt$plot_std_error, opt$insect, opt$location, latitude, start_date, end_date, days, maxval,
                opt$replications, life_stage, group=total_adults, group_std_error=total_adults.std_error, group2=total_nymphs, group2_std_error=total_nymphs.std_error, group3=eggs,
                group3_std_error=eggs.std_error);
            # Turn off device driver to flush output.
            dev.off();
        }
    }
}
