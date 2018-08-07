#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option(c("--adult_mortality"), action="store", dest="adult_mortality", type="integer", help="Adjustment rate for adult mortality"),
    make_option(c("--adult_accumulation"), action="store", dest="adult_accumulation", type="integer", help="Adjustment of degree-days accumulation (old nymph->adult)"),
    make_option(c("--egg_mortality"), action="store", dest="egg_mortality", type="integer", help="Adjustment rate for egg mortality"),
    make_option(c("--input_norm"), action="store", dest="input_norm", help="30 year normals temperature data for selected station"),
    make_option(c("--input_ytd"), action="store", dest="input_ytd", default=NULL, help="Year-to-date temperature data for selected location"),
    make_option(c("--insect"), action="store", dest="insect", help="Insect name"),
    make_option(c("--insects_per_replication"), action="store", dest="insects_per_replication", type="integer", help="Number of insects with which to start each replication"),
    make_option(c("--life_stages"), action="store", dest="life_stages", help="Selected life stages for plotting"),
    make_option(c("--life_stages_adult"), action="store", dest="life_stages_adult", default=NULL, help="Adult life stages for plotting"),
    make_option(c("--life_stages_nymph"), action="store", dest="life_stages_nymph", default=NULL, help="Nymph life stages for plotting"),
    make_option(c("--location"), action="store", dest="location", default=NULL, help="Selected location"),
    make_option(c("--min_clutch_size"), action="store", dest="min_clutch_size", type="integer", help="Adjustment of minimum clutch size"),
    make_option(c("--max_clutch_size"), action="store", dest="max_clutch_size", type="integer", help="Adjustment of maximum clutch size"),
    make_option(c("--num_days_ytd"), action="store", dest="num_days_ytd", default=NULL, type="integer", help="Total number of days in the year-to-date temperature dataset"),
    make_option(c("--nymph_mortality"), action="store", dest="nymph_mortality", type="integer", help="Adjustment rate for nymph mortality"),
    make_option(c("--old_nymph_accumulation"), action="store", dest="old_nymph_accumulation", type="integer", help="Adjustment of degree-days accumulation (young nymph->old nymph)"),
    make_option(c("--oviposition"), action="store", dest="oviposition", type="integer", help="Adjustment for oviposition rate"),
    make_option(c("--photoperiod"), action="store", dest="photoperiod", type="double", help="Critical photoperiod for diapause induction/termination"),
    make_option(c("--plot_generations_separately"), action="store", dest="plot_generations_separately", help="Plot Plot P, F1 and F2 as separate lines or pool across them"),
    make_option(c("--plot_std_error"), action="store", dest="plot_std_error", help="Plot Standard error"),
    make_option(c("--replications"), action="store", dest="replications", type="integer", help="Number of replications"),
    make_option(c("--script_dir"), action="store", dest="script_dir", help="R script source directory"),
    make_option(c("--young_nymph_accumulation"), action="store", dest="young_nymph_accumulation", type="integer", help="Adjustment of degree-days accumulation (egg->young nymph)")
)

parser <- OptionParser(usage="%prog [options] file", option_list=option_list);
args <- parse_args(parser, positional_arguments=TRUE);
opt <- args$options;

add_daylight_length = function(temperature_data_frame) {
    # Return temperature_data_frame with an added column
    # of daylight length (photoperiod profile).
    num_rows = dim(temperature_data_frame)[1];
    # From Forsythe 1995.
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

from_30_year_normals = function(norm_data_frame, start_date_doy, end_date_doy, year) {
    # The data we want is fully contained within the 30 year normals data.
    first_norm_row = which(norm_data_frame$DOY==start_date_doy);
    last_norm_row = which(norm_data_frame$DOY==end_date_doy);
    # Add 1 to the number of rows to ensure that the end date is included.
    tmp_data_frame_rows = last_norm_row - first_norm_row + 1;
    tmp_data_frame = get_new_temperature_data_frame(nrow=tmp_data_frame_rows);
    j = 0;
    for (i in first_norm_row:last_norm_row) {
        j = j + 1;
        tmp_data_frame[j,] = get_next_normals_row(norm_data_frame, year, i);
    }
    return (tmp_data_frame);
}

get_new_norm_data_frame = function(is_leap_year, input_norm=NULL, nrow=0) {
    # The input_norm data has the following 10 columns:
    # STATIONID, LATITUDE, LONGITUDE, ELEV_M, NAME, ST, MMDD, DOY, TMIN, TMAX
    column_names = c("STATIONID", "LATITUDE","LONGITUDE", "ELEV_M", "NAME", "ST", "MMDD", "DOY", "TMIN", "TMAX");
    if (is.null(input_norm)) {
        norm_data_frame = data.frame(matrix(ncol=10, nrow));
        # Set the norm_data_frame column names for access.
        colnames(norm_data_frame) = column_names;
    } else {
        norm_data_frame = read.csv(file=input_norm, header=T, strip.white=TRUE, stringsAsFactors=FALSE, sep=",");
        # Set the norm_data_frame column names for access.
        colnames(norm_data_frame) = column_names;
        if (!is_leap_year) {
            # All normals data includes Feb 29 which is row 60 in
            # the data, so delete that row if we're not in a leap year.
            norm_data_frame = norm_data_frame[-c(60),];
            # Since we've removed row 60, we need to subtract 1 from
            # each value in the DOY column of the data frame starting
            # with the 60th row.
            num_rows = dim(norm_data_frame)[1];
            for (i in 60:num_rows) {
                leap_year_doy = norm_data_frame$DOY[i];
                non_leap_year_doy = leap_year_doy - 1;
                norm_data_frame$DOY[i] = non_leap_year_doy;
            }
        }
    }
    return (norm_data_frame);
}

get_new_temperature_data_frame = function(input_ytd=NULL, nrow=0) {
    # The input_ytd data has the following 6 columns:
    # LATITUDE, LONGITUDE, DATE, DOY, TMIN, TMAX
    if (is.null(input_ytd)) {
        temperature_data_frame = data.frame(matrix(ncol=6, nrow));
    } else {
        temperature_data_frame = read.csv(file=input_ytd, header=T, strip.white=TRUE, stringsAsFactors=FALSE, sep=",");
    }
    colnames(temperature_data_frame) = c("LATITUDE", "LONGITUDE", "DATE", "DOY", "TMIN", "TMAX");
    return(temperature_data_frame);
}

get_next_normals_row = function(norm_data_frame, year, index) {
    # Return the next 30 year normals row formatted
    # appropriately for the year-to-date data.
    latitude = norm_data_frame[index,"LATITUDE"][1];
    longitude = norm_data_frame[index,"LONGITUDE"][1];
    # Format the date.
    mmdd = norm_data_frame[index,"MMDD"][1];
    date_str = paste(year, mmdd, sep="-");
    doy = norm_data_frame[index,"DOY"][1];
    tmin = norm_data_frame[index,"TMIN"][1];
    tmax = norm_data_frame[index,"TMAX"][1];
    return(list(latitude, longitude, date_str, doy, tmin, tmax));
}

get_temperature_at_hour = function(latitude, temperature_data_frame, row) {
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

is_leap_year = function(date_str) {
    # Extract the year from the date_str.
    date = format(date_str);
    items = strsplit(date, "-")[[1]];
    year = as.integer(items[1]);
    if (((year %% 4 == 0) & (year %% 100 != 0)) | (year %% 400 == 0)) {
        return(TRUE);
    } else {
        return(FALSE);
    }
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

parse_input_data = function(input_ytd, input_norm, location, start_date, end_date) {
    # The end DOY for norm data prepended to ytd data.
    prepend_end_doy_norm = 0;
    # The start DOY for norm data appended to ytd data.
    append_start_doy_norm = 0;
    if (is.null(start_date) && is.null(end_date)) {
        # We're not dealing with a date interval.
        date_interval = FALSE;
        if (is.null(input_ytd)) {
            # Base all dates on the current date since 30 year
            # normals data does not include any dates.
            year = format(Sys.Date(), "%Y");
        }
    } else {
        date_interval = TRUE;
        year = get_year_from_date(start_date);
        # Get the DOY for start_date and end_date.
        start_date_doy = as.integer(strftime(start_date, format="%j"));
        end_date_doy = as.integer(strftime(end_date, format="%j"));
    }
    if (is.null(input_ytd)) {
        # We're processing only the 30 year normals data.
        processing_year_to_date_data = FALSE;
        if (is.null(start_date) && is.null(end_date)) {
            # We're processing the entire year, so we can
            # set the start_date to Jan 1.
            start_date = paste(year, "01", "01", sep="-");
        }
    } else {
        processing_year_to_date_data = TRUE;
        # Read the input_ytd temperature data file into a data frame.
        temperature_data_frame = get_new_temperature_data_frame(input_ytd=input_ytd);
        num_ytd_rows = dim(temperature_data_frame)[1];
        if (!date_interval) {
            start_date = temperature_data_frame$DATE[1];
            year = get_year_from_date(start_date);
        }
    }
    # See if we're in a leap year.
    is_leap_year = is_leap_year(start_date);
    # Read the input_norm temperature datafile into a data frame.
    norm_data_frame = get_new_norm_data_frame(is_leap_year, input_norm=input_norm);
    if (processing_year_to_date_data) {
        if (date_interval) {
            # We're plotting a date interval.
            start_date_ytd_row = which(temperature_data_frame$DATE==start_date);
            if (length(start_date_ytd_row) > 0) {
                # The start date is contained within the input_ytd data.
                start_date_ytd_row = start_date_ytd_row[1];
                start_doy_ytd = as.integer(temperature_data_frame$DOY[start_date_ytd_row]);
            } else {
                # The start date is contained within the input_norm data.
                start_date_ytd_row = 0;
                start_date_norm_row = which(norm_data_frame$DOY==start_date_doy);
            }
            end_date_ytd_row = which(temperature_data_frame$DATE==end_date);
            if (length(end_date_ytd_row) > 0) {
                end_date_ytd_row = end_date_ytd_row[1];
                # The end date is contained within the input_ytd data.
                end_doy_ytd = as.integer(temperature_data_frame$DOY[end_date_ytd_row]);
            } else {
                end_date_ytd_row = 0;
            }
        } else {
            # We're plotting an entire year.
            # Get the start date and end date from temperature_data_frame.
            start_date_ytd_row = 1;
            # Temporarily set start_date to get the year.
            start_date = temperature_data_frame$DATE[1];
            end_date_ytd_row = num_ytd_rows;
            end_date = temperature_data_frame$DATE[num_ytd_rows];
            date_str = format(start_date);
            # Extract the year from the start date.
            date_str_items = strsplit(date_str, "-")[[1]];
            # Get the year.
            year = date_str_items[1];
            # Properly set the start_date to be Jan 1 of the year.
            start_date = paste(year, "01", "01", sep="-");
            # Properly set the end_date to be Dec 31 of the year.
            end_date = paste(year, "12", "31", sep="-");
            # Save the first DOY to later check if start_date is Jan 1.
            start_doy_ytd = as.integer(temperature_data_frame$DOY[1]);
            end_doy_ytd = as.integer(temperature_data_frame$DOY[num_ytd_rows]);
        }
    } else {
        # We're processing only the 30 year normals data, so create an empty
        # data frame for containing temperature data after it is converted
        # from the 30 year normals format to the year-to-date format.
        temperature_data_frame = get_new_temperature_data_frame();
        if (date_interval) {
            # We're plotting a date interval.
            # Extract the year, month and day from the start date.
            start_date_str = format(start_date);
            start_date_str_items = strsplit(start_date_str, "-")[[1]];
            year = start_date_str_items[1];
            start_date_month = start_date_str_items[2];
            start_date_day = start_date_str_items[3];
            start_date = paste(year, start_date_month, start_date_day, sep="-");
            # Extract the month and day from the end date.
            end_date_str = format(start_date);
            end_date_str_items = strsplit(end_date_str, "-")[[1]];
            end_date_month = end_date_str_items[2];
            end_date_day = end_date_str_items[3];
            end_date = paste(year, end_date_month, end_date_day, sep="-");
        } else {
            # We're plotting an entire year.
            start_date = paste(year, "01", "01", sep="-");
            end_date = paste(year, "12", "31", sep="-");
        }
    }
    # Set the location to be the station name if the user elected not to enter it.
    if (is.null(location) | length(location) == 0) {
        location = norm_data_frame$NAME[1];
    }
    if (processing_year_to_date_data) {
        # Merge the year-to-date data with the 30 year normals data.
        if (date_interval) {
            # The values of start_date_ytd_row and end_date_ytd_row were set above.
            if (start_date_ytd_row > 0 & end_date_ytd_row > 0) {
                # The date interval is contained within the input_ytd
                # data, so we don't need to merge the 30 year normals data.
                temperature_data_frame = temperature_data_frame[start_date_ytd_row:end_date_ytd_row,];
            } else if (start_date_ytd_row == 0 & end_date_ytd_row > 0) {
                # The date interval starts in input_norm and ends in
                # input_ytd, so prepend appropriate rows from input_norm
                # to appropriate rows from input_ytd.
                first_norm_row = which(norm_data_frame$DOY==start_date_doy);
                # Get the first DOY from temperature_data_frame.
                first_ytd_doy = temperature_data_frame$DOY[1];
                # End DOY of input_norm data prepended to input_ytd.
                prepend_end_doy_norm = first_ytd_doy - 1;
                # Get the number of rows for the restricted date interval
                # that are contained in temperature_data_frame.
                num_temperature_data_frame_rows = end_date_ytd_row;
                # Get the last row needed from the 30 year normals data.
                last_norm_row = which(norm_data_frame$DOY==prepend_end_doy_norm);
                # Get the number of rows for the restricted date interval
                # that are contained in norm_data_frame.
                num_norm_data_frame_rows = last_norm_row - first_norm_row;
                # Create a temporary data frame to contain the 30 year normals
                # data from the start date to the date immediately prior to the
                # first row of the input_ytd data.
                tmp_norm_data_frame = get_new_temperature_data_frame(nrow=num_temperature_data_frame_rows+num_norm_data_frame_rows);
                j = 1;
                for (i in first_norm_row:last_norm_row) {
                    # Populate the temp_data_frame row with
                    # values from norm_data_frame.
                    tmp_norm_data_frame[j,] = get_next_normals_row(norm_data_frame, year, i);
                    j = j + 1;
                }
                # Create a second temporary data frame containing the
                # appropriate rows from temperature_data_frame.
                tmp_temperature_data_frame = temperature_data_frame[1:num_temperature_data_frame_rows,];
                # Merge the 2 temporary data frames.
                temperature_data_frame = rbind(tmp_norm_data_frame, tmp_temperature_data_frame);
            } else if (start_date_ytd_row > 0 & end_date_ytd_row == 0) {
                # The date interval starts in input_ytd and ends in input_norm,
                # so append appropriate rows from input_norm to appropriate rows
                # from input_ytd. First, get the number of rows for the restricted
                # date interval that are contained in temperature_data_frame.
                num_temperature_data_frame_rows = num_ytd_rows - start_date_ytd_row + 1;
                # Get the DOY of the last row in the input_ytd data.
                last_ytd_doy = temperature_data_frame$DOY[num_ytd_rows];
                # Get the DOYs for the first and last rows from norm_data_frame
                # that will be appended to temperature_data_frame.
                append_start_doy_norm = last_ytd_doy + 1;
                # Get the row from norm_data_frame containing first_norm_doy.
                first_norm_row = which(norm_data_frame$DOY == append_start_doy_norm);
                # Get the row from norm_data_frame containing end_date_doy.
                last_norm_row = which(norm_data_frame$DOY == end_date_doy);
                # Get the number of rows for the restricted date interval
                # that are contained in norm_data_frame.
                num_norm_data_frame_rows = last_norm_row - first_norm_row;
                # Create a temporary data frame to contain the data
                # taken from both temperature_data_frame and norm_data_frame
                # for the date interval.
                tmp_data_frame = get_new_temperature_data_frame(nrow=num_temperature_data_frame_rows+num_norm_data_frame_rows);
                # Populate tmp_data_frame with the appropriate rows from temperature_data_frame.
                j = start_date_ytd_row;
                for (i in 1:num_temperature_data_frame_rows) {
                    tmp_data_frame[i,] = temperature_data_frame[j,];
                    j = j + 1;
                }
                # Apppend the appropriate rows from norm_data_frame to tmp_data_frame.
                current_iteration = num_temperature_data_frame_rows + 1;
                num_iterations = current_iteration + num_norm_data_frame_rows;
                j = first_norm_row;
                for (i in current_iteration:num_iterations) {
                    tmp_data_frame[i,] = get_next_normals_row(norm_data_frame, year, j);
                    j = j + 1;
                }
                temperature_data_frame = tmp_data_frame[,];
            } else if (start_date_ytd_row == 0 & end_date_ytd_row == 0) {
                # The date interval is contained witin input_norm.
                temperature_data_frame = from_30_year_normals(norm_data_frame, start_date_doy, end_date_doy, year);
            }
        } else {
            # We're plotting an entire year.
            if (start_doy_ytd > 1) {
                # The input_ytd data starts after Jan 1, so prepend
                # appropriate rows from input_norm to temperature_data_frame.
                prepend_end_doy_norm = start_doy_ytd - 1;
                last_norm_row = which(norm_data_frame$DOY == prepend_end_doy_norm);
                # Create a temporary data frame to contain the input_norm data
                # from Jan 1 to the date immediately prior to start_date.
                tmp_data_frame = temperature_data_frame[FALSE,];
                # Populate tmp_data_frame with appropriate rows from norm_data_frame.
                for (i in 1:last_norm_row) {
                    tmp_data_frame[i,] = get_next_normals_row(norm_data_frame, year, i);
                }
                # Merge the temporary data frame with temperature_data_frame.
                temperature_data_frame = rbind(tmp_data_frame, temperature_data_frame);
            }
            # Set the value of total_days.
            total_days = get_total_days(is_leap_year);
            if (end_doy_ytd < total_days) {
                # Define the next row for the year-to-date data from the 30 year normals data.
                append_start_doy_norm = end_doy_ytd + 1;
                first_norm_row = which(norm_data_frame$DOY == append_start_doy_norm);
                # Append the 30 year normals data to the year-to-date data.
                for (i in first_norm_row:total_days) {
                    temperature_data_frame[i,] = get_next_normals_row(norm_data_frame, year, i);
                }
            }
        }
    } else {
        # We're processing only the 30 year normals data.
        if (date_interval) {
            # Populate temperature_data_frame from norm_data_frame.
            temperature_data_frame = from_30_year_normals(norm_data_frame, start_date_doy, end_date_doy, year);
        } else {
            total_days = get_total_days(is_leap_year);
            for (i in 1:total_days) {
                temperature_data_frame[i,] = get_next_normals_row(norm_data_frame, year, i);
            }
        }
    }
    # Add a column containing the daylight length for each day.
    temperature_data_frame = add_daylight_length(temperature_data_frame);
    return(list(temperature_data_frame, start_date, end_date, prepend_end_doy_norm, append_start_doy_norm, is_leap_year, location));
}

# Import the shared utility functions.
utils_path <- paste(opt$script_dir, "utils.R", sep="/");
source(utils_path);

if (is.null(opt$input_ytd)) {
    processing_year_to_date_data = FALSE;
} else {
    processing_year_to_date_data = TRUE;
}
# Determine if we're plotting generations separately.
if (opt$plot_generations_separately=="yes") {
    plot_generations_separately = TRUE;
} else {
    plot_generations_separately = FALSE;
}
# Parse the inputs.
data_list = parse_input_data(opt$input_ytd, opt$input_norm, opt$location, opt$start_date, opt$end_date);
temperature_data_frame = data_list[[1]];
# Information needed for plots, some of these values are
# being reset here since in some case they were set above.
start_date = data_list[[2]];
end_date = data_list[[3]];
prepend_end_doy_norm = data_list[[4]];
append_start_doy_norm = data_list[[5]];
is_leap_year = data_list[[6]];
location = data_list[[7]];

# We're plotting an entire year.
# Display the total number of days in the Galaxy history item blurb.
if (processing_year_to_date_data) {
    cat("Number of days year-to-date: ", opt$num_days_ytd, "\n");
} else {
    if (is_leap_year) {
        num_days = 366;
    } else {
        num_days = 365;
    }
    cat("Number of days in year: ", num_days, "\n");
}

# Create copies of the temperature data for generations P, F1 and F2 if we're plotting generations separately.
if (plot_generations_separately) {
    temperature_data_frame_P = data.frame(temperature_data_frame);
    temperature_data_frame_F1 = data.frame(temperature_data_frame);
    temperature_data_frame_F2 = data.frame(temperature_data_frame);
}

# Get the ticks date labels for plots.
ticks_and_labels = get_x_axis_ticks_and_labels(temperature_data_frame, prepend_end_doy_norm, append_start_doy_norm);
ticks = c(unlist(ticks_and_labels[1]));
date_labels = c(unlist(ticks_and_labels[2]));
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
total_days = dim(temperature_data_frame)[1];
if (process_eggs) {
    Eggs.replications = matrix(rep(0, total_days*opt$replications), ncol=opt$replications);
}
if (process_young_nymphs | process_total_nymphs) {
    YoungNymphs.replications = matrix(rep(0, total_days*opt$replications), ncol=opt$replications);
}
if (process_old_nymphs | process_total_nymphs) {
    OldNymphs.replications = matrix(rep(0, total_days*opt$replications), ncol=opt$replications);
}
if (process_previttelogenic_adults | process_total_adults) {
    Previttelogenic.replications = matrix(rep(0, total_days*opt$replications), ncol=opt$replications);
}
if (process_vittelogenic_adults | process_total_adults) {
    Vittelogenic.replications = matrix(rep(0, total_days*opt$replications), ncol=opt$replications);
}
if (process_diapausing_adults | process_total_adults) {
    Diapausing.replications = matrix(rep(0, total_days*opt$replications), ncol=opt$replications);
}
newborn.replications = matrix(rep(0, total_days*opt$replications), ncol=opt$replications);
adult.replications = matrix(rep(0, total_days*opt$replications), ncol=opt$replications);
death.replications = matrix(rep(0, total_days*opt$replications), ncol=opt$replications);
if (plot_generations_separately) {
    # P is Parental, or overwintered adults.
    P.replications = matrix(rep(0, total_days*opt$replications), ncol=opt$replications);
    # F1 is the first field-produced generation.
    F1.replications = matrix(rep(0, total_days*opt$replications), ncol=opt$replications);
    # F2 is the second field-produced generation.
    F2.replications = matrix(rep(0, total_days*opt$replications), ncol=opt$replications);
    if (process_eggs) {
        P_eggs.replications = matrix(rep(0, total_days*opt$replications), ncol=opt$replications);
        F1_eggs.replications = matrix(rep(0, total_days*opt$replications), ncol=opt$replications);
        F2_eggs.replications = matrix(rep(0, total_days*opt$replications), ncol=opt$replications);
    }
    if (process_young_nymphs) {
        P_young_nymphs.replications = matrix(rep(0, total_days*opt$replications), ncol=opt$replications);
        F1_young_nymphs.replications = matrix(rep(0, total_days*opt$replications), ncol=opt$replications);
        F2_young_nymphs.replications = matrix(rep(0, total_days*opt$replications), ncol=opt$replications);
    }
    if (process_old_nymphs) {
        P_old_nymphs.replications = matrix(rep(0, total_days*opt$replications), ncol=opt$replications);
        F1_old_nymphs.replications = matrix(rep(0, total_days*opt$replications), ncol=opt$replications);
        F2_old_nymphs.replications = matrix(rep(0, total_days*opt$replications), ncol=opt$replications);
    }
    if (process_total_nymphs) {
        P_total_nymphs.replications = matrix(rep(0, total_days*opt$replications), ncol=opt$replications);
        F1_total_nymphs.replications = matrix(rep(0, total_days*opt$replications), ncol=opt$replications);
        F2_total_nymphs.replications = matrix(rep(0, total_days*opt$replications), ncol=opt$replications);
    }
    if (process_previttelogenic_adults) {
        P_previttelogenic_adults.replications = matrix(rep(0, total_days*opt$replications), ncol=opt$replications);
        F1_previttelogenic_adults.replications = matrix(rep(0, total_days*opt$replications), ncol=opt$replications);
        F2_previttelogenic_adults.replications = matrix(rep(0, total_days*opt$replications), ncol=opt$replications);
    }
    if (process_vittelogenic_adults) {
        P_vittelogenic_adults.replications = matrix(rep(0, total_days*opt$replications), ncol=opt$replications);
        F1_vittelogenic_adults.replications = matrix(rep(0, total_days*opt$replications), ncol=opt$replications);
        F2_vittelogenic_adults.replications = matrix(rep(0, total_days*opt$replications), ncol=opt$replications);
    }
    if (process_diapausing_adults) {
        P_diapausing_adults.replications = matrix(rep(0, total_days*opt$replications), ncol=opt$replications);
        F1_diapausing_adults.replications = matrix(rep(0, total_days*opt$replications), ncol=opt$replications);
        F2_diapausing_adults.replications = matrix(rep(0, total_days*opt$replications), ncol=opt$replications);
    }
    if (process_total_adults) {
        P_total_adults.replications = matrix(rep(0, total_days*opt$replications), ncol=opt$replications);
        F1_total_adults.replications = matrix(rep(0, total_days*opt$replications), ncol=opt$replications);
        F2_total_adults.replications = matrix(rep(0, total_days*opt$replications), ncol=opt$replications);
    }
}
# Total population.
population.replications = matrix(rep(0, total_days*opt$replications), ncol=opt$replications);

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
        Eggs = rep(0, total_days);
    }
    if (process_young_nymphs | process_total_nymphs) {
        YoungNymphs = rep(0, total_days);
    }
    if (process_old_nymphs | process_total_nymphs) {
        OldNymphs = rep(0, total_days);
    }
    if (process_previttelogenic_adults | process_total_adults) {
        Previttelogenic = rep(0, total_days);
    }
    if (process_vittelogenic_adults | process_total_adults) {
        Vittelogenic = rep(0, total_days);
    }
    if (process_diapausing_adults | process_total_adults) {
        Diapausing = rep(0, total_days);
    }
    N.newborn = rep(0, total_days);
    N.adult = rep(0, total_days);
    N.death = rep(0, total_days);
    overwintering_adult.population = rep(0, total_days);
    first_generation.population = rep(0, total_days);
    second_generation.population = rep(0, total_days);
    if (plot_generations_separately) {
        # P is Parental, or overwintered adults.
        # F1 is the first field-produced generation.
        # F2 is the second field-produced generation.
        if (process_eggs) {
            P.egg = rep(0, total_days);
            F1.egg = rep(0, total_days);
            F2.egg = rep(0, total_days);
        }
        if (process_young_nymphs) {
            P.young_nymph = rep(0, total_days);
            F1.young_nymph = rep(0, total_days);
            F2.young_nymph = rep(0, total_days);
        }
        if (process_old_nymphs) {
            P.old_nymph = rep(0, total_days);
            F1.old_nymph = rep(0, total_days);
            F2.old_nymph = rep(0, total_days);
        }
        if (process_total_nymphs) {
            P.total_nymph = rep(0, total_days);
            F1.total_nymph = rep(0, total_days);
            F2.total_nymph = rep(0, total_days);
        }
        if (process_previttelogenic_adults) {
            P.previttelogenic_adult = rep(0, total_days);
            F1.previttelogenic_adult = rep(0, total_days);
            F2.previttelogenic_adult = rep(0, total_days);
        }
        if (process_vittelogenic_adults) {
            P.vittelogenic_adult = rep(0, total_days);
            F1.vittelogenic_adult = rep(0, total_days);
            F2.vittelogenic_adult = rep(0, total_days);
        }
        if (process_diapausing_adults) {
            P.diapausing_adult = rep(0, total_days);
            F1.diapausing_adult = rep(0, total_days);
            F2.diapausing_adult = rep(0, total_days);
        }
        if (process_total_adults) {
            P.total_adult = rep(0, total_days);
            F1.total_adult = rep(0, total_days);
            F2.total_adult = rep(0, total_days);
        }
    }
    total.population = NULL;
    averages.day = rep(0, total_days);
    # All the days included in the input_ytd temperature dataset.
    for (row in 1:total_days) {
        # Get the integer day of the year for the current row.
        doy = temperature_data_frame$DOY[row];
        # Photoperiod in the day.
        photoperiod = temperature_data_frame$DAYLEN[row];
        temp.profile = get_temperature_at_hour(latitude, temperature_data_frame, row);
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
    }   # End of days specified in the input_ytd temperature data.

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
        if (life_stage_nymph=="Young") {
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
        } else if (life_stage_nymph=="Total") {
            # Mean value for all nymphs.
            total_nymphs = apply((YoungNymphs.replications+OldNymphs.replications), 1, mean);
            temperature_data_frame = append_vector(temperature_data_frame, total_nymphs, "TOTALNYMPH");
            # Standard error for all nymphs.
            total_nymphs.std_error = apply((YoungNymphs.replications+OldNymphs.replications) / sqrt(opt$replications), 1, sd);
            temperature_data_frame = append_vector(temperature_data_frame, total_nymphs.std_error, "TOTALNYMPHSE");
        }
    }
}
if (process_adults) {
    # Calculate adult populations for selected life stage.
    for (life_stage_adult in life_stages_adult) {
        if (life_stage_adult == "Pre-vittelogenic") {
            # Mean value for previttelogenic adults.
            previttelogenic_adults = apply(Previttelogenic.replications, 1, mean);
            temperature_data_frame = append_vector(temperature_data_frame, previttelogenic_adults, "PRE.VITADULT");
            # Standard error for previttelogenic adults.
            previttelogenic_adults.std_error = apply(Previttelogenic.replications, 1, sd) / sqrt(opt$replications);
            temperature_data_frame = append_vector(temperature_data_frame, previttelogenic_adults.std_error, "PRE.VITADULTSE");
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
        } else if (life_stage_adult=="Total") {
            # Mean value for all adults.
            total_adults = apply((Previttelogenic.replications+Vittelogenic.replications+Diapausing.replications), 1, mean);
            temperature_data_frame = append_vector(temperature_data_frame, total_adults, "TOTALADULT");
            # Standard error for all adults.
            total_adults.std_error = apply((Previttelogenic.replications+Vittelogenic.replications+Diapausing.replications), 1, sd) / sqrt(opt$replications);
            temperature_data_frame = append_vector(temperature_data_frame, total_adults.std_error, "TOTALADULTSE");
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
        temperature_data_frame_P = append_vector(temperature_data_frame_P, P_eggs, "EGG.P");
        temperature_data_frame_P = append_vector(temperature_data_frame_P, P_eggs.std_error, "EGG.P.SE");
        F1_eggs = m_se[[3]];
        F1_eggs.std_error = m_se[[4]];
        temperature_data_frame_F1 = append_vector(temperature_data_frame_F1, F1_eggs, "EGG.F1");
        temperature_data_frame_F1 = append_vector(temperature_data_frame_F1, F1_eggs.std_error, "EGG.F1.SE");
        F2_eggs = m_se[[5]];
        F2_eggs.std_error = m_se[[6]];
        temperature_data_frame_F2 = append_vector(temperature_data_frame_F2, F2_eggs, "EGG.F2");
        temperature_data_frame_F2 = append_vector(temperature_data_frame_F2, F2_eggs.std_error, "EGG.F2.SE");
    }
    if (process_young_nymphs) {
        m_se = get_mean_and_std_error(P_young_nymphs.replications, F1_young_nymphs.replications, F2_young_nymphs.replications);
        P_young_nymphs = m_se[[1]];
        P_young_nymphs.std_error = m_se[[2]];
        temperature_data_frame_P = append_vector(temperature_data_frame_P, P_young_nymphs, "YOUNGNYMPH.P");
        temperature_data_frame_P = append_vector(temperature_data_frame_P, P_young_nymphs.std_error, "YOUNGNYMPH.P.SE");
        F1_young_nymphs = m_se[[3]];
        F1_young_nymphs.std_error = m_se[[4]];
        temperature_data_frame_F1 = append_vector(temperature_data_frame_F1, F1_young_nymphs, "YOUNGNYMPH.F1");
        temperature_data_frame_F1 = append_vector(temperature_data_frame_F1, F1_young_nymphs.std_error, "YOUNGNYMPH.F1.SE");
        F2_young_nymphs = m_se[[5]];
        F2_young_nymphs.std_error = m_se[[6]];
        temperature_data_frame_F2 = append_vector(temperature_data_frame_F2, F2_young_nymphs, "YOUNGNYMPH.F2");
        temperature_data_frame_F2 = append_vector(temperature_data_frame_F2, F2_young_nymphs.std_error, "YOUNGNYMPH.F2.SE");
    }
    if (process_old_nymphs) {
        m_se = get_mean_and_std_error(P_old_nymphs.replications, F1_old_nymphs.replications, F2_old_nymphs.replications);
        P_old_nymphs = m_se[[1]];
        P_old_nymphs.std_error = m_se[[2]];
        temperature_data_frame_P = append_vector(temperature_data_frame_P, P_old_nymphs, "OLDNYMPH.P");
        temperature_data_frame_P = append_vector(temperature_data_frame_P, P_old_nymphs.std_error, "OLDNYMPH.P.SE");
        F1_old_nymphs = m_se[[3]];
        F1_old_nymphs.std_error = m_se[[4]];
        temperature_data_frame_F1 = append_vector(temperature_data_frame_F1, F1_old_nymphs, "OLDNYMPH.F1");
        temperature_data_frame_F1 = append_vector(temperature_data_frame_F1, F1_old_nymphs.std_error, "OLDNYMPH.F1.SE");
        F2_old_nymphs = m_se[[5]];
        F2_old_nymphs.std_error = m_se[[6]];
        temperature_data_frame_F2 = append_vector(temperature_data_frame_F2, F2_old_nymphs, "OLDNYMPH.F2");
        temperature_data_frame_F2 = append_vector(temperature_data_frame_F2, F2_old_nymphs.std_error, "OLDNYMPH.F2.SE");
    }
    if (process_total_nymphs) {
        m_se = get_mean_and_std_error(P_total_nymphs.replications, F1_total_nymphs.replications, F2_total_nymphs.replications);
        P_total_nymphs = m_se[[1]];
        P_total_nymphs.std_error = m_se[[2]];
        temperature_data_frame_P = append_vector(temperature_data_frame_P, P_total_nymphs, "TOTALNYMPH.P");
        temperature_data_frame_P = append_vector(temperature_data_frame_P, P_total_nymphs.std_error, "TOTALNYMPH.P.SE");
        F1_total_nymphs = m_se[[3]];
        F1_total_nymphs.std_error = m_se[[4]];
        temperature_data_frame_F1 = append_vector(temperature_data_frame_F1, F1_total_nymphs, "TOTALNYMPH.F1");
        temperature_data_frame_F1 = append_vector(temperature_data_frame_F1, F1_total_nymphs.std_error, "TOTALNYMPH.F1.SE");
        F2_total_nymphs = m_se[[5]];
        F2_total_nymphs.std_error = m_se[[6]];
        temperature_data_frame_F2 = append_vector(temperature_data_frame_F2, F2_total_nymphs, "TOTALNYMPH.F2");
        temperature_data_frame_F2 = append_vector(temperature_data_frame_F2, F2_total_nymphs.std_error, "TOTALNYMPH.F2.SE");
    }
    if (process_previttelogenic_adults) {
        m_se = get_mean_and_std_error(P_previttelogenic_adults.replications, F1_previttelogenic_adults.replications, F2_previttelogenic_adults.replications);
        P_previttelogenic_adults = m_se[[1]];
        P_previttelogenic_adults.std_error = m_se[[2]];
        temperature_data_frame_P = append_vector(temperature_data_frame_P, P_previttelogenic_adults, "PRE.VITADULT.P");
        temperature_data_frame_P = append_vector(temperature_data_frame_P, P_previttelogenic_adults.std_error, "PRE.VITADULT.P.SE");
        F1_previttelogenic_adults = m_se[[3]];
        F1_previttelogenic_adults.std_error = m_se[[4]];
        temperature_data_frame_F1 = append_vector(temperature_data_frame_F1, F1_previttelogenic_adults, "PRE.VITADULT.F1");
        temperature_data_frame_F1 = append_vector(temperature_data_frame_F1, F1_previttelogenic_adults.std_error, "PRE.VITADULT.F1.SE");
        F2_previttelogenic_adults = m_se[[5]];
        F2_previttelogenic_adults.std_error = m_se[[6]];
        temperature_data_frame_F2 = append_vector(temperature_data_frame_F2, F2_previttelogenic_adults, "PRE.VITADULT.F2");
        temperature_data_frame_F2 = append_vector(temperature_data_frame_F2, F2_previttelogenic_adults.std_error, "PRE.VITADULT.F2.SE");
    }
    if (process_vittelogenic_adults) {
        m_se = get_mean_and_std_error(P_vittelogenic_adults.replications, F1_vittelogenic_adults.replications, F2_vittelogenic_adults.replications);
        P_vittelogenic_adults = m_se[[1]];
        P_vittelogenic_adults.std_error = m_se[[2]];
        temperature_data_frame_P = append_vector(temperature_data_frame_P, P_vittelogenic_adults, "VITADULT.P");
        temperature_data_frame_P = append_vector(temperature_data_frame_P, P_vittelogenic_adults.std_error, "VITADULT.P.SE");
        F1_vittelogenic_adults = m_se[[3]];
        F1_vittelogenic_adults.std_error = m_se[[4]];
        temperature_data_frame_F1 = append_vector(temperature_data_frame_F1, F1_vittelogenic_adults, "VITADULT.F1");
        temperature_data_frame_F1 = append_vector(temperature_data_frame_F1, F1_vittelogenic_adults.std_error, "VITADULT.F1.SE");
        F2_vittelogenic_adults = m_se[[5]];
        F2_vittelogenic_adults.std_error = m_se[[6]];
        temperature_data_frame_F2 = append_vector(temperature_data_frame_F2, F2_vittelogenic_adults, "VITADULT.F2");
        temperature_data_frame_F2 = append_vector(temperature_data_frame_F2, F2_vittelogenic_adults.std_error, "VITADULT.F2.SE");
    }
    if (process_diapausing_adults) {
        m_se = get_mean_and_std_error(P_diapausing_adults.replications, F1_diapausing_adults.replications, F2_diapausing_adults.replications);
        P_diapausing_adults = m_se[[1]];
        P_diapausing_adults.std_error = m_se[[2]];
        temperature_data_frame_P = append_vector(temperature_data_frame_P, P_diapausing_adults, "DIAPAUSINGADULT.P");
        temperature_data_frame_P = append_vector(temperature_data_frame_P, P_diapausing_adults.std_error, "DIAPAUSINGADULT.P.SE");
        F1_diapausing_adults = m_se[[3]];
        F1_diapausing_adults.std_error = m_se[[4]];
        temperature_data_frame_F1 = append_vector(temperature_data_frame_F1, F1_diapausing_adults, "DIAPAUSINGADULT.F1");
        temperature_data_frame_F1 = append_vector(temperature_data_frame_F1, F1_diapausing_adults.std_error, "DIAPAUSINGADULT.F1.SE");
        F2_diapausing_adults = m_se[[5]];
        F2_diapausing_adults.std_error = m_se[[6]];
        temperature_data_frame_F2 = append_vector(temperature_data_frame_F2, F2_diapausing_adults, "DIAPAUSINGADULT.F2");
        temperature_data_frame_F2 = append_vector(temperature_data_frame_F2, F2_diapausing_adults.std_error, "DIAPAUSINGADULT.F2.SE");
    }
    if (process_total_adults) {
        m_se = get_mean_and_std_error(P_total_adults.replications, F1_total_adults.replications, F2_total_adults.replications);
        P_total_adults = m_se[[1]];
        P_total_adults.std_error = m_se[[2]];
        temperature_data_frame_P = append_vector(temperature_data_frame_P, P_total_adults, "TOTALADULT.P");
        temperature_data_frame_P = append_vector(temperature_data_frame_P, P_total_adults.std_error, "TOTALADULT.P.SE");
        F1_total_adults = m_se[[3]];
        F1_total_adults.std_error = m_se[[4]];
        temperature_data_frame_F1 = append_vector(temperature_data_frame_F1, F1_total_adults, "TOTALADULT.F1");
        temperature_data_frame_F1 = append_vector(temperature_data_frame_F1, F1_total_adults.std_error, "TOTALADULT.F1.SE");
        F2_total_adults = m_se[[5]];
        F2_total_adults.std_error = m_se[[6]];
        temperature_data_frame_F2 = append_vector(temperature_data_frame_F2, F2_total_adults, "TOTALADULT.F2");
        temperature_data_frame_F2 = append_vector(temperature_data_frame_F2, F2_total_adults.std_error, "TOTALADULT.F2.SE");
    }
}

# Save the analyzed data for combined generations.
file_path = paste("output_data_dir", "04_combined_generations.csv", sep="/");
write.csv(temperature_data_frame, file=file_path, row.names=F);
if (plot_generations_separately) {
    # Save the analyzed data for generation P.
    file_path = paste("output_data_dir", "01_generation_P.csv", sep="/");
    write.csv(temperature_data_frame_P, file=file_path, row.names=F);
    # Save the analyzed data for generation F1.
    file_path = paste("output_data_dir", "02_generation_F1.csv", sep="/");
    write.csv(temperature_data_frame_F1, file=file_path, row.names=F);
    # Save the analyzed data for generation F2.
    file_path = paste("output_data_dir", "03_generation_F2.csv", sep="/");
    write.csv(temperature_data_frame_F2, file=file_path, row.names=F);
}

total_days_vector = c(1:dim(temperature_data_frame)[1]);
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
            render_chart(ticks, date_labels, "pop_size_by_generation", opt$plot_std_error, opt$insect, location, latitude,
                start_date, end_date, total_days_vector, maxval, opt$replications, life_stage, group=P_eggs, group_std_error=P_eggs.std_error,
                group2=F1_eggs, group2_std_error=F1_eggs.std_error, group3=F2_eggs, group3_std_error=F2_eggs.std_error);
            # Turn off device driver to flush output.
            dev.off();
        } else if (life_stage == "Nymph") {
            for (life_stage_nymph in life_stages_nymph) {
                # Start PDF device driver.
                dev.new(width=20, height=30);
                file_path = get_file_path(life_stage, "nymph_pop_by_generation.pdf", sub_life_stage=life_stage_nymph)
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
                render_chart(ticks, date_labels, "pop_size_by_generation", opt$plot_std_error, opt$insect, location, latitude,
                    start_date, end_date, total_days_vector, maxval, opt$replications, life_stage, group=group, group_std_error=group_std_error,
                    group2=group2, group2_std_error=group2_std_error, group3=group3, group3_std_error=group3_std_error, sub_life_stage=life_stage_nymph);
                # Turn off device driver to flush output.
                dev.off();
            }
        } else if (life_stage == "Adult") {
            for (life_stage_adult in life_stages_adult) {
                # Start PDF device driver.
                dev.new(width=20, height=30);
                file_path = get_file_path(life_stage, "adult_pop_by_generation.pdf", sub_life_stage=life_stage_adult)
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
                render_chart(ticks, date_labels, "pop_size_by_generation", opt$plot_std_error, opt$insect, location, latitude,
                    start_date, end_date, total_days_vector, maxval, opt$replications, life_stage, group=group, group_std_error=group_std_error,
                    group2=group2, group2_std_error=group2_std_error, group3=group3, group3_std_error=group3_std_error, sub_life_stage=life_stage_adult);
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
            render_chart(ticks, date_labels, "pop_size_by_generation", opt$plot_std_error, opt$insect, location, latitude,
                start_date, end_date, total_days_vector, maxval, opt$replications, life_stage, group=P, group_std_error=P.std_error,
                group2=F1, group2_std_error=F1.std_error, group3=F2, group3_std_error=F2.std_error);
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
            render_chart(ticks, date_labels, "pop_size_by_life_stage", opt$plot_std_error, opt$insect, location, latitude,
                start_date, end_date, total_days_vector, maxval, opt$replications, life_stage, group=eggs, group_std_error=eggs.std_error);
            # Turn off device driver to flush output.
            dev.off();
        } else if (life_stage == "Nymph") {
            for (life_stage_nymph in life_stages_nymph) {
                # Start PDF device driver.
                dev.new(width=20, height=30);
                file_path = get_file_path(life_stage, "nymph_pop.pdf", sub_life_stage=life_stage_nymph)
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
                render_chart(ticks, date_labels, "pop_size_by_life_stage", opt$plot_std_error, opt$insect, location, latitude,
                    start_date, end_date, total_days_vector, maxval, opt$replications, life_stage, group=group, group_std_error=group_std_error,
                    sub_life_stage=life_stage_nymph);
                # Turn off device driver to flush output.
                dev.off();
            }
        } else if (life_stage == "Adult") {
            for (life_stage_adult in life_stages_adult) {
                # Start PDF device driver.
                dev.new(width=20, height=30);
                file_path = get_file_path(life_stage, "adult_pop.pdf", sub_life_stage=life_stage_adult)
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
                render_chart(ticks, date_labels, "pop_size_by_life_stage", opt$plot_std_error, opt$insect, location, latitude,
                    start_date, end_date, total_days_vector, maxval, opt$replications, life_stage, group=group, group_std_error=group_std_error,
                    sub_life_stage=life_stage_adult);
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
            render_chart(ticks, date_labels, "pop_size_by_life_stage", opt$plot_std_error, opt$insect, location, latitude,
                start_date, end_date, total_days_vector, maxval, opt$replications, life_stage, group=total_adults, group_std_error=total_adults.std_error,
                group2=total_nymphs, group2_std_error=total_nymphs.std_error, group3=eggs, group3_std_error=eggs.std_error);
            # Turn off device driver to flush output.
            dev.off();
        }
    }
}
