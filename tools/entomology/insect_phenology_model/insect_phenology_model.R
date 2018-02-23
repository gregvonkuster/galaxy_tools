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
    make_option(c("--nymph_mortality"), action="store", dest="nymph_mortality", type="integer", help="Adjustment rate for nymph mortality"),
    make_option(c("--old_nymph_accumulation"), action="store", dest="old_nymph_accumulation", type="integer", help="Adjustment of degree-days accumulation (young nymph->old nymph)"),
    make_option(c("--num_days"), action="store", dest="num_days", type="integer", help="Total number of days in the temperature dataset"),
    make_option(c("--oviposition"), action="store", dest="oviposition", type="integer", help="Adjustment for oviposition rate"),
    make_option(c("--photoperiod"), action="store", dest="photoperiod", type="double", help="Critical photoperiod for diapause induction/termination"),
    make_option(c("--replications"), action="store", dest="replications", type="integer", help="Number of replications"),
    make_option(c("--plot_generations_separately"), action="store", dest="plot_generations_separately", help="Plot Plot P, F1 and F2 as separate lines or pool across them"),
    make_option(c("--plot_std_error"), action="store", dest="plot_std_error", help="Plot Standard error"),
    make_option(c("--young_nymph_accumulation"), action="store", dest="young_nymph_accumulation", type="integer", help="Adjustment of degree-days accumulation (egg->young nymph)")
)

parser <- OptionParser(usage="%prog [options] file", option_list=option_list);
args <- parse_args(parser, positional_arguments=TRUE);
opt <- args$options;

add_daylight_length = function(temperature_data_frame, num_columns, num_rows) {
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
    temperature_data_frame[, num_columns+1] = daylight_length_vector;
    return(temperature_data_frame);
}

dev.egg = function(temperature) {
    dev.rate = -0.9843 * temperature + 33.438;
    return(dev.rate);
}

dev.emerg = function(temperature) {
    emerg.rate = -0.5332 * temperature + 24.147;
    return(emerg.rate);
}

dev.old = function(temperature) {
    n34 = -0.6119 * temperature + 17.602;
    n45 = -0.4408 * temperature + 19.036;
    dev.rate = mean(n34 + n45);
    return(dev.rate);
}

dev.young = function(temperature) {
    n12 = -0.3728 * temperature + 14.68;
    n23 = -0.6119 * temperature + 25.249;
    dev.rate = mean(n12 + n23);
    return(dev.rate);
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
        # Add a column containing the daylight length for each day.
        temperature_data_frame = add_daylight_length(temperature_data_frame, num_columns, num_rows);
        # Reset the column names with the additional column for later access.
        colnames(temperature_data_frame) = c("LATITUDE","LONGITUDE", "DATE", "DOY", "TMIN", "TMAX", "DAYLEN");
    }
    return(temperature_data_frame);
}


render_chart = function(date_labels, chart_type, plot_std_error, insect, location, latitude, start_date, end_date, days, maxval,
    life_stage, group, group_std_error, group2=NULL, group2_std_error=NULL, group3=NULL, group3_std_error=NULL,
    life_stages_adult=NULL, life_stages_nymph=NULL) {
    cat("RRR chart_type: ", chart_type, "\n");
    if (chart_type=="pop_size_by_life_stage") {
        if (life_stage=="All") {
            title = paste(insect, ": Total Pop :", location, ": Lat:", latitude, ":", start_date, "-", end_date, sep=" ");
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
        } else if (life_stage=="Egg") {
            title = paste(insect, ": Egg Pop :", location, ": Lat:", latitude, ":", start_date, "-", end_date, sep=" ");
            legend_text = c("Egg");
            columns = c(4); 
            plot(days, group, main=title, type="l", ylim=c(0, maxval), axes=F, lwd=2, xlab="", ylab="", cex=3, cex.lab=3, cex.axis=3, cex.main=3);
            legend("topleft", legend_text, lty=c(1, 1, 1), col=columns, cex=3);
            axis(1, at=c(1:length(date_labels)) * 30 - 15, cex.axis=3, labels=date_labels);
            axis(2, cex.axis=3);
            if (plot_std_error=="yes") {
                # Standard error for group.
                lines(days, group+group_std_error, lty=2);
                lines(days, group-group_std_error, lty=2);
            }
        } else if (life_stage=="Nymph") {
            stage = paste(life_stages_nymph, "Nymph Pop", sep=" ");
            title = paste(insect, ": ", stage, location, ": Lat:", latitude, ":", start_date, "-", end_date, sep=" ");
            legend_text = c("Nymph");
            columns = c(2);
            plot(days, group, main=title, type="l", ylim=c(0, maxval), axes=F, lwd=2, xlab="", ylab="", cex=3, cex.lab=3, cex.axis=3, cex.main=3);
            legend("topleft", legend_text, lty=c(1, 1, 1), col=columns, cex=3);
            axis(1, at=c(1:length(date_labels)) * 30 - 15, cex.axis=3, labels=date_labels);
            axis(2, cex.axis=3);
            if (plot_std_error=="yes") {
                # Standard error for group.
                lines(days, group+group_std_error, lty=2);
                lines(days, group-group_std_error, lty=2);
            }
        } else if (life_stage=="Adult") {
            stage = paste(life_stages_adult, "Adult Pop", sep=" ");
            title = paste(insect, ": ", stage, location, ": Lat:", latitude, ":", start_date, "-", end_date, sep=" ");
            legend_text = c("Adult");
            columns = c(1);
            plot(days, group, main=title, type="l", ylim=c(0, maxval), axes=F, lwd=2, xlab="", ylab="", cex=3, cex.lab=3, cex.axis=3, cex.main=3);
            legend("topleft", legend_text, lty=c(1, 1, 1), col=columns, cex=3);
            axis(1, at=c(1:length(date_labels)) * 30 - 15, cex.axis=3, labels=date_labels);
            axis(2, cex.axis=3);
            if (plot_std_error=="yes") {
                # Standard error for group.
                lines(days, group+group_std_error, lty=2);
                lines(days, group-group_std_error, lty=2);
            }
        }
    } else if (chart_type=="pop_size_by_generation") {
        if (life_stage=="All") {
            title = paste(insect, ": Total Pop by Generation :", location, ": Lat:", latitude, ":", start_date, "-", end_date, sep=" ");
        } else if (life_stage=="Egg") {
            cat("Not implemented!");
            quit(save="no", status=1);
        } else if (life_stage=="Nymph") {
            cat("Not implemented!");
            quit(save="no", status=1);
        } else if (life_stage=="Adult") {
            title = paste(insect, ": Adult Pop by Generation :", location, ": Lat:", latitude, ":", start_date, "-", end_date, sep=" ");
        }
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

# Read the temperature data into a data frame.
temperature_data_frame = parse_input_data(opt$input, opt$num_days);
output_dir = "output_dir";
# Split life_stages into a list of strings for plots.
life_stages_str = as.character(opt$life_stages);
life_stages = strsplit(life_stages_str, ",")[[1]];
# Get the date labels for plots.
date_labels = get_date_labels(temperature_data_frame, opt$num_days);
# All latitude values are the same, so get the value for plots from the first row.
latitude = temperature_data_frame$LATITUDE[1];
# Get the number of days for plots.
num_columns = dim(temperature_data_frame)[2];

# Initialize matrices.
Eggs.replications = matrix(rep(0, opt$num_days*opt$replications), ncol=opt$replications);
YoungNymphs.replications = matrix(rep(0, opt$num_days*opt$replications), ncol=opt$replications);
OldNymphs.replications = matrix(rep(0, opt$num_days*opt$replications), ncol=opt$replications);
Previtellogenic.replications = matrix(rep(0, opt$num_days*opt$replications), ncol=opt$replications);
Vitellogenic.replications = matrix(rep(0, opt$num_days*opt$replications), ncol=opt$replications);
Diapausing.replications = matrix(rep(0, opt$num_days*opt$replications), ncol=opt$replications);

newborn.replications = matrix(rep(0, opt$num_days*opt$replications), ncol=opt$replications);
adult.replications = matrix(rep(0, opt$num_days*opt$replications), ncol=opt$replications);
death.replications = matrix(rep(0, opt$num_days*opt$replications), ncol=opt$replications);

# P is Parental, or overwintered adults.
P.replications = matrix(rep(0, opt$num_days*opt$replications), ncol=opt$replications);
P_adults.replications = matrix(rep(0, opt$num_days*opt$replications), ncol=opt$replications);
# F1 is the first field-produced generation.
F1.replications = matrix(rep(0, opt$num_days*opt$replications), ncol=opt$replications);
F1_adults.replications = matrix(rep(0, opt$num_days*opt$replications), ncol=opt$replications);
# F2 is the second field-produced generation.
F2.replications = matrix(rep(0, opt$num_days*opt$replications), ncol=opt$replications);
F2_adults.replications = matrix(rep(0, opt$num_days*opt$replications), ncol=opt$replications);

population.replications = matrix(rep(0, opt$num_days*opt$replications), ncol=opt$replications);

# Process replications.
for (N.replications in 1:opt$replications) {
    # Start with the user-defined number of insects per replication.
    num_insects = opt$insects_per_replication;
    # Generation, Stage, degree-days, T, Diapause.
    vector.ini = c(0, 3, 0, 0, 0);
    # Overwintering, previttelogenic, degree-days=0, T=0, no-diapause.
    vector.matrix = rep(vector.ini, num_insects);
    # Complete matrix for the population.
    vector.matrix = base::t(matrix(vector.matrix, nrow=5));
    # Time series of population size.
    Eggs = rep(0, opt$num_days);
    YoungNymphs = rep(0, opt$num_days);
    OldNymphs = rep(0, opt$num_days);
    Previtellogenic = rep(0, opt$num_days);
    Vitellogenic = rep(0, opt$num_days);
    Diapausing = rep(0, opt$num_days);

    N.newborn = rep(0, opt$num_days);
    N.adult = rep(0, opt$num_days);
    N.death = rep(0, opt$num_days);

    overwintering_adult.population = rep(0, opt$num_days);
    first_generation.population = rep(0, opt$num_days);
    second_generation.population = rep(0, opt$num_days);

    # P is Parental, or overwintered adults.
    P.adult = rep(0, opt$num_days);
    # F1 is the first field-produced generation.
    F1.adult = rep(0, opt$num_days);
    # F2 is the second field-produced generation.
    F2.adult = rep(0, opt$num_days);

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
                    # Overwintering adult (previttelogenic).
                    if (photoperiod > opt$photoperiod && vector.individual[3] > 68 && doy < 180) {
                        # Add 68C to become fully reproductively matured.
                        # Transfer to vittelogenic.
                        vector.individual = c(0, 4, 0, 0, 0);
                        vector.matrix[i,] = vector.individual;
                    }
                    else {
                        # Add to # Add average temperature for current day.
                        vector.individual[3] = vector.individual[3] + averages.temp;
                        # Add 1 day in current stage.
                        vector.individual[4] = vector.individual[4] + 1;
                        vector.matrix[i,] = vector.individual;
                    }
                }
                if (vector.individual[1] != 0 && vector.individual[2] == 3) {
                    # Not overwintering adult (previttelogenic).
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
                # Old nymph to adult: previttelogenic or diapausing?
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

        # Aggregate results by day.
        # Egg population size.
        Eggs[row] = sum(vector.matrix[,2]==0);
        # Young nymph population size.
        YoungNymphs[row] = sum(vector.matrix[,2]==1);
        # Old nymph population size.
        OldNymphs[row] = sum(vector.matrix[,2]==2);
        # Previtellogenic population size.
        Previtellogenic[row] = sum(vector.matrix[,2]==3);
        # Vitellogenic population size.
        Vitellogenic[row] = sum(vector.matrix[,2]==4);
        # Diapausing population size.
        Diapausing[row] = sum(vector.matrix[,2]==5);

        # Newborn population size.
        N.newborn[row] = num_insects.newborn;
        # Adult population size.
        N.adult[row] = sum(vector.matrix[,2]==3) + sum(vector.matrix[,2]==4) + sum(vector.matrix[,2]==5);
        # Dead population size.
        N.death[row] = num_insects.death;

        total.population = c(total.population, num_insects);

        # Overwintering adult population size.
        overwintering_adult.population[row] = sum(vector.matrix[,1]==0);
        # First generation population size.
        first_generation.population[row] = sum(vector.matrix[,1]==1);
        # Second generation population size.
        second_generation.population[row] = sum(vector.matrix[,1]==2);

        # P is Parental, or overwintered adults.
        P.adult[row] = sum(vector.matrix[,1]==0);
        # F1 is the first field-produced generation.
        F1.adult[row] = sum((vector.matrix[,1]==1 & vector.matrix[,2]==3) | (vector.matrix[,1]==1 & vector.matrix[,2]==4) | (vector.matrix[,1]==1 & vector.matrix[,2]==5));
        # F2 adult population size
        F2.adult[row] = sum((vector.matrix[,1]==2 & vector.matrix[,2]==3) | (vector.matrix[,1]==2 & vector.matrix[,2]==4) | (vector.matrix[,1]==2 & vector.matrix[,2]==5));
    }   # End of days specified in the input temperature data.

    averages.cum = cumsum(averages.day);

    # Define the output values.
    Eggs.replications[,N.replications] = Eggs;
    YoungNymphs.replications[,N.replications] = YoungNymphs;
    OldNymphs.replications[,N.replications] = OldNymphs;
    Previtellogenic.replications[,N.replications] = Previtellogenic;
    Vitellogenic.replications[,N.replications] = Vitellogenic;
    Diapausing.replications[,N.replications] = Diapausing;

    newborn.replications[,N.replications] = N.newborn;
    adult.replications[,N.replications] = N.adult;
    death.replications[,N.replications] = N.death;

    # P is Parental, or overwintered adults.
    P.replications[,N.replications] = overwintering_adult.population;
    P_adults.replications[,N.replications] = P.adult;
    # F1 is the first field-produced generation.
    F1.replications[,N.replications] = first_generation.population;
    F1_adults.replications[,N.replications] = F1.adult;
    # F2 is the second field-produced generation.
    F2.replications[,N.replications] = second_generation.population;
    F2_adults.replications[,N.replications] = F2.adult;

    population.replications[,N.replications] = total.population;
}

# Mean value for eggs.
eggs = apply(Eggs.replications, 1, mean);
# Standard error for eggs.
eggs.std_error = apply(Eggs.replications, 1, sd) / sqrt(opt$replications);

# Calculate nymph populations for selected life stage.
if (opt$life_stages_nymph=="All") {
    # Mean value for all nymphs.
    nymphs = apply((YoungNymphs.replications+OldNymphs.replications), 1, mean);
    # Standard error for all nymphs.
    nymphs.std_error = apply((YoungNymphs.replications+OldNymphs.replications) / sqrt(opt$replications), 1, sd);
} else if (opt$life_stages_nymph=="Young") {
    # Mean value for young nymphs.
    nymphs = apply(YoungNymphs.replications, 1, mean);
    # Standard error for young nymphs.
    nymphs.std_error = apply(YoungNymphs.replications / sqrt(opt$replications), 1, sd);
} else if (opt$life_stages_nymph=="Old") {
    # Mean value for old nymphs.
    nymphs = apply(OldNymphs.replications, 1, mean);
    # Standard error for old nymphs.
    nymphs.std_error = apply(OldNymphs.replications / sqrt(opt$replications), 1, sd);
}

# Calculate adult populations for selected life stage.
if (opt$life_stages_adult=="All") {
    # Mean value for all adults.
    adults = apply((Previtellogenic.replications+Vitellogenic.replications+Diapausing.replications), 1, mean);
    # Standard error for all adults.
    adults.std_error = apply((Previtellogenic.replications+Vitellogenic.replications+Diapausing.replications), 1, sd) / sqrt(opt$replications);
} else if (opt$life_stages_adult == "Pre-vittelogenic") {
    # Mean value for previtellogenic adults.
    adults = apply(Previtellogenic.replications, 1, mean);
    # Standard error for previtellogenic adults.
    adults.std_error = apply(Previtellogenic.replications, 1, sd) / sqrt(opt$replications);
} else if (opt$life_stages_adult == "Vittelogenic") {
    # Mean value for vitellogenic adults.
    adults = apply(Vitellogenic.replications, 1, mean);
    # Standard error for vitellogenic adults.
    adults.std_error = apply(Vitellogenic.replications, 1, sd) / sqrt(opt$replications);
} else if (opt$life_stages_adult == "Diapausing") {
    # Mean value for vitellogenic adults.
    adults = apply(Diapausing.replications, 1, mean);
    # Standard error for vitellogenic adults.
    adults.std_error = apply(Diapausing.replications, 1, sd) / sqrt(opt$replications);
}

# Mean value for P which is Parental, or overwintered adults.
P = apply(P.replications, 1, mean);
# Standard error for P.
P.std_error = apply(P.replications, 1, sd) / sqrt(opt$replications);

# Mean value for P adults.
P_adults = apply(P_adults.replications, 1, mean);
# Standard error for P_adult.
P_adults.std_error = apply(P_adults.replications, 1, sd) / sqrt(opt$replications);

# Mean value for F1, which is the first field-produced generation.
F1 = apply(F1.replications, 1, mean);
# Standard error for F1.
F1.std_error = apply(F1.replications, 1, sd) / sqrt(opt$replications);

# Mean value for F1 adults.
F1_adults = apply(F1_adults.replications, 1, mean);
# Standard error for F1 adult.
F1_adults.std_error = apply(F1_adults.replications, 1, sd) / sqrt(opt$replications);

# Mean value for F2, which is the second field-produced generation.
F2 = apply(F2.replications, 1, mean);
# Standard error for F2.
F2.std_error = apply(F2.replications, 1, sd) / sqrt(opt$replications);

# Mean value for F2 adults.
F2_adults = apply(F2_adults.replications, 1, mean);
# Standard error for F2 adult.
F2_adults.std_error = apply(F2_adults.replications, 1, sd) / sqrt(opt$replications);

# Display the total number of days in the Galaxy history item blurb.
cat("Number of days: ", opt$num_days, "\n");

# Information needed for plots plots.
days = c(1:opt$num_days);
start_date = temperature_data_frame$DATE[1];
end_date = temperature_data_frame$DATE[opt$num_days];

cat("\n plot_generations_separately: ", opt$plot_generations_separately, "\n");

if (opt$plot_generations_separately=="yes") {
    for (life_stage in life_stages) {}
        if (life_stage == "Egg") {
            cat("Not implemented!");
            quit(save="no", status=1);
        } else if (life_stage == "Nymph") {
            cat("Not implemented!");
            quit(save="no", status=1);
        } else if (life_stage == "Adult") {
            # Start PDF device driver.
            dev.new(width=20, height=30);
            file_path = paste(output_dir, "adult_pop_size_by_generation.pdf", sep="/");
            pdf(file=file_path, width=20, height=30, bg="white");
            par(mar=c(5, 6, 4, 4), mfrow=c(3, 1));
            # Adult population size by generation.
            maxval = max(F2_adults) + 100;
            render_chart(date_labels, "pop_size_by_generation", opt$plot_std_error, opt$insect, opt$location, latitude, start_date, end_date, days, maxval,
                life_stage, group=P_adults, group_std_error=P_adults.std_error, group2=F1_adults, group2_std_error=F1_adults.std_error,
                group3=F2_adults, group3_std_error=F2_adults.std_error);
            # Turn off device driver to flush output.
            dev.off();
        } else if (life_stage == "All") {
            # Start PDF device driver.
            dev.new(width=20, height=30);
            file_path = paste(output_dir, "total_pop_size_by_generation.pdf", sep="/");
            pdf(file=file_path, width=20, height=30, bg="white");
            par(mar=c(5, 6, 4, 4), mfrow=c(3, 1));
            # Total population size by generation.
            maxval = max(F2);
            render_chart(date_labels, "pop_size_by_generation", opt$plot_std_error, opt$insect, opt$location, latitude, start_date, end_date, days, maxval,
                life_stage, group=P, group_std_error=P.std_error, group2=F1, group2_std_error=F1.std_error, group3=F2, group3_std_error=F2.std_error);
            # Turn off device driver to flush output.
            dev.off();
        }
} else {
    for (life_stage in life_stages) {}
        if (life_stage == "Egg") {
            # Start PDF device driver.
            dev.new(width=20, height=30);
            file_path = paste(output_dir, "egg_pop_size.pdf", sep="/");
            pdf(file=file_path, width=20, height=30, bg="white");
            par(mar=c(5, 6, 4, 4), mfrow=c(3, 1));
            # Egg population size.
            maxval = max(eggs+eggs.std_error);
            render_chart(date_labels, "pop_size_by_life_stage", opt$plot_std_error, opt$insect, opt$location, latitude, start_date, end_date, days, maxval,
                life_stage, group=eggs, group_std_error=eggs.std_error);
            # Turn off device driver to flush output.
            dev.off();
        } else if (life_stage == "Nymph") {
            # Start PDF device driver.
            dev.new(width=20, height=30);
            file_name = paste(opt$life_stages_nymph, "nymph_pop_size.pdf", sep="_");
            file_path = paste(output_dir, file_name, sep="/");
            pdf(file=file_path, width=20, height=30, bg="white");
            par(mar=c(5, 6, 4, 4), mfrow=c(3, 1));
            # Nymph population size.
            maxval = max(nymphs+nymphs.std_error);
            render_chart(date_labels, "pop_size_by_life_stage", opt$plot_std_error, opt$insect, opt$location, latitude, start_date, end_date, days, maxval,
                life_stage, group=nymphs, group_std_error=nymphs.std_error, life_stages_nymph=opt$life_stages_nymph);
            # Turn off device driver to flush output.
            dev.off();
        } else if (life_stage == "Adult") {
            # Start PDF device driver.
            dev.new(width=20, height=30);
            file_name = paste(opt$life_stages_adult, "adult_pop_size.pdf", sep="_");
            file_path = paste(output_dir, file_name, sep="/");
            pdf(file=file_path, width=20, height=30, bg="white");
            par(mar=c(5, 6, 4, 4), mfrow=c(3, 1));
            # Adult population size.
            maxval = max(adults+adults.std_error);
            render_chart(date_labels, "pop_size_by_life_stage", opt$plot_std_error, opt$insect, opt$location, latitude, start_date, end_date, days, maxval,
                life_stage, group=adults, group_std_error=adults.std_error, life_stages_adult=opt$life_stages_adult);
            # Turn off device driver to flush output.
            dev.off();
        } else if (life_stage == "All") {
            # Start PDF device driver.
            dev.new(width=20, height=30);
            file_path = paste(output_dir, "total_pop_size.pdf", sep="/");
            pdf(file=file_path, width=20, height=30, bg="white");
            par(mar=c(5, 6, 4, 4), mfrow=c(3, 1));
            # Total population size.
            maxval = max(eggs+eggs.std_error, nymphs+nymphs.std_error, adults+adults.std_error);
            render_chart(date_labels, "pop_size_by_life_stage", opt$plot_std_error, opt$insect, opt$location, latitude, start_date, end_date, days, maxval,
                life_stage, group=adults, group_std_error=adults.std_error, group2=nymphs, group2_std_error=nymphs.std_error, group3=eggs,
                group3_std_error=eggs.std_error);
            # Turn off device driver to flush output.
            dev.off();
        }
}
