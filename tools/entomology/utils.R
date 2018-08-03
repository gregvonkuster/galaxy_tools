#!/usr/bin/env Rscript

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
    file_path = paste("output_plots_dir", file_name, sep="/");
    return(file_path);
}

get_year_from_date = function(date_str) {
    date_str_items = strsplit(date_str, "-")[[1]];
    return (date_str_items[1]);
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

get_tick_index = function(index, last_tick, ticks, tick_labels, tick_sep) {
    # The R code tries hard not to draw overlapping tick labels, and so
    # will omit labels where they would abut or overlap previously drawn
    # labels. This can result in, for example, every other tick being
    # labelled.  We'll keep track of the last tick to make sure all of
    # the month labels are displayed, and missing ticks are restricted
    # to Sundays which have no labels anyway.
    if (last_tick==0) {
        return(length(ticks)+1);
    }
    last_saved_tick = ticks[[length(ticks)]];
    if (index-last_saved_tick<tick_sep) {
        last_saved_month = tick_labels[[length(tick_labels)]];
        if (last_saved_month=="") {
            # We're safe overwriting a tick
            # with no label (i.e., a Sunday tick).
            return(length(ticks));
        } else {
            # Don't eliminate a Month label.
            return(NULL);
        }
    }
    return(length(ticks)+1);
}

get_total_days = function(is_leap_year) {
    # Get the total number of days in the current year.
    if (is_leap_year) {
        return(366);
    } else {
        return(365);
    }
}

get_x_axis_ticks_and_labels = function(temperature_data_frame, prepend_end_doy_norm=0, append_start_doy_norm=0, date_interval=FALSE) {
    # Generate a list of ticks and labels for plotting the x axis.
    if (prepend_end_doy_norm > 0) {
        prepend_end_norm_row = which(temperature_data_frame$DOY==prepend_end_doy_norm);
    } else {
        prepend_end_norm_row = 0;
    }
    if (append_start_doy_norm > 0) {
        append_start_norm_row = which(temperature_data_frame$DOY==append_start_doy_norm);
    } else {
        append_start_norm_row = 0;
    }
    num_rows = dim(temperature_data_frame)[1];
    tick_labels = list();
    ticks = list();
    current_month_label = NULL;
    last_tick = 0;
    if (date_interval) {
        tick_sep = 0;
    } else {
        tick_sep = 3;
    }
    for (i in 1:num_rows) {
        # Get the year and month from the date which
        # has the format YYYY-MM-DD.
        date = format(temperature_data_frame$DATE[i]);
        # Get the month label.
        items = strsplit(date, "-")[[1]];
        month = items[2];
        month_label = month.abb[as.integer(month)];
        day = as.integer(items[3]);
        doy = as.integer(temperature_data_frame$DOY[i]);
        # We're plotting the entire year, so ticks will
        # occur on Sundays and the first of each month.
        if (i == prepend_end_norm_row) {
            # Add a tick for the end of the 30 year normnals data
            # that was prepended to the year-to-date data.
            label_str = "End prepended 30 year normals";
            tick_index = get_tick_index(i, last_tick, ticks, tick_labels, tick_sep)
            ticks[tick_index] = i;
            if (date_interval) {
                # Append the day to label_str
                tick_labels[tick_index] = paste(label_str, day, sep=" ");
            } else {
                tick_labels[tick_index] = label_str;
            }
            last_tick = i;
        } else if (doy == append_start_doy_norm) {
            # Add a tick for the start of the 30 year normnals data
            # that was appended to the year-to-date data.
            label_str = "Start appended 30 year normals";
            tick_index = get_tick_index(i, last_tick, ticks, tick_labels, tick_sep)
            ticks[tick_index] = i;
            if (!identical(current_month_label, month_label)) {
                # Append the month to label_str.
                label_str = paste(label_str, month_label, spe=" ");
                current_month_label = month_label;
            }
            if (date_interval) {
                # Append the day to label_str
                label_str = paste(label_str, day, sep=" ");
            }
            tick_labels[tick_index] = label_str;
            last_tick = i;
        } else if (i==num_rows) {
            # Add a tick for the last day of the year.
            label_str = "";
            tick_index = get_tick_index(i, last_tick, ticks, tick_labels, tick_sep)
            ticks[tick_index] = i;
            if (!identical(current_month_label, month_label)) {
                # Append the month to label_str.
                label_str = month_label;
                current_month_label = month_label;
            }
            if (date_interval) {
                # Append the day to label_str
                label_str = paste(label_str, day, sep=" ");
            }
            tick_labels[tick_index] = label_str;
        } else {
            if (!identical(current_month_label, month_label)) {
                # Add a tick for the month.
                tick_index = get_tick_index(i, last_tick, ticks, tick_labels, tick_sep)
                ticks[tick_index] = i;
                if (date_interval) {
                    # Append the day to the month.
                    tick_labels[tick_index] = paste(month_label, day, sep=" ");
                } else {
                    tick_labels[tick_index] = month_label;
                }
                current_month_label = month_label;
                last_tick = i;
            }
            tick_index = get_tick_index(i, last_tick, ticks, tick_labels, tick_sep)
            if (!is.null(tick_index)) {
                if (date_interval) {
                    # Add a tick for every day. The first tick is the
                    # month label, so add a tick only if i is not 1
                    if (i>1 & day>1) {
                        tick_index = get_tick_index(i, last_tick, ticks, tick_labels, tick_sep)
                        ticks[tick_index] = i;
                        # Add the day as the label.
                        tick_labels[tick_index] = day;
                        last_tick = i;
                    }
                } else {
                    # Get the day.
                    day = weekdays(as.Date(date));
                    if (day=="Sunday") {
                        # Add a tick if we're on a Sunday.
                        ticks[tick_index] = i;
                        # Add a blank month label so it is not displayed.
                        tick_labels[tick_index] = "";
                        last_tick = i;
                    }
                }
            }
        }
    }
    return(list(ticks, tick_labels));
}

render_chart = function(ticks, date_labels, chart_type, plot_std_error, insect, location, latitude, start_date, end_date, days, maxval,
        replications, life_stage, group, group_std_error, group2=NULL, group2_std_error=NULL, group3=NULL, group3_std_error=NULL,
        life_stages_adult=NULL, life_stages_nymph=NULL) {
    if (chart_type=="pop_size_by_life_stage") {
        if (life_stage=="Total") {
            title = paste(insect, ": Reps", replications, ":", life_stage, "Pop :", location, ": Lat", latitude, ":", start_date, "-", end_date, sep=" ");
            legend_text = c("Egg", "Nymph", "Adult");
            columns = c(4, 2, 1);
            plot(days, group, main=title, type="l", ylim=c(0, maxval), axes=FALSE, lwd=2, xlab="", ylab="", cex=3, cex.lab=3, cex.axis=3, cex.main=3);
            legend("topleft", legend_text, lty=c(1, 1, 1), col=columns, cex=3);
            lines(days, group2, lwd=2, lty=1, col=2);
            lines(days, group3, lwd=2, lty=1, col=4);
            axis(side=1, at=ticks, labels=date_labels, las=2, font.axis=3, xpd=TRUE, cex=3, cex.lab=3, cex.axis=3, cex.main=3);
            axis(side=2, font.axis=3, xpd=TRUE, cex=3, cex.lab=3, cex.axis=3, cex.main=3);
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
            plot(days, group, main=title, type="l", ylim=c(0, maxval), axes=FALSE, lwd=2, xlab="", ylab="", cex=3, cex.lab=3, cex.axis=3, cex.main=3);
            legend("topleft", legend_text, lty=c(1), col="black", cex=3);
            axis(side=1, at=ticks, labels=date_labels, las=2, font.axis=3, xpd=TRUE, cex=3, cex.lab=3, cex.axis=3, cex.main=3);
            axis(side=2, font.axis=3, xpd=TRUE, cex=3, cex.lab=3, cex.axis=3, cex.main=3);
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
        plot(days, group, main=title, type="l", ylim=c(0, maxval), axes=FALSE, lwd=2, xlab="", ylab="", cex=3, cex.lab=3, cex.axis=3, cex.main=3);
        legend("topleft", legend_text, lty=c(1, 1, 1), col=columns, cex=3);
        lines(days, group2, lwd=2, lty=1, col=2);
        lines(days, group3, lwd=2, lty=1, col=4);
        axis(side=1, at=ticks, labels=date_labels, las=2, font.axis=3, xpd=TRUE, cex=3, cex.lab=3, cex.axis=3, cex.main=3);
        axis(side=2, font.axis=3, xpd=TRUE, cex=3, cex.lab=3, cex.axis=3, cex.main=3);
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

stop_err = function(msg) {
    cat(msg, file=stderr());
    quit(save="no", status=1);
}

validate_date = function(date_str) {
    valid_date = as.Date(date_str, format="%Y-%m-%d");
    if( class(valid_date)=="try-error" || is.na(valid_date)) {
        msg = paste("Invalid date: ", date_str, ", valid date format is yyyy-mm-dd.", sep="");
        stop_err(msg);
    }
    return(valid_date);
}