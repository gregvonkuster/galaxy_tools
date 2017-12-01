#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option(c("-a", "--adult_mort"), action="store", dest="adult_mort", type="integer", help="Adjustment rate for adult mortality"),
    make_option(c("-b", "--adult_accum"), action="store", dest="adult_accum", type="integer", help="Adjustment of degree-days accumulation (old nymph->adult)"),
    make_option(c("-c", "--egg_mort"), action="store", dest="egg_mort", type="integer", help="Adjustment rate for egg mortality"),
    make_option(c("-e", "--location"), action="store", dest="location", help="Selected location"),
    make_option(c("-f", "--min_clutch_size"), action="store", dest="min_clutch_size", type="integer", help="Adjustment of minimum clutch size"),
    make_option(c("-i", "--max_clutch_size"), action="store", dest="max_clutch_size", type="integer", help="Adjustment of maximum clutch size"),
    make_option(c("-j", "--nymph_mort"), action="store", dest="nymph_mort", type="integer", help="Adjustment rate for nymph mortality"),
    make_option(c("-k", "--old_nymph_accum"), action="store", dest="old_nymph_accum", type="integer", help="Adjustment of degree-days accumulation (young nymph->old nymph)"),
    make_option(c("-n", "--num_days"), action="store", dest="num_days", type="integer", help="Total number of days in the temperature dataset"),
    make_option(c("-o", "--output"), action="store", dest="output", help="Output dataset"),
    make_option(c("-p", "--oviposition"), action="store", dest="oviposition", type="integer", help="Adjustment for oviposition rate"),
    make_option(c("-q", "--photoperiod"), action="store", dest="photoperiod", type="double", help="Critical photoperiod for diapause induction/termination"),
    make_option(c("-s", "--replications"), action="store", dest="replications", type="integer", help="Number of replications"),
    make_option(c("-t", "--std_error_plot"), action="store", dest="std_error_plot", help="Plot Standard error"),
    make_option(c("-v", "--input"), action="store", dest="input", help="Temperature data for selected location"),
    make_option(c("-y", "--young_nymph_accum"), action="store", dest="young_nymph_accum", type="integer", help="Adjustment of degree-days accumulation (egg->young nymph)"),
    make_option(c("-x", "--insect"), action="store", dest="insect", help="Insect name")
)

parser <- OptionParser(usage="%prog [options] file", option_list=option_list)
args <- parse_args(parser, positional_arguments=TRUE)
opt <- args$options

parse_input_data = function(input_file, num_rows) {
    # Read in the input temperature datafile into a data frame.
    temperature_data_frame <- read.csv(file=input_file, header=T, strip.white=TRUE, sep=",")
    num_columns <- dim(temperature_data_frame)[2]
    if (num_columns == 6) {
        # The input data has the following 6 columns:
        # LATITUDE, LONGITUDE, DATE, DOY, TMIN, TMAX
        # Set the column names for access when adding daylight length..
        colnames(temperature_data_frame) <- c("LATITUDE","LONGITUDE", "DATE", "DOY", "TMIN", "TMAX")
        # Add a column containing the daylight length for each day.
        temperature_data_frame <- add_daylight_length(temperature_data_frame, num_columns, num_rows)
        # Reset the column names with the additional column for later access.
        colnames(temperature_data_frame) <- c("LATITUDE","LONGITUDE", "DATE", "DOY", "TMIN", "TMAX", "DAYLEN")
    }
    return(temperature_data_frame)
}

add_daylight_length = function(temperature_data_frame, num_columns, num_rows) {
    # Return a vector of daylight length (photoperido profile) for
    # the number of days specified in the input temperature data
    # (from Forsythe 1995).
    p = 0.8333
    latitude <- temperature_data_frame$LATITUDE[1]
    daylight_length_vector <- NULL
    for (i in 1:num_rows) {
        # Get the day of the year from the current row
        # of the temperature data for computation.
        doy <- temperature_data_frame$DOY[i]
        theta <- 0.2163108 + 2 * atan(0.9671396 * tan(0.00860 * (doy - 186)))
        phi <- asin(0.39795 * cos(theta))
        # Compute the length of daylight for the day of the year.
        daylight_length_vector[i] <- 24 - (24 / pi * acos((sin(p * pi / 180) + sin(latitude * pi / 180) * sin(phi)) / (cos(latitude * pi / 180) * cos(phi))))
    }
    # Append daylight_length_vector as a new column to temperature_data_frame.
    temperature_data_frame[, num_columns+1] <- daylight_length_vector
    return(temperature_data_frame)
}

get_temperature_at_hour = function(latitude, temperature_data_frame, row, num_days) {
    # Base development threshold for Brown Marmolated Stink Bug
    # insect phenology model.
    # TODO: Pass insect on the command line to accomodate more
    # the just the Brown Marmolated Stink Bub.
    threshold <- 14.17
    # Minimum temperature for current row.
    dnp <- temperature_data_frame$TMIN[row]
    # Maximum temperature for current row.
    dxp <- temperature_data_frame$TMAX[row]
    # Mean temperature for current row.
    dmean <- 0.5 * (dnp + dxp)
    # Initialize degree day accumulation
    degree_days <- 0
    if (dxp < threshold) {
        degree_days <- 0
    }
    else {
        # Initialize hourly temperature.
        T <- NULL
        # Initialize degree hour vector.
        dh <- NULL
        # Daylight length for current row.
        y <- temperature_data_frame$DAYLEN[row]
        # Darkness length.
        z <- 24 - y
        # Lag coefficient.
        a <- 1.86
        # Darkness coefficient.
        b <- 2.20
        # Sunrise time.
        risetime <- 12 - y / 2
        # Sunset time.
        settime <- 12 + y / 2
        ts <- (dxp - dnp) * sin(pi * (settime - 5) / (y + 2 * a)) + dnp
        for (i in 1:24) {
            if (i > risetime && i < settime) {
                # Number of hours after Tmin until sunset.
                m <- i - 5
                T[i] = (dxp - dnp) * sin(pi * m / (y + 2 * a)) + dnp
                if (T[i] < 8.4) {
                    dh[i] <- 0
                }
                else {
                    dh[i] <- T[i] - 8.4
                }
            }
            else if (i > settime) { 
                n <- i - settime
                T[i] = dnp + (ts - dnp) * exp( - b * n / z)
                if (T[i] < 8.4) {
                    dh[i] <- 0
                }
                else {
                    dh[i] <- T[i] - 8.4
                }
            }
            else {
                n <- i + 24 - settime
                T[i]=dnp + (ts - dnp) * exp( - b * n / z)
                if (T[i] < 8.4) {
                    dh[i] <- 0
                }
                else {
                    dh[i] <- T[i] - 8.4
                }
            }
        }
        degree_days <- sum(dh) / 24
    }
    return(c(dmean, degree_days))
}

dev.egg = function(temperature) {
    dev.rate= -0.9843 * temperature + 33.438
    return(dev.rate)
}

dev.young = function(temperature) {
    n12 <- -0.3728 * temperature + 14.68
    n23 <- -0.6119 * temperature + 25.249
    dev.rate = mean(n12 + n23)
    return(dev.rate)
}

dev.old = function(temperature) {
    n34 <- -0.6119 * temperature + 17.602
    n45 <- -0.4408 * temperature + 19.036
    dev.rate = mean(n34 + n45)
    return(dev.rate)
}

dev.emerg = function(temperature) {
    emerg.rate <- -0.5332 * temperature + 24.147
    return(emerg.rate)
}

mortality.egg = function(temperature) {
    if (temperature < 12.7) {
        mort.prob = 0.8
    }
    else {
        mort.prob = 0.8 - temperature / 40.0
        if (mort.prob < 0) {
            mort.prob = 0.01
        }
    }
    return(mort.prob)
}

mortality.nymph = function(temperature) {
    if (temperature < 12.7) {
        mort.prob = 0.03
    }
    else {
        mort.prob = temperature * 0.0008 + 0.03
    }
    return(mort.prob)
}

mortality.adult = function(temperature) {
    if (temperature < 12.7) {
        mort.prob = 0.002
    }
    else {
        mort.prob = temperature * 0.0005 + 0.02
    }
    return(mort.prob)
}

temperature_data_frame <- parse_input_data(opt$input, opt$num_days)
# All latitude values are the same,
# so get the value from the first row.
latitude <- temperature_data_frame$LATITUDE[1]

cat("Number of days: ", opt$num_days, "\n")

# Initialize matrix for results from all replications.
S0.rep <- S1.rep <- S2.rep <- S3.rep <- S4.rep <- S5.rep <- matrix(rep(0, opt$num_days * opt$replications), ncol = opt$replications)
newborn.rep <- death.rep <- adult.rep <- pop.rep <- g0.rep <- g1.rep <- g2.rep <- g0a.rep <- g1a.rep <- g2a.rep <- matrix(rep(0, opt$num_days * opt$replications), ncol=opt$replications)

# Loop through replications.
for (N.rep in 1:opt$replications) {
    # During each replication start with 1000 individuals.
    # TODO: user definable as well?
    num_insects <- 1000
    # Generation, Stage, degree-days, T, Diapause.
    vec.ini <- c(0, 3, 0, 0, 0)
    # Overwintering, previttelogenic, degree-days=0, T=0, no-diapause.
    vec.mat <- rep(vec.ini, num_insects)
    # Complete matrix for the population.
    vec.mat <- base::t(matrix(vec.mat, nrow=5))
    # Time series of population size.
    tot.pop <- NULL
    gen0.pop <- rep(0, opt$num_days)
    gen1.pop <- rep(0, opt$num_days)
    gen2.pop <- rep(0, opt$num_days)
    S0 <- S1 <- S2 <- S3 <- S4 <- S5 <- rep(0, opt$num_days)
    g0.adult <- g1.adult <- g2.adult <- rep(0, opt$num_days)
    N.newborn <- N.death <- N.adult <- rep(0, opt$num_days)
    degree_days.day <- rep(0, opt$num_days)
    # All the days included in the input temperature dataset.
    for (row in 1:opt$num_days) {
        # Get the integer day of the year for the current row.
        doy <- temperature_data_frame$DOY[row]
        # Photoperiod in the day.
        photoperiod <- temperature_data_frame$DAYLEN[row]
        temp.profile <- get_temperature_at_hour(latitude, temperature_data_frame, row, opt$num_days)
        mean.temp <- temp.profile[1]
        degree_days.temp <- temp.profile[2]
        degree_days.day[row] <- degree_days.temp
        # Trash bin for death.
        death.vec <- NULL
        # Newborn.
        birth.vec <- NULL
        # All individuals.
        for (i in 1:num_insects) {
            # Find individual record.
            vec.ind <- vec.mat[i,]
            # First of all, still alive?
            # Adjustment for late season mortality rate.
            if (latitude < 40.0) {
                post.mort <- 1
                day.kill <- 300
            }
            else {
                post.mort <- 2
                day.kill <- 250
            }
            if (vec.ind[2] == 0) {
                # Egg.
                death.prob = opt$egg_mort * mortality.egg(mean.temp)
            }
            else if (vec.ind[2] == 1 | vec.ind[2] == 2) {
                death.prob = opt$nymph_mort * mortality.nymph(mean.temp)
            }
            else if (vec.ind[2] == 3 | vec.ind[2] == 4 | vec.ind[2] == 5) {
                # For adult.
                if (doy < day.kill) {
                    death.prob = opt$adult_mort * mortality.adult(mean.temp)
                }
                else {
                    # Increase adult mortality after fall equinox.
                    death.prob = opt$adult_mort * post.mort * mortality.adult(mean.temp)
                }
            }
            # (or dependent on temperature and life stage?)
            u.d <- runif(1)
            if (u.d < death.prob) {
                death.vec <- c(death.vec, i)
            } 
            else {
                # Aggregrate index of dead bug.
                # Event 1 end of diapause.
                if (vec.ind[1] == 0 && vec.ind[2] == 3) {
                    # Overwintering adult (previttelogenic).
                    if (photoperiod > opt$photoperiod && vec.ind[3] > 68 && doy < 180) {
                        # Add 68C to become fully reproductively matured.
                        # Transfer to vittelogenic.
                        vec.ind <- c(0, 4, 0, 0, 0)
                        vec.mat[i,] <- vec.ind
                    }
                    else {
                        # Add to degree_days.
                        vec.ind[3] <- vec.ind[3] + degree_days.temp
                        # Add 1 day in current stage.
                        vec.ind[4] <- vec.ind[4] + 1
                        vec.mat[i,] <- vec.ind
                    }
                }
                if (vec.ind[1] != 0 && vec.ind[2] == 3) {
                    # Not overwintering adult (previttelogenic).
                    current.gen <- vec.ind[1]
                    if (vec.ind[3] > 68) {
                        # Add 68C to become fully reproductively matured.
                        # Transfer to vittelogenic.
                        vec.ind <- c(current.gen, 4, 0, 0, 0)
                        vec.mat[i,] <- vec.ind
                    }
                    else {
                        # Add to degree_days.
                        vec.ind[3] <- vec.ind[3] + degree_days.temp
                        # Add 1 day in current stage.
                        vec.ind[4] <- vec.ind[4] + 1
                        vec.mat[i,] <- vec.ind
                    }
                }
                # Event 2 oviposition -- where population dynamics comes from.
                if (vec.ind[2] == 4 && vec.ind[1] == 0 && mean.temp > 10) {
                    # Vittelogenic stage, overwintering generation.
                    if (vec.ind[4] == 0) {
                        # Just turned in vittelogenic stage.
                        num_insects.birth = round(runif(1, 2 + opt$min_clutch_size, 8 + opt$max_clutch_size))
                    }
                    else {
                        # Daily probability of birth.
                        p.birth = opt$oviposition * 0.01
                        u1 <- runif(1)
                        if (u1 < p.birth) {
                            num_insects.birth = round(runif(1, 2, 8))
                        }
                    }
                    # Add to degree_days.
                    vec.ind[3] <- vec.ind[3] + degree_days.temp
                    # Add 1 day in current stage.
                    vec.ind[4] <- vec.ind[4] + 1
                    vec.mat[i,] <- vec.ind
                    if (num_insects.birth > 0) {
                        # Add new birth -- might be in different generations.
                        new.gen <- vec.ind[1] + 1
                        # Egg profile.
                        new.ind <- c(new.gen, 0, 0, 0, 0)
                        new.vec <- rep(new.ind, num_insects.birth)
                        # Update batch of egg profile.
                        new.vec <- t(matrix(new.vec, nrow=5))
                        # Group with total eggs laid in that day.
                        birth.vec <- rbind(birth.vec, new.vec)
                    }
                }
                # Event 2 oviposition -- for generation 1.
                if (vec.ind[2] == 4 && vec.ind[1] == 1 && mean.temp > 12.5 && doy < 222) {
                    # Vittelogenic stage, 1st generation
                    if (vec.ind[4] == 0) {
                        # Just turned in vittelogenic stage.
                        num_insects.birth=round(runif(1, 2 + opt$min_clutch_size, 8 + opt$max_clutch_size))
                    }
                    else {
                        # Daily probability of birth.
                        p.birth = opt$oviposition * 0.01
                        u1 <- runif(1)
                        if (u1 < p.birth) {
                            num_insects.birth = round(runif(1, 2, 8))
                        }
                    }
                    # Add to degree_days.
                    vec.ind[3] <- vec.ind[3] + degree_days.temp
                    # Add 1 day in current stage.
                    vec.ind[4] <- vec.ind[4] + 1
                    vec.mat[i,] <- vec.ind
                    if (num_insects.birth > 0) {
                        # Add new birth -- might be in different generations.
                        new.gen <- vec.ind[1] + 1
                        # Egg profile.
                        new.ind <- c(new.gen, 0, 0, 0, 0)
                        new.vec <- rep(new.ind, num_insects.birth)
                        # Update batch of egg profile.
                        new.vec <- t(matrix(new.vec, nrow=5))
                        # Group with total eggs laid in that day.
                        birth.vec <- rbind(birth.vec, new.vec)
                    }
                }
                # Event 3 development (with diapause determination).
                # Event 3.1 egg development to young nymph (vec.ind[2]=0 -> egg).
                if (vec.ind[2] == 0) {
                    # Egg stage.
                    # Add to degree_days.
                    vec.ind[3] <- vec.ind[3] + degree_days.temp
                    if (vec.ind[3] >= (68 + opt$young_nymph_accum)) {
                        # From egg to young nymph, degree-days requirement met.
                        current.gen <- vec.ind[1]
                        # Transfer to young nymph stage.
                        vec.ind <- c(current.gen, 1, 0, 0, 0)
                    }
                    else {
                        # Add 1 day in current stage.
                        vec.ind[4] <- vec.ind[4] + 1
                    }
                    vec.mat[i,] <- vec.ind
                }
                # Event 3.2 young nymph to old nymph (vec.ind[2]=1 -> young nymph: determines diapause).
                if (vec.ind[2] == 1) {
                    # Young nymph stage.
                    # Add to degree_days.
                    vec.ind[3] <- vec.ind[3] + degree_days.temp
                    if (vec.ind[3] >= (250 + opt$old_nymph_accum)) {
                        # From young to old nymph, degree_days requirement met.
                        current.gen <- vec.ind[1]
                        # Transfer to old nym stage.
                        vec.ind <- c(current.gen, 2, 0, 0, 0)
                        if (photoperiod < opt$photoperiod && doy > 180) {
                            vec.ind[5] <- 1
                        } # Prepare for diapausing.
                    }
                    else {
                        # Add 1 day in current stage.
                        vec.ind[4] <- vec.ind[4] + 1
                    }
                    vec.mat[i,] <- vec.ind
                }  
                # Event 3.3 old nymph to adult: previttelogenic or diapausing?
                if (vec.ind[2] == 2) {
                    # Old nymph stage.
                    # Add to degree_days.
                    vec.ind[3] <- vec.ind[3] + degree_days.temp
                    if (vec.ind[3] >= (200 + opt$adult_accum)) {
                        # From old to adult, degree_days requirement met.
                        current.gen <- vec.ind[1]
                        if (vec.ind[5] == 0) {
                            # Non-diapausing adult -- previttelogenic.
                            vec.ind <- c(current.gen, 3, 0, 0, 0)
                        }
                        else {
                            # Diapausing.
                            vec.ind <- c(current.gen, 5, 0, 0, 1)
                        }
                    }
                    else {
                        # Add 1 day in current stage.
                        vec.ind[4] <- vec.ind[4] + 1
                    }
                    vec.mat[i,] <- vec.ind
                }
                # Event 4 growing of diapausing adult (unimportant, but still necessary).
                if (vec.ind[2] == 5) {
                    vec.ind[3] <- vec.ind[3] + degree_days.temp
                    vec.ind[4] <- vec.ind[4] + 1
                    vec.mat[i,] <- vec.ind
                }
            } # Else if it is still alive.
        } # End of the individual bug loop.
        # Find how many died.
        num_insects.death <- length(death.vec)
        if (num_insects.death > 0) {
            vec.mat <- vec.mat[-death.vec, ]
        }
        # Remove record of dead.
        # Find how many new born.
        num_insects.newborn <- length(birth.vec[,1])
        vec.mat <- rbind(vec.mat, birth.vec)
        # Update population size for the next day.
        num_insects <- num_insects - num_insects.death + num_insects.newborn 

        # Aggregate results by day.
        tot.pop <- c(tot.pop, num_insects) 
        # Egg.
        s0 <- sum(vec.mat[,2] == 0)
        # Young nymph.
        s1 <- sum(vec.mat[,2] == 1)
        # Old nymph.
        s2 <- sum(vec.mat[,2] == 2)
        # Previtellogenic.
        s3 <- sum(vec.mat[,2] == 3)
        # Vitellogenic.
        s4 <- sum(vec.mat[,2] == 4)
        # Diapausing.
        s5 <- sum(vec.mat[,2] == 5)
        # Overwintering adult.
        gen0 <- sum(vec.mat[,1] == 0)
        # First generation.
        gen1 <- sum(vec.mat[,1] == 1)
        # Second generation.
        gen2 <- sum(vec.mat[,1] == 2)
        # Sum of all adults.
        num_insects.adult <- sum(vec.mat[,2] == 3) + sum(vec.mat[,2] == 4) + sum(vec.mat[,2] == 5)

        # Generation 0 pop size.
        gen0.pop[row] <- gen0
        gen1.pop[row] <- gen1
        gen2.pop[row] <- gen2

        S0[row] <- s0
        S1[row] <- s1
        S2[row] <- s2
        S3[row] <- s3
        S4[row] <- s4
        S5[row] <- s5

        g0.adult[row] <- sum(vec.mat[,1] == 0)
        g1.adult[row] <- sum((vec.mat[,1] == 1 & vec.mat[,2] == 3) | (vec.mat[,1] == 1 & vec.mat[,2] == 4) | (vec.mat[,1] == 1 & vec.mat[,2] == 5))
        g2.adult[row] <- sum((vec.mat[,1]== 2 & vec.mat[,2] == 3) | (vec.mat[,1] == 2 & vec.mat[,2] == 4) | (vec.mat[,1] == 2 & vec.mat[,2] == 5))

        N.newborn[row] <- num_insects.newborn
        N.death[row] <- num_insects.death
        N.adult[row] <- num_insects.adult
    }   # end of days specified in the input temperature data

    degree_days.cum <- cumsum(degree_days.day)

    # Collect all the outputs.
    S0.rep[,N.rep] <- S0
    S1.rep[,N.rep] <- S1
    S2.rep[,N.rep] <- S2
    S3.rep[,N.rep] <- S3
    S4.rep[,N.rep] <- S4
    S5.rep[,N.rep] <- S5
    newborn.rep[,N.rep] <- N.newborn
    death.rep[,N.rep] <- N.death
    adult.rep[,N.rep] <- N.adult
    pop.rep[,N.rep] <- tot.pop
    g0.rep[,N.rep] <- gen0.pop
    g1.rep[,N.rep] <- gen1.pop
    g2.rep[,N.rep] <- gen2.pop
    g0a.rep[,N.rep] <- g0.adult
    g1a.rep[,N.rep] <- g1.adult
    g2a.rep[,N.rep] <- g2.adult
}

# Mean value for adults
mean_value_adult <- apply((S3.rep + S4.rep + S5.rep), 1, mean)
# Mean value for nymphs
mean_value_nymph <- apply((S1.rep + S2.rep), 1, mean)
# Mean value for eggs
mean_value_egg <- apply(S0.rep, 1, mean)
# Mean value for P
g0 <- apply(g0.rep, 1, mean)
# Mean value for F1
g1 <- apply(g1.rep, 1, mean)
# Mean value for F2
g2 <- apply(g2.rep, 1, mean)
# Mean value for P adult
g0a <- apply(g0a.rep, 1, mean)
# Mean value for F1 adult
g1a <- apply(g1a.rep, 1, mean)
# Mean value for F2 adult
g2a <- apply(g2a.rep, 1, mean)

# Standard error for adults
mean_value_adult.std_error <- apply((S3.rep + S4.rep + S5.rep), 1, sd) / sqrt(opt$replications)
# Standard error for nymphs
mean_value_nymph.std_error <- apply((S1.rep + S2.rep) / sqrt(opt$replications), 1, sd)
# Standard error for eggs
mean_value_egg.std_error <- apply(S0.rep, 1, sd) / sqrt(opt$replications)
# Standard error value for P
g0.std_error <- apply(g0.rep, 1, sd) / sqrt(opt$replications)
# Standard error for F1
g1.std_error <- apply(g1.rep, 1, sd) / sqrt(opt$replications)
# Standard error for F2
g2.std_error <- apply(g2.rep, 1, sd) / sqrt(opt$replications)
# Standard error for P adult
g0a.std_error <- apply(g0a.rep, 1, sd) / sqrt(opt$replications)
# Standard error for F1 adult
g1a.std_error <- apply(g1a.rep, 1, sd) / sqrt(opt$replications)
# Standard error for F2 adult
g2a.std_error <- apply(g2a.rep, 1, sd) / sqrt(opt$replications)

dev.new(width=20, height=30)

# Start PDF device driver to save charts to output.
pdf(file=opt$output, width=20, height=30, bg="white")

par(mar=c(5, 6, 4, 4), mfrow=c(3, 1))

# Data analysis and visualization plots
# only within a single calendar year.
day.all <- c(1:opt$num_days)
start_date <- temperature_data_frame$DATE[1]
end_date <- temperature_data_frame$DATE[opt$num_days]

# Subfigure 1: population size by life stage
title <- paste(opt$insect, ": Total pop. by life stage :", opt$location, ": Lat:", latitude, ":", start_date, "-", end_date, sep=" ")
plot(day.all, mean_value_adult, main=title, type="l", ylim=c(0, max(mean_value_egg + mean_value_egg.std_error, mean_value_nymph + mean_value_nymph.std_error, mean_value_adult + mean_value_adult.std_error)), axes=F, lwd=2, xlab="", ylab="", cex=3, cex.lab=3, cex.axis=3, cex.main=3)
# Young and old nymphs.
lines(day.all, mean_value_nymph, lwd=2, lty=1, col=2)
# Eggs
lines(day.all, mean_value_egg, lwd=2, lty=1, col=4)
axis(1, at=c(1:12) * 30 - 15, cex.axis=3, labels=c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
axis(2, cex.axis=3)
leg.text <- c("Egg", "Nymph", "Adult")
legend("topleft", leg.text, lty=c(1, 1, 1), col=c(4, 2, 1), cex=3)
if (opt$se_plot == 1) {
    # Add Standard error lines to plot
    # Standard error for adults
    lines (day.all, mean_value_adult+mean_value_adult.std_error, lty=2)
    lines (day.all, mean_value_adult-mean_value_adult.std_error, lty=2) 
    # Standard error for nymphs
    lines (day.all, mean_value_nymph+mean_value_nymph.std_error, col=2, lty=2)
    lines (day.all, mean_value_nymph-mean_value_nymph.std_error, col=2, lty=2)
    # Standard error for eggs
    lines (day.all, mean_value_egg+mean_value_egg.std_error, col=4, lty=2)
    lines (day.all, mean_value_egg-mean_value_egg.std_error, col=4, lty=2)
}

# Subfigure 2: population size by generation
title <- paste(opt$insect, ": Total pop. by generation :", opt$location, ": Lat:", latitude, ":", start_date, "-", end_date, sep=" ")
plot(day.all, g0, main=title, type="l", ylim=c(0, max(g2)), axes=F, lwd=2, xlab="", ylab="", cex=3, cex.lab=3, cex.axis=3, cex.main=3)
lines(day.all, g1, lwd = 2, lty = 1, col=2)
lines(day.all, g2, lwd = 2, lty = 1, col=4)
axis(1, at=c(1:12) * 30 - 15, cex.axis=3, labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
axis(2, cex.axis=3)
leg.text <- c("P", "F1", "F2")
legend("topleft", leg.text, lty=c(1, 1, 1), col=c(1, 2, 4), cex=3)
if (opt$se_plot == 1) {
    # Add Standard error lines to plot
    # Standard error for adults
    lines (day.all, g0+g0.std_error, lty=2)
    lines (day.all, g0-g0.std_error, lty=2) 
    # Standard error for nymphs
    lines (day.all, g1+g1.std_error, col=2, lty=2)
    lines (day.all, g1-g1.std_error, col=2, lty=2)
    # Standard error for eggs
    lines (day.all, g2+g2.std_error, col=4, lty=2)
    lines (day.all, g2-g2.std_error, col=4, lty=2)
}

# Subfigure 3: adult population size by generation
title <- paste(opt$insect, ": Adult pop. by generation :", opt$location, ": Lat:", latitude, ":", start_date, "-", end_date, sep=" ")
plot(day.all, g0a, ylim=c(0, max(g2a) + 100), main=title, type="l", axes=F, lwd=2, xlab="", ylab="", cex=3, cex.lab=3, cex.axis=3, cex.main=3)
lines(day.all, g1a, lwd = 2, lty = 1, col=2)
lines(day.all, g2a, lwd = 2, lty = 1, col=4)
axis(1, at=c(1:12) * 30 - 15, cex.axis=3, labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
axis(2, cex.axis=3)
leg.text <- c("P", "F1", "F2")
legend("topleft", leg.text, lty=c(1, 1, 1), col=c(1, 2, 4), cex=3)
if (opt$std_error_plot == 1) {
    # Add Standard error lines to plot
    # Standard error for adults
    lines (day.all, g0a+g0a.std_error, lty=2)
    lines (day.all, g0a-g0a.std_error, lty=2) 
    # Standard error for nymphs
    lines (day.all, g1a+g1a.std_error, col=2, lty=2)
    lines (day.all, g1a-g1a.std_error, col=2, lty=2)
    # Standard error for eggs
    lines (day.all, g2a+g2a.std_error, col=4, lty=2)
    lines (day.all, g2a-g2a.std_error, col=4, lty=2)
}

# Turn off device driver to flush output.
dev.off()
