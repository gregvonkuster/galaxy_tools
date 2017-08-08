#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option(c("-a", "--adult_mort"), action="store", dest="adult_mort", type="integer", help="Adjustment rate for adult mortality"),
    make_option(c("-b", "--adult_accum"), action="store", dest="adult_accum", type="integer", help="Adjustment of DD accumulation (old nymph->adult)"),
    make_option(c("-c", "--egg_mort"), action="store", dest="egg_mort", type="integer", help="Adjustment rate for egg mortality"),
    make_option(c("-d", "--latitude"), action="store", dest="latitude", type="double", help="Latitude of selected location"),
    make_option(c("-e", "--location"), action="store", dest="location", help="Selected location"),
    make_option(c("-f", "--min_clutch_size"), action="store", dest="min_clutch_size", type="integer", help="Adjustment of minimum clutch size"),
    make_option(c("-i", "--max_clutch_size"), action="store", dest="max_clutch_size", type="integer", help="Adjustment of maximum clutch size"),
    make_option(c("-j", "--nymph_mort"), action="store", dest="nymph_mort", type="integer", help="Adjustment rate for nymph mortality"),
    make_option(c("-k", "--old_nymph_accum"), action="store", dest="old_nymph_accum", type="integer", help="Adjustment of DD accumulation (young nymph->old nymph)"),
    make_option(c("-o", "--output"), action="store", dest="output", help="Output dataset"),
    make_option(c("-p", "--oviposition"), action="store", dest="oviposition", type="integer", help="Adjustment for oviposition rate"),
    make_option(c("-q", "--photoperiod"), action="store", dest="photoperiod", type="double", help="Critical photoperiod for diapause induction/termination"),
    make_option(c("-s", "--replications"), action="store", dest="replications", type="integer", help="Number of replications"),
    make_option(c("-t", "--se_plot"), action="store", dest="se_plot", help="Plot SE"),
    make_option(c("-u", "--year"), action="store", dest="year", type="integer", help="Starting year"),
    make_option(c("-v", "--temperature_dataset"), action="store", dest="temperature_dataset", help="Temperature data for selected location"),
    make_option(c("-y", "--young_nymph_accum"), action="store", dest="young_nymph_accum", type="integer", help="Adjustment of DD accumulation (egg->young nymph)")
)

parser <- OptionParser(usage="%prog [options] file", option_list=option_list)
args <- parse_args(parser, positional_arguments=TRUE)
opt <- args$options

data.input=function(loc, year, temperature.dataset)
{
    expdata <- matrix(rep(0, 365 * 3), nrow=365)
    namedat <- paste(loc,  year, ".Rdat", sep="")
    temp.data <- read.csv(file=temperature.dataset, header=T)

    expdata[,1] <- c(1:365)
    # Minimum
    expdata[,2] <- temp.data[c(1:365), 3]
    # Maximum
    expdata[,3] <- temp.data[c(1:365), 2]
    save(expdata, file=namedat)
    namedat
}

daylength=function(latitude)
{
    # from Forsythe 1995
    p=0.8333
    dl <- NULL
    for (i in 1:365) {
        theta <- 0.2163108 + 2 * atan(0.9671396 * tan(0.00860 * (i - 186)))
        phi <- asin(0.39795 * cos(theta))
        dl[i] <- 24 - 24 / pi * acos((sin(p * pi / 180) + sin(latitude * pi / 180) * sin(phi)) / (cos(latitude * pi / 180) * cos(phi)))
    }
    dl   # return a vector of daylength in 365 days
}

hourtemp=function(latitude, date, temperature_file_path)
{
    load(temperature_file_path)
    threshold <- 14.17  # base development threshold for BMSB
    dnp <- expdata[date, 2]  # daily minimum
    dxp <- expdata[date, 3]  # daily maximum
    dmean <- 0.5 * (dnp + dxp)
    dd <- 0  # initialize degree day accumulation

    if (dxp<threshold) {
        dd <- 0
    }
    else {
        dlprofile <- daylength(latitude)  # extract daylength data for entire year
        T <- NULL  # initialize hourly temperature
        dh <- NULL #initialize degree hour vector
        # date <- 200
        y <- dlprofile[date]  # calculate daylength in given date
        z <- 24 - y     # night length
        a <- 1.86     # lag coefficient
        b <- 2.20     # night coefficient
        #tempdata <- read.csv("tempdata.csv") #import raw data set
        # Should be outside function otherwise its redundant
        risetime <- 12 - y / 2      # sunrise time
        settime <- 12 + y / 2       # sunset time
        ts <- (dxp - dnp) * sin(pi * (settime - 5) / (y + 2 * a)) + dnp
        for (i in 1:24) {
            if (i > risetime && i<settime) {
                m <- i - 5  # number of hours after Tmin until sunset
                T[i]=(dxp - dnp) * sin(pi * m / (y + 2 * a)) + dnp
                if (T[i]<8.4) {
                    dh[i] <- 0
                }
                else {
                    dh[i] <- T[i] - 8.4
                }
            }
            else if (i > settime) { 
                n <- i - settime
                T[i]=dnp + (ts - dnp) * exp( - b * n / z)
                if (T[i]<8.4) {
                    dh[i] <- 0
                }
                else {
                    dh[i] <- T[i] - 8.4
                }
            }
            else {
                n <- i + 24 - settime
                T[i]=dnp + (ts - dnp) * exp( - b * n / z)
                if (T[i]<8.4) {
                    dh[i] <- 0
                }
                else {
                    dh[i] <- T[i] - 8.4
                }
            }
        }
        dd <- sum(dh) / 24
    }
    return=c(dmean, dd)
    return
}

dev.egg = function(temperature)
{
    dev.rate= -0.9843 * temperature + 33.438
    return = dev.rate
    return
}

dev.young = function(temperature)
{
    n12 <- -0.3728 * temperature + 14.68
    n23 <- -0.6119 * temperature + 25.249
    dev.rate = mean(n12 + n23)
    return = dev.rate
    return
}

dev.old = function(temperature)
{
    n34 <- -0.6119 * temperature + 17.602
    n45 <- -0.4408 * temperature + 19.036
    dev.rate = mean(n34 + n45)
    return = dev.rate
    return
}

dev.emerg = function(temperature)
{
    emerg.rate <- -0.5332 * temperature + 24.147
    return = emerg.rate
    return
}

mortality.egg = function(temperature)
{
    if (temperature < 12.7) {
        mort.prob = 0.8
    }
    else {
        mort.prob = 0.8 - temperature / 40.0
        if (mort.prob < 0) {
            mort.prob = 0.01
        }
    }
    return = mort.prob
    return
}

mortality.nymph = function(temperature)
{
    if (temperature < 12.7) {
        mort.prob = 0.03
    }
    else {
        mort.prob = temperature * 0.0008 + 0.03
    }
    return = mort.prob
    return
}

mortality.adult = function(temperature)
{
    if (temperature < 12.7) {
        mort.prob = 0.002
    }
    else {
        mort.prob = temperature * 0.0005 + 0.02
    }
    return = mort.prob
    return
}

cat("Replications: ", opt$replications, "\n")
cat("Photoperiod: ", opt$photoperiod, "\n")
cat("Oviposition rate: ", opt$oviposition, "\n")
cat("Egg mortality rate: ", opt$egg_mort, "\n")
cat("Nymph mortality rate: ", opt$nymph_mort, "\n")
cat("Adult mortality rate: ", opt$adult_mort, "\n")
cat("Min clutch size: ", opt$min_clutch_size, "\n")
cat("Max clutch size: ", opt$max_clutch_size, "\n")
cat("(egg->young nymph): ", opt$young_nymph_accum, "\n")
cat("(young nymph->old nymph): ", opt$old_nymph_accum, "\n")
cat("(old nymph->adult): ", opt$adult_accum, "\n")

# Read in the input temperature datafile
temperature_file_path <- data.input(opt$location, opt$year, opt$temperature_dataset)

# Initialize matrix for results from all replications
S0.rep <- S1.rep <- S2.rep <- S3.rep <- S4.rep <- S5.rep <- matrix(rep(0, 365 * opt$replications), ncol = opt$replications)
newborn.rep <- death.rep <- adult.rep <- pop.rep <- g0.rep <- g1.rep <- g2.rep <- g0a.rep <- g1a.rep <- g2a.rep <- matrix(rep(0, 365 * opt$replications), ncol=opt$replications)

# loop through replications
for (N.rep in 1:opt$replications) {
    # during each replication
    # start with 1000 individuals -- user definable as well?
    n <- 1000
    # Generation, Stage, DD, T, Diapause
    vec.ini <- c(0, 3, 0, 0, 0)
    # overwintering, previttelogenic, DD=0, T=0, no-diapause
    vec.mat <- rep(vec.ini, n)
    # complete matrix for the population
    vec.mat <- t(matrix(vec.mat, nrow=5))
    # complete photoperiod profile in a year, requires daylength function
    ph.p <- daylength(opt$latitude)

    # time series of population size
    tot.pop <- NULL
    # gen.0 pop size
    gen0.pop <- rep(0, 365)
    gen1.pop <- rep(0, 365)
    gen2.pop <- rep(0, 365)
    S0 <- S1 <- S2 <- S3 <- S4 <- S5 <- rep(0, 365)
    g0.adult <- g1.adult <- g2.adult <- rep(0, 365)
    N.newborn <- N.death <- N.adult <- rep(0, 365)
    dd.day <- rep(0, 365)

    # start tick
    ptm <- proc.time()

    # all the days
    for (day in 1:365) {
        # photoperiod in the day
        photoperiod <- ph.p[day]
        temp.profile <- hourtemp(opt$latitude, day, temperature_file_path)
        mean.temp <- temp.profile[1]
        dd.temp <- temp.profile[2]
        dd.day[day] <- dd.temp
        # trash bin for death
        death.vec <- NULL
        # new born
        birth.vec <- NULL

        # all individuals
        for (i in 1:n) {
            # find individual record
            vec.ind <- vec.mat[i,]
            # first of all, still alive?  
            # adjustment for late season mortality rate
            if (opt$latitude < 40.0) {
                post.mort <- 1
                day.kill <- 300
            }
            else {
                post.mort <- 2
                day.kill <- 250
            }
            if (vec.ind[2] == 0) {
                # egg
                death.prob = opt$egg_mort * mortality.egg(mean.temp)
            }
            else if (vec.ind[2] == 1 | vec.ind[2] == 2) {
                death.prob = opt$nymph_mort * mortality.nymph(mean.temp)
            }
            else if (vec.ind[2] == 3 | vec.ind[2] == 4 | vec.ind[2] == 5) {
                # for adult
                if (day < day.kill) {
                    death.prob = opt$adult_mort * mortality.adult(mean.temp)
                }
                else {
                    # increase adult mortality after fall equinox
                    death.prob = opt$adult_mort * post.mort * mortality.adult(mean.temp)
                }
            }
            # (or dependent on temperature and life stage?)
            u.d <- runif(1)
            if (u.d < death.prob) {
                death.vec <- c(death.vec, i)
            } 
            else {
                # aggregrate index of dead bug
                # event 1 end of diapause
                if (vec.ind[1] == 0 && vec.ind[2] == 3) {
                    # overwintering adult (previttelogenic)
                    if (photoperiod > opt$photoperiod && vec.ind[3] > 68 && day < 180) {
                        # add 68C to become fully reproductively matured
                        # transfer to vittelogenic
                        vec.ind <- c(0, 4, 0, 0, 0)
                        vec.mat[i,] <- vec.ind
                    }
                    else {
                        # add to DD
                        vec.ind[3] <- vec.ind[3] + dd.temp
                        # add 1 day in current stage
                        vec.ind[4] <- vec.ind[4] + 1
                        vec.mat[i,] <- vec.ind
                    }
                }
                if (vec.ind[1] != 0 && vec.ind[2] == 3) {
                    # NOT overwintering adult (previttelogenic)
                    current.gen <- vec.ind[1]
                    if (vec.ind[3] > 68) {
                        # add 68C to become fully reproductively matured
                        # transfer to vittelogenic
                        vec.ind <- c(current.gen, 4, 0, 0, 0)
                        vec.mat[i,] <- vec.ind
                    }
                    else {
                        # add to DD
                        vec.ind[3] <- vec.ind[3] + dd.temp
                        # add 1 day in current stage
                        vec.ind[4] <- vec.ind[4] + 1
                        vec.mat[i,] <- vec.ind
                    }
                }

                # event 2 oviposition -- where population dynamics comes from
                if (vec.ind[2] == 4 && vec.ind[1] == 0 && mean.temp > 10) {
                    # vittelogenic stage, overwintering generation
                    if (vec.ind[4] == 0) {
                        # just turned in vittelogenic stage
                        n.birth=round(runif(1, 2 + opt$min_clutch_size, 8 + opt$max_clutch_size))
                    }
                    else {
                        # daily probability of birth
                        p.birth = opt$oviposition * 0.01
                        u1 <- runif(1)
                        if (u1 < p.birth) {
                            n.birth=round(runif(1, 2, 8))
                        }
                    }
                    # add to DD
                    vec.ind[3] <- vec.ind[3] + dd.temp
                    # add 1 day in current stage
                    vec.ind[4] <- vec.ind[4] + 1
                    vec.mat[i,] <- vec.ind
                    if (n.birth > 0) {
                        # add new birth -- might be in different generations
                        # generation + 1
                        new.gen <- vec.ind[1] + 1
                        # egg profile
                        new.ind <- c(new.gen, 0, 0, 0, 0)
                        new.vec <- rep(new.ind, n.birth)
                        # update batch of egg profile
                        new.vec <- t(matrix(new.vec, nrow=5))
                        # group with total eggs laid in that day
                        birth.vec <- rbind(birth.vec, new.vec)
                    }
                }

                # event 2 oviposition -- for gen 1.
                if (vec.ind[2] == 4 && vec.ind[1] == 1 && mean.temp > 12.5 && day < 222) {
                    # vittelogenic stage, 1st generation
                    if (vec.ind[4] == 0) {
                        # just turned in vittelogenic stage
                        n.birth=round(runif(1, 2 + opt$min_clutch_size, 8 + opt$max_clutch_size))
                    }
                    else {
                        # daily probability of birth
                        p.birth = opt$oviposition * 0.01
                        u1 <- runif(1)
                        if (u1 < p.birth) {
                            n.birth = round(runif(1, 2, 8))
                        }
                    }
                    # add to DD
                    vec.ind[3] <- vec.ind[3] + dd.temp
                    # add 1 day in current stage
                    vec.ind[4] <- vec.ind[4] + 1
                    vec.mat[i,] <- vec.ind
                    if (n.birth > 0) {
                        # add new birth -- might be in different generations
                        # generation + 1
                        new.gen <- vec.ind[1] + 1
                        # egg profile
                        new.ind <- c(new.gen, 0, 0, 0, 0)
                        new.vec <- rep(new.ind, n.birth)
                        # update batch of egg profile
                        new.vec <- t(matrix(new.vec, nrow=5))
                        # group with total eggs laid in that day
                        birth.vec <- rbind(birth.vec, new.vec)
                    }
                }

                # event 3 development (with diapause determination)
                # event 3.1 egg development to young nymph (vec.ind[2]=0 -> egg)
                if (vec.ind[2] == 0) {
                    # egg stage
                    # add to DD
                    vec.ind[3] <- vec.ind[3] + dd.temp
                    if (vec.ind[3] >= (68 + opt$young_nymph_accum)) {
                        # from egg to young nymph, DD requirement met
                        current.gen <- vec.ind[1]
                        # transfer to young nym stage
                        vec.ind <- c(current.gen, 1, 0, 0, 0)
                    }
                    else {
                        # add 1 day in current stage
                        vec.ind[4] <- vec.ind[4] + 1
                    }
                    vec.mat[i,] <- vec.ind
                }

                # event 3.2 young nymph to old nymph (vec.ind[2]=1 -> young nymph: determines diapause)
                if (vec.ind[2] == 1) {
                    # young nymph stage
                    # add to DD
                    vec.ind[3] <- vec.ind[3] + dd.temp
                    if (vec.ind[3] >= (250 + opt$old_nymph_accum)) {
                        # from young to old nymph, DD requirement met
                        current.gen <- vec.ind[1]
                        # transfer to old nym stage
                        vec.ind <- c(current.gen, 2, 0, 0, 0)
                        if (photoperiod < opt$photoperiod && day > 180) {
                            vec.ind[5] <- 1
                        } # prepare for diapausing
                    }
                    else {
                        # add 1 day in current stage
                        vec.ind[4] <- vec.ind[4] + 1
                    }
                    vec.mat[i,] <- vec.ind
                }  

                # event 3.3 old nymph to adult: previttelogenic or diapausing?
                if (vec.ind[2] == 2) {
                    # old nymph stage
                    # add to DD
                    vec.ind[3] <- vec.ind[3] + dd.temp
                    if (vec.ind[3] >= (200 + opt$adult_accum)) {
                        # from old to adult, DD requirement met
                        current.gen <- vec.ind[1]
                        if (vec.ind[5] == 0) {
                            # non-diapausing adult -- previttelogenic
                            vec.ind <- c(current.gen, 3, 0, 0, 0)
                        }
                        else {
                            # diapausing 
                            vec.ind <- c(current.gen, 5, 0, 0, 1)
                        }
                    }
                    else {
                        # add 1 day in current stage
                        vec.ind[4] <- vec.ind[4] + 1
                    }
                    vec.mat[i,] <- vec.ind
                }

                # event 4 growing of diapausing adult (unimportant, but still necessary)## 
                if (vec.ind[2] == 5) {
                    vec.ind[3] <- vec.ind[3] + dd.temp
                    vec.ind[4] <- vec.ind[4] + 1
                    vec.mat[i,] <- vec.ind
                }
            } # else if it is still alive
        } # end of the individual bug loop

        # find how many died
        n.death <- length(death.vec)
        if (n.death > 0) {
            vec.mat <- vec.mat[-death.vec, ]
        }
        # remove record of dead
        # find how many new born  
        n.newborn <- length(birth.vec[,1])
        vec.mat <- rbind(vec.mat, birth.vec)
        # update population size for the next day
        n <- n - n.death + n.newborn 

        # aggregate results by day
        tot.pop <- c(tot.pop, n) 
        # egg
        s0 <- sum(vec.mat[,2] == 0)
        # young nymph
        s1 <- sum(vec.mat[,2] == 1)
        # old nymph
        s2 <- sum(vec.mat[,2] == 2)
        # previtellogenic
        s3 <- sum(vec.mat[,2] == 3)
        # vitellogenic
        s4 <- sum(vec.mat[,2] == 4)
        # diapausing
        s5 <- sum(vec.mat[,2] == 5)
        # overwintering adult
        gen0 <- sum(vec.mat[,1] == 0)
        # first generation
        gen1 <- sum(vec.mat[,1] == 1)
        # second generation
        gen2 <- sum(vec.mat[,1] == 2)
        # sum of all adults
        n.adult <- sum(vec.mat[,2] == 3) + sum(vec.mat[,2] == 4) + sum(vec.mat[,2] == 5)
        # gen.0 pop size
        gen0.pop[day] <- gen0
        gen1.pop[day] <- gen1
        gen2.pop[day] <- gen2
        S0[day] <- s0
        S1[day] <- s1
        S2[day] <- s2
        S3[day] <- s3
        S4[day] <- s4
        S5[day] <- s5
        g0.adult[day] <- sum(vec.mat[,1] == 0)
        g1.adult[day] <- sum((vec.mat[,1] == 1 & vec.mat[,2] == 3) | (vec.mat[,1] == 1 & vec.mat[,2] == 4) | (vec.mat[,1] == 1 & vec.mat[,2] == 5))
        g2.adult[day] <- sum((vec.mat[,1]== 2 & vec.mat[,2] == 3) | (vec.mat[,1] == 2 & vec.mat[,2] == 4) | (vec.mat[,1] == 2 & vec.mat[,2] == 5))

        N.newborn[day] <- n.newborn
        N.death[day] <- n.death
        N.adult[day] <- n.adult
        #print(c(N.rep, day, n, n.adult))
    }   # end of 365 days

    dd.cum <- cumsum(dd.day)
    # collect all the outputs
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

# save(dd.day, dd.cum, S0.rep, S1.rep, S2.rep, S3.rep, S4.rep, S5.rep, newborn.rep, death.rep, adult.rep, pop.rep, g0.rep, g1.rep, g2.rep, g0a.rep, g1a.rep, g2a.rep, file=opt$output)
# maybe do not need to export this bit, but for now just leave it as-is
# do we need to export this Rdat file? 

# Data analysis and visualization
# default: plot 1 year of result
# but can be expanded to accommodate multiple years
n.yr <- 1
day.all <- c(1:365 * n.yr)

# mean value for adults
sa <- apply((S3.rep + S4.rep + S5.rep), 1, mean)
# mean value for nymphs
sn <- apply((S1.rep + S2.rep), 1,mean)
# mean value for eggs
se <- apply(S0.rep, 1, mean)
# mean value for P
g0 <- apply(g0.rep, 1, mean)
# mean value for F1
g1 <- apply(g1.rep, 1, mean)
# mean value for F2
g2 <- apply(g2.rep, 1, mean)
# mean value for P adult
g0a <- apply(g0a.rep, 1, mean)
# mean value for F1 adult
g1a <- apply(g1a.rep, 1, mean)
# mean value for F2 adult
g2a <- apply(g2a.rep, 1, mean)

# SE for adults
sa.se <- apply((S3.rep + S4.rep + S5.rep), 1, sd) / sqrt(opt$replications)
# SE for nymphs
sn.se <- apply((S1.rep + S2.rep) / sqrt(opt$replications), 1, sd)
# SE for eggs
se.se <- apply(S0.rep, 1, sd) / sqrt(opt$replications)
# SE value for P
g0.se <- apply(g0.rep, 1, sd) / sqrt(opt$replications)
# SE for F1
g1.se <- apply(g1.rep, 1, sd) / sqrt(opt$replications)
# SE for F2
g2.se <- apply(g2.rep, 1, sd) / sqrt(opt$replications)
# SE for P adult
g0a.se <- apply(g0a.rep, 1, sd) / sqrt(opt$replications)
# SE for F1 adult
g1a.se <- apply(g1a.rep, 1, sd) / sqrt(opt$replications)
# SE for F2 adult
g2a.se <- apply(g2a.rep, 1, sd) / sqrt(opt$replications)

dev.new(width=20, height=20)

# Start PDF device driver to save charts to output.
pdf(file=opt$output, height=20, width=20, bg="white")

par(mar = c(5, 6, 4, 4), mfrow=c(3, 1))

# Subfigure 2: population size by life stage
plot(day.all, sa, main = "BSMB Total Population Size by Life Stage", type = "l", ylim = c(0, max(se + se.se, sn + sn.se, sa + sa.se)), axes = F, lwd = 2, xlab = "", ylab = "Number", cex = 2, cex.lab = 2, cex.axis = 2, cex.main = 2)
# Young and old nymphs
lines(day.all, sn, lwd = 2, lty = 1, col = 2)
# Eggs
lines(day.all, se, lwd = 2, lty = 1, col = 4)
axis(1, at = c(1:12) * 30 - 15, cex.axis = 2, labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
axis(2, cex.axis = 2)
leg.text <- c("Egg", "Nymph", "Adult")
legend("topleft", leg.text, lty = c(1, 1, 1), col = c(4, 2, 1), cex = 2)
if (opt$se_plot == 1) {
    # add SE lines to plot
    # SE for adults
    lines (day.all, sa + sa.se, lty = 2)
    lines (day.all, sa - sa.se, lty = 2) 
    # SE for nymphs
    lines (day.all, sn + sn.se, col = 2, lty = 2)
    lines (day.all, sn - sn.se, col = 2, lty = 2) 
    # SE for eggs
    lines (day.all, se + se.se, col = 4, lty = 2)
    lines (day.all, se - se.se, col = 4, lty = 2) 
}

# Subfigure 3: population size by generation
plot(day.all, g0, main = "BSMB Total Population Size by Generation", type = "l", ylim = c(0, max(g2)), axes = F, lwd = 2, xlab = "", ylab = "Number", cex = 2, cex.lab = 2, cex.axis = 2, cex.main = 2)
lines(day.all, g1, lwd = 2, lty = 1, col = 2)
lines(day.all, g2, lwd = 2, lty = 1, col = 4)
axis(1, at = c(1:12) * 30 - 15, cex.axis = 2, labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
axis(2, cex.axis = 2)
leg.text <- c("P", "F1", "F2")
legend("topleft", leg.text, lty = c(1, 1, 1), col =c(1, 2, 4), cex = 2)
if (opt$se_plot == 1) {
    # add SE lines to plot
    # SE for adults
    lines (day.all, g0 + g0.se, lty = 2)
    lines (day.all, g0 - g0.se, lty = 2) 
    # SE for nymphs
    lines (day.all, g1 + g1.se, col = 2, lty = 2)
    lines (day.all, g1 - g1.se, col = 2, lty = 2) 
    # SE for eggs
    lines (day.all, g2 + g2.se, col = 4, lty = 2)
    lines (day.all, g2 - g2.se, col = 4, lty = 2) 
}

# Subfigure 4: adult population size by generation
plot(day.all, g0a, ylim = c(0, max(g2a) + 100), main = "BSMB Adult Population Size by Generation", type = "l", axes = F, lwd = 2, xlab = "Year", ylab = "Number", cex = 2, cex.lab = 2, cex.axis = 2, cex.main = 2)
lines(day.all, g1a, lwd = 2, lty = 1, col = 2)
lines(day.all, g2a, lwd = 2, lty = 1, col = 4)
axis(1, at = c(1:12) * 30 - 15, cex.axis = 2, labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
axis(2, cex.axis = 2)
leg.text <- c("P", "F1", "F2")
legend("topleft", leg.text, lty = c(1, 1, 1), col = c(1, 2, 4), cex = 2)
if (opt$se_plot == 1) {
    # add SE lines to plot
    # SE for adults
    lines (day.all, g0a + g0a.se, lty = 2)
    lines (day.all, g0a - g0a.se, lty = 2) 
    # SE for nymphs
    lines (day.all, g1a + g1a.se, col = 2, lty = 2)
    lines (day.all, g1a - g1a.se, col = 2, lty = 2) 
    # SE for eggs
    lines (day.all, g2a + g2a.se, col = 4, lty = 2)
    lines (day.all, g2a - g2a.se, col = 4, lty = 2) 
}

# Turn off device driver to flush output.
dev.off()
