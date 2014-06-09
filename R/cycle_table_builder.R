script_version <- 9
# This script is based on calculations in microcosm_datalogger_visualisatoin_v8.Rmd

# expt : the numerical index of the experinent - this identifies the data files to use as input
# expt_name : the name of the experiment - which will define the output directory of the cycle tables

#######################################################################################
# these variables need to be set by the calling script, or uncommented to run directly:
# expt <- "044"
# expt_name = "30C"
# data_dir <- paste("C:\\sync\\lab_book\\kalahari\\files\\microcosms",sep="")
# max_cycle_number <- 200
# min_cycle_number <- 0
# verbose=TRUE
#######################################################################################

save_cycle_table <- TRUE

day_start <- as.POSIXlt("01/01/01-14:00","UTC","%d/%m/%Y-%H:%M")
day_end <- as.POSIXlt("01/01/01-07:00","UTC","%d/%m/%Y-%H:%M")

#day_start <- as.POSIXlt("01/01/01-06:00","UTC","%d/%m/%Y-%H:%M")
#day_end <- as.POSIXlt("01/01/01-18:00","UTC","%d/%m/%Y-%H:%M")
tz="UTC"
# adjustment for british summer time (error was made in setup of expt 36). Expts 36 on were run on BST
# No time adjustment should be needed in this script because the Gallenkamp chamber program was adjusted

# NOTE - this commented out when converting script for general use in repository "fluxan"
#if(as.numeric(expt)==36) {
#  day_start <- as.POSIXlt("01/01/01-07:00","UTC","%d/%m/%Y-%H:%M")
#  day_end <- as.POSIXlt("01/01/01-19:00","UTC","%d/%m/%Y-%H:%M")
#  tz="UTC+1"
#}
#if(as.numeric(expt)>36) {
#  # day_start <- as.POSIXlt("01/01/01-05:00","UTC","%d/%m/%Y-%H:%M")
#  # day_end <- as.POSIXlt("01/01/01-17:00","UTC","%d/%m/%Y-%H:%M")
#  tz="UTC+1"
#}

expt_dir <- paste(data_dir,"\\",expt,sep="")
# setwd(expt_dir)
cycle_table_output_dir <- paste(output_dir,"\\cycle_tables_",expt_name,"_v",script_version,sep="")
simple_cycle_table_dir <- paste(output_dir,"\\cycle_tables(simple)_",expt_name,"_v",script_version,sep="")

cycle_table_output_file <- paste(cycle_table_output_dir,"\\",expt,"-cycle_table.txt",sep="")
simple_cycle_table_output_file <- paste(simple_cycle_table_dir,"\\",expt,"-cycle_table.txt",sep="")

functions_dir <- "R"


# define data file locations
# infile <- paste(expt_dir,"\\","2012-12-14_10,18,24_datalogger.txt",sep="")
infiles <- Sys.glob(file.path(expt_dir, "*datalogger.txt"))
core_mass_file <- Sys.glob(file.path(expt_dir, "*mass.txt"))

# load libraries
library("reshape2")
source("R/CO2_flux_calc_v6.R")

# read in data files
d <- data.frame()
for(i in infiles) {
  if(verbose) { print(i)}
  x <- read.table(i,header=TRUE,stringsAsFactors=TRUE,na.strings=c("-1","-1.0","NA"))
  d <- rbind(d,x)
}

# Work out which core is connected to each chamber, and also make a table relating the temp probes
chamber_temp_probe_map <- as.data.frame(matrix( 
  c("0","t57240332","1","tC2B20342","2","tBCB203AB","3","tA8B2039D","4","t40C3032B"), 
  nrow=5, ncol=2,byrow = TRUE)  ,row.names=c("0","1","2","3","4") )  
colnames(chamber_temp_probe_map) <- c("chamber_t","probe")

core_chamber_map <- subset(unique(d[c("IRGA_conn","chamber")]),chamber!="NA")
row.names(core_chamber_map) <- core_chamber_map$chamber
colnames(core_chamber_map) <- c("core","chamber")

core_chamber_map <- cbind(core_chamber_map,chamber_temp_probe_map)[c("chamber","core","probe")]

cols <- colnames(d)
core_names <- levels(d$IRGA_conn)
core_names <- subset(core_names, (core_names!="ambient")&(core_names!="ref_gas")&(core_names!="ND"))

for (i in 1:5){
  find <- as.character(core_chamber_map$probe)[i]
  replace <- paste(as.character(core_chamber_map$core)[i],"_t",sep="")
  #replace <- i
  cols <- sub(find,replace,cols)
}

colnames(d) <- cols
temp_sensor_col_names <- paste(core_names,"_t",sep="")


if(length(core_mass_file)>0) {
  core_mass <- read.table(core_mass_file,header=TRUE,stringsAsFactors=FALSE,check.names=FALSE)
  core_mass$time <- as.POSIXct(row.names(core_mass),"UTC","%d/%m/%Y_%H:%M")
  # make a table of relative core masses.
  # i.e. subtract the first mass from each core
  core_mass_names <- names(subset(core_mass, select = -c(time)))
  core_mass_relative <- core_mass
  
  for(this in core_mass_names) {
    core_mass_relative[this] <- core_mass[this] - core_mass[1,this]
  }
  
  
  # add elapsed time in seconds to the main dataframe and core mass tables, to help time-based calculations
  core_mass$time_elapsed <- as.numeric(difftime(core_mass$time,core_mass$time[1],units="sec"))
  core_mass_relative$time_elapsed <- as.numeric(difftime(core_mass_relative$time,core_mass_relative$time[1],units="sec"))
  
  # add some convenience time columns:
  core_mass$time_elapsed_min <- core_mass$time_elapsed / 60
  
  # Find out when water was added:
  max_row <- which.max(rowSums(core_mass_relative[1:5]))[[1]]
  time_of_water_addition <- core_mass_relative[max_row,"time"]
}

# convert time string into time object for main dataframe and mass table
d$time <- as.POSIXct(d$time,"UTC","%d/%m/%Y-%H:%M:%S")

# sort the dataframe. Although it should be in order, when multiple data files are 
# combined they might not be cobmined in the right order.
d <- d[order(d$time),]

# add elapsed time in seconds to the main dataframe and core mass tables, to help time-based calculations
d$time_elapsed <- as.numeric(difftime(d$time,d$time[1],units="sec"))

# add some convenience time columns:
d$conn_timer_sec <- d$conn_timer / 1000
d$time_elapsed_min <- d$time_elapsed /60
d$time_elapsed_hour <- d$time_elapsed /3600

d$time_UTC_decimal_hour <- unclass(as.POSIXlt(d$time))$hour + (unclass(as.POSIXlt(d$time))$min/60) + (unclass(as.POSIXlt(d$time))$sec/3600)

# set our day start/end times to the first day of the data to help with calculations later
data_time <- min(d$time)
day_end$year <- as.POSIXlt(data_time)$year
day_end$mon <- as.POSIXlt(data_time)$mon
day_end$mday <- as.POSIXlt(data_time)$mday
day_start$year <- as.POSIXlt(data_time)$year
day_start$mon <- as.POSIXlt(data_time)$mon
day_start$mday <- as.POSIXlt(data_time)$mday

# work out the experiment elapsed time that corresponds to when the lights are switched off
day_end_time_elapsed <- min(d$time_elapsed_hour[d$time>day_end])

# make sure we don't look for non-existent cycles
max_cycle_number <- min(max_cycle_number,max(d$cycle,na.rm=TRUE))
min_cycle_number <- max(min_cycle_number,min(d$cycle,na.rm=TRUE))

# Only work with the specified cycles 
# (taking care to retain NA cycles during the same period because these are useful for environmental data ), and making sure we don't look for non-existent cycles
max_cycle_number <- min(max_cycle_number,max(d$cycle,na.rm=TRUE))
data_start_time <- min(d$time_elapsed[d$cycle==min_cycle_number],na.rm=TRUE)
data_end_time <-   max(d$time_elapsed[d$cycle==max_cycle_number],na.rm=TRUE)

# work out the ref gas concentration
ref_gas_conc <- min(d$IRGA_CO2[d$time_elapsed<data_start_time], na.rm=TRUE)

d <- d[(d$time_elapsed>=data_start_time)&(d$time_elapsed<=data_end_time),]


# output some summary info
if(verbose) {
cat("Script version: ",script_version,
    ".\nCore mass file: ", core_mass_file,
    ".\nSave cycle table: ", save_cycle_table,
    ".\nMax cycle: ", max_cycle_number,
    ".\nMin cycle: ",min_cycle_number)
cat("Start time: ",as.character(d$time[1]),
    "\nEnd time: ",as.character(d$time[nrow(d)]),
    "\nLights on:",day_start$hour,":",day_start$min,
    "\nLights off:",day_end$hour,":",day_end$min," (",day_end_time_elapsed," hours into experiment)",
    "\nRef gas conc:",ref_gas_conc,
    "\nTime zone: ",tz,
    "\nInfiles: ", infiles)
}
if(length(core_mass_file)>0) {
  if(verbose) { 
    print(core_mass_relative[,1:5])
  cat("water added at ",as.character(time_of_water_addition))
  }
}


if(length(core_mass_file)>0) {
  # core mass calculations, and moisture estimation
  
  d.m <- melt(core_mass_relative,measure.vars=core_mass_names)
  
  for(thisCore in core_names) {
    # get data for 1 core at a time, and only data after water was added
    d.m.this <- subset(d.m,(variable==thisCore)&value>0)
    
    mass.lm <- lm(formula=value~time_elapsed,data=d.m.this)
    c <- coef(mass.lm)[1]
    m <- coef(mass.lm)[2]
    calculated_mass <- (subset(d,IRGA_conn==thisCore)["time_elapsed"] * m) + c
    # find out which rows we need to update for this core
    row <- which(d$IRGA_conn==thisCore)
    d[row,"water_content_estimate"] <- (d$time_elapsed[row] * m) + c
  }
  
  # Now set all the times before the water addition to 0. 
  # find out which rows we need to update 
  row <- which(d$time<time_of_water_addition)
  # and set those to 0
  d[row,"water_content_estimate"] <- 0
} else {
  time_of_water_addition <- FALSE
  d$water_content_estimate <- -1
  # cycle table will not contain moisture data, so save it in a different folder
  cycle_table_output_file <- simple_cycle_table_output_file
  cycle_table_output_dir <- simple_cycle_table_dir
}
d.chambers <- subset(d,IRGA_conn!="ND"&IRGA_conn!="ambient"&IRGA_conn!="ref_gas"&cycle!=0)
d.chambers_ambient <- subset(d,IRGA_conn!="ND"&IRGA_conn!="ref_gas")

day_start_decimal <- day_start$hour + (day_start$min/60)
day_end_decimal <- day_end$hour + (day_end$min/60)

cycle_numbers <- unique(d$cycle)[!is.na(unique(d$cycle))]
# never process cycle 0, or the last cycle. They are meaningless and can cause errors
# also don't process cycles 1-3 because they are part of the equilibration time
cycle_numbers <- cycle_numbers[cycle_numbers>2]
cycle_numbers <- cycle_numbers[cycle_numbers!=max(cycle_numbers)]

cycle_table <- data.frame()

# cycle_numbers[9] in expt 38 causes NAN error -> cycle 11

if(length(cycle_numbers)>0) {
  for(this_cycle in cycle_numbers) {
    if(verbose) { print(this_cycle) }
    this_cycle_IRGA_conn.unique <- unique(d[d$cycle==this_cycle,"IRGA_conn"])
    this_cycle_core_names <- this_cycle_IRGA_conn.unique[this_cycle_IRGA_conn.unique%in%core_names]
    next_cycle <- this_cycle+1
    for(core in this_cycle_core_names) {
      #core <- this_cycle_core_names[5]
      if(verbose) { print(core) }
      d.this_core <- d[d$IRGA_conn==core,]
      this_core_temperature_name <- paste(core,"_t",sep="")
      
      # get the relevant data:
      # IRGA_conn indicates when each chamber is vented.
      # The relevant data for each chamber is:
      # 1. the period when the chamber is closed
      # 2. the period when the chamber is open
      # the logger identifies when the chamber is open by naming it in the IRGA_conn field
      
      chamber_closed_time <- max(d.this_core$time_elapsed[d.this_core$cycle==this_cycle])
      chamber_open_time <- min(d.this_core$time_elapsed[d.this_core$cycle==next_cycle])
      
      cycle_length_s <- chamber_open_time - chamber_closed_time
      
      wanted_fields <- c("time","time_elapsed","conn_timer","IRGA_conn","IRGA_CO2","IRGA_Pressure","ref_pressure","ref_temp","IRGA_Humidity",this_core_temperature_name,"time_UTC_decimal_hour","water_content_estimate") # 

      
      d.this_cycle.closed <- d[(d$time_elapsed>chamber_closed_time)&(d$time_elapsed<chamber_open_time),wanted_fields]
      d.this_cycle.open <- subset(d,(cycle==this_cycle)&(IRGA_conn==core))[,colnames(d.this_cycle.closed)]
      cycle_start_hour <- d$time_elapsed_hour[d$time_elapsed==chamber_closed_time]
      
      # define factor for wet or dry (only if we have a mass table)
      if(time_of_water_addition) {
        if (max(d.this_cycle.closed$time)>time_of_water_addition) {
          factor1 <- "wet"
        } else {
          factor1 <- "dry"
        }
      } else {
        factor1 <- "none"
      }
      
      # cycle time of day
      cycle_start_time <- d$time_UTC_decimal_hour[d$time_elapsed==chamber_closed_time]
      cycle_end_time <- d$time_UTC_decimal_hour[d$time_elapsed==chamber_closed_time]
      cycle_start_time_UTC <- d$time[d$time_elapsed==chamber_closed_time]
      cycle_start_time_elapsed <- d$time_elapsed_hour[d$time_elapsed==chamber_closed_time]
      
      if ((cycle_start_time>day_end_decimal)|(cycle_end_time<day_start_decimal)) {
        period <- "night"
      } else if ((cycle_start_time>day_end_decimal)&(cycle_end_time>day_end_decimal)) {
        # this condition is just to catch the occasional case that happens near midnight
        period <- "night"
      }else if ((cycle_start_time<day_end_decimal)&(cycle_end_time>day_start_decimal)){
        period <- "day"
      } else {
        period <- "dusk"   
      }
      
      # identify if any results are from the second day (i.e. the morning when the experiment is stopped)
      cycle_start_time_within_expt <- d$time_elapsed_hour[d$time_elapsed==chamber_closed_time]
      if ((period=="day")&(cycle_start_time_within_expt>14)) {
        period <- "day2"
      }
      
      # NOTE: sub-period calc below was added after the period calc (above). It is probably better.
      # consider removing original period calc (above), and use code below instead to identify period
      
      # identify the cycle time more accurately (define sub_period)
      # time split into 3 hour blocks for comparisons: day1, night1-4, day2
      # This calculation relies on the assumption that lights come on 12 hours after lights off!
      # time_to_lights_off <- as.numeric(day_end - cycle_start_time_UTC)
      time_to_lights_off <- day_end_time_elapsed - cycle_start_time_elapsed
      
      if (time_to_lights_off > 3) {
        sub_period <- "pre_day1"
      }else if ((time_to_lights_off <= 3)&(time_to_lights_off > 0)) {
        sub_period <- "day1"
      }else if ((time_to_lights_off <= 0)&(time_to_lights_off > -3)) {
        sub_period <- "night1"
      }else if ((time_to_lights_off <= -3)&(time_to_lights_off > -6)) {
        sub_period <- "night2"
      }else if ((time_to_lights_off <= -6)&(time_to_lights_off > -9)) {
        sub_period <- "night3"
      }else if ((time_to_lights_off <= -9)&(time_to_lights_off > -12)) {
        sub_period <- "night4"
      }else if ((time_to_lights_off <= -12)&(time_to_lights_off > -15)) {
        sub_period <- "day2"
      }else if (time_to_lights_off <= -15) {
        sub_period <- "post_day2"
      }else {
        sub_period = NA
      }
      
      # temperature during cycle
      cycle_temperature <- d.this_cycle.closed[,c(this_core_temperature_name)]
      
      # ambient pressure during cycle (ref gas system - should exactly match chamber pressure)
      cycle_ambient_pressure <- subset(d.this_cycle.closed)[,c("ref_pressure")]
      
      # ambient pressure during cycle (IRGA)
      cycle_ambient_pressure_IRGA <- subset(d.this_cycle.closed,IRGA_conn=="ambient")[,c("IRGA_Pressure")]
      colnames(d.this_cycle.closed)
      
      # pressure during measurement (IRGA)
      cycle_measurement_pressure_IRGA <- d.this_cycle.open[,c("IRGA_Pressure")]
      
      # humidity during measurement
      cycle_measurement_humidity <- d.this_cycle.open[,c("IRGA_Humidity")]
      
      # CO2 during measurement
      cycle_measurement_CO2 <- d.this_cycle.open[,c("IRGA_CO2")]
      measurements <- d.this_cycle.open[,c("conn_timer","IRGA_CO2")]
      
      # time since water addition
      # first work out the experiment elapsed time that corresponds to when water was added
      water_added_time_elapsed <- min(d$time_elapsed_hour[d$time>time_of_water_addition])
      time_since_water_addition_h <- cycle_start_time_elapsed - water_added_time_elapsed
      
      # time since water addition original calculation - doesn't work as expedcted (usually works but not always)
      #time_since_water_addition_h <- as.numeric(cycle_start_time_UTC - time_of_water_addition)
      
      # amount of water in the core
      water_content_ml <- mean(d.this_cycle.closed$water_content_estimate,na.rm=TRUE)
      
      mean_temperature <- mean(cycle_temperature)
      if(is.na(mean_temperature)) {
        cycle_temperature <- d.this_cycle.closed$ref_temp
        mean_temperature <- mean(cycle_temperature)
      }
      mean_pressure <- mean(cycle_ambient_pressure,na.rm=TRUE)
      CO2_change_sign <- sign(mean(cycle_measurement_CO2[5:15]) - mean(cycle_measurement_CO2[50:60]))
      CO2_change <- CO2_change_sign * (max(cycle_measurement_CO2) - min(cycle_measurement_CO2))
      measurement_pressure_change <- mean(cycle_measurement_pressure_IRGA,na.rm=TRUE) - mean(cycle_ambient_pressure_IRGA,na.rm=TRUE)
      
      # find out the true max CO2 during the cycle:
      testing <- peak_CO2(measurements)
      if (verbose) {
        cat("MAX CO2 detected: ",max(cycle_measurement_CO2),". TRUE MAX CO2 calculated: ",testing)
      }
      
      # flux(cell_temp, CO2_change_ppm, time_period_s, expt_pressure_kPa,verbose)
      # Flux calculation
      # Note that this step can generate NaN errors, but it doesn't matter.
      # The flux function calculates both diff_corrected_low and diff_corrected_high, then selects the best one.
      # Sometimes the other one is NaN
      CO2_flux <- flux(mean_temperature, CO2_change, cycle_length_s, mean_pressure,verbose=FALSE)
      
      cycle <- this_cycle # for naming consistency in cycle table
      cycle_table <- rbind(cycle_table,data.frame(expt,core,cycle,cycle_start_hour,cycle_start_time,cycle_end_time,period,sub_period,time_to_lights_off,factor1,mean_temperature,mean_pressure,CO2_change,measurement_pressure_change,ref_gas_conc,time_since_water_addition_h,water_content_ml,CO2_flux))
      
      #
    }
  }
  
  # Sorting the factors determines the order in ggplots
  #cycle_table$factor1 <- factor(cycle_table$factor1, levels = sort(unique(cycle_table$factor1),decreasing=TRUE))
  
  if (save_cycle_table) {
    
    if (!file.exists(cycle_table_output_dir)){
      dir.create(cycle_table_output_dir)
    } 
    
    write.table(cycle_table,cycle_table_output_file)
    cat("The generated cycle table was saved to ",cycle_table_output_file)
  } else {
    cat("Cycle table was not saved.") 
  }

}
