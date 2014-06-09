# These calculations are based on a spreadsheet made by Steve Hoon
# setwd("C:\\Users\\david\\mmu\\lab_book\\kalahari\\files\\microcosms")
# 0.183333333333333 * 3600
# example: flux(28.5,5.42,660,101325,TRUE)

# v.5 returns NA if given CO2_change = NA
# functions dir added in v.4. - calling script needs to define it
# if using standalone then define functions dir by uncommenting and adjusting following line:
# functions_dir <- "C:\\sync\\lab_book\\kalahari\\files\\microcosms\\R"

# v.6 adding in new function to determine the peak CO2 conc
# based on analysis provided by SRH 
# Growth Cabinet Respiration Chamber Modelling and analysis.docx
# MicrocosmModelling.xlsx


const_file <- "data/CO2_flux_calc_constants.txt"
dimensions_file <- "data/microcosm_dimensions.txt"



#Constants
constants <- read.table(const_file,header=T)
microcosm <- read.table(dimensions_file,header=T,stringsAsFactors=FALSE)

# Constants
#constants <- read.table("CO2_flux_calc_constants.txt",header=T)
#microcosm <- read.table("microcosm_dimensions.txt",header=T,stringsAsFactors=FALSE)
#microcosm <- read.table("field_chamber_dimensions.txt",header=T,stringsAsFactors=FALSE)

# derived variables
crust_area <- pi*(microcosm["crust_diameter",1]/2)^2
cell_volume <- (pi*(microcosm["cell_diameter",1]/2)^2)*(microcosm["cell_height",1])
microcosm["crust_area",] <- list(crust_area,"mm^2 (derived)")
microcosm["crust_area_metre",] <- list(crust_area*10^-6,"m^2 (derived)")
microcosm["cell_volume",] <- list(cell_volume,"mm^3 (derived)")
microcosm["cell_volume_litre",] <- list(cell_volume*10^-6,"litre (derived)")

flux_consts <- function(x) {
  print(constants)
  print("")
  print(format(microcosm, scientific = FALSE,digits=2))
}

flux <- function(cell_temp, CO2_change_ppm, time_period_s, expt_pressure_kPa,verbose) {
  if (is.na(CO2_change_ppm)) {return(NA)}
  time_period_h <- time_period_s / 3600
  input <- c(cell_temp, CO2_change_ppm, time_period_h, expt_pressure_kPa)
  input <- data.frame(input,row.names=c("cell temperature","CO2 change ppm","time period (h)","pressure (kPa)"))
  CO2_change_v_pc <- CO2_change_ppm / 10000 # need to check this with Steve -> DONE, it is right now.
  
  # calculation 
  #Observed Two Vial [a_obs] (v/v CO2 hr-1)  
  # CO2_change / interval / 100
  obs_corrected_2vial_vv <- CO2_change_v_pc / time_period_h / 100
  
  #Observed Two Vial C mass flux (mg C m-2 hr-1)  
  # =AtWt_C/(NTP / STP vol 1 g-mole of gas V0) *obs_corrected_2vial_vv *Cell Volume (litres)/Cell (soil) Area = 
  obs_corrected_2vial_mg <- constants["AtWt_C",1] / constants["NTP/STP_vol_1_g-mole_of_gas_V0",1] * obs_corrected_2vial_vv * microcosm["cell_volume_litre",1] / microcosm["crust_area_metre",1]
  
  #STP corrected Observed Two Vial C mass flux (mg C m-2 hr-1)	
  #=obs_corrected_2vial_mg * (expt_pressure_kPa/Standard_Pressure_(P0))   *   (cell_temp/Standard_Temperature_(T0)+1)
  STP_corrected_2vial <- obs_corrected_2vial_mg * (expt_pressure_kPa/constants["Standard_Pressure_(P0)",1])  *   (cell_temp/constants["Standard_Temperature_(T0)",1]+1)
  
  #Diff. Corrected [a_obs] low rates C mass flux (mg C m-2 hr-1)	
  # sign(obs_corrected_2vial_vv) * 
  # (-1 * (1 + mq * time_period_h) +
  # sqrt((1 + mq * time_period_h)^2 + 4 *
  # nq * time_period_h * abs(obs_corrected_2vial_vv))) / 
  # (2* nq * time_period_h) / 
  # NTP--- * (cell_vol/cell_area) * atWtC
  
  diff_corrected_low <- 
    sign(obs_corrected_2vial_vv) * 
    (-1 * ( 1 + constants["mq",1] * time_period_h) +
    sqrt(( 1 + constants["mq",1] * time_period_h)^2 + 4 * 
    constants["nq",1] * time_period_h * abs(obs_corrected_2vial_vv)))   /
    (2 * constants["nq",1] * time_period_h) / 
    constants["NTP/STP_vol_1_g-mole_of_gas_V0",1] * (microcosm["cell_volume_litre",1]/microcosm["crust_area_metre",1]) * constants["AtWt_C",1]
  
  #Diff. Corrected [a_obs] high rates C mass flux (mg C m-2 hr-1)
  # = (H25 -$B$12*G25)/(1+$B$13*G25)/$B$7*(D25/$B$5)*$B$8
  diff_corrected_high <- (obs_corrected_2vial_vv - constants["d",1] * time_period_h) /
    (1 + constants["m",1] * time_period_h) /
    constants["NTP/STP_vol_1_g-mole_of_gas_V0",1] * (microcosm["cell_volume_litre",1] / microcosm["crust_area_metre",1]) * constants["AtWt_C",1]
  
  #Diffusion Corrected CO2 Flux mg C m-2 hr-1
  # = IF(ABS(H25)<Limit,K25,L25)
  if (abs(obs_corrected_2vial_vv)<constants["Limit",1]) {
    diff_corrected <- diff_corrected_low
  } else   
  {
    diff_corrected <- diff_corrected_high
  }
  
  # STP & Diffusion  Corrected CO2 Flux mg C m-2 hr-2
  #=M25*(B$19/B$9)*(C25/B$10+1)
  STP_diff_corrected <- diff_corrected *
    (expt_pressure_kPa/constants["Standard_Pressure_(P0)",1]) *
    (cell_temp / constants["Standard_Temperature_(T0)",1]+1)
  
  results_rows <- c("obs_corrected_2vial_vv","obs_corrected_2vial_mg","STP_corrected_2vial","diff_corrected_low","diff_corrected_high","diff_corrected","STP_diff_corrected")
  results <- c(obs_corrected_2vial_vv,obs_corrected_2vial_mg,STP_corrected_2vial,diff_corrected_low,diff_corrected_high,diff_corrected,STP_diff_corrected)
  
  results <- data.frame(results,row.names=results_rows)
  #results <- data.frame(row.names=c(2:8),c(obs_corrected_2vial_vv,obs_corrected_2vial_mg,STP_corrected_2vial,diff_corrected_low,diff_corrected_high,diff_corrected,STP_diff_corrected))
  #results <- data.frame(rows=c(2:8),c(obs_corrected_2vial_vv,obs_corrected_2vial_mg,STP_corrected_2vial,diff_corrected_low,diff_corrected_high,diff_corrected,STP_diff_corrected))
  
  if(verbose) {
    cat("Input\n")
    print(input,digits=3)
    cat("\nChamber details\n")
    print(microcosm,digits=3)
    cat("\nConstants\n")
    print(constants,digits=3)
    cat("\nResults\n")
    print(results,digits=3)
  }
  return (STP_diff_corrected)
}


peak_CO2 <- function(measurements) {
  plot_min <- min(measurements$IRGA_CO2)-20
  plot_max <- max(measurements$IRGA_CO2)+20
  
  # fitting parameters
  tclsd <- 3.3
  
  # Physical parameters:
  S <- 3
  Vc <- 108
  Ve <- 3.4
  dt <- 0.5
  
  # Derived Model Constants:
  dV <- S * dt
  Alpha <- S/Vc
  Beta <- S/Ve
  B_fac <- 1/(1-Alpha/Beta)
  
  # Analytical Model
  t_peak <- tclsd+(Ve/S)*log(1+Vc/Ve)
  t_peak_millis <- t_peak * 1000
  
  # calculate peak and ref gas:
  # calculations use 2 points on the data curve

  # peak
  # =IF(Q48<(t_Peak+1),"",((T48*R48-T83*R83)/(T48-T83)))
  
  # first we discard the points before the peak: (IF(Q48<(t_Peak+1))
  measurements <- measurements[measurements$conn_timer>(t_peak_millis+1000),]
  # work out the middle measurement index, with a bit of leeway in-case of odd number of measurements
  middle_measurement <- round((nrow(measurements)-2)/2)
  # now work through the first half of measurements, and calculate results in conjunction with second-half measurements:
  #1:middle_measurement
  

  
  # Work out Cr
  for(i in 1:middle_measurement){
    j <- i+middle_measurement
    t1 <- (measurements[i,"conn_timer"])/1000
    t2 <- (measurements[j,"conn_timer"])/1000
    c1 <- measurements[i,"IRGA_CO2"]
    c2 <- measurements[j,"IRGA_CO2"]
    # fjfk=IF(Q49<(tclsd+1),"NA",((EXP(Alpha*(Q49-tclsd)))/(1-EXP(-Beta*(Q49-tclsd)))))
    fj <- (exp(Alpha*(t1-tclsd)))/(1-exp(-Beta*(t1-tclsd)))
    fk <- (exp(Alpha*(t2-tclsd)))/(1-exp(-Beta*(t2-tclsd)))
    #((T52*R52-T87*R87)/(T52-T87)))
    Cr <- ((fj*c1-fk*c2)/(fj-fk))
    measurements[i,"Cr"] <- Cr
    #cat(i,j,t1,t2,c1,c2,"\n")  
  }
  
  Cr.mean <- mean(measurements$Cr,na.rm=TRUE)
  
  # Work out Cci
  for(i in 1:middle_measurement){
    j <- i+middle_measurement
    t1 <- (measurements[i,"conn_timer"])/1000
    t2 <- (measurements[j,"conn_timer"])/1000
    c1 <- measurements[i,"IRGA_CO2"]
    c2 <- measurements[j,"IRGA_CO2"]
    # fjfk=IF(Q49<(tclsd+1),"NA",((EXP(Alpha*(Q49-tclsd)))/(1-EXP(-Beta*(Q49-tclsd)))))
    fj <- (exp(Alpha*(t1-tclsd)))/(1-exp(-Beta*(t1-tclsd)))
    fk <- (exp(Alpha*(t2-tclsd)))/(1-exp(-Beta*(t2-tclsd)))
    # (Cr_Mean+(R48-Cr_Mean)*B_fac*T48))
    Cci <- Cr.mean+(c1-Cr.mean)*B_fac*fj
    measurements[i,"Cci"] <- Cci
  }
  
  if(verbose) {
    plot(measurements$conn_timer,measurements$IRGA_CO2,ylim=c(plot_min,plot_max))
    # plot the detected peak position:
    abline(v=t_peak_millis)
    # now plot the remaining points in green so we can see which ones are used for fitting
    points(measurements$conn_timer,measurements$IRGA_CO2,ylim=c(plot_min,plot_max),col="green",pch=16)
    points(measurements$conn_timer,measurements$Cr,col="blue")
    points(measurements$conn_timer,measurements$Cci,col="red")
  }
  
  # ref gas
  Cci.mean <- mean(measurements$Cci,na.rm=TRUE)
  return (Cci.mean)
}

#measurements <- d.this_cycle.open[,c("conn_timer","IRGA_CO2")]
#peak_CO2(measurements)
