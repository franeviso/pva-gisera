setwd("~/Documents/Programming/C++/Population_viability_simulations")

seedbank_data  <- read.table('demo_sb_fire.dat', skip = 10, header = TRUE, sep ='\t', stringsAsFactors = FALSE)
View(seedbank_data)
length(names(seedbank_data))
colnames(seedbank_data) <- c("Gen", "Nrun", "Veg","Veg_err", "Repro", "Repro_err", "Age", "Age_err", "Empty", "Empty_err", "Lrep", "Lrep_err", "Mate", "Mate_err")

seedbank_data_num <- data.frame(sapply(seedbank_data, function(x) as.numeric(as.character(x))))
View(seedbank_data_num)

seedbank_data_num$Total_pop <- seedbank_data_num$Veg + seedbank_data_num$Repro
cleaned_total <- seedbank_data_num$Total_pop[which(!is.na(seedbank_data_num$Total_pop))]


log_stochastic_growth <- function(data_sims){
   log_lambda <- array(0,dim=c(length(data_sims),0))
   for(i in 1:length(data_sims)-1){
     log_lambda[i] = log(data_sims[i+1]/data_sims[i])
   }
   return(log_lambda)
}

maxt <- length(cleaned_total)
log_lambda_s <- mean(log_stochastic_growth(cleaned_total))
log_lambda_s_se <- 1.96*sqrt(var(log_lambda)/maxt)
CI_log_lambda = c(log_lambda_s  - log_lambda_s_se, log_lambda_s + log_lambda_s_se)


quasi_ext_threshold <- function(data_sims, threshold){
  quasi <- rep(0,length(data_sims))
  for(t in 1:length(data_sims)){
    if(data_sims[t] < threshold){
      quasi[t] = quasi[t] + 1
      break;
    }
  }
  return(quasi)
}
