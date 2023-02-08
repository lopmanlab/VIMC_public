# Code to multiply predicted incidence by Montagu Pop estimates
#############################################################################

library(reshape2)
library(stringr)
library(plyr)
library(dplyr)
library(foreach)
library(data.table)

#Montagu estimated population sizes
#pops_data=read.csv(paste0(pathi,pvm.file),header=T)


#parammat=read.csv(paste0(pathi,pvm.file),header=T)
parammat <- foreach(listi = 1:length(runs), .combine='rbind') %do% { 
  as.data.frame(runs[[listi]], stringsAsFactors=FALSE)
}

# Create dataset with population data only in same format as PVM output data
year_ends <- c("00", paste0("0",1:9), "10", paste0("1", 1:9), "20", 
  paste0("2", 1:9), "30", paste0("3",1:5))

pops <- foreach(iter = 1:length(year_ends), .combine = "rbind") %do% {
  dat <- cbind(parammat$country, year = paste0("20", year_ends[iter]),
    unname(select(parammat, starts_with(paste0("p20", year_ends[iter])), 
      -ends_with("over5"))))

# note this is summing the pop sizes for 0-5 to get <5 pop size
  cbind(dat, rowSums(dat[,3:7]))
}
colnames(pops) <- c("country", "year", "pop_u1", "pop_y1", "pop_y2", "pop_y3", 
  "pop_y4","pop_u5")

# empty dataframe
results_output <- setNames(data.frame(matrix(ncol = 10, nrow = 0)),
  c("run_id",	"year",	"age", "country",	"country_name",	"cohort_size", 
    "deaths", "cases", "dalys"))
temp.cohort_size <- melt(pops, id.vars = c("country", "year"), 
  variable.name = "age", value.name = "cohort_size")
names(temp.cohort_size)[names(temp.cohort_size) == 'country'] <- 'country_name'
temp.cohort_size$age <- ifelse(temp.cohort_size$age == "pop_u1",0,
                             ifelse(temp.cohort_size$age == "pop_y1",1,
                                    ifelse(temp.cohort_size$age == "pop_y2",2,
                                           ifelse(temp.cohort_size$age == "pop_y3",3,
                                                  ifelse(temp.cohort_size$age == "pop_y4",4,
                                                         ifelse(temp.cohort_size$age == "pop_u5",5,temp.cohort_size$age))))))
                               

for (lnum in 1:nrow(parammat)){
  #Country names
  country <- parammat[lnum, "country"]

  #Country names
  c_code <- parammat[lnum, "country_code"]
  
  #iterate over stochastic runs for each country
  #for (i in 1:length(set_run_ids)){
  #Run ID
  #run_idlab=set_run_ids[i]
  run_idlab <- 'central'
  
  infile_case_vac <- paste0(patho,'total_cases_3q2r_2mab_',vax.scenario,c_code,run_idlab,"_",mod_time_step,".csv")
  infile_sevcase_vac <- paste0(patho,'SevCases_3q2r_2mab_',vax.scenario,c_code,run_idlab,"_",mod_time_step,".csv")
  infile_deaths_vac <- paste0(patho,'deaths_3q2r_2mab_',vax.scenario,c_code,run_idlab,"_",mod_time_step,".csv")
  infile_DALYs_vac <- paste0(patho,'DALYs_3q2r_2mab_',vax.scenario,c_code,run_idlab,"_",mod_time_step,".csv")
  
  cases_vac <- try(as.data.frame(fread(infile_case_vac,header = T))) 
  if(class(cases_vac) == "try-error") {cat('!!file does not exist!!')}
  
  sevcases_vac <- try(as.data.frame(fread(infile_sevcase_vac,header = T)))
  if(class(sevcases_vac) == "try-error") {cat('!!file does not exist!!')}

  deaths_vac <- try(as.data.frame(fread(infile_deaths_vac,header = T)))
  if(class(deaths_vac) == "try-error") {cat('!!file does not exist!!')}

  DALYs_vac <- try(as.data.frame(fread(infile_DALYs_vac,header = T)))
  if(class(DALYs_vac) == "try-error") {cat('!!file does not exist!!')}
  
  #with vaccination
  if(class(cases_vac)=="data.frame") {cases_vacnum=cbind(cases_vac[,c('cases_0_11','cases_12_23','cases_24_35','cases_36_47','cases_48_59','cases_0_59')]*pops[pops$country==country,c('pop_u1','pop_y1','pop_y2','pop_y3','pop_y4','pop_u5')],run_id=rep(run_idlab,nrow(cases_vac)), year=seq(2000,2035))
  outfile_case_vac=paste0(pathr,'num_total_cases_3q2r_2mab_',vax.scenario,c_code,run_idlab,"_",mod_time_step,'_', as.numeric(Sys.time()), ".csv")
  write.csv(cases_vacnum,outfile_case_vac)}
  
  if(class(sevcases_vac)=="data.frame") {sevcases_vacnum=cbind(sevcases_vac[,c('sevcases_0_11','sevcases_12_23','sevcases_24_35','sevcases_36_47','sevcases_48_59','sevcases_0_59')]*pops[pops$country==country,c('pop_u1','pop_y1','pop_y2','pop_y3','pop_y4','pop_u5')],run_id=rep(run_idlab,nrow(sevcases_vac)), year=seq(2000,2035))
  outfile_sevcase_vac=paste0(pathr,'num_SevCases_3q2r_2mab_',vax.scenario,c_code,run_idlab,"_",mod_time_step,'_', as.numeric(Sys.time()), ".csv")
  write.csv(sevcases_vacnum,outfile_sevcase_vac)}
  
  if(class(deaths_vac)=="data.frame") {deaths_vacnum=cbind(deaths_vac[,c('deaths_0_11','deaths_12_23','deaths_24_35','deaths_36_47','deaths_48_59','deaths_0_59')]*pops[pops$country==country,c('pop_u1','pop_y1','pop_y2','pop_y3','pop_y4','pop_u5')],run_id=rep(run_idlab,nrow(deaths_vac)), year=seq(2000,2035))
  outfile_deaths_vac=paste0(pathr,'num_deaths_3q2r_2mab_',vax.scenario,c_code,run_idlab,"_",mod_time_step,'_', as.numeric(Sys.time()), ".csv")
  write.csv(deaths_vacnum,outfile_deaths_vac)}
  
  if(class(DALYs_vac)=="data.frame") {DALYs_vacnum=cbind(DALYs_vac[,c('DALYs_0_11','DALYs_12_23','DALYs_24_35','DALYs_36_47','DALYs_48_59','DALYs_0_59')]*pops[pops$country==country,c('pop_u1','pop_y1','pop_y2','pop_y3','pop_y4','pop_u5')],run_id=rep(run_idlab,nrow(DALYs_vac)), year=seq(2000,2035))
  outfile_DALYs_vac=paste0(pathr,'num_DALYs_3q2r_2mab_',vax.scenario,c_code,run_idlab,"_",mod_time_step,'_', as.numeric(Sys.time()), ".csv")
  write.csv(DALYs_vacnum,outfile_DALYs_vac)}
  
  temp.long.cases<-melt(cases_vacnum, id.vars=c("run_id", "year"),variable.name = "age",value.name = "cases")
  temp.long.cases$age<-str_sub(temp.long.cases$age, start= -5)
  
  temp.long.sevcases<-melt(sevcases_vacnum, id.vars=c("run_id", "year"),variable.name = "age",value.name = "severe_cases")
  temp.long.sevcases$age<-str_sub(temp.long.sevcases$age, start= -5)
  
  temp.long.deaths<-melt(deaths_vacnum, id.vars=c("run_id", "year"),variable.name = "age",value.name = "deaths")
  temp.long.deaths$age<-str_sub(temp.long.deaths$age, start= -5)
  
  temp.long.dalys<-melt(DALYs_vacnum, id.vars=c("run_id", "year"),variable.name = "age",value.name = "dalys")
  temp.long.dalys$age<-str_sub(temp.long.dalys$age, start= -5)

                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
  # merge temporary data frames and add additional columns
  temp.long<-merge(temp.long.cases,merge(temp.long.deaths,merge(temp.long.sevcases,temp.long.dalys,by=c("run_id","year","age")),by=c("run_id","year","age")),by=c("run_id","year","age"))
  temp.long$disease<-"Rota"
  temp.long$country<-c_code
  temp.long$country_name<-country
  
  #append to formatted dataframe
  results_output<-rbind.data.frame(results_output, temp.long)
}
#  }

# format age values in dataframe
results_output$age<-ifelse(results_output$age=="_0_11",0,
                           ifelse(results_output$age=="12_23",1,
                                  ifelse(results_output$age=="24_35",2,
                                         ifelse(results_output$age=="36_47",3,
                                                ifelse(results_output$age=="48_59",4,
                                                       ifelse(results_output$age=="_0_59",5,results_output$age))))))

# merge in cohort size
complete_results_output<-merge(results_output, temp.cohort_size, by=c("country_name","year","age"), all.x = TRUE)
# Reorder columns and sort
complete_results_output<-complete_results_output[order(complete_results_output$run_id, complete_results_output$age,complete_results_output$country,complete_results_output$year )
                                                 ,c("run_id",	"year",	"age",	"country",	"country_name",	"cohort_size",	"deaths",	"cases",	"dalys", "severe_cases")]


# save final output
write.csv(complete_results_output,paste0(pathr,checkout.file), row.names = FALSE)

# Create file that only contains VIMC requested output and save
vimc_results_output <- complete_results_output[complete_results_output$age!=5,]
vimc_results_output <- vimc_results_output[,c("run_id",	"year",	"age",	"country",	"country_name",	"cohort_size",	"deaths",	"cases",	"dalys")]

write.csv(vimc_results_output,paste0(pathr,vimcout.file), row.names = FALSE)
