#Local runs for PINE countries without cluster use

#############################################################################
#############################################################################

# user defined inputs and environments

## user defined paths
### input data
pathi <- '/home/a-m/akraay/vimc/Chaney/demog/' 
### modeling code
pathm <- '/home/a-m/akraay/vimc/Chaney/code/'  
### output files
patho <-'/home/a-m/akraay/vimc/Chaney/results/'
### compilations
#pathr <- '/Users/akraay/Dropbox/Chaney Thesis/RotaDCM-legacy/compile/' 

## define input data files
input_data <- "VIMC_Parameters_contdictcbr_2019-11-25.csv"
### life expectancy file (as .csv) NOTE needs to be LONG format
le.file <- "201910gavi-1_dds-201910_2_life_ex_both.csv"
### coverage data
coverage_data <- "coverage_201910gavi-4_rota-routine-default.csv"

## user defined arguments ##
## vax: {'yes', 'no'}
## vax.scenario: {'novax', 'bestcase', 'routine'}
## doses: {'single', 'two'}
## effect_type: {'overalleffects', 'directeffects'}
## run_type: {'central', 'stochastic'}
## set_run_ids
vax <- "yes"
vax.scenario <- "routine"
doses <-'two'
effect_type <-'overalleffects'
run_type <-"central"
set_run_ids <- 2
mod_time_step <- 1
## this output format includes severe cases and data for all under 5
checkout.file <- paste0('Checkformat_ModelValidation_', vax.scenario, 
	effect_type, 'Doses', doses, "_", Sys.Date(), "_",mod_time_step,'.csv')
## this output file follows VIMC specified output format
vimcout.file <- paste0('VIMCformat_ModelValidation_', vax.scenario, effect_type, 
	'Doses', doses, "_", Sys.Date(), "_",mod_time_step,'.csv')
## set time resolution/step in days
#############################################################################
#############################################################################

# data cleaning

## read in datasets
country_params<-read.csv(paste0(pathi,input_data))
coverage_params<-read.csv(paste0(pathi,coverage_data))

#coverage_params$ID<-seq(from=1, to=nrow(coverage_params), by=1)
#NE<-subset(coverage_params, coverage_params$country_code %in% c('FJI', 'FSM', 'KIR', 'MHL', 'PHL', 'TON', 'TUV', 'VUT', 'WSM', 
#                                                            'CIV', 'GIN', 'GNB', 'NER', 'SEN', 'TCD', 'GTM', 'VEN', 'SYR'))
#rows with negative indirect effects: 17, 31, 32, 35, 37, 38, 51, 59, 68, 74, 81, 90, 91, 97, 99, 104, 106, 107
country_params[,paste0('ihme_under1.', 2000:2010)]
## merge datasets
#country indices with negative IEs: 31 already done and 17
#coverage_params<-coverage_params[c(1, 3:16, 18:30, 33, 34, 36, 39:50, 52:58, 60:67, 69:73, 75:80, 
#                                   82:89, 92:96, 98, 100:103, 105, 108:112),]
coverage_params<-coverage_params[c(23:30, 33, 34, 36, 39:50, 52:58, 60:67, 69:73, 75:80, 
                                   82:89, 92:96, 98, 100:103, 105, 108:112),]
country_params2 <- merge(country_params,coverage_params, 
	by = c("country","country_code"))
param_names <- names(country_params2)
country_params.list <- setNames(split(country_params2, 
	seq(nrow(country_params2))), rownames(country_params2))

#############################################################################
#############################################################################

# run model and compile
tstart <- proc.time()
## run model
source(paste0(pathm, 'model3.R'))
## run code to compile results
