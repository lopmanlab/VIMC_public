# Fit regression model that predicts median (and mean) age of severe infection for all countries
# Last updated: 3/28/2021

# Assign file path
#data<-"C:/Users/ewhall/Box Sync/VIMC/Age of Infection Regression/data/"
data<-"C:/Users/swimf/Documents/GitHub/VIMC/Age of Infection Regression/data/"
outpath<-"C:/Users/swimf/Documents/GitHub/VIMC/VIMC Model Runs/201910gavi_touchstone/Demographic and vaccine coverage data from Montagu/"

##### Predict Median Age of Severe Infection #####
      # Select year of predictor data to use for predicting median age of severe infection
      year<-2015
      
      
      ## Read in outcome data ##
        library(readxl)
        library(openxlsx)
        library(dplyr)
        library(foreach)
        library(reshape2)
        library(tidyr)
        
        files <- list.files(path = paste0(data,"infection data/"), pattern = "*.xlsx", full.names = T)
        
        
        # Read only the columns of interest from each file
        allxlsx.files <- list()  # create a list to populate with xlsx data
        count <- 1
        for (file in files) {
          sheet.names<-getSheetNames(file)
          
          # Some files have a sheet titled "Pre" and some are "Sheet1"
          if("Pre" %in% sheet.names){
            dat <- readxl::read_xlsx(file, sheet = "Pre")[c("RVCount","age_wk", "ISO3_code")]
          }else if ("Sheet1" %in% sheet.names){
            dat <- readxl::read_xlsx(file, sheet = "Sheet1")[c("RVCount","age_wk", "ISO3_code")]
          }else{
            dat <- readxl::read_xlsx(file, sheet = 1)[c("RVCount","age_wk", "ISO3_code")]
            dat$RVCount<-NA
            dat$age_wk<-NA
          }
          allxlsx.files[[count]] <-dat # creat a list of rows from xls files
          count <- count + 1
        }
        
        
        # Loop to calcualte median age of infection for each country
        outcome<-list()
        for (i in 1:length(allxlsx.files)){
          temp<-list()
          
          df<-as.data.frame(allxlsx.files[i])
          
          temp[1]<-df$ISO3_code[1]
          
          # Calculate median age
          df<-df[order(df$age_wk),]
          df$cumsump<-cumsum(df$RVCount)/sum(df$RVCount)
          
         temp[2]<-df[df$cumsump==min(df[df$cumsump>0.5,]$cumsump),"age_wk"]
         temp[3]<-weighted.mean(df$age_wk,df$RVCount)
        
         
         outcome<<-rbind(outcome, temp[1:3])
         rm(df)
         rm(temp)
        }
        outcome.df<-as.data.frame(outcome)
        outcome.df<-setNames(outcome.df, c("country_code", "median_age_wk", "mean_age_wk"))
        outcome.df <- as.data.frame(lapply(outcome.df, unlist))
        outcome.df<- outcome.df[!is.na(outcome.df$median_age_wk),]
      
      ## Read in predictor data ##
      
        # Chandresh's dataset with sub-region values imputed for Eritrea, Iran, Libya, 
        # Papua New Guinea, Syria Venezuela, Mauritania, 
        country_data<-read.csv(paste0(data,"wb_wdi_3_23_19_imputed.csv"))
      
        country_data_year<-country_data[-(1:329),c("Country.Name", "Country.Code", "Indicator.Code",paste0("X",year))]
        country_data_year_wide<-reshape(country_data_year, 
                timevar = "Indicator.Code",
                idvar = c("Country.Name", "Country.Code"),
                direction = "wide")
        colnames(country_data_year_wide)<-c("Country Name", "country_code",paste0("gdp_",year),paste0("gdp_pcap_", year),"gini","u5_mort_1000","pop_density","percent_pop_urban","percent_pop_rural")
        country_data_year_wide<-subset(country_data_year_wide, select=-c(gini,percent_pop_urban))
      #  # WHO Region and Subregion codes
      #  region_data<-read.csv(paste0(data,"who_regions.csv"))
      #  # change regions that are listed as Oceania to Southern Asia
      #  region_data$sub_region[region_data$region=='Oceania'] <- 'Southern Asia'
      #    region_data$region[region_data$region=='Oceania'] <- 'Asia'
      
        # Crude Birth Rate (from Montagu) NOTE CODE REQUIRES WIDE FORMAT
        crude_birth_rate<-read.csv(paste0(data,"201910gavi-1_dds-201910_cbr_both.csv"))
          crude_birth_rate_year<-crude_birth_rate[,c("country_code",paste0("X",year) )]
          colnames(crude_birth_rate_year)<-c("country_code","crude_birth_rate")
        
        # Life expectance (from Montagu) NOTE CODE REQUIRES WIDE FORMAT
        life_expectancy<-read.csv(paste0(data,"201910gavi-1_dds-201910_lx0_both.csv"))
        life_expectancy_year<-life_expectancy[,c("country_code",paste0("X",year) )]
          colnames(life_expectancy_year)<-c("country_code","life_expectancy")
        
        # List of 112 countries we are responsible for modeling
        template <- read.csv(paste0(data,"central-burden-template.201910gavi-1.Rota_Emory-Lopman_standard.csv"))
        vimc_countries <- unique(template$country)
          montagu_data_pre<-merge(crude_birth_rate_year,life_expectancy_year)
          montagu_data <- montagu_data_pre[which(montagu_data_pre$country_code %in% vimc_countries),]
        
        # UNWPP full dataset
        unwpp_data<-read.csv(paste0(data,"unwpp_data.csv"))
        
        # Region and codes
        country_codes<-read.csv(paste0(data,"country-codes_csv.csv"))
        
        
      ## Merge dataset ##
        # Full list of countries
        full.df<-merge(outcome.df, 
                       merge(unwpp_data,
                             merge(montagu_data,
                                   merge(country_codes, country_data_year_wide, by="country_code", all.y=TRUE ), by="country_code", all.y=TRUE), by="country_code_numeric", all.y=TRUE), by="country_code", all.y=TRUE)
        
            # Impute missing data that is in VIMC Montagu dataset
            full.df$crude_birth_rate_2015<-ifelse(is.na(full.df$crude_birth_rate_2015), full.df$crude_birth_rate,full.df$crude_birth_rate_2015)
            full.df$life_expectancy_2015<-ifelse(is.na(full.df$life_expectancy_2015), full.df$life_expectancy,full.df$life_expectancy_2015)
            
            # Change micronesia, melanesia, polynesia, australia and New Zealand to Oceana
            full.df$sub_region<-ifelse(full.df$region=="Oceania", as.character(full.df$region), as.character(full.df$sub_region))
            # Change Southern Europe and Eastern Europe to South/Eastern Europe
            full.df$sub_region<-ifelse(full.df$sub_region=="Eastern Europe"|full.df$sub_region=="Southern Europe", "South/East Europe", full.df$sub_region)
            
            # Categorize crude birth rate into a dichotomous variable
            full.df$crude_birth_rate_dict<-ifelse(full.df$crude_birth_rate_2015<0.02,0,1)

        # Only VIMC Countries
        vimc.countries.df<-full.df[full.df$country_code %in% montagu_data$country_code,]
          # Fix Kosovo country code
          vimc.countries.df$country_code_numeric<-ifelse(vimc.countries.df$country_code=="XK",999,vimc.countries.df$country_code_numeric)
        # Only countries with outcome data
        outcome.countries.df<-full.df[full.df$country_code %in% outcome.df$country_code,]
          

      ## Meadian age - Run model with full oucome dataset.  Make predictions on vimc dataset.
        BestModMedian<-lm(log(median_age_wk)~u5_mort_1000+crude_birth_rate+crude_birth_rate_dict +life_expectancy_2015+percent_pop_rural+gdp_pcap_2015, data=full.df)
        VIMCPredMedAgeWk <- exp(predict(BestModMedian, vimc.countries.df))
        
        ## Mean age - Run model with full oucome dataset.  Make predictions on vimc dataset.
        BestModMean<-lm(log(mean_age_wk)~u5_mort_1000+crude_birth_rate+crude_birth_rate_dict +life_expectancy_2015+percent_pop_rural+gdp_pcap_2015, data=full.df)
        VIMCPredMeanAgeWk <- exp(predict(BestModMean, vimc.countries.df))
        
        ## Predicted Mean age with confidence interval
        VIMCPredMeanAgeWk_interval <- as.data.frame(exp(predict(BestModMean, vimc.countries.df,interval = "confidence")))
        
        ## Sample 200 values of mean age based on uniform distribution for stochastic runs
        stochastic_meanage <- foreach(row=1:nrow(VIMCPredMeanAgeWk_interval),.combine="rbind")%do%{
          runif(200, min=VIMCPredMeanAgeWk_interval$lwr,max=VIMCPredMeanAgeWk_interval$upr)
        }
        colnames(stochastic_meanage) <- paste0("meanage_runid",1:200)
        
        premelt_meanage <- cbind.data.frame(country_code=vimc_countries,stochastic_meanage)
        ## Add predicted to total VIMC dataset
        #select columns to drop
        drops<-c("wuwpp_name", "crude_birth_rate_2015", "life_expectancy_2015",
                 "official_name_en","Country Name", "region", "sub_region")
        vimc.countries.df<-cbind(vimc.countries.df[, !(names(vimc.countries.df) %in% drops)], VIMCPredMedAgeWk, VIMCPredMeanAgeWk,stochastic_meanage)
        
                          
##### Add population sizes from Montagu #####
          # read in data and limit to single year totals (<5) from 2000-2045
          montagu.pops<-read.csv(paste0(data,"201910gavi-1_dds-201910_2_int_pop_both.csv"))
          montagu.pops.2000.2045<-montagu.pops[which(montagu.pops$year>1999 & montagu.pops$year<2046),-7]
          montagu.pops.2000.2045$age<-ifelse(montagu.pops.2000.2045$age_to==0,"u1",
                                             ifelse(montagu.pops.2000.2045$age_to==1,"y1",
                                                    ifelse(montagu.pops.2000.2045$age_to==2,"y2",
                                                           ifelse(montagu.pops.2000.2045$age_to==3,"y3",
                                                                  ifelse(montagu.pops.2000.2045$age_to==4,"y4","over5")))))
          
          montagu.agg<-aggregate(value ~ country_code_numeric+country_code+country+year+age, data=montagu.pops.2000.2045, FUN=sum)
          
          colnames(montagu.agg)[colnames(montagu.agg)=="value"] <- "p"
          # Convert to wide data
          pops.wide.year<-reshape(montagu.agg[,!names(montagu.agg) %in% c("age_to","age_from")],
                             timevar = "year",
                             idvar = c("country_code_numeric", "country_code", "country", "age"),
                             sep = "",
                             direction = "wide"
                             )
          pops.wide<-reshape(pops.wide.year,
                             timevar = "age",
                             idvar = c("country_code_numeric", "country_code", "country"),
                             sep = "",
                             direction = "wide"
          )

                    
 ##### Add IHME deaths data #####
          # Read in IHME GBD Financing dataset to get ISO and IHME country codes
          ihme_iso<-read.csv(paste0(data,"IHME_ISO_codes.csv"))
          
          # Read in IHME mortality data
         # ihme_mortality<-read.csv(paste0(data,"IHME-GBD_2017_DATA-87937164-1.csv"))
           ihme_mortality_pre<-read.csv(paste0(data,"IHME-GBD_2017_DATA-f214422c-1.csv"))
          
          ihme_mortality<-ihme_mortality_pre[ihme_mortality_pre$metric_name=="Number",c("location_id","location_name", "age_name", "year", "val")]
          
          # Convert to wide 
          w1 <- reshape(ihme_mortality, 
                        timevar = "age_name",
                        idvar = c("location_id","location_name","year"),
                        direction = "wide")
          w1$ihme_under1<-w1$`val.<1 year`
          w1$ihme_1_4<-w1$`val.1 to 4`
          w1$ihme_5over<-w1$`val.All Ages`-w1$ihme_1_4-w1$ihme_under1

          
          ihme_wide <- reshape(w1[order(w1$year),c("location_id","location_name", "year", "ihme_under1", "ihme_1_4","ihme_5over")], 
                               timevar = "year",
                               idvar = c("location_id","location_name"),
                               direction = "wide")
          
          ihme_full<-merge(ihme_wide,ihme_iso, by="location_id", all.x=TRUE)
          ihme_rota_deaths<-ihme_full[,-1]
          
          #Using proxy countries to estimate number of rota deaths in Kosovo and Tuvalu
          ihme_mort_proxy<-ihme_mortality_pre[ihme_mortality_pre$metric_name=="Rate",]
          ihme_mort_proxy<-ihme_mort_proxy[ihme_mort_proxy$location_name=="Albania"|
                                           ihme_mort_proxy$location_name=="Montenegro"|
                                           ihme_mort_proxy$location_name=="Serbia"|
                                           ihme_mort_proxy$location_name=="Macedonia"|
                                           ihme_mort_proxy$location_name=="Tonga",c("location_id","location_name", "age_name", "year", "val")]
          #covert proxy mortality rates df to wide 
          proxy1 <- reshape(ihme_mort_proxy, 
                        timevar = "age_name",
                        idvar = c("location_id","location_name","year"),
                        direction = "wide")
          proxy1$ihme_under1_rate<-proxy1$`val.<1 year`/100000
          proxy1$ihme_1_4_rate<-proxy1$`val.1 to 4`/100000
          proxy1$ihme_5over_rate<-(rowMeans(cbind(proxy1$`val.5-14 years`,proxy1$`val.15-49 years`,proxy1$`val.50-69 years`,proxy1$`val.70+ years`)))/100000
          
          #average across SE European countries for Kosovo
          proxy2 <- proxy1[proxy1$location_name=="Albania"|
                             proxy1$location_name=="Montenegro"|
                             proxy1$location_name=="Serbia"|
                             proxy1$location_name=="Macedonia",]
          
          age_names=c('ihme_under1_rate', 'ihme_1_4_rate', 'ihme_5over_rate')
          
          proxy3 <-  foreach(bla=1:length(age_names),.combine='cbind')%do%{
          df=aggregate(proxy2[,age_names[bla]]~proxy2$year,data=proxy2,FUN=mean)
          df$`proxy2[, age_names[bla]]`
          }
          colnames(proxy3) <- age_names
           
         Kosovo_proxy <- cbind.data.frame(year=seq(2000,2017,1),proxy3)
         Tuvalu_proxy <- proxy1[proxy1$location_name=="Tonga",c('year','ihme_under1_rate', 'ihme_1_4_rate', 'ihme_5over_rate')]
         
         Kosovo_mort <- Kosovo_proxy[,c('ihme_under1_rate', 'ihme_1_4_rate', 'ihme_5over_rate')]*cbind.data.frame(montagu.agg[montagu.agg$country=="Kosovo"& montagu.agg$year<=2017& montagu.agg$age=="u1","p"],
                                                                                                                  rowSums(cbind.data.frame(montagu.agg[montagu.agg$country=="Kosovo"& montagu.agg$year<=2017& montagu.agg$age=="y1","p"],
                                                                                                                                           montagu.agg[montagu.agg$country=="Kosovo"& montagu.agg$year<=2017& montagu.agg$age=="y2","p"],
                                                                                                                                           montagu.agg[montagu.agg$country=="Kosovo"& montagu.agg$year<=2017& montagu.agg$age=="y3","p"],
                                                                                                                                           montagu.agg[montagu.agg$country=="Kosovo"& montagu.agg$year<=2017& montagu.agg$age=="y4","p"])),
                                                                                                                  montagu.agg[montagu.agg$country=="Kosovo"& montagu.agg$year<=2017& montagu.agg$age=="over5","p"])
         names(Kosovo_mort)
         Kosovo_mort$year=2000:2017
         Kosovo_mort$location_name="Kosovo"
         Tuvalu_mort <- Tuvalu_proxy[,c('ihme_under1_rate', 'ihme_1_4_rate', 'ihme_5over_rate')]*cbind.data.frame(montagu.agg[montagu.agg$country=="Tuvalu"& montagu.agg$year<=2017& montagu.agg$age=="u1","p"],
                                                                                                                  rowSums(cbind.data.frame(montagu.agg[montagu.agg$country=="Tuvalu"& montagu.agg$year<=2017& montagu.agg$age=="y1","p"],
                                                                                                                                           montagu.agg[montagu.agg$country=="Tuvalu"& montagu.agg$year<=2017& montagu.agg$age=="y2","p"],
                                                                                                                                           montagu.agg[montagu.agg$country=="Tuvalu"& montagu.agg$year<=2017& montagu.agg$age=="y3","p"],
                                                                                                                                           montagu.agg[montagu.agg$country=="Tuvalu"& montagu.agg$year<=2017& montagu.agg$age=="y4","p"])),
                                                                                                                  montagu.agg[montagu.agg$country=="Tuvalu"& montagu.agg$year<=2017& montagu.agg$age=="over5","p"])
         Tuvalu_mort$year=2000:2017
         Tuvalu_mort$location_name="Tuvalu"
         
         
         proxy_mort <- rbind(Kosovo_mort,Tuvalu_mort)
         names(proxy_mort) <- c('ihme_under1','ihme_1_4','ihme_5over', 'year','location_name')
         ihme_proxy_mort <- reshape(proxy_mort[order(proxy_mort$year),c("location_name", "year", "ihme_under1", "ihme_1_4","ihme_5over")], 
                 timevar = "year",
                 idvar = c("location_name"),
                 direction = "wide")
        for(i in 1:2){
          if(ihme_proxy_mort$location_name[i]=="Kosovo"){
           ihme_proxy_mort$country_code[i]="XK"} else {
             ihme_proxy_mort$country_code[i]="TUV"
           }}
         ihme_rota_deaths_all <- rbind(ihme_rota_deaths,ihme_proxy_mort)
##### Add country sepcific income classification #####
          # Read in world bank income classifications
          wbank_inc<-read.csv(paste0(data,"income_classifications.csv"))
          levels(wbank_inc$income_group)
          vac_imm_interval=matrix(NA,nrow=nrow(wbank_inc),ncol=200)
          for(m in 1:nrow(wbank_inc)){
          if(wbank_inc$income_group[m]=="Low income"){
            wbank_inc$vac_imm[m]=0.63;vac_imm_interval[m,]=runif(200, min=0.58,max=0.67)
          } else if(wbank_inc$income_group[m]=="Lower middle income" |wbank_inc$income_group[m]=="Upper middle income"){
            wbank_inc$vac_imm[m]=0.74;vac_imm_interval[m,]=runif(200, min=0.63,max=0.86)
          } else {
           wbank_inc$vac_imm[m]=0.86;vac_imm_interval[m,]=runif(200, min=0.67,max=0.98)}
          }
          

          colnames(vac_imm_interval) <- paste0("vac_imm_runid",1:200)
          
          wbank_inc2 <- cbind.data.frame(wbank_inc,vac_imm_interval)
          
##### Merge together #####  
    parameters_contcbr.df<-merge(pops.wide,
                         merge(ihme_rota_deaths_all, 
                         merge(vimc.countries.df,wbank_inc2, by="country_code",all.x=TRUE), 
                         by="country_code",all.x=TRUE), 
                         by=c("country_code","country_code_numeric"), all.y = TRUE)
    parameters_contcbr.df2<-parameters_contcbr.df[which(parameters_contcbr.df$country_code%in%vimc_countries),]

    save(parameters_contcbr.df2, file=paste0(outpath,"VIMC_Parameters_contdictcbr_",Sys.Date(),".RData"))
    write.csv(parameters_contcbr.df2,file=paste0(outpath,"VIMC_Parameters_contdictcbr_",Sys.Date(),".csv"),row.names = FALSE)
    
    
    parnames_stoch <- names(select(parameters_contcbr.df2,contains("vac_imm_runid")))
    parnames <- names(select(parameters_contcbr.df2,-contains("vac_imm_runid")))
    parnames_stoch2 <- names(select(parameters_contcbr.df2,contains("meanage_runid")))
    parnames2 <- names(select(parameters_contcbr.df2,-contains("meanage_runid")))
    
    parameters_contcbr.df3 <-gather(parameters_contcbr.df2, key = "run_id", value = "vac_imm_stoch",
           parnames_stoch)
    parameters_contcbr.df4 <-gather(parameters_contcbr.df2, key = "run_id2", value = "meanage_stoch",
                                    parnames_stoch2)
    
    parameters_contcbr_long_pre <- cbind.data.frame(parameters_contcbr.df3,meanage_stoch=parameters_contcbr.df4$meanage_stoch)
    parameters_contcbr_long_pre2 <- parameters_contcbr_long_pre[ , -which(names(parameters_contcbr_long_pre) %in% parnames_stoch2)]
    parameters_contcbr_long <- separate(parameters_contcbr_long_pre2,run_id,into=c("tag","run_id"),sep='_imm_runid',remove = TRUE)
    parameters_contcbr_long <- parameters_contcbr_long[,-which(names(parameters_contcbr_long)=="tag")]
    
    save(parameters_contcbr_long, file=paste0(outpath,"VIMC_Parameters_contdictcbr_long",Sys.Date(),".RData"))
    write.csv(parameters_contcbr_long,file=paste0(outpath,"VIMC_Parameters_contdictcbr_long",Sys.Date(),".csv"),row.names = FALSE)
  
##### Write separate csv file for mean age uncertainty data #####
    #write.csv(stochastic_meanage.df,file=paste0(outpath,"VIMC_stochmeanage_contdictcbr",Sys.Date(),".csv"),row.names = FALSE)
