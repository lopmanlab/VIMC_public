library(tidyverse)
library(readxl)


params <- read.csv("/Users/daniellechaney/Desktop/MPH Documents/Thesis/Data Analysis/Demographic and vaccine coverage data from Montagu/VIMC_Parameters_contdictcbr_2021-10-15.csv")
coverage <- read_excel("/Users/daniellechaney/Desktop/MPH Documents/Thesis/Data Analysis/coverage_estimates_gavi.xlsx")
ie_by_year<- read.csv("/Users/daniellechaney/Desktop/MPH Documents/Thesis/Data Analysis/country_ie_by_year.csv")
who_region<- data.frame(read_xlsx("/Users/daniellechaney/Desktop/MPH Documents/Thesis/Data Analysis/WHO Region.xlsx"))
vax_intro <- data.frame(read_xlsx("/Users/daniellechaney/Desktop/MPH Documents/Thesis/Data Analysis/Year of Vaccine Intro.xlsx"))


coverage <- coverage %>% select(country_code, year,coverage) %>% rename(country=country_code) %>% 
  left_join(vax_intro)
params <- params %>% select(country_code, crude_birth_rate, u5_mort_1000) %>% rename(country=country_code) %>% 
  left_join(who_region)
ie_by_year <- ie_by_year %>% select(-X) %>% left_join(vax_intro)

###### Prepping Variables #######

## Recode regions (SEAR is reference)
region1 <- 99
region1[params$region=="EMR"] <- 1
region1[params$region=="AFR" | params$region=="EUR" | params$region=="SEAR" | params$region=="AMR" | params$region=="WPR"] <- 0
params <-cbind(params, data.frame(region1))

region2 <- 99
region2[params$region=="AFR"] <- 1
region2[params$region=="EMR" | params$region=="EUR" | params$region=="SEAR" | params$region=="AMR" | params$region=="WPR"] <- 0
params <-cbind(params, data.frame(region2))

region3 <- 99
region3[params$region=="EUR"] <- 1
region3[params$region=="AFR" | params$region=="EMR" | params$region=="SEAR" | params$region=="AMR" | params$region=="WPR"] <- 0
params <-cbind(params, data.frame(region3))

region4 <- 99
region4[params$region=="AMR"] <- 1
region4[params$region=="AFR" | params$region=="EUR" | params$region=="SEAR" | params$region=="EMR" | params$region=="WPR"] <- 0
params <-cbind(params, data.frame(region4))

region5 <- 99
region5[params$region=="WPR"] <- 1
region5[params$region=="AFR" | params$region=="EUR" | params$region=="SEAR" | params$region=="AMR" | params$region=="EMR"] <- 0
params <-cbind(params, data.frame(region5))

## Select years for vaccine coverage
year0 <- coverage %>% filter(year == intro_year)
year5 <- coverage %>% filter(year == intro_year+5)
year8 <- coverage %>% filter(year == intro_year+8)

## Joining vax coverage into param data frame
params<- cbind(params, data.frame(year0$coverage, year5$coverage, year8$coverage))

## Select IE for each year
year0_ie <- ie_by_year %>% filter(year == intro_year)
year5_ie <- ie_by_year %>% filter(year == intro_year+5)
year8_ie <- ie_by_year %>% filter(year == intro_year+8)

## Joining ie into param data frame
params<- cbind(params, data.frame(year0_ie$ie, year5_ie$ie, year8_ie$ie))


###### Assumptions check #######
## Independence (birth rate and under 5 mortality are highly correlated, but that was expected)
cor(params$crude_birth_rate, params$u5_mort_1000)
cor(params$crude_birth_rate, params$year0.coverage)
cor(params$crude_birth_rate, params$year5.coverage)
cor(params$crude_birth_rate, params$year8.coverage)
cor(params$u5_mort_1000, params$year0.coverage)
cor(params$u5_mort_1000, params$year5.coverage)
cor(params$u5_mort_1000, params$year8.coverage)

#Normality
hist(params$year0_ie.ie)
hist(params$year5_ie.ie)
hist(params$year8_ie.ie)

#Linearity
plot(params$year0_ie.ie, params$crude_birth_rate)
plot(params$year0_ie.ie, params$u5_mort_1000)
plot(params$year0_ie.ie, params$year0.coverage)

plot(params$year5_ie.ie, params$crude_birth_rate)
plot(params$year5_ie.ie, params$u5_mort_1000)
plot(params$year5_ie.ie, params$year5.coverage)

plot(params$year8_ie.ie, params$crude_birth_rate)
plot(params$year8_ie.ie, params$u5_mort_1000)
plot(params$year8_ie.ie, params$year8.coverage)

###### Linear Regression #######
year0_lm <- lm(year0_ie.ie ~ crude_birth_rate + u5_mort_1000 + year0.coverage + region1 + region2 + region3 + region4 + region5, data=params)
summary(year0_lm)
confint(year0_lm, level=0.95)

par(mfrow = c(2,2))
plot(year0_lm)
par(mfrow=c(1,1))

year5_lm <- lm(year5_ie.ie ~ crude_birth_rate + u5_mort_1000 + year5.coverage + region1 + region2 + region3 + region4 + region5, data=params)
summary(year5_lm)
confint(year5_lm, level=0.95)

par(mfrow = c(2,2))
plot(year5_lm)
par(mfrow=c(1,1))

year8_lm <- lm(year8_ie.ie ~ crude_birth_rate + u5_mort_1000 + year8.coverage + region1 + region2 + region3 + region4 + region5, data=params)
summary(year8_lm)
confint(year8_lm, level=0.95)

par(mfrow = c(2,2))
plot(year8_lm)
par(mfrow=c(1,1))

###### "Table 1" of Params by Region #######

wpr_meds <- params %>% filter(region == "WPR")
median(wpr_meds$crude_birth_rate)
median(wpr_meds$u5_mort_1000)
median(wpr_meds$year0.coverage)
median(wpr_meds$year5.coverage)
median(wpr_meds$year8.coverage)

params %>% count(region == "WPR")

###### Countries with Negative Indirect Effects #######
neg_effect <- subset(ie_by_year, ie < 0) %>% left_join(who_region)
neg_effect2 <- neg_effect %>% left_join(params)

###### Negative Effect Logistic Regression #######
log_coverage <- coverage %>% filter(year == intro_year|
                                    year == intro_year+1|
                                    year == intro_year+2|
                                    year == intro_year+3|
                                    year == intro_year+4|
                                    year == intro_year+5|
                                    year == intro_year+6|
                                    year == intro_year+7|
                                    year == intro_year+8|
                                    year == intro_year+9|
                                    year == intro_year+10|
                                    year == intro_year+11|
                                    year == intro_year+12)
log_df <- ie_by_year %>% left_join(log_coverage) %>% drop_na()

logreg_data <- log_df %>% left_join(params)
logreg_data <- logreg_data %>% select(-c(16:21))

neg <- 99
neg[logreg_data$ie < 0] <- 1
neg[logreg_data$ie >= 0] <- 0
logreg_data <- cbind(logreg_data, data.frame(neg))

logreg_data <- logreg_data %>% select(-c(11:15))

region_AFR <- 99
region_AFR[logreg_data$region=="AFR"]<- 1
region_AFR[logreg_data$region == "AMR" | logreg_data$region == "EMR" | logreg_data$region == "EUR" | logreg_data$region == "SEAR" | logreg_data$region == "WPR"] <- 0 

logreg_data <- cbind(logreg_data, data.frame(region_AFR))

region_WPR <- 99
region_WPR[logreg_data$region=="WPR"]<- 1
region_WPR[logreg_data$region == "AMR" | logreg_data$region == "EMR" | logreg_data$region == "EUR" | logreg_data$region == "SEAR" | logreg_data$region == "AFR"] <- 0 

logreg_data <- cbind(logreg_data, data.frame(region_WPR))

region_AMR <- 99
region_AMR[logreg_data$region=="AMR"]<- 1
region_AMR[logreg_data$region == "WPR" | logreg_data$region == "EMR" | logreg_data$region == "EUR" | logreg_data$region == "SEAR" | logreg_data$region == "AFR"] <- 0 

logreg_data <- cbind(logreg_data, data.frame(region_AMR))

region_EMR <- 99
region_EMR[logreg_data$region=="EMR"]<- 1
region_EMR[logreg_data$region == "WPR" | logreg_data$region == "AMR" | logreg_data$region == "EUR" | logreg_data$region == "SEAR" | logreg_data$region == "AFR"] <- 0 

logreg_data <- cbind(logreg_data, data.frame(region_EMR))



#logreg <- glm(formula = neg ~ crude_birth_rate+  u5_mort_1000 + region_AFR + region_WPR + coverage, family= "binomial", data=logreg_data)
#summary(logreg)
#cbind(exp(coef(logreg)), exp(confint.default(logreg)))

#full <- glm(neg ~ region_AFR + region_WPR + crude_birth_rate + u5_mort_1000 + coverage, data = logreg_data, family = "binomial")
#reduced <- glm(neg ~ crude_birth_rate + u5_mort_1000 + coverage, data = logreg_data, family = "binomial")
#anova(full,reduced, test="Chisq")

#drop1(full, test="Chisq")


chi_data <- logreg_data %>% select(country, neg, region, region_AFR, region_WPR, region_EMR, region_AMR) #%>% distinct(country, .keep_all =TRUE)

ever_neg <- 0
ever_neg[logreg_data$country == "CIV"|logreg_data$country =="FJI"|
           logreg_data$country =="FSM"|logreg_data$country =="GIN"|
           logreg_data$country =="GNB"|logreg_data$country =="GTM"|
           logreg_data$country =="KIR"|logreg_data$country =="MHL"|
           logreg_data$country =="NER"|logreg_data$country =="PHL"|
           logreg_data$country =="SEN"|logreg_data$country =="SYR"|
           logreg_data$country =="TCD"|logreg_data$country =="TON"|
           logreg_data$country =="TUV"|logreg_data$country =="VEN"|
           logreg_data$country =="VUT"|logreg_data$country =="WSM"]<- 1

chi_data <- cbind(chi_data, data.frame(ever_neg))
chi_data[is.na(chi_data)] = 0 
chi_data_final <- chi_data %>% distinct(country, .keep_all =TRUE)


test1 <- fisher.test(chi_data_final$ever_neg, chi_data_final$region_AFR)
test1

test2 <- fisher.test(chi_data_final$ever_neg, chi_data_final$region_WPR)
test2

test3 <- fisher.test(chi_data_final$ever_neg, chi_data_final$region_AMR)
test3

test4 <- fisher.test(chi_data_final$ever_neg, chi_data_final$region_EMR)
test4



