library(tidyverse)
library(readxl)
library(data.table)

######################################## DIRECT EFFECTS ##########################################

##Read in datasets, remove unneeded columns##
novax_og <- data.frame(read_xlsx("/Users/daniellechaney/Desktop/MPH Documents/Thesis/Original Data/novax.xlsx"))
novax1 <- select(novax_og, -1, -9)

routine_og <- data.frame(read_xlsx("/Users/daniellechaney/Desktop/MPH Documents/Thesis/Original Data/routine_direct.xlsx"))
routine1 <- select(routine_og, -1, -9)

vax_intro <- data.frame(read_xlsx("/Users/daniellechaney/Desktop/MPH Documents/Thesis/Data Analysis/Year of Vaccine Intro.xlsx"))

##Join vax_intro to other datasets##

routine <- routine1 %>% left_join(vax_intro)
novax <- novax1 %>% left_join(vax_intro) %>% mutate(novax_death = deaths)
novax_clean <- novax %>% select(-5, -6, -7, -8)
routine_clean <- routine %>% left_join(novax_clean)

## Direct effect = num deaths in default / num deaths in no vax after 10 years (By age?) ##
## Gives PERCENT REDUCTION in deaths ##

#Loop to select time points 
countries <- list(routine_clean$country) 
age <- list(routine_clean$age)

result1 <- list()
for (x in countries) {
  
df <- routine_clean %>% filter(country == x) %>% filter(year == intro_year-5 |year == intro_year+0 |year == intro_year+1 
  | year == intro_year+2 | year == intro_year+3 |year == intro_year+4 |year == intro_year+5 | year == intro_year+6 | year == intro_year+7 |
  year == intro_year+8 | year == intro_year+9 | year == intro_year+10 | year == intro_year+11 | year == intro_year+12)
  
  result1[[1 + length(result1)]] <- data.frame(df)

}

result1 <- as.data.frame(do.call(rbind, result1))

#Loop to calculate DEs#
countries2 <- list(result1$country)
age2 <- list(result1$age)

result_de_age <- list()
for (x in countries2) {
  
  for (y in age2) {
  
    de = 1-((result1 %>% filter(country == x, age == y) %>% select(deaths))/(result1 %>% filter(country == x, age == y) %>% select(novax_death)))
    
    result_de_age[[1 + length(result_de_age)]] <- data.frame(result1$year,result1$country, result1$age, result1$deaths,result1$novax_death,de)
  } 
}

result_de_age <- as.data.frame(do.call(rbind, result_de_age))
result_de_age <- result_de_age %>% rename(de = deaths, year = result1.year, country = result1.country, age = result1.age, vax_deaths = result1.deaths, novax_deaths = result1.novax_death)

## Total Direct Effect per Country (not split in age groups) ##

#Loop to select time points (don't want to include 5 year pre mark like above) 
countries2 <- list(routine_clean$country) 
age <- list(routine_clean$age)


result2 <- list()
for (x in countries2) {
  
  df <- routine_clean %>% filter(country == x) %>% filter(year == intro_year+0 |year == intro_year+1 
                                                          | year == intro_year+2 | year == intro_year+3 |year == intro_year+4 |year == intro_year+5 | year == intro_year+6 | year == intro_year+7 |
                                                            year == intro_year+8 | year == intro_year+9 | year == intro_year+10 | year == intro_year+11 | year == intro_year+12)
  
  result2[[1 + length(result2)]] <- data.frame(df)
  
}

result2 <- as.data.frame(do.call(rbind, result2))

#Find country DE (without age groups) (gives one DE per country from time of vax intro + 12 years)

vax_deaths <- aggregate(result2$deaths, by=list(country=result2$country), FUN=sum) %>% mutate(vax_deaths = x) %>% select(-2)
novax_deaths <- aggregate(result2$novax_death, by=list(country=result2$country), FUN=sum) %>% mutate(novax_deaths = x) %>% select(-2)
death_sums <- vax_deaths %>% left_join(novax_deaths)
country_de_noage <- death_sums %>% mutate(country_de = (1-(death_sums$vax_deaths/death_sums$novax_deaths)))

#Find country DE (without age groups)
df4 <- result1 %>% select(-2,-4,-5,-7,-8) %>% rename(vax_deaths= deaths, novax_deaths = novax_death)


by_year <- as.data.table(df4)[, sum(vax_deaths), by = .(country, year)]
by_year1 <- by_year %>% select(country, year, V1) %>% rename(vax_deaths = V1)
by_year2 <- as.data.table(df4)[, sum(novax_deaths), by = .(country, year)]
by_year3 <- by_year2 %>% select(country, year, V1) %>% rename(novax_deaths = V1)

country_de_by_year <- by_year1 %>% left_join(by_year3)
country_de_by_year <- country_de_by_year %>% mutate(de = 1-(vax_deaths/novax_deaths))

## Writing main datatables to csv files
write.csv(result_de_age, file="result_de_age.csv")
write.csv(country_de_noage, file="country_de_noage.csv")
write.csv(country_de_by_year, file="country_de_by_year.csv")

######################################## OVERALL EFFECTS ##########################################

##Read in datasets, remove unneeded columns##

routine_overall_og <- data.frame(read_xlsx("/Users/daniellechaney/Desktop/MPH Documents/Thesis/Original Data/routine_overall_double.xlsx"))
routine_overall1 <- select(routine_overall_og, -1, -9)

##Join vax_intro to other datasets##

routine_overall <- routine_overall1 %>% left_join(vax_intro)
routine_overall_clean <- routine_overall %>% left_join(novax_clean)


## Gives PERCENT REDUCTION in deaths ##

#Loop to select time points 
countries_oe <- list(routine_overall_clean$country) 
age_oe <- list(routine_overall_clean$age)

result1_oe <- list()
for (x in countries_oe) {
  
  df <- routine_overall_clean %>% filter(country == x) %>% filter(year == intro_year-5 |year == intro_year+0 |year == intro_year+1 
                                                          | year == intro_year+2 | year == intro_year+3 |year == intro_year+4 |year == intro_year+5 | year == intro_year+6 | year == intro_year+7 |
                                                            year == intro_year+8 | year == intro_year+9 | year == intro_year+10 | year == intro_year+11 | year == intro_year+12)
  
  result1_oe[[1 + length(result1_oe)]] <- data.frame(df)
  
}

result1_oe <- as.data.frame(do.call(rbind, result1_oe))

#Loop to calculate OEs#
countries2_oe <- list(result1_oe$country)
age2_oe <- list(result1_oe$age)

result_oe_age <- list()
for (x in countries2_oe) {
  
  for (y in age2_oe) {
    
    oe = 1-((result1_oe %>% filter(country == x, age == y) %>% select(deaths))/(result1_oe %>% filter(country == x, age == y) %>% select(novax_death)))
    
    result_oe_age[[1 + length(result_oe_age)]] <- data.frame(result1_oe$year,result1_oe$country, result1_oe$age, result1_oe$deaths,result1_oe$novax_death,oe)
  } 
}

result_oe_age <- as.data.frame(do.call(rbind, result_oe_age))
result_oe_age <- result_oe_age %>% rename(oe = deaths, year = result1_oe.year, country = result1_oe.country, age = result1_oe.age, vax_deaths = result1_oe.deaths, novax_deaths = result1_oe.novax_death)

## Total Overall Effect per Country (not split in age groups) ##

#Loop to select time points (don't want to include 5 year pre mark like above) 
countries2_oe <- list(routine_overall_clean$country) 
age_oe <- list(routine_overall_clean$age)


result2_oe <- list()
for (x in countries2_oe) {
  
  df <- routine_overall_clean %>% filter(country == x) %>% filter(year == intro_year+0 |year == intro_year+1 
                                                          | year == intro_year+2 | year == intro_year+3 |year == intro_year+4 |year == intro_year+5 | year == intro_year+6 | year == intro_year+7 |
                                                            year == intro_year+8 | year == intro_year+9 | year == intro_year+10 | year == intro_year+11 | year == intro_year+12)
  
  result2_oe[[1 + length(result2_oe)]] <- data.frame(df)
  
}

result2_oe <- as.data.frame(do.call(rbind, result2_oe))

#Find country OE (without age groups) (gives one OE per country from time of vax intro + 12 years)

vax_deaths_oe <- aggregate(result2_oe$deaths, by=list(country=result2_oe$country), FUN=sum) %>% mutate(vax_deaths = x) %>% select(-2)
novax_deaths_oe <- aggregate(result2_oe$novax_death, by=list(country=result2_oe$country), FUN=sum) %>% mutate(novax_deaths = x) %>% select(-2)
death_sums_oe <- vax_deaths_oe %>% left_join(novax_deaths_oe)
country_oe_noage <- death_sums_oe %>% mutate(country_oe = 1-(death_sums_oe$vax_deaths/death_sums_oe$novax_deaths))

#Find country OE (without age groups)
df4_oe <- result1_oe %>% select(-2,-4,-5,-7,-8) %>% rename(vax_deaths= deaths, novax_deaths = novax_death)


by_year_oe <- as.data.table(df4_oe)[, sum(vax_deaths), by = .(country, year)]
by_year1_oe <- by_year_oe %>% select(country, year, V1) %>% rename(vax_deaths = V1)
by_year2_oe <- as.data.table(df4_oe)[, sum(novax_deaths), by = .(country, year)]
by_year3_oe <- by_year2_oe %>% select(country, year, V1) %>% rename(novax_deaths = V1)

country_oe_by_year <- by_year1_oe %>% left_join(by_year3_oe)
country_oe_by_year <- country_oe_by_year %>% mutate(oe = 1-(vax_deaths/novax_deaths))

## Writing main datatables to csv files
write.csv(result_oe_age, file="result_oe_age.csv")
write.csv(country_oe_noage, file="country_oe_noage.csv")
write.csv(country_oe_by_year, file="country_oe_by_year.csv")

######################################## INDIRECT EFFECTS ##########################################


a <- result_de_age %>% select(year, country, age, de)
b <- result_oe_age %>% select(year, country, age, oe)
result_ie_age <- a %>% left_join(b) %>% mutate(ie = oe-de)


c <- country_de_noage %>% select(country, country_de)
d <- country_oe_noage %>% select(country, country_oe)
country_ie_noage <- c %>% left_join(d) %>% mutate(ie = country_oe-country_de)

e <- country_de_by_year %>% select(country, year, de)
f <- country_oe_by_year %>% select(country, year, oe)
country_ie_by_year <- e %>% left_join(f) %>% mutate(ie = oe-de)

#Writing main datatables to csv files
write.csv(result_ie_age, file="result_ie_age.csv")
write.csv(country_ie_noage, file="country_ie_noage.csv")
write.csv(country_ie_by_year, file="country_ie_by_year.csv")


