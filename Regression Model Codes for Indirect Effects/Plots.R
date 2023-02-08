library(tidyverse)
library(readxl)
library(patchwork)
library(RColorBrewer)
library(Hmisc)
library(ggridges)

##Read in datasets
de_by_year<- read.csv("/Users/daniellechaney/Desktop/MPH Documents/Thesis/Data Analysis/country_de_by_year.csv")
de_by_year<- select(de_by_year, -1)
oe_by_year<- read.csv("/Users/daniellechaney/Desktop/MPH Documents/Thesis/Data Analysis/country_oe_by_year.csv")
oe_by_year<- select(oe_by_year, -1)
ie_by_year<- read.csv("/Users/daniellechaney/Desktop/MPH Documents/Thesis/Data Analysis/country_ie_by_year.csv")
ie_by_year<- select(ie_by_year, -1)

who_region<- data.frame(read_xlsx("/Users/daniellechaney/Desktop/MPH Documents/Thesis/Data Analysis/WHO Region.xlsx"))

##Join datasets with region dataset
de_by_year <- de_by_year %>% left_join(who_region)
oe_by_year <- oe_by_year %>% left_join(who_region)
ie_by_year <- ie_by_year %>% left_join(who_region)

##Distribution of indirect effects
hist(ie_by_year$ie)

##Add timepoint column
ie_by_year <- ie_by_year %>% mutate(timepoint = rep(c(-5, 0, 1:12),times=112))

#EMR
EMR <- ie_by_year %>% filter(region == "EMR")
p1 <- ggplot(EMR, aes(x = timepoint, y = ie, color=country)) + 
  geom_line() +
  xlab("Years Since Vaccine Introduction")+
  ylab("% Reduction")+
  ggtitle("Reduction in Deaths due to IE (EMR)")+
  theme (legend.position = "none",
         legend.title= element_blank())

ggplot(EMR, aes(x=timepoint))+
  facet_wrap(~country)+
  geom_line(aes(y=oe), color = "black")+
  geom_line(aes(y=de), color="black", linetype="twodash")+
  xlab("Years Since Vaccine Introduction")+
  ylab("% Reduction in Deaths")



#AFR
AFR <- ie_by_year %>% filter(region == "AFR")
p2 <- ggplot(AFR, aes(x = timepoint, y = ie, color=country)) + 
  geom_line() +
  xlab("Years Since Vaccine Introduction")+
  ylab("% Reduction")+
  ggtitle("Reduction in Deaths due to IE (AFR)")+
  theme (legend.position = "none",
         legend.title= element_blank())

ggplot(AFR, aes(x=timepoint))+
  facet_wrap(~country)+
  geom_line(aes(y=oe), color = "black")+
  geom_line(aes(y=de), color="black", linetype="twodash")+
  xlab("Years Since Vaccine Introduction")+
  ylab("% Reduction in Deaths")

#EUR
EUR <- ie_by_year %>% filter(region == "EUR")
p3 <- ggplot(EUR, aes(x = timepoint, y = ie, color=country)) + 
  geom_line() +
  xlab("Years Since Vaccine Introduction")+
  ylab("% Reduction")+
  ggtitle("Reduction in Deaths due to IE (EUR)")+
  theme (legend.position = "none",
         legend.title= element_blank())

ggplot(EUR, aes(x=timepoint))+
  facet_wrap(~country)+
  geom_line(aes(y=oe), color = "black")+
  geom_line(aes(y=de), color="black", linetype="twodash")+
  xlab("Years Since Vaccine Introduction")+
  ylab("% Reduction in Deaths")

#SEAR
SEAR <- ie_by_year %>% filter(region == "SEAR")
p4 <- ggplot(SEAR, aes(x = timepoint, y = ie, color=country)) + 
  geom_line() +
  xlab("Years Since Vaccine Introduction")+
  ylab("% Reduction")+
  ggtitle("Reduction in Deaths due to IE (SEAR)")+
  theme (legend.position = "none",
         legend.title= element_blank())

ggplot(SEAR, aes(x=timepoint))+
  facet_wrap(~country)+
  geom_line(aes(y=oe), color = "black")+
  geom_line(aes(y=de), color="black", linetype="twodash")+
  xlab("Years Since Vaccine Introduction")+
  ylab("% Reduction in Deaths")

#AMR
AMR <- ie_by_year %>% filter(region == "AMR")
p5 <- ggplot(AMR, aes(x = timepoint, y = ie, color=country)) + 
  geom_line() +
  xlab("Years Since Vaccine Introduction")+
  ylab("% Reduction")+
  ggtitle("Reduction in Deaths due to IE (AMR)")+
  theme (legend.position = "none",
         legend.title= element_blank())

ggplot(AMR, aes(x=timepoint))+
  facet_wrap(~country)+
  geom_line(aes(y=oe), color = "black")+
  geom_line(aes(y=de), color="black", linetype="twodash")+
  xlab("Years Since Vaccine Introduction")+
  ylab("% Reduction in Deaths")

#WPR
WPR <- ie_by_year %>% filter(region == "WPR")
p6 <- ggplot(WPR, aes(x = timepoint, y = ie, color=country)) + 
  geom_line() +
  xlab("Years Since Vaccine Introduction")+
  ylab("% Reduction")+
  ggtitle("Reduction in Deaths due to IE (WPR)")+
  theme (legend.position = "none",
         legend.title= element_blank())

ggplot(WPR, aes(x=timepoint))+
  facet_wrap(~country)+
  geom_line(aes(y=oe), color = "black")+
  geom_line(aes(y=de), color="black", linetype="twodash")+
  xlab("Years Since Vaccine Introduction")+
  ylab("% Reduction in Deaths")

#All regions, IEs
(p1+p2+p3)/(p4+p5+p6)


##### PLOTS FOR MANUSCRIPT ######
year0_ie <- ie_by_year %>% filter(timepoint==0)
year5_ie <- ie_by_year %>% filter(timepoint==5)
year8_ie <- ie_by_year %>% filter(timepoint==8)

##Violin Plots by region
p <- ggplot(year8_ie) +
  geom_density_ridges(aes(x = ie,
                          y = region, group = region, fill = region), alpha = 0.5,
                      scale = 0.9) +
  theme_bw() +
  xlab('Proportion of Deaths Adverted') +
  ylab('Region')+
  theme(legend.position="none")
  



##Direct, Overall Curves
median_oe <- aggregate(oe ~ region + timepoint, data=ie_by_year, FUN = median)
median_de <- aggregate(de ~ region + timepoint, data=ie_by_year, FUN = median)
median_ie <- aggregate(ie ~ region + timepoint, data=ie_by_year, FUN = median)
median_effect<- median_oe %>% left_join(median_de) 

afr_point<- data.frame("AFR", -0.1, 0, 0)
names(afr_point)<- c("region", "timepoint", "oe", "de")

amr_point<- data.frame("AMR", -0.1, 0, 0)
names(amr_point)<- c("region", "timepoint", "oe", "de")

emr_point<- data.frame("EMR", -0.1, 0, 0)
names(emr_point)<- c("region", "timepoint", "oe", "de")

eur_point<- data.frame("EUR", -0.1, 0, 0)
names(eur_point)<- c("region", "timepoint", "oe", "de")

sear_point<- data.frame("SEAR", -0.1, 0, 0)
names(sear_point)<- c("region", "timepoint", "oe", "de")

wpr_point<- data.frame("WPR", -0.1, 0, 0)
names(wpr_point)<- c("region", "timepoint", "oe", "de")

df <- rbind(median_effect, afr_point, amr_point, emr_point, eur_point, sear_point, wpr_point)

region.labs <- c("African", "Americas", "Eastern Mediterranean", "European", "South-East Asian", "Western Pacific")
names(region.labs) <- c("AFR", "AMR", "EMR", "EUR", "SEAR", "WPR")

ggplot(df, aes(x=timepoint))+
  facet_wrap(~region, labeller=labeller(region=region.labs))+
  geom_line(aes(y=oe), color = "black")+
  geom_line(aes(y=de), color="black", linetype="twodash")+
  xlab("Years Since Vaccine Introduction")+
  ylab("Proportion of Deaths Adverted") + 
  theme(strip.background = element_blank())+
  theme(strip.text=element_text(size=8, face="bold"))+
  theme(axis.title=element_text(size=8))+
  theme_bw()
 

median_de %>% filter(region == "WPR", timepoint == 8)
median_oe %>% filter(region == "WPR", timepoint == 8)
median_ie %>% filter(region == "WPR", timepoint == 8)

max <- median_ie %>% filter(region == "WPR")
which.max(max$ie)

which.max(ie_by_year$ie)


