###proportion symtpomatic from Velazquez
alpha1 = 0.30 	#proportion symptomatic at 1st infection (0.30)
alpha2 = 0.28		#...2nd infection (0.28)
alpha3 = 0.18		#...3rd infection (0.18)
alpha4 = 0.21		#...4th infection (0.21)

rel_inf = 0.25		

pathi <- '/Users/akraay/Dropbox/VIMC/demog/' 
input_data <- "VIMC_Parameters_contdictcbr_2019-11-25.csv"
coverage_data <- "coverage_201910gavi-4_rota-routine-default.csv"
coverage_data2 <- "country_ie_by_year.csv"
pathi2<-'/Users/akraay/Dropbox/ChaneyThesis/'
country_params<-read.csv(paste0(pathi,input_data))
coverage_params<-read.csv(paste0(pathi,coverage_data))
ie_dat<-read.csv(paste0(pathi2,coverage_data2))

#Get year of introduction
startv<-character(112)
for(i in 1:112){
  startv[i]<-names(coverage_params[i,12:ncol(coverage_params)][which(coverage_params[i,12:ncol(coverage_params)]>0)])[1]
}
startc<-substr(startv, start=10, stop=13)

startn<-as.numeric(as.character(startc))
#change to numeric and make into data frame, then merge below to make figures as years since vaccine introduction

intro.dat<-data.frame(country_code=coverage_params$country_code, vaccine_intro=startn)
negIE<-c(c('FJI', 'FSM', 'KIR', 'MHL', 'PHL', 'TON', 'TUV', 'VUT', 'WSM', 'CIV', 'GIN', 'GNB', 'NER', 'SEN', 'TCD', 'GTM', 'VEN', 'SYR'))
plotlists<-list(length(negIE))

IEvals<-read.csv("/Users/akraay/Dropbox/ChaneyThesis/country_ie_by_year.csv")
IEneg<-subset(IEvals, IEvals$country %in% unique(negIE))
#names(IEneg)[2]<-'country_code'
#IEneg2<-merge(IEneg, intro.dat, by='country_code')
intro.dat2<-intro.dat
names(intro.dat2)[1]<-'country'
IEneg2<-merge(IEneg, intro.dat2, by='country')
IEneg2$rel_yr<-IEneg2$year-IEneg2$vaccine_intro
IEneg2<-subset(IEneg2, IEneg2$rel_yr>=0)
IEneg2<-subset(IEneg2, IEneg2$ie<0)
names(IEneg2)[1]<-'country_code'
IEneg3<-merge(IEneg2, coverage_params[,c(6, 7)], by='country_code')

IEneg3$ordervec<-NA
for(i in 1:length(negIE)){
  IEneg3$ordervec[which(IEneg3$country_code==negIE[i])]<-i
}

IEneg4<-IEneg3[order(IEneg3$ordervec),]
IEpick1<-ddply(IEneg4, .(ordervec), function(x) head(x,1))
for (i in 1:length(negIE)){
  patho<-'/Users/akraay/Dropbox/ChaneyThesis/RotaDCM-legacy/results/'
  type<-'ODE_res_routine'
  est<-'central.csv'
  dat<-read.csv(paste0(patho,type,negIE[i],est))
  dat$time_yrs<-dat$time/365
  dat.end<-subset(dat, dat$time_yrs>=2000)
  dat.end$N<-dat.end$M_0_1+dat.end$S1_0_1+dat.end$I1_0_1+dat.end$S2_0_1+dat.end$I2_0_1+dat.end$S3_0_1+dat.end$I3_0_1+dat.end$S4_0_1+dat.end$I4_0_1+dat.end$R_0_1+
    dat.end$M_2_3+dat.end$S1_2_3+dat.end$I1_2_3+dat.end$S2_2_3+dat.end$I2_2_3+dat.end$S3_2_3+dat.end$I3_2_3+dat.end$S4_2_3+dat.end$I4_2_3+dat.end$R_2_3+dat.end$V_2_3+
    dat.end$VS1_2_3+dat.end$VI1_2_3+dat.end$VS2_2_3+dat.end$VI2_2_3+dat.end$VS3_2_3+dat.end$VI3_2_3+dat.end$VS4_2_3+dat.end$VI4_2_3+dat.end$VR_2_3+
    dat.end$M_4_5+dat.end$S1_4_5+dat.end$I1_4_5+dat.end$S2_4_5+dat.end$I2_4_5+dat.end$S3_4_5+dat.end$I3_4_5+dat.end$S4_4_5+dat.end$I4_4_5+dat.end$R_4_5+dat.end$V_4_5+
    dat.end$VS1_4_5+dat.end$VI1_4_5+dat.end$VS2_4_5+dat.end$VI2_4_5+dat.end$VS3_4_5+dat.end$VI3_4_5+dat.end$VS4_4_5+dat.end$VI4_4_5+dat.end$VR_4_5+
    dat.end$M_6_11+dat.end$S1_6_11+dat.end$I1_6_11+dat.end$S2_6_11+dat.end$I2_6_11+dat.end$S3_6_11+dat.end$I3_6_11+dat.end$S4_6_11+dat.end$I4_6_11+dat.end$R_6_11+dat.end$V_6_11+
    dat.end$VS1_6_11+dat.end$VI1_6_11+dat.end$VS2_6_11+dat.end$VI2_6_11+dat.end$VS3_6_11+dat.end$VI3_6_11+dat.end$VS4_6_11+dat.end$VI4_6_11+dat.end$VR_6_11+
    dat.end$S1_12_23+dat.end$I1_12_23+dat.end$S2_12_23+dat.end$I2_12_23+dat.end$S3_12_23+dat.end$I3_12_23+dat.end$S4_12_23+dat.end$I4_12_23+dat.end$R_12_23+dat.end$V_12_23+
    dat.end$VS1_12_23+dat.end$VI1_12_23+dat.end$VS2_12_23+dat.end$VI2_12_23+dat.end$VS3_12_23+dat.end$VI3_12_23+dat.end$VS4_12_23+dat.end$VI4_12_23+dat.end$VR_12_23+
    dat.end$S1_24_35+dat.end$I1_24_35+dat.end$S2_24_35+dat.end$I2_24_35+dat.end$S3_24_35+dat.end$I3_24_35+dat.end$S4_24_35+dat.end$I4_24_35+dat.end$R_24_35+dat.end$V_24_35+
    dat.end$VS1_24_35+dat.end$VI1_24_35+dat.end$VS2_24_35+dat.end$VI2_24_35+dat.end$VS3_24_35+dat.end$VI3_24_35+dat.end$VS4_24_35+dat.end$VI4_24_35+dat.end$VR_24_35+
    dat.end$S1_36_47+dat.end$I1_36_47+dat.end$S2_36_47+dat.end$I2_36_47+dat.end$S3_36_47+dat.end$I3_36_47+dat.end$S4_36_47+dat.end$I4_36_47+dat.end$R_36_47+dat.end$V_36_47+
    dat.end$VS1_36_47+dat.end$VI1_36_47+dat.end$VS2_36_47+dat.end$VI2_36_47+dat.end$VS3_36_47+dat.end$VI3_36_47+dat.end$VS4_36_47+dat.end$VI4_36_47+dat.end$VR_36_47+
    dat.end$S1_48_59+dat.end$I1_48_59+dat.end$S2_48_59+dat.end$I2_48_59+dat.end$S3_48_59+dat.end$I3_48_59+dat.end$S4_48_59+dat.end$I4_48_59+dat.end$R_48_59+dat.end$V_48_59+
    dat.end$VS1_48_59+dat.end$VI1_48_59+dat.end$VS2_48_59+dat.end$VI2_48_59+dat.end$VS3_48_59+dat.end$VI3_48_59+dat.end$VS4_48_59+dat.end$VI4_48_59+dat.end$VR_48_59+
    dat.end$S1_60_up+dat.end$I1_60_up+dat.end$S2_60_up+dat.end$I2_60_up+dat.end$S3_60_up+dat.end$I3_60_up+dat.end$S4_60_up+dat.end$I4_60_up+dat.end$R_60_up+dat.end$V_60_up+
    dat.end$VS1_60_up+dat.end$VI1_60_up+dat.end$VS2_60_up+dat.end$VI2_60_up+dat.end$VS3_60_up+dat.end$VI3_60_up+dat.end$VS4_60_up+dat.end$VI4_60_up+dat.end$VR_60_up

  dat.end$infec<-(alpha1*dat.end$I1_0_1 + (alpha2*dat.end$I2_0_1+alpha3*dat.end$I3_0_1+alpha4*dat.end$I4_0_1)*rel_inf)+
    (alpha1*dat.end$I1_2_3 + (alpha2*dat.end$I2_2_3+alpha3*dat.end$I3_2_3+alpha4*dat.end$I4_2_3)*rel_inf)+
    (alpha1*dat.end$VI1_2_3 + (alpha2*dat.end$VI2_2_3+alpha3*dat.end$VI3_2_3+alpha4*dat.end$VI4_2_3)*rel_inf)+
    (alpha1*dat.end$I1_4_5 + (alpha2*dat.end$I2_4_5+alpha3*dat.end$I3_4_5+alpha4*dat.end$I4_4_5)*rel_inf)+
    (alpha1*dat.end$VI1_4_5 + (alpha2*dat.end$VI2_4_5+alpha3*dat.end$VI3_4_5+alpha4*dat.end$VI4_4_5)*rel_inf)+
    (alpha1*dat.end$I1_6_11 + (alpha2*dat.end$I2_6_11+alpha3*dat.end$I3_6_11+alpha4*dat.end$I4_6_11)*rel_inf)+
    (alpha1*dat.end$VI1_6_11 + (alpha2*dat.end$VI2_6_11+alpha3*dat.end$VI3_6_11+alpha4*dat.end$VI4_6_11)*rel_inf)+
    (alpha1*dat.end$I1_12_23 + (alpha2*dat.end$I2_12_23+alpha3*dat.end$I3_12_23+alpha4*dat.end$I4_12_23)*rel_inf)+
    (alpha1*dat.end$VI1_12_23 + (alpha2*dat.end$VI2_12_23+alpha3*dat.end$VI3_12_23+alpha4*dat.end$VI4_12_23)*rel_inf)+
    (alpha1*dat.end$I1_24_35 + (alpha2*dat.end$I2_24_35+alpha3*dat.end$I3_24_35+alpha4*dat.end$I4_24_35)*rel_inf)+
    (alpha1*dat.end$VI1_24_35 + (alpha2*dat.end$VI2_24_35+alpha3*dat.end$VI3_24_35+alpha4*dat.end$VI4_24_35)*rel_inf)+
    (alpha1*dat.end$I1_36_47 + (alpha2*dat.end$I2_36_47+alpha3*dat.end$I3_36_47+alpha4*dat.end$I4_36_47)*rel_inf)+
    (alpha1*dat.end$VI1_36_47 + (alpha2*dat.end$VI2_36_47+alpha3*dat.end$VI3_36_47+alpha4*dat.end$VI4_36_47)*rel_inf)+
    (alpha1*dat.end$I1_48_59 + (alpha2*dat.end$I2_48_59+alpha3*dat.end$I3_48_59+alpha4*dat.end$I4_48_59)*rel_inf)+
    (alpha1*dat.end$VI1_48_59 + (alpha2*dat.end$VI2_48_59+alpha3*dat.end$VI3_48_59+alpha4*dat.end$VI4_48_59)*rel_inf)+
    (alpha1*dat.end$I1_60_up + (alpha2*dat.end$I2_60_up+alpha3*dat.end$I3_60_up+alpha4*dat.end$I4_60_up)*rel_inf)+
    (alpha1*dat.end$VI1_60_up + (alpha2*dat.end$VI2_60_up+alpha3*dat.end$VI3_60_up+alpha4*dat.end$VI4_60_up)*rel_inf)
  
  dat.end$susceptible<-dat.end$S1_0_1+dat.end$S2_0_1+dat.end$S3_0_1+dat.end$S4_0_1+
    dat.end$S1_2_3+dat.end$S2_2_3+dat.end$S3_2_3+dat.end$S4_2_3+dat.end$VS1_2_3+dat.end$VS2_2_3+dat.end$VS3_2_3+dat.end$VS4_2_3+
    dat.end$S1_4_5+dat.end$S2_4_5+dat.end$S3_4_5+dat.end$S4_4_5+dat.end$VS1_4_5+dat.end$VS2_4_5+dat.end$VS3_4_5+dat.end$VS4_4_5+
    dat.end$S1_6_11+dat.end$S2_6_11+dat.end$S3_6_11+dat.end$S4_6_11+dat.end$VS1_6_11+dat.end$VS2_6_11+dat.end$VS3_6_11+dat.end$VS4_6_11+
    dat.end$S1_12_23+dat.end$S2_12_23+dat.end$S3_12_23+dat.end$S4_12_23+
    dat.end$VS1_12_23+dat.end$VS2_12_23+dat.end$VS3_12_23+dat.end$VS4_12_23+
    dat.end$S1_24_35+dat.end$S2_24_35+dat.end$S3_24_35+dat.end$S4_24_35+
    dat.end$VS1_24_35+dat.end$VS2_24_35+dat.end$VS3_24_35+dat.end$VS4_24_35+
    dat.end$S1_36_47+dat.end$S2_36_47+dat.end$S3_36_47+dat.end$S4_36_47+
    dat.end$VS1_36_47+dat.end$VS2_36_47+dat.end$VS3_36_47+dat.end$VS4_36_47+
    dat.end$S1_48_59+dat.end$S2_48_59+dat.end$S3_48_59+dat.end$S4_48_59+
    dat.end$VS1_48_59+dat.end$VS2_48_59+dat.end$VS3_48_59+dat.end$VS4_48_59+
    dat.end$S1_60_up+dat.end$S2_60_up+dat.end$S3_60_up+dat.end$S4_60_up+
    dat.end$VS1_60_up+dat.end$VS2_60_up+dat.end$VS3_60_up+dat.end$VS4_60_up
  
  param.temp<-subset(country_params, country_params$country_code==negIE[i])
  
  q1<-param.temp$life_expectancy/(param.temp$VIMCPredMeanAgeWk/52)
  
  dat.end$FOI<-q1*(dat.end$infec/dat.end$N)
  scale.var1<-max(dat.end$susceptible)/max(dat.end$FOI)
  dat.end$country_code<-negIE[i]
  dat.end2<-merge(dat.end, intro.dat, by='country_code')
  dat.end2$rel.yr<-dat.end2$time_yrs-dat.end2$vaccine_intro
  
  dat.end3<-subset(dat.end2, dat.end2$rel.yr>-5 & dat.end2$rel.yr<=12)
  #p1<-ggplot(data=dat.end3, aes(x=rel.yr))+geom_line(aes(y=FOI))+geom_line(aes(y=susceptible/scale.var1), linetype='dashed')+
  #  scale_y_continuous(name="Force of Infection", 
  #                     sec.axis=sec_axis(~.*scale.var1, name='Susceptible population'))+theme_bw()+labs(x='year')+
  #ggtitle(negIE[i])
  
  scale.var2<-max(dat.end3$susceptible*100/dat.end3$N)/max(dat.end3$FOI)
  negbar<-subset(IEneg2, IEneg2$country==negIE[i])
  p2<-ggplot(data=dat.end3, aes(x=rel.yr))+geom_line(aes(y=FOI))+geom_line(aes(y=(susceptible*100/N)/scale.var2), linetype='dashed')+
    scale_y_continuous(name="Force of Infection", 
                       sec.axis=sec_axis(~.*scale.var2, name='Susceptible population (%)'))+
    theme_bw()+labs(x='year')+#geom_rect(aes(xmin=min(negbar$rel_yr), xmax=max(negbar$rel_yr)), ymin=min(dat.end3$FOI-0.001), ymax=max(dat.end3$FOI+0.01), alpha=0.01)+
    ggtitle(IEpick1$country[i])
  
  plotlists[[i]]<-p2
  
  #pathf<-'/Users/akraay/Dropbox/Chaney Thesis/RotaDCM-legacy/figs/'
  
  #ggsave(paste0(pathf, negIE[i], '.pdf'), p1)
  
  #ggsave(paste0(pathf, negIE[i], 'proportionS', '.pdf'), p2)
  
}

require(ggpubr)
require(gridExtra)
pcombine<-grid.arrange(plotlists[[1]], plotlists[[2]], plotlists[[3]], plotlists[[4]],
             plotlists[[5]], plotlists[[6]], plotlists[[7]], plotlists[[8]],
             plotlists[[9]], plotlists[[10]], plotlists[[11]], plotlists[[12]],
             plotlists[[13]], plotlists[[14]], plotlists[[15]], plotlists[[16]],
             plotlists[[17]], plotlists[[18]], ncol=4)

ggsave('/Users/akraay/Dropbox/ChaneyThesis/RotaDCM-legacy/figs/IE_combined_jan2023.pdf', pcombine, 
       w=14, h=12)



  patho<-'/Users/akraay/Dropbox/ChaneyThesis/RotaDCM-legacy/results/'
  type<-'ODE_res_routine'
  est<-'central.csv'
  dat<-read.csv(paste0(patho,type,'AGO',est))
  dat$time_yrs<-dat$time/365
  dat.end<-subset(dat, dat$time_yrs>=2000)
  dat.end$N<-dat.end$M_0_1+dat.end$S1_0_1+dat.end$I1_0_1+dat.end$S2_0_1+dat.end$I2_0_1+dat.end$S3_0_1+dat.end$I3_0_1+dat.end$S4_0_1+dat.end$I4_0_1+dat.end$R_0_1+
    dat.end$M_2_3+dat.end$S1_2_3+dat.end$I1_2_3+dat.end$S2_2_3+dat.end$I2_2_3+dat.end$S3_2_3+dat.end$I3_2_3+dat.end$S4_2_3+dat.end$I4_2_3+dat.end$R_2_3+dat.end$V_2_3+
    dat.end$VS1_2_3+dat.end$VI1_2_3+dat.end$VS2_2_3+dat.end$VI2_2_3+dat.end$VS3_2_3+dat.end$VI3_2_3+dat.end$VS4_2_3+dat.end$VI4_2_3+dat.end$VR_2_3+
    dat.end$M_4_5+dat.end$S1_4_5+dat.end$I1_4_5+dat.end$S2_4_5+dat.end$I2_4_5+dat.end$S3_4_5+dat.end$I3_4_5+dat.end$S4_4_5+dat.end$I4_4_5+dat.end$R_4_5+dat.end$V_4_5+
    dat.end$VS1_4_5+dat.end$VI1_4_5+dat.end$VS2_4_5+dat.end$VI2_4_5+dat.end$VS3_4_5+dat.end$VI3_4_5+dat.end$VS4_4_5+dat.end$VI4_4_5+dat.end$VR_4_5+
    dat.end$M_6_11+dat.end$S1_6_11+dat.end$I1_6_11+dat.end$S2_6_11+dat.end$I2_6_11+dat.end$S3_6_11+dat.end$I3_6_11+dat.end$S4_6_11+dat.end$I4_6_11+dat.end$R_6_11+dat.end$V_6_11+
    dat.end$VS1_6_11+dat.end$VI1_6_11+dat.end$VS2_6_11+dat.end$VI2_6_11+dat.end$VS3_6_11+dat.end$VI3_6_11+dat.end$VS4_6_11+dat.end$VI4_6_11+dat.end$VR_6_11+
    dat.end$S1_12_23+dat.end$I1_12_23+dat.end$S2_12_23+dat.end$I2_12_23+dat.end$S3_12_23+dat.end$I3_12_23+dat.end$S4_12_23+dat.end$I4_12_23+dat.end$R_12_23+dat.end$V_12_23+
    dat.end$VS1_12_23+dat.end$VI1_12_23+dat.end$VS2_12_23+dat.end$VI2_12_23+dat.end$VS3_12_23+dat.end$VI3_12_23+dat.end$VS4_12_23+dat.end$VI4_12_23+dat.end$VR_12_23+
    dat.end$S1_24_35+dat.end$I1_24_35+dat.end$S2_24_35+dat.end$I2_24_35+dat.end$S3_24_35+dat.end$I3_24_35+dat.end$S4_24_35+dat.end$I4_24_35+dat.end$R_24_35+dat.end$V_24_35+
    dat.end$VS1_24_35+dat.end$VI1_24_35+dat.end$VS2_24_35+dat.end$VI2_24_35+dat.end$VS3_24_35+dat.end$VI3_24_35+dat.end$VS4_24_35+dat.end$VI4_24_35+dat.end$VR_24_35+
    dat.end$S1_36_47+dat.end$I1_36_47+dat.end$S2_36_47+dat.end$I2_36_47+dat.end$S3_36_47+dat.end$I3_36_47+dat.end$S4_36_47+dat.end$I4_36_47+dat.end$R_36_47+dat.end$V_36_47+
    dat.end$VS1_36_47+dat.end$VI1_36_47+dat.end$VS2_36_47+dat.end$VI2_36_47+dat.end$VS3_36_47+dat.end$VI3_36_47+dat.end$VS4_36_47+dat.end$VI4_36_47+dat.end$VR_36_47+
    dat.end$S1_48_59+dat.end$I1_48_59+dat.end$S2_48_59+dat.end$I2_48_59+dat.end$S3_48_59+dat.end$I3_48_59+dat.end$S4_48_59+dat.end$I4_48_59+dat.end$R_48_59+dat.end$V_48_59+
    dat.end$VS1_48_59+dat.end$VI1_48_59+dat.end$VS2_48_59+dat.end$VI2_48_59+dat.end$VS3_48_59+dat.end$VI3_48_59+dat.end$VS4_48_59+dat.end$VI4_48_59+dat.end$VR_48_59+
    dat.end$S1_60_up+dat.end$I1_60_up+dat.end$S2_60_up+dat.end$I2_60_up+dat.end$S3_60_up+dat.end$I3_60_up+dat.end$S4_60_up+dat.end$I4_60_up+dat.end$R_60_up+dat.end$V_60_up+
    dat.end$VS1_60_up+dat.end$VI1_60_up+dat.end$VS2_60_up+dat.end$VI2_60_up+dat.end$VS3_60_up+dat.end$VI3_60_up+dat.end$VS4_60_up+dat.end$VI4_60_up+dat.end$VR_60_up
  
  dat.end$infec<-(alpha1*dat.end$I1_0_1 + (alpha2*dat.end$I2_0_1+alpha3*dat.end$I3_0_1+alpha4*dat.end$I4_0_1)*rel_inf)+
    (alpha1*dat.end$I1_2_3 + (alpha2*dat.end$I2_2_3+alpha3*dat.end$I3_2_3+alpha4*dat.end$I4_2_3)*rel_inf)+
    (alpha1*dat.end$VI1_2_3 + (alpha2*dat.end$VI2_2_3+alpha3*dat.end$VI3_2_3+alpha4*dat.end$VI4_2_3)*rel_inf)+
    (alpha1*dat.end$I1_4_5 + (alpha2*dat.end$I2_4_5+alpha3*dat.end$I3_4_5+alpha4*dat.end$I4_4_5)*rel_inf)+
    (alpha1*dat.end$VI1_4_5 + (alpha2*dat.end$VI2_4_5+alpha3*dat.end$VI3_4_5+alpha4*dat.end$VI4_4_5)*rel_inf)+
    (alpha1*dat.end$I1_6_11 + (alpha2*dat.end$I2_6_11+alpha3*dat.end$I3_6_11+alpha4*dat.end$I4_6_11)*rel_inf)+
    (alpha1*dat.end$VI1_6_11 + (alpha2*dat.end$VI2_6_11+alpha3*dat.end$VI3_6_11+alpha4*dat.end$VI4_6_11)*rel_inf)+
    (alpha1*dat.end$I1_12_23 + (alpha2*dat.end$I2_12_23+alpha3*dat.end$I3_12_23+alpha4*dat.end$I4_12_23)*rel_inf)+
    (alpha1*dat.end$VI1_12_23 + (alpha2*dat.end$VI2_12_23+alpha3*dat.end$VI3_12_23+alpha4*dat.end$VI4_12_23)*rel_inf)+
    (alpha1*dat.end$I1_24_35 + (alpha2*dat.end$I2_24_35+alpha3*dat.end$I3_24_35+alpha4*dat.end$I4_24_35)*rel_inf)+
    (alpha1*dat.end$VI1_24_35 + (alpha2*dat.end$VI2_24_35+alpha3*dat.end$VI3_24_35+alpha4*dat.end$VI4_24_35)*rel_inf)+
    (alpha1*dat.end$I1_36_47 + (alpha2*dat.end$I2_36_47+alpha3*dat.end$I3_36_47+alpha4*dat.end$I4_36_47)*rel_inf)+
    (alpha1*dat.end$VI1_36_47 + (alpha2*dat.end$VI2_36_47+alpha3*dat.end$VI3_36_47+alpha4*dat.end$VI4_36_47)*rel_inf)+
    (alpha1*dat.end$I1_48_59 + (alpha2*dat.end$I2_48_59+alpha3*dat.end$I3_48_59+alpha4*dat.end$I4_48_59)*rel_inf)+
    (alpha1*dat.end$VI1_48_59 + (alpha2*dat.end$VI2_48_59+alpha3*dat.end$VI3_48_59+alpha4*dat.end$VI4_48_59)*rel_inf)+
    (alpha1*dat.end$I1_60_up + (alpha2*dat.end$I2_60_up+alpha3*dat.end$I3_60_up+alpha4*dat.end$I4_60_up)*rel_inf)+
    (alpha1*dat.end$VI1_60_up + (alpha2*dat.end$VI2_60_up+alpha3*dat.end$VI3_60_up+alpha4*dat.end$VI4_60_up)*rel_inf)
  
  dat.end$susceptible<-dat.end$S1_0_1+dat.end$S2_0_1+dat.end$S3_0_1+dat.end$S4_0_1+
    dat.end$S1_2_3+dat.end$S2_2_3+dat.end$S3_2_3+dat.end$S4_2_3+dat.end$VS1_2_3+dat.end$VS2_2_3+dat.end$VS3_2_3+dat.end$VS4_2_3+
    dat.end$S1_4_5+dat.end$S2_4_5+dat.end$S3_4_5+dat.end$S4_4_5+dat.end$VS1_4_5+dat.end$VS2_4_5+dat.end$VS3_4_5+dat.end$VS4_4_5+
    dat.end$S1_6_11+dat.end$S2_6_11+dat.end$S3_6_11+dat.end$S4_6_11+dat.end$VS1_6_11+dat.end$VS2_6_11+dat.end$VS3_6_11+dat.end$VS4_6_11+
    dat.end$S1_12_23+dat.end$S2_12_23+dat.end$S3_12_23+dat.end$S4_12_23+
    dat.end$VS1_12_23+dat.end$VS2_12_23+dat.end$VS3_12_23+dat.end$VS4_12_23+
    dat.end$S1_24_35+dat.end$S2_24_35+dat.end$S3_24_35+dat.end$S4_24_35+
    dat.end$VS1_24_35+dat.end$VS2_24_35+dat.end$VS3_24_35+dat.end$VS4_24_35+
    dat.end$S1_36_47+dat.end$S2_36_47+dat.end$S3_36_47+dat.end$S4_36_47+
    dat.end$VS1_36_47+dat.end$VS2_36_47+dat.end$VS3_36_47+dat.end$VS4_36_47+
    dat.end$S1_48_59+dat.end$S2_48_59+dat.end$S3_48_59+dat.end$S4_48_59+
    dat.end$VS1_48_59+dat.end$VS2_48_59+dat.end$VS3_48_59+dat.end$VS4_48_59+
    dat.end$S1_60_up+dat.end$S2_60_up+dat.end$S3_60_up+dat.end$S4_60_up+
    dat.end$VS1_60_up+dat.end$VS2_60_up+dat.end$VS3_60_up+dat.end$VS4_60_up
  
  param.temp<-subset(country_params, country_params$country_code=='AGO')
  
  q1<-param.temp$life_expectancy/(param.temp$VIMCPredMeanAgeWk/52)
  
  dat.end$FOI<-q1*(dat.end$infec/dat.end$N)
  scale.var1<-max(dat.end$susceptible)/max(dat.end$FOI)
  dat.end$country_code<-negIE[i]
  dat.end2<-merge(dat.end, intro.dat, by='country_code')
  dat.end2$rel.yr<-dat.end2$time_yrs-dat.end2$vaccine_intro
  
  dat.end3<-subset(dat.end2, dat.end2$rel.yr>-5 & dat.end2$rel.yr<=12)
  #p1<-ggplot(data=dat.end3, aes(x=rel.yr))+geom_line(aes(y=FOI))+geom_line(aes(y=susceptible/scale.var1), linetype='dashed')+
  #  scale_y_continuous(name="Force of Infection", 
  #                     sec.axis=sec_axis(~.*scale.var1, name='Susceptible population'))+theme_bw()+labs(x='year')+
  #ggtitle(negIE[i])
  
  scale.var2<-max(dat.end3$susceptible*100/dat.end3$N)/max(dat.end3$FOI)
  
ggplot(data=dat.end3, aes(x=rel.yr))+geom_line(aes(y=FOI))+geom_line(aes(y=(susceptible*100/N)/scale.var2), linetype='dashed')+
    scale_y_continuous(name="Force of Infection", 
                       sec.axis=sec_axis(~.*scale.var2, name='Susceptible population (%)'))+
    theme_bw()+labs(x='year')+
    ggtitle('AGO')
  
