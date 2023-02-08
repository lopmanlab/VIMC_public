#################################################################################
# Rotavirus model for VIMC. Code produces output under vaccination
# Things to add:
# Tracking number of vaccinated individuals for each year
# Fixed FOI post-vaccination
# Options for one dose only
# Coded by: Molly Steele 
# Last updated: 8/1/2019 by Molly Steele
##################################################################################


stat="mean"
wane=250

#designate which RData file to use
#input_data <-"VIMC_Parameters_contcbr2019-08-19.RData";cbrvar="contcbr"
#input_data <-"VIMC_Parameters_catcbr2019-08-19.RData";cbrvar="cattcbr"
#input_data <-"VIMC_Parameters_contdictcbr2019-08-19.RData";cbrvar="contdictcbr"


library(deSolve)  #loads ODE solver package
library(foreach)
library(ggplot2)
library(plyr)
library(dplyr)

#Fixed Parameters   
#demography
#time period is per day

#age rates
age_rate_0_5 = 1/(60)       	#daily rate at which infants age (0-1, 2-3, 4-5 month age groups)
age_rate_6_11 = 1/(30*6)      #daily rate at which infants age (6-11 month age group)
age_rate_12_60 = 1/(365)      #daily rate at which young children age (ages 1-4)



#natural history parameters
#nh_names=c("eta","gamma","rec","rel_inf","P1","P2","P3","alpha1","alpha2","alpha3","alpha4","sev_alpha1","sev_alpha2","sev_alpha3","sev_alpha4","vac_imm")
eta = (1/(26*7))   #rate of loss of maternal antibody protection (1/26 weeks)
gamma = (1/7)  #rate of loss of infectiousness (1/7 days)
rec = (1/(wane*365))       #rate of waning immunity (potentially fit?)

P1 = 0.61				#proportion who DO NOT become immune after 1st infection (0.61)
P2 = 0.48				#...2nd infection  (0.48)
P3 = 0.33				#...3rd infection  (0.33)

rel_inf = 0.25				#relative infectiousness of non-primary infections (25%)

###proportion symtpomatic from Velazquez
alpha1 = 0.30 	#proportion symptomatic at 1st infection (0.30)
alpha2 = 0.28		#...2nd infection (0.28)
alpha3 = 0.18		#...3rd infection (0.18)
alpha4 = 0.21		#...4th infection (0.21)

#severe infections by order of infection# from Velazquez
sev_alpha1 = 0.17	#proportion of symptomatic infections associated with severe disease at nth infection 
sev_alpha2 = 0.23  
sev_alpha3 = 0.24 	
sev_alpha4 = 0.18 



################################################################################################
#
# FUNCTION WITH MODEL: q1 for < 4 mo, q2 for >4 mo, Rep2 for 4-11mo, Rep for all but 4-11 mo
#
###############################################################################################
ystart=1975;yend=2036
times = seq(ystart*365,yend*365, by=mod_time_step) #defining daily timesteps from 1975 (to allow sufficient burn-in) to 2036
times_yr <- floor(times/365)

atolv=1e-6; rtolv=1e-6; #tolerances for ODE solver
#Function for transmission model
RM = function(t, x, p) {
  M_0_1 =x[1];S1_0_1 = x[2];I1_0_1 = x[3];S2_0_1 = x[4];I2_0_1 =x[5];S3_0_1 =x[6];I3_0_1 =x[7];S4_0_1 =x[8];I4_0_1 =x[9];R_0_1 =x[10]
  M_2_3 = x[11];S1_2_3 = x[12];I1_2_3 = x[13];S2_2_3 =x[14];I2_2_3 =x[15];S3_2_3 =x[16];I3_2_3 =x[17]; S4_2_3 =x[18];I4_2_3 =x[19];R_2_3 =x[20];V_2_3=x[21]
  VS1_2_3 = x[22];VI1_2_3 = x[23];VS2_2_3 = x[24];VI2_2_3 = x[25];VS3_2_3 = x[26];VI3_2_3 = x[27];VS4_2_3 = x[28];VI4_2_3 = x[29];VR_2_3 = x[30];
  M_4_5 = x[31];S1_4_5 = x[32];I1_4_5 = x[33];S2_4_5 =x[34];I2_4_5 =x[35];S3_4_5 =x[36];I3_4_5 =x[37];S4_4_5 =x[38];I4_4_5 =x[39];R_4_5 =x[40];V_4_5=x[41]
  VS1_4_5 = x[42];VI1_4_5 = x[43];VS2_4_5 = x[44];VI2_4_5 = x[45];VS3_4_5 = x[46];VI3_4_5 = x[47];VS4_4_5 = x[48];VI4_4_5 = x[49];VR_4_5 = x[50];
  M_6_11 = x[51];S1_6_11 = x[52];I1_6_11 = x[53];S2_6_11 =x[54];I2_6_11 =x[55];S3_6_11 =x[56];I3_6_11 =x[57];S4_6_11 =x[58];I4_6_11 =x[59];R_6_11 =x[60];V_6_11=x[61]
  VS1_6_11 = x[62];VI1_6_11 = x[63];VS2_6_11 = x[64];VI2_6_11 = x[65];VS3_6_11 = x[66];VI3_6_11 = x[67];VS4_6_11 = x[68];VI4_6_11 = x[69];VR_6_11 = x[70];
  S1_12_23 = x[71];I1_12_23 = x[72];S2_12_23 =x[73];I2_12_23 =x[74];S3_12_23 =x[75];I3_12_23 =x[76];S4_12_23 =x[77];I4_12_23 =x[78];R_12_23 =x[79];V_12_23=x[80]
  VS1_12_23 = x[81];VI1_12_23 = x[82];VS2_12_23 = x[83];VI2_12_23 = x[84];VS3_12_23 = x[85];VI3_12_23 = x[86];VS4_12_23 = x[87];VI4_12_23 = x[88];VR_12_23 = x[89];
  S1_24_35 = x[90];I1_24_35 = x[91];S2_24_35 =x[92];I2_24_35 =x[93];S3_24_35 =x[94];I3_24_35 =x[95];S4_24_35 =x[96];I4_24_35 =x[97];R_24_35 =x[98];V_24_35=x[99]
  VS1_24_35 = x[100];VI1_24_35 = x[101];VS2_24_35 = x[102];VI2_24_35 = x[103];VS3_24_35 = x[104];VI3_24_35 = x[105];VS4_24_35 = x[106];VI4_24_35 = x[107];VR_24_35 = x[108];
  S1_36_47 = x[109];I1_36_47 = x[110];S2_36_47 =x[111];I2_36_47 =x[112];S3_36_47 =x[113];I3_36_47 =x[114];S4_36_47 =x[115];I4_36_47 =x[116];R_36_47 =x[117];V_36_47=x[118]
  VS1_36_47 = x[119];VI1_36_47 = x[120];VS2_36_47 = x[121];VI2_36_47 = x[122];VS3_36_47 = x[123];VI3_36_47 = x[124];VS4_36_47 = x[125];VI4_36_47 = x[126];VR_36_47 = x[127];
  S1_48_59 = x[128];I1_48_59 = x[129];S2_48_59 =x[130];I2_48_59 =x[131];S3_48_59 =x[132];I3_48_59 =x[133];S4_48_59 =x[134];I4_48_59 =x[135];R_48_59 =x[136];V_48_59=x[137]
  VS1_48_59 = x[138];VI1_48_59 = x[139];VS2_48_59 = x[140];VI2_48_59 = x[141];VS3_48_59 = x[142];VI3_48_59 = x[143];VS4_48_59 = x[144];VI4_48_59 = x[145];VR_48_59 = x[146];
  S1_60_up = x[147];I1_60_up = x[148];S2_60_up =x[149];I2_60_up =x[150];S3_60_up =x[151];I3_60_up =x[152];S4_60_up =x[153];I4_60_up =x[154];R_60_up =x[155];V_60_up=x[156]
  VS1_60_up = x[157];VI1_60_up = x[158];VS2_60_up = x[159];VI2_60_up = x[160];VS3_60_up = x[161];VI3_60_up = x[162];VS4_60_up = x[163];VI4_60_up = x[164];VR_60_up = x[165];
  
  
  N = M_0_1+S1_0_1+I1_0_1+S2_0_1+I2_0_1+S3_0_1+I3_0_1+S4_0_1+I4_0_1+R_0_1+
    M_2_3+S1_2_3+I1_2_3+S2_2_3+I2_2_3+S3_2_3+I3_2_3+S4_2_3+I4_2_3+R_2_3+V_2_3+
    VS1_2_3+VI1_2_3+VS2_2_3+VI2_2_3+VS3_2_3+VI3_2_3+VS4_2_3+VI4_2_3+VR_2_3+
    M_4_5+S1_4_5+I1_4_5+S2_4_5+I2_4_5+S3_4_5+I3_4_5+S4_4_5+I4_4_5+R_4_5+V_4_5+
    VS1_4_5+VI1_4_5+VS2_4_5+VI2_4_5+VS3_4_5+VI3_4_5+VS4_4_5+VI4_4_5+VR_4_5+
    M_6_11+S1_6_11+I1_6_11+S2_6_11+I2_6_11+S3_6_11+I3_6_11+S4_6_11+I4_6_11+R_6_11+V_6_11+
    VS1_6_11+VI1_6_11+VS2_6_11+VI2_6_11+VS3_6_11+VI3_6_11+VS4_6_11+VI4_6_11+VR_6_11+
    S1_12_23+I1_12_23+S2_12_23+I2_12_23+S3_12_23+I3_12_23+S4_12_23+I4_12_23+R_12_23+V_12_23+
    VS1_12_23+VI1_12_23+VS2_12_23+VI2_12_23+VS3_12_23+VI3_12_23+VS4_12_23+VI4_12_23+VR_12_23+
    S1_24_35+I1_24_35+S2_24_35+I2_24_35+S3_24_35+I3_24_35+S4_24_35+I4_24_35+R_24_35+V_24_35+
    VS1_24_35+VI1_24_35+VS2_24_35+VI2_24_35+VS3_24_35+VI3_24_35+VS4_24_35+VI4_24_35+VR_24_35+
    S1_36_47+I1_36_47+S2_36_47+I2_36_47+S3_36_47+I3_36_47+S4_36_47+I4_36_47+R_36_47+V_36_47+
    VS1_36_47+VI1_36_47+VS2_36_47+VI2_36_47+VS3_36_47+VI3_36_47+VS4_36_47+VI4_36_47+VR_36_47+
    S1_48_59+I1_48_59+S2_48_59+I2_48_59+S3_48_59+I3_48_59+S4_48_59+I4_48_59+R_48_59+V_48_59+
    VS1_48_59+VI1_48_59+VS2_48_59+VI2_48_59+VS3_48_59+VI3_48_59+VS4_48_59+VI4_48_59+VR_48_59+
    S1_60_up+I1_60_up+S2_60_up+I2_60_up+S3_60_up+I3_60_up+S4_60_up+I4_60_up+R_60_up+V_60_up+
    VS1_60_up+VI1_60_up+VS2_60_up+VI2_60_up+VS3_60_up+VI3_60_up+VS4_60_up+VI4_60_up+VR_60_up
  
  
  
  
  
  #indicator variables for coverage years
  vacc_era_2006=ifelse(trunc(t)>=2006*365 & trunc(t)<2007*365,1,0);
  vacc_era_2007=ifelse(trunc(t)>=2007*365 & trunc(t)<2008*365,1,0);
  vacc_era_2008=ifelse(trunc(t)>=2008*365 & trunc(t)<2009*365,1,0);
  vacc_era_2009=ifelse(trunc(t)>=2009*365 & trunc(t)<2010*365,1,0);
  vacc_era_2010=ifelse(trunc(t)>=2010*365 & trunc(t)<2011*365,1,0);
  vacc_era_2011=ifelse(trunc(t)>=2011*365 & trunc(t)<2012*365,1,0);
  vacc_era_2012=ifelse(trunc(t)>=2012*365 & trunc(t)<2013*365,1,0);
  vacc_era_2013=ifelse(trunc(t)>=2013*365 & trunc(t)<2014*365,1,0);
  vacc_era_2014=ifelse(trunc(t)>=2014*365 & trunc(t)<2015*365,1,0);vacc_era_2015=ifelse(trunc(t)>=2015*365 & trunc(t)<2016*365,1,0);vacc_era_2016=ifelse(trunc(t)>=2016*365 & trunc(t)<2017*365,1,0);vacc_era_2017=ifelse(trunc(t)>=2017*365 & trunc(t)<2018*365,1,0);vacc_era_2018=ifelse(trunc(t)>=2018*365 & trunc(t)<2019*365,1,0);vacc_era_2019=ifelse(trunc(t)>=2019*365 & trunc(t)<2020*365,1,0); vacc_era_2020=ifelse(trunc(t)>=2020*365 & trunc(t)<2021*365,1,0) ;vacc_era_2021=ifelse(trunc(t)>=2021*365 & trunc(t)<2022*365,1,0) ;vacc_era_2022=ifelse(trunc(t)>=2022*365 & trunc(t)<2023*365,1,0) ;vacc_era_2023=ifelse(trunc(t)>=2023*365 & trunc(t)<2024*365,1,0) ;vacc_era_2024=ifelse(trunc(t)>=2024*365 & trunc(t)<2025*365,1,0) ;vacc_era_2025=ifelse(trunc(t)>=2025*365 & trunc(t)<2026*365,1,0) ;vacc_era_2026=ifelse(trunc(t)>=2026*365 & trunc(t)<2027*365,1,0) ;vacc_era_2027=ifelse(trunc(t)>=2027*365 & trunc(t)<2028*365,1,0) ;vacc_era_2028=ifelse(trunc(t)>=2028*365 & trunc(t)<2029*365,1,0) ;vacc_era_2029=ifelse(trunc(t)>=2029*365 & trunc(t)<2030*365,1,0) ;vacc_era_2030=ifelse(trunc(t)>=2030*365 & trunc(t)<2031*365,1,0) ;vacc_era_2031=ifelse(trunc(t)>=2031*365 & trunc(t)<2032*365,1,0) ;vacc_era_2032=ifelse(trunc(t)>=2032*365 & trunc(t)<2033*365,1,0) ;vacc_era_2033=ifelse(trunc(t)>=2033*365 & trunc(t)<2034*365,1,0) ;vacc_era_2034=ifelse(trunc(t)>=2034*365 & trunc(t)<2035*365,1,0) ;vacc_era_2035=ifelse(trunc(t)>=2035*365 & trunc(t)<2036*365,1,0) ;vacc_era_2036=ifelse(trunc(t)>=2036*365 & trunc(t)<2037*365,1,0) 
  
  #ONLY 2 AND 4 MONTH OLDS GET VACCINATED
  vac1 = vac_imm*((p_vac1[1]*vacc_era_2006)+(p_vac1[2]*vacc_era_2007)+(p_vac1[3]*vacc_era_2008)+(p_vac1[4]*vacc_era_2009)+(p_vac1[5]*vacc_era_2010)+(p_vac1[6]*vacc_era_2011)+
                      (p_vac1[7]*vacc_era_2012)+(p_vac1[8]*vacc_era_2013)+(p_vac1[9]*vacc_era_2014)+(p_vac1[10]*vacc_era_2015)+(p_vac1[11]*vacc_era_2016)+
                      (p_vac1[12]*vacc_era_2017)+(p_vac1[13]*vacc_era_2018)+(p_vac1[14]*vacc_era_2019)+
                      (p_vac1[15]*vacc_era_2020)+(p_vac1[16]*vacc_era_2021)+(p_vac1[17]*vacc_era_2022)+(p_vac1[18]*vacc_era_2023)+(p_vac1[19]*vacc_era_2024)+(p_vac1[20]*vacc_era_2025)+
                      (p_vac1[21]*vacc_era_2026)+(p_vac1[22]*vacc_era_2027)+(p_vac1[23]*vacc_era_2028)+(p_vac1[24]*vacc_era_2029)+(p_vac1[24]*vacc_era_2030)+(p_vac1[25]*vacc_era_2031)+
                      (p_vac1[27]*vacc_era_2032)+(p_vac1[28]*vacc_era_2033)+(p_vac1[29]*vacc_era_2034)+(p_vac1[30]*vacc_era_2035)+(p_vac1[31]*vacc_era_2036))           #effectively immunised after full vaccine course of 1 doses

  vac2 = vac_imm*((p_vac2[1]*vacc_era_2006)+(p_vac2[2]*vacc_era_2007)+(p_vac2[3]*vacc_era_2008)+(p_vac2[4]*vacc_era_2009)+(p_vac2[5]*vacc_era_2010)+(p_vac2[6]*vacc_era_2011)+
                    (p_vac2[7]*vacc_era_2012)+(p_vac2[8]*vacc_era_2013)+(p_vac2[9]*vacc_era_2014)+(p_vac2[10]*vacc_era_2015)+(p_vac2[11]*vacc_era_2016)+
                    (p_vac2[12]*vacc_era_2017)+(p_vac2[13]*vacc_era_2018)+(p_vac2[14]*vacc_era_2019)+
                    (p_vac2[15]*vacc_era_2020)+(p_vac2[16]*vacc_era_2021)+(p_vac2[17]*vacc_era_2022)+(p_vac2[18]*vacc_era_2023)+(p_vac2[19]*vacc_era_2024)+(p_vac2[20]*vacc_era_2025)+
                    (p_vac2[21]*vacc_era_2026)+(p_vac2[22]*vacc_era_2027)+(p_vac2[23]*vacc_era_2028)+(p_vac2[24]*vacc_era_2029)+(p_vac2[24]*vacc_era_2030)+(p_vac2[25]*vacc_era_2031)+
                    (p_vac2[27]*vacc_era_2032)+(p_vac2[28]*vacc_era_2033)+(p_vac2[29]*vacc_era_2034)+(p_vac2[30]*vacc_era_2035)+(p_vac2[31]*vacc_era_2036))           #effectively immunised after full vaccine course of 1 doses
  
  #!!!NOTE: births currently set equal to deaths so population does not fluctuate
  births=brate*(N)
  
  #Individuals in each age group who are infectious, we assume that only symptomatic individuals (alpha) are infectious  AND that non-primary infections are less infections as per rel_inf 
  Infec = (alpha1*I1_0_1 + (alpha2*I2_0_1+alpha3*I3_0_1+alpha4*I4_0_1)*rel_inf)+
    (alpha1*I1_2_3 + (alpha2*I2_2_3+alpha3*I3_2_3+alpha4*I4_2_3)*rel_inf)+
    (alpha1*VI1_2_3 + (alpha2*VI2_2_3+alpha3*VI3_2_3+alpha4*VI4_2_3)*rel_inf)+
    (alpha1*I1_4_5 + (alpha2*I2_4_5+alpha3*I3_4_5+alpha4*I4_4_5)*rel_inf)+
    (alpha1*VI1_4_5 + (alpha2*VI2_4_5+alpha3*VI3_4_5+alpha4*VI4_4_5)*rel_inf)+
    (alpha1*I1_6_11 + (alpha2*I2_6_11+alpha3*I3_6_11+alpha4*I4_6_11)*rel_inf)+
    (alpha1*VI1_6_11 + (alpha2*VI2_6_11+alpha3*VI3_6_11+alpha4*VI4_6_11)*rel_inf)+
    (alpha1*I1_12_23 + (alpha2*I2_12_23+alpha3*I3_12_23+alpha4*I4_12_23)*rel_inf)+
    (alpha1*VI1_12_23 + (alpha2*VI2_12_23+alpha3*VI3_12_23+alpha4*VI4_12_23)*rel_inf)+
    (alpha1*I1_24_35 + (alpha2*I2_24_35+alpha3*I3_24_35+alpha4*I4_24_35)*rel_inf)+
    (alpha1*VI1_24_35 + (alpha2*VI2_24_35+alpha3*VI3_24_35+alpha4*VI4_24_35)*rel_inf)+
    (alpha1*I1_36_47 + (alpha2*I2_36_47+alpha3*I3_36_47+alpha4*I4_36_47)*rel_inf)+
    (alpha1*VI1_36_47 + (alpha2*VI2_36_47+alpha3*VI3_36_47+alpha4*VI4_36_47)*rel_inf)+
    (alpha1*I1_48_59 + (alpha2*I2_48_59+alpha3*I3_48_59+alpha4*I4_48_59)*rel_inf)+
    (alpha1*VI1_48_59 + (alpha2*VI2_48_59+alpha3*VI3_48_59+alpha4*VI4_48_59)*rel_inf)+
    (alpha1*I1_60_up + (alpha2*I2_60_up+alpha3*I3_60_up+alpha4*I4_60_up)*rel_inf)+
    (alpha1*VI1_60_up + (alpha2*VI2_60_up+alpha3*VI3_60_up+alpha4*VI4_60_up)*rel_inf)
  
#Force of Infection equations
FOI1=q1*(Infec/N)
FOI2=q2*(Infec/N)
FOI3=q3*(Infec/N)


#Save FOI numbers for year before vaccine introduction
#Need the double headed arrows so that this can be saved outside the function
if(trunc(t)==(yintro-1)){
  FOI1.pre<<-FOI1
  FOI2.pre<<-FOI2
  FOI3.pre<<-FOI3  
}

if(DE_FOI=='yes'){
  FOI1=ifelse(trunc(t)<yintro,q1*(Infec/N),FOI1.pre)
  FOI2=ifelse(trunc(t)<yintro,q2*(Infec/N),FOI2.pre)
  FOI3=ifelse(trunc(t)<yintro,q3*(Infec/N),FOI3.pre)
}


  #EQUATIONS
  
  #Infant age group:  0 -1 months old
  dM_0_1 = births - eta*M_0_1 - age_rate_0_5*M_0_1 - m_0_59*M_0_1		                     					
  dS1_0_1 = eta*M_0_1 - FOI1*S1_0_1 - age_rate_0_5*S1_0_1 + rec*R_0_1	 - m_0_59*S1_0_1                   					
  dI1_0_1 = FOI1*S1_0_1 - gamma*I1_0_1 - age_rate_0_5*I1_0_1 	  - m_0_59*I1_0_1                  				
  dS2_0_1 =  P1*gamma*I1_0_1 - FOI1*S2_0_1 - age_rate_0_5*S2_0_1	- m_0_59*S2_0_1         							
  dI2_0_1 = FOI1*S2_0_1 - gamma*I2_0_1 - age_rate_0_5*I2_0_1 	  - m_0_59*I2_0_1         					
  dS3_0_1 = gamma*P2*I2_0_1 - FOI1*S3_0_1 - age_rate_0_5*S3_0_1  - m_0_59*S3_0_1                					
  dI3_0_1 = FOI1*S3_0_1 - gamma*I3_0_1 - age_rate_0_5*I3_0_1    - m_0_59*I3_0_1                   				
  dS4_0_1 = gamma*P3*I3_0_1 - FOI1*S4_0_1 - age_rate_0_5*S4_0_1 	- m_0_59*S4_0_1	 				
  dI4_0_1 = FOI1*S4_0_1 - gamma*I4_0_1 	- age_rate_0_5*I4_0_1 	   - m_0_59*I4_0_1   				
  dR_0_1 = gamma*(1-P1)*I1_0_1 + gamma*(1-P2)*I2_0_1 + gamma*(1-P3)*I3_0_1 + gamma*I4_0_1 - age_rate_0_5*R_0_1 - rec*R_0_1    - m_0_59*R_0_1        	#immune to infection
  
  #Infant age group:  2 - 3 months old
  dM_2_3 =  age_rate_0_5*M_0_1*(1-vac1) - M_2_3*(age_rate_0_5 + eta + m_0_59)    
  dS1_2_3 = eta*M_2_3 + age_rate_0_5*S1_0_1*(1-vac1) + rec*(R_2_3) - S1_2_3*(age_rate_0_5  + m_0_59 + FOI1)
  dI1_2_3 = FOI1*S1_2_3 + age_rate_0_5*I1_0_1*(1-vac1) - I1_2_3*(gamma  + age_rate_0_5 + m_0_59)
  dS2_2_3 = P1*gamma*I1_2_3 + age_rate_0_5*S2_0_1*(1-vac1) - S2_2_3*(FOI1 + age_rate_0_5 + m_0_59)
  dI2_2_3 = FOI1*S2_2_3 + age_rate_0_5*I2_0_1*(1-vac1) - I2_2_3*(gamma  + age_rate_0_5 + m_0_59)
  dS3_2_3 = gamma*P2*I2_2_3 + age_rate_0_5*S3_0_1*(1-vac1) - S3_2_3*(FOI1 + age_rate_0_5 + m_0_59) 
  dI3_2_3 = FOI1*S3_2_3 + age_rate_0_5*I3_0_1*(1-vac1) - I3_2_3*(gamma  + age_rate_0_5 + m_0_59)           
  dS4_2_3 = gamma*P3*I3_2_3 + age_rate_0_5*S4_0_1*(1-vac1) - S4_2_3*(FOI1 + age_rate_0_5 + m_0_59)
  dI4_2_3 = FOI1*S4_2_3 + age_rate_0_5*I4_0_1*(1-vac1) - I4_2_3*(gamma  + age_rate_0_5 + m_0_59)  
  dR_2_3 = gamma*((1-P1)*I1_2_3 + (1-P2)*I2_2_3 + (1-P3)*I3_2_3 + I4_2_3) + age_rate_0_5*R_0_1*(1-vac1)  - R_2_3*(rec + age_rate_0_5 + m_0_59)
  
  #Successful vaccination
  dV_2_3 = (age_rate_0_5*vac1)*(M_0_1*(1-P1) + S1_0_1*(1-P1) + I1_0_1*(1-P1) +
                                  S2_0_1*(1-P2) + I2_0_1*(1-P2) +
                                  S3_0_1*(1-P3) + I3_0_1*(1-P3)) + age_rate_0_5*vac1*(S4_0_1 + I4_0_1 + R_0_1) - V_2_3*(rec + age_rate_0_5 + m_0_59)
  
  #Unsuccessful vaccination
  dVS1_2_3 =  rec*(VR_2_3 + V_2_3) - VS1_2_3*(age_rate_0_5  + m_0_59 + FOI1)
  dVI1_2_3 = FOI1*VS1_2_3 - VI1_2_3*(gamma  + age_rate_0_5 + m_0_59)
  dVS2_2_3 = P1*gamma*VI1_2_3 + age_rate_0_5*M_0_1*vac1*P1 + age_rate_0_5*S1_0_1*vac1*P1  - VS2_2_3*(FOI1 + age_rate_0_5 + m_0_59)
  dVI2_2_3 = age_rate_0_5*I1_0_1*vac1*P1 + FOI1*VS2_2_3 - VI2_2_3*(gamma  + age_rate_0_5 + m_0_59)
  dVS3_2_3 = gamma*P2*VI2_2_3 + age_rate_0_5*S2_0_1*vac1*P2 - VS3_2_3*(FOI1 + age_rate_0_5 + m_0_59) 
  dVI3_2_3 = age_rate_0_5*I2_0_1*vac1*P2 + FOI1*VS3_2_3 - VI3_2_3*(gamma  + age_rate_0_5 + m_0_59)           
  dVS4_2_3 = gamma*P3*VI3_2_3 + age_rate_0_5*S3_0_1*vac1*P3 - VS4_2_3*(FOI1 + age_rate_0_5 + m_0_59)
  dVI4_2_3 = age_rate_0_5*I3_0_1*vac1*P3 + FOI1*VS4_2_3 - VI4_2_3*(gamma  + age_rate_0_5 + m_0_59)  
  dVR_2_3 = gamma*((1-P1)*VI1_2_3 + (1-P2)*VI2_2_3 + (1-P3)*VI3_2_3 + VI4_2_3) - VR_2_3*(rec + age_rate_0_5 + m_0_59)
  
  
  #Infant age group:  4- 5 months old
  dM_4_5 =  age_rate_0_5*M_2_3*(1-vac2) - M_4_5*(age_rate_0_5 + eta + m_0_59)    
  dS1_4_5 = eta*M_4_5 + age_rate_0_5*S1_2_3*(1-vac2) + rec*(R_4_5) - S1_4_5*(age_rate_0_5  + m_0_59 + FOI2)
  dI1_4_5 = FOI2*S1_4_5 + age_rate_0_5*I1_2_3*(1-vac2) - I1_4_5*(gamma  + age_rate_0_5 + m_0_59)
  dS2_4_5 = P1*gamma*I1_4_5 + age_rate_0_5*S2_2_3*(1-vac2) - S2_4_5*(FOI2 + age_rate_0_5+ m_0_59)
  dI2_4_5 = FOI2*S2_4_5 + age_rate_0_5*I2_2_3*(1-vac2) - I2_4_5*(gamma  + age_rate_0_5 + m_0_59)
  dS3_4_5 = gamma*P2*I2_4_5 + age_rate_0_5*S3_2_3*(1-vac2) - S3_4_5*(FOI2 + age_rate_0_5 + m_0_59) 
  dI3_4_5 = FOI2*S3_4_5 + age_rate_0_5*I3_2_3*(1-vac2) - I3_4_5*(gamma  + age_rate_0_5 + m_0_59)           
  dS4_4_5 = gamma*P3*I3_4_5 + age_rate_0_5*S4_2_3*(1-vac2) - S4_4_5*(FOI2 + age_rate_0_5 + m_0_59)
  dI4_4_5 = FOI2*S4_4_5 + age_rate_0_5*I4_2_3*(1-vac2) - I4_4_5*(gamma  + age_rate_0_5 + m_0_59)  
  dR_4_5 = gamma*((1-P1)*I1_4_5 + (1-P2)*I2_4_5 + (1-P3)*I3_4_5 + I4_4_5) + age_rate_0_5*R_2_3*(1-vac2) - R_4_5*(rec + age_rate_0_5 + m_0_59)
  
  #Successfully vaccinated
  dV_4_5 = age_rate_0_5*V_2_3 + 
    (age_rate_0_5*vac2)*(M_2_3*(1-P1) + S1_2_3*(1-P1) + VS1_2_3*(1-P1) + VI1_2_3*(1-P1) + I1_2_3*(1-P1) +
                           S2_2_3*(1-P2) + VS2_2_3*(1-P2) + I2_2_3*(1-P2) + VI2_2_3*(1-P2) +
                           S3_2_3*(1-P3) + VS3_2_3*(1-P3) + I3_2_3*(1-P3) + VI3_2_3*(1-P3)) + 
    age_rate_0_5*vac2*(S4_2_3 + I4_2_3 + R_2_3 + VS4_2_3 + VI4_2_3 + VR_2_3)  - V_4_5*(rec + age_rate_0_5 + m_0_59)
  
  #Unsuccessful vaccination
  dVS1_4_5 = age_rate_0_5*(VS1_2_3)*(1-vac2) + rec*(VR_4_5 + V_4_5) - VS1_4_5*(age_rate_0_5  + m_0_59 + FOI1)
  dVI1_4_5 = age_rate_0_5*(VI1_2_3)*(1-vac2) + FOI1*VS1_4_5 - VI1_4_5*(gamma  + age_rate_0_5 + m_0_59)
  dVS2_4_5 = age_rate_0_5*(VS2_2_3)*(1-vac2) + P1*gamma*VI1_4_5 + vac2*P1*(age_rate_0_5*M_2_3 + age_rate_0_5*S1_2_3  + age_rate_0_5*VS1_2_3) - VS2_4_5*(FOI1 + age_rate_0_5 + m_0_59)
  dVI2_4_5 = age_rate_0_5*(VI2_2_3)*(1-vac2) + age_rate_0_5*(I1_2_3+VI1_2_3)*vac2*P1 + FOI1*VS2_4_5 - VI2_4_5*(gamma  + age_rate_0_5 + m_0_59)
  dVS3_4_5 = age_rate_0_5*(VS3_2_3)*(1-vac2) + gamma*P2*VI2_4_5 + age_rate_0_5*(S2_2_3+VS2_2_3)*vac2*P2 - VS3_4_5*(FOI1 + age_rate_0_5 + m_0_59) 
  dVI3_4_5 = age_rate_0_5*(VI3_2_3)*(1-vac2) + age_rate_0_5*(I2_2_3+VI2_2_3)*vac2*P2 + FOI1*VS3_4_5 - VI3_4_5*(gamma  + age_rate_0_5 + m_0_59)           
  dVS4_4_5 = age_rate_0_5*(VS4_2_3)*(1-vac2) + gamma*P3*VI3_4_5 + age_rate_0_5*(S3_2_3+VS3_2_3)*vac2*P3 - VS4_4_5*(FOI1 + age_rate_0_5 + m_0_59)
  dVI4_4_5 = age_rate_0_5*(VI4_2_3)*(1-vac2) + age_rate_0_5*(I3_2_3+VI3_2_3)*vac2*P3 + FOI1*VS4_4_5 - VI4_4_5*(gamma  + age_rate_0_5 + m_0_59)  
  dVR_4_5 = age_rate_0_5*(VR_2_3)*(1-vac2) + gamma*((1-P1)*VI1_4_5 + (1-P2)*VI2_4_5 + (1-P3)*VI3_4_5 + VI4_4_5) - VR_4_5*(rec + age_rate_0_5 + m_0_59)
  
  
  #Infant age group:  6- 7 months old
  dM_6_11 =  age_rate_0_5*M_4_5 - M_6_11*(age_rate_6_11 + eta + m_0_59)    
  dS1_6_11 = eta*M_6_11 + age_rate_0_5*S1_4_5 + rec*(R_6_11) - S1_6_11*(age_rate_6_11  + m_0_59 + FOI2)
  dI1_6_11 = FOI2*S1_6_11 + age_rate_0_5*I1_4_5 - I1_6_11*(gamma  + age_rate_6_11 + m_0_59)
  dS2_6_11 = P1*gamma*I1_6_11 + age_rate_0_5*S2_4_5 - S2_6_11*(FOI2 + age_rate_6_11+ m_0_59)
  dI2_6_11 = FOI2*S2_6_11 + age_rate_0_5*I2_4_5 - I2_6_11*(gamma  + age_rate_6_11 + m_0_59)
  dS3_6_11 = gamma*P2*I2_6_11 + age_rate_0_5*S3_4_5 - S3_6_11*(FOI2 + age_rate_6_11 + m_0_59) 
  dI3_6_11 = FOI2*S3_6_11 + age_rate_0_5*I3_4_5 - I3_6_11*(gamma + age_rate_6_11 + m_0_59)           
  dS4_6_11 = gamma*P3*I3_6_11 + age_rate_0_5*S4_4_5 - S4_6_11*(FOI2 + age_rate_6_11 + m_0_59)
  dI4_6_11 = FOI2*S4_6_11 + age_rate_0_5*I4_4_5 - I4_6_11*(gamma  + age_rate_6_11 + m_0_59)  
  dR_6_11 = gamma*((1-P1)*I1_6_11 + (1-P2)*I2_6_11 + (1-P3)*I3_6_11 + I4_6_11) + age_rate_0_5*R_4_5  - R_6_11*(rec + age_rate_6_11 + m_0_59)
  
  #Successful vaccination 
  dV_6_11 = age_rate_0_5*V_4_5 - V_6_11*(rec + age_rate_6_11 + m_0_59)
  
  #Unsuccessful vaccination
  dVS1_6_11 = age_rate_0_5*(VS1_4_5) + rec*(VR_6_11 + V_6_11) - VS1_6_11*(age_rate_6_11  + m_0_59 + FOI1)
  dVI1_6_11 = age_rate_0_5*(VI1_4_5) + FOI1*VS1_6_11 - VI1_6_11*(gamma  + age_rate_6_11 + m_0_59)
  dVS2_6_11 = age_rate_0_5*(VS2_4_5) + P1*gamma*VI1_6_11 - VS2_6_11*(FOI1 + age_rate_6_11 + m_0_59)
  dVI2_6_11 = age_rate_0_5*(VI2_4_5) + FOI1*VS2_6_11 - VI2_6_11*(gamma  + age_rate_6_11 + m_0_59)
  dVS3_6_11 = age_rate_0_5*(VS3_4_5) + gamma*P2*VI2_6_11 - VS3_6_11*(FOI1 + age_rate_6_11 + m_0_59) 
  dVI3_6_11 = age_rate_0_5*(VI3_4_5) + FOI1*VS3_6_11 - VI3_6_11*(gamma  + age_rate_6_11 + m_0_59)           
  dVS4_6_11 = age_rate_0_5*(VS4_4_5) + gamma*P3*VI3_6_11 - VS4_6_11*(FOI1 + age_rate_6_11 + m_0_59)
  dVI4_6_11 = age_rate_0_5*(VI4_4_5) + FOI1*VS4_6_11 - VI4_6_11*(gamma  + age_rate_6_11 + m_0_59)  
  dVR_6_11 = age_rate_0_5*(VR_4_5) + gamma*((1-P1)*VI1_6_11 + (1-P2)*VI2_6_11 + (1-P3)*VI3_6_11 + VI4_6_11) - VR_6_11*(rec + age_rate_6_11 + m_0_59)
  
  
  #Infant age groups: 12 - 23 months old
  dS1_12_23 = age_rate_6_11*(M_6_11+S1_6_11) + rec*(R_12_23) - S1_12_23*(FOI2 + age_rate_12_60 + m_0_59)                      
  dI1_12_23 = FOI2*S1_12_23 + age_rate_6_11*I1_6_11 - I1_12_23*(gamma + age_rate_12_60 + m_0_59)
  dS2_12_23 = P1*gamma*I1_12_23 + age_rate_6_11*S2_6_11 - S2_12_23*(FOI2 + age_rate_12_60 + m_0_59)
  dI2_12_23 = FOI2*S2_12_23 + age_rate_6_11*I2_6_11 - I2_12_23*(gamma + age_rate_12_60 + m_0_59)
  dS3_12_23 = gamma*P2*I2_12_23 + age_rate_6_11*S3_6_11 - S3_12_23*(FOI2 + age_rate_12_60 + m_0_59)
  dI3_12_23 = FOI2*S3_12_23 + age_rate_6_11*I3_6_11 - I3_12_23*(gamma + age_rate_12_60 + m_0_59)          
  dS4_12_23 = gamma*P3*I3_12_23 + age_rate_6_11*S4_6_11 - S4_12_23*(FOI2 + age_rate_12_60 + m_0_59) 
  dI4_12_23 = FOI2*S4_12_23 + age_rate_6_11*I4_6_11 - I4_12_23*(gamma + age_rate_12_60 + m_0_59)          
  dR_12_23 = gamma*(1-P1)*I1_12_23 + (gamma*(1-P2)*I2_12_23) + (gamma*(1-P3)*I3_12_23) + (gamma*I4_12_23) + age_rate_6_11*R_6_11 - R_12_23*(age_rate_12_60 + rec + m_0_59)       
  
  #Successful vaccination
  dV_12_23 = age_rate_6_11*V_6_11 - V_12_23*(rec + age_rate_12_60 + m_0_59)
  
  #Unsuccessful vaccination
  dVS1_12_23 = age_rate_6_11*(VS1_6_11) + rec*(VR_12_23 + V_12_23) - VS1_12_23*(age_rate_12_60  + m_0_59 + FOI1)
  dVI1_12_23 = age_rate_6_11*(VI1_6_11) + FOI1*VS1_12_23 - VI1_12_23*(gamma  + age_rate_12_60 + m_0_59)
  dVS2_12_23 = age_rate_6_11*(VS2_6_11) + P1*gamma*VI1_12_23 - VS2_12_23*(FOI1 + age_rate_12_60 + m_0_59)
  dVI2_12_23 = age_rate_6_11*(VI2_6_11) + FOI1*VS2_12_23 - VI2_12_23*(gamma  + age_rate_12_60 + m_0_59)
  dVS3_12_23 = age_rate_6_11*(VS3_6_11) + gamma*P2*VI2_12_23 - VS3_12_23*(FOI1 + age_rate_12_60 + m_0_59) 
  dVI3_12_23 = age_rate_6_11*(VI3_6_11) + FOI1*VS3_12_23 - VI3_12_23*(gamma  + age_rate_12_60 + m_0_59)           
  dVS4_12_23 = age_rate_6_11*(VS4_6_11) + gamma*P3*VI3_12_23 - VS4_12_23*(FOI1 + age_rate_12_60 + m_0_59)
  dVI4_12_23 = age_rate_6_11*(VI4_6_11) + FOI1*VS4_12_23 - VI4_12_23*(gamma  + age_rate_12_60 + m_0_59)  
  dVR_12_23 = age_rate_6_11*(VR_6_11) + gamma*((1-P1)*VI1_12_23 + (1-P2)*VI2_12_23 + (1-P3)*VI3_12_23 + VI4_12_23) - VR_12_23*(rec + age_rate_12_60 + m_0_59)
  
  #Infant age groups: 24 - 35 months old
  dS1_24_35 = age_rate_12_60*S1_12_23 + rec*(R_24_35) - S1_24_35*(FOI3 + age_rate_12_60 + m_0_59) 
  dI1_24_35 = FOI3*S1_24_35 + age_rate_12_60*I1_12_23 - I1_24_35*(gamma + age_rate_12_60 + m_0_59)
  dS2_24_35 = P1*gamma*I1_24_35 + age_rate_12_60*S2_12_23 - S2_24_35*(FOI3 + age_rate_12_60 + m_0_59)
  dI2_24_35 = FOI3*S2_24_35 + age_rate_12_60*I2_12_23 - I2_24_35*(gamma + age_rate_12_60 + m_0_59)
  dS3_24_35 = gamma*P2*I2_24_35 + age_rate_12_60*S3_12_23 - S3_24_35*(FOI3 + age_rate_12_60 + m_0_59)
  dI3_24_35 = FOI3*S3_24_35 + age_rate_12_60*I3_12_23 - I3_24_35*(gamma + age_rate_12_60 + m_0_59)  
  dS4_24_35 = gamma*P3*I3_24_35 + age_rate_12_60*S4_12_23 - S4_24_35*(FOI3 + age_rate_12_60 + m_0_59) 
  dI4_24_35 = FOI3*S4_24_35 + age_rate_12_60*I4_12_23 - I4_24_35*(gamma + age_rate_12_60 + m_0_59)          
  dR_24_35 = gamma*(1-P1)*I1_24_35 + (gamma*(1-P2)*I2_24_35) + (gamma*(1-P3)*I3_24_35) + (gamma*I4_24_35) + age_rate_12_60*R_12_23 - R_24_35*(age_rate_12_60 + rec + m_0_59)       
  
  #Successful vaccination
  dV_24_35 = age_rate_12_60*V_12_23 - V_24_35*(rec + age_rate_12_60 + m_0_59)
  
  #Unsuccessful vaccination
  dVS1_24_35 = age_rate_12_60*(VS1_12_23) + rec*(VR_24_35 + V_24_35) - VS1_24_35*(age_rate_12_60  + m_0_59 + FOI1)
  dVI1_24_35 = age_rate_12_60*(VI1_12_23) + FOI1*VS1_24_35 - VI1_24_35*(gamma  + age_rate_12_60 + m_0_59)
  dVS2_24_35 = age_rate_12_60*(VS2_12_23) + P1*gamma*VI1_24_35 - VS2_24_35*(FOI1 + age_rate_12_60 + m_0_59)
  dVI2_24_35 = age_rate_12_60*(VI2_12_23) + FOI1*VS2_24_35 - VI2_24_35*(gamma  + age_rate_12_60 + m_0_59)
  dVS3_24_35 = age_rate_12_60*(VS3_12_23) + gamma*P2*VI2_24_35 - VS3_24_35*(FOI1 + age_rate_12_60 + m_0_59) 
  dVI3_24_35 = age_rate_12_60*(VI3_12_23) + FOI1*VS3_24_35 - VI3_24_35*(gamma  + age_rate_12_60 + m_0_59)           
  dVS4_24_35 = age_rate_12_60*(VS4_12_23) + gamma*P3*VI3_24_35 - VS4_24_35*(FOI1 + age_rate_12_60 + m_0_59)
  dVI4_24_35 = age_rate_12_60*(VI4_12_23) + FOI1*VS4_24_35 - VI4_24_35*(gamma  + age_rate_12_60 + m_0_59)  
  dVR_24_35 = age_rate_12_60*(VR_12_23) + gamma*((1-P1)*VI1_24_35 + (1-P2)*VI2_24_35 + (1-P3)*VI3_24_35 + VI4_24_35) - VR_24_35*(rec + age_rate_12_60 + m_0_59)
  
  #Infant age groups: 36 - 47 months old
  dS1_36_47 = age_rate_12_60*S1_24_35 + rec*(R_36_47) - S1_36_47*(FOI3 + age_rate_12_60 + m_0_59) 
  dI1_36_47 = FOI3*S1_36_47 + age_rate_12_60*I1_24_35 - I1_36_47*(gamma + age_rate_12_60 + m_0_59)
  dS2_36_47 = P1*gamma*I1_36_47 + age_rate_12_60*S2_24_35 - S2_36_47*(FOI3 + age_rate_12_60 + m_0_59)
  dI2_36_47 = FOI3*S2_36_47 + age_rate_12_60*I2_24_35 - I2_36_47*(gamma + age_rate_12_60 + m_0_59)
  dS3_36_47 = gamma*P2*I2_36_47 + age_rate_12_60*S3_24_35 - S3_36_47*(FOI3 + age_rate_12_60 + m_0_59)
  dI3_36_47 = FOI3*S3_36_47 + age_rate_12_60*I3_24_35 - I3_36_47*(gamma + age_rate_12_60 + m_0_59)  
  dS4_36_47 = gamma*P3*I3_36_47 + age_rate_12_60*S4_24_35 - S4_36_47*(FOI3 + age_rate_12_60 + m_0_59) 
  dI4_36_47 = FOI3*S4_36_47 + age_rate_12_60*I4_24_35 - I4_36_47*(gamma + age_rate_12_60 + m_0_59)          
  dR_36_47 = gamma*(1-P1)*I1_36_47 + (gamma*(1-P2)*I2_36_47) + (gamma*(1-P3)*I3_36_47) + (gamma*I4_36_47) + age_rate_12_60*R_24_35 - R_36_47*(age_rate_12_60 + rec + m_0_59)       
  
  #Succesful vaccination
  dV_36_47 = age_rate_12_60*V_24_35 - V_36_47*(rec + age_rate_12_60 + m_0_59)
  
  #Unsuccessful vaccination
  dVS1_36_47 = age_rate_12_60*(VS1_24_35) + rec*(VR_36_47 + V_36_47) - VS1_36_47*(age_rate_12_60  + m_0_59 + FOI1)
  dVI1_36_47 = age_rate_12_60*(VI1_24_35) + FOI1*VS1_36_47 - VI1_36_47*(gamma  + age_rate_12_60 + m_0_59)
  dVS2_36_47 = age_rate_12_60*(VS2_24_35) + P1*gamma*VI1_36_47 - VS2_36_47*(FOI1 + age_rate_12_60 + m_0_59)
  dVI2_36_47 = age_rate_12_60*(VI2_24_35) + FOI1*VS2_36_47 - VI2_36_47*(gamma  + age_rate_12_60 + m_0_59)
  dVS3_36_47 = age_rate_12_60*(VS3_24_35) + gamma*P2*VI2_36_47 - VS3_36_47*(FOI1 + age_rate_12_60 + m_0_59) 
  dVI3_36_47 = age_rate_12_60*(VI3_24_35) + FOI1*VS3_36_47 - VI3_36_47*(gamma  + age_rate_12_60 + m_0_59)           
  dVS4_36_47 = age_rate_12_60*(VS4_24_35) + gamma*P3*VI3_36_47 - VS4_36_47*(FOI1 + age_rate_12_60 + m_0_59)
  dVI4_36_47 = age_rate_12_60*(VI4_24_35) + FOI1*VS4_36_47 - VI4_36_47*(gamma  + age_rate_12_60 + m_0_59)  
  dVR_36_47 = age_rate_12_60*(VR_24_35) + gamma*((1-P1)*VI1_36_47 + (1-P2)*VI2_36_47 + (1-P3)*VI3_36_47 + VI4_36_47) - VR_36_47*(rec + age_rate_12_60 + m_0_59)
  
  
  #Infant age groups: 48 - 59 months old
  dS1_48_59 = age_rate_12_60*S1_36_47 + rec*(R_48_59) - S1_48_59*(FOI3 + age_rate_12_60 + m_0_59) 
  dI1_48_59 = FOI3*S1_48_59 + age_rate_12_60*I1_36_47 - I1_48_59*(gamma + age_rate_12_60 + m_0_59)
  dS2_48_59 = P1*gamma*I1_48_59 + age_rate_12_60*S2_36_47 - S2_48_59*(FOI3 + age_rate_12_60 + m_0_59)
  dI2_48_59 = FOI3*S2_48_59 + age_rate_12_60*I2_36_47 - I2_48_59*(gamma + age_rate_12_60 + m_0_59)
  dS3_48_59 = gamma*P2*I2_48_59 + age_rate_12_60*S3_36_47 - S3_48_59*(FOI3 + age_rate_12_60 + m_0_59)
  dI3_48_59 = FOI3*S3_48_59 + age_rate_12_60*I3_36_47 - I3_48_59*(gamma + age_rate_12_60 + m_0_59)  
  dS4_48_59 = gamma*P3*I3_48_59 + age_rate_12_60*S4_36_47 - S4_48_59*(FOI3 + age_rate_12_60 + m_0_59) 
  dI4_48_59 = FOI3*S4_48_59 + age_rate_12_60*I4_36_47 - I4_48_59*(gamma + age_rate_12_60 + m_0_59)          
  dR_48_59 = gamma*(1-P1)*I1_48_59 + (gamma*(1-P2)*I2_48_59) + (gamma*(1-P3)*I3_48_59) + (gamma*I4_48_59) + age_rate_12_60*R_36_47 - R_48_59*(age_rate_12_60 + rec + m_0_59)       
  
  #Succesful vaccination
  dV_48_59 = age_rate_12_60*V_36_47 -  V_48_59*(rec + age_rate_12_60 + m_0_59)
  
  
  #Unsuccessful vaccination
  dVS1_48_59 = age_rate_12_60*(VS1_36_47) + rec*(VR_48_59 + V_48_59) - VS1_48_59*(age_rate_12_60  + m_0_59 + FOI1)
  dVI1_48_59 = age_rate_12_60*(VI1_36_47) + FOI1*VS1_48_59 - VI1_48_59*(gamma  + age_rate_12_60 + m_0_59)
  dVS2_48_59 = age_rate_12_60*(VS2_36_47) + P1*gamma*VI1_48_59 - VS2_48_59*(FOI1 + age_rate_12_60 + m_0_59)
  dVI2_48_59 = age_rate_12_60*(VI2_36_47) + FOI1*VS2_48_59 - VI2_48_59*(gamma  + age_rate_12_60 + m_0_59)
  dVS3_48_59 = age_rate_12_60*(VS3_36_47) + gamma*P2*VI2_48_59 - VS3_48_59*(FOI1 + age_rate_12_60 + m_0_59) 
  dVI3_48_59 = age_rate_12_60*(VI3_36_47) + FOI1*VS3_48_59 - VI3_48_59*(gamma  + age_rate_12_60 + m_0_59)           
  dVS4_48_59 = age_rate_12_60*(VS4_36_47) + gamma*P3*VI3_48_59 - VS4_48_59*(FOI1 + age_rate_12_60 + m_0_59)
  dVI4_48_59 = age_rate_12_60*(VI4_36_47) + FOI1*VS4_48_59 - VI4_48_59*(gamma  + age_rate_12_60 + m_0_59)  
  dVR_48_59 = age_rate_12_60*(VR_36_47) + gamma*((1-P1)*VI1_48_59 + (1-P2)*VI2_48_59 + (1-P3)*VI3_48_59 + VI4_48_59) - VR_48_59*(rec + age_rate_12_60 + m_0_59)
  
  #Age group: child = 5-29 years old
  dS1_60_up = age_rate_12_60*S1_48_59 - FOI3*S1_60_up + rec*(R_60_up)          - m_60_up*S1_60_up	  	
  dI1_60_up = age_rate_12_60*I1_48_59 + FOI3*S1_60_up - gamma*I1_60_up       - m_60_up*I1_60_up                 	
  dS2_60_up =  P1*gamma*I1_60_up + age_rate_12_60*S2_48_59 - FOI3*S2_60_up   - m_60_up*S2_60_up         	
  dI2_60_up = FOI3*S2_60_up + age_rate_12_60*I2_48_59 - gamma*I2_60_up      - m_60_up*I2_60_up                	
  dS3_60_up = gamma*P2*I2_60_up + age_rate_12_60*S3_48_59 - FOI3*S3_60_up   - m_60_up*S3_60_up         	
  dI3_60_up = FOI3*S3_60_up + age_rate_12_60*I3_48_59 - gamma*I3_60_up     - m_60_up*I3_60_up               	
  dS4_60_up = gamma*P3*I3_60_up + age_rate_12_60*S4_48_59  - FOI3*S4_60_up  - m_60_up*S4_60_up
  dI4_60_up = FOI3*S4_60_up + age_rate_12_60*I4_48_59 - gamma*I4_60_up      - m_60_up*I4_60_up               	
  dR_60_up = gamma*(1-P1)*I1_60_up + gamma*(1-P2)*I2_60_up + gamma*(1-P3)*I3_60_up + gamma*I4_60_up + age_rate_12_60*R_48_59 - rec*R_60_up   - m_60_up*R_60_up     
  
  #Successful vaccination
  dV_60_up = age_rate_12_60*V_48_59 - V_60_up*(rec+ m_60_up)
  
  #Unsuccessful vaccination
  dVS1_60_up = age_rate_12_60*(VS1_48_59) + rec*(VR_60_up + V_60_up) - VS1_60_up*(m_60_up + FOI1)
  dVI1_60_up = age_rate_12_60*(VI1_48_59) + FOI1*VS1_60_up - VI1_60_up*(gamma  + m_60_up)
  dVS2_60_up = age_rate_12_60*(VS2_48_59) + P1*gamma*VI1_60_up - VS2_60_up*(FOI1 + m_60_up)
  dVI2_60_up = age_rate_12_60*(VI2_48_59) + FOI1*VS2_60_up - VI2_60_up*(gamma  + m_60_up)
  dVS3_60_up = age_rate_12_60*(VS3_48_59) + gamma*P2*VI2_60_up - VS3_60_up*(FOI1 + m_60_up) 
  dVI3_60_up = age_rate_12_60*(VI3_48_59) + FOI1*VS3_60_up - VI3_60_up*(gamma  + m_60_up)           
  dVS4_60_up = age_rate_12_60*(VS4_48_59) + gamma*P3*VI3_60_up - VS4_60_up*(FOI1 + m_60_up)
  dVI4_60_up = age_rate_12_60*(VI4_48_59) + FOI1*VS4_60_up - VI4_60_up*(gamma  + m_60_up)  
  dVR_60_up = age_rate_12_60*(VR_48_59) + gamma*((1-P1)*VI1_60_up + (1-P2)*VI2_60_up + (1-P3)*VI3_60_up + VI4_60_up) - VR_60_up*(rec + m_60_up)
  
  res = c(dM_0_1,dS1_0_1, dI1_0_1,dS2_0_1,dI2_0_1,dS3_0_1,dI3_0_1,dS4_0_1,dI4_0_1,dR_0_1,
          dM_2_3,dS1_2_3, dI1_2_3,dS2_2_3,dI2_2_3,dS3_2_3,dI3_2_3,dS4_2_3,dI4_2_3,dR_2_3,dV_2_3,
          dVS1_2_3, dVI1_2_3,dVS2_2_3,dVI2_2_3,dVS3_2_3,dVI3_2_3,dVS4_2_3,dVI4_2_3,dVR_2_3,
          dM_4_5,dS1_4_5,dI1_4_5,dS2_4_5,dI2_4_5,dS3_4_5,dI3_4_5,dS4_4_5,dI4_4_5,dR_4_5,dV_4_5,
          dVS1_4_5, dVI1_4_5,dVS2_4_5,dVI2_4_5,dVS3_4_5,dVI3_4_5,dVS4_4_5,dVI4_4_5,dVR_4_5,
          dM_6_11,dS1_6_11,dI1_6_11,dS2_6_11,dI2_6_11,dS3_6_11,dI3_6_11,dS4_6_11,dI4_6_11,dR_6_11,dV_6_11,
          dVS1_6_11, dVI1_6_11,dVS2_6_11,dVI2_6_11,dVS3_6_11,dVI3_6_11,dVS4_6_11,dVI4_6_11,dVR_6_11,
          dS1_12_23,dI1_12_23,dS2_12_23,dI2_12_23,dS3_12_23,dI3_12_23,dS4_12_23,dI4_12_23,dR_12_23,dV_12_23,
          dVS1_12_23, dVI1_12_23,dVS2_12_23,dVI2_12_23,dVS3_12_23,dVI3_12_23,dVS4_12_23,dVI4_12_23,dVR_12_23,
          dS1_24_35,dI1_24_35,dS2_24_35,dI2_24_35,dS3_24_35,dI3_24_35,dS4_24_35,dI4_24_35,dR_24_35,dV_24_35,
          dVS1_24_35, dVI1_24_35,dVS2_24_35,dVI2_24_35,dVS3_24_35,dVI3_24_35,dVS4_24_35,dVI4_24_35,dVR_24_35,
          dS1_36_47,dI1_36_47,dS2_36_47,dI2_36_47,dS3_36_47,dI3_36_47,dS4_36_47,dI4_36_47,dR_36_47,dV_36_47,
          dVS1_36_47, dVI1_36_47,dVS2_36_47,dVI2_36_47,dVS3_36_47,dVI3_36_47,dVS4_36_47,dVI4_36_47,dVR_36_47,
          dS1_48_59,dI1_48_59,dS2_48_59,dI2_48_59,dS3_48_59,dI3_48_59,dS4_48_59,dI4_48_59,dR_48_59,dV_48_59,
          dVS1_48_59, dVI1_48_59,dVS2_48_59,dVI2_48_59,dVS3_48_59,dVI3_48_59,dVS4_48_59,dVI4_48_59,dVR_48_59,
          dS1_60_up,dI1_60_up,dS2_60_up,dI2_60_up, dS3_60_up,dI3_60_up,dS4_60_up, dI4_60_up,dR_60_up,dV_60_up,
          dVS1_60_up, dVI1_60_up,dVS2_60_up,dVI2_60_up,dVS3_60_up,dVI3_60_up,dVS4_60_up,dVI4_60_up,dVR_60_up)
  list(res)
}

#Read in dataset that contains all model parameter values and country specific demographics
#country_params<-get(load(paste0(pathi,input_data)))
#country_params<-read.csv(paste0(pathi,input_data))
#coverage_params<-read.csv(paste0(pathi,coverage_data))

#country_params2 <- merge(country_params,coverage_params,by=c("country","country_code"))
#param_names=names(country_params2)


#country_params.list <- setNames(split(country_params2, seq(nrow(country_params2))), rownames(country_params2))


country_runs.fn=function(data_init){
  #for(i in 1:length(set_run_ids)){
  data_in <- as.data.frame(data_init)
  names(data_in) <- param_names
  
  tstart=proc.time()
  
  p_vac2i<<-c(data_in[["coverage_2006"]],data_in[["coverage_2007"]],data_in[["coverage_2008"]],data_in[["coverage_2009"]],data_in[["coverage_2010"]],data_in[["coverage_2011"]],data_in[["coverage_2012"]],data_in[["coverage_2013"]],
             data_in[["coverage_2014"]],data_in[["coverage_2015"]],data_in[["coverage_2016"]],data_in[["coverage_2017"]],
             data_in[["coverage_2018"]],data_in[["coverage_2019"]],data_in[["coverage_2020"]],data_in[["coverage_2021"]],
             data_in[["coverage_2022"]],data_in[["coverage_2023"]],data_in[["coverage_2024"]],data_in[["coverage_2025"]],
             data_in[["coverage_2026"]],data_in[["coverage_2027"]],data_in[["coverage_2028"]],data_in[["coverage_2029"]],
             data_in[["coverage_2030"]],data_in[["coverage_2031"]],data_in[["coverage_2032"]],data_in[["coverage_2033"]],
             data_in[["coverage_2034"]],data_in[["coverage_2035"]],data_in[["coverage_2036"]],data_in[["coverage_2037"]],
             data_in[["coverage_2038"]],data_in[["coverage_2039"]],data_in[["coverage_2040"]],data_in[["coverage_2041"]],
             data_in[["coverage_2042"]],data_in[["coverage_2043"]]) #proportion receiving TWO doses
  
  p_vac1<<-sqrt(p_vac2i)     #proportion receiving ONE OR MORE doses

  yintro<<-(min(which(p_vac2i>0))+2005)*365

  #recode p_vac2 to zero if considering just one dose  
  if(doses=='single'){
    p_vac2<<-c(rep(0,38))
  } else{p_vac2<<-p_vac2i}
  
  #country and run identifiers
  incom_rank<<-data_in[["income_group"]]

  iso <<- as.character(data_in[['country_code']])
  
  #identify the values of the fitted parameters THIS IS CHANGED
  e0<<- data_in[['life_expectancy']]  #lnum=1, 53.43
  
  #identify the values of parameters that are varide in uncertainty/stochastic runs
  if(run_type=="stochastic"){
    run_id <<- set_run_ids[i]
    #set_run_id<<-set_run_ids[i]
    vac_imm<<-data_in[[paste0("vac_imm_runid",run_id)]]
    mean_age_regress <<- data_in[[paste0('meanage_runid',run_id)]]
  } else {
  run_id<<-"central"
  mean_age_regress <<- data_in[['VIMCPredMeanAgeWk']]  #lnum=1 ,35.59205
  vac_imm<<-data_in[["vac_imm"]]
  }
  
  R0 <<- e0/(mean_age_regress/52)
  q1 <<- R0
  q2<<-q1
  q3<<-q1
  
  if(effect_type=="directeffects"){
    DE_FOI <<- as.character('yes')
  } else {
    DE_FOI <<- as.character('no')
  }
  
  FOI1.pre<<-NA
  FOI2.pre<<-NA
  FOI3.pre<<-NA
  
  #mortality rates-NOTE this the birth rate from 2015 data from the Montagu data
  m_0_59 <<- data_in[["crude_birth_rate"]]/365       # mortality rate under5 
  m_60_up <<- data_in[["crude_birth_rate"]]/365      # mortality rate 5+
  brate <<- data_in[["crude_birth_rate"]]/365        # birth rate
  u1_pop<<-data_in[["p2018u1"]]
  y1_pop<<-data_in[["p2018y1"]]
  y2_pop<<-data_in[["p2018y2"]]
  y3_pop<<-data_in[["p2018y3"]]
  y4_pop<<-data_in[["p2018y4"]]
  pop_5up<<-data_in[["p2018over5"]]
  
  FOI1.pre<<-numeric(length(times))
  FOI2.pre<<-numeric(length(times))
  FOI3.pre<<-numeric(length(times))
  
  #Defining initial conditions of the model
  xstart = c( M_0_1 =1000,S1_0_1 = (u1_pop/12)*2-(2001),I1_0_1 = 1000,S2_0_1 =  0,I2_0_1 = 0,S3_0_1 = 0,I3_0_1 = 0,S4_0_1 = 0,I4_0_1 = 0,R_0_1 = 0,
              #Infant age group:  2 - 3 months old
              M_2_3 = 1000,S1_2_3 = (u1_pop/12)*2-(2001),I1_2_3 = 1000,S2_2_3 = 0,I2_2_3 = 0,S3_2_3 = 0,I3_2_3 = 0, S4_2_3 = 0,I4_2_3 = 0,R_2_3 = 0,V_2_3 = 0,
              #unsuccessfully vaccinated 2 - 3 month olds
              VS1_2_3 = 0,VI1_2_3 = 0,VS2_2_3 = 0,VI2_2_3 = 0,VS3_2_3 = 0,VI3_2_3 = 0, VS4_2_3 = 0,VI4_2_3 = 0,VR_2_3 = 0,
              #Infant age group:  4 - 5 months old
              M_4_5 = 1000,S1_4_5 = (u1_pop/12)*2-(2001),I1_4_5 = 1000,S2_4_5 = 0,I2_4_5 = 0,S3_4_5 = 0,I3_4_5 = 0,S4_4_5 = 0,I4_4_5 = 0,R_4_5 = 0,V_4_5 = 0,
              #unsuccessfully vaccinated 4 - 5 month olds
              VS1_4_5 = 0,VI1_4_5 = 0,VS2_4_5 = 0,VI2_4_5 = 0,VS3_4_5 = 0,VI3_4_5 = 0, VS4_4_5 = 0,VI4_4_5 = 0,VR_4_5 = 0,
              #Infant age group:  6 - 11 months old
              M_6_11 = 1000,S1_6_11 = (u1_pop/12)*6-(2001),I1_6_11 = 1000,S2_6_11 = 0,I2_6_11 = 0,S3_6_11 = 0,I3_6_11 = 0,S4_6_11 = 0,I4_6_11 = 0,R_6_11 = 0,V_6_11 = 0,
              #unsuccessfully vaccinated 6 - 11 month olds
              VS1_6_11 = 0,VI1_6_11 = 0,VS2_6_11 = 0,VI2_6_11 = 0,VS3_6_11 = 0,VI3_6_11 = 0, VS4_6_11 = 0,VI4_6_11 = 0,VR_6_11 = 0,
              #Infant age groups: 12 - 23 months old
              S1_12_23 = y1_pop-1000,I1_12_23 = 1000,S2_12_23 = 0,I2_12_23 = 0,S3_12_23 = 0,I3_12_23 = 0,S4_12_23 = 0,I4_12_23 = 0,R_12_23 = 0,V_12_23 = 0,
              #unsuccessfully vaccinated 6 - 11 month olds
              VS1_12_23 = 0,VI1_12_23 = 0,VS2_12_23 = 0,VI2_12_23 = 0,VS3_12_23 = 0,VI3_12_23 = 0, VS4_12_23 = 0,VI4_12_23 = 0,VR_12_23 = 0,
              #Infant age groups: 24 - 35 months old
              S1_24_35 = y2_pop-1000,I1_24_35 = 1000,S2_24_35 = 0,I2_24_35 = 0,S3_24_35 = 0,I3_24_35 = 0,S4_24_35 = 0,I4_24_35 = 0,R_24_35 = 0,V_24_35 = 0,
              #unsuccessfully vaccinated 6 - 11 month olds
              VS1_24_35 = 0,VI1_24_35 = 0,VS2_24_35 = 0,VI2_24_35 = 0,VS3_24_35 = 0,VI3_24_35 = 0, VS4_24_35 = 0,VI4_24_35 = 0,VR_24_35 = 0,
              #Infant age groups: 36 - 47 months old
              S1_36_47 = y3_pop-1000,I1_36_47 = 1000,S2_36_47 = 0,I2_36_47 = 0,S3_36_47 = 0,I3_36_47 = 0,S4_36_47 = 0,I4_36_47 = 0,R_36_47 = 0,V_36_47 = 0,
              #unsuccessfully vaccinated 6 - 11 month olds
              VS1_36_47 = 0,VI1_36_47 = 0,VS2_36_47 = 0,VI2_36_47 = 0,VS3_36_47 = 0,VI3_36_47 = 0, VS4_36_47 = 0,VI4_36_47 = 0,VR_36_47 = 0,
              #Infant age groups: 48 - 59 months old
              S1_48_59 =y4_pop-1000,I1_48_59 = 1000,S2_48_59 = 0,I2_48_59 = 0,S3_48_59 = 0,I3_48_59 = 0,S4_48_59 = 0,I4_48_59 = 0,R_48_59 = 0, V_48_59 = 0,
              #unsuccessfully vaccinated 6 - 11 month olds
              VS1_48_59 = 0,VI1_48_59 = 0,VS2_48_59 = 0,VI2_48_59 = 0,VS3_48_59 = 0,VI3_48_59 = 0, VS4_48_59 = 0,VI4_48_59 = 0,VR_48_59 = 0,
              #Age group: child = 5-29 years old
              S1_60_up = pop_5up-1000,I1_60_up = 1000,S2_60_up = 0,I2_60_up = 0,S3_60_up = 0,I3_60_up = 0,S4_60_up = 0,I4_60_up = 0,R_60_up = 0, V_60_up = 0,
              #unsuccessfully vaccinated 6 - 11 month olds
              VS1_60_up = 0,VI1_60_up = 0,VS2_60_up = 0,VI2_60_up = 0,VS3_60_up = 0,VI3_60_up = 0, VS4_60_up = 0,VI4_60_up = 0,VR_60_up = 0)
  
  #Line below runs the model
  out_PVM=try(lsoda(xstart, times, RM, parameters, atol=atolv,rtol=rtolv)) #
  #try command catches error from lsoda. If error occurs and lsoda "breaks", we exit the whole optimizer routine with a high objective function value
  if (dim(out_PVM)[1]<length(times)) {cat('!!unresolvable ODE Solver error - triggering early return from optimizer!!'); return(1e10) }
  out_PVM = as.data.frame(out_PVM)

  #saveRDS(out_PVM, paste0('/Users/ani/Documents/phd/research/lopman/vimc_rota/data/03_compile/mod_', mod_time_step,'_', as.numeric(Sys.time()), 
  #  '.rds'))
  
  #saveRDS(out_PVM, paste0('/Users/akraay/Dropbox/RotaDCM-legacy/mod_', mod_time_step,'_', as.numeric(Sys.time()), 
  #                        '.rds'))
  outfile_out_PVM=paste0(patho,'ODE_res_',vax.scenario,iso,run_id,".csv")
  write.csv(out_PVM,outfile_out_PVM)
  
  mean_age_obs_regress <<- data_in[['VIMCPredMeanAgeWk']]
  mean_age_obs_wks <<- data_in[['mean_age_wk']]
  #incom<<-data_in[['Income.group']]
  
  
  #mean_age_pred = mean(SevCases_mean_age, which(times == ystart*365))
  #out_data <- data.frame(mean_age_pred,mean_age_obs_regress,mean_age_obs_wks,iso)
  
  #return(out_data)
  tend=proc.time()
  tdiff=tend-tstart;
  print(sprintf("model run took %f minutes",tdiff[3]/60))
#Closing bracket to loop over run ids
    }
#}

runs=country_params.list[1:73]

tstart=proc.time()
result_df=lapply(runs, FUN= country_runs.fn)
tend=proc.time(); tdiff=tend-tstart



#result_df2 <- ldply(result_df, data.frame)
#write.csv(result_df2,paste0(patho,"pred_mean_age_",vax.scenario,set_run_id,".csv"))
