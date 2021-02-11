#DO THESE FIRST!!!
#Import dataset into R from your directory
#Here is my directory
dataSG_seed_surv <- read.csv("D:/MERDS_2018/merds/Switchgrass/R_data/SG_surv_Seed_germ.csv")
dataSG_trans_surv <- read.csv("D:/MERDS_2018/merds/Switchgrass/R_data/SG_surv_transplant_replace.csv")
data_SG_biomass <- read.csv("D:/MERDS_2018/merds/Switchgrass/R_data/SG_TotalBiomass.csv")
SG_trt <- read.csv("D:/MERDS_2018/merds/Switchgrass/R_data/SG_trt_sheet.csv")
SG_height_combin<- read.csv("D:/MERDS_2018/merds/Switchgrass/R_data/SG_height_time_combine.csv", header = T)
SG_potweight_combin<- read.csv("D:/MERDS_2018/merds/Switchgrass/R_data/SG_potweight_time_combine.csv", header = T)
#let's calculated the percentage change in pot weight from intial potweights
SG_potweight_combin_int=subset(SG_potweight_combin, date=="6_13")
colnames(SG_potweight_combin_int)[3]="int_potweight_g"
SG_potweight_combin_int=SG_potweight_combin_int[,c("Plant_Number","int_potweight_g")]
nrow(SG_potweight_combin_int)

dataSG_delta_potweight=merge(SG_potweight_combin, SG_potweight_combin_int, by="Plant_Number", all.x = T)
summary(dataSG_delta_potweight)
dataSG_delta_potweight$per_ch_potweight=(dataSG_delta_potweight$potweight_g-dataSG_delta_potweight$int_potweight_g)/dataSG_delta_potweight$int_potweight_g
dataSG_delta_potweight_bio=merge(dataSG_delta_potweight, data_SG_biomass, by="Plant_Number", all.x = T)
dataSG_delta_potweight_bio$total_biomass_0=dataSG_delta_potweight_bio$total_biomass
dataSG_delta_potweight_bio$total_biomass_0[is.na(dataSG_delta_potweight_bio$total_biomass)]=0
dataSG_delta_potweight_bio$potweight_g_bio_crt=dataSG_delta_potweight_bio$potweight_g-dataSG_delta_potweight_bio$total_biomass_0
dataSG_delta_potweight_bio$per_ch_potweight_bio_crt=(dataSG_delta_potweight_bio$potweight_g_bio_crt-dataSG_delta_potweight_bio$int_potweight_g)/dataSG_delta_potweight_bio$int_potweight_g
dataSG_delta_potweight2=dataSG_delta_potweight_bio[,!(names(dataSG_delta_potweight_bio) %in% c( "total_biomass_0","root_weight_g","shoot_weight_g","total_biomass","notes_roots","notes_shoots","sandy_root"))]

SG_roottraits<- read.csv("D:/MERDS_2018/merds/Switchgrass/R_data/SG_roottraits.csv", header = T)
summary(SG_roottraits)

SG_inorg_N<- read.csv("D:/MERDS_2018/merds/Switchgrass/R_data/SG_NO3NH4_ug_gdrysoil.csv", header = T)
summary(SG_inorg_N)
SG_start_inorg_N<-read.csv("D:/MERDS_2018/merds/Switchgrass/R_data/SG_start_NO3NH4_ug_gdrysoil.csv", header = T)

#now we need to merge the data file with the trt file


summary(dataSG_seed_surv)
dataSG_seed_surv[,8:14]=NULL
summary(dataSG_trans_surv)
dataSG_trans_surv[,9:14]=NULL
summary(dataSG_delta_potweight)


dataSG_trans_surv_trt=merge(dataSG_trans_surv, SG_trt, by="Plant_Number", all.x = T)
dataSG_seed_surv_trt=merge(dataSG_seed_surv, SG_trt, by="Plant_Number", all.x = T)
dataSG_potweight_trt=merge(dataSG_delta_potweight2, SG_trt, by="Plant_Number", all.x = T)
SG_roottraits_trt=merge(SG_roottraits, SG_trt, by="Plant_Number", all.x = T)
SG_inorg_N_trt=merge(SG_inorg_N, SG_trt, by="Plant_Number", all.x = T)

summary(dataSG_trans_surv_trt)
summary(dataSG_seed_surv_trt)
summary(dataSG_potweight_trt)

library(MASS)
library(reshape2)
library(reshape)
library(multcomp)
library(ggplot2)
library(car)
library(MASS)
library(dplyr)
library(emmeans)
library(lme4)
library(ggpubr)
options(contrasts=c("contr.sum", "contr.poly"))
"%w/o%" <- function(x,y)!('%in%'(x,y))

#####SEED GERMINATION PORTION OF THE EXPERIEMENT BEGINING####

#let's look at germination and number of germinates

#let's look at prop of pots with germinates through time
#THIS ALSO INCLUDES SURVIVAL

dataSG_seed_surv_trt$pot_w_germ=dataSG_seed_surv_trt$num_germinates
dataSG_seed_surv_trt$pot_w_germ[dataSG_seed_surv_trt$pot_w_germ>0]=1
summary(dataSG_seed_surv_trt)
hist(dataSG_seed_surv_trt$pot_w_germ)
dataSG_seed_surv_trt$soil_root=with(dataSG_seed_surv_trt, interaction(soil_status,root_association))





#Time to first germination in a given pot 
dataSG_seed_surv_trt_1=subset(dataSG_seed_surv_trt, pot_w_germ==1)
nrow(dataSG_seed_surv_trt_1)
#764
dataSG_seed_surv_trt_1_pot_g=group_by(dataSG_seed_surv_trt_1, Plant_Number)

dataSG_seed_surv_trt_1_pot_first_germ=summarise_at(dataSG_seed_surv_trt_1_pot_g, "exp_days", min)
colnames(dataSG_seed_surv_trt_1_pot_first_germ)[2]="t0_germ"

dataSG_seed_1st_germ_trt=merge(dataSG_seed_surv_trt_1_pot_first_germ, SG_trt, by="Plant_Number", all.y = T)
dataSG_seed_1st_germ_trt$soil_root=with(dataSG_seed_1st_germ_trt, interaction(soil_status,root_association))
dataSG_seed_1st_germ_trt=subset(dataSG_seed_1st_germ_trt, life_stage=="S")
nrow(dataSG_seed_1st_germ_trt)
summary(dataSG_seed_1st_germ_trt)


day_t0_germ_model= lmer((t0_germ)^-1~precip*soil_root+(1|block), data= dataSG_seed_1st_germ_trt)
qqPlot(resid(day_t0_germ_model))
hist(resid(day_t0_germ_model))
shapiro.test(resid(day_t0_germ_model))
#0.1514


anova(day_t0_germ_model, type=3)
#precip             3.5356  1   0.060065 .  
#soil_root          9.2636  2   0.009737 ** 
emmeans(day_t0_germ_model, ~precip)
emmeans(day_t0_germ_model, pairwise~soil_root)
emmeans(day_t0_germ_model, pairwise~soil_root|precip)


dataSG_seed_1st_germ_trt_na=subset(dataSG_seed_1st_germ_trt, !is.na(t0_germ))
dataSG_seed_1st_germ_trt_precip_soil_g=dataSG_seed_1st_germ_trt_na %>% group_by(soil_root,precip)
day_t0_soil_root_precip=summarise_at(dataSG_seed_1st_germ_trt_precip_soil_g, 
                                   "t0_germ", funs(n(),mean,sd,se=sd(.)/sqrt(n())))

ggplot(day_t0_soil_root_precip, aes(x=precip,y=mean, ymin = mean-se, ymax= mean+se,fill=soil_root))+
  geom_bar(position="dodge",stat="identity", color="black",aes(fill=soil_root))+
  geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+
  scale_fill_manual(values = c("gray", "white", "black"),labels=c("Bulk","Sterile","Rhizo"))+
  scale_x_discrete(labels=c("Ambient","Drought"))+ylab("Days to first germination")+
  geom_text(aes(y=mean+se+1, label=n),position=position_dodge(width=0.9))+theme_bw()

treatment_order=c("S.B","L.B","L.R")
(rate_germination_seedling_p=ggplot(day_t0_soil_root_precip, aes(x=precip,y=mean, ymin = mean-se, ymax= mean+se,fill=factor(soil_root,levels = treatment_order)))+
    geom_bar(position="dodge",stat="identity", color="black",aes(fill=factor(soil_root,levels = treatment_order)))+
    geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+
    scale_fill_manual(values = c( "white","lightgray", "darkgrey"),labels=c("Sterile","Bulk","Rhizo"), limits=treatment_order, name=NULL)+
    scale_x_discrete(labels=c("Ambient","Drought"), name=NULL)+ylab("Days to first germination")+
    geom_text(aes(y=1, label=n),size=10,position=position_dodge(width=0.9))+theme_bw()+scale_y_continuous(limits = c(0,15), breaks = c(3,6,9,12,15))+
    theme(axis.title = element_text(size = 23), axis.text = element_text(size = 23),
          legend.position = "none", panel.grid.major = element_blank(),panel.grid.minor = element_blank()))
#9X8
#800x750


#Quick analyses of germination rate diff
dataSG_seed_1st_germ_trt %>% group_by(soil_root)  %>% summarise_at("t0_germ", list(~mean(.,na.rm = TRUE),~var(.,na.rm = TRUE)))
#   soil_root exp_days
#<fct>        <dbl>
#  1 L.B           5.77
#2 S.B           9.5 
#3 L.R           6.22
9.5-5.77
#3.73
9.5-6.22

dataSG_seed_1st_germ_trt %>% group_by(precip,soil_root)  %>% summarise_at("t0_germ", list(~mean(.,na.rm = TRUE),~sd(.,na.rm = TRUE),~n()))
#drought presence 
7.83-5.2

#ambient presence 
10.8-6.25


dataSG_seed_1st_germ_trt_na=subset(dataSG_seed_1st_germ_trt, t0_germ>0)
nrow(dataSG_seed_1st_germ_trt)
nrow(dataSG_seed_1st_germ_trt_na)
dataSG_seed_1st_germ_trt_na %>% group_by(soil_root,precip)  %>% summarise_at("t0_germ", ~n())



#####Presence Days to first Germination#####
nrow(dataSG_seed_1st_germ_trt)

dataSG_seed_1st_germ_pres=subset(dataSG_seed_1st_germ_trt, root_association =="B")
unique(dataSG_seed_1st_germ_pres$soil_root)
nrow(dataSG_seed_1st_germ_pres)
#60


day_t0_germ_model_pres= lmer((t0_germ)^-1~precip*soil_root+(1|block), data= dataSG_seed_1st_germ_pres)
qqPlot(resid(day_t0_germ_model_pres))
hist(resid(day_t0_germ_model_pres))
shapiro.test(resid(day_t0_germ_model_pres))
#0.8234

item=anova(day_t0_germ_model_pres, type=3)
item$`Pr(>Chisq)`

anova(day_t0_germ_model_pres, type=3)
#precip           0.0122308 0.0122308     1 29.562  4.2813 0.047376 * 
#soil_root        0.0252094 0.0252094     1 28.581  8.8243 0.005968 **



#####Origin Days to first Germination#####
nrow(dataSG_seed_1st_germ_trt)

dataSG_seed_1st_germ_orig=subset(dataSG_seed_1st_germ_trt, soil_status =="L")
unique(dataSG_seed_1st_germ_orig$soil_root)
nrow(dataSG_seed_1st_germ_orig)
#60


day_t0_germ_model_orig= lm((t0_germ)^-2~precip*soil_root+(1|block), data= dataSG_seed_1st_germ_orig)
qqPlot(resid(day_t0_germ_model_orig))
hist(resid(day_t0_germ_model_orig))
shapiro.test(resid(day_t0_germ_model_orig))
#0.1707
boxCox(day_t0_germ_model_orig)

anova(day_t0_germ_model_orig, type=3)

emmeans(day_t0_germ_model_orig, pairwise~soil_root)
emmeans(day_t0_germ_model_orig, pairwise~soil_root|precip)


#let's look at the final number of germinates
#harvested transplants plus final survival

dataSG_seed_surv_trt_g=group_by(dataSG_seed_surv_trt, Plant_Number)
total_transs_harvested=summarise_at(dataSG_seed_surv_trt_g, "num_removed", sum)
colnames(total_transs_harvested)[2]="tot_num_removed"





fin_dataSG_seed_surv_trt=subset(dataSG_seed_surv_trt, date=="7_09")

fin_dataSG_seed_surv_trt$tot_num_germ_harvest=fin_dataSG_seed_surv_trt$num_germinates+total_transs_harvested$tot_num_removed
#this gives the number of plants that may have survived to the end of the experiment but there were some seedlings 
#that died before the experiment ended and were not harvested

#let's make a germinated seeds value
fin_dataSG_seed_surv_trt=merge(fin_dataSG_seed_surv_trt, dataSG_seed_1st_germ_trt[,c("t0_germ","Plant_Number")], by="Plant_Number", all.y = T)


fin_dataSG_seed_surv_trt$germinated=fin_dataSG_seed_surv_trt$t0_germ
fin_dataSG_seed_surv_trt$germinated[fin_dataSG_seed_surv_trt$germinated>0]=1
fin_dataSG_seed_surv_trt$germinated[is.na(fin_dataSG_seed_surv_trt$germinated)]=0
fin_dataSG_seed_surv_trt$soil_root=with(fin_dataSG_seed_surv_trt, interaction(soil_status,root_association))



fin_dataSG_seed_surv_trt %>% group_by(soil_root, precip)  %>% summarise_at("germinated", ~sum(.))

#I am going to make a column with the total number of seedlings (estimated since we did not mark seedlings
#as they germinated) based off if the "germinated" and "tot_num_germ_harvest"


fin_dataSG_seed_surv_trt_sub1=subset(fin_dataSG_seed_surv_trt, germinated>0&tot_num_germ_harvest==0)
nrow(fin_dataSG_seed_surv_trt_sub1)
#12
fin_dataSG_seed_surv_trt_sub1$tot_num_germ=rep(1)
fin_dataSG_seed_surv_trt_sub2=subset(fin_dataSG_seed_surv_trt, germinated==0&tot_num_germ_harvest==0)
nrow(fin_dataSG_seed_surv_trt_sub2)
#31
fin_dataSG_seed_surv_trt_sub2$tot_num_germ=rep(0)
fin_dataSG_seed_surv_trt_sub3=subset(fin_dataSG_seed_surv_trt, germinated==1&tot_num_germ_harvest>0)
nrow(fin_dataSG_seed_surv_trt_sub3)
#47
fin_dataSG_seed_surv_trt_sub3$tot_num_germ=fin_dataSG_seed_surv_trt_sub3$tot_num_germ_harvest
fin_dataSG_seed_surv_trt_tot_seedling=rbind(fin_dataSG_seed_surv_trt_sub1,fin_dataSG_seed_surv_trt_sub2,fin_dataSG_seed_surv_trt_sub3)
nrow(fin_dataSG_seed_surv_trt_tot_seedling)
summary(fin_dataSG_seed_surv_trt_tot_seedling)

write.csv(fin_dataSG_seed_surv_trt_tot_seedling, "D:/MERDS_2018/merds/Switchgrass/R_data/fin_dataSG_seed_surv_trt_tot_seedling.csv")


summary(fin_dataSG_seed_surv_trt)


hist(fin_dataSG_seed_surv_trt$germinated)

#####IS THIS ANALYSIS WRONG####

#pots with germination

seed_germ_model= glm((germinated)~precip*soil_root+(1|block), data= fin_dataSG_seed_surv_trt,family="binomial")
qqPlot(resid(seed_germ_model))
hist(resid(seed_germ_model))
shapiro.test(resid(seed_germ_model))
#p-value = 8.29e-05
#plot(seed_germ_model)


aov(seed_germ_model)
summary(seed_germ_model)
anova(seed_germ_model, type=3)
#precip            10.6579  1   0.001096 **
#soil_root         11.9973  2   0.002482 **
#precip:soil_root   6.6558  2   0.035868 * 
fin_dataSG_seed_surv_trt %>% group_by(precip)  %>% summarise_at("germinated", ~mean(.))
fin_dataSG_seed_surv_trt %>% group_by(soil_root)  %>% summarise_at("germinated", ~mean(.))
fin_dataSG_seed_surv_trt %>% group_by(soil_root,precip)  %>% summarise_at("germinated", ~mean(.))
fin_dataSG_seed_surv_trt   %>% summarise_at("germinated", ~mean(.))
fin_dataSG_seed_surv_trt$soil_root_precip=with(fin_dataSG_seed_surv_trt,interaction(soil_root,precip))
seed_germ_modelP_test= glm(germinated~soil_root_precip+(1|block), data= fin_dataSG_seed_surv_trt)
anova(seed_germ_modelP_test, type=3)
mtl=glht(seed_germ_modelP_test, linfct = mcp(soil_root_precip = "Tukey"))
summary(mtl)

anova(seed_germ_model, type=2)
emmeans(seed_germ_model, pairwise~precip)
emmeans(seed_germ_model, pairwise~soil_root)
emmeans(seed_germ_model, pairwise~soil_root|precip)
"$contrasts
precip = A:
 contrast  estimate       SE  df z.ratio p.value
 L.B - S.B    1.368    0.866 Inf  1.580  0.2544 
 L.B - L.R  -17.153 1599.114 Inf -0.011  0.9999 
 S.B - L.R  -18.520 1599.114 Inf -0.012  0.9999 

precip = D:
 contrast  estimate       SE  df z.ratio p.value
 L.B - S.B    1.220    0.804 Inf  1.517  0.2830 
 L.B - L.R    0.619    0.793 Inf  0.780  0.7152 
 S.B - L.R   -0.601    0.781 Inf -0.770  0.7217 

Results are averaged over the levels of: block 
Results are given on the log odds ratio (not the response) scale. 
P value adjustment: tukey method for comparing a family of 3 estimates "

emmeans(seed_germ_model, pairwise~soil_root*precip, adjust="tukey")
"$contrasts
 contrast      estimate       SE  df z.ratio p.value
 L.B,A - S.B,A    1.368    0.866 Inf  1.580  0.6121 
 L.B,A - L.R,A  -17.153 1599.114 Inf -0.011  1.0000 
 L.B,A - L.B,D    0.749    0.880 Inf  0.851  0.9578 
 L.B,A - S.B,D    1.969    0.878 Inf  2.242  0.2187 
 L.B,A - L.R,D    1.368    0.866 Inf  1.580  0.6121 
 S.B,A - L.R,A  -18.520 1599.114 Inf -0.012  1.0000 
 S.B,A - L.B,D   -0.619    0.793 Inf -0.780  0.9709 
 S.B,A - S.B,D    0.601    0.781 Inf  0.770  0.9726 
 S.B,A - L.R,D    0.000    0.772 Inf  0.000  1.0000 
 L.R,A - L.B,D   17.902 1599.114 Inf  0.011  1.0000 
 L.R,A - S.B,D   19.122 1599.114 Inf  0.012  1.0000 
 L.R,A - L.R,D   18.520 1599.114 Inf  0.012  1.0000 
 L.B,D - S.B,D    1.220    0.804 Inf  1.517  0.6535 
 L.B,D - L.R,D    0.619    0.793 Inf  0.780  0.9709 
 S.B,D - L.R,D   -0.601    0.781 Inf -0.770  0.9726 

Results are averaged over the levels of: block 
Results are given on the log odds ratio (not the response) scale. 
P value adjustment: tukey method for comparing a family of 6 estimates"

seed_germ_model_no_block= glm((germinated)~precip*soil_root, data= fin_dataSG_seed_surv_trt,family="binomial")
qqPlot(resid(seed_germ_model_no_block))
hist(resid(seed_germ_model_no_block))
shapiro.test(resid(seed_germ_model_no_block))
#p-value = 1.366e-07
#plot(seed_germ_model)


aov(seed_germ_model_no_block)
summary(seed_germ_model_no_block)
anova(seed_germ_model_no_block, type=3)
#precip            10.1493  1   0.001444 **
#soil_root         11.3950  2   0.003354 **
#precip:soil_root   6.4784  2   0.039196 *

emmeans(seed_germ_model_no_block, pairwise~soil_root|precip)

AIC(seed_germ_model,seed_germ_model_no_block)

#negative binomial 

seed_germ_model_nb= glmer.nb((germinated)~precip*soil_root+(1|block), data= fin_dataSG_seed_surv_trt)
qqPlot(resid(seed_germ_model_nb))
hist(resid(seed_germ_model_nb))
shapiro.test(resid(seed_germ_model_nb))
#p-value = 6.285e-07

anova(seed_germ_model_nb, type=3)

#Let's try a mixed effects model

#pots with germination
fin_dataSG_seed_surv_trt$germinated
seed_germ_model_mix= glmer((germinated)~precip*soil_root+(1|block), data= fin_dataSG_seed_surv_trt,family="binomial")
qqPlot(resid(seed_germ_model_mix))
hist(resid(seed_germ_model_mix))
shapiro.test(resid(seed_germ_model_mix))
#p-value = 5.294e-07

Anova(seed_germ_model_mix, type=3)

#####IS THIS ANALYSIS WRONG####


#Bars by soil_root

fin_dataSG_seed_surv_trt_soil_root_g=group_by(fin_dataSG_seed_surv_trt, soil_root)
overal_germ_soil_root=summarise_at(fin_dataSG_seed_surv_trt_soil_root_g, "germinated", mean)

ggplot(overal_germ_soil_root, aes(soil_root,germinated))+geom_bar(stat="identity", color="black",aes(fill=soil_root))+
  scale_fill_manual(values = c("gray", "white", "black"),labels=c("Bulk","Sterile","Rhizo"))+
  scale_x_discrete(limits=c("S.B", "L.B", "L.R"),labels=c("Sterile","Bulk","Rhizo"))+ylab("prop Pots with at least one germinate")



fin_dataSG_seed_surv_trt_precip_g=group_by(fin_dataSG_seed_surv_trt, precip)
overal_germ_precip=summarise_at(fin_dataSG_seed_surv_trt_precip_g, "germinated", mean)

ggplot(overal_germ_precip, aes(precip,germinated))+geom_bar(stat="identity", color="black",aes(fill=precip))+
  scale_fill_manual(values = c("gray", "white"),labels=c("Ambient","Drought"))+
  scale_x_discrete(labels=c("Ambient","Drought"))+ylab("prop Pots with at least one germinate")


fin_dataSG_seed_surv_trt_precip_soil_g=fin_dataSG_seed_surv_trt %>% group_by(soil_root,precip)
surv_soil_root_precip=summarise_at(fin_dataSG_seed_surv_trt_precip_soil_g, 
                                   "germinated", funs(n(),mean,sd,se=sd(.)/sqrt(n())))

ggplot(surv_soil_root_precip, aes(x=precip,y=mean,fill=soil_root))+
  geom_bar(position="dodge",stat="identity", color="black",aes(fill=soil_root))+
  scale_fill_manual(values = c("gray", "white", "black"),labels=c("Bulk","Sterile","Rhizo"))+
  scale_x_discrete(labels=c("Ambient","Drought"))+ylab("prop Pots with at least one germinate")+
  theme_bw()

#Let's use CHi-Sq to analyses the proportional difference

seed_germ_matrix_chi_sq=fin_dataSG_seed_surv_trt %>% group_by(soil_root,precip)  %>% summarise_at("germinated", ~sum(.))
head(seed_germ_matrix_chi_sq)

seed_germ_matrix_chi_sq_c=cast(seed_germ_matrix_chi_sq, precip~ soil_root)
row.names(seed_germ_matrix_chi_sq_c)=seed_germ_matrix_chi_sq_c$precip
seed_germ_matrix_chi_sq_c$precip=NULL
seed_germ_chisq_model <- chisq.test(seed_germ_matrix_chi_sq_c)
seed_germ_chisq_model$observed


#combined seed germination

fin_dataSG_seed_surv_trt %>% group_by(precip) %>% summarise_at("germinated", list(~n(),~mean(.),~sd(.),se=~sd(.)/sqrt(n())))
(0.533-0.778)/0.778

fin_dataSG_seed_surv_trt %>% group_by(soil_root) %>% summarise_at("germinated", list(~n(),~mean(.),~sd(.),se=~sd(.)/sqrt(n())))
#presence
(0.467-0.733)/0.733

fin_dataSG_seed_surv_trt %>% group_by(precip,soil_root) %>% summarise_at("germinated", list(~n(),~mean(.),~sd(.),se=~sd(.)/sqrt(n())))
#drought presence
(0.4-0.667)/0.667
#ambient presence
(0.533-0.8)/0.8

fin_dataSG_seed_surv_trt_g=fin_dataSG_seed_surv_trt %>% group_by(soil_root,precip)
surv_soil_root_precip=summarise_at(fin_dataSG_seed_surv_trt_g, 
                                               "germinated", list(~n(),~mean(.),~sd(.),se=~sd(.)/sqrt(n())))

treatment_order=c("S.B","L.B","L.R")
(pot_germination_seedling_p=ggplot(surv_soil_root_precip, aes(x=precip,y=mean,fill=factor(soil_root,levels = treatment_order)))+
  geom_bar(position="dodge",stat="identity", color="black",aes(fill=factor(soil_root,levels = treatment_order)))+
  scale_fill_manual(values = c( "white","lightgray", "darkgrey"),labels=c("Sterile","Bulk","Rhizo"), limits=treatment_order, name=NULL)+
  scale_x_discrete(labels=c("Ambient","Drought"), name=NULL)+ylab("Pots with at least one germinate")+
  geom_text(aes(y=0.07, label=n),size=10,position=position_dodge(width=0.9))+theme_bw()+scale_y_continuous(limits = c(0,1.1), breaks = c(0,0.2,0.4,0.6,0.8,1))+
  theme(axis.title = element_text(size = 23), axis.text = element_text(size = 23),
        legend.position = c(0.85,.9), legend.text=element_text(size=20),
        legend.background = element_rect(size=0.5,linetype="solid",colour ="black")))

surv_soil_root_precip$mean_per=surv_soil_root_precip$mean*100

(pot_germination_per_seedling_p=ggplot(surv_soil_root_precip, aes(x=precip,y=mean_per,fill=factor(soil_root,levels = treatment_order)))+
    geom_bar(position="dodge",stat="identity", color="black",aes(fill=factor(soil_root,levels = treatment_order)))+
    scale_fill_manual(values = c( "white","lightgray", "darkgrey"),labels=c("Sterile","Bulk","Rhizo"), limits=treatment_order, name=NULL)+
    scale_x_discrete(labels=c("Ambient","Drought"), name=NULL)+ylab("% pots with at least one germinate")+
    geom_text(aes(y=7, label=n),size=10,position=position_dodge(width=0.9))+theme_bw()+scale_y_continuous(limits = c(0,120), breaks = c(0,20,40,60,80,100))+
    theme(axis.title = element_text(size = 23), axis.text = element_text(size = 23),
          legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank()))
#9X8
#800x750
#####Combined germination rate of germ graph#####
ggarrange(pot_germination_seedling_p,rate_germination_seedling_p,ncol = 2,  legend = "none")
#15x7.38


#####Presence pots with germination#####
fin_dataSG_seed_surv_pres=subset(fin_dataSG_seed_surv_trt, root_association =="B")
nrow(fin_dataSG_seed_surv_pres)

seed_germ_model_pres= glm((germinated)~precip*soil_root+(1|block), data= fin_dataSG_seed_surv_pres,family="binomial")
qqPlot(resid(seed_germ_model_pres))
hist(resid(seed_germ_model_pres))
shapiro.test(resid(seed_germ_model_pres))
#p-value = 8.202e-05
#plot(seed_germ_model)


aov(seed_germ_model_pres)
summary(seed_germ_model_pres)
anova(seed_germ_model_pres, type=3)
#soil_root          5.2465  1    0.02199 *
fin_dataSG_seed_surv_trt %>% group_by(precip)  %>% summarise_at("germinated", ~mean(.))
fin_dataSG_seed_surv_trt %>% group_by(soil_root)  %>% summarise_at("germinated", ~mean(.))
fin_dataSG_seed_surv_trt %>% group_by(soil_root,precip)  %>% summarise_at("germinated", ~mean(.))
fin_dataSG_seed_surv_trt   %>% summarise_at("germinated", ~mean(.))
fin_dataSG_seed_surv_trt$soil_root_precip=with(fin_dataSG_seed_surv_trt,interaction(soil_root,precip))


seed_germ_model_pres_no_block= glm((germinated)~precip*soil_root, data= fin_dataSG_seed_surv_pres,family="binomial")
qqPlot(resid(seed_germ_model_pres_no_block))
hist(resid(seed_germ_model_pres_no_block))
shapiro.test(resid(seed_germ_model_pres_no_block))
#p-value = 2.145e-07
#plot(seed_germ_model)


aov(seed_germ_model_pres_no_block)
summary(seed_germ_model_pres_no_block)
anova(seed_germ_model_pres_no_block, type=3)
#soil_root          4.6207  1    0.03159 *

emmeans(seed_germ_model_pres_no_block, pairwise~soil_root|precip)

emmeans(seed_germ_model_pres_no_block, pairwise~soil_root*precip,adjust="none")

AIC(seed_germ_model_pres,seed_germ_model_pres_no_block)


#####Origin pots with germination#####
fin_dataSG_seed_surv_orig=subset(fin_dataSG_seed_surv_trt, soil_status =="L")
nrow(fin_dataSG_seed_surv_orig)

seed_germ_model_orig= glm((germinated)~precip*soil_root+(1|block), data= fin_dataSG_seed_surv_orig,family="binomial")
qqPlot(resid(seed_germ_model_orig))
hist(resid(seed_germ_model_orig))
shapiro.test(resid(seed_germ_model_orig))
#p-value = 0.001212
#plot(seed_germ_model)


aov(seed_germ_model_orig)
summary(seed_germ_model_orig)
anova(seed_germ_model_orig, type=3)
#precip:soil_root   5.1778  1  0.0228776 *  
fin_dataSG_seed_surv_trt %>% group_by(precip)  %>% summarise_at("germinated", ~mean(.))
fin_dataSG_seed_surv_trt %>% group_by(soil_root)  %>% summarise_at("germinated", ~mean(.))
fin_dataSG_seed_surv_trt %>% group_by(soil_root,precip)  %>% summarise_at("germinated", ~mean(.))
fin_dataSG_seed_surv_trt   %>% summarise_at("germinated", ~mean(.))
fin_dataSG_seed_surv_trt$soil_root_precip=with(fin_dataSG_seed_surv_trt,interaction(soil_root,precip))


seed_germ_model_orig_no_block= glm((germinated)~precip*soil_root, data= fin_dataSG_seed_surv_orig,family="binomial")
qqPlot(resid(seed_germ_model_orig_no_block))
hist(resid(seed_germ_model_orig_no_block))
shapiro.test(resid(seed_germ_model_orig_no_block))
#p-value = 1.479e-06
#plot(seed_germ_model)


aov(seed_germ_model_orig_no_block)
summary(seed_germ_model_orig_no_block)
anova(seed_germ_model_orig_no_block, type=3)
#precip:soil_root   4.9494  1   0.026100 * 

emmeans(seed_germ_model_orig_no_block, pairwise~soil_root*precip,adjust="none")

AIC(seed_germ_model_orig,seed_germ_model_orig_no_block)


#now lets look at the total number of germinates
#####Diagnostic plots for probablity disturbutions####

#https://ase.tufts.edu/gsc/gradresources/guidetomixedmodelsinr/mixed%20model%20guide.html
require(MASS)

#I want no zeros
fin_dataSG_seed_surv_trt_tot_seedling$tot_num_germ.t <- fin_dataSG_seed_surv_trt_tot_seedling$tot_num_germ + 1
min(fin_dataSG_seed_surv_trt_tot_seedling$tot_num_germ.t)

# lnorm means lognormal
qqp(fin_dataSG_seed_surv_trt_tot_seedling$tot_num_germ.t, "lnorm")

# qqp requires estimates of the parameters of the negative binomial, Poisson
# and gamma distributions. You can generate estimates using the fitdistr
# function. Save the output and extract the estimates of each parameter as I
# have shown below.
nbinom <- fitdistr(fin_dataSG_seed_surv_trt_tot_seedling$tot_num_germ.t, "Negative Binomial")
qqp(fin_dataSG_seed_surv_trt_tot_seedling$tot_num_germ.t, "nbinom", size = nbinom$estimate[[1]], mu = nbinom$estimate[[2]])

poisson <- fitdistr(fin_dataSG_seed_surv_trt_tot_seedling$tot_num_germ.t, "Poisson")
qqp(fin_dataSG_seed_surv_trt_tot_seedling$tot_num_germ.t, "pois", lambda=poisson$estimate)

gamma <- fitdistr(fin_dataSG_seed_surv_trt_tot_seedling$tot_num_germ.t, "gamma")
qqp(fin_dataSG_seed_surv_trt_tot_seedling$tot_num_germ.t, "gamma", shape = gamma$estimate[[1]], rate = gamma$estimate[[2]])

hist(fin_dataSG_seed_surv_trt_tot_seedling$tot_num_germ)

seed_num_germ_model= lm((tot_num_germ+1)^(1/3)~precip*soil_root+(1|block), data= fin_dataSG_seed_surv_trt_tot_seedling)
qqPlot(resid(seed_num_germ_model))
hist(resid(seed_num_germ_model))
shapiro.test(resid(seed_num_germ_model))
#p-value = 0.01155

anova(seed_num_germ_model, type=3)
#precip             0.245  1    8.8857  0.003806 ** 
#soil_root          0.273  2    4.9524  0.009382 ** 
emmeans(seed_num_germ_model, ~precip)
emmeans(seed_num_germ_model, pairwise~soil_root)

emmeans(seed_num_germ_model, pairwise~soil_root|precip)


seed_num_germ_model_no_block= lm((tot_num_germ+1)^(1/3)~precip*soil_root, data= fin_dataSG_seed_surv_trt_tot_seedling)
qqPlot(resid(seed_num_germ_model_no_block))
hist(resid(seed_num_germ_model_no_block))
shapiro.test(resid(seed_num_germ_model_no_block))
#p-value = 0.0005655

anova(seed_num_germ_model_no_block, type=3)
#precip             0.245  1    8.5377  0.004467 ** 
#soil_root          0.273  2    4.7585  0.011024 * 

AIC(seed_num_germ_model,seed_num_germ_model_no_block)

seed_num_germ_model_no_block.2= lm(log(tot_num_germ+1)~precip*soil_root, data= fin_dataSG_seed_surv_trt_tot_seedling)
qqPlot(resid(seed_num_germ_model_no_block.2))
hist(resid(seed_num_germ_model_no_block.2))
shapiro.test(resid(seed_num_germ_model_no_block.2))
#p-value = 0.0006849

anova(seed_num_germ_model_no_block.2, type=3)


emmeans(seed_num_germ_model_no_block, pairwise~soil_root|precip)

fin_dataSG_seed_surv_trt_tot_seedling_g=fin_dataSG_seed_surv_trt_tot_seedling %>% group_by(soil_root,precip)
germination_tot_soil_root_precip=summarise_at(fin_dataSG_seed_surv_trt_tot_seedling_g, 
                                   "tot_num_germ", list(~n(),~mean(.),~sd(.),se=~sd(.)/sqrt(n())))
germination_tot_soil_root_precip$mean_per=(germination_tot_soil_root_precip$mean/10)*100
germination_tot_soil_root_precip$se_per=(germination_tot_soil_root_precip$se/10)*100
treatment_order=c("S.B","L.B","L.R")
(germination_tot_seedling_p=ggplot(germination_tot_soil_root_precip, aes(x=precip,y=mean_per, ymin = mean_per-se_per, ymax= mean_per+se_per,fill=factor(soil_root,levels = treatment_order)))+
    geom_bar(position="dodge",stat="identity", color="black",aes(fill=factor(soil_root,levels = treatment_order)))+
    geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+
    scale_fill_manual(values = c( "white","lightgray", "darkgrey"),labels=c("Sterile","Bulk","Rhizo"), limits=treatment_order, name=NULL)+
    scale_x_discrete(labels=c("Ambient","Drought"), name=NULL)+ylab("Percentage germination per pot")+theme_bw()+ylim(c(0,21))+
    theme(axis.title = element_text(size = 23), axis.text = element_text(size = 23),
          legend.position = c(0.85,.9), legend.text=element_text(size=20),
          legend.background = element_rect(size=0.5,linetype="solid",colour ="black")))


#####Combined germination rate and number of germ graph#####
ggarrange(germination_tot_seedling_p,rate_germination_seedling_p,ncol = 2,  legend = "none")
#15x7.38



#####Presence total number of germinates####
fin_dataSG_seed_surv_pres_tot_seedling=subset(fin_dataSG_seed_surv_trt_tot_seedling, root_association =="B")
nrow(fin_dataSG_seed_surv_pres_tot_seedling)
#60
hist(fin_dataSG_seed_surv_trt_tot_seedling$tot_num_germ)

seed_num_germ_model_pres= lm((tot_num_germ)~precip*soil_root+(1|block), data= fin_dataSG_seed_surv_pres_tot_seedling)
qqPlot(resid(seed_num_germ_model_pres))
hist(resid(seed_num_germ_model_pres))
shapiro.test(resid(seed_num_germ_model_pres))
#p-value = 0.1838

anova(seed_num_germ_model_pres, type=3)
#soil_root         3.750  1  6.6402  0.01285 *  
#(1|block)  5.567  4  2.4642  0.05639 .  
emmeans(seed_num_germ_model_pres, ~precip)
emmeans(seed_num_germ_model_pres, pairwise~soil_root)

emmeans(seed_num_germ_model_pres, pairwise~soil_root*precip,adjust="none")


seed_num_germ_model_pres_no_block= lm((tot_num_germ)~precip*soil_root, data= fin_dataSG_seed_surv_pres_tot_seedling)
qqPlot(resid(seed_num_germ_model_pres_no_block))
hist(resid(seed_num_germ_model_pres_no_block))
shapiro.test(resid(seed_num_germ_model_pres_no_block))
#p-value = 0.0245
boxCox(seed_num_germ_model_pres_no_block)


anova(seed_num_germ_model_pres_no_block, type=3)
#soil_root         3.750  1  6.0115  0.01736 * 

AIC(seed_num_germ_model_pres,seed_num_germ_model_pres_no_block)

seed_num_germ_model_no_block.2= lm(log(tot_num_germ+1)~precip*soil_root, data= fin_dataSG_seed_surv_trt_tot_seedling)
qqPlot(resid(seed_num_germ_model_no_block.2))
hist(resid(seed_num_germ_model_no_block.2))
shapiro.test(resid(seed_num_germ_model_no_block.2))
#p-value = 0.0006849

anova(seed_num_germ_model_no_block.2, type=3)


emmeans(seed_num_germ_model_no_block, pairwise~soil_root|precip)


#####Origin total number of germinates####
fin_dataSG_seed_surv_orig_tot_seedling=subset(fin_dataSG_seed_surv_trt_tot_seedling, soil_status =="L")
nrow(fin_dataSG_seed_surv_orig_tot_seedling)
#60
hist(fin_dataSG_seed_surv_orig_tot_seedling$tot_num_germ)

seed_num_germ_model_orig= lm(sqrt(tot_num_germ+1)~precip*soil_root+(1|block), data= fin_dataSG_seed_surv_orig_tot_seedling)
qqPlot(resid(seed_num_germ_model_orig))
hist(resid(seed_num_germ_model_orig))
shapiro.test(resid(seed_num_germ_model_orig))
#p-value = 0.2825

anova(seed_num_germ_model_orig, type=3)
#precip             0.815  1   10.3158  0.002264 ** 
emmeans(seed_num_germ_model_orig, ~precip)
emmeans(seed_num_germ_model_orig, pairwise~soil_root)

emmeans(seed_num_germ_model_orig, pairwise~soil_root*precip,adjust="none")


seed_num_germ_model_orig_no_block= lm(sqrt(tot_num_germ+1)~precip*soil_root, data= fin_dataSG_seed_surv_orig_tot_seedling)
qqPlot(resid(seed_num_germ_model_orig_no_block))
hist(resid(seed_num_germ_model_orig_no_block))
shapiro.test(resid(seed_num_germ_model_orig_no_block))
#p-value = 0.08012
boxCox(seed_num_germ_model_orig_no_block)


anova(seed_num_germ_model_orig_no_block, type=3)
#precip             0.815  1   10.1178  0.002395 ** 

emmeans(seed_num_germ_model_orig_no_block, pairwise~soil_root*precip,adjust="none")


AIC(seed_num_germ_model_orig,seed_num_germ_model_orig_no_block)

seed_num_germ_model_no_block.2= lm(log(tot_num_germ+1)~precip*soil_root, data= fin_dataSG_seed_surv_trt_tot_seedling)
qqPlot(resid(seed_num_germ_model_no_block.2))
hist(resid(seed_num_germ_model_no_block.2))
shapiro.test(resid(seed_num_germ_model_no_block.2))
#p-value = 0.0006849

anova(seed_num_germ_model_no_block.2, type=3)


emmeans(seed_num_germ_model_no_block, pairwise~soil_root|precip)




#now lets look at the total number of germinates given there was at least one germinate

fin_dataSG_seed_surv_trt_tot_seedling_w_germ=subset(fin_dataSG_seed_surv_trt_tot_seedling, germinated>0)

hist(fin_dataSG_seed_surv_trt_tot_seedling_w_germ$tot_num_germ)

seed_num_germ_w_germ_model= glm((tot_num_germ)~precip*soil_root+(1|block), data= fin_dataSG_seed_surv_trt_tot_seedling_w_germ,family = "poisson")
qqPlot(resid(seed_num_germ_w_germ_model))
hist(resid(seed_num_germ_w_germ_model))
shapiro.test(resid(seed_num_germ_w_germ_model))
#p-value = 0.01552
anova(seed_num_germ_w_germ_model, type=3)
#nada sig
emmeans(seed_num_germ_w_germ_model, ~precip)
emmeans(seed_num_germ_w_germ_model, pairwise~soil_root)

fin_dataSG_seed_surv_trt_tot_seedling_w_germ_soil_root_g=group_by(fin_dataSG_seed_surv_trt_tot_seedling_w_germ, soil_root)
overal_num_germ_soil_root=summarise_at(fin_dataSG_seed_surv_trt_tot_seedling_w_germ_soil_root_g, "tot_num_germ", funs(mean,sd))

ggplot(overal_num_germ_soil_root, aes(soil_root,mean))+geom_bar(stat="identity", color="black",aes(fill=soil_root))+
  scale_fill_manual(values = c("gray", "white", "black"),labels=c("Bulk","Sterile","Rhizo"))+
  scale_x_discrete(limits=c("S.B", "L.B", "L.R"),labels=c("Sterile","Bulk","Rhizo"))+ylab("Num germ with zeros")



fin_dataSG_seed_surv_w_germ_trt_precip_g=group_by(fin_dataSG_seed_surv_trt_tot_seedling_w_germ, precip)
overal_w_germ_precip=summarise_at(fin_dataSG_seed_surv_w_germ_trt_precip_g, "germinated", mean)

ggplot(overal_w_germ_precip, aes(precip,germinated))+geom_bar(stat="identity", color="black",aes(fill=precip))+
  scale_fill_manual(values = c("gray", "white"),labels=c("Ambient","Drought"))+
  scale_x_discrete(labels=c("Ambient","Drought"))+ylab("At least one germinate")


###


data_SG_biomass_seed_surv_trt=merge(fin_dataSG_seed_surv_trt, data_SG_biomass, by="Plant_Number", all.x = T)

#Survival given germination
fin_dataSG_seed_surv_trt
fin_dataSG_biomass_seed_surv_trt_w_germ=subset(data_SG_biomass_seed_surv_trt, germinated>0)
summary(fin_dataSG_biomass_seed_surv_trt_w_germ)



fin_dataSG_biomass_seed_surv_trt_w_germ$surv_germ=fin_dataSG_biomass_seed_surv_trt_w_germ$total_biomass
fin_dataSG_biomass_seed_surv_trt_w_germ$surv_germ[is.na(fin_dataSG_biomass_seed_surv_trt_w_germ$total_biomass)]=0
fin_dataSG_biomass_seed_surv_trt_w_germ$surv_germ[fin_dataSG_biomass_seed_surv_trt_w_germ$total_biomass>0]=1
hist(fin_dataSG_biomass_seed_surv_trt_w_germ$surv_germ)


seed_surv_germ_model= glm(surv_germ~precip*soil_root+(1|block), data= fin_dataSG_biomass_seed_surv_trt_w_germ, family="binomial")
qqPlot(resid(seed_surv_germ_model))
hist(resid(seed_surv_germ_model))
shapiro.test(resid(seed_surv_germ_model))
#4.631e-06
#plot(seed_surv_germ_model)

anova(seed_surv_germ_model, type=3)
#soil_root          9.8880  2   0.007126 **
#(1|block)  11.6063  4   0.020533 * 
#precip:soil_root   4.3070  2   0.116077 
emmeans(seed_surv_germ_model, ~precip)
emmeans(seed_surv_germ_model, pairwise~soil_root)
emmeans(seed_surv_germ_model, pairwise~soil_root|precip)

seed_surv_germ_model_no_block= glm(surv_germ~precip*soil_root, data= fin_dataSG_biomass_seed_surv_trt_w_germ, family="binomial")
qqPlot(resid(seed_surv_germ_model_no_block))
hist(resid(seed_surv_germ_model_no_block))
shapiro.test(resid(seed_surv_germ_model_no_block))
#6.334e-09
#plot(seed_surv_germ_model)

anova(seed_surv_germ_model_no_block, type=3)
#soil_root          7.1039  2    0.02867 *

emmeans(seed_surv_germ_model_no_block, pairwise~soil_root|precip)

AIC(seed_surv_germ_model,seed_surv_germ_model_no_block)

fin_dataSG_biomass_seed_surv_trt_w_germ %>% group_by(soil_root) %>% summarise_at("surv_germ", list(~n(),~mean(.),~sd(.),se=~sd(.)/sqrt(n())))

#Live bulk 
(0.5-0.864)/0.864
#-0.4212963


fin_dataSG_biomass_seed_surv_trt_w_germ %>% group_by(precip,soil_root) %>% summarise_at("surv_germ", list(~n(),~mean(.),~sd(.),se=~sd(.)/sqrt(n())))

#Drought Presence
(0.333-0.9)/0.9

#combined seed survival

fin_dataSG_biomass_seed_surv_trt_w_germ_g=fin_dataSG_biomass_seed_surv_trt_w_germ %>% group_by(soil_root,precip)
seed_surv_soil_root_precip=summarise_at(fin_dataSG_biomass_seed_surv_trt_w_germ_g, 
                                   "surv_germ", list(~n(),~mean(.),~sd(.),se=~sd(.)/sqrt(n())))

#Live bulk 
(0.9-0.333)/0.333
#Live Rhizosphere
(0.875-0.333)/0.333

(seedling_suv_p=ggplot(seed_surv_soil_root_precip, aes(x=precip,y=mean,fill=factor(soil_root,levels = treatment_order)))+
  geom_bar(position="dodge",stat="identity", color="black",aes(fill=factor(soil_root,levels = treatment_order)))+
  scale_fill_manual(values = c( "white","light gray", "dark grey"),labels=c("Sterile","Bulk","Rhizo"), limits=treatment_order, name=NULL)+
  scale_x_discrete(labels=c("Ambient","Drought"), name=NULL)+ylab("prop Seedling survival \ngiven gemination")+ylim(c(0,1))+
  geom_text(aes(y=0.07, label=n),size=10,position=position_dodge(width=0.9))+theme_bw()+
  theme(axis.title = element_text(size = 23), axis.text = element_text(size = 23),
        legend.position = c(0.85,.9), legend.text=element_text(size=20),
        legend.background = element_rect(size=0.5,linetype="solid",colour ="black")))


seed_surv_soil_root_precip$mean_per=seed_surv_soil_root_precip$mean*100
(seedling_suv_p2=ggplot(seed_surv_soil_root_precip, aes(x=precip,y=mean_per,fill=factor(soil_root,levels = treatment_order)))+
    geom_bar(position="dodge",stat="identity", color="black",aes(fill=factor(soil_root,levels = treatment_order)))+
    scale_fill_manual(values = c( "white","lightgray", "darkgrey"),labels=c("Sterile","Bulk","Rhizo"), limits=treatment_order, name=NULL)+
    scale_x_discrete(labels=c("Ambient","Drought"), name=NULL)+ylab("% seedling survival \ngiven gemination")+
    scale_y_continuous(limits = c(0,120), breaks = c(0,20,40,60,80,100))+
    geom_text(aes(y=7, label=n),size=10,position=position_dodge(width=0.9))+theme_bw()+
    theme(axis.title = element_text(size = 23), axis.text = element_text(size = 23),
          legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank()))
#800x750
####Presence survival given germination####
fin_dataSG_biomass_seed_surv_pres_w_germ=subset(fin_dataSG_biomass_seed_surv_trt_w_germ, root_association =="B")
nrow(fin_dataSG_biomass_seed_surv_pres_w_germ)
#36
unique(fin_dataSG_biomass_seed_surv_pres_w_germ$soil_root)

seed_surv_germ_model_pres= glm(surv_germ~precip*soil_root+(1|block), data= fin_dataSG_biomass_seed_surv_pres_w_germ, family="binomial")
qqPlot(resid(seed_surv_germ_model_pres))
hist(resid(seed_surv_germ_model_pres))
shapiro.test(resid(seed_surv_germ_model_pres))
#0.0001507
#plot(seed_surv_germ_model_pres)

anova(seed_surv_germ_model_pres, type=3)
#soil_root          6.3748  1    0.01158 *
emmeans(seed_surv_germ_model, ~precip)
emmeans(seed_surv_germ_model, pairwise~soil_root)
emmeans(seed_surv_germ_model, pairwise~soil_root|precip)

seed_surv_germ_model_pres_no_block= glm(surv_germ~precip*soil_root, data= fin_dataSG_biomass_seed_surv_pres_w_germ, family="binomial")
qqPlot(resid(seed_surv_germ_model_pres_no_block))
hist(resid(seed_surv_germ_model_pres_no_block))
shapiro.test(resid(seed_surv_germ_model_pres_no_block))
#6.334e-09
#plot(seed_surv_germ_model)

anova(seed_surv_germ_model_pres_no_block, type=3)
#soil_root          6.0809  1    0.01366 *

emmeans(seed_surv_germ_model_pres_no_block, pairwise~soil_root|precip)
emmeans(seed_surv_germ_model_pres_no_block, pairwise~soil_root*precip, adjust="none")


AIC(seed_surv_germ_model_pres,seed_surv_germ_model_pres_no_block)


####Origin survival given germination####
fin_dataSG_biomass_seed_surv_orig_w_germ=subset(fin_dataSG_biomass_seed_surv_trt_w_germ, soil_status =="L")
nrow(fin_dataSG_biomass_seed_surv_orig_w_germ)
#45
unique(fin_dataSG_biomass_seed_surv_orig_w_germ$soil_root)

seed_surv_germ_model_orig= glm(surv_germ~precip*soil_root+(1|block), data= fin_dataSG_biomass_seed_surv_orig_w_germ, family="binomial")
qqPlot(resid(seed_surv_germ_model_orig))
hist(resid(seed_surv_germ_model_orig))
shapiro.test(resid(seed_surv_germ_model_orig))
#2.128e-06
#plot(seed_surv_germ_model_pres)

anova(seed_surv_germ_model_orig, type=3)
#nada sig
emmeans(seed_surv_germ_model_orig, ~precip)
emmeans(seed_surv_germ_model_orig, pairwise~soil_root)
emmeans(seed_surv_germ_model_orig, pairwise~soil_root|precip)

seed_surv_germ_model_orig_no_block= glm(surv_germ~precip*soil_root, data= fin_dataSG_biomass_seed_surv_orig_w_germ, family="binomial")
qqPlot(resid(seed_surv_germ_model_orig_no_block))
hist(resid(seed_surv_germ_model_orig_no_block))
shapiro.test(resid(seed_surv_germ_model_orig_no_block))
#6.334e-09
#plot(seed_surv_germ_model)

anova(seed_surv_germ_model_orig_no_block, type=3)
#nada sig

emmeans(seed_surv_germ_model_orig_no_block, pairwise~soil_root|precip)
emmeans(seed_surv_germ_model_orig_no_block, pairwise~soil_root*precip, adjust="none")


AIC(seed_surv_germ_model_orig,seed_surv_germ_model_orig_no_block)



#now lets look at total biomass of germinated 
#fin_dataSG_biomass_seed_surv_trt_w_germ_test=subset(data_SG_biomass_seed_surv_trt, !is.na(total_biomass))
#nrow(fin_dataSG_biomass_seed_surv_trt_w_germ_test)
#45
fin_dataSG_biomass_seed_surv_trt_w_germ=subset(fin_dataSG_biomass_seed_surv_trt_w_germ, germinated>0)
fin_dataSG_biomass_seed_surv_trt_given_germ=subset(fin_dataSG_biomass_seed_surv_trt_w_germ, !is.na(total_biomass))
nrow(fin_dataSG_biomass_seed_surv_trt_given_germ)
#45
summary(fin_dataSG_biomass_seed_surv_trt_given_germ)
fin_dataSG_biomass_seed_surv_trt_given_germ%>%group_by(precip,soil_root)%>%summarise_at("total_biomass",~n())
seed_biomas_model= lm(log(total_biomass)~precip*soil_root+(1|block), data= fin_dataSG_biomass_seed_surv_trt_given_germ)
qqPlot(resid(seed_biomas_model))
hist(resid(seed_biomas_model))
boxCox(seed_biomas_model)
shapiro.test(resid(seed_biomas_model))
#p-value = 0.3121
anova(seed_biomas_model, type=3)
#precip:soil_root  10.855  2   5.7748  0.006802 ** 
#soil_root          4.776  2   2.5409  0.093244 . 

emmeans(seed_biomas_model, ~precip)
emmeans(seed_biomas_model, pairwise~soil_root)
emmeans(seed_biomas_model, pairwise~soil_root|precip)
#$contrasts
#precip = A:
#  contrast  estimate    SE df t.ratio p.value
#L.B - S.B    0.462 0.542 35  0.853  0.6732 
#L.B - L.R   -0.396 0.418 35 -0.947  0.6148 
#S.B - L.R   -0.858 0.520 35 -1.650  0.2385 

#precip = D:
#  contrast  estimate    SE df t.ratio p.value
#L.B - S.B   -2.627 0.813 35 -3.231  0.0074 
#L.B - L.R   -0.135 0.496 35 -0.273  0.9598 
#S.B - L.R    2.491 0.864 35  2.883  0.0179 

emmeans(seed_biomas_model, pairwise~soil_root*precip, adjust="fdr")

"$contrasts
 contrast      estimate    SE df t.ratio p.value
L.B,A - S.B,A    0.529 0.482 34  1.097  0.3823 
L.B,A - L.R,A   -0.364 0.372 34 -0.980  0.4177 
L.B,A - L.B,D    0.775 0.403 34  1.924  0.1176 
L.B,A - S.B,D   -1.870 0.708 34 -2.642  0.0371 
L.B,A - L.R,D    0.207 0.459 34  0.450  0.6552 
S.B,A - L.R,A   -0.893 0.462 34 -1.932  0.1176 
S.B,A - L.B,D    0.247 0.484 34  0.509  0.6552 
S.B,A - S.B,D   -2.399 0.771 34 -3.113  0.0258 
S.B,A - L.R,D   -0.322 0.532 34 -0.606  0.6330 
L.R,A - L.B,D    1.140 0.381 34  2.989  0.0258 
L.R,A - S.B,D   -1.506 0.704 34 -2.140  0.0990 
L.R,A - L.R,D    0.571 0.440 34  1.299  0.3379 
L.B,D - S.B,D   -2.646 0.722 34 -3.663  0.0126 
L.B,D - L.R,D   -0.569 0.461 34 -1.234  0.3383 
S.B,D - L.R,D    2.077 0.779 34  2.667  0.0371 

Results are averaged over the levels of: block 
Results are given on the log (not the response) scale. 
P value adjustment: fdr method for 15 tests 
"
seed_biomas_model_no_block= lm(log(total_biomass)~precip*soil_root, data= fin_dataSG_biomass_seed_surv_trt_given_germ)
qqPlot(resid(seed_biomas_model_no_block))
hist(resid(seed_biomas_model_no_block))
boxCox(seed_biomas_model_no_block)
shapiro.test(resid(seed_biomas_model_no_block))
#p-value = 0.457
anova(seed_biomas_model_no_block, type=3)
#soil_root          7.220  2   3.9802 0.0267300 *  
#precip:soil_root  18.153  2  10.0075 0.0003104 ***

emmeans(seed_biomas_model_no_block, pairwise~soil_root|precip)

AIC(seed_biomas_model,seed_biomas_model_no_block)

fin_dataSG_biomass_seed_surv_trt_given_germ %>% group_by(precip) %>% summarise_at("total_biomass", funs(n(),mean,sd,se=sd(.)/sqrt(n())))
## A tibble: 2 x 5
#precip     n  mean    sd     se
#<fct>  <int> <dbl> <dbl>  <dbl>
#1 A         27 0.136 0.235 0.0452
#2 D         18 0.122 0.223 0.0526

(0.122-0.136)/0.136
#-0.1029412

fin_dataSG_biomass_seed_surv_trt_given_germ %>% group_by(precip,soil_root) %>% summarise_at("total_biomass", funs(n(),mean,sd,se=sd(.)/sqrt(n())))
(0.0514-0.0863)/0.0863



fin_dataSG_biomass_seed_surv_trt_given_germ_soil_root_precip_g=fin_dataSG_biomass_seed_surv_trt_given_germ %>% group_by(soil_root,precip)
total_bio_soil_root_precip=summarise_at(fin_dataSG_biomass_seed_surv_trt_given_germ_soil_root_precip_g, 
                                       "total_biomass", funs(n(),mean,sd,se=sd(.)/sqrt(n())))

ggplot(total_bio_soil_root_precip, aes(x=precip,y=mean, ymin = mean-se, ymax= mean+se,fill=soil_root))+
  geom_bar(position="dodge",stat="identity", color="black",aes(fill=soil_root))+
  geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+
  scale_fill_manual(values = c("gray", "white", "black"),labels=c("Bulk","Sterile","Rhizo"))+
  scale_x_discrete(labels=c("Ambient","Drought"))+ylab("Total biomass given gemination")+
  geom_text(aes(y=mean+se+0.02, label=n),position=position_dodge(width=0.9))+theme_bw()


#####Seed Bio Reformated bar graph #####

(seed_biomass_p=ggplot(total_bio_soil_root_precip, aes(x=precip,y=mean, ymin = mean-se, ymax= mean+se,fill=factor(soil_root,levels = treatment_order)))+
    geom_bar(position="dodge",stat="identity", color="black",aes(fill=factor(soil_root,levels = treatment_order)))+
    geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+ylim(c(0,.95))+
    scale_fill_manual(values = c( "white","lightgray", "darkgrey"),labels=c("Sterile","Bulk","Rhizo"), limits=treatment_order, name=NULL)+
    scale_x_discrete(labels=c("Ambient","Drought"), name=NULL)+ylab("Total biomass \ngiven gemination (g)")+
    geom_text(aes(y=mean+se+0.05, label=n),size=10,position=position_dodge(width=0.9))+theme_bw()+
    theme(axis.title = element_text(size = 23), axis.text = element_text(size = 23),
          legend.position = c(0.85,.9), legend.text=element_text(size=20),
          legend.background = element_rect(size=0.5,linetype="solid",colour ="black")))

(seed_biomass_p2=ggplot(total_bio_soil_root_precip, aes(x=precip,y=mean, ymin = mean-se, ymax= mean+se,fill=factor(soil_root,levels = treatment_order)))+
    geom_bar(position="dodge",stat="identity", color="black",aes(fill=factor(soil_root,levels = treatment_order)))+
    geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+ylim(c(0,.95))+
    scale_fill_manual(values = c( "white","lightgray", "darkgrey"),labels=c("Sterile","Bulk","Rhizo"), limits=treatment_order, name=NULL)+
    scale_x_discrete(labels=c("Ambient","Drought"), name=NULL)+ylab("Biomass \ngiven gemination (g)")+
    geom_text(aes(y=mean+se+0.05, label=n),size=10,position=position_dodge(width=0.9))+theme_bw()+
    theme(axis.title = element_text(size = 23), axis.text = element_text(size = 23),
          legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

#####Presence total biomass of germinated#####

fin_dataSG_biomass_seed_surv_pres_w_germ=subset(fin_dataSG_biomass_seed_surv_trt_given_germ, root_association =="B")
nrow(fin_dataSG_biomass_seed_surv_pres_w_germ)
#26
unique(fin_dataSG_biomass_seed_surv_pres_w_germ$soil_root)

seed_biomas_model_pres= lm((total_biomass)~precip*soil_root+(1|block), data= fin_dataSG_biomass_seed_surv_pres_w_germ)
qqPlot(resid(seed_biomas_model_pres))
hist(resid(seed_biomas_model_pres))
boxCox(seed_biomas_model_pres)
shapiro.test(resid(seed_biomas_model_pres))
#p-value = 0.1883
anova(seed_biomas_model_pres, type=3)
#precip           0.36711  1 112.9819 3.464e-09 ***
#soil_root        0.37063  1 114.0662 3.214e-09 ***
#precip:soil_root 0.44304  1 136.3510 7.821e-10 ***

emmeans(seed_biomas_model_pres, ~precip)
emmeans(seed_biomas_model_pres, pairwise~soil_root)
emmeans(seed_biomas_model_pres, pairwise~soil_root|precip)


emmeans(seed_biomas_model_pres, pairwise~soil_root*precip, adjust="fdr")

"$contrasts
 contrast      estimate    SE df t.ratio p.value
 L.B,A - S.B,A  0.043308 0.0325 18   1.332 0.2395 
 L.B,A - L.B,D  0.043786 0.0271 18   1.618 0.1846 
 L.B,A - S.B,D -0.617603 0.0482 18 -12.826 <.0001 
 S.B,A - L.B,D  0.000478 0.0323 18   0.015 0.9884 
 S.B,A - S.B,D -0.660912 0.0527 18 -12.548 <.0001 
 L.B,D - S.B,D -0.661389 0.0494 18 -13.395 <.0001 
 
 Results are averaged over the levels of: block 
Note: contrasts are still on the ( scale 
P value adjustment: fdr method for 6 tests 
"
seed_biomas_model_pres_no_block= lm((total_biomass)~precip*soil_root, data= fin_dataSG_biomass_seed_surv_pres_w_germ)
qqPlot(resid(seed_biomas_model_pres_no_block))
hist(resid(seed_biomas_model_pres_no_block))
boxCox(seed_biomas_model_pres_no_block)
shapiro.test(resid(seed_biomas_model_pres_no_block))
#p-value = 0.1224
anova(seed_biomas_model_pres_no_block, type=3)
#precip           0.43448  1  147.58 3.139e-11 ***
#soil_root        0.44799  1  152.17 2.336e-11 ***
#precip:soil_root 0.55107  1  187.19 3.080e-12 ***

emmeans(seed_biomas_model_pres_no_block, pairwise~soil_root*precip,adjust="none")

AIC(seed_biomas_model_pres,seed_biomas_model_pres_no_block)


#####Origin total biomass of germinated#####

fin_dataSG_biomass_seed_surv_orig_w_germ=subset(fin_dataSG_biomass_seed_surv_trt_given_germ, soil_status =="L")
nrow(fin_dataSG_biomass_seed_surv_orig_w_germ)
#38
unique(fin_dataSG_biomass_seed_surv_orig_w_germ$soil_root)

seed_biomas_model_orig= lm(log(total_biomass)~precip*soil_root+(1|block), data= fin_dataSG_biomass_seed_surv_orig_w_germ)
qqPlot(resid(seed_biomas_model_orig))
hist(resid(seed_biomas_model_orig))
boxCox(seed_biomas_model_orig)
shapiro.test(resid(seed_biomas_model_orig))
#p-value = 0.2373
anova(seed_biomas_model_orig, type=3)
#precip             5.795  1   5.6648 0.02386 *  

emmeans(seed_biomas_model_orig, ~precip)
emmeans(seed_biomas_model_orig, pairwise~soil_root)
emmeans(seed_biomas_model_orig, pairwise~soil_root|precip)


emmeans(seed_biomas_model_orig, pairwise~soil_root*precip, adjust="fdr")

"$contrasts
 contrast      estimate    SE df t.ratio p.value
 L.B,A - L.R,A   -0.421 0.437 30 -0.963  0.4122 
 L.B,A - L.B,D    0.685 0.474 30  1.443  0.3187 
 L.B,A - L.R,D    0.554 0.524 30  1.057  0.4122 
 L.R,A - L.B,D    1.105 0.448 30  2.467  0.1173 
 L.R,A - L.R,D    0.975 0.495 30  1.969  0.1749 
 L.B,D - L.R,D   -0.131 0.518 30 -0.252  0.8028 

Results are averaged over the levels of: block 
Results are given on the log (not the response) scale. 
P value adjustment: fdr method for 6 tests  
"
seed_biomas_model_orig_no_block= lm(log(total_biomass)~precip*soil_root, data= fin_dataSG_biomass_seed_surv_orig_w_germ)
qqPlot(resid(seed_biomas_model_orig_no_block))
hist(resid(seed_biomas_model_orig_no_block))
boxCox(seed_biomas_model_orig_no_block)
shapiro.test(resid(seed_biomas_model_orig_no_block))
#p-value = 0.5891
anova(seed_biomas_model_orig_no_block, type=3)
#precip             8.416  1   8.7072  0.005705 ** 

emmeans(seed_biomas_model_orig_no_block, pairwise~soil_root*precip,adjust="none")

AIC(seed_biomas_model_orig,seed_biomas_model_orig_no_block)

#How does the time to germination affect the biomass
#Time to first germination in a given pot 
dataSG_seed_surv_trt_1=subset(dataSG_seed_surv_trt, pot_w_germ==1)

dataSG_seed_surv_trt_1_pot_g=group_by(dataSG_seed_surv_trt_1, Plant_Number)

dataSG_seed_surv_trt_1_pot_first_germ=summarise_at(dataSG_seed_surv_trt_1_pot_g, "exp_days", min)


data_SG_biomass_seed_surv_to_germ_trt=merge(fin_dataSG_biomass_seed_surv_trt_given_germ, dataSG_seed_surv_trt_1_pot_first_germ, by="Plant_Number")
summary(data_SG_biomass_seed_surv_to_germ_trt)
days_to_germ_seed_biomas_model= lm(log(total_biomass)~exp_days.y, data= data_SG_biomass_seed_surv_to_germ_trt)
qqPlot(resid(days_to_germ_seed_biomas_model))
hist(resid(days_to_germ_seed_biomas_model))

anova(days_to_germ_seed_biomas_model, type=3)
#exp_days.y   2.663  1   2.490    0.1221 


ggplot(data_SG_biomass_seed_surv_to_germ_trt, aes(x=exp_days.y,y=log(total_biomass)))+geom_point(aes(fill=soil_root, shape=precip), size=4)+
  theme_bw()+scale_fill_manual(values = c("gray", "white", "black"),labels=c("Bulk","Sterile","Rhizo"))+scale_shape_manual(values=c(21,25))+
  geom_text(aes(y=log(total_biomass)+.1, x= exp_days.y,label=Plant_Number))

#Roots alone
summary(fin_dataSG_biomass_seed_surv_trt_given_germ)
seed_root_model= lm(log(root_weight_g)~precip*soil_root+(1|block), data= fin_dataSG_biomass_seed_surv_trt_given_germ)
qqPlot(resid(seed_root_model))
hist(resid(seed_root_model))

anova(seed_root_model, type=3)
#precip:soil_root  10.773  2   6.1438  0.005275 ** 
#soil_root          5.526  2   3.1512  0.055528 .  

emmeans(seed_root_model, pairwise~soil_root)
emmeans(seed_root_model, pairwise~soil_root|precip)
#precip = A:
#contrast    estimate        SE df t.ratio p.value
#L.B - S.B  0.5289687 0.4822125 34   1.097  0.5225
#L.B - L.R -0.3642337 0.3717826 34  -0.980  0.5945
#S.B - L.R -0.8932024 0.4622551 34  -1.932  0.1451

#precip = D:
#  contrast    estimate        SE df t.ratio p.value
#L.B - S.B -2.6455949 0.7222441 34  -3.663  0.0024
#L.B - L.R -0.5686425 0.4606463 34  -1.234  0.4415
#S.B - L.R  2.0769524 0.7786334 34   2.667  0.0304



fin_dataSG_biomass_seed_surv_trt_given_germ_soil_root_precip_g=fin_dataSG_biomass_seed_surv_trt_given_germ %>% group_by(soil_root,precip)
root_soil_root_precip=summarise_at(fin_dataSG_biomass_seed_surv_trt_given_germ_soil_root_precip_g, 
                                        "root_weight_g", funs(n(),mean,sd,se=sd(.)/sqrt(n())))

ggplot(root_soil_root_precip, aes(x=precip,y=mean, ymin = mean-se, ymax= mean+se,fill=soil_root))+
  geom_bar(position="dodge",stat="identity", color="black",aes(fill=soil_root))+
  geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+
  scale_fill_manual(values = c("gray", "white", "black"),labels=c("Bulk","Sterile","Rhizo"))+
  scale_x_discrete(labels=c("Ambient","Drought"))+ylab("Root biomass given gemination")+
  geom_text(aes(y=mean+se+0.02, label=n),position=position_dodge(width=0.9))+theme_bw()

#now root:shoot ratio
fin_dataSG_biomass_seed_surv_trt_given_germ$root_shoot=fin_dataSG_biomass_seed_surv_trt_given_germ$root_weight_g/fin_dataSG_biomass_seed_surv_trt_given_germ$shoot_weight_g
summary(fin_dataSG_biomass_seed_surv_trt_given_germ)
root_shoot_seed_model= lm(log(root_shoot)~precip*soil_root+(1|block), data= fin_dataSG_biomass_seed_surv_trt_given_germ)
qqPlot(resid(root_shoot_seed_model))
hist(resid(root_shoot_seed_model))
shapiro.test(resid(root_shoot_seed_model))
#0.4825

anova(root_shoot_seed_model, type=3)
#precip             6.246  1   8.2080  0.007008 ** 

fin_dataSG_biomass_seed_surv_trt_given_germ %>% group_by(precip)  %>% summarise_at("root_shoot", ~mean(.))

emmeans(root_shoot_seed_model, pairwise~precip)
emmeans(root_shoot_seed_model, pairwise~soil_root)
emmeans(root_shoot_seed_model, pairwise~soil_root|precip)
emmeans(root_shoot_seed_model, pairwise~soil_root*precip, adjust ="fdr")


root_shoot_seed_model_no_block= lm(log(root_shoot)~precip*soil_root, data= fin_dataSG_biomass_seed_surv_trt_given_germ)
qqPlot(resid(root_shoot_seed_model_no_block))
hist(resid(root_shoot_seed_model_no_block))
shapiro.test(resid(root_shoot_seed_model_no_block))
#0.2474

anova(root_shoot_seed_model_no_block, type=3)
#precip             6.323  1   7.8877  0.007737 ** 
#precip:soil_root   4.828  2   3.0113  0.060792 .

emmeans(root_shoot_seed_model_no_block, pairwise~soil_root|precip)

AIC(root_shoot_seed_model,root_shoot_seed_model_no_block)

fin_dataSG_biomass_seed_surv_trt_given_germ %>% group_by(precip) %>% summarise_at("root_shoot", funs(n(),mean,sd,se=sd(.)/sqrt(n())))
(18.3-11.6)/11.6

fin_dataSG_biomass_seed_surv_trt_given_germ %>% group_by(precip,soil_root) %>% summarise_at("root_shoot", funs(n(),mean,sd,se=sd(.)/sqrt(n())))
#drought presence
(4.88-10.4)/10.4

fin_dataSG_biomass_seed_surv_trt_given_germ_soil_root_precip_g=fin_dataSG_biomass_seed_surv_trt_given_germ %>% group_by(soil_root,precip)
root_shoot_soil_root_precip=summarise_at(fin_dataSG_biomass_seed_surv_trt_given_germ_soil_root_precip_g, 
                                   "root_shoot", funs(n(),mean,sd,se=sd(.)/sqrt(n())))

ggplot(root_shoot_soil_root_precip, aes(x=precip,y=mean, ymin = mean-se, ymax= mean+se,fill=soil_root))+
  geom_bar(position="dodge",stat="identity", color="black",aes(fill=soil_root))+
  geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+
  scale_fill_manual(values = c("gray", "white", "black"),labels=c("Bulk","Sterile","Rhizo"))+
  scale_x_discrete(labels=c("Ambient","Drought"))+ylab("Root:shoot given gemination")+
  geom_text(aes(y=mean+se+1, label=n),position=position_dodge(width=0.9))+theme_bw()


#####Seed root:shoot reformated graph####

(seed_root_t_shoot_p=ggplot(root_shoot_soil_root_precip, aes(x=precip,y=mean, ymin = mean-se, ymax= mean+se,fill=factor(soil_root,levels = treatment_order)))+
    geom_bar(position="dodge",stat="identity", color="black",aes(fill=factor(soil_root,levels = treatment_order)))+
    geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+ylim(c(0,45))+
    scale_fill_manual(values = c( "white","light gray", "dark grey"),labels=c("Sterile","Bulk","Rhizo"), limits=treatment_order, name=NULL)+
    scale_x_discrete(labels=c("Ambient","Drought"), name=NULL)+ylab("Root:shoot given gemination")+
    geom_text(aes(y=2, label=n),size=10,position=position_dodge(width=0.9))+theme_bw()+
    theme(axis.title = element_text(size = 23), axis.text = element_text(size = 23),
          legend.position = c(0.85,.9), legend.text=element_text(size=20),
          legend.background = element_rect(size=0.5,linetype="solid",colour ="black")))

(seed_root_t_shoot_p=ggplot(root_shoot_soil_root_precip, aes(x=precip,y=mean, ymin = mean-se, ymax= mean+se,fill=factor(soil_root,levels = treatment_order)))+
    geom_bar(position="dodge",stat="identity", color="black",aes(fill=factor(soil_root,levels = treatment_order)))+
    geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+ylim(c(0,50))+
    scale_fill_manual(values = c( "white","lightgray", "darkgrey"),labels=c("Sterile","Bulk","Rhizo"), limits=treatment_order, name=NULL)+
    scale_x_discrete(labels=c("Ambient","Drought"), name=NULL)+ylab("Root:shoot given gemination")+
    geom_text(aes(y=2, label=n),size=10,position=position_dodge(width=0.9))+theme_bw()+
    theme(axis.title = element_text(size = 23), axis.text = element_text(size = 23),
          legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

#####Presence root:shoot ratio#####
fin_dataSG_biomass_seed_surv_pres_w_germ=subset(fin_dataSG_biomass_seed_surv_trt_given_germ, root_association =="B")

summary(fin_dataSG_biomass_seed_surv_pres_w_germ)
root_shoot_seed_model_pres= lm(log(root_shoot)~precip*soil_root+(1|block), data= fin_dataSG_biomass_seed_surv_pres_w_germ)
qqPlot(resid(root_shoot_seed_model_pres))
hist(resid(root_shoot_seed_model_pres))
shapiro.test(resid(root_shoot_seed_model_pres))
#0.5222

anova(root_shoot_seed_model_pres, type=3)
#precip            4.595  1   7.4445   0.01379 * 
#precip:soil_root  3.299  1   5.3456   0.03282 * 

fin_dataSG_biomass_seed_surv_trt_given_germ %>% group_by(precip)  %>% summarise_at("root_shoot", ~mean(.))

emmeans(root_shoot_seed_model_pres, pairwise~precip)
emmeans(root_shoot_seed_model_pres, pairwise~soil_root)
emmeans(root_shoot_seed_model_pres, pairwise~soil_root|precip)
emmeans(root_shoot_seed_model_pres, pairwise~soil_root*precip, adjust ="fdr")


root_shoot_seed_model_pres_no_block= lm(log(root_shoot)~precip*soil_root, data= fin_dataSG_biomass_seed_surv_pres_w_germ)
qqPlot(resid(root_shoot_seed_model_pres_no_block))
hist(resid(root_shoot_seed_model_pres_no_block))
shapiro.test(resid(root_shoot_seed_model_pres_no_block))
#0.6602

anova(root_shoot_seed_model_pres_no_block, type=3)
#precip            5.277  1   9.7046  0.005042 ** 
#precip:soil_root  4.823  1   8.8704  0.006934 ** 

emmeans(root_shoot_seed_model_pres_no_block, pairwise~soil_root*precip,adjust="none")

AIC(root_shoot_seed_model_pres,root_shoot_seed_model_pres_no_block)


#####Origin root:shoot ratio#####
fin_dataSG_biomass_seed_surv_orig_w_germ=subset(fin_dataSG_biomass_seed_surv_trt_given_germ, soil_status =="L")

summary(fin_dataSG_biomass_seed_surv_orig_w_germ)
root_shoot_seed_model_orig= lm(log(root_shoot)~precip*soil_root+(1|block), data= fin_dataSG_biomass_seed_surv_orig_w_germ)
qqPlot(resid(root_shoot_seed_model_orig))
hist(resid(root_shoot_seed_model_orig))
shapiro.test(resid(root_shoot_seed_model_orig))
#0.5222

anova(root_shoot_seed_model_orig, type=3)
#nada sig

fin_dataSG_biomass_seed_surv_trt_given_germ %>% group_by(precip)  %>% summarise_at("root_shoot", ~mean(.))

emmeans(root_shoot_seed_model_orig, pairwise~precip)
emmeans(root_shoot_seed_model_orig, pairwise~soil_root)
emmeans(root_shoot_seed_model_orig, pairwise~soil_root|precip)
emmeans(root_shoot_seed_model_orig, pairwise~soil_root*precip, adjust ="fdr")


root_shoot_seed_model_orig_no_block= lm(log(root_shoot)~precip*soil_root, data= fin_dataSG_biomass_seed_surv_orig_w_germ)
qqPlot(resid(root_shoot_seed_model_orig_no_block))
hist(resid(root_shoot_seed_model_orig_no_block))
shapiro.test(resid(root_shoot_seed_model_orig_no_block))
#0.2467

anova(root_shoot_seed_model_orig_no_block, type=3)
#nada sig

emmeans(root_shoot_seed_model_orig_no_block, pairwise~soil_root*precip,adjust="none")

AIC(root_shoot_seed_model_orig,root_shoot_seed_model_orig_no_block)





#let's remove the sandy roots

#now lets look at total biomass of germinated transs
fin_dataSG_biomass_seed_surv_trt_given_germ=subset(data_SG_biomass_seed_surv_trt, germinated>0)
fin_dataSG_biomass_seed_surv_trt_given_germ_N_sandy=subset(fin_dataSG_biomass_seed_surv_trt_given_germ, sandy_root=="N")
summary(fin_dataSG_biomass_seed_surv_trt_given_germ_N_sandy)

seed_biomas_model_N_sandy= lm(log(total_biomass)~precip*soil_root+(1|block), data= fin_dataSG_biomass_seed_surv_trt_given_germ_N_sandy)
qqPlot(resid(seed_biomas_model_N_sandy))
hist(resid(seed_biomas_model_N_sandy))

anova(seed_biomas_model_N_sandy, type=3)
#precip             1.470  1   2.9040  0.098696 .  
#soil_root          2.761  2   2.7274  0.081602 .  
#precip:soil_root   7.339  2   7.2505  0.002699 ** 

emmeans(seed_biomas_model_N_sandy, ~precip)
emmeans(seed_biomas_model_N_sandy, pairwise~soil_root|precip)
#$contrasts
#precip = A:
#  contrast  estimate    SE df t.ratio p.value
#L.B - S.B   0.7512 0.441 30  1.704  0.2203 
#L.B - L.R   0.0864 0.323 30  0.267  0.9614 
#S.B - L.R  -0.6647 0.429 30 -1.549  0.2830 

#precip = D:
#  contrast  estimate    SE df t.ratio p.value
#L.B - S.B  -2.7719 0.788 30 -3.516  0.0039 
#L.B - L.R  -0.5114 0.382 30 -1.340  0.3849 
#S.B - L.R   2.2606 0.829 30  2.727  0.0277 



fin_dataSG_biomass_seed_surv_trt_N_sandy_w_germ_soil_root_precip_g=fin_dataSG_biomass_seed_surv_trt_given_germ_N_sandy %>% group_by(soil_root,precip)
total_bio_N_sandy_soil_root_precip=summarise_at(fin_dataSG_biomass_seed_surv_trt_N_sandy_w_germ_soil_root_precip_g, 
                                        "total_biomass", funs(n(),mean,sd,se=sd(.)/sqrt(n())))

ggplot(total_bio_N_sandy_soil_root_precip, aes(x=precip,y=mean, ymin = mean-se, ymax= mean+se,fill=soil_root))+
  geom_bar(position="dodge",stat="identity", color="black",aes(fill=soil_root))+
  geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+
  scale_fill_manual(values = c("gray", "white", "black"),labels=c("Bulk","Sterile","Rhizo"))+
  scale_x_discrete(labels=c("Ambient","Drought"))+ylab("Total biomass given gemination (with sandy samples removed)")+
  geom_text(aes(y=mean+0.04, label=n),position=position_dodge(width=0.9))+theme_bw()



#Roots alone
summary(fin_dataSG_biomass_seed_surv_trt_given_germ_N_sandy)
seed_root_model_N_sandy= lm(log(root_weight_g)~precip*soil_root+(1|block), data= fin_dataSG_biomass_seed_surv_trt_given_germ_N_sandy)
qqPlot(resid(seed_root_model_N_sandy))
hist(resid(seed_root_model_N_sandy))

anova(seed_root_model_N_sandy, type=3)
#precip:soil_root   8.064  2   6.4893   0.00455 ** 
#precip             2.007  1   3.2309   0.08233 .   

emmeans(seed_root_model_N_sandy, pairwise~soil_root)
emmeans(seed_root_model_N_sandy, pairwise~soil_root|precip)
#precip = A:
#contrast  estimate    SE df t.ratio p.value
#L.B - S.B    0.811 0.488 30  1.660  0.2369 
#L.B - L.R    0.105 0.358 30  0.292  0.9542 
#S.B - L.R   -0.706 0.475 30 -1.486  0.3118 

#precip = D:
#  contrast  estimate    SE df t.ratio p.value
#L.B - S.B   -2.873 0.873 30 -3.289  0.0070 
#L.B - L.R   -0.563 0.423 30 -1.331  0.3897 
#S.B - L.R    2.310 0.919 30  2.515  0.0448 



fin_dataSG_biomass_seed_surv_trt_N_sandy_w_germ_soil_root_precip_g=fin_dataSG_biomass_seed_surv_trt_given_germ_N_sandy %>% group_by(soil_root,precip)
root_N_sandy_soil_root_precip=summarise_at(fin_dataSG_biomass_seed_surv_trt_N_sandy_w_germ_soil_root_precip_g, 
                                   "root_weight_g", funs(n(),mean,sd,se=sd(.)/sqrt(n())))

ggplot(root_N_sandy_soil_root_precip, aes(x=precip,y=mean, ymin = mean-se, ymax= mean+se,fill=soil_root))+
  geom_bar(position="dodge",stat="identity", color="black",aes(fill=soil_root))+
  geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+
  scale_fill_manual(values = c("gray", "white", "black"),labels=c("Bulk","Sterile","Rhizo"))+
  scale_x_discrete(labels=c("Ambient","Drought"))+ylab("Root biomass given gemination (with sandy samples removed)")+
  geom_text(aes(y=mean+0.04, label=n),position=position_dodge(width=0.9))+theme_bw()




#Shoots alone
summary(fin_dataSG_biomass_seed_surv_trt_given_germ)
seed_shoot_model= lm(log(shoot_weight_g)~precip*soil_root+(1|block), data= fin_dataSG_biomass_seed_surv_trt_given_germ)
qqPlot(resid(seed_shoot_model))
hist(resid(seed_shoot_model))

anova(seed_shoot_model, type=3)
#precip             2.06  1   10.0215  0.003257 ** 
#soil_root          1.63  2    3.9621  0.028395 *  
#(1|block)   4.70  4    5.7232  0.001240 ** 
#precip:soil_root   3.28  2    7.9750  0.001446 ** 

emmeans(seed_shoot_model, ~precip)
emmeans(seed_shoot_model, pairwise~soil_root)
emmeans(seed_shoot_model, pairwise~soil_root|precip)
#precip = A:
# L.B - S.B  0.14946329 0.2537568 34   0.589  0.8269
#L.B - L.R -0.25671123 0.1956448 34  -1.312  0.3983
#S.B - L.R -0.40617452 0.2432545 34  -1.670  0.2314

#precip = D:
#  contrast     estimate        SE df t.ratio p.value
#L.B - S.B -1.41470888 0.3800696 34  -3.722  0.0020
#L.B - L.R  0.08175893 0.2424079 34   0.337  0.9393
#S.B - L.R  1.49646782 0.4097436 34   3.652  0.0024



fin_dataSG_biomass_seed_surv_trt_given_germ_soil_root_precip_g=fin_dataSG_biomass_seed_surv_trt_given_germ %>% group_by(soil_root,precip)
shoot_soil_root_precip=summarise_at(fin_dataSG_biomass_seed_surv_trt_given_germ_soil_root_precip_g, 
                                   "shoot_weight_g", funs(n(),mean,sd,se=sd(.)/sqrt(n())))

ggplot(shoot_soil_root_precip, aes(x=precip,y=mean, ymin = mean-se, ymax= mean+se,fill=soil_root))+
  geom_bar(position="dodge",stat="identity", color="black",aes(fill=soil_root))+
  geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+
  scale_fill_manual(values = c("gray", "white", "black"),labels=c("Bulk","Sterile","Rhizo"))+
  scale_x_discrete(labels=c("Ambient","Drought"))+ylab("Shoot biomass given gemination")+
  geom_text(aes(y=mean+se+0.001, label=n),position=position_dodge(width=0.9))+theme_bw()

#pot weight of seed germination pots

dataSG_potweight_trt_seed=subset(dataSG_potweight_trt, life_stage=="S")
nrow(dataSG_potweight_trt_seed)
#450
dataSG_potweight_trt_seed$soil_root=with(dataSG_potweight_trt_seed, interaction(soil_status,root_association))
dataSG_potweight_trt_seed_sub=subset(dataSG_potweight_trt_seed, date!="6_13")
summary(dataSG_potweight_trt_seed_sub)

min(dataSG_potweight_trt_seed_sub$per_ch_potweight)
delta_pot_weight_seed_model= lmer((per_ch_potweight+0.3913044)^2~precip*soil_root*date+(1|block)+(1|Plant_Number), data= dataSG_potweight_trt_seed_sub)
plot(delta_pot_weight_seed_model)
qqPlot(resid(delta_pot_weight_seed_model))
hist(resid(delta_pot_weight_seed_model))
shapiro.test(resid(delta_pot_weight_seed_model))

anova(delta_pot_weight_seed_model, type=3)
#date                  190.6915  3  < 2.2e-16 ***
#precip:date            39.7663  3  1.194e-08 ***
#soil_root:date         14.7113  6    0.02262 *  
emmeans(delta_pot_weight_seed_model, ~precip)
emmeans(delta_pot_weight_seed_model, pairwise~soil_root)

emmeans(delta_pot_weight_seed_model, pairwise~soil_root|date|precip)
emmeans(delta_pot_weight_seed_model, pairwise~soil_root|date)
emmeans(delta_pot_weight_seed_model, pairwise~precip|date)

seedSG_delta_potweight_trt_soil_root_time_g=dataSG_potweight_trt_seed_sub %>% group_by(date,soil_root)
pot_delta_potweight_seed_soil_root_time=summarise_at(seedSG_delta_potweight_trt_soil_root_time_g, "per_ch_potweight", mean)
#dot graph
ggplot(pot_delta_potweight_seed_soil_root_time, aes(x=date,y=per_ch_potweight))+geom_point(size=4,aes(color=soil_root))+
  scale_color_manual(values = c("gray", "white", "black"),labels=c("Bulk","Sterile","Rhizo"))+ylab("percentage change in pot weight from int")+
  xlab("date")
#boxplot
ggplot(dataSG_potweight_trt_seed_sub, aes(x=date,y=per_ch_potweight))+geom_boxplot(aes(fill=soil_root))+
  scale_fill_manual(values = c("gray", "white", "black"),labels=c("Bulk","Sterile","Rhizo"))+ylab("percentage change in pot weight from int")+
  xlab("date")
#dot graph
seedSG_delta_potweight_trt_soil_root_time_drought_g=dataSG_potweight_trt_seed_sub %>% group_by(date,soil_root,precip)
pot_delta_potweight_seed_soil_root_time_drought=summarise_at(seedSG_delta_potweight_trt_soil_root_time_drought_g, "per_ch_potweight", mean)

ggplot(pot_delta_potweight_seed_soil_root_time_drought, aes(x=date,y=per_ch_potweight))+geom_point(shape=21,size=4,aes(color=precip,fill=soil_root))+
  ylab("percentage change in pot weight from int")+scale_fill_manual(values = c("gray", "white", "black"),labels=c("Bulk","Sterile","Rhizo"))+
  scale_color_manual(values = c("blue", "red"),labels=c("Ambient","Drought"))+ xlab("Date")

#boxplot
ggplot(dataSG_potweight_trt_seed_sub, aes(x=date,y=per_ch_potweight))+geom_boxplot(aes(color=precip,fill=soil_root))+
  ylab("percentage change in pot weight from int")+scale_fill_manual(values = c("gray", "white", "black"),labels=c("Bulk","Sterile","Rhizo"))+
  scale_color_manual(values = c("blue", "red"),labels=c("Ambient","Drought"))+ xlab("Date")


#Extremely large outliers in bulk treatment
rank(dataSG_potweight_trt_seed_sub$per_ch_potweight)
#I am going to remove pot 20 and 71 since their intial weight is suspicious and the very lowest values

dataSG_potweight_trt_seed_sub_out=subset(dataSG_potweight_trt_seed_sub, per_ch_potweight<.50)
dataSG_potweight_trt_seed_sub_out=subset(dataSG_potweight_trt_seed_sub_out, per_ch_potweight>-.3)
nrow(dataSG_potweight_trt_seed_sub_out)
#351


#####Diagnostic plots for probablity disturbutions####

#https://ase.tufts.edu/gsc/gradresources/guidetomixedmodelsinr/mixed%20model%20guide.html
require(MASS)

#I want no zeros
dataSG_potweight_trt_seed_sub_out$per_ch_potweight_bio_crt.t <- dataSG_potweight_trt_seed_sub_out$per_ch_potweight_bio_crt + 0.178
min(dataSG_potweight_trt_seed_sub_out$per_ch_potweight_bio_crt.t)

# lnorm means lognormal
qqp(dataSG_potweight_trt_seed_sub_out$per_ch_potweight_bio_crt.t, "lnorm")

# qqp requires estimates of the parameters of the negative binomial, Poisson
# and gamma distributions. You can generate estimates using the fitdistr
# function. Save the output and extract the estimates of each parameter as I
# have shown below.
nbinom <- fitdistr(dataSG_potweight_trt_seed_sub_out$per_ch_potweight_bio_crt.t, "Negative Binomial")
qqp(dataSG_potweight_trt_seed_sub_out$per_ch_potweight.t, "nbinom", size = nbinom$estimate[[1]], mu = nbinom$estimate[[2]])

poisson <- fitdistr(dataSG_potweight_trt_seed_sub_out$per_ch_potweight_bio_crt.t, "Poisson")
qqp(dataSG_potweight_trt_seed_sub_out$per_ch_potweight.t, "pois", lambda=poisson$estimate)

gamma <- fitdistr(dataSG_potweight_trt_seed_sub_out$per_ch_potweight_bio_crt.t, "gamma")
qqp(dataSG_potweight_trt_seed_sub_out$per_ch_potweight_bio_crt.t, "gamma", shape = gamma$estimate[[1]], rate = gamma$estimate[[2]])

colnames(dataSG_potweight_trt_seed_sub_out)
min(dataSG_potweight_trt_seed_sub_out$per_ch_potweight_bio_crt)
delta_pot_weight_seed__out_model= lmer(sqrt(per_ch_potweight_bio_crt+0.178)~precip*soil_root*date+(1|block)+(1|Plant_Number), data= dataSG_potweight_trt_seed_sub_out)
plot(delta_pot_weight_seed__out_model)
qqPlot(resid(delta_pot_weight_seed__out_model))
hist(resid(delta_pot_weight_seed__out_model))
shapiro.test(resid(delta_pot_weight_seed__out_model))
# 2.574e-10



anova(delta_pot_weight_seed__out_model, type=3)
#precip                  69.5599  1  < 2.2e-16 ***
#soil_root               75.4637  2  < 2.2e-16 ***
#date                   821.8939  3  < 2.2e-16 ***
#precip:date             41.4204  3  5.325e-09 ***
#soil_root:date          65.7602  6  3.017e-12 ***
emmeans(delta_pot_weight_seed__out_model, ~precip)
emmeans(delta_pot_weight_seed__out_model, pairwise~soil_root)

emmeans(delta_pot_weight_seed__out_model, pairwise~soil_root|date|precip)
emmeans(delta_pot_weight_seed__out_model, pairwise~soil_root|date)
emmeans(delta_pot_weight_seed__out_model, pairwise~precip|date)

seedSG_delta_potweight_trt_soil_root_time_g=dataSG_potweight_trt_seed_sub_out %>% group_by(date,soil_root)
pot_delta_potweight_seed_soil_root_time=summarise_at(seedSG_delta_potweight_trt_soil_root_time_g, "per_ch_potweight", mean)
#dot graph
ggplot(pot_delta_potweight_seed_soil_root_time, aes(x=date,y=per_ch_potweight))+geom_point(size=4,aes(color=soil_root))+
  scale_color_manual(values = c("gray", "white", "black"),labels=c("Bulk","Sterile","Rhizo"))+ylab("percentage change in pot weight from int")+
  xlab("date")
#boxplot
ggplot(dataSG_potweight_trt_seed_sub_out, aes(x=date,y=per_ch_potweight))+geom_boxplot(aes(fill=soil_root))+
  scale_fill_manual(values = c("gray", "white", "black"),labels=c("Bulk","Sterile","Rhizo"))+ylab("percentage change in pot weight from int")+
  xlab("date")
#dot graph
seedSG_delta_potweight_trt_soil_root_time_drought_g=dataSG_potweight_trt_seed_sub_out %>% group_by(date,soil_root,precip)
pot_delta_potweight_seed_soil_root_time_drought=summarise_at(seedSG_delta_potweight_trt_soil_root_time_drought_g, "per_ch_potweight", mean)

ggplot(pot_delta_potweight_seed_soil_root_time_drought, aes(x=date,y=per_ch_potweight))+geom_point(shape=21,size=4,aes(color=precip,fill=soil_root))+
  ylab("percentage change in pot weight from int")+scale_fill_manual(values = c("gray", "white", "black"),labels=c("Bulk","Sterile","Rhizo"))+
  scale_color_manual(values = c("blue", "red"),labels=c("Ambient","Drought"))+ xlab("Date")

#boxplot
ggplot(dataSG_potweight_trt_seed_sub_out, aes(x=date,y=per_ch_potweight))+geom_boxplot(aes(color=precip,fill=soil_root))+
  ylab("percentage change in pot weight from int")+scale_fill_manual(values = c("gray", "white", "black"),labels=c("Bulk","Sterile","Rhizo"))+
  scale_color_manual(values = c("blue", "red"),labels=c("Ambient","Drought"))+ xlab("Date")

#let's look at the final VWC alone
unique(dataSG_potweight_trt_seed_sub_out$date)
dataSG_potweight_trt_seed_sub_out_end=subset(dataSG_potweight_trt_seed_sub_out, date=="7_7")

delta_pot_weight_seed_end_out_model= lm((per_ch_potweight_bio_crt)~precip*soil_root+(1|block), data= dataSG_potweight_trt_seed_sub_out_end)
#plot(delta_pot_weight_seed_end_out_model)
qqPlot(resid(delta_pot_weight_seed_end_out_model))
hist(resid(delta_pot_weight_seed_end_out_model))
shapiro.test(resid(delta_pot_weight_seed_end_out_model))
# 0.8819



anova(delta_pot_weight_seed_end_out_model, type=3)
#precip           0.02520  1   29.5266 6.087e-07 ***
#soil_root        0.07180  2   42.0700 4.036e-13 ***




#####VWC as a covariate#####

dataSG_potweight_trt_seed=subset(dataSG_potweight_trt, life_stage=="S")
nrow(dataSG_potweight_trt_seed)
summary(dataSG_potweight_trt_seed)
dataSG_potweight_trt_seed$soil_root=with(dataSG_potweight_trt_seed, interaction(soil_status,root_association))
dataSG_potweight_trt_seed_sub=subset(dataSG_potweight_trt_seed, date!="6_13")
summary(dataSG_potweight_trt_seed_sub)
#Extremely large outliers in bulk treatment

#I am going to remove pot 20 and 71 since their intial weight is suspicious and the very lowest values

dataSG_potweight_trt_seed_sub_out=subset(dataSG_potweight_trt_seed_sub, per_ch_potweight<.50)
dataSG_potweight_trt_seed_sub_out=subset(dataSG_potweight_trt_seed_sub_out, per_ch_potweight>-.3)
nrow(dataSG_potweight_trt_seed_sub_out)
#351

#I am going to use the end chnage in pot weight
dataSG_potweight_trt_seed_sub_out_end=subset(dataSG_potweight_trt_seed_sub_out,date=="7_7")
nrow(dataSG_potweight_trt_seed_sub_out_end)
#88
dataSG_potweight_end=dataSG_potweight_trt_seed_sub_out_end[,c("Plant_Number","per_ch_potweight","per_ch_potweight_bio_crt")]
colnames(dataSG_potweight_end)[2]="end_per_ch_potweight"
colnames(dataSG_potweight_end)[3]="end_per_ch_potweight_bio_crt"
summary(dataSG_potweight_trans_end)

#let's look at germination and number of germinates

#let's look at prop of pots with germinates through time
summary(dataSG_potweight_end)
summary(dataSG_seed_surv_trt)


dataSG_seed_surv_trt_pw=merge(dataSG_seed_surv_trt,dataSG_potweight_end, by="Plant_Number", all.y = T)
summary(dataSG_seed_surv_trt_pw)
nrow(dataSG_seed_surv_trt_pw)
#1584
dataSG_seed_surv_trt_pw$pot_w_germ=dataSG_seed_surv_trt_pw$num_germinates
dataSG_seed_surv_trt_pw$pot_w_germ[dataSG_seed_surv_trt_pw$pot_w_germ>0]=1
summary(dataSG_seed_surv_trt_pw)
hist(dataSG_seed_surv_trt_pw$pot_w_germ)
dataSG_seed_surv_trt_pw$soil_root=with(dataSG_seed_surv_trt_pw, interaction(soil_status,root_association))

pot_germ_VWC_model= glmer(pot_w_germ~end_per_ch_potweight*soil_root*exp_days+(1|block)+(1|Plant_Number), data= dataSG_seed_surv_trt_pw, family = binomial)
#Warning messages:
#1: In checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,  :
#                  Model failed to converge with max|grad| = 0.254375 (tol = 0.001, component 1)
#                2: In checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,  :
#                                  Model is nearly unidentifiable: large eigenvalue ratio
#                                - Rescale variables?
qqPlot(resid(pot_germ_VWC_model))
hist(resid(pot_germ_VWC_model))
plot(pot_germ_VWC_model)

anova(pot_germ_VWC_model, type=3)
#exp_days                                 5.5701  1   0.018269 * 
#(1|block)                        14.2177  4   0.006632 **


emmeans(pot_germ_VWC_model, pairwise~soil_root)


#UPDATE BELOW
dataSG_seed_surv_trt_soil_root_time_g=dataSG_seed_surv_trt %>% group_by(exp_days,soil_root)
pot_germ_soil_root_time=summarise_at(dataSG_seed_surv_trt_soil_root_time_g, "pot_w_germ", mean)

ggplot(pot_germ_soil_root_time, aes(x=exp_days,y=pot_w_germ))+geom_point(size=4,aes(color=soil_root))+
  scale_color_manual(values = c("gray", "white", "black"),labels=c("Bulk","Sterile","Rhizo"))+ylab("prop Pots with at least one germinate")+
  xlab("Days from start of exp")


dataSG_seed_surv_trt_soil_root_time_drought_g=dataSG_seed_surv_trt %>% group_by(exp_days,soil_root,precip)
pot_germ_soil_root_time_drought=summarise_at(dataSG_seed_surv_trt_soil_root_time_drought_g, "pot_w_germ", mean)

ggplot(pot_germ_soil_root_time_drought, aes(x=exp_days,y=pot_w_germ))+geom_point(shape=21,size=4,aes(color=precip,fill=soil_root))+
  ylab("prop Pots with at least one germinate")+scale_fill_manual(values = c("gray", "white", "black"),labels=c("Bulk","Sterile","Rhizo"))+
  scale_color_manual(values = c("blue", "red"),labels=c("Ambient","Drought"))+ xlab("Days from start of exp")




#let's look at the final number of germinates
#harvested transs plus final survival

dataSG_seed_surv_trt_pw_g=group_by(dataSG_seed_surv_trt_pw, Plant_Number)
total_transs_harvested_pw=summarise_at(dataSG_seed_surv_trt_pw_g, "num_removed", sum)
summary(total_transs_harvested_pw)
nrow(total_transs_harvested_pw)
colnames(total_transs_harvested_pw)[2]="tot_num_removed"

fin_dataSG_seed_surv_trt_pw=subset(dataSG_seed_surv_trt_pw, date=="7_09")

fin_dataSG_seed_surv_trt_pw$tot_num_germ=fin_dataSG_seed_surv_trt_pw$num_germinates+total_transs_harvested_pw$tot_num_removed

fin_dataSG_seed_surv_trt_pw$germinated=fin_dataSG_seed_surv_trt_pw$tot_num_germ
fin_dataSG_seed_surv_trt_pw$germinated[fin_dataSG_seed_surv_trt_pw$germinated>0]=1
fin_dataSG_seed_surv_trt_pw$soil_root=with(fin_dataSG_seed_surv_trt_pw, interaction(soil_status,root_association))

summary(fin_dataSG_seed_surv_trt_pw)


hist(fin_dataSG_seed_surv_trt_pw$tot_num_germ)


#pots with germination

seed_germ_model_pw= glm(germinated~end_per_ch_potweight_bio_crt*soil_root+(1|block), data= fin_dataSG_seed_surv_trt_pw, family = binomial)
qqPlot(resid(seed_germ_model_pw))
hist(resid(seed_germ_model_pw))

anova(seed_germ_model_pw, type=3)
#(1|block)                 9.5179  4    0.04938 *
anova(seed_germ_model_pw, type=2)
aov(seed_germ_model_pw)
emmeans(seed_germ_model_pw, pairwise~soil_root)
emmeans(seed_germ_model_pw, pairwise~soil_root|precip)
"$contrasts
contrast  estimate    SE  df z.ratio p.value
L.B - S.B   1.6140 0.789 Inf  2.045  0.1015 
L.B - L.R   0.0633 0.750 Inf  0.084  0.9961 
S.B - L.R  -1.5508 0.775 Inf -2.002  0.1117  "


#Bars by soil_root

fin_dataSG_seed_surv_trt_soil_root_g=group_by(fin_dataSG_seed_surv_trt, soil_root)
overal_germ_soil_root=summarise_at(fin_dataSG_seed_surv_trt_soil_root_g, "germinated", mean)

ggplot(overal_germ_soil_root, aes(soil_root,germinated))+geom_bar(stat="identity", color="black",aes(fill=soil_root))+
  scale_fill_manual(values = c("gray", "white", "black"),labels=c("Bulk","Sterile","Rhizo"))+
  scale_x_discrete(limits=c("S.B", "L.B", "L.R"),labels=c("Sterile","Bulk","Rhizo"))+ylab("prop Pots with at least one germinate")



fin_dataSG_seed_surv_trt_precip_g=group_by(fin_dataSG_seed_surv_trt, precip)
overal_germ_precip=summarise_at(fin_dataSG_seed_surv_trt_precip_g, "germinated", mean)

ggplot(overal_germ_precip, aes(precip,germinated))+geom_bar(stat="identity", color="black",aes(fill=precip))+
  scale_fill_manual(values = c("gray", "white"),labels=c("Ambient","Drought"))+
  scale_x_discrete(labels=c("Ambient","Drought"))+ylab("prop Pots with at least one germinate")


fin_dataSG_seed_surv_trt_precip_soil_g=fin_dataSG_seed_surv_trt %>% group_by(soil_root,precip)
surv_soil_root_precip=summarise_at(fin_dataSG_seed_surv_trt_precip_soil_g, 
                                   "germinated", funs(n(),mean,sd,se=sd(.)/sqrt(n())))

ggplot(surv_soil_root_precip, aes(x=precip,y=mean,fill=soil_root))+
  geom_bar(position="dodge",stat="identity", color="black",aes(fill=soil_root))+
  scale_fill_manual(values = c("gray", "white", "black"),labels=c("Bulk","Sterile","Rhizo"))+
  scale_x_discrete(labels=c("Ambient","Drought"))+ylab("prop Pots with at least one germinate")+
  theme_bw()

####UPDATE BELOW
#now lets look at the total number of germinates given there was at least one germinate


fin_dataSG_seed_surv_trt_w_germ=subset(fin_dataSG_seed_surv_trt, germinated>0)

hist(fin_dataSG_seed_surv_trt_w_germ$tot_num_germ)

seed_num_germ_model= lm(tot_num_germ~precip*soil_root+(1|block), data= fin_dataSG_seed_surv_trt_w_germ)
qqPlot(resid(seed_num_germ_model))
hist(resid(seed_num_germ_model))

anova(seed_num_germ_model, type=3)
emmeans(seed_num_germ_model, ~precip)
emmeans(seed_num_germ_model, pairwise~soil_root)

fin_dataSG_seed_surv_trt_w_germ_soil_root_g=group_by(fin_dataSG_seed_surv_trt_w_germ, soil_root)
overal_num_germ_soil_root=summarise_at(fin_dataSG_seed_surv_trt_w_germ_soil_root_g, "tot_num_germ", funs(mean,sd))

ggplot(overal_num_germ_soil_root, aes(soil_root,mean))+geom_bar(stat="identity", color="black",aes(fill=soil_root))+
  scale_fill_manual(values = c("gray", "white", "black"),labels=c("Bulk","Sterile","Rhizo"))+
  scale_x_discrete(limits=c("S.B", "L.B", "L.R"),labels=c("Sterile","Bulk","Rhizo"))+ylab("Num germ given one germ")



fin_dataSG_seed_surv_trt_precip_g=group_by(fin_dataSG_seed_surv_trt, precip)
overal_germ_precip=summarise_at(fin_dataSG_seed_surv_trt_precip_g, "germinated", mean)

ggplot(overal_germ_precip, aes(precip,germinated))+geom_bar(stat="identity", color="black",aes(fill=precip))+
  scale_fill_manual(values = c("gray", "white"),labels=c("Ambient","Drought"))+
  scale_x_discrete(labels=c("Ambient","Drought"))+ylab("At least one germinate")


###UPDATE ABOVE


summary(fin_dataSG_seed_surv_trt_pw)

data_SG_biomass_seed_surv_trt_pw=merge(fin_dataSG_seed_surv_trt_pw, data_SG_biomass, by="Plant_Number", all.x = T)
summary(data_SG_biomass_seed_surv_trt_pw)


#now lets look at total biomass of germinated 




seed_biomas_model_pw= lm(log(total_biomass)~end_per_ch_potweight_bio_crt*soil_root+(1|block), data= data_SG_biomass_seed_surv_trt_pw)
qqPlot(resid(seed_biomas_model_pw))
hist(resid(seed_biomas_model_pw))
shapiro.test(resid(seed_biomas_model_pw))
anova(seed_biomas_model_pw, type=3)
#nada sig





fin_dataSG_biomass_seed_surv_trt_given_germ_soil_root_precip_g=fin_dataSG_biomass_seed_surv_trt_given_germ %>% group_by(soil_root,precip)
total_bio_soil_root_precip=summarise_at(fin_dataSG_biomass_seed_surv_trt_given_germ_soil_root_precip_g, 
                                        "total_biomass", funs(n(),mean,sd,se=sd(.)/sqrt(n())))

ggplot(total_bio_soil_root_precip, aes(x=precip,y=mean, ymin = mean-se, ymax= mean+se,fill=soil_root))+
  geom_bar(position="dodge",stat="identity", color="black",aes(fill=soil_root))+
  geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+
  scale_fill_manual(values = c("gray", "white", "black"),labels=c("Bulk","Sterile","Rhizo"))+
  scale_x_discrete(labels=c("Ambient","Drought"))+ylab("Total biomass given gemination")+
  geom_text(aes(y=mean+se+0.02, label=n),position=position_dodge(width=0.9))+theme_bw()


#Root:shoot ratio

data_SG_biomass_seed_surv_trt_pw$root_shoot=data_SG_biomass_seed_surv_trt_pw$root_weight_g/data_SG_biomass_seed_surv_trt_pw$shoot_weight_g
summary(data_SG_biomass_seed_surv_trt_pw)


seed_root_shoot_model_pw= lm(log(root_shoot)~end_per_ch_potweight_bio_crt*soil_root+(1|block), data= data_SG_biomass_seed_surv_trt_pw)
qqPlot(resid(seed_root_shoot_model_pw))
hist(resid(seed_root_shoot_model_pw))

anova(seed_root_shoot_model_pw, type=3)
#end_per_ch_potweight_bio_crt            2.3982  1  3.0157 0.09152 .


ggplot(data_SG_biomass_seed_surv_trt_pw, aes(x=end_per_ch_potweight_bio_crt,y=log(root_shoot)))+geom_point(aes(fill=soil_root), size=4, shape=21)+
  geom_smooth(method='lm')+theme_bw()+scale_fill_manual(values = c("gray", "white", "black"),labels=c("Bulk","Sterile","Rhizo"))

#Roots alone
summary(fin_dataSG_biomass_seed_surv_trt_given_germ)
seed_root_model= lm(log(root_weight_g)~precip*soil_root+(1|block), data= fin_dataSG_biomass_seed_surv_trt_given_germ)
qqPlot(resid(seed_root_model))
hist(resid(seed_root_model))

anova(seed_root_model, type=3)
#precip:soil_root  10.773  2   6.1438  0.005275 ** 
#soil_root          5.526  2   3.1512  0.055528 .  

emmeans(seed_root_model, pairwise~soil_root)
emmeans(seed_root_model, pairwise~soil_root|precip)
#precip = A:
#contrast    estimate        SE df t.ratio p.value
#L.B - S.B  0.5289687 0.4822125 34   1.097  0.5225
#L.B - L.R -0.3642337 0.3717826 34  -0.980  0.5945
#S.B - L.R -0.8932024 0.4622551 34  -1.932  0.1451

#precip = D:
#  contrast    estimate        SE df t.ratio p.value
#L.B - S.B -2.6455949 0.7222441 34  -3.663  0.0024
#L.B - L.R -0.5686425 0.4606463 34  -1.234  0.4415
#S.B - L.R  2.0769524 0.7786334 34   2.667  0.0304



fin_dataSG_biomass_seed_surv_trt_given_germ_soil_root_precip_g=fin_dataSG_biomass_seed_surv_trt_given_germ %>% group_by(soil_root,precip)
root_soil_root_precip=summarise_at(fin_dataSG_biomass_seed_surv_trt_given_germ_soil_root_precip_g, 
                                   "root_weight_g", funs(n(),mean,sd,se=sd(.)/sqrt(n())))

ggplot(root_soil_root_precip, aes(x=precip,y=mean, ymin = mean-se, ymax= mean+se,fill=soil_root))+
  geom_bar(position="dodge",stat="identity", color="black",aes(fill=soil_root))+
  geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+
  scale_fill_manual(values = c("gray", "white", "black"),labels=c("Bulk","Sterile","Rhizo"))+
  scale_x_discrete(labels=c("Ambient","Drought"))+ylab("Root biomass given gemination")+
  geom_text(aes(y=mean+se+0.02, label=n),position=position_dodge(width=0.9))+theme_bw()




#let's remove the sandy roots

#now lets look at total biomass of germinated transs
fin_dataSG_biomass_seed_surv_trt_given_germ=subset(data_SG_biomass_seed_surv_trt, germinated>0)
fin_dataSG_biomass_seed_surv_trt_given_germ_N_sandy=subset(fin_dataSG_biomass_seed_surv_trt_given_germ, sandy_root=="N")
summary(fin_dataSG_biomass_seed_surv_trt_given_germ_N_sandy)

seed_biomas_model_N_sandy= lm(log(total_biomass)~precip*soil_root+(1|block), data= fin_dataSG_biomass_seed_surv_trt_given_germ_N_sandy)
qqPlot(resid(seed_biomas_model_N_sandy))
hist(resid(seed_biomas_model_N_sandy))

anova(seed_biomas_model_N_sandy, type=3)
#precip             1.470  1   2.9040  0.098696 .  
#soil_root          2.761  2   2.7274  0.081602 .  
#precip:soil_root   7.339  2   7.2505  0.002699 ** 

emmeans(seed_biomas_model_N_sandy, ~precip)
emmeans(seed_biomas_model_N_sandy, pairwise~soil_root|precip)
#$contrasts
#precip = A:
#  contrast  estimate    SE df t.ratio p.value
#L.B - S.B   0.7512 0.441 30  1.704  0.2203 
#L.B - L.R   0.0864 0.323 30  0.267  0.9614 
#S.B - L.R  -0.6647 0.429 30 -1.549  0.2830 

#precip = D:
#  contrast  estimate    SE df t.ratio p.value
#L.B - S.B  -2.7719 0.788 30 -3.516  0.0039 
#L.B - L.R  -0.5114 0.382 30 -1.340  0.3849 
#S.B - L.R   2.2606 0.829 30  2.727  0.0277 



fin_dataSG_biomass_seed_surv_trt_N_sandy_w_germ_soil_root_precip_g=fin_dataSG_biomass_seed_surv_trt_given_germ_N_sandy %>% group_by(soil_root,precip)
total_bio_N_sandy_soil_root_precip=summarise_at(fin_dataSG_biomass_seed_surv_trt_N_sandy_w_germ_soil_root_precip_g, 
                                                "total_biomass", funs(n(),mean,sd,se=sd(.)/sqrt(n())))

ggplot(total_bio_N_sandy_soil_root_precip, aes(x=precip,y=mean, ymin = mean-se, ymax= mean+se,fill=soil_root))+
  geom_bar(position="dodge",stat="identity", color="black",aes(fill=soil_root))+
  geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+
  scale_fill_manual(values = c("gray", "white", "black"),labels=c("Bulk","Sterile","Rhizo"))+
  scale_x_discrete(labels=c("Ambient","Drought"))+ylab("Total biomass given gemination (with sandy samples removed)")+
  geom_text(aes(y=mean+0.04, label=n),position=position_dodge(width=0.9))+theme_bw()



#Roots alone
summary(fin_dataSG_biomass_seed_surv_trt_given_germ_N_sandy)
seed_root_model_N_sandy= lm(log(root_weight_g)~precip*soil_root+(1|block), data= fin_dataSG_biomass_seed_surv_trt_given_germ_N_sandy)
qqPlot(resid(seed_root_model_N_sandy))
hist(resid(seed_root_model_N_sandy))

anova(seed_root_model_N_sandy, type=3)
#precip:soil_root   8.064  2   6.4893   0.00455 ** 
#precip             2.007  1   3.2309   0.08233 .   

emmeans(seed_root_model_N_sandy, pairwise~soil_root)
emmeans(seed_root_model_N_sandy, pairwise~soil_root|precip)
#precip = A:
#contrast  estimate    SE df t.ratio p.value
#L.B - S.B    0.811 0.488 30  1.660  0.2369 
#L.B - L.R    0.105 0.358 30  0.292  0.9542 
#S.B - L.R   -0.706 0.475 30 -1.486  0.3118 

#precip = D:
#  contrast  estimate    SE df t.ratio p.value
#L.B - S.B   -2.873 0.873 30 -3.289  0.0070 
#L.B - L.R   -0.563 0.423 30 -1.331  0.3897 
#S.B - L.R    2.310 0.919 30  2.515  0.0448 



fin_dataSG_biomass_seed_surv_trt_N_sandy_w_germ_soil_root_precip_g=fin_dataSG_biomass_seed_surv_trt_given_germ_N_sandy %>% group_by(soil_root,precip)
root_N_sandy_soil_root_precip=summarise_at(fin_dataSG_biomass_seed_surv_trt_N_sandy_w_germ_soil_root_precip_g, 
                                           "root_weight_g", funs(n(),mean,sd,se=sd(.)/sqrt(n())))

ggplot(root_N_sandy_soil_root_precip, aes(x=precip,y=mean, ymin = mean-se, ymax= mean+se,fill=soil_root))+
  geom_bar(position="dodge",stat="identity", color="black",aes(fill=soil_root))+
  geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+
  scale_fill_manual(values = c("gray", "white", "black"),labels=c("Bulk","Sterile","Rhizo"))+
  scale_x_discrete(labels=c("Ambient","Drought"))+ylab("Root biomass given gemination (with sandy samples removed)")+
  geom_text(aes(y=mean+0.04, label=n),position=position_dodge(width=0.9))+theme_bw()




#Shoots alone
summary(fin_dataSG_biomass_seed_surv_trt_given_germ)
seed_shoot_model= lm(log(shoot_weight_g)~precip*soil_root+(1|block), data= fin_dataSG_biomass_seed_surv_trt_given_germ)
qqPlot(resid(seed_shoot_model))
hist(resid(seed_shoot_model))

anova(seed_shoot_model, type=3)
#precip             2.06  1   10.0215  0.003257 ** 
#soil_root          1.63  2    3.9621  0.028395 *  
#(1|block)   4.70  4    5.7232  0.001240 ** 
#precip:soil_root   3.28  2    7.9750  0.001446 ** 

emmeans(seed_shoot_model, ~precip)
emmeans(seed_shoot_model, pairwise~soil_root)
emmeans(seed_shoot_model, pairwise~soil_root|precip)
#precip = A:
# L.B - S.B  0.14946329 0.2537568 34   0.589  0.8269
#L.B - L.R -0.25671123 0.1956448 34  -1.312  0.3983
#S.B - L.R -0.40617452 0.2432545 34  -1.670  0.2314

#precip = D:
#  contrast     estimate        SE df t.ratio p.value
#L.B - S.B -1.41470888 0.3800696 34  -3.722  0.0020
#L.B - L.R  0.08175893 0.2424079 34   0.337  0.9393
#S.B - L.R  1.49646782 0.4097436 34   3.652  0.0024



fin_dataSG_biomass_seed_surv_trt_given_germ_soil_root_precip_g=fin_dataSG_biomass_seed_surv_trt_given_germ %>% group_by(soil_root,precip)
shoot_soil_root_precip=summarise_at(fin_dataSG_biomass_seed_surv_trt_given_germ_soil_root_precip_g, 
                                    "shoot_weight_g", funs(n(),mean,sd,se=sd(.)/sqrt(n())))

ggplot(shoot_soil_root_precip, aes(x=precip,y=mean, ymin = mean-se, ymax= mean+se,fill=soil_root))+
  geom_bar(position="dodge",stat="identity", color="black",aes(fill=soil_root))+
  geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+
  scale_fill_manual(values = c("gray", "white", "black"),labels=c("Bulk","Sterile","Rhizo"))+
  scale_x_discrete(labels=c("Ambient","Drought"))+ylab("Shoot biomass given gemination")+
  geom_text(aes(y=mean+se+0.001, label=n),position=position_dodge(width=0.9))+theme_bw()


#####Soil nitrogen####

summary(SG_inorg_N_trt)
SG_inorg_N_trt$soil_root=with(SG_inorg_N_trt, interaction(soil_status,root_association))

SG_inorg_N_seed_trt= subset(SG_inorg_N_trt, life_stage=="S")
summary(SG_inorg_N_seed_trt)

#nitrate
SG_soil_nit_model= lm((ug_N_NO3_g_dry_soil)~precip*soil_root+(1|block), data= SG_inorg_N_seed_trt)
qqPlot(resid(SG_soil_nit_model))
hist(resid(SG_soil_nit_model))
shapiro.test(resid(SG_soil_nit_model))
#0.9866
anova(SG_soil_nit_model, type=3)
#precip             9.606  1  19.7580   0.01129 *  
#soil_root          8.119  2   8.3489   0.03735 * 
#precip:soil_root   4.949  2   5.0893   0.07959 .
SG_inorg_N_seed_trt %>% group_by(precip)  %>% summarise_at("ug_N_NO3_g_dry_soil", ~mean(.))
emmeans(SG_soil_nit_model, ~precip)
emmeans(SG_soil_nit_model, pairwise~soil_root)
emmeans(SG_soil_nit_model, pairwise~soil_root|precip)
#$contrasts
#precip = A:
#  contrast  estimate    SE df t.ratio p.value
#L.B - S.B    0.364 0.750  4  0.485  0.8821 
#L.B - L.R   -0.283 0.551  4 -0.513  0.8695 
#S.B - L.R   -0.646 0.664  4 -0.973  0.6293 

#precip = D:
#  contrast  estimate    SE df t.ratio p.value
#L.B - S.B    4.743 1.239  4  3.829  0.0399 
#L.B - L.R    2.832 0.897  4  3.155  0.0721 
#S.B - L.R   -1.911 1.075  4 -1.778  0.2868 

SG_soil_nit_model_no_block= lm((ug_N_NO3_g_dry_soil)~precip*soil_root, data= SG_inorg_N_seed_trt)
qqPlot(resid(SG_soil_nit_model_no_block))
hist(resid(SG_soil_nit_model_no_block))
shapiro.test(resid(SG_soil_nit_model_no_block))
#0.7402
anova(SG_soil_nit_model_no_block, type=3)
#precip            11.466  1  38.967 0.0004273 ***
#soil_root          9.913  2  16.845 0.0021117 ** 
#precip:soil_root   8.826  2  14.998 0.0029466 ** 

emmeans(SG_soil_nit_model_no_block, pairwise~soil_root|precip)

AIC(SG_soil_nit_model,SG_soil_nit_model_no_block)

SG_inorg_N_seed_trt_soil_root_precip_g=SG_inorg_N_seed_trt %>% group_by(soil_root,precip)
soil_nit_seed_trt=summarise_at(SG_inorg_N_seed_trt_soil_root_precip_g, 
                                    "ug_N_NO3_g_dry_soil", funs(n(),mean,sd,se=sd(.)/sqrt(n())))

ggplot(soil_nit_seed_trt, aes(x=precip,y=mean, ymin = mean-se, ymax= mean+se,fill=soil_root))+
  geom_bar(position="dodge",stat="identity", color="black",aes(fill=soil_root))+
  geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+
  scale_fill_manual(values = c("gray", "white", "black"),labels=c("Bulk","Sterile","Rhizo"))+
  scale_x_discrete(labels=c("Ambient","Drought"))+ylab("Soil nitrate (ug per g soil)")+
  geom_text(aes(y=mean+se+0.3, label=n),position=position_dodge(width=0.9))+theme_bw()

#####Seed nitrate reformated graph####

(seed_soil_nitrate_p2=ggplot(soil_nit_seed_trt, aes(x=precip,y=mean, ymin = mean-se, ymax= mean+se,fill=factor(soil_root,levels = treatment_order)))+
   geom_bar(position="dodge",stat="identity", color="black",aes(fill=factor(soil_root,levels = treatment_order)))+
   geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+ylim(c(0,9))+
   scale_fill_manual(values = c( "white","lightgray", "darkgrey"),labels=c("Sterile","Bulk","Rhizo"), limits=treatment_order, name=NULL)+
   scale_x_discrete(labels=c("Ambient","Drought"), name=NULL)+ylab("Soil nitrate (ug per g soil)")+
   geom_text(aes(y=0.5, label=n),size=10,position=position_dodge(width=0.9))+theme_bw()+
   theme(axis.title = element_text(size = 23), axis.text = element_text(size = 23),
         legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank()))
#800x750
(seed_soil_nitrate_p=ggplot(soil_nit_seed_trt, aes(x=precip,y=mean, ymin = mean-se, ymax= mean+se,fill=factor(soil_root,levels = treatment_order)))+
   geom_bar(position="dodge",stat="identity", color="black",aes(fill=factor(soil_root,levels = treatment_order)))+
   geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+ylim(c(0,9))+
   scale_fill_manual(values = c( "white","light gray", "dark grey"),labels=c("Sterile","Bulk","Rhizo"), limits=treatment_order, name=NULL)+
   scale_x_discrete(labels=c("Ambient","Drought"), name=NULL)+ylab("Soil nitrate (ug per g soil)")+
   geom_text(aes(y=0.5, label=n),size=10,position=position_dodge(width=0.9))+theme_bw()+
   theme(axis.title.x = element_text(size = 23), axis.text.x = element_text(size = 23),
         axis.title.y = element_blank(), axis.text.y = element_blank(),
         legend.position = c(0.85,.9), legend.text=element_text(size=20),
         legend.background = element_rect(size=0.5,linetype="solid",colour ="black")))

#Amonnium
SG_soil_amm_model= lm((ug_N_NH4_g_dry_soil_negto0)~precip*soil_root+(1|block), data= SG_inorg_N_seed_trt)
qqPlot(resid(SG_soil_amm_model))
hist(resid(SG_soil_amm_model))
#plot(SG_soil_amm_model)
shapiro.test(resid(SG_soil_amm_model))
#0.8105
anova(SG_soil_amm_model, type=3)
#Nada sig



emmeans(SG_soil_amm_model, ~precip)
emmeans(SG_soil_amm_model, pairwise~soil_root)
emmeans(SG_soil_amm_model, pairwise~soil_root|precip)


SG_soil_amm_model_no_block= lm((ug_N_NH4_g_dry_soil_negto0)~precip*soil_root, data= SG_inorg_N_seed_trt)
qqPlot(resid(SG_soil_amm_model_no_block))
hist(resid(SG_soil_amm_model_no_block))
#plot(SG_soil_amm_model)
shapiro.test(resid(SG_soil_amm_model_no_block))
#7.899e-05
anova(SG_soil_amm_model_no_block, type=3)
#Nada sig

emmeans(SG_soil_amm_model, pairwise~soil_root|precip)

AIC(SG_soil_amm_model,SG_soil_amm_model_no_block)

SG_inorg_N_seed_trt_soil_root_precip_g=SG_inorg_N_seed_trt %>% group_by(soil_root,precip)
soil_amm_seed_trt=summarise_at(SG_inorg_N_seed_trt_soil_root_precip_g, 
                               "ug_N_NH4_g_dry_soil_negto0", funs(n(),mean,sd,se=sd(.)/sqrt(n())))

ggplot(soil_amm_seed_trt, aes(x=precip,y=mean, ymin = mean-se, ymax= mean+se,fill=soil_root))+
  geom_bar(position="dodge",stat="identity", color="black",aes(fill=soil_root))+
  geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+
  scale_fill_manual(values = c("gray", "white", "black"),labels=c("Bulk","Sterile","Rhizo"))+
  scale_x_discrete(labels=c("Ambient","Drought"))+ylab("Soil ammonium (ug per g soil)")+
  geom_text(aes(y=mean+se+0.3, label=n),position=position_dodge(width=0.9))+theme_bw()


#####Seed Ammonium reformated graph####
(seed_soil_amm_p2=ggplot(soil_amm_seed_trt, aes(x=precip,y=mean, ymin = mean-se, ymax= mean+se,fill=factor(soil_root,levels = treatment_order)))+
   geom_bar(position="dodge",stat="identity", color="black",aes(fill=factor(soil_root,levels = treatment_order)))+
   geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+ylim(c(0,4.5))+
   scale_fill_manual(values = c( "white","lightgray", "darkgrey"),labels=c("Sterile","Bulk","Rhizo"), limits=treatment_order, name=NULL)+
   scale_x_discrete(labels=c("Ambient","Drought"), name=NULL)+ylab("Soil ammonium (ug per g soil)")+
   geom_text(aes(y=mean+se+0.3, label=n),size=10,position=position_dodge(width=0.9))+theme_bw()+
   theme(axis.title = element_text(size = 23), axis.text = element_text(size = 23),
         legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

#soil_ammonium_germination 800x750
(seed_soil_amm_p=ggplot(soil_amm_seed_trt, aes(x=precip,y=mean, ymin = mean-se, ymax= mean+se,fill=factor(soil_root,levels = treatment_order)))+
   geom_bar(position="dodge",stat="identity", color="black",aes(fill=factor(soil_root,levels = treatment_order)))+
   geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+ylim(c(0,4.5))+
   scale_fill_manual(values = c( "white","light gray", "dark grey"),labels=c("Sterile","Bulk","Rhizo"), limits=treatment_order, name=NULL)+
   scale_x_discrete(labels=c("Ambient","Drought"), name=NULL)+ylab("Soil ammonium (ug per g soil)")+
   geom_text(aes(y=mean+se+0.3, label=n),size=10,position=position_dodge(width=0.9))+theme_bw()+
   theme(axis.title.x = element_text(size = 23), axis.text.x = element_text(size = 23),
         axis.title.y = element_blank(), axis.text.y = element_blank(),
         legend.position = c(0.85,.9), legend.text=element_text(size=20),
         legend.background = element_rect(size=0.5,linetype="solid",colour ="black")))


#Using soil moisture instead of the precip treat 

dataSG_potweight_end=dataSG_potweight_trt_seed_sub_out_end[,c("Plant_Number","per_ch_potweight","per_ch_potweight_bio_crt")]
colnames(dataSG_potweight_end)[2]="end_per_ch_potweight"
colnames(dataSG_potweight_end)[3]="end_per_ch_potweight_bio_crt"
nrow(dataSG_potweight_end)

SG_inorg_N_seed_trt_pw=merge(SG_inorg_N_seed_trt,dataSG_potweight_end, by="Plant_Number", all.x = T)

#nitrate
SG_soil_nit_model_pw= lm((ug_N_NO3_g_dry_soil)~end_per_ch_potweight_bio_crt*soil_root+(1|block), data= SG_inorg_N_seed_trt_pw)
qqPlot(resid(SG_soil_nit_model_pw))
hist(resid(SG_soil_nit_model_pw))

anova(SG_soil_nit_model_pw, type=3)
#end_per_ch_potweight_bio_crt           7.8061  1  5.1416 0.08596 .

ggplot(SG_inorg_N_seed_trt_pw, aes(x=end_per_ch_potweight,y=ug_N_NO3_g_dry_soil))+geom_point(aes(fill=soil_root), size=4, shape=21)+
  geom_smooth(method='lm')+theme_bw()+scale_fill_manual(values = c("gray", "white", "black"),labels=c("Bulk","Sterile","Rhizo"))

#Amonnium
SG_soil_amm_model_pw= lm((ug_N_NH4_g_dry_soil_negto0)~end_per_ch_potweight_bio_crt*soil_root+(1|block), data= SG_inorg_N_seed_trt_pw)
qqPlot(resid(SG_soil_amm_model_pw))
hist(resid(SG_soil_amm_model_pw))
#plot(SG_soil_amm_model_pw)
anova(SG_soil_amm_model_pw, type=3)
#soil_root                              4.1723  2  5.5067 0.07098 .

emmeans(SG_soil_amm_model_pw, pairwise~soil_root)


#Is change in pot weight correlated with gravimetric water

seed_gravimetric_mod=lm(percent_soil_moisture_dry_weight~end_per_ch_potweight_bio_crt,data=SG_inorg_N_seed_trt_pw)
qqPlot(resid(seed_gravimetric_mod))
hist(resid(seed_gravimetric_mod))
#plot(seed_gravimetric_mod)
summary(seed_gravimetric_mod)
plot(SG_inorg_N_seed_trt_pw$end_per_ch_potweight,SG_inorg_N_seed_trt_pw$percent_soil_moisture_dry_weight)

ggplot(SG_inorg_N_seed_trt_pw, aes(x=end_per_ch_potweight_bio_crt,y=percent_soil_moisture_dry_weight))+geom_point(aes(fill=soil_root), size=4, shape=21)+
  theme_bw()+scale_fill_manual(values = c("gray", "white", "black"),labels=c("Bulk","Sterile","Rhizo"))+
  geom_text(aes(y=percent_soil_moisture_dry_weight+.7, x= end_per_ch_potweight,label=Plant_Number))


#Treat effects on the gravimetric 



SG_soil_VWC_model= lm((percent_soil_moisture_dry_weight)~precip*soil_root+(1|block), data= SG_inorg_N_seed_trt)
qqPlot(resid(SG_soil_VWC_model))
hist(resid(SG_soil_VWC_model))
plot(SG_soil_VWC_model)
anova(SG_soil_VWC_model, type=3)
#Nada sig



emmeans(SG_soil_VWC_model, ~precip)
emmeans(SG_soil_VWC_model, pairwise~soil_root)
emmeans(SG_soil_VWC_model, pairwise~soil_root|precip)


#####SEED GERMINATION PORTION OF THE EXPERIEMENT ENDING####

#####TRANSPLANT PORTION OF THE EXPERIEMENT BEGINING####
#let's look at the transplanted transs

#now I want to see if surv in the first week of the experiement is affect by treatments
summary(dataSG_trans_surv_trt)
dataSG_trans_surv_trt_g=group_by(dataSG_trans_surv_trt, Plant_Number)
total_transs_transplanted=summarise_at(dataSG_trans_surv_trt_g, "Transplant_height_cm", sum,na.rm = TRUE)
colnames(total_transs_transplanted)[2]="tot_Transplant_height_cm"
summary(total_transs_transplanted)
nrow(total_transs_transplanted)

fin_dataSG_trans_surv_trt=subset(dataSG_trans_surv_trt, date=="7_9")
nrow(fin_dataSG_trans_surv_trt)
fin_dataSG_trans_surv_trt_merge=merge(fin_dataSG_trans_surv_trt,total_transs_transplanted, by="Plant_Number")
nrow(fin_dataSG_trans_surv_trt_merge)
fin_dataSG_trans_surv_trt_merge$int_surv=fin_dataSG_trans_surv_trt_merge$tot_Transplant_height_cm
fin_dataSG_trans_surv_trt_merge$int_surv[fin_dataSG_trans_surv_trt_merge$int_surv==0]=1
fin_dataSG_trans_surv_trt_merge$int_surv[fin_dataSG_trans_surv_trt_merge$int_surv!=1]=0
fin_dataSG_trans_surv_trt_merge$soil_root=with(fin_dataSG_trans_surv_trt_merge, interaction(soil_status,root_association))

summary(fin_dataSG_trans_surv_trt_merge)


hist(fin_dataSG_trans_surv_trt_merge$int_surv)

int_surv_model= glm(int_surv~precip*soil_root+(1|block), data= fin_dataSG_trans_surv_trt_merge, family = binomial)
qqPlot(resid(int_surv_model))
hist(resid(int_surv_model))

anova(int_surv_model, type=3)
emmeans(int_surv_model, ~precip)
emmeans(int_surv_model, pairwise~soil_root)
emmeans(int_surv_model, pairwise~precip|soil_root)

fin_fin_dataSG_trans_surv_trt_soil_root_precip_g=fin_dataSG_trans_surv_trt_merge %>% group_by(soil_root,precip)
int_surv_soil_root_precip=summarise_at(fin_fin_dataSG_trans_surv_trt_soil_root_precip_g, "int_surv", mean)

ggplot(int_surv_soil_root_precip, aes(precip,int_surv))+geom_bar(position="dodge",stat="identity", color="black",aes(fill=soil_root))+
  scale_fill_manual(values = c("gray", "white", "black"),labels=c("Bulk","Sterile","Rhizo"))+
  scale_x_discrete(labels=c("Ambient","Drought"))+ylab("First week Transplant Survival")


summary(fin_dataSG_trans_surv_trt)
summary(data_SG_biomass)
data_SG_biomass_trans_surv_trt=merge(fin_dataSG_trans_surv_trt, data_SG_biomass, by="Plant_Number", all.x = T)
data_SG_biomass_trans_surv_trt$soil_root=with(data_SG_biomass_trans_surv_trt, interaction(soil_status,root_association))
summary(data_SG_biomass_trans_surv_trt)

#now lets look at biomass with dead as zero

data_SG_biomass_trans_surv_trt$total_biomass[data_SG_biomass_trans_surv_trt$total_biomass==0]=NA
data_SG_biomass_trans_surv_trt$total_biomass_0=data_SG_biomass_trans_surv_trt$total_biomass
data_SG_biomass_trans_surv_trt$total_biomass_0[is.na(data_SG_biomass_trans_surv_trt$total_biomass_0)]=0

data_SG_biomass_trans_surv_trt$root_weight_g_0=data_SG_biomass_trans_surv_trt$root_weight_g
data_SG_biomass_trans_surv_trt$root_weight_g_0[is.na(data_SG_biomass_trans_surv_trt$root_weight_g_0)]=0

data_SG_biomass_trans_surv_trt$shoot_weight_g_0=data_SG_biomass_trans_surv_trt$shoot_weight_g
data_SG_biomass_trans_surv_trt$shoot_weight_g_0[is.na(data_SG_biomass_trans_surv_trt$shoot_weight_g_0)]=0

summary(data_SG_biomass_trans_surv_trt)


#now lets look at total biomass



trans_biomas_model= lm(log(total_biomass_0+.001)~precip*soil_root+(1|block), data= data_SG_biomass_trans_surv_trt)
qqPlot(resid(trans_biomas_model))
hist(resid(trans_biomas_model))
shapiro.test(resid(trans_biomas_model))
#0.348
anova(trans_biomas_model, type=3)
#precip            31.95  1  13.4201 0.0004457 ***
#soil_root         23.30  2   4.8939 0.0098833 ** 
#precip:soil_root  15.88  2   3.3353 0.0406192 *  

emmeans(trans_biomas_model, pairwise~soil_root)
data_SG_biomass_trans_surv_trt %>% group_by(precip) %>% summarise_at("total_biomass_0", funs(n(),mean,sd,se=sd(.)/sqrt(n())))

emmeans(trans_biomas_model, pairwise~soil_root|precip)
#precip = A:
#contrast     estimate        SE df t.ratio p.value
#L.B - S.B -1.86324726 0.5633964 80  -3.307  0.0040
#L.B - L.R  0.14030501 0.5633964 80   0.249  0.9664
#S.B - L.R  2.00355226 0.5633964 80   3.556  0.0018

#precip = D:
#  contrast     estimate        SE df t.ratio p.value
#L.B - S.B  0.04493185 0.5633964 80   0.080  0.9965
#L.B - L.R  0.42719232 0.5633964 80   0.758  0.7295
#S.B - L.R  0.38226047 0.5633964 80   0.678  0.7767

emmeans(trans_biomas_model, pairwise~soil_root*precip, adjust="fdr")
"$contrasts
 contrast      estimate    SE df t.ratio p.value
L.B,A - S.B,A  -1.8632 0.563 80 -3.307  0.0042 
L.B,A - L.R,A   0.1403 0.563 80  0.249  0.8614 
L.B,A - L.B,D   0.4599 0.563 80  0.816  0.6493 
L.B,A - S.B,D   0.5048 0.563 80  0.896  0.6493 
L.B,A - L.R,D   0.8871 0.563 80  1.575  0.2983 
S.B,A - L.R,A   2.0036 0.563 80  3.556  0.0024 
S.B,A - L.B,D   2.3232 0.563 80  4.123  0.0005 
S.B,A - S.B,D   2.3681 0.563 80  4.203  0.0005 
S.B,A - L.R,D   2.7504 0.563 80  4.882  0.0001 
L.R,A - L.B,D   0.3196 0.563 80  0.567  0.6601 
L.R,A - S.B,D   0.3645 0.563 80  0.647  0.6493 
L.R,A - L.R,D   0.7468 0.563 80  1.326  0.4045 
L.B,D - S.B,D   0.0449 0.563 80  0.080  0.9366 
L.B,D - L.R,D   0.4272 0.563 80  0.758  0.6493 
S.B,D - L.R,D   0.3823 0.563 80  0.678  0.6493 

Results are averaged over the levels of: block 
Results are given on the log (not the response) scale. 
P value adjustment: fdr method for 15 tests "


trans_biomas_model_no_block= lm(log(total_biomass_0+.001)~precip*soil_root, data= data_SG_biomass_trans_surv_trt)
qqPlot(resid(trans_biomas_model_no_block))
hist(resid(trans_biomas_model_no_block))
shapiro.test(resid(trans_biomas_model_no_block))
#0.001249
anova(trans_biomas_model_no_block, type=3)
#precip            31.95  1  12.8430 0.0005671 ***
#soil_root         23.30  2   4.6835 0.0117928 *  
#precip:soil_root  15.88  2   3.1919 0.0461250 * 

emmeans(trans_biomas_model_no_block, pairwise~soil_root|precip)

AIC(trans_biomas_model,trans_biomas_model_no_block)

data_SG_biomass_trans_surv_trt %>% group_by(precip) %>% summarise_at("total_biomass_0", funs(n(),mean,sd,se=sd(.)/sqrt(n())))

(0.0717-0.174)/0.174

data_SG_biomass_trans_surv_trt %>% group_by(precip,soil_root) %>% summarise_at("total_biomass_0", funs(n(),mean,sd,se=sd(.)/sqrt(n())))


data_SG_biomass_trans_surv_trt_soil_root_precip_g=data_SG_biomass_trans_surv_trt %>% group_by(soil_root,precip)
trans_total_bio_soil_root_precip=summarise_at(data_SG_biomass_trans_surv_trt_soil_root_precip_g, 
                                        "total_biomass_0", funs(n(),mean,sd,se=sd(.)/sqrt(n())))

#Sterile drought versus ambient
(0.0876-0.403)/0.403
#Live Bulk drought versus ambient
(0.0712- 0.0653)/ 0.0653
#Live Bulk drought versus ambient
(0.0562-0.0548)/0.0548

ggplot(trans_total_bio_soil_root_precip, aes(x=precip,y=mean, ymin = mean-se, ymax= mean+se,fill=soil_root))+
  geom_bar(position="dodge",stat="identity", color="black",aes(fill=soil_root))+
  geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+
  scale_fill_manual(values = c("gray", "white", "black"),labels=c("Bulk","Sterile","Rhizo"))+
  scale_x_discrete(labels=c("Ambient","Drought"))+ylab("Total biomass transplants")+
  geom_text(aes(y=mean+se+0.02, label=n),position=position_dodge(width=0.9))+theme_bw()



#####Trans Bio Reformated bar graph #####

(trans_biomass_p=ggplot(trans_total_bio_soil_root_precip, aes(x=precip,y=mean, ymin = mean-se, ymax= mean+se,fill=factor(soil_root,levels = treatment_order)))+
   geom_bar(position="dodge",stat="identity", color="black",aes(fill=factor(soil_root,levels = treatment_order)))+
   geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+ylim(c(0,.95))+
   scale_fill_manual(values = c( "white","light gray", "dark grey"),labels=c("Sterile","Bulk","Rhizo"), limits=treatment_order, name=NULL)+
   scale_x_discrete(labels=c("Ambient","Drought"), name=NULL)+ylab("Total biomass (g)")+
   geom_text(aes(y=mean+se+0.04, label=n),size=10,position=position_dodge(width=0.9))+theme_bw()+
   theme(axis.title = element_text(size = 23), axis.text = element_text(size = 23),
         legend.position = c(0.85,.9), legend.text=element_text(size=20),
         legend.background = element_rect(size=0.5,linetype="solid",colour ="black")))


(trans_biomass_p=ggplot(trans_total_bio_soil_root_precip, aes(x=precip,y=mean, ymin = mean-se, ymax= mean+se,fill=factor(soil_root,levels = treatment_order)))+
    geom_bar(position="dodge",stat="identity", color="black",aes(fill=factor(soil_root,levels = treatment_order)))+
    geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+ylim(c(0,.95))+
    scale_fill_manual(values = c( "white","lightgray", "darkgrey"),labels=c("Sterile","Bulk","Rhizo"), limits=treatment_order, name=NULL)+
    scale_x_discrete(labels=c("Ambient","Drought"), name=NULL)+ylab("Seedling biomass (g)")+theme_bw()+
    theme(axis.title = element_text(size = 23), axis.text = element_text(size = 23),
          legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank()))
#800x750
#####Combined biomass graph#####
ggarrange(trans_biomass_p,seed_biomass_p,ncol = 2,  legend = "none")
#15x7.38

#####Presence at total biomass#####

data_SG_biomass_trans_surv_pres=subset(data_SG_biomass_trans_surv_trt, root_association =="B")

trans_biomas_model_pres= lm((total_biomass_0+.001)^(1/4)~precip*soil_root+(1|block), data= data_SG_biomass_trans_surv_pres)
qqPlot(resid(trans_biomas_model_pres))
hist(resid(trans_biomas_model_pres))
shapiro.test(resid(trans_biomas_model_pres))
#0.4747
anova(trans_biomas_model_pres, type=3)
#precip            0.3824  1  13.1022 0.0006683 ***
#soil_root         0.2985  1  10.2276 0.0023565 ** 
#precip:soil_root  0.2559  1   8.7667 0.0046136 ** 

emmeans(trans_biomas_model_pres, pairwise~soil_root)
data_SG_biomass_trans_surv_pres %>% group_by(precip) %>% summarise_at("total_biomass_0", funs(n(),mean,sd,se=sd(.)/sqrt(n())))

emmeans(trans_biomas_model_pres, pairwise~soil_root|precip)


emmeans(trans_biomas_model_pres, pairwise~soil_root*precip, adjust="fdr")



trans_biomas_model_pres_no_block= lm((total_biomass_0+.001)^(1/4)~precip*soil_root, data= data_SG_biomass_trans_surv_pres)
qqPlot(resid(trans_biomas_model_pres_no_block))
hist(resid(trans_biomas_model_pres_no_block))
shapiro.test(resid(trans_biomas_model_pres_no_block))
#0.254
boxCox(trans_biomas_model_pres_no_block)
anova(trans_biomas_model_pres_no_block, type=3)
#precip            0.3824  1  13.2458 0.0005969 ***
#soil_root         0.2985  1  10.3396 0.0021637 ** 
#precip:soil_root  0.2559  1   8.8628 0.0042933 ** 

emmeans(trans_biomas_model_pres_no_block, pairwise~soil_root*precip,adjust="none")

AIC(trans_biomas_model_pres,trans_biomas_model_pres_no_block)


#####Origin at total biomass#####

data_SG_biomass_trans_surv_orig=subset(data_SG_biomass_trans_surv_trt, soil_status =="L")
nrow(data_SG_biomass_trans_surv_orig)

trans_biomas_model_orig= lm((total_biomass_0+0.001)^(1/4)~precip*soil_root+(1|block), data= data_SG_biomass_trans_surv_orig)
qqPlot(resid(trans_biomas_model_orig))
hist(resid(trans_biomas_model_orig))
shapiro.test(resid(trans_biomas_model_orig))
#0.1223
anova(trans_biomas_model_orig, type=3)
#nada sig

emmeans(trans_biomas_model_orig, pairwise~soil_root)
data_SG_biomass_trans_surv_orig %>% group_by(precip) %>% summarise_at("total_biomass_0", funs(n(),mean,sd,se=sd(.)/sqrt(n())))

emmeans(trans_biomas_model_pres, pairwise~soil_root|precip)


emmeans(trans_biomas_model_pres, pairwise~soil_root*precip, adjust="fdr")



trans_biomas_model_orig_no_block= lm((total_biomass_0+.001)^(1/4)~precip*soil_root, data= data_SG_biomass_trans_surv_orig)
qqPlot(resid(trans_biomas_model_orig_no_block))
hist(resid(trans_biomas_model_orig_no_block))
shapiro.test(resid(trans_biomas_model_orig_no_block))
#0.0864
boxCox(trans_biomas_model_orig_no_block)
anova(trans_biomas_model_orig_no_block, type=3)
#nada sig

emmeans(trans_biomas_model_orig_no_block, pairwise~soil_root*precip,adjust="none")

AIC(trans_biomas_model_orig,trans_biomas_model_orig_no_block)

#Roots alone
summary(data_SG_biomass_trans_surv_trt)
trans_root_model= lm(log(root_weight_g_0+.001)~precip*soil_root+(1|block), data= data_SG_biomass_trans_surv_trt)
qqPlot(resid(trans_root_model))
hist(resid(trans_root_model))

anova(trans_root_model, type=3)
#precip             24.14  1   8.6958 0.00418 ** 
#soil_root          20.25  2   3.6482 0.03046 * 
#precip:soil_root   16.45  2   2.9639 0.05731 .  


emmeans(trans_root_model, pairwise~soil_root)
emmeans(trans_root_model, pairwise~soil_root|precip)
#precip = A:
#contrast    estimate        SE df t.ratio p.value
#L.B - S.B -1.8156429 0.6083664 80  -2.984  0.0104
#L.B - L.R  0.1377609 0.6083664 80   0.226  0.9721
#S.B - L.R  1.9534037 0.6083664 80   3.211  0.0054

#precip = D:
#  contrast    estimate        SE df t.ratio p.value
#L.B - S.B  0.1003581 0.6083664 80   0.165  0.9851
#L.B - L.R  0.3625702 0.6083664 80   0.596  0.8227
#S.B - L.R  0.2622120 0.6083664 80   0.431  0.9028



data_SG_biomass_trans_surv_trt_soil_root_precip_g=data_SG_biomass_trans_surv_trt %>% group_by(soil_root,precip)
trans_root_soil_root_precip=summarise_at(data_SG_biomass_trans_surv_trt_soil_root_precip_g, 
                                   "root_weight_g_0", funs(n(),mean,sd,se=sd(.)/sqrt(n())))

ggplot(trans_root_soil_root_precip, aes(x=precip,y=mean, ymin = mean-se, ymax= mean+se,fill=soil_root))+
  geom_bar(position="dodge",stat="identity", color="black",aes(fill=soil_root))+
  geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+
  scale_fill_manual(values = c("gray", "white", "black"),labels=c("Bulk","Sterile","Rhizo"))+
  scale_x_discrete(labels=c("Ambient","Drought"))+ylab("Root biomass transplants")+
  geom_text(aes(y=mean+se+0.02, label=n),position=position_dodge(width=0.9))+theme_bw()

#with sandy roots removed
colnames(data_SG_biomass_trans_surv_trt)
data_SG_biomass_trans_surv_trt_N_sandy=subset(data_SG_biomass_trans_surv_trt, sandy_root=="N")

#now lets look at total biomass



trans_biomas_model_N_sandy= lm(log(total_biomass_0+0.001)~precip*soil_root+(1|block), data= data_SG_biomass_trans_surv_trt_N_sandy)
qqPlot(resid(trans_biomas_model_N_sandy))
hist(resid(trans_biomas_model_N_sandy))

anova(trans_biomas_model_N_sandy, type=3)
#precip            19.52  1   9.0147  0.003623 ** 
#soil_root         26.57  2   6.1367  0.003384 **  

emmeans(trans_biomas_model_N_sandy, pairwise~soil_root)
emmeans(trans_biomas_model_N_sandy, pairwise~soil_root|precip)
#$contrasts
#precip = A:
#  contrast  estimate    SE df t.ratio p.value
#L.B - S.B   -1.791 0.547 76 -3.273  0.0045 
#L.B - L.R    0.140 0.537 76  0.261  0.9631 
#S.B - L.R    1.931 0.547 76  3.529  0.0020 

#precip = D:
#  contrast  estimate    SE df t.ratio p.value
#L.B - S.B   -0.399 0.558 76 -0.714  0.7560 
#L.B - L.R    0.236 0.547 76  0.432  0.9025 
#S.B - L.R    0.635 0.567 76  1.119  0.5052 



data_SG_biomass_trans_surv_trt_N_sandy_soil_root_precip_g=data_SG_biomass_trans_surv_trt_N_sandy %>% group_by(soil_root,precip)
trans_total_bio_N_sandy_soil_root_precip=summarise_at(data_SG_biomass_trans_surv_trt_N_sandy_soil_root_precip_g, 
                                                 "total_biomass_0", funs(n(),mean,sd,se=sd(.)/sqrt(n())))

ggplot(trans_total_bio_N_sandy_soil_root_precip, aes(x=precip,y=mean, ymin = mean-se, ymax= mean+se,fill=soil_root))+
  geom_bar(position="dodge",stat="identity", color="black",aes(fill=soil_root))+
  geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+
  scale_fill_manual(values = c("gray", "white", "black"),labels=c("Bulk","Sterile","Rhizo"))+
  scale_x_discrete(labels=c("Ambient","Drought"))+ylab("Total biomass transplants (with sandy samples removed)")+
  geom_text(aes(y=mean+se+0.02, label=n),position=position_dodge(width=0.9))+theme_bw()



#Roots alone
summary(data_SG_biomass_trans_surv_trt_N_sandy)
trans_root_model_N_sandy= lm(log(root_weight_g_0+.001)~precip*soil_root+(1|block), data= data_SG_biomass_trans_surv_trt_N_sandy)
qqPlot(resid(trans_root_model_N_sandy))
hist(resid(trans_root_model_N_sandy))

anova(trans_root_model_N_sandy, type=3)
#precip             14.10  1   5.3558 0.02336 *  
#soil_root          22.21  2   4.2191 0.01830 *  


emmeans(trans_root_model_N_sandy, pairwise~soil_root)
emmeans(trans_root_model_N_sandy, pairwise~soil_root|precip)
#$contrasts
#precip = A:
#  contrast  estimate    SE df t.ratio p.value
#L.B - S.B   -1.719 0.603 76 -2.849  0.0154 
#L.B - L.R    0.138 0.592 76  0.233  0.9706 
#S.B - L.R    1.857 0.603 76  3.077  0.0081 

#precip = D:
#  contrast  estimate    SE df t.ratio p.value
#L.B - S.B   -0.296 0.616 76 -0.481  0.8807 
#L.B - L.R    0.186 0.603 76  0.308  0.9491 
#S.B - L.R    0.482 0.625 76  0.770  0.7223 



data_SG_biomass_trans_surv_trt_N_sandy_soil_root_precip_g=data_SG_biomass_trans_surv_trt_N_sandy %>% group_by(soil_root,precip)
trans_root_N_sandy_soil_root_precip=summarise_at(data_SG_biomass_trans_surv_trt_N_sandy_soil_root_precip_g, 
                                            "root_weight_g_0", funs(n(),mean,sd,se=sd(.)/sqrt(n())))

ggplot(trans_root_N_sandy_soil_root_precip, aes(x=precip,y=mean, ymin = mean-se, ymax= mean+se,fill=soil_root))+
  geom_bar(position="dodge",stat="identity", color="black",aes(fill=soil_root))+
  geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+
  scale_fill_manual(values = c("gray", "white", "black"),labels=c("Bulk","Sterile","Rhizo"))+
  scale_x_discrete(labels=c("Ambient","Drought"))+ylab("Root biomass transplants (with sandy samples removed)")+
  geom_text(aes(y=mean+se+0.02, label=n),position=position_dodge(width=0.9))+theme_bw()




#Shoots alone
summary(data_SG_biomass_trans_surv_trt)
trans_shoot_model= lm(sqrt(shoot_weight_g_0+.001)~precip*soil_root+(1|block), data= data_SG_biomass_trans_surv_trt)
qqPlot(resid(trans_shoot_model))
hist(resid(trans_shoot_model))

anova(trans_shoot_model, type=3)
#precip           0.08441  1   67.5994 2.970e-12 ***
#soil_root        0.15209  2   60.8972 < 2.2e-16 ***
#precip:soil_root 0.05837  2   23.3708 1.016e-08 ***

emmeans(trans_shoot_model, ~precip)
emmeans(trans_shoot_model, pairwise~soil_root)
emmeans(trans_shoot_model, pairwise~soil_root|precip)
#precip = A:
# L.B - S.B  0.14946329 0.2537568 34   0.589  0.8269
#L.B - L.R -0.25671123 0.1956448 34  -1.312  0.3983
#S.B - L.R -0.40617452 0.2432545 34  -1.670  0.2314

#precip = D:
#  contrast     estimate        SE df t.ratio p.value
#L.B - S.B -1.41470888 0.3800696 34  -3.722  0.0020
#L.B - L.R  0.08175893 0.2424079 34   0.337  0.9393
#S.B - L.R  1.49646782 0.4097436 34   3.652  0.0024



data_SG_biomass_trans_surv_trt_soil_root_precip_g=data_SG_biomass_trans_surv_trt %>% group_by(soil_root,precip)
trans_shoot_soil_root_precip=summarise_at(data_SG_biomass_trans_surv_trt_soil_root_precip_g, 
                                    "shoot_weight_g_0", funs(n(),mean,sd,se=sd(.)/sqrt(n())))

ggplot(trans_shoot_soil_root_precip, aes(x=precip,y=mean, ymin = mean-se, ymax= mean+se,fill=soil_root))+
  geom_bar(position="dodge",stat="identity", color="black",aes(fill=soil_root))+
  geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+
  scale_fill_manual(values = c("gray", "white", "black"),labels=c("Bulk","Sterile","Rhizo"))+
  scale_x_discrete(labels=c("Ambient","Drought"))+ylab("Shoot biomass transplants")+
  geom_text(aes(y=mean+se+0.005, label=n),position=position_dodge(width=0.9))+theme_bw()


#filter out the dead transs

data_SG_biomass_trans_surv_trt_na=subset(data_SG_biomass_trans_surv_trt,total_biomass!=0)
summary(data_SG_biomass_trans_surv_trt_na)

#now lets look at total biomass



trans_biomas_model_na= lm(log(total_biomass_0)~precip*soil_root+(1|block), data= data_SG_biomass_trans_surv_trt_na)
qqPlot(resid(trans_biomas_model_na))
hist(resid(trans_biomas_model_na))

anova(trans_biomas_model_na, type=3)
#soil_root         27.89  2  11.2867 5.642e-05 ***
#precip:soil_root   8.41  2   3.4050    0.0388 *  

emmeans(trans_biomas_model_na, pairwise~soil_root)
emmeans(trans_biomas_model_na, pairwise~soil_root|precip)
#precip = A:
#contrast    estimate        SE df t.ratio p.value
#L.B - S.B -1.9058971 0.4058580 70  -4.696  <.0001
#L.B - L.R  0.1350397 0.4058580 70   0.333  0.9409
#S.B - L.R  2.0409369 0.4058580 70   5.029  <.0001

#precip = D:
#  contrast    estimate        SE df t.ratio p.value
#L.B - S.B -0.6390409 0.4574494 70  -1.397  0.3480
#L.B - L.R -0.1156812 0.4621815 70  -0.250  0.9661
#S.B - L.R  0.5233597 0.4785246 70   1.094  0.5211



data_SG_biomass_trans_surv_trt_na_soil_root_precip_g=data_SG_biomass_trans_surv_trt_na %>% group_by(soil_root,precip)
trans_total_bio_soil_root_precip_na=summarise_at(data_SG_biomass_trans_surv_trt_na_soil_root_precip_g, 
                                                 "total_biomass_0", funs(n(),mean,sd,se=sd(.)/sqrt(n())))

ggplot(trans_total_bio_soil_root_precip_na, aes(x=precip,y=mean, ymin = mean-se, ymax= mean+se,fill=soil_root))+
  geom_bar(position="dodge",stat="identity", color="black",aes(fill=soil_root))+
  geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+
  scale_fill_manual(values = c("gray", "white", "black"),labels=c("Bulk","Sterile","Rhizo"))+
  scale_x_discrete(labels=c("Ambient","Drought"))+ylab("Total biomass transplants (na removed)")+
  geom_text(aes(y=mean+se+0.02, label=n),position=position_dodge(width=0.9))+theme_bw()



#Roots alone
summary(data_SG_biomass_trans_surv_trt_na)
trans_root_model_na= lm(sqrt(root_weight_g_0)~precip*soil_root+(1|block), data= data_SG_biomass_trans_surv_trt_na)
qqPlot(resid(trans_root_model_na))
hist(resid(trans_root_model_na))

anova(trans_root_model_na, type=3)
#soil_root        0.5631  2   8.1806 0.0006416 ***
#precip:soil_root 0.2938  2   4.2685 0.0178163 *  


emmeans(trans_root_model_na, pairwise~soil_root)
emmeans(trans_root_model_na, pairwise~soil_root|precip)
#precip = A:
#contrast      estimate         SE df t.ratio p.value
#L.B - S.B -0.305804122 0.06774070 70  -4.514  0.0001
#L.B - L.R  0.013737110 0.06774070 70   0.203  0.9776
#S.B - L.R  0.319541232 0.06774070 70   4.717  <.0001

#precip = D:
#contrast      estimate         SE df t.ratio p.value
#L.B - S.B -0.050628188 0.07635168 70  -0.663  0.7856
#L.B - L.R -0.000403386 0.07714151 70  -0.005  1.0000
#S.B - L.R  0.050224802 0.07986929 70   0.629  0.8048



data_SG_biomass_trans_surv_trt_na_soil_root_precip_g=data_SG_biomass_trans_surv_trt_na %>% group_by(soil_root,precip)
trans_root_soil_root_precip_na=summarise_at(data_SG_biomass_trans_surv_trt_na_soil_root_precip_g, 
                                            "root_weight_g_0", funs(n(),mean,sd,se=sd(.)/sqrt(n())))

ggplot(trans_root_soil_root_precip_na, aes(x=precip,y=mean, ymin = mean-se, ymax= mean+se,fill=soil_root))+
  geom_bar(position="dodge",stat="identity", color="black",aes(fill=soil_root))+
  geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+
  scale_fill_manual(values = c("gray", "white", "black"),labels=c("Bulk","Sterile","Rhizo"))+
  scale_x_discrete(labels=c("Ambient","Drought"))+ylab("Root biomass transplants (na removed)")+
  geom_text(aes(y=mean+se+0.02, label=n),position=position_dodge(width=0.9))+theme_bw()



#Shoots alone
summary(data_SG_biomass_trans_surv_trt_na)
trans_shoot_model_na= lm(sqrt(shoot_weight_g_0)~precip*soil_root+(1|block), data= data_SG_biomass_trans_surv_trt_na)
qqPlot(resid(trans_shoot_model_na))
hist(resid(trans_shoot_model_na))

anova(trans_shoot_model_na, type=3)
#precip           0.03748  1   46.7229 2.489e-09 ***
#soil_root        0.16183  2  100.8740 < 2.2e-16 ***
#precip:soil_root 0.03762  2   23.4478 1.605e-08 ***

emmeans(trans_shoot_model_na, ~precip)
emmeans(trans_shoot_model_na, pairwise~soil_root)
emmeans(trans_shoot_model_na, pairwise~soil_root|precip)
#precip = A:
# L.B - S.B  0.14946329 0.2537568 34   0.589  0.8269
#L.B - L.R -0.25671123 0.1956448 34  -1.312  0.3983
#S.B - L.R -0.40617452 0.2432545 34  -1.670  0.2314

#precip = D:
#  contrast     estimate        SE df t.ratio p.value
#L.B - S.B -1.41470888 0.3800696 34  -3.722  0.0020
#L.B - L.R  0.08175893 0.2424079 34   0.337  0.9393
#S.B - L.R  1.49646782 0.4097436 34   3.652  0.0024



data_SG_biomass_trans_surv_trt_na_soil_root_precip_g=data_SG_biomass_trans_surv_trt_na %>% group_by(soil_root,precip)
trans_shoot_soil_root_precip_na=summarise_at(data_SG_biomass_trans_surv_trt_na_soil_root_precip_g, 
                                             "shoot_weight_g_0", funs(n(),mean,sd,se=sd(.)/sqrt(n())))

ggplot(trans_shoot_soil_root_precip_na, aes(x=precip,y=mean, ymin = mean-se, ymax= mean+se,fill=soil_root))+
  geom_bar(position="dodge",stat="identity", color="black",aes(fill=soil_root))+
  geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+
  scale_fill_manual(values = c("gray", "white", "black"),labels=c("Bulk","Sterile","Rhizo"))+
  scale_x_discrete(labels=c("Ambient","Drought"))+ylab("Shoot biomass transplants (na removed)")+
  geom_text(aes(y=mean+se+0.005, label=n),position=position_dodge(width=0.9))+theme_bw()


#Root:shoot ratio
data_SG_biomass_trans_surv_trt_na=subset(data_SG_biomass_trans_surv_trt,total_biomass!=0)
summary(data_SG_biomass_trans_surv_trt_na)
data_SG_biomass_trans_surv_trt_na$root_shoot=data_SG_biomass_trans_surv_trt_na$root_weight_g/data_SG_biomass_trans_surv_trt_na$shoot_weight_g


summary(data_SG_biomass_trans_surv_trt_na)

#####Diagnostic plots for probablity disturbutions####

#https://ase.tufts.edu/gsc/gradresources/guidetomixedmodelsinr/mixed%20model%20guide.html
require(MASS)



# lnorm means lognormal
qqp(data_SG_biomass_trans_surv_trt_na$root_shoot, "lnorm")

# qqp requires estimates of the parameters of the negative binomial, Poisson
# and gamma distributions. You can generate estimates using the fitdistr
# function. Save the output and extract the estimates of each parameter as I
# have shown below.
nbinom <- fitdistr(data_SG_biomass_trans_surv_trt_na$root_shoot, "Negative Binomial")
qqp(data_SG_biomass_trans_surv_trt_na$root_shoot, "nbinom", size = nbinom$estimate[[1]], mu = nbinom$estimate[[2]])

poisson <- fitdistr(data_SG_biomass_trans_surv_trt_na$root_shoot, "Poisson")
qqp(data_SG_biomass_trans_surv_trt_na$root_shoot, "pois", lambda=poisson$estimate)

gamma <- fitdistr(data_SG_biomass_trans_surv_trt_na$root_shoot, "gamma")
qqp(data_SG_biomass_trans_surv_trt_na$root_shoot, "gamma", shape = gamma$estimate[[1]], rate = gamma$estimate[[2]])


root_shoot_model_na= lm(log(root_shoot)~precip*soil_root+(1|block), data= data_SG_biomass_trans_surv_trt_na)
qqPlot(resid(root_shoot_model_na))
hist(resid(root_shoot_model_na))
shapiro.test(resid(root_shoot_model_na))
#0.0005919
boxCox(root_shoot_model_na)

anova(root_shoot_model_na, type=3)
#nada sig

emmeans(root_shoot_model_na, ~precip)
emmeans(root_shoot_model_na, pairwise~soil_root)
emmeans(root_shoot_model_na, pairwise~soil_root|precip)
#precip = A:
# L.B - S.B  0.14946329 0.2537568 34   0.589  0.8269
#L.B - L.R -0.25671123 0.1956448 34  -1.312  0.3983
#S.B - L.R -0.40617452 0.2432545 34  -1.670  0.2314

#precip = D:
#  contrast     estimate        SE df t.ratio p.value
#L.B - S.B -1.41470888 0.3800696 34  -3.722  0.0020
#L.B - L.R  0.08175893 0.2424079 34   0.337  0.9393
#S.B - L.R  1.49646782 0.4097436 34   3.652  0.0024

root_shoot_model_na_no_block= lm(log(root_shoot)~precip*soil_root, data= data_SG_biomass_trans_surv_trt_na)
qqPlot(resid(root_shoot_model_na_no_block))
hist(resid(root_shoot_model_na_no_block))
shapiro.test(resid(root_shoot_model_na_no_block))
#0.001287
boxCox(root_shoot_model_na_no_block)

anova(root_shoot_model_na_no_block, type=3)
#nada sig

emmeans(root_shoot_model_na_no_block, pairwise~soil_root|precip)

AIC(root_shoot_model_na,root_shoot_model_na_no_block)

data_SG_biomass_trans_surv_trt_na %>% group_by(soil_root,precip) %>% summarise_at("root_shoot", funs(n(),mean,sd,se=sd(.)/sqrt(n())))



data_SG_biomass_trans_surv_trt_na_soil_root_precip_g=data_SG_biomass_trans_surv_trt_na %>% group_by(soil_root,precip)
trans_root_shoot_soil_root_precip_na=summarise_at(data_SG_biomass_trans_surv_trt_na_soil_root_precip_g, 
                                                "root_shoot", funs(n(),mean,sd,se=sd(.)/sqrt(n())))

ggplot(trans_root_shoot_soil_root_precip_na, aes(x=precip,y=mean, ymin = mean-se, ymax= mean+se,fill=soil_root))+
  geom_bar(position="dodge",stat="identity", color="black",aes(fill=soil_root))+
  geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+
  scale_fill_manual(values = c("gray", "white", "black"),labels=c("Bulk","Sterile","Rhizo"))+
  scale_x_discrete(labels=c("Ambient","Drought"))+ylab("Root:shoot (na removed)")+
  geom_text(aes(y=mean+se+1, label=n),position=position_dodge(width=0.9))+theme_bw()


#####Trans root:shoot reformated graph####

(trans_root_t_shoot_p=ggplot(trans_root_shoot_soil_root_precip_na, aes(x=precip,y=mean, ymin = mean-se, ymax= mean+se,fill=factor(soil_root,levels = treatment_order)))+
   geom_bar(position="dodge",stat="identity", color="black",aes(fill=factor(soil_root,levels = treatment_order)))+
   geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+ylim(c(0,13))+
   scale_fill_manual(values = c( "white","lightgray", "darkgrey"),labels=c("Sterile","Bulk","Rhizo"), limits=treatment_order, name=NULL)+
   scale_x_discrete(labels=c("Ambient","Drought"), name=NULL)+ylab("Root:shoot")+
   geom_text(aes(y=1, label=n),size=10,position=position_dodge(width=0.9))+theme_bw()+
   theme(axis.title = element_text(size = 23), axis.text = element_text(size = 23),
         legend.position = c(0.85,.9), legend.text=element_text(size=20),
         legend.background = element_rect(size=0.5,linetype="solid",colour ="black")))


(trans_root_t_shoot_p2=ggplot(trans_root_shoot_soil_root_precip_na, aes(x=precip,y=mean, ymin = mean-se, ymax= mean+se,fill=factor(soil_root,levels = treatment_order)))+
    geom_bar(position="dodge",stat="identity", color="black",aes(fill=factor(soil_root,levels = treatment_order)))+
    geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+ylim(c(0,14))+
    scale_fill_manual(values = c( "white","lightgray", "darkgrey"),labels=c("Sterile","Bulk","Rhizo"), limits=treatment_order, name=NULL)+
    scale_x_discrete(labels=c("Ambient","Drought"), name=NULL)+ylab("Root:shoot")+
    geom_text(aes(y=1, label=n),size=10,position=position_dodge(width=0.9))+theme_bw()+
    theme(axis.title = element_text(size = 23), axis.text = element_text(size = 23),
          legend.position = "none"))



#####Combined root shoot graph#####
ggarrange(trans_root_t_shoot_p,seed_root_t_shoot_p,ncol = 2,  legend = "none")
#15x7.38


#remove the outlier
data_SG_biomass_trans_surv_trt_na_out=subset(data_SG_biomass_trans_surv_trt_na, root_shoot<50)
summary(data_SG_biomass_trans_surv_trt_na_out)

root_shoot_model_na_out= lm(sqrt(root_shoot)~precip*soil_root+(1|block), data= data_SG_biomass_trans_surv_trt_na_out)
qqPlot(resid(root_shoot_model_na_out))
hist(resid(root_shoot_model_na_out))

anova(root_shoot_model_na_out, type=3)
#precip             4.33  1   3.0869 0.08336 . 

emmeans(root_shoot_model_na_out, ~precip)
emmeans(root_shoot_model_na_out, pairwise~soil_root)
emmeans(root_shoot_model_na_out, pairwise~soil_root|precip)




data_SG_biomass_trans_surv_trt_na_soil_root_precip_out_g=data_SG_biomass_trans_surv_trt_na_out %>% group_by(soil_root,precip)
trans_root_shoot_soil_root_precip_na_out=summarise_at(data_SG_biomass_trans_surv_trt_na_soil_root_precip_out_g, 
                                                     "root_shoot", funs(n(),mean,sd,se=sd(.)/sqrt(n())))

ggplot(trans_root_shoot_soil_root_precip_na_out, aes(x=precip,y=mean, ymin = mean-se, ymax= mean+se,fill=soil_root))+
  geom_bar(position="dodge",stat="identity", color="black",aes(fill=soil_root))+
  geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+
  scale_fill_manual(values = c("gray", "white", "black"),labels=c("Bulk","Sterile","Rhizo"))+
  scale_x_discrete(labels=c("Ambient","Drought"))+ylab("Root:shoot (na removed)")+
  geom_text(aes(y=mean+se+1, label=n),position=position_dodge(width=0.9))+theme_bw()



#####Presence Root:shoot ratio####
data_SG_biomass_trans_surv_pres_na=subset(data_SG_biomass_trans_surv_trt_na, root_association =="B")

root_shoot_model_pres_na= lm((root_shoot)^(1/4)~precip*soil_root+(1|block), data= data_SG_biomass_trans_surv_pres_na)
qqPlot(resid(root_shoot_model_pres_na))
hist(resid(root_shoot_model_pres_na))
shapiro.test(resid(root_shoot_model_pres_na))
#0.0005919
boxCox(root_shoot_model_pres_na)

anova(root_shoot_model_pres_na, type=3)
#nada sig

emmeans(root_shoot_model_pres_na, ~precip)
emmeans(root_shoot_model_pres_na, pairwise~soil_root)
emmeans(root_shoot_model_pres_na, pairwise~soil_root|precip)


root_shoot_model_pres_na_no_block= lm((root_shoot)^(1/4)~precip*soil_root, data= data_SG_biomass_trans_surv_pres_na)
qqPlot(resid(root_shoot_model_pres_na_no_block))
hist(resid(root_shoot_model_pres_na_no_block))
shapiro.test(resid(root_shoot_model_pres_na_no_block))
#0.01335
boxCox(root_shoot_model_pres_na_no_block)

anova(root_shoot_model_pres_na_no_block, type=3)
#nada sig

emmeans(root_shoot_model_pres_na_no_block, pairwise~soil_root|precip)

AIC(root_shoot_model_pres_na,root_shoot_model_pres_na_no_block)

#####Origin Root:shoot ratio####
data_SG_biomass_trans_surv_orig_na=subset(data_SG_biomass_trans_surv_trt_na,soil_status =="L")

root_shoot_model_orig_na= lm(log(root_shoot)~precip*soil_root+(1|block), data= data_SG_biomass_trans_surv_orig_na)
qqPlot(resid(root_shoot_model_orig_na))
hist(resid(root_shoot_model_orig_na))
shapiro.test(resid(root_shoot_model_orig_na))
#0.03771
boxCox(root_shoot_model_orig_na)

anova(root_shoot_model_orig_na, type=3)
#nada sig

emmeans(root_shoot_model_orig_na, ~precip)
emmeans(root_shoot_model_orig_na, pairwise~soil_root)
emmeans(root_shoot_model_orig_na, pairwise~soil_root|precip)


root_shoot_model_orig_na_no_block= lm(log(root_shoot)~precip*soil_root, data= data_SG_biomass_trans_surv_orig_na)
qqPlot(resid(root_shoot_model_orig_na_no_block))
hist(resid(root_shoot_model_orig_na_no_block))
shapiro.test(resid(root_shoot_model_orig_na_no_block))
#0.02169
boxCox(root_shoot_model_orig_na_no_block)

anova(root_shoot_model_orig_na_no_block, type=3)
#nada sig

emmeans(root_shoot_model_orig_na_no_block, pairwise~soil_root|precip)

AIC(root_shoot_model_orig_na,root_shoot_model_orig_na_no_block)

#Final survival


summary(data_SG_biomass_trans_surv_trt)


data_SG_biomass_trans_surv_trt$surv=data_SG_biomass_trans_surv_trt$total_biomass
data_SG_biomass_trans_surv_trt$surv[is.na(data_SG_biomass_trans_surv_trt$surv)]=0
data_SG_biomass_trans_surv_trt$surv[data_SG_biomass_trans_surv_trt$surv!=0]=1

data_SG_biomass_trans_surv_trt %>% group_by(precip) %>% summarise_at("surv", funs(n(),mean))
(0.778-1)/1

data_SG_biomass_trans_surv_trt %>% group_by(precip, soil_root) %>% summarise_at("surv", funs(n(),mean))


#####Trans survival Reformated bar graph #####
data_SG_biomass_trans_surv_trt_soil_root_precip_out_g=data_SG_biomass_trans_surv_trt %>% group_by(soil_root,precip)
trans_survial_soil_root_precip=summarise_at(data_SG_biomass_trans_surv_trt_soil_root_precip_out_g, 
                                                      "surv", funs(n(),mean))


(trans_surv_p=ggplot(trans_survial_soil_root_precip, aes(x=precip,y=mean,fill=factor(soil_root,levels = treatment_order)))+
   geom_bar(position="dodge",stat="identity", color="black",aes(fill=factor(soil_root,levels = treatment_order)))+
   scale_fill_manual(values = c( "white","light gray", "dark grey"),labels=c("Sterile","Bulk","Rhizo"), limits=treatment_order, name=NULL)+
   scale_x_discrete(labels=c("Ambient","Drought"), name=NULL)+ylab("prop Transplant survival)")+ylim(c(0,1))+
   geom_text(aes(y=0.07, label=n),size=10,position=position_dodge(width=0.9))+theme_bw()+
   theme(axis.title = element_text(size = 23), axis.text = element_text(size = 23),
         legend.position = c(0.85,.9), legend.text=element_text(size=20),
         legend.background = element_rect(size=0.5,linetype="solid",colour ="black")))

trans_survial_soil_root_precip$mean_per=trans_survial_soil_root_precip$mean*100
(trans_surv_per_p=ggplot(trans_survial_soil_root_precip, aes(x=precip,y=mean_per,fill=factor(soil_root,levels = treatment_order)))+
    geom_bar(position="dodge",stat="identity", color="black",aes(fill=factor(soil_root,levels = treatment_order)))+
    scale_fill_manual(values = c( "white","lightgray", "darkgrey"),labels=c("Sterile","Bulk","Rhizo"), limits=treatment_order, name=NULL)+
    scale_x_discrete(labels=c("Ambient","Drought"), name=NULL)+ylab("% seedling survival")+scale_y_continuous(limits = c(0,120), breaks = c(0,20,40,60,80,100))+
    geom_text(aes(y=7, label=n),size=10,position=position_dodge(width=0.9))+theme_bw()+
    theme(axis.title = element_text(size = 23), axis.text = element_text(size = 23),
          legend.position = "none"))


#####Combined survial graph#####
ggarrange(trans_surv_p,seedling_suv_p,ncol = 2,  legend = "none")
#15x7.38
#Final survival in the drought treatment only since all pots had live transs in the ambient treatments

data_SG_biomass_trans_surv_drought=subset(data_SG_biomass_trans_surv_trt, precip=="D")

#survival of the transplant

drought_surv_model= glm(surv~soil_root+(1|block), data= data_SG_biomass_trans_surv_drought, family = binomial)
qqPlot(resid(drought_surv_model))
hist(resid(drought_surv_model))

anova(drought_surv_model, type=3)


emmeans(drought_surv_model, pairwise~soil_root)

drought_surv_model_no_block= glm(surv~soil_root, data= data_SG_biomass_trans_surv_drought, family = binomial)
qqPlot(resid(drought_surv_model_no_block))
hist(resid(drought_surv_model_no_block))

anova(drought_surv_model_no_block, type=3)
#soil_root   1.0984  2     0.5774

AIC(drought_surv_model,drought_surv_model_no_block)
#> AIC(drought_surv_model,drought_surv_model_no_block)
#                           df      AIC
#drought_surv_model           7 49.93455
#drought_surv_model_no_block  3 52.57514
#Bars by soil_root

data_SG_biomass_trans_surv_drought_g=group_by(data_SG_biomass_trans_surv_drought, soil_root)
drought_surv_soil_root=summarise_at(data_SG_biomass_trans_surv_drought_g, "surv", mean)

ggplot(drought_surv_soil_root, aes(soil_root,surv))+geom_bar(stat="identity", color="black",aes(fill=soil_root))+
  scale_fill_manual(values = c("gray", "white", "black"),labels=c("Bulk","Sterile","Rhizo"))+
  scale_x_discrete(limits=c("S.B", "L.B", "L.R"),labels=c("Sterile","Bulk","Rhizo"))+ylab("Survival under drought treatment")



####Presence survival of the transplant####

data_SG_biomass_trans_surv_drought_pres=subset(data_SG_biomass_trans_surv_drought, root_association =="B")
unique(data_SG_biomass_trans_surv_drought_pres$soil_root)

drought_surv_model_pres= glm(surv~soil_root+(1|block), data= data_SG_biomass_trans_surv_drought_pres, family = binomial)
qqPlot(resid(drought_surv_model_pres))
hist(resid(drought_surv_model_pres))

anova(drought_surv_model_pres, type=3)


emmeans(drought_surv_model_pres, pairwise~soil_root)

drought_surv_model_pres_no_block= glm(surv~soil_root, data= data_SG_biomass_trans_surv_drought_pres, family = binomial)
qqPlot(resid(drought_surv_model_pres_no_block))
hist(resid(drought_surv_model_pres_no_block))

anova(drought_surv_model_pres_no_block, type=3)
#soil_root  0.84646  1     0.3576

AIC(drought_surv_model_pres,drought_surv_model_pres_no_block)

####Origin survival of the transplant####

data_SG_biomass_trans_surv_drought_orig=subset(data_SG_biomass_trans_surv_drought, soil_status =="L")
unique(data_SG_biomass_trans_surv_drought_orig$soil_root)
drought_surv_model_orig= glm(surv~soil_root+(1|block), data= data_SG_biomass_trans_surv_drought_orig, family = binomial)
qqPlot(resid(drought_surv_model_orig))
hist(resid(drought_surv_model_orig))

anova(drought_surv_model_orig, type=3)
#(1|block)   8.9261  4    0.06297 .

emmeans(drought_surv_model_orig, pairwise~soil_root)

drought_surv_model_orig_no_block= glm(surv~soil_root, data= data_SG_biomass_trans_surv_drought_orig, family = binomial)
qqPlot(resid(drought_surv_model_orig_no_block))
hist(resid(drought_surv_model_orig_no_block))

anova(drought_surv_model_orig_no_block, type=3)
#soil_root  0.84646  1     0.3576

AIC(drought_surv_model_orig,drought_surv_model_orig_no_block)


#####Change in height####

#let's look at height from intial planting to the end of the experiment

summary(SG_height_combin)

SG_height_combin_planting=subset(SG_height_combin, exp_period=="planting")

summary(SG_height_combin_planting)


#remove any plants that had to be replaced in first week
replacement_no=as.numeric(c(SG_height_combin$Plant_Number[SG_height_combin$exp_period=="replacement"]))
summary(replacement_no)
SG_height_combin_int_H<-SG_height_combin_planting[SG_height_combin_planting$Plant_Number %w/o% replacement_no,]
nrow(SG_height_combin_int_H)

SG_height_combin_int_height<-rbind(SG_height_combin_int_H, SG_height_combin[SG_height_combin$exp_period=="replacement",])
nrow(SG_height_combin_int_height)
anyDuplicated(SG_height_combin_int_height$Plant_Number)
SG_height_harvest=subset(SG_height_combin,exp_period=="harvest")
nrow(SG_height_harvest)

SG_height_plant_harvest=merge(SG_height_combin_int_height,SG_height_harvest, by="Plant_Number")
nrow(SG_height_plant_harvest)
SG_height_plant_harvest$int_height=SG_height_plant_harvest$height_cm.x
SG_height_plant_harvest$int_height[is.na(SG_height_plant_harvest$height_cm.x)]=0
SG_height_plant_harvest$har_height=SG_height_plant_harvest$height_cm.y
SG_height_plant_harvest$har_height[is.na(SG_height_plant_harvest$height_cm.y)]=0
colnames(SG_height_plant_harvest)[19]="surv"
nrow(SG_height_plant_harvest)

SG_height_plant_harvest$delta_height=SG_height_plant_harvest$har_height-SG_height_plant_harvest$int_height
nrow(SG_height_plant_harvest)
SG_height_plant_harvest_trt=merge(SG_height_plant_harvest, SG_trt, by="Plant_Number", all.x = T)
nrow(SG_height_plant_harvest_trt)
SG_height_plant_harvest_trt$soil_root=with(SG_height_plant_harvest_trt, interaction(soil_status,root_association))
nrow(SG_height_plant_harvest_trt)
#Let's look at the transplants
transplant_SG_height_plant_harvest_trt=subset(SG_height_plant_harvest_trt,life_stage=="G")


summary(SG_height_plant_harvest_trt)

delta_height_trans_model= lm((delta_height)~precip*soil_root+(1|block), data= transplant_SG_height_plant_harvest_trt)
qqPlot(resid(delta_height_trans_model))
hist(resid(delta_height_trans_model))

anova(delta_height_trans_model, type=3)
#precip            770.3  1  64.4430 7.140e-12 ***
#soil_root        1074.7  2  44.9547 8.220e-14 ***
#(1|block)  151.0  4   3.1587    0.0183 *  
#precip:soil_root  463.1  2  19.3710 1.379e-07 ***

emmeans(delta_height_trans_model, ~precip)
emmeans(delta_height_trans_model, pairwise~soil_root)
emmeans(delta_height_trans_model, pairwise~soil_root|precip)
#precip = A:
#contrast  estimate   SE  df t.ratio p.value
#L.B - S.B   -4.410 1.31 170 -3.357  0.0028 
#L.B - L.R   -0.587 1.31 170 -0.447  0.8960 
#S.B - L.R    3.823 1.31 170  2.910  0.0114 

#precip = D:
#  contrast  estimate   SE  df t.ratio p.value
#L.B - S.B   -0.827 1.31 170 -0.629  0.8042 
#L.B - L.R    1.253 1.31 170  0.954  0.6069 
#S.B - L.R    2.080 1.31 170  1.583  0.2556 

emmeans(delta_height_trans_model, pairwise~soil_root*precip, adjust="fdr")
"$contrasts
 contrast      estimate   SE df t.ratio p.value
 L.B,A - S.B,A  -11.467 1.26 80 -9.083  <.0001 
 L.B,A - L.R,A    1.213 1.26 80  0.961  0.4242 
 L.B,A - L.B,D    2.700 1.26 80  2.139  0.0666 
 L.B,A - S.B,D    0.800 1.26 80  0.634  0.5658 
 L.B,A - L.R,D    3.800 1.26 80  3.010  0.0087 
 S.B,A - L.R,A   12.680 1.26 80 10.044  <.0001 
 S.B,A - L.B,D   14.167 1.26 80 11.222  <.0001 
 S.B,A - S.B,D   12.267 1.26 80  9.717  <.0001 
 S.B,A - L.R,D   15.267 1.26 80 12.093  <.0001 
 L.R,A - L.B,D    1.487 1.26 80  1.178  0.3306 
 L.R,A - S.B,D   -0.413 1.26 80 -0.327  0.7442 
 L.R,A - L.R,D    2.587 1.26 80  2.049  0.0729 
 L.B,D - S.B,D   -1.900 1.26 80 -1.505  0.2044 
 L.B,D - L.R,D    1.100 1.26 80  0.871  0.4456 
 S.B,D - L.R,D    3.000 1.26 80  2.376  0.0426 

Results are averaged over the levels of: block 
P value adjustment: fdr method for 15 tests "
transplant_SG_height_plant_harvest_trt_soil_root_precip_g=transplant_SG_height_plant_harvest_trt %>% group_by(soil_root,precip)
delta_height_trans_soil_root_precip=summarise_at(transplant_SG_height_plant_harvest_trt_soil_root_precip_g, 
                                                "delta_height", funs(n(),mean,sd,se=sd(.)/sqrt(n())))

ggplot(delta_height_trans_soil_root_precip, aes(x=precip,y=mean, ymin = mean-se, ymax= mean+se,fill=soil_root))+
  geom_bar(position="dodge",stat="identity", color="black",aes(fill=soil_root))+
  geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+
  scale_fill_manual(values = c("gray", "white", "black"),labels=c("Bulk","Sterile","Rhizo"))+
  scale_x_discrete(labels=c("Ambient","Drought"))+ylab("Change in height")+
  geom_text(aes(y=mean+se+0.5, label=n),position=position_dodge(width=0.9))+theme_bw()


#with dead plants removed

summary(transplant_SG_height_plant_harvest_trt)

transplant_SG_height_plant_harvest_trt_0=subset(transplant_SG_height_plant_harvest_trt, surv==1)


summary(transplant_SG_height_plant_harvest_trt_0)

delta_height_trans_model_0= lm((delta_height)~precip*soil_root+(1|block), data= transplant_SG_height_plant_harvest_trt_0)
qqPlot(resid(delta_height_trans_model_0))
hist(resid(delta_height_trans_model_0))

anova(delta_height_trans_model_0, type=3)
#precip            339.9  1  47.7230 2.211e-09 ***
#soil_root        1210.5  2  84.9733 < 2.2e-16 ***
#(1|block)   66.5  4   2.3352   0.06437 .  
#precip:soil_root  185.3  2  13.0044 1.690e-05 ***

emmeans(delta_height_trans_model_0, ~precip)
emmeans(delta_height_trans_model_0, pairwise~soil_root)
emmeans(delta_height_trans_model_0, pairwise~soil_root|precip)

#precip = A:
#contrast  estimate    SE df t.ratio p.value
#L.B - S.B   -11.17 0.993 67 -11.252 <.0001 
#L.B - L.R     1.51 0.993 67   1.519 0.2887 
#S.B - L.R    12.68 0.975 67  13.011 <.0001 

#precip = D:
#  contrast  estimate    SE df t.ratio p.value
#L.B - S.B    -4.60 1.145 67  -4.016 0.0004 
#L.B - L.R     1.15 1.134 67   1.018 0.5680 
#S.B - L.R     5.75 1.179 67   4.881 <.0001 



transplant_SG_height_plant_harvest_trt_0_soil_root_precip_g=transplant_SG_height_plant_harvest_trt_0 %>% group_by(soil_root,precip)
delta_height_trans_0_soil_root_precip=summarise_at(transplant_SG_height_plant_harvest_trt_0_soil_root_precip_g, 
                                                 "delta_height", funs(n(),mean,sd,se=sd(.)/sqrt(n())))

ggplot(delta_height_trans_0_soil_root_precip, aes(x=precip,y=mean, ymin = mean-se, ymax= mean+se,fill=soil_root))+
  geom_bar(position="dodge",stat="identity", color="black",aes(fill=soil_root))+
  geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+
  scale_fill_manual(values = c("gray", "white", "black"),labels=c("Bulk","Sterile","Rhizo"))+
  scale_x_discrete(labels=c("Ambient","Drought"))+ylab("Change in height (with dead plants removed)")+
  geom_text(aes(y=mean+se+0.5, label=n),position=position_dodge(width=0.9))+theme_bw()



#pot weight of transplants pots

dataSG_potweight_trt_trans=subset(dataSG_potweight_trt, life_stage=="G")
dataSG_potweight_trt_trans$soil_root=with(dataSG_potweight_trt_trans, interaction(soil_status,root_association))
dataSG_potweight_trt_trans_sub=subset(dataSG_potweight_trt_trans, date!="6_13")
min(dataSG_potweight_trt_trans_sub$per_ch_potweight)
delta_pot_weight_model= lmer((per_ch_potweight+0.234105)^2~precip*soil_root*date+(1|block)+(1|Plant_Number), data= dataSG_potweight_trt_trans_sub)
plot(delta_pot_weight_model)
qqPlot(resid(delta_pot_weight_model))
hist(resid(delta_pot_weight_model))

anova(delta_pot_weight_model, type=3)
#precip:date             84.4040  3  < 2.2e-16 ***
#soil_root:date          38.6782  6  8.276e-07 ***
#precip:soil_root:date   16.1963  6    0.01274 *  
#precip                 106.2803  1  < 2.2e-16 ***
#soil_root              107.5332  2  < 2.2e-16 ***
#date                   801.5201  3  < 2.2e-16 ***
emmeans(delta_pot_weight_model, ~precip)
emmeans(delta_pot_weight_model, pairwise~soil_root)

emmeans(delta_pot_weight_model, pairwise~soil_root|date|precip)

dataSG_delta_potweight_trt_soil_root_time_g=dataSG_potweight_trt_trans_sub %>% group_by(date,soil_root)
pot_delta_potweight_soil_root_time=summarise_at(dataSG_delta_potweight_trt_soil_root_time_g, "per_ch_potweight", mean)
#dot graph
ggplot(pot_delta_potweight_soil_root_time, aes(x=date,y=per_ch_potweight))+geom_point(size=4,aes(color=soil_root))+
  scale_color_manual(values = c("gray", "white", "black"),labels=c("Bulk","Sterile","Rhizo"))+ylab("percentage change in pot weight from int")+
  xlab("date")
#boxplot
ggplot(dataSG_potweight_trt_trans_sub, aes(x=date,y=per_ch_potweight))+geom_boxplot(aes(fill=soil_root))+
  scale_fill_manual(values = c("gray", "white", "black"),labels=c("Bulk","Sterile","Rhizo"))+ylab("percentage change in pot weight from int")+
  xlab("date")
#dot graph
dataSG_delta_potweight_trt_soil_root_time_drought_g=dataSG_potweight_trt_trans_sub %>% group_by(date,soil_root,precip)
pot_delta_potweight_soil_root_time_drought=summarise_at(dataSG_delta_potweight_trt_soil_root_time_drought_g, "per_ch_potweight", mean)
#boxplot
ggplot(pot_delta_potweight_soil_root_time_drought, aes(x=date,y=per_ch_potweight))+geom_point(shape=21,size=4,aes(color=precip,fill=soil_root))+
  ylab("percentage change in pot weight from int")+scale_fill_manual(values = c("gray", "white", "black"),labels=c("Bulk","Sterile","Rhizo"))+
  scale_color_manual(values = c("blue", "red"),labels=c("Ambient","Drought"))+ xlab("Date")


ggplot(dataSG_potweight_trt_trans_sub, aes(x=date,y=per_ch_potweight))+geom_boxplot(aes(color=precip,fill=soil_root))+
  ylab("percentage change in pot weight from int")+scale_fill_manual(values = c("gray", "white", "black"),labels=c("Bulk","Sterile","Rhizo"))+
  scale_color_manual(values = c("blue", "red"),labels=c("Ambient","Drought"))+ xlab("Date")


#remove the min and max value is a next step here




#####VWC as a covariate#####
dataSG_potweight_trt_trans=subset(dataSG_potweight_trt, life_stage=="G")
nrow(dataSG_potweight_trt_trans)
#450
summary(dataSG_potweight_trt_trans)
dataSG_potweight_trt_trans$soil_root=with(dataSG_potweight_trt_trans, interaction(soil_status,root_association))
dataSG_potweight_trt_trans_sub=subset(dataSG_potweight_trt_trans, date!="6_13")
summary(dataSG_potweight_trt_trans_sub)
#Extremely large outliers in bulk treatment


#I am going to use the end change in pot weight
dataSG_potweight_trt_trans_end=subset(dataSG_potweight_trt_trans_sub,date=="7_7")
nrow(dataSG_potweight_trt_trans_end)
#90
summary(dataSG_potweight_trt_trans_end)
dataSG_potweight_trans_end=dataSG_potweight_trt_trans_end[,c("Plant_Number","per_ch_potweight","per_ch_potweight_bio_crt")]
summary(dataSG_potweight_trans_end)
colnames(dataSG_potweight_trans_end)[2]="end_per_ch_potweight"
colnames(dataSG_potweight_trans_end)[3]="end_per_ch_potweight_bio_crt"
summary(dataSG_potweight_trans_end)


summary(dataSG_trans_surv_trt)
nrow(dataSG_trans_surv_trt)
dataSG_trans_surv_trt_pw=merge(dataSG_trans_surv_trt,dataSG_potweight_trans_end, by="Plant_Number")
nrow(dataSG_trans_surv_trt_pw)
dataSG_trans_surv_trt_pw_g=group_by(dataSG_trans_surv_trt_pw, Plant_Number)
total_transs_transplanted_pw=summarise_at(dataSG_trans_surv_trt_pw_g, "Transplant_height_cm", sum,na.rm = TRUE)
colnames(total_transs_transplanted_pw)[2]="tot_Transplant_height_cm"
summary(total_transs_transplanted_pw)
nrow(total_transs_transplanted_pw)

fin_dataSG_trans_surv_trt_pw=subset(dataSG_trans_surv_trt_pw, date=="7_9")
nrow(fin_dataSG_trans_surv_trt_pw)
fin_dataSG_trans_surv_trt_pw_merge=merge(fin_dataSG_trans_surv_trt_pw,total_transs_transplanted_pw, by="Plant_Number")
nrow(fin_dataSG_trans_surv_trt_pw_merge)
fin_dataSG_trans_surv_trt_pw_merge$int_surv=fin_dataSG_trans_surv_trt_pw_merge$tot_Transplant_height_cm
fin_dataSG_trans_surv_trt_pw_merge$int_surv[fin_dataSG_trans_surv_trt_pw_merge$int_surv==0]=1
fin_dataSG_trans_surv_trt_pw_merge$int_surv[fin_dataSG_trans_surv_trt_pw_merge$int_surv!=1]=0
fin_dataSG_trans_surv_trt_pw_merge$soil_root=with(fin_dataSG_trans_surv_trt_pw_merge, interaction(soil_status,root_association))


#now lets look at biomass with dead as zero
summary(data_SG_biomass)
data_SG_biomass_trans_surv_trt_pw=merge(fin_dataSG_trans_surv_trt_pw, data_SG_biomass, by="Plant_Number", all.x = T)
data_SG_biomass_trans_surv_trt_pw$soil_root=with(data_SG_biomass_trans_surv_trt_pw, interaction(soil_status,root_association))
summary(data_SG_biomass_trans_surv_trt_pw)


data_SG_biomass_trans_surv_trt_pw$total_biomass[data_SG_biomass_trans_surv_trt_pw$total_biomass==0]=NA
data_SG_biomass_trans_surv_trt_pw$total_biomass_0=data_SG_biomass_trans_surv_trt_pw$total_biomass
data_SG_biomass_trans_surv_trt_pw$total_biomass_0[is.na(data_SG_biomass_trans_surv_trt_pw$total_biomass_0)]=0

data_SG_biomass_trans_surv_trt_pw$root_weight_g_0=data_SG_biomass_trans_surv_trt_pw$root_weight_g
data_SG_biomass_trans_surv_trt_pw$root_weight_g_0[is.na(data_SG_biomass_trans_surv_trt_pw$root_weight_g_0)]=0

data_SG_biomass_trans_surv_trt_pw$shoot_weight_g_0=data_SG_biomass_trans_surv_trt_pw$shoot_weight_g
data_SG_biomass_trans_surv_trt_pw$shoot_weight_g_0[is.na(data_SG_biomass_trans_surv_trt_pw$shoot_weight_g_0)]=0

summary(data_SG_biomass_trans_surv_trt_pw)


#now lets look at total biomass



trans_biomas_model_pw= lm(log(total_biomass_0+0.001)~end_per_ch_potweight_bio_crt*soil_root+(1|block), data= data_SG_biomass_trans_surv_trt_pw)
qqPlot(resid(trans_biomas_model_pw))
hist(resid(trans_biomas_model_pw))
shapiro.test(resid(trans_biomas_model_pw))
#0.02558
anova(trans_biomas_model_pw, type=3)
#nada sig

emmeans(trans_biomas_model_pw, pairwise~soil_root)
emmeans(trans_biomas_model_pw, pairwise~soil_root|precip)
#precip = A:
#contrast     estimate        SE df t.ratio p.value
#L.B - S.B -1.86324726 0.5633964 80  -3.307  0.0040
#L.B - L.R  0.14030501 0.5633964 80   0.249  0.9664
#S.B - L.R  2.00355226 0.5633964 80   3.556  0.0018

#precip = D:
#  contrast     estimate        SE df t.ratio p.value
#L.B - S.B  0.04493185 0.5633964 80   0.080  0.9965
#L.B - L.R  0.42719232 0.5633964 80   0.758  0.7295
#S.B - L.R  0.38226047 0.5633964 80   0.678  0.7767

emmeans(trans_biomas_model_pw, pairwise~soil_root*precip, adjust="fdr")
"$contrasts
contrast      estimate    SE df t.ratio p.value
L.B,A - S.B,A  -1.8632 0.563 80 -3.307  0.0042 
L.B,A - L.R,A   0.1403 0.563 80  0.249  0.8614 
L.B,A - L.B,D   0.4599 0.563 80  0.816  0.6493 
L.B,A - S.B,D   0.5048 0.563 80  0.896  0.6493 
L.B,A - L.R,D   0.8871 0.563 80  1.575  0.2983 
S.B,A - L.R,A   2.0036 0.563 80  3.556  0.0024 
S.B,A - L.B,D   2.3232 0.563 80  4.123  0.0005 
S.B,A - S.B,D   2.3681 0.563 80  4.203  0.0005 
S.B,A - L.R,D   2.7504 0.563 80  4.882  0.0001 
L.R,A - L.B,D   0.3196 0.563 80  0.567  0.6601 
L.R,A - S.B,D   0.3645 0.563 80  0.647  0.6493 
L.R,A - L.R,D   0.7468 0.563 80  1.326  0.4045 
L.B,D - S.B,D   0.0449 0.563 80  0.080  0.9366 
L.B,D - L.R,D   0.4272 0.563 80  0.758  0.6493 
S.B,D - L.R,D   0.3823 0.563 80  0.678  0.6493 

Results are averaged over the levels of: block 
Results are given on the log (not the response) scale. 
P value adjustment: fdr method for 15 tests "
data_SG_biomass_trans_surv_trt_soil_root_precip_g=data_SG_biomass_trans_surv_trt %>% group_by(soil_root,precip)
trans_total_bio_soil_root_precip=summarise_at(data_SG_biomass_trans_surv_trt_soil_root_precip_g, 
                                        "total_biomass_0", funs(n(),mean,sd,se=sd(.)/sqrt(n())))

ggplot(trans_total_bio_soil_root_precip, aes(x=precip,y=mean, ymin = mean-se, ymax= mean+se,fill=soil_root))+
  geom_bar(position="dodge",stat="identity", color="black",aes(fill=soil_root))+
  geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+
  scale_fill_manual(values = c("gray", "white", "black"),labels=c("Bulk","Sterile","Rhizo"))+
  scale_x_discrete(labels=c("Ambient","Drought"))+ylab("Total biomass transplants")+
  geom_text(aes(y=mean+se+0.02, label=n),position=position_dodge(width=0.9))+theme_bw()


#let's look at the root:shoot ratio
data_SG_biomass_trans_surv_trt_pw$root_shoot=data_SG_biomass_trans_surv_trt_pw$root_weight_g/data_SG_biomass_trans_surv_trt_pw$shoot_weight_g
summary(data_SG_biomass_trans_surv_trt_pw)


trans_root_shoot_model_pw= lm(log(root_shoot)~end_per_ch_potweight_bio_crt*soil_root+(1|block), data= data_SG_biomass_trans_surv_trt_pw)
qqPlot(resid(trans_root_shoot_model_pw))
hist(resid(trans_root_shoot_model_pw))
boxCox(trans_root_shoot_model_pw)
shapiro.test(resid(trans_root_shoot_model_pw))
#0.000549
anova(trans_root_shoot_model_pw, type=3)
#nada sig




#let's look at height from intial planting to the end of the experiment

summary(SG_height_combin)

SG_height_combin_planting=subset(SG_height_combin, exp_period=="planting")

summary(SG_height_combin_planting)


#remove any plants that had to be replaced in first week
replacement_no=as.numeric(c(SG_height_combin$Plant_Number[SG_height_combin$exp_period=="replacement"]))
summary(replacement_no)
SG_height_combin_int_H<-SG_height_combin_planting[SG_height_combin_planting$Plant_Number %w/o% replacement_no,]
nrow(SG_height_combin_int_H)

SG_height_combin_int_height<-rbind(SG_height_combin_int_H, SG_height_combin[SG_height_combin$exp_period=="replacement",])
nrow(SG_height_combin_int_height)
anyDuplicated(SG_height_combin_int_height$Plant_Number)
SG_height_harvest=subset(SG_height_combin,exp_period=="harvest")
nrow(SG_height_harvest)

SG_height_plant_harvest=merge(SG_height_combin_int_height,SG_height_harvest, by="Plant_Number")
summary(SG_height_plant_harvest)
nrow(SG_height_plant_harvest)
SG_height_plant_harvest$int_height=SG_height_plant_harvest$height_cm.x
SG_height_plant_harvest$int_height[is.na(SG_height_plant_harvest$height_cm.x)]=0
SG_height_plant_harvest$har_height=SG_height_plant_harvest$height_cm.y
SG_height_plant_harvest$har_height[is.na(SG_height_plant_harvest$height_cm.y)]=0
colnames(SG_height_plant_harvest)[19]="surv"
nrow(SG_height_plant_harvest)

SG_height_plant_harvest$delta_height=SG_height_plant_harvest$har_height-SG_height_plant_harvest$int_height
nrow(SG_height_plant_harvest)
summary(SG_height_plant_harvest)
SG_height_plant_harvest_trt=merge(SG_height_plant_harvest, SG_trt, by="Plant_Number", all.x = T)
nrow(SG_height_plant_harvest_trt)
dataSG_potweight_trans_end=dataSG_potweight_trt_trans_end[,c("Plant_Number","per_ch_potweight","per_ch_potweight_bio_crt")]
colnames(dataSG_potweight_trans_end)[2]="end_per_ch_potweight"
colnames(dataSG_potweight_trans_end)[3]="end_per_ch_potweight_bio_crt"

summary(dataSG_potweight_trans_end)
nrow(dataSG_potweight_trans_end)
SG_height_plant_harvest_trt_pw=merge(SG_height_plant_harvest_trt,dataSG_potweight_trans_end, by="Plant_Number")
nrow(SG_height_plant_harvest_trt_pw)
summary(SG_height_plant_harvest_trt_pw)

SG_height_plant_harvest_trt_pw$soil_root=with(SG_height_plant_harvest_trt_pw, interaction(soil_status,root_association))
nrow(SG_height_plant_harvest_trt_pw)
#Let's look at the transplants
transplant_SG_height_plant_harvest_trt_pw=subset(SG_height_plant_harvest_trt_pw,life_stage=="G")


summary(transplant_SG_height_plant_harvest_trt_pw)
nrow(transplant_SG_height_plant_harvest_trt_pw)
delta_height_trans_model_pw= lm(delta_height~end_per_ch_potweight_bio_crt*soil_root+(1|block), data= transplant_SG_height_plant_harvest_trt_pw)
qqPlot(resid(delta_height_trans_model_pw))
hist(resid(delta_height_trans_model_pw))
boxCox(delta_height_trans_model_pw)
shapiro.test(resid(delta_height_trans_model_pw))
#0.0005005


anova(delta_height_trans_model_pw, type=3)
#end_per_ch_potweight_bio_crt            118.45  1  4.5840   0.03532 *  


emmeans(delta_height_trans_model_pw, pairwise~soil_root)


ggplot(transplant_SG_height_plant_harvest_trt_pw, aes(x=end_per_ch_potweight_bio_crt,y=delta_height))+geom_point(aes(fill=soil_root), size=4, shape=21)+
  geom_smooth(method='lm')+theme_bw()+scale_fill_manual(values = c("gray", "white", "black"),labels=c("Bulk","Sterile","Rhizo"))


transplant_SG_height_plant_harvest_trt_soil_root_precip_g=transplant_SG_height_plant_harvest_trt %>% group_by(soil_root,precip)
delta_height_trans_soil_root_precip=summarise_at(transplant_SG_height_plant_harvest_trt_soil_root_precip_g, 
                                                 "delta_height", funs(n(),mean,sd,se=sd(.)/sqrt(n())))

ggplot(delta_height_trans_soil_root_precip, aes(x=precip,y=mean, ymin = mean-se, ymax= mean+se,fill=soil_root))+
  geom_bar(position="dodge",stat="identity", color="black",aes(fill=soil_root))+
  geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+
  scale_fill_manual(values = c("gray", "white", "black"),labels=c("Bulk","Sterile","Rhizo"))+
  scale_x_discrete(labels=c("Ambient","Drought"))+ylab("Change in height")+
  geom_text(aes(y=mean+se+0.5, label=n),position=position_dodge(width=0.9))+theme_bw()


#with dead plants removed

summary(transplant_SG_height_plant_harvest_trt_pw)

transplant_SG_height_plant_harvest_trt_pw_0=subset(transplant_SG_height_plant_harvest_trt_pw, surv==1)


summary(transplant_SG_height_plant_harvest_trt_pw_0)

delta_height_trans_model_0_pw= lm((delta_height)~end_per_ch_potweight_bio_crt*soil_root+(1|block), data= transplant_SG_height_plant_harvest_trt_pw_0)
qqPlot(resid(delta_height_trans_model_0_pw))
hist(resid(delta_height_trans_model_0_pw))
#plot(delta_height_trans_model_0_pw)
boxCox(delta_height_trans_model_0_pw)

anova(delta_height_trans_model_0_pw, type=3)
#end_per_ch_potweight_bio_crt            101.96  1  7.8684  0.006578 ** 
#soil_root                               191.91  2  7.4048  0.001243 ** 


emmeans(delta_height_trans_model_0, pairwise~soil_root)



ggplot(transplant_SG_height_plant_harvest_trt_pw_0, aes(x=end_per_ch_potweight_bio_crt,y=delta_height))+geom_point(aes(fill=soil_root), size=4, shape=21)+
  geom_smooth(method='lm')+theme_bw()+scale_fill_manual(values = c("gray", "white", "black"),labels=c("Bulk","Sterile","Rhizo"))


#####Root traits####

summary(SG_roottraits_trt)
SG_roottraits_trt$soil_root=with(SG_roottraits_trt, interaction(soil_status,root_association))
#Specific root length
SRL_trans_model= lm((SRL_length_drymass)~precip*soil_root+(1|block), data= SG_roottraits_trt)
qqPlot(resid(SRL_trans_model))
hist(resid(SRL_trans_model))
boxCox(SRL_trans_model)
#plot(SRL_trans_model_pw)


anova(SRL_trans_model, type=3)
#precip            140431454  1   6.5482  0.02033 *  
#precip:soil_root  200016455  2   4.6633  0.02429 *  

emmeans(SRL_trans_model, pairwise~soil_root|precip)

SRL_trans_model_no_block= lm((SRL_length_drymass)~precip*soil_root, data= SG_roottraits_trt)
qqPlot(resid(SRL_trans_model_no_block))
hist(resid(SRL_trans_model_no_block))
boxCox(SRL_trans_model_no_block)
shapiro.test(resid(SRL_trans_model_no_block))
#0.6666
#plot(SRL_trans_model_pw)


anova(SRL_trans_model_no_block, type=3)
#precip            112669528  1   4.8980   0.03810 *  
#precip:soil_root  166002900  2   3.6082   0.04498 *  

emmeans(SRL_trans_model_no_block, pairwise~soil_root|precip)

AIC(SRL_trans_model,SRL_trans_model_no_block)

SG_roottraits_trt %>% group_by(precip) %>% summarise_at("SRL_length_drymass", funs(n(),mean,sd,se=sd(.)/sqrt(n())))

(15522-12067)/12067


SG_roottraits_trt %>% group_by(precip,soil_root) %>% summarise_at("SRL_length_drymass", funs(n(),mean,sd,se=sd(.)/sqrt(n())))

treatment_order=c("S.B","L.B","L.R")
SG_roottraits_trt_soil_root_precip_g=SG_roottraits_trt %>% group_by(soil_root,precip)
SRL_trans_soil_root_precip=summarise_at(SG_roottraits_trt_soil_root_precip_g, 
                                        "SRL_length_drymass", funs(n(),mean,sd,se=sd(.)/sqrt(n())))
#Live Bulk
(12476-20949)/20949
#Live Rhizosphere
(15260-20949)/20949

(SRL_trans_p=ggplot(SRL_trans_soil_root_precip, aes(x=precip,y=mean, ymin = mean-se, ymax= mean+se,fill=factor(soil_root,levels = treatment_order)))+
  geom_bar(position="dodge",stat="identity", color="black",aes(fill=factor(soil_root,levels = treatment_order)))+
  geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+ylim(c(0,27500))+
  scale_fill_manual(values = c( "white","light gray", "dark grey"),labels=c("Sterile","Bulk","Rhizo"), limits=treatment_order, name=NULL)+
  scale_x_discrete(labels=c("Ambient","Drought"), name=NULL)+ylab("Specific root length (cm/g)")+
  geom_text(aes(y=2000, label=n),size=10,position=position_dodge(width=0.9))+theme_bw()+
  theme(axis.title = element_text(size = 23), axis.text = element_text(size = 23),
        legend.position = c(0.85,.9), legend.text=element_text(size=20),
        legend.background = element_rect(size=0.5,linetype="solid",colour ="black")))

(SRL_trans_p2=ggplot(SRL_trans_soil_root_precip, aes(x=precip,y=mean, ymin = mean-se, ymax= mean+se,fill=factor(soil_root,levels = treatment_order)))+
    geom_bar(position="dodge",stat="identity", color="black",aes(fill=factor(soil_root,levels = treatment_order)))+
    geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+ylim(c(0,29000))+
    scale_fill_manual(values = c( "white","lightgray", "darkgrey"),labels=c("Sterile","Bulk","Rhizo"), limits=treatment_order, name=NULL)+
    scale_x_discrete(labels=c("Ambient","Drought"), name=NULL)+ylab("Specific root length (cm/g)")+
    geom_text(aes(y=2000, label=n),size=10,position=position_dodge(width=0.9))+theme_bw()+
    theme(axis.title = element_text(size = 23), axis.text = element_text(size = 23),
          legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank()))  

####Presence specific root length####

SG_roottraits_pres=subset(SG_roottraits_trt, root_association =="B")
unique(SG_roottraits_pres$soil_root)

SRL_trans_model_pres= lm((SRL_length_drymass)~precip*soil_root+(1|block), data= SG_roottraits_pres)
qqPlot(resid(SRL_trans_model_pres))
hist(resid(SRL_trans_model_pres))
boxCox(SRL_trans_model_pres)
#plot(SRL_trans_model_pw)


anova(SRL_trans_model_pres, type=3)
#precip            140431454  1   6.5482  0.02033 *  
#precip:soil_root  200016455  2   4.6633  0.02429 *  

emmeans(SRL_trans_model_pres, pairwise~soil_root|precip)

emmeans(SRL_trans_model_pres, pairwise~soil_root*precip,adjust="none")

SRL_trans_model_pres_no_block= lm((SRL_length_drymass)~precip*soil_root, data= SG_roottraits_pres)
qqPlot(resid(SRL_trans_model_pres_no_block))
hist(resid(SRL_trans_model_pres_no_block))
boxCox(SRL_trans_model_pres_no_block)
shapiro.test(resid(SRL_trans_model_pres_no_block))
#0.6666
#plot(SRL_trans_model_pw)


anova(SRL_trans_model_pres_no_block, type=3)
#precip            112669528  1   4.8980   0.03810 *  
#precip:soil_root  166002900  2   3.6082   0.04498 *  

emmeans(SRL_trans_model_pres_no_block, pairwise~soil_root|precip)

AIC(SRL_trans_model_pres,SRL_trans_model_pres_no_block)


####Origin specific root length####

SG_roottraits_orig=subset(SG_roottraits_trt,soil_status =="L")
unique(SG_roottraits_orig$soil_root)

SRL_trans_model_orig= lm((SRL_length_drymass)~precip*soil_root+(1|block), data= SG_roottraits_orig)
qqPlot(resid(SRL_trans_model_orig))
hist(resid(SRL_trans_model_orig))
boxCox(SRL_trans_model_orig)
#plot(SRL_trans_model_pw)


anova(SRL_trans_model_orig, type=3)
#nada sig 

emmeans(SRL_trans_model_orig, pairwise~soil_root|precip)

emmeans(SRL_trans_model_orig, pairwise~soil_root*precip,adjust="none")

SRL_trans_model_orig_no_block= lm((SRL_length_drymass)~precip*soil_root, data= SG_roottraits_orig)
qqPlot(resid(SRL_trans_model_orig_no_block))
hist(resid(SRL_trans_model_orig_no_block))
boxCox(SRL_trans_model_orig_no_block)
shapiro.test(resid(SRL_trans_model_orig_no_block))
#0.6444
#plot(SRL_trans_model_pw)


anova(SRL_trans_model_orig_no_block, type=3)
#nada sig 

emmeans(SRL_trans_model_orig_no_block, pairwise~soil_root|precip)

AIC(SRL_trans_model_orig,SRL_trans_model_orig_no_block)


#root length density 
RLD_trans_model= lm((RLD_length_volume)~precip*soil_root+(1|block), data= SG_roottraits_trt)
qqPlot(resid(RLD_trans_model))
hist(resid(RLD_trans_model))
#plot(RLD_trans_model_pw)
boxCox(RLD_trans_model)

anova(RLD_trans_model, type=3)
#precip            0.3138  1   4.0265  0.060968 .  
#soil_root         6.2268  2  39.9434 3.763e-07 ***
#precip:soil_root  1.2581  2   8.0705  0.003433 ** 

emmeans(RLD_trans_model, pairwise~precip|soil_root)
emmeans(RLD_trans_model, pairwise~soil_root|precip)

RLD_trans_model_no_block= lm((RLD_length_volume)~precip*soil_root, data= SG_roottraits_trt)
qqPlot(resid(RLD_trans_model_no_block))
hist(resid(RLD_trans_model_no_block))
#plot(RLD_trans_model_pw)
boxCox(RLD_trans_model_no_block)
shapiro.test(resid(RLD_trans_model_no_block))
#0.5027

anova(RLD_trans_model_no_block, type=3)
#precip            0.4261  1   5.0492  0.035515 *  
#soil_root         6.3855  2  37.8319 1.092e-07 ***
#precip:soil_root  1.2362  2   7.3238  0.003864 **

emmeans(RLD_trans_model_no_block, pairwise~soil_root|precip)

AIC(RLD_trans_model,RLD_trans_model_no_block)

SG_roottraits_trt %>% group_by(precip) %>% summarise_at("RLD_length_volume", funs(n(),mean,sd,se=sd(.)/sqrt(n())))
(0.719-1.04)/1.04

SG_roottraits_trt %>% group_by(precip, soil_root) %>% summarise_at("RLD_length_volume", funs(n(),mean,sd,se=sd(.)/sqrt(n())))

treatment_order=c("S.B","L.B","L.R")
SG_roottraits_trt_soil_root_precip_g=SG_roottraits_trt %>% group_by(soil_root,precip)
RLD_trans_soil_root_precip=summarise_at(SG_roottraits_trt_soil_root_precip_g, 
                                                 "RLD_length_volume", funs(n(),mean,sd,se=sd(.)/sqrt(n())))

ggplot(RLD_trans_soil_root_precip, aes(x=precip,y=mean, ymin = mean-se, ymax= mean+se,fill=factor(soil_root,levels = treatment_order)))+
  geom_bar(position="dodge",stat="identity", color="black",aes(fill=factor(soil_root,levels = treatment_order)))+
  geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+
  scale_fill_manual(values = c( "white","light gray", "dark grey"),labels=c("Sterile","Bulk","Rhizo"), limits=treatment_order, name=NULL)+
  scale_x_discrete(labels=c("Ambient","Drought"), name=NULL)+ylab("Root length density (cm/mL)")+
  geom_text(aes(y=.2, label=n),position=position_dodge(width=0.9))+theme_bw()+theme(axis.title = element_text(size = 23), axis.text = element_text(size = 23),
                                                                                    legend.position = c(0.85,.9), legend.text=element_text(size=20),
                                                                                    legend.background = element_rect(size=0.5,linetype="solid",colour ="black"))

(RLD_trans_p=ggplot(RLD_trans_soil_root_precip, aes(x=precip,y=mean, ymin = mean-se, ymax= mean+se,fill=factor(soil_root,levels = treatment_order)))+
    geom_bar(position="dodge",stat="identity", color="black",aes(fill=factor(soil_root,levels = treatment_order)))+
    geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+ylim(c(0,2.5))+
    scale_fill_manual(values = c( "white","lightgray", "darkgrey"),labels=c("Sterile","Bulk","Rhizo"), limits=treatment_order, name=NULL)+
    scale_x_discrete(labels=c("Ambient","Drought"), name=NULL)+ylab("Root length density (cm/mL)")+
    geom_text(aes(y=0.2, label=n),size=10,position=position_dodge(width=0.9))+theme_bw()+
    theme(axis.title = element_text(size = 23), axis.text = element_text(size = 23),
          legend.position = c(0.85,.9), legend.text=element_text(size=20),
          legend.background = element_rect(size=0.5,linetype="solid",colour ="black")))

(RLD_trans_p2=ggplot(RLD_trans_soil_root_precip, aes(x=precip,y=mean, ymin = mean-se, ymax= mean+se,fill=factor(soil_root,levels = treatment_order)))+
    geom_bar(position="dodge",stat="identity", color="black",aes(fill=factor(soil_root,levels = treatment_order)))+
    geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+ylim(c(0,3))+
    scale_fill_manual(values = c( "white","lightgray", "darkgrey"),labels=c("Sterile","Bulk","Rhizo"), limits=treatment_order, name=NULL)+
    scale_x_discrete(labels=c("Ambient","Drought"), name=NULL)+ylab("Root length density (cm/mL)")+
    geom_text(aes(y=0.2, label=n),size=10,position=position_dodge(width=0.9))+theme_bw()+
    theme(axis.title = element_text(size = 23), axis.text = element_text(size = 23),
          legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank()))
#800x750

(legend_RLD_trans_p2=ggplot(RLD_trans_soil_root_precip, aes(x=precip,y=mean, ymin = mean-se, ymax= mean+se,fill=factor(soil_root,levels = treatment_order)))+
    geom_bar(position="dodge",stat="identity", color="black",aes(fill=factor(soil_root,levels = treatment_order)))+
    geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+ylim(c(0,3))+
    scale_fill_manual(values = c( "white","lightgray", "darkgrey"),labels=c("Sterile","Bulk","Rhizo"), limits=treatment_order, name=NULL)+
    scale_x_discrete(labels=c("Ambient","Drought"), name=NULL)+ylab("Root length density (cm/mL)")+
    geom_text(aes(y=0.2, label=n),size=10,position=position_dodge(width=0.9))+theme_bw()+
    theme(axis.title = element_text(size = 23), axis.text = element_text(size = 23),
          legend.position = c(0.85,.9), legend.text=element_text(size=20),
          legend.background = element_rect(size=0.5,linetype="solid",colour ="black")
          ,panel.grid.major = element_blank(), panel.grid.minor = element_blank()))


#####Combined root traits graph#####
ggarrange(SRL_trans_p,RLD_trans_p,ncol = 2,  legend = "none")
#15x7.38

####Presence root length density####
SG_roottraits_pres=subset(SG_roottraits_trt, root_association =="B")
unique(SG_roottraits_pres$soil_root)


RLD_trans_model_pres= lm((RLD_length_volume)~precip*soil_root+(1|block), data= SG_roottraits_pres)
qqPlot(resid(RLD_trans_model_pres))
hist(resid(RLD_trans_model_pres))
#plot(RLD_trans_model_pw)
boxCox(RLD_trans_model_pres)

anova(RLD_trans_model_pres, type=3)
#precip            0.9260  1  13.2637  0.004522 ** 
#soil_root         5.1649  1  73.9829 6.207e-06 ***
#precip:soil_root  0.6550  1   9.3830  0.011980 *  

emmeans(RLD_trans_model_pres, pairwise~precip|soil_root)
emmeans(RLD_trans_model_pres, pairwise~soil_root|precip)
emmeans(RLD_trans_model_pres, pairwise~soil_root*precip,adjust="none")



RLD_trans_model_pres_no_block= lm((RLD_length_volume)~precip*soil_root, data= SG_roottraits_pres)
qqPlot(resid(RLD_trans_model_pres_no_block))
hist(resid(RLD_trans_model_pres_no_block))
#plot(RLD_trans_model_pw)
boxCox(RLD_trans_model_pres_no_block)
shapiro.test(resid(RLD_trans_model_pres_no_block))
#0.4136

anova(RLD_trans_model_pres_no_block, type=3)
#precip            0.9774  1  10.9665  0.005143 ** 
#soil_root         5.3936  1  60.5195 1.893e-06 ***
#precip:soil_root  0.6924  1   7.7689  0.014540 *  

emmeans(RLD_trans_model_no_block, pairwise~soil_root|precip)

AIC(RLD_trans_model_pres,RLD_trans_model_pres_no_block)

####Origin root length density####
SG_roottraits_orig=subset(SG_roottraits_trt,soil_status =="L")
unique(SG_roottraits_orig$soil_root)


RLD_trans_model_orig= lm((RLD_length_volume)~precip*soil_root+(1|block), data= SG_roottraits_orig)
qqPlot(resid(RLD_trans_model_orig))
hist(resid(RLD_trans_model_orig))
#plot(RLD_trans_model_pw)
boxCox(RLD_trans_model_orig)

anova(RLD_trans_model_orig, type=3)
#nada sig

emmeans(RLD_trans_model_orig, pairwise~precip|soil_root)
emmeans(RLD_trans_model_orig, pairwise~soil_root|precip)
emmeans(RLD_trans_model_orig, pairwise~soil_root*precip,adjust="none")



RLD_trans_model_orig_no_block= lm((RLD_length_volume)~precip*soil_root, data= SG_roottraits_orig)
qqPlot(resid(RLD_trans_model_orig_no_block))
hist(resid(RLD_trans_model_orig_no_block))
#plot(RLD_trans_model_pw)
boxCox(RLD_trans_model_orig_no_block)
shapiro.test(resid(RLD_trans_model_orig_no_block))
#0.6776

anova(RLD_trans_model_orig_no_block, type=3)
#nada sig 

emmeans(RLD_trans_model_orig_no_block, pairwise~soil_root|precip)
emmeans(RLD_trans_model_orig_no_block, pairwise~soil_root*precip,adjust="none")

AIC(RLD_trans_model_orig,RLD_trans_model_orig_no_block)


#Rhizosheath RhizosheathSoil_DryRoots_g 
min(SG_roottraits_trt$RhizosheathSoil_DryRoots_g)
rhizosheath_model= lm(log(RhizosheathSoil_DryRoots_g +1)~precip*soil_root+(1|block), data= SG_roottraits_trt)
qqPlot(resid(rhizosheath_model))
hist(resid(rhizosheath_model))
shapiro.test(resid(rhizosheath_model))
#p-value = 0.5554
#plot(rhizosheath_model_pw)
boxCox(rhizosheath_model)

anova(rhizosheath_model, type=3)
#nada sig


rhizosheath_model_no_block= lm(log(RhizosheathSoil_DryRoots_g +1)~precip*soil_root, data= SG_roottraits_trt)
qqPlot(resid(rhizosheath_model_no_block))
hist(resid(rhizosheath_model_no_block))
shapiro.test(resid(rhizosheath_model_no_block))
#p-value = 0.6068
#plot(rhizosheath_model_pw)
boxCox(rhizosheath_model_no_block)

anova(rhizosheath_model_no_block, type=3)
#soil_root          5.148  2   2.5945   0.09841 .

emmeans(rhizosheath_model_no_block, pairwise~soil_root|precip)

AIC(rhizosheath_model,rhizosheath_model_no_block)

emmeans(rhizosheath_model_no_block, pairwise~soil_root)

SG_roottraits_trt %>% group_by(soil_root,precip) %>% summarise_at("RhizosheathSoil_DryRoots_g", funs(n(),mean,sd,se=sd(.)/sqrt(n())))

treatment_order=c("S.B","L.B","L.R")
SG_roottraits_trt_soil_root_precip_g=SG_roottraits_trt %>% group_by(soil_root,precip)
RhizoSheath_trans_soil_root_precip=summarise_at(SG_roottraits_trt_soil_root_precip_g, 
                                        "RhizosheathSoil_DryRoots_g", funs(n(),mean,sd,se=sd(.)/sqrt(n())))

(rhizosheath_trans_p=ggplot(RhizoSheath_trans_soil_root_precip, aes(x=precip,y=mean, ymin = mean-se, ymax= mean+se,fill=factor(soil_root,levels = treatment_order)))+
    geom_bar(position="dodge",stat="identity", color="black",aes(fill=factor(soil_root,levels = treatment_order)))+
    geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+
    scale_fill_manual(values = c( "white","lightgray", "darkgrey"),labels=c("Sterile","Bulk","Rhizo"), limits=treatment_order, name=NULL)+
    scale_x_discrete(labels=c("Ambient","Drought"), name=NULL)+ylab("Rhizosheath soil \nper dry roots (g/g)")+
    geom_text(aes(y=1.8, label=n),size=10,position=position_dodge(width=0.9))+theme_bw()+
    theme(axis.title = element_text(size = 23), axis.text = element_text(size = 23),
          legend.position = "none", legend.text=element_text(size=20),
          legend.background = element_rect(size=0.5,linetype="solid",colour ="black")))


(rhizosheath_trans_p2=ggplot(RhizoSheath_trans_soil_root_precip, aes(x=precip,y=mean, ymin = mean-se, ymax= mean+se,fill=factor(soil_root,levels = treatment_order)))+
    geom_bar(position="dodge",stat="identity", color="black",aes(fill=factor(soil_root,levels = treatment_order)))+
    geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+
    scale_fill_manual(values = c( "white","lightgray", "darkgrey"),labels=c("Sterile","Bulk","Rhizo"), limits=treatment_order, name=NULL)+
    scale_x_discrete(labels=c("Ambient","Drought"), name=NULL)+ylab("Rhizosheath soil \nper dry roots (g/g)")+
    geom_text(aes(y=1.8, label=n),size=10,position=position_dodge(width=0.9))+theme_bw()+ylim(0,35)+
    theme(axis.title = element_text(size = 23), axis.text = element_text(size = 23),
          legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank()))
#800x750
#7x5
#####Combined root traits graph#####
ggarrange(RLD_trans_p,rhizosheath_trans_p,ncol = 2,  legend = "none")
#15x7.38


####Presence Rhizosheath RhizosheathSoil_DryRoots_g####
SG_roottraits_pres=subset(SG_roottraits_trt, root_association =="B")
unique(SG_roottraits_pres$soil_root)


min(SG_roottraits_pres$RhizosheathSoil_DryRoots_g)
rhizosheath_model_pres= lm((RhizosheathSoil_DryRoots_g)~precip*soil_root+(1|block), data= SG_roottraits_pres)
qqPlot(resid(rhizosheath_model_pres))
hist(resid(rhizosheath_model_pres))
shapiro.test(resid(rhizosheath_model_pres))
#p-value = 0.9993
#plot(rhizosheath_model_pw)
boxCox(rhizosheath_model_pres)


anova(rhizosheath_model_pres, type=3)
#soil_root         429.4  1  4.5279  0.05923 .  
#(1|block) 1127.3  4  2.9720  0.07396 .  


rhizosheath_model_pres_no_block= lm((RhizosheathSoil_DryRoots_g)~precip*soil_root, data= SG_roottraits_pres)
qqPlot(resid(rhizosheath_model_pres_no_block))
hist(resid(rhizosheath_model_pres_no_block))
shapiro.test(resid(rhizosheath_model_pres_no_block))
#p-value = 0.6068
#plot(rhizosheath_model_pw)
boxCox(rhizosheath_model_pres_no_block)

anova(rhizosheath_model_pres_no_block, type=3)
#nada sig

emmeans(rhizosheath_model_pres_no_block, pairwise~soil_root|precip)

AIC(rhizosheath_model_pres,rhizosheath_model_pres_no_block)

####Origin Rhizosheath RhizosheathSoil_DryRoots_g####
SG_roottraits_orig=subset(SG_roottraits_trt,soil_status =="L")
unique(SG_roottraits_orig$soil_root)


min(SG_roottraits_pres$RhizosheathSoil_DryRoots_g)
rhizosheath_model_orig= lm(log(RhizosheathSoil_DryRoots_g+1)~precip*soil_root+(1|block), data= SG_roottraits_orig)
qqPlot(resid(rhizosheath_model_orig))
hist(resid(rhizosheath_model_orig))
shapiro.test(resid(rhizosheath_model_orig))
#p-value = 0.05013
#plot(rhizosheath_model_pw)
boxCox(rhizosheath_model_orig)


anova(rhizosheath_model_orig, type=3)
#nada sif


rhizosheath_model_orig_no_block= lm(log(RhizosheathSoil_DryRoots_g+1)~precip*soil_root, data= SG_roottraits_orig)
qqPlot(resid(rhizosheath_model_orig_no_block))
hist(resid(rhizosheath_model_orig_no_block))
shapiro.test(resid(rhizosheath_model_orig_no_block))
#p-value = 0.7577
#plot(rhizosheath_model_pw)
boxCox(rhizosheath_model_orig_no_block)

anova(rhizosheath_model_orig_no_block, type=3)
#nada sig

emmeans(rhizosheath_model_orig_no_block, pairwise~soil_root|precip)

AIC(rhizosheath_model_orig,rhizosheath_model_orig_no_block)



#Potweight
dataSG_potweight_trt_trans=subset(dataSG_potweight_trt, life_stage=="G")
nrow(dataSG_potweight_trt_trans)
#450
summary(dataSG_potweight_trt_trans)
dataSG_potweight_trt_trans$soil_root=with(dataSG_potweight_trt_trans, interaction(soil_status,root_association))
dataSG_potweight_trt_trans_sub=subset(dataSG_potweight_trt_trans, date!="6_13")
summary(dataSG_potweight_trt_trans_sub)
dataSG_potweight_trt_trans_end=subset(dataSG_potweight_trt_trans_sub,date=="7_7")
nrow(dataSG_potweight_trt_trans_end)
dataSG_potweight_trans_end=dataSG_potweight_trt_trans_end[,c("Plant_Number","per_ch_potweight","per_ch_potweight_bio_crt")]
colnames(dataSG_potweight_trans_end)[2]="end_per_ch_potweight"
colnames(dataSG_potweight_trans_end)[3]="end_per_ch_potweight_bio_crt"



SG_roottraits_trt_pw=merge(SG_roottraits_trt,dataSG_potweight_trans_end, by="Plant_Number")
nrow(SG_roottraits_trt_pw)
SG_roottraits_trt_pw$soil_root=with(SG_roottraits_trt_pw, interaction(soil_status,root_association))
head(SG_roottraits_trt_pw)

#Specific root length
SRL_trans_model_pw= lm((SRL_length_drymass)~end_per_ch_potweight_bio_crt*soil_root+(1|block), data= SG_roottraits_trt_pw)
qqPlot(resid(SRL_trans_model_pw))
hist(resid(SRL_trans_model_pw))
#plot(SRL_trans_model_pw)


anova(SRL_trans_model_pw, type=3)
#end_per_ch_potweight           141733838  1  5.2270 0.03534 *
summary(SRL_trans_model_pw)

SRL_trans_model_pw_reg= lm((SRL_length_drymass)~end_per_ch_potweight_bio_crt, data= SG_roottraits_trt_pw)
qqPlot(resid(SRL_trans_model_pw_reg))
hist(resid(SRL_trans_model_pw_reg))
#plot(SRL_trans_model_pw)


anova(SRL_trans_model_pw_reg, type=3)

summary(SRL_trans_model_pw_reg)

ggplot(SG_roottraits_trt_pw, aes(x=end_per_ch_potweight,y=SRL_length_drymass))+geom_point(aes(fill=soil_root), size=4, shape=21)+
  geom_smooth(method='lm')+theme_bw()+scale_fill_manual(values = c("gray", "white", "black"),labels=c("Bulk","Sterile","Rhizo"))

#Sterile
SG_roottraits_trt_pw_sterile=subset(SG_roottraits_trt_pw, soil_root=="S.B")
SRL_trans_model_pw_sterile= lm((SRL_length_drymass)~end_per_ch_potweight_bio_crt, data= SG_roottraits_trt_pw_sterile)
qqPlot(resid(SRL_trans_model_pw))
hist(resid(SRL_trans_model_pw))
#plot(SRL_trans_model_pw)


anova(SRL_trans_model_pw_sterile, type=3)
#end_per_ch_potweight           141733838  1  5.2270 0.03534 *
summary(SRL_trans_model_pw_sterile)


#Bulk live
SG_roottraits_trt_pw_Live_bulk=subset(SG_roottraits_trt_pw, soil_root=="L.B")
SRL_trans_model_pw_Live_bulk= lm((SRL_length_drymass)~end_per_ch_potweight_bio_crt, data= SG_roottraits_trt_pw_Live_bulk)
qqPlot(resid(SRL_trans_model_pw_Live_bulk))
hist(resid(SRL_trans_model_pw_Live_bulk))
#plot(SRL_trans_model_pw)


anova(SRL_trans_model_pw_Live_bulk, type=3)
#end_per_ch_potweight           141733838  1  5.2270 0.03534 *
summary(SRL_trans_model_pw_Live_bulk)


#Rhizosphere live
SG_roottraits_trt_pw_Live_rhizo=subset(SG_roottraits_trt_pw, soil_root=="L.R")
SRL_trans_model_pw_Live_rhizo= lm((SRL_length_drymass)~end_per_ch_potweight_bio_crt, data= SG_roottraits_trt_pw_Live_rhizo)
qqPlot(resid(SRL_trans_model_pw_Live_rhizo))
hist(resid(SRL_trans_model_pw_Live_rhizo))
#plot(SRL_trans_model_pw)


anova(SRL_trans_model_pw_Live_rhizo, type=3)
#end_per_ch_potweight           141733838  1  5.2270 0.03534 *
summary(SRL_trans_model_pw_Live_rhizo)


ggplot(SG_roottraits_trt_pw, aes(x=end_per_ch_potweight,y=SRL_length_drymass))+geom_point(aes(fill=soil_root), size=4, shape=21)+
  geom_smooth(method='lm', aes(linetype=soil_root))+theme_bw()+scale_fill_manual(values = c("gray", "white", "black"),labels=c("Bulk","Sterile","Rhizo"))


#root length density 
RLD_trans_model_pw= lm((RLD_length_volume)~end_per_ch_potweight_bio_crt*soil_root+(1|block), data= SG_roottraits_trt_pw)
qqPlot(resid(RLD_trans_model_pw))
hist(resid(RLD_trans_model_pw))
#plot(RLD_trans_model_pw)


anova(RLD_trans_model_pw, type=3)
#end_per_ch_potweight           0.49256  1  6.6694  0.019367 *  
#soil_root                      2.76879  2 18.7449 5.013e-05 ***
#end_per_ch_potweight:soil_root 1.12406  2  7.6100  0.004363 ** 


ggplot(SG_roottraits_trt_pw, aes(x=end_per_ch_potweight_bio_crt,y=RLD_length_volume))+geom_point(aes(fill=soil_root), size=4, shape=21)+
  geom_smooth(method='lm')+theme_bw()+scale_fill_manual(values = c("gray", "white", "black"),labels=c("Bulk","Sterile","Rhizo"))

ggplot(SG_roottraits_trt_pw, aes(x=end_per_ch_potweight_bio_crt,y=RLD_length_volume))+geom_point(aes(fill=soil_root), size=4, shape=21)+
  geom_smooth(method='lm',aes(linetype=soil_root))+theme_bw()+scale_fill_manual(values = c("gray", "white", "black"),labels=c("Bulk","Sterile","Rhizo"))


#Rhizosheath RhizosheathSoil_DryRoots_g 
rhizosheath_model_pw= lm((RhizosheathSoil_DryRoots_g )~end_per_ch_potweight_bio_crt*soil_root+(1|block), data= SG_roottraits_trt_pw)
qqPlot(resid(rhizosheath_model_pw))
hist(resid(rhizosheath_model_pw))
shapiro.test(resid(rhizosheath_model_pw))
#plot(rhizosheath_model_pw)


anova(rhizosheath_model_pw, type=3)
#end_per_ch_potweight           0.49256  1  6.6694  0.019367 *  
#soil_root                      2.76879  2 18.7449 5.013e-05 ***
#end_per_ch_potweight:soil_root 1.12406  2  7.6100  0.004363 ** 


ggplot(SG_roottraits_trt_pw, aes(x=end_per_ch_potweight_bio_crt,y=RhizosheathSoil_DryRoots_g))+geom_point(aes(fill=soil_root), size=4, shape=21)+
  geom_smooth(method='lm')+theme_bw()+scale_fill_manual(values = c("gray", "white", "black"),labels=c("Bulk","Sterile","Rhizo"))

ggplot(SG_roottraits_trt_pw, aes(x=end_per_ch_potweight,y=RLD_length_volume))+geom_point(aes(fill=soil_root), size=4, shape=21)+
  geom_smooth(method='lm',aes(linetype=soil_root))+theme_bw()+scale_fill_manual(values = c("gray", "white", "black"),labels=c("Bulk","Sterile","Rhizo"))


#soil moisture in block c which we used for the root traits


soil_moisture_root_traits_model= lm((end_per_ch_potweight_bio_crt)~precip*soil_root+(1|block), data= SG_roottraits_trt_pw)
qqPlot(resid(soil_moisture_root_traits_model))
hist(resid(soil_moisture_root_traits_model))
#plot(soil_moisture_root_traits_model)


anova(soil_moisture_root_traits_model, type=3)
emmeans(soil_moisture_root_traits_model, pairwise~soil_root|precip)

SG_roottraits_trt_pw_precip_soil_g=SG_roottraits_trt_pw %>% group_by(soil_root,precip)
soil_moist_soil_root_precip=summarise_at(SG_roottraits_trt_pw_precip_soil_g, 
                                     "end_per_ch_potweight_bio_crt", funs(n(),mean,sd,se=sd(.)/sqrt(n())))

ggplot(soil_moist_soil_root_precip, aes(x=precip,y=mean, ymin = mean-se, ymax= mean+se,fill=soil_root))+
  geom_bar(position="dodge",stat="identity", color="black",aes(fill=soil_root))+
  geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+
  scale_fill_manual(values = c("gray", "white", "black"),labels=c("Bulk","Sterile","Rhizo"))+
  scale_x_discrete(labels=c("Ambient","Drought"))+ylab("Soil moisture in block C")+
  geom_text(aes(y=mean-se-0.01, label=n),position=position_dodge(width=0.9))+theme_bw()





#####Soil nitrogen####

summary(SG_inorg_N_trt)
SG_inorg_N_trt$soil_root=with(SG_inorg_N_trt, interaction(soil_status,root_association))

SG_inorg_N_trans_trt= subset(SG_inorg_N_trt, life_stage=="G")
summary(SG_inorg_N_trans_trt)

#nitrate
SG_soil_nit_trans_model= lm((ug_N_NO3_g_dry_soil)~precip*soil_root+(1|block), data= SG_inorg_N_trans_trt)
qqPlot(resid(SG_soil_nit_trans_model))
hist(resid(SG_soil_nit_trans_model))

anova(SG_soil_nit_trans_model, type=3)
#precip             5.73  1   15.9870  0.001036 ** 
#soil_root          6.79  2    9.4836  0.001922 ** 
#precip:soil_root   3.57  2    4.9799  0.020824 *  
SG_inorg_N_trans_trt %>% group_by(precip) %>% summarise_at("ug_N_NO3_g_dry_soil", funs(n(),mean,sd,se=sd(.)/sqrt(n())))

(4.42-3.65)/3.65

SG_inorg_N_trans_trt %>% group_by(precip,soil_root) %>% summarise_at("ug_N_NO3_g_dry_soil", funs(n(),mean,sd,se=sd(.)/sqrt(n())))

emmeans(SG_soil_nit_trans_model, ~precip)
emmeans(SG_soil_nit_trans_model, pairwise~soil_root)
emmeans(SG_soil_nit_trans_model, pairwise~soil_root|precip)
#$contrasts
#precip = A:
#  contrast  estimate    SE df t.ratio p.value
#L.B - S.B    0.458 0.409 16  1.120  0.5162 
#L.B - L.R    0.216 0.409 16  0.528  0.8589 
#S.B - L.R   -0.242 0.378 16 -0.639  0.8010 

#precip = D:
#  contrast  estimate    SE df t.ratio p.value
#L.B - S.B    1.609 0.409 16  3.935  0.0032 
#L.B - L.R   -0.431 0.464 16 -0.929  0.6305 
#S.B - L.R   -2.039 0.450 16 -4.526  0.0009 


SG_soil_nit_trans_model_no_block= lm((ug_N_NO3_g_dry_soil)~precip*soil_root, data= SG_inorg_N_trans_trt)
qqPlot(resid(SG_soil_nit_trans_model_no_block))
hist(resid(SG_soil_nit_trans_model_no_block))
shapiro.test(resid(SG_soil_nit_trans_model_no_block))
#0.2901
anova(SG_soil_nit_trans_model_no_block, type=3)
#precip             5.28  1   13.8979  0.001328 ** 
#soil_root          6.34  2    8.3557  0.002303 ** 
#precip:soil_root   3.18  2    4.1948  0.030110 *  

emmeans(SG_soil_nit_trans_model_no_block, pairwise~soil_root|precip)

AIC(SG_soil_nit_trans_model,SG_soil_nit_trans_model_no_block)

SG_inorg_N_trans_trt %>% group_by(soil_root,precip) %>% summarise_at("ug_N_NO3_g_dry_soil", funs(n(),mean,sd,se=sd(.)/sqrt(n())))


trans_inorg_N_seed_trt_soil_root_precip_g=SG_inorg_N_trans_trt %>% group_by(soil_root,precip)
soil_nit_trans_trt=summarise_at(trans_inorg_N_seed_trt_soil_root_precip_g, 
                               "ug_N_NO3_g_dry_soil", funs(n(),mean,sd,se=sd(.)/sqrt(n())))

ggplot(soil_nit_trans_trt, aes(x=precip,y=mean, ymin = mean-se, ymax= mean+se,fill=soil_root))+
  geom_bar(position="dodge",stat="identity", color="black",aes(fill=soil_root))+
  geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+
  scale_fill_manual(values = c("gray", "white", "black"),labels=c("Bulk","Sterile","Rhizo"))+
  scale_x_discrete(labels=c("Ambient","Drought"))+ylab("Soil nitrate (ug per g soil)")+
  geom_text(aes(y=mean+se+0.3, label=n),position=position_dodge(width=0.9))+theme_bw()

SG_inorg_N_trans_trt$ug_N_NH4_g_dry_soil_negto0


#####Transplant nitrate reformated graph####

(trans_soil_nitrate_p=ggplot(soil_nit_trans_trt, aes(x=precip,y=mean, ymin = mean-se, ymax= mean+se,fill=factor(soil_root,levels = treatment_order)))+
   geom_bar(position="dodge",stat="identity", color="black",aes(fill=factor(soil_root,levels = treatment_order)))+
   geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+ylim(c(0,9))+
   scale_fill_manual(values = c( "white","lightgray", "darkgrey"),labels=c("Sterile","Bulk","Rhizo"), limits=treatment_order, name=NULL)+
   scale_x_discrete(labels=c("Ambient","Drought"), name=NULL)+ylab("Soil nitrate (ug per g soil)")+
   geom_text(aes(y=.5, label=n),size=10,position=position_dodge(width=0.9))+theme_bw()+
   theme(axis.title = element_text(size = 23), axis.text = element_text(size = 23),
         legend.position = c(0.85,.9), legend.text=element_text(size=20),
         legend.background = element_rect(size=0.5,linetype="solid",colour ="black"),
         panel.grid.major = element_blank(), panel.grid.minor = element_blank()))
(trans_soil_nitrate_p2=ggplot(soil_nit_trans_trt, aes(x=precip,y=mean, ymin = mean-se, ymax= mean+se,fill=factor(soil_root,levels = treatment_order)))+
    geom_bar(position="dodge",stat="identity", color="black",aes(fill=factor(soil_root,levels = treatment_order)))+
    geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+ylim(c(0,9))+
    scale_fill_manual(values = c( "white","lightgray", "darkgrey"),labels=c("Sterile","Bulk","Rhizo"), limits=treatment_order, name=NULL)+
    scale_x_discrete(labels=c("Ambient","Drought"), name=NULL)+ylab("Soil nitrate (ug per g soil)")+
    geom_text(aes(y=.5, label=n),size=10,position=position_dodge(width=0.9))+theme_bw()+
    theme(axis.title = element_text(size = 23), axis.text = element_text(size = 23),
          legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank()))
#####Combined germination nitrate graph#####
ggarrange(trans_soil_nitrate_p,seed_soil_nitrate_p,ncol = 2,  legend = "none",widths = c(1,.90))
#15x7.38
attributes(trans_soil_nitrate_p)
attributes(seed_soil_nitrate_p)
summary(trans_soil_nitrate_p)
summary(seed_soil_nitrate_p)




####Presence nitrate transplant####
SG_inorg_N_trans_pres=subset(SG_inorg_N_trans_trt, root_association =="B")
unique(SG_inorg_N_trans_pres$soil_root)

SG_soil_nit_trans_model_pres= lm((ug_N_NO3_g_dry_soil)~precip*soil_root+(1|block), data= SG_inorg_N_trans_pres)
qqPlot(resid(SG_soil_nit_trans_model_pres))
hist(resid(SG_soil_nit_trans_model_pres))
boxCox(SG_soil_nit_trans_model_pres)
shapiro.test(resid(SG_soil_nit_trans_model_pres))
#0.2007

anova(SG_soil_nit_trans_model_pres, type=3)
#precip             1.387  1   5.3310 0.043589 *  
#soil_root          3.910  1  15.0244 0.003078 ** 
#precip:soil_root   1.471  1   5.6548 0.038745 *

SG_inorg_N_trans_trt %>% group_by(precip) %>% summarise_at("ug_N_NO3_g_dry_soil", funs(n(),mean,sd,se=sd(.)/sqrt(n())))

emmeans(SG_soil_nit_trans_model_pres, ~precip)
emmeans(SG_soil_nit_trans_model_pres, pairwise~soil_root)
emmeans(SG_soil_nit_trans_model_pres, pairwise~soil_root|precip)
emmeans(SG_soil_nit_trans_model_pres, pairwise~soil_root*precip,adjust="none")


SG_soil_nit_trans_model_pres_no_block= lm((ug_N_NO3_g_dry_soil)~precip*soil_root, data= SG_inorg_N_trans_pres)
qqPlot(resid(SG_soil_nit_trans_model_pres_no_block))
hist(resid(SG_soil_nit_trans_model_pres_no_block))
shapiro.test(resid(SG_soil_nit_trans_model_pres_no_block))
#0.04276
anova(SG_soil_nit_trans_model_pres_no_block, type=3)
#precip             1.387  1    5.5927  0.033012 *  
#soil_root          4.219  1   17.0078  0.001032 ** 
#precip:soil_root   1.471  1    5.9324  0.028828 *  

emmeans(SG_soil_nit_trans_model_pres_no_block, pairwise~soil_root|precip)

AIC(SG_soil_nit_trans_model_pres,SG_soil_nit_trans_model_pres_no_block)


####Origin nitrate transplant####


SG_inorg_N_trans_orig=subset(SG_inorg_N_trans_trt,soil_status =="L")
unique(SG_inorg_N_trans_orig$soil_root)

SG_soil_nit_trans_model_orig= lm((ug_N_NO3_g_dry_soil)~precip*soil_root+(1|block), data= SG_inorg_N_trans_orig)
qqPlot(resid(SG_soil_nit_trans_model_orig))
hist(resid(SG_soil_nit_trans_model_orig))
boxCox(SG_soil_nit_trans_model_orig)
shapiro.test(resid(SG_soil_nit_trans_model_orig))
#0.7787

anova(SG_soil_nit_trans_model_orig, type=3)
#precip             8.274  1  14.8091   0.00489 ** 

SG_inorg_N_trans_trt %>% group_by(precip) %>% summarise_at("ug_N_NO3_g_dry_soil", funs(n(),mean,sd,se=sd(.)/sqrt(n())))

emmeans(SG_soil_nit_trans_model_orig, ~precip)
emmeans(SG_soil_nit_trans_model_orig, pairwise~soil_root)
emmeans(SG_soil_nit_trans_model_orig, pairwise~soil_root|precip)
emmeans(SG_soil_nit_trans_model_orig, pairwise~soil_root*precip,adjust="none")


SG_soil_nit_trans_model_orig_no_block= lm((ug_N_NO3_g_dry_soil)~precip*soil_root, data= SG_inorg_N_trans_orig)
qqPlot(resid(SG_soil_nit_trans_model_orig_no_block))
hist(resid(SG_soil_nit_trans_model_orig_no_block))
shapiro.test(resid(SG_soil_nit_trans_model_orig_no_block))
#0.3
anova(SG_soil_nit_trans_model_orig_no_block, type=3)
#precip             7.407  1  12.0132  0.004665 ** 

emmeans(SG_soil_nit_trans_model_orig_no_block, pairwise~soil_root*precip,adjust="none")

AIC(SG_soil_nit_trans_model_orig,SG_soil_nit_trans_model_orig_no_block)




#Amonnium
SG_soil_amm_trans_model= lm(log(ug_N_NH4_g_dry_soil_negto0+0.007)~precip*soil_root+(1|block), data= SG_inorg_N_trans_trt)
qqPlot(resid(SG_soil_amm_trans_model))
hist(resid(SG_soil_amm_trans_model))
#plot(SG_soil_amm_trans_model)
shapiro.test(resid(SG_soil_amm_trans_model))
#0.2646

anova(SG_soil_amm_trans_model, type=3)
#precip            8.035  1  5.4444 0.0330058 *  
#soil_root        37.549  2 12.7215 0.0004936 ***
#precip:soil_root  9.350  2  3.1677 0.0693432 . 

SG_inorg_N_trans_trt %>% group_by(precip) %>% summarise_at("ug_N_NH4_g_dry_soil_negto0", funs(n(),mean,sd,se=sd(.)/sqrt(n())))
emmeans(SG_soil_amm_trans_model, ~precip)
emmeans(SG_soil_amm_trans_model, pairwise~soil_root)
emmeans(SG_soil_amm_trans_model, pairwise~soil_root|precip)


SG_soil_amm_trans_model_no_block= lm(log(ug_N_NH4_g_dry_soil_negto0+0.007)~precip*soil_root, data= SG_inorg_N_trans_trt)
qqPlot(resid(SG_soil_amm_trans_model_no_block))
hist(resid(SG_soil_amm_trans_model_no_block))
#plot(SG_soil_amm_trans_model)
shapiro.test(resid(SG_soil_amm_trans_model_no_block))
#0.3981

anova(SG_soil_amm_trans_model_no_block, type=3)
#precip            7.010  1  5.3624   0.03131 *  
#soil_root        43.257  2 16.5458 5.755e-05 ***
#precip:soil_root  8.962  2  3.4281   0.05246 . 

emmeans(SG_soil_amm_trans_model_no_block, pairwise~soil_root|precip)

AIC(SG_soil_amm_trans_model,SG_soil_amm_trans_model_no_block)

SG_inorg_N_trans_trt %>% group_by(precip) %>% summarise_at("ug_N_NH4_g_dry_soil_negto0", funs(n(),mean,sd,se=sd(.)/sqrt(n())))
(1.46-0.215)/0.215


SG_inorg_N_trans_trt %>% group_by(soil_root) %>% summarise_at("ug_N_NH4_g_dry_soil_negto0", funs(n(),mean,sd,se=sd(.)/sqrt(n())))
(1.86-0.117)/0.117

SG_inorg_N_trans_trt %>% group_by(precip,soil_root) %>% summarise_at("ug_N_NH4_g_dry_soil_negto0", funs(n(),mean,sd,se=sd(.)/sqrt(n())))
(3.30-0.0910)/0.0910
#35.26374


trans_inorg_N_seed_trt_soil_root_precip_g=SG_inorg_N_trans_trt %>% group_by(soil_root,precip)
soil_amm_trans_trt=summarise_at(trans_inorg_N_seed_trt_soil_root_precip_g, 
                               "ug_N_NH4_g_dry_soil_negto0", funs(n(),mean,sd,se=sd(.)/sqrt(n())))

ggplot(soil_amm_trans_trt, aes(x=precip,y=mean, ymin = mean-se, ymax= mean+se,fill=soil_root))+
  geom_bar(position="dodge",stat="identity", color="black",aes(fill=soil_root))+
  geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+
  scale_fill_manual(values = c("gray", "white", "black"),labels=c("Bulk","Sterile","Rhizo"))+
  scale_x_discrete(labels=c("Ambient","Drought"))+ylab("Soil ammonium (ug per g soil)")+
  geom_text(aes(y=mean+se+0.3, label=n),position=position_dodge(width=0.9))+theme_bw()


#####Transplant ammonium reformated graph####

(trans_soil_amm_p=ggplot(soil_amm_trans_trt, aes(x=precip,y=mean, ymin = mean-se, ymax= mean+se,fill=factor(soil_root,levels = treatment_order)))+
   geom_bar(position="dodge",stat="identity", color="black",aes(fill=factor(soil_root,levels = treatment_order)))+
   geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+ylim(c(0,4.5))+
   scale_fill_manual(values = c( "white","lightgray", "darkgrey"),labels=c("Sterile","Bulk","Rhizo"), limits=treatment_order, name=NULL)+
   scale_x_discrete(labels=c("Ambient","Drought"), name=NULL)+ylab("Soil ammonium (ug per g soil)")+
   geom_text(aes(y=mean+se+0.3, label=n),size=10,position=position_dodge(width=0.9))+theme_bw()+
   theme(axis.title = element_text(size = 23), axis.text = element_text(size = 23),
         legend.position = c(0.85,.9), legend.text=element_text(size=20),
         legend.background = element_rect(size=0.5,linetype="solid",colour ="black")))

(trans_soil_amm_p=ggplot(soil_amm_trans_trt, aes(x=precip,y=mean, ymin = mean-se, ymax= mean+se,fill=factor(soil_root,levels = treatment_order)))+
    geom_bar(position="dodge",stat="identity", color="black",aes(fill=factor(soil_root,levels = treatment_order)))+
    geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+ylim(c(0,4.5))+
    scale_fill_manual(values = c( "white","lightgray", "darkgrey"),labels=c("Sterile","Bulk","Rhizo"), limits=treatment_order, name=NULL)+
    scale_x_discrete(labels=c("Ambient","Drought"), name=NULL)+ylab("Soil ammonium (ug per g soil)")+
    geom_text(aes(y=mean+se+0.3, label=n),size=10,position=position_dodge(width=0.9))+theme_bw()+
    theme(axis.title = element_text(size = 23), axis.text = element_text(size = 23),
          legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

#
#####Combined germination ammonium graph#####
ggarrange(trans_soil_amm_p,seed_soil_amm_p,ncol = 2,  legend = "none",widths = c(1,.90))
#15x7.38


####Presence Amonnium####
SG_inorg_N_trans_pres=subset(SG_inorg_N_trans_trt, root_association =="B")
unique(SG_inorg_N_trans_pres$soil_root)


SG_soil_amm_trans_model_pres= lm((ug_N_NH4_g_dry_soil_negto0)~precip*soil_root+(1|block), data= SG_inorg_N_trans_pres)
qqPlot(resid(SG_soil_amm_trans_model_pres))
hist(resid(SG_soil_amm_trans_model_pres))
#plot(SG_soil_amm_trans_model)
shapiro.test(resid(SG_soil_amm_trans_model_pres))
#0.1107
boxCox(SG_soil_amm_trans_model_pres)
anova(SG_soil_amm_trans_model_pres, type=3)
#precip            8.9092  1 15.8691 0.0025859 ** 
#soil_root        12.6845  1 22.5938 0.0007764 ***
#precip:soil_root  9.5648  1 17.0370 0.0020521 ** 


SG_inorg_N_trans_trt %>% group_by(precip) %>% summarise_at("ug_N_NH4_g_dry_soil_negto0", funs(n(),mean,sd,se=sd(.)/sqrt(n())))
emmeans(SG_soil_amm_trans_model_pres, ~precip)
emmeans(SG_soil_amm_trans_model_pres, pairwise~soil_root)
emmeans(SG_soil_amm_trans_model_pres, pairwise~soil_root|precip)


SG_soil_amm_trans_model_pres_no_block= lm((ug_N_NH4_g_dry_soil_negto0)~precip*soil_root, data= SG_inorg_N_trans_pres)
qqPlot(resid(SG_soil_amm_trans_model_pres_no_block))
hist(resid(SG_soil_amm_trans_model_pres_no_block))
#plot(SG_soil_amm_trans_model)
shapiro.test(resid(SG_soil_amm_trans_model_pres_no_block))
#0.001322

anova(SG_soil_amm_trans_model_pres_no_block, type=3)
#precip            8.9092  1  17.509 0.0009180 ***
#soil_root        13.5020  1  26.535 0.0001472 ***
# precip:soil_root  9.5648  1  18.798 0.0006847 ***

emmeans(SG_soil_amm_trans_model_pres_no_block, pairwise~soil_root|precip)
emmeans(SG_soil_amm_trans_model_pres_no_block, pairwise~soil_root*precip,adjust="none")


AIC(SG_soil_amm_trans_model_pres,SG_soil_amm_trans_model_pres_no_block)


####Origin ammonium transplant####


SG_inorg_N_trans_orig=subset(SG_inorg_N_trans_trt,soil_status =="L")
unique(SG_inorg_N_trans_orig$soil_root)


SG_soil_amm_trans_model_orig= lm((ug_N_NH4_g_dry_soil_negto0)~precip*soil_root+(1|block), data= SG_inorg_N_trans_orig)
qqPlot(resid(SG_soil_amm_trans_model_orig))
hist(resid(SG_soil_amm_trans_model_orig))
#plot(SG_soil_amm_trans_model)
shapiro.test(resid(SG_soil_amm_trans_model_orig))
#0.8792
boxCox(SG_soil_amm_trans_model_orig)
anova(SG_soil_amm_trans_model_orig, type=3)
#nada


SG_inorg_N_trans_trt %>% group_by(precip) %>% summarise_at("ug_N_NH4_g_dry_soil_negto0", funs(n(),mean,sd,se=sd(.)/sqrt(n())))
emmeans(SG_soil_amm_trans_model_orig, ~precip)
emmeans(SG_soil_amm_trans_model_orig, pairwise~soil_root)
emmeans(SG_soil_amm_trans_model_orig, pairwise~soil_root|precip)


SG_soil_amm_trans_model_orig_no_block= lm((ug_N_NH4_g_dry_soil_negto0)~precip*soil_root, data= SG_inorg_N_trans_orig)
qqPlot(resid(SG_soil_amm_trans_model_orig_no_block))
hist(resid(SG_soil_amm_trans_model_orig_no_block))
#plot(SG_soil_amm_trans_model)
shapiro.test(resid(SG_soil_amm_trans_model_orig_no_block))
#0.1014

anova(SG_soil_amm_trans_model_orig_no_block, type=3)
#nada

emmeans(SG_soil_amm_trans_model_orig_no_block, pairwise~soil_root|precip)
emmeans(SG_soil_amm_trans_model_orig_no_block, pairwise~soil_root*precip,adjust="none")


AIC(SG_soil_amm_trans_model_orig,SG_soil_amm_trans_model_orig_no_block)





attributes(seed_soil_amm_p)
attributes(trans_soil_amm_p)
#Using soil moisture instead of the precip treat 

dataSG_potweight_trans_end=dataSG_potweight_trt_trans_end[,c("Plant_Number","per_ch_potweight","per_ch_potweight_bio_crt")]
colnames(dataSG_potweight_trans_end)[2]="end_per_ch_potweight"
colnames(dataSG_potweight_trans_end)[3]="end_per_ch_potweight_bio_crt"
nrow(dataSG_potweight_trans_end)

SG_inorg_N_trans_trt_pw=merge(SG_inorg_N_trans_trt,dataSG_potweight_trans_end, by="Plant_Number", all.x = T)
summary(SG_inorg_N_trans_trt_pw)
#nitrate
SG_soil_nit_trans_model_pw= lm((ug_N_NO3_g_dry_soil)^-1~end_per_ch_potweight_bio_crt*soil_root+(1|block), data= SG_inorg_N_trans_trt_pw)
qqPlot(resid(SG_soil_nit_trans_model_pw))
hist(resid(SG_soil_nit_trans_model_pw))
boxCox(SG_soil_nit_trans_model_pw)
anova(SG_soil_nit_trans_model_pw, type=3)
#nada sig

ggplot(SG_inorg_N_trans_trt_pw, aes(x=end_per_ch_potweight_bio_crt,y=ug_N_NO3_g_dry_soil))+geom_point(aes(fill=soil_root), size=4, shape=21)+
  theme_bw()+scale_fill_manual(values = c("gray", "white", "black"),labels=c("Bulk","Sterile","Rhizo"))

#Amonnium
SG_soil_amm_trans_model_pw= lm(log(ug_N_NH4_g_dry_soil_negto0+.007)~end_per_ch_potweight_bio_crt*soil_root+(1|block), data= SG_inorg_N_trans_trt_pw)
qqPlot(resid(SG_soil_amm_trans_model_pw))
hist(resid(SG_soil_amm_trans_model_pw))
shapiro.test(resid(SG_soil_amm_trans_model_pw))
anova(SG_soil_amm_trans_model_pw, type=3)
#nada sig

emmeans(SG_soil_amm_trans_model_pw, pairwise~soil_root)

ggplot(SG_inorg_N_trans_trt_pw, aes(x=end_per_ch_potweight,y=ug_N_NH4_g_dry_soil_negto0))+geom_point(aes(fill=soil_root), size=4, shape=21)+
  theme_bw()+scale_fill_manual(values = c("gray", "white", "black"),labels=c("Bulk","Sterile","Rhizo"))

#Is change in pot weight correlated with gravimetric water

trans_gravimetric_mod=lm(log(percent_soil_moisture_dry_weight)~end_per_ch_potweight_bio_crt,data=SG_inorg_N_trans_trt_pw)
qqPlot(resid(trans_gravimetric_mod))
hist(resid(trans_gravimetric_mod))
#plot(trans_gravimetric_mod)
summary(trans_gravimetric_mod)
plot(SG_inorg_N_trans_trt_pw$end_per_ch_potweight,SG_inorg_N_trans_trt_pw$percent_soil_moisture_dry_weight)

ggplot(SG_inorg_N_trans_trt_pw, aes(x=end_per_ch_potweight_bio_crt,y=percent_soil_moisture_dry_weight))+geom_point(aes(fill=soil_root), size=4, shape=21)+
  theme_bw()+scale_fill_manual(values = c("gray", "white", "black"),labels=c("Bulk","Sterile","Rhizo"))+
  geom_text(aes(y=percent_soil_moisture_dry_weight+.7, x= end_per_ch_potweight,label=Plant_Number))


#Treat effects on the gravimetric 



SG_soil_VWC_trans_model= lm((percent_soil_moisture_dry_weight)~precip*soil_root+(1|block), data= SG_inorg_N_trans_trt_pw)
qqPlot(resid(SG_soil_VWC_trans_model))
hist(resid(SG_soil_VWC_trans_model))
#plot(SG_soil_VWC_trans_model)
anova(SG_soil_VWC_trans_model, type=3)
#precip            281.7  1  5.2726   0.03551 *  



emmeans(SG_soil_VWC_trans_model, ~precip)
emmeans(SG_soil_VWC_trans_model, pairwise~soil_root)
emmeans(SG_soil_VWC_trans_model, pairwise~soil_root|precip)


#####TRANSPLANT PORTION OF THE EXPERIEMENT ENDING####


#####Soil nitrogen####

summary(SG_inorg_N_trt)
SG_inorg_N_trt$soil_root=with(SG_inorg_N_trt, interaction(soil_status,root_association))

#SOil nitrate
SG_Soil_nit_combined= lm((ug_N_NO3_g_dry_soil)~precip*soil_root*life_stage+(1|block), data= SG_inorg_N_trt)
qqPlot(resid(SG_Soil_nit_combined))
hist(resid(SG_Soil_nit_combined))
shapiro.test(resid(SG_Soil_nit_combined))
#0.9184
anova(SG_Soil_nit_combined, type=3)
#precip                       15.96  1   45.9200 6.533e-07 ***
#soil_root                    13.48  2   19.3870 1.163e-05 ***
#life_stage                    1.73  1    4.9640  0.035945 *  
#precip:soil_root              6.09  2    8.7622  0.001483 ** 
#precip:life_stage             2.29  1    6.5836  0.017279 *  
#soil_root:life_stage          3.12  2    4.4886  0.022606 *  
#precip:soil_root:life_stage   4.35  2    6.2588  0.006757 ** 


#SOil nitrate
SG_Soil_nit_combined_no_block= lm((ug_N_NO3_g_dry_soil)~precip*soil_root*life_stage, data= SG_inorg_N_trt)
qqPlot(resid(SG_Soil_nit_combined_no_block))
hist(resid(SG_Soil_nit_combined_no_block))
shapiro.test(resid(SG_Soil_nit_combined_no_block))
#0.1548
anova(SG_Soil_nit_combined_no_block, type=3)
#precip                       16.72  1   46.7861 2.391e-07 ***
#soil_root                    14.65  2   20.4851 3.865e-06 ***
#life_stage                    2.14  1    5.9838  0.021235 *  
#precip:soil_root              9.46  2   13.2326 9.873e-05 ***
#precip:life_stage             2.67  1    7.4703  0.010931 *  
#soil_root:life_stage          3.60  2    5.0322  0.013883 *  
#precip:soil_root:life_stage   5.34  2    7.4752  0.002609 ** 

emmeans(SG_Soil_nit_combined_no_block, pairwise~soil_root|precip)


AIC(SG_Soil_nit_combined,SG_Soil_nit_combined_no_block)

#SOil ammonium
SG_Soil_amm_combined= lm(log(ug_N_NH4_g_dry_soil_negto0+.001)~precip*soil_root*life_stage+(1|block), data= SG_inorg_N_trt)
qqPlot(resid(SG_Soil_amm_combined))
hist(resid(SG_Soil_amm_combined))
shapiro.test(resid(SG_Soil_amm_combined))
#0.7772
anova(SG_Soil_amm_combined, type=3)
#precip                       37.947  1  9.5370  0.005191 ** 
#soil_root                    36.168  2  4.5450  0.021708 *  

SG_Soil_amm_combined_no_block= lm(log(ug_N_NH4_g_dry_soil_negto0+.001)~precip*soil_root*life_stage, data= SG_inorg_N_trt)
qqPlot(resid(SG_Soil_amm_combined_no_block))
hist(resid(SG_Soil_amm_combined_no_block))
shapiro.test(resid(SG_Soil_amm_combined_no_block))
#0.4174
anova(SG_Soil_amm_combined_no_block, type=3)
#precip                       44.035  1 11.5159  0.002146 ** 
#soil_root                    36.989  2  4.8366  0.016021 *   
#precip:life_stage            13.906  1  3.6367  0.067211 .  

emmeans(SG_Soil_amm_combined_no_block, pairwise~soil_root|precip)

AIC(SG_Soil_amm_combined,SG_Soil_amm_combined_no_block)

#percentage change in soil moisture
nrow(dataSG_potweight_trt)
#900
summary(dataSG_potweight_trt)
dataSG_potweight_trt$soil_root=with(dataSG_potweight_trt, interaction(soil_status,root_association))
dataSG_potweight_trt_sub=subset(dataSG_potweight_trt, date!="6_13")
summary(dataSG_potweight_trt_sub)
#Extremely large outliers in bulk treatment

#I am going to remove pot 20 and 71 since their intial weight is suspicious and the very lowest values

dataSG_potweight_trt_sub_out=subset(dataSG_potweight_trt_sub, per_ch_potweight<.50)
dataSG_potweight_trt_sub_out=subset(dataSG_potweight_trt_sub_out, per_ch_potweight>-.3)
nrow(dataSG_potweight_trt_sub_out)
#711

#I am going to use the end chnage in pot weight
dataSG_potweight_trt_sub_out_end=subset(dataSG_potweight_trt_sub_out,date=="7_7")
nrow(dataSG_potweight_trt_sub_out_end)
#178
dataSG_potweight_end=dataSG_potweight_trt_sub_out_end[,c("Plant_Number","per_ch_potweight", "per_ch_potweight_bio_crt")]
colnames(dataSG_potweight_end)[2]="end_per_ch_potweight"
colnames(dataSG_potweight_end)[3]="end_per_ch_potweight_bio_crt"


#####VWC Analyses####

SG_Soil_potweight_combined= lm((per_ch_potweight_bio_crt)~precip*soil_root*life_stage+(1|block), data= dataSG_potweight_trt_sub_out_end)
qqPlot(resid(SG_Soil_potweight_combined))
hist(resid(SG_Soil_potweight_combined))
shapiro.test(resid(SG_Soil_potweight_combined))
#0.8658

anova(SG_Soil_potweight_combined, type=3)
#precip                      0.05518   1   74.0995 6.202e-15 ***
#soil_root                   0.10599   2   71.1692 < 2.2e-16 ***
#life_stage                  0.00647   1    8.6841  0.003684 ** 
#soil_root:life_stage        0.00436   2    2.9291  0.056282 .  

emmeans(SG_Soil_potweight_combined, pairwise~soil_root|life_stage)
emmeans(SG_Soil_potweight_combined, pairwise~life_stage|soil_root)

SG_Soil_potweight_combined_no_block= lm((per_ch_potweight_bio_crt)~precip*soil_root*life_stage, data= dataSG_potweight_trt_sub_out_end)
qqPlot(resid(SG_Soil_potweight_combined_no_block))
hist(resid(SG_Soil_potweight_combined_no_block))
shapiro.test(resid(SG_Soil_potweight_combined_no_block))
#0.9003

anova(SG_Soil_potweight_combined_no_block, type=3)
#precip                      0.05504   1   73.9038  5.85e-15 ***
#soil_root                   0.10631   2   71.3708 < 2.2e-16 ***
#life_stage                  0.00656   1    8.8063  0.003446 ** 
#soil_root:life_stage        0.00436   2    2.9276  0.056293 .  

AIC(SG_Soil_potweight_combined,SG_Soil_potweight_combined_no_block)

emmeans(SG_Soil_potweight_combined_no_block, pairwise~soil_root|precip)


dataSG_potweight_trt_sub_out_end %>% group_by(precip)  %>% summarise_at("per_ch_potweight_bio_crt", list(~mean(.),~var(.)))

#  precip    mean     var
#1 A      -0.0835 0.00153
#2 D      -0.119  0.00124

dataSG_potweight_trt_sub_out_end %>% group_by(soil_root)  %>% summarise_at("per_ch_potweight_bio_crt", list(~mean(.),~var(.)))
#soil_root    mean      var
#1 L.B       -0.122  0.000947
#2 S.B       -0.0672 0.00157 
#3 L.R       -0.115  0.000789

#S.B to L.R

(-0.115--0.0672)/-0.0672


dataSG_potweight_trt_sub_out_end %>% group_by(soil_root,precip)  %>% summarise_at("per_ch_potweight_bio_crt", list(~mean(.),~var(.)))
#L.B
(-0.103--0.141)/-0.141
#-0.2695035
#S.B 
(-0.0480--0.0864)/-0.0864
#-0.4444444
#L.R
(-0.100--0.130)/-0.130
#-0.2307692

emmeans(SG_Soil_potweight_combined_no_block, pairwise~soil_root|precip)

combin_pot_weight_seed_trt_soil_root_precip_g=dataSG_potweight_trt_sub_out_end %>% group_by(soil_root,precip)
soil_pot_weight_out_trt=summarise_at(combin_pot_weight_seed_trt_soil_root_precip_g, 
                                "per_ch_potweight_bio_crt", funs(n(),mean,sd,se=sd(.)/sqrt(n())))

treatment_order=c("S.B","L.B","L.R")
(combined_soil_moist_p=ggplot(soil_pot_weight_out_trt, aes(x=precip,y=mean, ymin = mean-se, ymax= mean+se,fill=factor(soil_root,levels = treatment_order)))+
    geom_bar(position="dodge",stat="identity", color="black",aes(fill=factor(soil_root,levels = treatment_order)))+
    geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+ylim(c(-0.16,0))+
    scale_fill_manual(values = c( "white","light gray", "dark grey"),labels=c("Sterile","Bulk","Rhizo"), limits=treatment_order, name=NULL)+
    scale_x_discrete(labels=c("Ambient","Drought"), name=NULL)+ylab("prop Change \nin soil moisture from start")+
    geom_text(aes(y=-.01, label=n),size=10,position=position_dodge(width=0.9))+theme_bw()+
    theme(axis.title = element_text(size = 23), axis.text = element_text(size = 23),
          legend.position = "none"))
#9x8
#percentage change in soil moisture
nrow(SG_Soil_potweight_combined)
#900
summary(dataSG_potweight_trt)
dataSG_potweight_trt$soil_root=with(dataSG_potweight_trt, interaction(soil_status,root_association))
dataSG_potweight_trt_sub=subset(dataSG_potweight_trt, date!="6_13")
summary(dataSG_potweight_trt_sub)


SG_inorg_N_trt_pw=merge(SG_inorg_N_trt,dataSG_potweight_end, by="Plant_Number", all.x = T)


#Is change in pot weight correlated with gravimetric water

soil_gravimetric_mod=lm(log(percent_soil_moisture_dry_weight)~end_per_ch_potweight,data=SG_inorg_N_trt_pw)
qqPlot(resid(soil_gravimetric_mod))
hist(resid(soil_gravimetric_mod))
#plot(soil_gravimetric_mod)
summary(soil_gravimetric_mod)
plot(SG_inorg_N_trt_pw$end_per_ch_potweight,log(SG_inorg_N_trt_pw$percent_soil_moisture_dry_weight))

ggplot(SG_inorg_N_trt_pw, aes(x=end_per_ch_potweight,y=percent_soil_moisture_dry_weight))+geom_point(aes(fill=soil_root), size=4, shape=21)+
  theme_bw()+scale_fill_manual(values = c("gray", "white", "black"),labels=c("Bulk","Sterile","Rhizo"))+
  geom_text(aes(y=percent_soil_moisture_dry_weight+.7, x= end_per_ch_potweight,label=Plant_Number))
ggplot(SG_inorg_N_trt_pw, aes(x=end_per_ch_potweight,y=log(percent_soil_moisture_dry_weight)))+geom_point(aes(fill=soil_root), size=4, shape=21)+
  theme_bw()+scale_fill_manual(values = c("gray", "white", "black"),labels=c("Bulk","Sterile","Rhizo"))+
  geom_text(aes(y=log(percent_soil_moisture_dry_weight)+.1, x= end_per_ch_potweight,label=Plant_Number))

#now using the corrected value to see if it improves the correlation
soil_gravimetric_mod_bio_crt=lm(log(percent_soil_moisture_dry_weight)~end_per_ch_potweight_bio_crt,data=SG_inorg_N_trt_pw)
qqPlot(resid(soil_gravimetric_mod_bio_crt))
hist(resid(soil_gravimetric_mod_bio_crt))
#plot(soil_gravimetric_mod)
summary(soil_gravimetric_mod_bio_crt)
plot(SG_inorg_N_trt_pw$end_per_ch_potweight_bio_crt,log(SG_inorg_N_trt_pw$percent_soil_moisture_dry_weight))

ggplot(SG_inorg_N_trt_pw, aes(x=end_per_ch_potweight_bio_crt,y=percent_soil_moisture_dry_weight))+geom_point(aes(fill=soil_root), size=4, shape=21)+
  theme_bw()+scale_fill_manual(values = c("gray", "white", "black"),labels=c("Bulk","Sterile","Rhizo"))+
  geom_text(aes(y=percent_soil_moisture_dry_weight+.7, x= end_per_ch_potweight,label=Plant_Number))
ggplot(SG_inorg_N_trt_pw, aes(x=end_per_ch_potweight_bio_crt,y=log(percent_soil_moisture_dry_weight)))+geom_point(aes(fill=soil_root), size=4, shape=21)+
  theme_bw()+scale_fill_manual(values = c("gray", "white", "black"),labels=c("Bulk","Sterile","Rhizo"))+
  geom_text(aes(y=log(percent_soil_moisture_dry_weight)+.1, x= end_per_ch_potweight,label=Plant_Number))


#####Gravimetric####

#Treat effects on the gravimetric 



SG_soil_VWC_model= lm((percent_soil_moisture_dry_weight)~precip*soil_root*life_stage+(1|block), data= SG_inorg_N_trt)
qqPlot(resid(SG_soil_VWC_model))
hist(resid(SG_soil_VWC_model))
#plot(SG_soil_VWC_trans_model)
shapiro.test(resid(SG_soil_VWC_model))
#0.2795

anova(SG_soil_VWC_model, type=3)
#precip                       259.11  1  6.4951  0.01796 *  



emmeans(SG_soil_VWC_model, ~precip)
emmeans(SG_soil_VWC_model, pairwise~soil_root)
emmeans(SG_soil_VWC_model, pairwise~soil_root|precip)

SG_soil_VWC_model_no_block= lm((percent_soil_moisture_dry_weight)~precip*soil_root*life_stage, data= SG_inorg_N_trt)
qqPlot(resid(SG_soil_VWC_model_no_block))
hist(resid(SG_soil_VWC_model_no_block))
#plot(SG_soil_VWC_trans_model)
shapiro.test(resid(SG_soil_VWC_model_no_block))
#0.3377

anova(SG_soil_VWC_model_no_block, type=3)
#precip                       307.12  1  7.9808  0.008782 ** 
#life_stage                   135.01  1  3.5084  0.071923 .  
emmeans(SG_soil_VWC_model_no_block, pairwise~soil_root|precip)


AIC(SG_soil_VWC_model,SG_soil_VWC_model_no_block)

combin_grav_moist_seed_trt_soil_root_precip_g=SG_inorg_N_trt %>% group_by(soil_root,precip)
soil_grav_moist_out_trt=summarise_at(combin_grav_moist_seed_trt_soil_root_precip_g, 
                                     "percent_soil_moisture_dry_weight", funs(n(),mean,sd,se=sd(.)/sqrt(n())))

treatment_order=c("S.B","L.B","L.R")
(combined_grav_soil_moist_p=ggplot(soil_grav_moist_out_trt, aes(x=precip,fill=factor(soil_root,levels = treatment_order)))+
    geom_bar(aes(y=mean),position="dodge",stat="identity", color="black")+
    geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1,aes(y=mean, ymin = mean-se, ymax= mean+se))+
    scale_fill_manual(values = c( "white","lightgray", "darkgrey"),labels=c("Sterile","Bulk","Rhizo"), limits=treatment_order, name=NULL)+
    scale_x_discrete(labels=c("Ambient","Drought"), name=NULL)+ylab("Gravimetric soil moisture (%)")+
    geom_text(aes(y=2, label=n),size=10,position=position_dodge(width=0.9))+theme_bw()+ylim(c(0,22))+
    theme(axis.title = element_text(size = 23), axis.text = element_text(size = 23),
          legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank()))
#9x8
soil_grav_moist_out_trt

#L.B A versus D
(10.1-15.9)/15.9
#-0.3647799
#S.B A versus D
(7.78-11.3)/11.3
#-0.3115044
#L.R A versus D
(6.77-15.2)/15.2
#-0.5546053


####Presence Gravimetric####
SG_inorg_N_trans_pres=subset(SG_inorg_N_trans_trt, root_association =="B")
unique(SG_inorg_N_trans_pres$soil_root)


####Origin Gravimetric####


SG_inorg_N_trans_orig=subset(SG_inorg_N_trans_trt,soil_status =="L")
unique(SG_inorg_N_trans_orig$soil_root)



#there is no significant effect of soil_root on soil moisture when using gravimetric
#Is this just because of the samples choosen to measure gravimetric moisture?

#let's look at the pot wieght of only these samples


SG_soil_VWC_pot_w_model= lm((end_per_ch_potweight_bio_crt)~precip*soil_root*life_stage+(1|block), data= SG_inorg_N_trt_pw)
qqPlot(resid(SG_soil_VWC_pot_w_model))
hist(resid(SG_soil_VWC_pot_w_model))
#plot(SG_soil_VWC_trans_model)
shapiro.test(resid(SG_soil_VWC_pot_w_model))
#0.7648

anova(SG_soil_VWC_pot_w_model, type=3)
#precip                      0.014801  1  14.3368 0.0009547 ***
#soil_root                   0.008510  2   4.1216 0.0295243 *  
#life_stage                  0.007848  1   7.6019 0.0112192 *



emmeans(SG_soil_VWC_pot_w_model, ~precip)
emmeans(SG_soil_VWC_pot_w_model, pairwise~soil_root)
emmeans(SG_soil_VWC_pot_w_model, pairwise~soil_root|precip)

SG_soil_VWC_pot_w_model_no_block= lm((end_per_ch_potweight_bio_crt)~precip*soil_root*life_stage, data= SG_inorg_N_trt_pw)
qqPlot(resid(SG_soil_VWC_pot_w_model_no_block))
hist(resid(SG_soil_VWC_pot_w_model_no_block))
#plot(SG_soil_VWC_trans_model)
shapiro.test(resid(SG_soil_VWC_pot_w_model_no_block))
#0.4286

anova(SG_soil_VWC_pot_w_model_no_block, type=3)
#precip                      0.01443  1  15.7850 0.0004753 ***
#soil_root                   0.00833  2   4.5563 0.0197257 *  
#life_stage                  0.00721  1   7.8899 0.0091285 ** 

emmeans(SG_soil_VWC_pot_w_model_no_block, pairwise~soil_root|precip)

AIC(SG_soil_VWC_pot_w_model,SG_soil_VWC_pot_w_model_no_block)




combin_grav_moist_PW_seed_trt_soil_root_precip_g=SG_inorg_N_trt_pw %>% group_by(soil_root,precip)
PW_soil_grav_moist_out_trt=summarise_at(combin_grav_moist_PW_seed_trt_soil_root_precip_g, 
                                     "end_per_ch_potweight_bio_crt", funs(n(),mean,sd,se=sd(.)/sqrt(n())))

treatment_order=c("S.B","L.B","L.R")
(combined_PW_grav_soil_moist_p=ggplot(PW_soil_grav_moist_out_trt, aes(x=precip,y=mean, ymin = mean-se, ymax= mean+se,fill=factor(soil_root,levels = treatment_order)))+
    geom_bar(position="dodge",stat="identity", color="black",aes(fill=factor(soil_root,levels = treatment_order)))+
    geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+
    scale_fill_manual(values = c( "white","light gray", "dark grey"),labels=c("Sterile","Bulk","Rhizo"), limits=treatment_order, name=NULL)+
    scale_x_discrete(labels=c("Ambient","Drought"), name=NULL)+ylab("prop Change \nin soil moisture from start")+
    geom_text(aes(y=-.012, label=n),size=10,position=position_dodge(width=0.9))+theme_bw()+
    theme(axis.title = element_text(size = 23), axis.text = element_text(size = 23),
          legend.position = "none"))


#Is plant biomass correlated with NO3 or NH4?
summary(data_SG_biomass)
SG_inorg_N_trt_bio_pw=merge(data_SG_biomass,SG_inorg_N_trt_pw, by="Plant_Number", all.y = T)
summary(SG_inorg_N_trt_bio_pw)

#nitrate

corrr::correlate(SG_inorg_N_trt_bio_pw[,c("root_weight_g","shoot_weight_g","total_biomass","end_per_ch_potweight_bio_crt","ug_N_NH4_g_dry_soil_negto0","ug_N_NO3_g_dry_soil")])


#####Start Soil nitrogen####
summary(SG_start_inorg_N)


#Soil nitrate
SG_start_inorg_N$soil_root=with(SG_start_inorg_N, interaction(soil_status,root_association))

#Start SOil nitrate
SG_Soil_nit_start= lm((ug_N_NO3_g_dry_soil)~soil_root+(1|block), data= SG_start_inorg_N)
qqPlot(resid(SG_Soil_nit_start))
hist(resid(SG_Soil_nit_start))
shapiro.test(resid(SG_Soil_nit_start))
#p-value = 0.1575
anova(SG_Soil_nit_start, type=3)

SG_Soil_nit_start_no_block= lm((ug_N_NO3_g_dry_soil)~soil_root, data= SG_start_inorg_N)
qqPlot(resid(SG_Soil_nit_start_no_block))
hist(resid(SG_Soil_nit_start_no_block))
shapiro.test(resid(SG_Soil_nit_start_no_block))
#p-value = 0.2278
anova(SG_Soil_nit_start_no_block, type=3)

AIC(SG_Soil_nit_start,SG_Soil_nit_start_no_block)

SG_start_inorg_N %>% group_by(soil_root) %>% summarise_at("ug_N_NO3_g_dry_soil", list(~n(),~mean(.),sd,se=~sd(.)/sqrt(n())))


SG_start_inorg_N_soil_root_g=SG_start_inorg_N %>% group_by(soil_root)
soil_nit_START_trt=summarise_at(SG_start_inorg_N_soil_root_g, 
                                "ug_N_NO3_g_dry_soil", list(~n(),~mean(.),sd,se=~sd(.)/sqrt(n())))
treatment_order=c("S.B","L.B","L.R")
(soil_nit_START_trt_p=ggplot(soil_nit_START_trt, aes(x=soil_root,y=mean, ymin = mean-se, ymax= mean+se,fill=factor(soil_root,levels = treatment_order)))+
    geom_bar(position="dodge",stat="identity", color="black",aes(fill=factor(soil_root,levels = treatment_order)))+
    geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+ylim(c(0,10.5))+
    scale_fill_manual(values = c( "white","lightgray", "darkgrey"),labels=NULL, limits=treatment_order, name=NULL)+
    scale_x_discrete(labels=c("Sterile","Bulk","Rhizo"), limits=treatment_order, name=NULL)+ylab("Starting soil nitrate \n(ug per g soil)")+
    geom_text(aes(y=.5, label=n),size=10,position=position_dodge(width=0.9))+theme_bw()+
    theme(axis.title = element_text(size = 23), axis.text = element_text(size = 23),
          legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank()))




####Presence Start SOil nitrate####
SG_start_inorg_N_pres=subset(SG_start_inorg_N, root_association =="B")
unique(SG_start_inorg_N_pres$soil_root)

SG_Soil_nit_start_pres= lm((ug_N_NO3_g_dry_soil)~soil_root+(1|block), data= SG_start_inorg_N_pres)
qqPlot(resid(SG_Soil_nit_start_pres))
hist(resid(SG_Soil_nit_start_pres))
shapiro.test(resid(SG_Soil_nit_start_pres))
#p-value = 0.5044
anova(SG_Soil_nit_start_pres, type=3)

SG_Soil_nit_start_pres_no_block= lm((ug_N_NO3_g_dry_soil)~soil_root, data= SG_start_inorg_N_pres)
qqPlot(resid(SG_Soil_nit_start_pres_no_block))
hist(resid(SG_Soil_nit_start_pres_no_block))
shapiro.test(resid(SG_Soil_nit_start_pres_no_block))
#p-value = 0.3099
anova(SG_Soil_nit_start_pres_no_block, type=3)

AIC(SG_Soil_nit_start_pres,SG_Soil_nit_start_pres_no_block)


####Origin Start SOil nitrate####
SG_start_inorg_N_orig=subset(SG_start_inorg_N, soil_status =="L")
unique(SG_start_inorg_N_orig$soil_root)

SG_Soil_nit_start_orig= lm((ug_N_NO3_g_dry_soil)~soil_root+(1|block), data= SG_start_inorg_N_orig)
qqPlot(resid(SG_Soil_nit_start_orig))
hist(resid(SG_Soil_nit_start_orig))
shapiro.test(resid(SG_Soil_nit_start_orig))
#p-value = 0.927
anova(SG_Soil_nit_start_orig, type=3)
#nada

SG_Soil_nit_start_orig_no_block= lm((ug_N_NO3_g_dry_soil)~soil_root, data= SG_start_inorg_N_orig)
qqPlot(resid(SG_Soil_nit_start_orig_no_block))
hist(resid(SG_Soil_nit_start_orig_no_block))
shapiro.test(resid(SG_Soil_nit_start_orig_no_block))
#p-value = 0.3099
anova(SG_Soil_nit_start_orig_no_block, type=3)

AIC(SG_Soil_nit_start_orig,SG_Soil_nit_start_orig_no_block)


#Soil ammonium
SG_start_inorg_N$soil_root=with(SG_start_inorg_N, interaction(soil_status,root_association))

#Start Soil ammonium
SG_Soil_amm_start= lm((ug_N_NH4_g_dry_soil)~soil_root+(1|block), data= SG_start_inorg_N)
qqPlot(resid(SG_Soil_amm_start))
hist(resid(SG_Soil_amm_start))
shapiro.test(resid(SG_Soil_amm_start))
#p-value = 0.2101
anova(SG_Soil_amm_start, type=3)
#soil_root        1911.17  2  65.7348 1.083e-05 ***

emmeans(SG_Soil_amm_start, pairwise~soil_root)

SG_Soil_amm_start_no_block= lm((ug_N_NH4_g_dry_soil)~soil_root, data= SG_start_inorg_N)
qqPlot(resid(SG_Soil_amm_start_no_block))
hist(resid(SG_Soil_amm_start_no_block))
shapiro.test(resid(SG_Soil_amm_start_no_block))
#p-value = 0.001853
anova(SG_Soil_amm_start_no_block, type=3)
#soil_root   1911.17  2  48.924 1.700e-06 ***

AIC(SG_Soil_amm_start,SG_Soil_amm_start_no_block)


SG_start_inorg_N %>% group_by(soil_root) %>% summarise_at("ug_N_NH4_g_dry_soil", list(~n(),~mean(.),~sd(.),se=~sd(.)/sqrt(n())))
#SB versus LB
(28.1-4.20)/4.20


SG_start_inorg_N_soil_root_g=SG_start_inorg_N %>% group_by(soil_root)
soil_amm_START_trt=summarise_at(SG_start_inorg_N_soil_root_g, 
                                "ug_N_NH4_g_dry_soil", list(~n(),~mean(.),sd,se=~sd(.)/sqrt(n())))
treatment_order=c("S.B","L.B","L.R")
(soil_amm_START_trt_p=ggplot(soil_amm_START_trt, aes(x=soil_root,y=mean, ymin = mean-se, ymax= mean+se,fill=factor(soil_root,levels = treatment_order)))+
    geom_bar(position="dodge",stat="identity", color="black",aes(fill=factor(soil_root,levels = treatment_order)))+
    geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+ylim(c(0,35))+
    scale_fill_manual(values = c( "white","lightgray", "darkgrey"),labels=NULL, limits=treatment_order, name=NULL)+
    scale_x_discrete(labels=c("Sterile","Bulk","Rhizo"), limits=treatment_order,name=NULL)+ylab("Starting soil ammonium \n(ug per g soil)")+
    geom_text(aes(y=1.5, label=n),size=10,position=position_dodge(width=0.9))+theme_bw()+
    theme(axis.title = element_text(size = 23), axis.text = element_text(size = 23),
          legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank()))
#ylim(c(0,9))+
#800x750
#SB versus LB
(28.1-4.20)/4.20
#SB versus LR
(28.1-4.04)/4.04
#5.955446

#####Combined start nitrate and ammonium graph#####
ggarrange(soil_nit_START_trt_p,soil_amm_START_trt_p,ncol = 2,  legend = "none")
#15x7.38

####Presence Soil ammonium####
SG_start_inorg_N_pres=subset(SG_start_inorg_N, root_association =="B")
unique(SG_start_inorg_N_pres$soil_root)

#Start Soil ammonium
SG_Soil_amm_start_pres= lm((ug_N_NH4_g_dry_soil)~soil_root+(1|block), data= SG_start_inorg_N_pres)
qqPlot(resid(SG_Soil_amm_start_pres))
hist(resid(SG_Soil_amm_start_pres))
shapiro.test(resid(SG_Soil_amm_start_pres))
#p-value = 0.9275
anova(SG_Soil_amm_start_pres, type=3)
#soil_root        1423.78  1  94.1106 0.0006320 ***

emmeans(SG_Soil_amm_start_pres, pairwise~soil_root)

SG_Soil_amm_start_no_block_pres= lm((ug_N_NH4_g_dry_soil)~soil_root, data= SG_start_inorg_N_pres)
qqPlot(resid(SG_Soil_amm_start_no_block_pres))
hist(resid(SG_Soil_amm_start_no_block_pres))
shapiro.test(resid(SG_Soil_amm_start_no_block_pres))
#p-value = 0.02782
boxCox(SG_Soil_amm_start_no_block_pres)
anova(SG_Soil_amm_start_no_block_pres, type=3)
#soil_root   1423.78  1  49.409 0.0001094 ***

AIC(SG_Soil_amm_start_pres,SG_Soil_amm_start_no_block_pres)

####Origin Soil ammonium####
SG_start_inorg_N_orig=subset(SG_start_inorg_N, soil_status =="L")
unique(SG_start_inorg_N_orig$soil_root)

#Start Soil ammonium
SG_Soil_amm_start_orig= lm((ug_N_NH4_g_dry_soil)^-1~soil_root+(1|block), data= SG_start_inorg_N_orig)
qqPlot(resid(SG_Soil_amm_start_orig))
hist(resid(SG_Soil_amm_start_orig))
shapiro.test(resid(SG_Soil_amm_start_orig))
#p-value = 0.9759
anova(SG_Soil_amm_start_orig, type=3)
#nada sig

emmeans(SG_Soil_amm_start_orig, pairwise~soil_root)

SG_Soil_amm_start_orig_no_block= lm((ug_N_NH4_g_dry_soil)^-1~soil_root, data= SG_start_inorg_N_orig)
qqPlot(resid(SG_Soil_amm_start_orig_no_block))
hist(resid(SG_Soil_amm_start_orig_no_block))
shapiro.test(resid(SG_Soil_amm_start_orig_no_block))
#p-value = 0.8879
boxCox(SG_Soil_amm_start_orig_no_block)
anova(SG_Soil_amm_start_orig_no_block, type=3)
#nada

AIC(SG_Soil_amm_start_orig,SG_Soil_amm_start_orig_no_block)

#####Change in soil nitrogen throughout the experiment####
summary(SG_start_inorg_N)
SG_start_inorg_N$soil_root_blck=with(SG_start_inorg_N,interaction(soil_root,block))
summary(SG_inorg_N_trt)
SG_inorg_N_trt$soil_root=with(SG_inorg_N_trt, interaction(soil_status,root_association))
SG_inorg_N_trt$soil_root_blck=with(SG_inorg_N_trt,interaction(soil_root,block))
nrow(SG_inorg_N_trt)
SG_inorg_N_trt_w_start=merge(SG_inorg_N_trt,SG_start_inorg_N[,c("ug_N_NO3_g_dry_soil",
                                                                "ug_N_NH4_g_dry_soil",
                                                                "soil_root_blck")], by="soil_root_blck", all.x=T)
nrow(SG_inorg_N_trt_w_start)
summary(SG_inorg_N_trt_w_start)
SG_inorg_N_trt_w_start$per_change_ug_N_NO3_g_dry_soil=(SG_inorg_N_trt_w_start$ug_N_NO3_g_dry_soil.x-
                                                         SG_inorg_N_trt_w_start$ug_N_NO3_g_dry_soil.y)/
  SG_inorg_N_trt_w_start$ug_N_NO3_g_dry_soil.y


SG_inorg_N_trt_w_start$per_change_ug_N_NH4_g_dry_soil=(SG_inorg_N_trt_w_start$ug_N_NH4_g_dry_soil_negto0-
                                                         SG_inorg_N_trt_w_start$ug_N_NH4_g_dry_soil.y)/
  SG_inorg_N_trt_w_start$ug_N_NH4_g_dry_soil.y


#SOil nitrate
SG_Soil_nit_per_chng_combined= lm((per_change_ug_N_NO3_g_dry_soil)~precip*soil_root*life_stage+(1|block), 
                                  data= SG_inorg_N_trt_w_start)
qqPlot(resid(SG_Soil_nit_per_chng_combined))
hist(resid(SG_Soil_nit_per_chng_combined))
shapiro.test(resid(SG_Soil_nit_per_chng_combined))
#p-value = 0.4422

anova(SG_Soil_nit_per_chng_combined, type=3)
#precip                      0.2432  1  21.1294 0.0001272 ***
#soil_root                   0.1779  2   7.7298 0.0027059 ** 
#(1|block)            0.2864  4   6.2214 0.0015176 ** 
#precip:soil_root            0.1028  2   4.4647 0.0229973 *  

emmeans(SG_Soil_nit_per_chng_combined, pairwise~soil_root|precip)

SG_per_chng_inorg_N_soil_root_g=SG_inorg_N_trt_w_start %>% group_by(soil_root,precip)
soil_per_chng_nit_START_trt=summarise_at(SG_per_chng_inorg_N_soil_root_g, 
                                "per_change_ug_N_NO3_g_dry_soil", list(~n(),~mean(.),sd,se=~sd(.)/sqrt(n())))
treatment_order=c("S.B","L.B","L.R")
(combine_per_chng_soil_nitrate_p=ggplot(soil_per_chng_nit_START_trt, aes(x=precip,y=mean, ymin = mean-se, ymax= mean+se,fill=factor(soil_root,levels = treatment_order)))+
    geom_bar(position="dodge",stat="identity", color="black",aes(fill=factor(soil_root,levels = treatment_order)))+
    geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+ylim(c(-0.65,0))+
    scale_fill_manual(values = c( "white","light gray", "dark grey"),labels=c("Sterile","Bulk","Rhizo"), limits=treatment_order, name=NULL)+
    scale_x_discrete(labels=c("Ambient","Drought"), name=NULL)+ylab("per Change in soil nitrate \n(Start to end)")+
    geom_text(aes(y=-.1, label=n),size=10,position=position_dodge(width=0.9))+theme_bw()+
    theme(axis.title = element_text(size = 23), axis.text = element_text(size = 23),
          legend.position = c(0.85,.9), legend.text=element_text(size=20),
          legend.background = element_rect(size=0.5,linetype="solid",colour ="black")))


#SOil ammonium
SG_Soil_amm_per_chng_combined= lm((per_change_ug_N_NH4_g_dry_soil)~precip*soil_root*life_stage+(1|block), 
                                  data= SG_inorg_N_trt_w_start)
qqPlot(resid(SG_Soil_amm_per_chng_combined))
hist(resid(SG_Soil_amm_per_chng_combined))
shapiro.test(resid(SG_Soil_amm_per_chng_combined))
#p-value = 0.707
anova(SG_Soil_amm_per_chng_combined, type=3)
#precip                       0.0210  1    16.8759 0.0004301 ***
#precip:soil_root             0.0099  2     3.9907 0.0325253 *  
#precip:soil_root:life_stage  0.0065  2     2.6043 0.0956028 .  

emmeans(SG_Soil_amm_per_chng_combined, pairwise~soil_root|precip)

SG_per_chng_inorg_N_soil_root_g=SG_inorg_N_trt_w_start %>% group_by(soil_root,precip)
soil_per_chng_amm_START_trt=summarise_at(SG_per_chng_inorg_N_soil_root_g, 
                                         "per_change_ug_N_NH4_g_dry_soil", list(~n(),~mean(.),sd,se=~sd(.)/sqrt(n())))
treatment_order=c("S.B","L.B","L.R")
(combine_per_chng_soil_ammonium_p=ggplot(soil_per_chng_amm_START_trt, aes(x=precip,y=mean, ymin = mean-se, ymax= mean+se,fill=factor(soil_root,levels = treatment_order)))+
    geom_bar(position="dodge",stat="identity", color="black",aes(fill=factor(soil_root,levels = treatment_order)))+
    geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+ylim(c(-1.1,0))+
    scale_fill_manual(values = c( "white","light gray", "dark grey"),labels=c("Sterile","Bulk","Rhizo"), limits=treatment_order, name=NULL)+
    scale_x_discrete(labels=c("Ambient","Drought"), name=NULL)+ylab("per Change in soil ammonium \n(Start to end)")+
    geom_text(aes(y=-.1, label=n),size=10,position=position_dodge(width=0.9))+theme_bw()+
    theme(axis.title = element_text(size = 23), axis.text = element_text(size = 23),
          legend.position = c(0.85,.9), legend.text=element_text(size=20),
          legend.background = element_rect(size=0.5,linetype="solid",colour ="black")))

#####Combined start nitrate and ammonium graph#####
ggarrange(combine_per_chng_soil_nitrate_p,combine_per_chng_soil_ammonium_p,ncol = 2,  legend = "none")
#15x7.38


#####Microbial Biomass####
SG_trt <- read.csv("D:/MERDS_2018/merds/Switchgrass/R_data/SG_trt_sheet.csv")
head(SG_trt)
SG_MBC_DOC <- read.csv("D:/MERDS_2018/merds/Switchgrass/R_data/SG_analyzed_MBC_DOC.csv")
head(SG_MBC_DOC)


SG_MBC_DOC_trt=merge(SG_MBC_DOC,SG_trt, by="Plant_Number",all.x = T)
SG_MBC_DOC_trt$soil_root=with(SG_MBC_DOC_trt, interaction(soil_status,root_association))



#S0il DOC
SG_Soil_DOC_mod= lm((unfumigated_DOC_concentration_ugC_gsoil)~soil_root+(1|block), 
                                  data= SG_MBC_DOC_trt)
qqPlot(resid(SG_Soil_DOC_mod))
hist(resid(SG_Soil_DOC_mod))
shapiro.test(resid(SG_Soil_DOC_mod))
#p-value =  0.265
anova(SG_Soil_DOC_mod, type=3)
#nada sig

#emmeans(SG_Soil_DOC_mod, pairwise~soil_root|precip)


soil_DOC_trt_SUM=SG_MBC_DOC_trt %>% group_by(soil_root,precip)%>%summarise_at(
                                         "unfumigated_DOC_concentration_ugC_gsoil", list(~n(),~mean(.),sd,se=~sd(.)/sqrt(n())))
treatment_order=c("S.B","L.B","L.R")
(soil_DOC_trt_p=ggplot(soil_DOC_trt_SUM, aes(x=factor(soil_root,levels = treatment_order),fill=factor(soil_root,levels = treatment_order)))+
    geom_bar(position="dodge",stat="identity", color="black",aes(y=mean))+
    geom_errorbar(aes(y=mean, ymin = mean-se, ymax= mean+se),position=position_dodge(width=0.9), width=0.2, size=1)+
    scale_fill_manual(values = c( "white","light gray", "dark grey"),labels=c("Sterile","Bulk","Rhizo"), limits=treatment_order, name=NULL)+
    scale_x_discrete(labels=c("Sterile","Bulk","Rhizo"), name=NULL)+ylab("Dissolved organic carbon (??g/g dry soil)")+
    geom_text(aes(y=3, label=n),size=10,position=position_dodge(width=0.9))+theme_bw()+
    theme(axis.title = element_text(size = 23), axis.text = element_text(size = 23),
          legend.position ="none",panel.grid.major = element_blank(), panel.grid.minor = element_blank()))



#S0il MBC
SG_Soil_MBC_mod= lm((ugC_MicBiomass_g_dry_soil)~soil_root+(1|block), 
                    data= SG_MBC_DOC_trt)
qqPlot(resid(SG_Soil_MBC_mod))
hist(resid(SG_Soil_MBC_mod))
shapiro.test(resid(SG_Soil_MBC_mod))
#p-value =  0.3948
anova(SG_Soil_MBC_mod, type=3)
#nada sig

#emmeans(SG_Soil_DOC_mod, pairwise~soil_root|precip)


soil_MBC_trt_SUM=SG_MBC_DOC_trt %>% group_by(soil_root,precip)%>%summarise_at(
  "ugC_MicBiomass_g_dry_soil", list(~n(),~mean(.),sd,se=~sd(.)/sqrt(n())))
treatment_order=c("S.B","L.B","L.R")
(soil_DOC_trt_p=ggplot(soil_MBC_trt_SUM, aes(x=factor(soil_root,levels = treatment_order),fill=factor(soil_root,levels = treatment_order)))+
    geom_bar(position="dodge",stat="identity", color="black",aes(y=mean))+
    geom_errorbar(aes(y=mean, ymin = mean-se, ymax= mean+se),position=position_dodge(width=0.9), width=0.2, size=1)+
    scale_fill_manual(values = c( "white","light gray", "dark grey"),labels=c("Sterile","Bulk","Rhizo"), limits=treatment_order, name=NULL)+
    scale_x_discrete(labels=c("Sterile","Bulk","Rhizo"), name=NULL)+ylab("Microbial biomass carbon (??g/g dry soil)")+
    geom_text(aes(y=20, label=n),size=10,position=position_dodge(width=0.9))+theme_bw()+
    theme(axis.title = element_text(size = 23), axis.text = element_text(size = 23),
          legend.position ="none",panel.grid.major = element_blank(), panel.grid.minor = element_blank()))




#######NEED TO UPDATE BELOW#########




library (emmeans)
emmeans(seed_survivor_model, pairwise ~soil_root)
boxplot(survivorship~soil_status, data= data_SG_seeds)
mean(data_SG_seeds["soil_status",]$survivorship)

biomass_summary =summary(data_SG_harvest_biomass)
str(biomass_summary)
write.csv(data_SG_harvest_biomass, "SG_biomass.csv")
#Exports data as a csv 

height_biomass = lm(total_biomass~height_cm+0, data=data_SG_germinates)

hist(resid(height_biomass))
qqPlot(resid(height_biomass))

summary(height_biomass)

plot(data_SG_germinates$height_cm, data_SG_germinates$total_biomass)
abline(height_biomass)

height_leafnum = data_SG_germinates$height_cm * data_SG_germinates$leaf_num
leafnum_height = lm(total_biomass~height_leafnum+0, data= data_SG_germinates)

hist(resid(leafnum_height))
qqPlot(resid(leafnum_height))

summary(leafnum_height)

data_SG_germinates = mutate(data_SG_germinates, Root_Shoot= root_weight_g / shoot_weight_g)

dataSG_height6_20 <- read.csv(file.choose())
colnames(dataSG_height6_20)
dataSG_height6_27 <- read.csv(file.choose())
colnames(dataSG_height6_27)
dataSG_height7_5 <- read.csv(file.choose())
colnames(dataSG_height7_5)
dataSG_height7_9 <- read.csv(file.choose())
colnames(dataSG_height7_9)

combined_heights = rbind(dataSG_height6_20, dataSG_height6_27, dataSG_height7_5, dataSG_height7_9)

data_SG_height_biomass = merge(dataSG_biomass, combined_heights, by = "Plant_Number", all = T)
data_SG_height_biomass = merge(data_SG_height_biomass, SG_trt, by = "Plant_Number", all = T)
data_SG_estbiomass = mutate(data_SG_height_biomass, est_biomass = 0.0109045*height_cm)

data_SG_estbiomass$soiltrt_rootass = with(data_SG_estbiomass, interaction(soil_status, root_association))
data_SG_estbiomass_germs = subset(data_SG_estbiomass, life_stage == "G")

biomass_time_model = lmer(est_biomass~date*soiltrt_rootass*precip+(1|block) + (1|Plant_Number), data = data_SG_estbiomass_germs)
anova(biomass_time_model, type=3)

write.csv(data_SG_estbiomass_germs, "SG_biomassovertime.csv")

#Potweight data

dataSG_potweight6_13 <- read.csv(file.choose())
colnames(dataSG_potweight6_13)
dataSG_potweight6_19 <- read.csv(file.choose())
colnames(dataSG_potweight6_19)
dataSG_potweight6_25 <- read.csv(file.choose())
colnames(dataSG_potweight6_25)
dataSG_potweight7_1 <- read.csv(file.choose())
colnames(dataSG_potweight7_1)
dataSG_potweight7_7 <- read.csv(file.choose())
colnames(dataSG_potweight7_7)

combined_potweight = rbind(dataSG_potweight6_13, dataSG_potweight6_19, dataSG_potweight6_25, dataSG_potweight7_1, dataSG_potweight7_7)

dataSG_potweight_trt = merge(combined_potweight, SG_trt, by = "Plant_Number", all = T)

dataSG_potweight_trt$soiltrt_rootass = with(dataSG_potweight_trt, interaction(soil_status, root_association))

potweight_time_model = lmer(potweight_g~date*precip*+(1|block) + (1|Plant_Number), data = dataSG_potweight_trt)
anova(potweight_time_model, type=3)



#MAKE A HISTOGRAM
#Displays the frequency distribution of any single variable you choose
#Example: What is the frequency distribution of fast plant height?
#You can alter the code by changing the name of objects
cover_scale<-data_SG_germinates$coverscale_1_5 #change name of the object here here
cover_scale

#this is not the same thing as calling the column as shown below
data_SG_germinates$coverscale_1_5
hist(cover_scale)
#What does the distribution look like to you?
#looks like a normalish distribution

#let look at QQplot, another way to check normality
qqPlot(cover_scale)
#"q-q plot is an exploratory graphical device used to check the validity of 
#a distributional assumption for a data set. In general, the basic idea is 
#to compute the theoretically expected value for each data point based on 
#the distribution in question." http://onlinestatbook.com/2/advanced_graphs/q-q_plots.html
#some of the samples fall outside the confidence intervals
#but it looks mostly normal

#Does transformation make it look more normal
#log transformation (natural log is the default)
hist(log(cover_scale))
qqPlot(log(height))
#NOT BETTER, SO MUCH WORSE

#Square-root transformation
hist(sqrt(cover_scale))
qqPlot(sqrt(height))
#better but I still think that the untransformed was better


#One way ANOVA = Analysis of Variance
#Example:Does fastplant height differ between our Live and sterile?
#In this case, treatments=soil_status (discrete)
#Our reponse=height_cm (continuous) and our data = dataFP_Tri_trt
#You can alter the code to fit your choice of treatments and response variables according to the questions

One_anova_Height<-aov(height_cm~soil_status, data = dataFP_Tri_trt)
#Look at the residuals to make sure they are normal
hist(resid(One_anova_Height))
qqPlot(resid(One_anova_Height))
#It all falls mostly within the confidence intervals so I think it is normal!

#you can also use Box Cox to determine what transformation you should use
boxCox(One_anova_Height)
#Since the confidence interval overlaps with 1, there is no need to transform


summary(One_anova_Height)#gives statistics for one-way anova

#Plot results as means with 95% confidence intervals
plot(dataFP_Tri_trt$soil_status,dataFP_Tri_trt$height_cm)

#Two-way ANOVA
#Example: How do soil status and inoculum affect fast plant height?
#In this case, treatments=soil_status (discrete) and inocul (discrete) 
#we are going include interactions
#why?
#Maybe, we think that the effects of SDS depend on the status of the soil 
#i.e. Live versus sterile
# "*"= signifies interactions while "+" is only main effects
#You can alter the code to fit your choice of treatments and response variables according to the questions


Two_anova_Height<-aov(height_cm~soil_status*inocul, data = dataFP_Tri_trt)
#Look at the residuals to make sure they are normal
hist(resid(Two_anova_Height))
qqPlot(resid(Two_anova_Height))
#It all falls mostly within the confidence intervals so I think it is normal!

#you can also use Box Cox to determine what transformation you should use
boxCox(Two_anova_Height)

summary(Two_anova_Height)#gives statistics for two-way anova
Two_anova_Height

#If we want to test the significance of the pairwise difference
TukeyHSD(Two_anova_Height)

#Plot results
#first we need to remove any blank cells
dataFP_Tri_trt_sub=subset(dataFP_Tri_trt, height_cm>0)

#plot using a boxplot
boxplot(height_cm~soil_status*inocul, data=dataFP_Tri_trt_sub, 
        ylab = "Fast plant height (cm)", xlab = "Inoculum x Soil Status")

#just inoculant effect 
boxplot(height_cm~inocul, data=dataFP_Tri_trt_sub)


#ANOVA = Analysis of Variance
#Example:Does fastplant height differ between our treatments?
#In this case, treatments=soil_status,block,precip,inocul (discrete)
#Our reponse=height_cm (continuous) and our data = dataFP_Tri_trt
#we are going include interactions
#why?
#Maybe, we think that the effects of SDS depend on the status of the soil 
#i.e. Live versus sterile
# "*"= signifies interactions while "+" is only main effects
#Our blocks are numbers but are they continous?
#They should probably be discrete? Use as.factor()

#You can alter the code to fit your choice of treatments and response variables according to the questions

anova_biomass<-aov(Root_Shoot~soil_status*precip*root_association+(1|block), data = data_SG_germinates)
#Look at the residuals to make sure they are normal
hist(resid(anova_biomass))
qqPlot(resid(anova_biomass))
#It mostly falls within the confidence intervals so I think it is normal!

#you can also use Box Cox to determine what transformation you should use
boxCox(anova_biomass)
#Since the confidence interval overlaps with 1, there is no need to transform


summary(anova_biomass)#gives statistics for one-way anova

#plot using a boxplot
boxplot(total_biomass~soil_status, data=data_SG_germinates)

#plot using a boxplot, but you can only choose two treatments
boxplot(total_biomass~precip*soil_status, data=data_SG_germinates)

#SIMPLE LINEAR REGRESSION
#Does the pot weight increase with Fast Plant height ?

#loading in the pot weights
FP_pot_weight <- read.csv(file.choose())
dataFP_Tri_PW=merge(dataFP_Tri_trt,FP_pot_weight, by="Plant_Number")

#In this case, treatments=potweight_g (continuous)
#Our reponse=height_cm (continuous) and our data = dataFP_Tri_PW

reg_Height<-lm(height_cm~potweight_g, data=dataFP_Tri_PW)
#Look at the residuals to make sure they are normal
hist(resid(reg_Height))
qqPlot(resid(reg_Height))
#It mostly falls within the confidence intervals so I think it is normal!

#you can also use Box Cox to determine what transformation you should use
boxCox(reg_Height)
#Since the confidence interval overlaps with 1, there is no need to transform

summary(reg_Height) #this will give you the statistics for the regression
#Run the next 2 lines together to plot the results
plot(dataFP_Tri_PW$potweight_g,dataFP_Tri_PW$height_cm) 
abline(reg_Height)
