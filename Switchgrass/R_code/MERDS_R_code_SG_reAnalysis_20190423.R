#DO THESE FIRST!!!
#Import dataset into R from your directory
#Here is my directory
dataSG_seed_surv <- read.csv("D:/MERDS_2018/merds/Switchgrass/R_data/SG_surv_Seed_germ.csv")
dataSG_seedling_surv <- read.csv("D:/MERDS_2018/merds/Switchgrass/R_data/SG_surv_Seedling_replace.csv")
data_SG_biomass <- read.csv("D:/MERDS_2018/merds/Switchgrass/R_data/SG_TotalBiomass.csv")
SG_trt <- read.csv("D:/MERDS_2018/merds/Switchgrass/R_data/SG_trt_sheet.csv")
SG_height_combin<- read.csv("D:/MERDS_2018/merds/Switchgrass/R_data/SG_height_time_combine.csv", header = T)

#now we need to merge the data file with the trt file


summary(dataSG_seed_surv)
dataSG_seed_surv[,8:14]=NULL
summary(dataSG_seedling_surv)
dataSG_seedling_surv[,9:14]=NULL



dataSG_seedling_surv_trt=merge(dataSG_seedling_surv, SG_trt, by="Plant_Number", all.x = T)
dataSG_seed_surv_trt=merge(dataSG_seed_surv, SG_trt, by="Plant_Number", all.x = T)


summary(dataSG_seedling_surv_trt)
summary(dataSG_seed_surv_trt)





library(ggplot2)
library(car)
library(MASS)
library(dplyr)
library(emmeans)
library(lme4)
options(contrasts=c("contr.sum", "contr.poly"))
"%w/o%" <- function(x,y)!('%in%'(x,y))
#let's look at germination and number of germinates

#let's look at prop of pots with germinates through time

dataSG_seed_surv_trt$pot_w_germ=dataSG_seed_surv_trt$num_germinates
dataSG_seed_surv_trt$pot_w_germ[dataSG_seed_surv_trt$pot_w_germ>0]=1
summary(dataSG_seed_surv_trt)
hist(dataSG_seed_surv_trt$pot_w_germ)
dataSG_seed_surv_trt$soil_root=with(dataSG_seed_surv_trt, interaction(soil_status,root_association))

pot_germ_model= glmer(pot_w_germ~precip*soil_root*exp_days+as.factor(block)+(1|Plant_Number), data= dataSG_seed_surv_trt, family = binomial)
#Warning message:
#  In checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,  :
#                 Model failed to converge with max|grad| = 0.613639 (tol = 0.001, component 1)
qqPlot(resid(pot_germ_model))
hist(resid(pot_germ_model))

Anova(pot_germ_model, type=3)
emmeans(pot_germ_model, ~precip)
emmeans(pot_germ_model, pairwise~soil_root)

emmeans(pot_germ_model, pairwise~soil_root|precip)

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
#harvested seedlings plus final survival

dataSG_seed_surv_trt_g=group_by(dataSG_seedling_surv_trt, Plant_Number)
total_seedlings_harvested=summarise_at(dataSG_seed_surv_trt_g, "num_removed", sum)
colnames(total_seedlings_harvested)[2]="tot_num_removed"

fin_dataSG_seed_surv_trt=subset(dataSG_seed_surv_trt, date=="7_09")

fin_dataSG_seed_surv_trt$tot_num_germ=fin_dataSG_seed_surv_trt$num_germinates+total_seedlings_harvested$tot_num_removed

fin_dataSG_seed_surv_trt$germinated=fin_dataSG_seed_surv_trt$tot_num_germ
fin_dataSG_seed_surv_trt$germinated[fin_dataSG_seed_surv_trt$germinated>0]=1
fin_dataSG_seed_surv_trt$soil_root=with(fin_dataSG_seed_surv_trt, interaction(soil_status,root_association))

summary(fin_dataSG_seed_surv_trt)


hist(fin_dataSG_seed_surv_trt$tot_num_germ)


#pots with germination

seed_germ_model= glm(germinated~precip*soil_root+as.factor(block), data= fin_dataSG_seed_surv_trt, family = binomial)
qqPlot(resid(seed_germ_model))
hist(resid(seed_germ_model))

Anova(seed_germ_model, type=3)
Anova(seed_germ_model, type=2)
emmeans(seed_germ_model, ~precip)
emmeans(seed_germ_model, pairwise~soil_root)
emmeans(seed_germ_model, pairwise~soil_root|precip)

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


#now lets look at the total number of germinates given there was at least one germinate


fin_dataSG_seed_surv_trt_w_germ=subset(fin_dataSG_seed_surv_trt, germinated>0)

hist(fin_dataSG_seed_surv_trt_w_germ$tot_num_germ)

seed_num_germ_model= lm(tot_num_germ~precip*soil_root+as.factor(block), data= fin_dataSG_seed_surv_trt_w_germ)
qqPlot(resid(seed_num_germ_model))
hist(resid(seed_num_germ_model))

Anova(seed_num_germ_model, type=3)
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







###


data_SG_biomass_seed_surv_trt=merge(fin_dataSG_seed_surv_trt, data_SG_biomass, by="Plant_Number")



#now lets look at total biomass of germinated seedlings
fin_dataSG_biomass_seed_surv_trt_w_germ=subset(data_SG_biomass_seed_surv_trt, germinated>0)

summary(fin_dataSG_biomass_seed_surv_trt_w_germ)

seed_biomas_model= lm(log(total_biomass)~precip*soil_root+as.factor(block), data= fin_dataSG_biomass_seed_surv_trt_w_germ)
qqPlot(resid(seed_biomas_model))
hist(resid(seed_biomas_model))

Anova(seed_biomas_model, type=3)
#precip:soil_root  10.020  2   6.7518  0.003395 ** 
#soil_root          5.193  2   3.4989  0.041514 * 

emmeans(seed_biomas_model, ~precip)
emmeans(seed_biomas_model, pairwise~soil_root|precip)
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



fin_dataSG_biomass_seed_surv_trt_w_germ_soil_root_precip_g=fin_dataSG_biomass_seed_surv_trt_w_germ %>% group_by(soil_root,precip)
total_bio_soil_root_precip=summarise_at(fin_dataSG_biomass_seed_surv_trt_w_germ_soil_root_precip_g, 
                                       "total_biomass", funs(n(),mean,sd,se=sd(.)/sqrt(n())))

ggplot(total_bio_soil_root_precip, aes(x=precip,y=mean, ymin = mean-se, ymax= mean+se,fill=soil_root))+
  geom_bar(position="dodge",stat="identity", color="black",aes(fill=soil_root))+
  geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+
  scale_fill_manual(values = c("gray", "white", "black"),labels=c("Bulk","Sterile","Rhizo"))+
  scale_x_discrete(labels=c("Ambient","Drought"))+ylab("Total biomass given gemination")+
  geom_text(aes(y=mean+se+0.02, label=n),position=position_dodge(width=0.9))+theme_bw()



#Roots alone
summary(fin_dataSG_biomass_seed_surv_trt_w_germ)
seed_root_model= lm(log(root_weight_g)~precip*soil_root+as.factor(block), data= fin_dataSG_biomass_seed_surv_trt_w_germ)
qqPlot(resid(seed_root_model))
hist(resid(seed_root_model))

Anova(seed_root_model, type=3)
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



fin_dataSG_biomass_seed_surv_trt_w_germ_soil_root_precip_g=fin_dataSG_biomass_seed_surv_trt_w_germ %>% group_by(soil_root,precip)
root_soil_root_precip=summarise_at(fin_dataSG_biomass_seed_surv_trt_w_germ_soil_root_precip_g, 
                                        "root_weight_g", funs(n(),mean,sd,se=sd(.)/sqrt(n())))

ggplot(root_soil_root_precip, aes(x=precip,y=mean, ymin = mean-se, ymax= mean+se,fill=soil_root))+
  geom_bar(position="dodge",stat="identity", color="black",aes(fill=soil_root))+
  geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+
  scale_fill_manual(values = c("gray", "white", "black"),labels=c("Bulk","Sterile","Rhizo"))+
  scale_x_discrete(labels=c("Ambient","Drought"))+ylab("Root biomass given gemination")+
  geom_text(aes(y=mean+se+0.02, label=n),position=position_dodge(width=0.9))+theme_bw()




#let's remove the sandy roots

#now lets look at total biomass of germinated seedlings
fin_dataSG_biomass_seed_surv_trt_w_germ=subset(data_SG_biomass_seed_surv_trt, germinated>0)
fin_dataSG_biomass_seed_surv_trt_w_germ_N_sandy=subset(fin_dataSG_biomass_seed_surv_trt_w_germ, sandy_root=="N")
summary(fin_dataSG_biomass_seed_surv_trt_w_germ_N_sandy)

seed_biomas_model_N_sandy= lm(log(total_biomass)~precip*soil_root+as.factor(block), data= fin_dataSG_biomass_seed_surv_trt_w_germ_N_sandy)
qqPlot(resid(seed_biomas_model_N_sandy))
hist(resid(seed_biomas_model_N_sandy))

Anova(seed_biomas_model_N_sandy, type=3)
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



fin_dataSG_biomass_seed_surv_trt_N_sandy_w_germ_soil_root_precip_g=fin_dataSG_biomass_seed_surv_trt_w_germ_N_sandy %>% group_by(soil_root,precip)
total_bio_N_sandy_soil_root_precip=summarise_at(fin_dataSG_biomass_seed_surv_trt_N_sandy_w_germ_soil_root_precip_g, 
                                        "total_biomass", funs(n(),mean,sd,se=sd(.)/sqrt(n())))

ggplot(total_bio_N_sandy_soil_root_precip, aes(x=precip,y=mean, ymin = mean-se, ymax= mean+se,fill=soil_root))+
  geom_bar(position="dodge",stat="identity", color="black",aes(fill=soil_root))+
  geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+
  scale_fill_manual(values = c("gray", "white", "black"),labels=c("Bulk","Sterile","Rhizo"))+
  scale_x_discrete(labels=c("Ambient","Drought"))+ylab("Total biomass given gemination (with sandy samples removed)")+
  geom_text(aes(y=mean+0.04, label=n),position=position_dodge(width=0.9))+theme_bw()



#Roots alone
summary(fin_dataSG_biomass_seed_surv_trt_w_germ_N_sandy)
seed_root_model_N_sandy= lm(log(root_weight_g)~precip*soil_root+as.factor(block), data= fin_dataSG_biomass_seed_surv_trt_w_germ_N_sandy)
qqPlot(resid(seed_root_model_N_sandy))
hist(resid(seed_root_model_N_sandy))

Anova(seed_root_model_N_sandy, type=3)
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



fin_dataSG_biomass_seed_surv_trt_N_sandy_w_germ_soil_root_precip_g=fin_dataSG_biomass_seed_surv_trt_w_germ_N_sandy %>% group_by(soil_root,precip)
root_N_sandy_soil_root_precip=summarise_at(fin_dataSG_biomass_seed_surv_trt_N_sandy_w_germ_soil_root_precip_g, 
                                   "root_weight_g", funs(n(),mean,sd,se=sd(.)/sqrt(n())))

ggplot(root_N_sandy_soil_root_precip, aes(x=precip,y=mean, ymin = mean-se, ymax= mean+se,fill=soil_root))+
  geom_bar(position="dodge",stat="identity", color="black",aes(fill=soil_root))+
  geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+
  scale_fill_manual(values = c("gray", "white", "black"),labels=c("Bulk","Sterile","Rhizo"))+
  scale_x_discrete(labels=c("Ambient","Drought"))+ylab("Root biomass given gemination (with sandy samples removed)")+
  geom_text(aes(y=mean+0.04, label=n),position=position_dodge(width=0.9))+theme_bw()




#Shoots alone
summary(fin_dataSG_biomass_seed_surv_trt_w_germ)
seed_shoot_model= lm(log(shoot_weight_g)~precip*soil_root+as.factor(block), data= fin_dataSG_biomass_seed_surv_trt_w_germ)
qqPlot(resid(seed_shoot_model))
hist(resid(seed_shoot_model))

Anova(seed_shoot_model, type=3)
#precip             2.06  1   10.0215  0.003257 ** 
#soil_root          1.63  2    3.9621  0.028395 *  
#as.factor(block)   4.70  4    5.7232  0.001240 ** 
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



fin_dataSG_biomass_seed_surv_trt_w_germ_soil_root_precip_g=fin_dataSG_biomass_seed_surv_trt_w_germ %>% group_by(soil_root,precip)
shoot_soil_root_precip=summarise_at(fin_dataSG_biomass_seed_surv_trt_w_germ_soil_root_precip_g, 
                                   "shoot_weight_g", funs(n(),mean,sd,se=sd(.)/sqrt(n())))

ggplot(shoot_soil_root_precip, aes(x=precip,y=mean, ymin = mean-se, ymax= mean+se,fill=soil_root))+
  geom_bar(position="dodge",stat="identity", color="black",aes(fill=soil_root))+
  geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+
  scale_fill_manual(values = c("gray", "white", "black"),labels=c("Bulk","Sterile","Rhizo"))+
  scale_x_discrete(labels=c("Ambient","Drought"))+ylab("Shoot biomass given gemination")+
  geom_text(aes(y=mean+se+0.001, label=n),position=position_dodge(width=0.9))+theme_bw()




#let's look at the transplanted seedlings

#now I want to see if surv in the first week of the experiement is affect by treatments
summary(dataSG_seedling_surv_trt)
dataSG_seedling_surv_trt_g=group_by(dataSG_seedling_surv_trt, Plant_Number)
total_seedlings_transplanted=summarise_at(dataSG_seedling_surv_trt_g, "Transplant_height_cm", sum,na.rm = TRUE)
colnames(total_seedlings_transplanted)[2]="tot_Transplant_height_cm"
summary(total_seedlings_transplanted)
nrow(total_seedlings_transplanted)

fin_dataSG_seedling_surv_trt=subset(dataSG_seedling_surv_trt, date=="7_9")
nrow(fin_dataSG_seedling_surv_trt)
fin_dataSG_seedling_surv_trt_merge=merge(fin_dataSG_seedling_surv_trt,total_seedlings_transplanted, by="Plant_Number")
nrow(fin_dataSG_seedling_surv_trt_merge)
fin_dataSG_seedling_surv_trt_merge$int_surv=fin_dataSG_seedling_surv_trt_merge$tot_Transplant_height_cm
fin_dataSG_seedling_surv_trt_merge$int_surv[fin_dataSG_seedling_surv_trt_merge$int_surv==0]=1
fin_dataSG_seedling_surv_trt_merge$int_surv[fin_dataSG_seedling_surv_trt_merge$int_surv!=1]=0
fin_dataSG_seedling_surv_trt_merge$soil_root=with(fin_dataSG_seedling_surv_trt_merge, interaction(soil_status,root_association))

summary(fin_dataSG_seedling_surv_trt_merge)


hist(fin_dataSG_seedling_surv_trt_merge$int_surv)

int_surv_model= glm(int_surv~precip*soil_root+as.factor(block), data= fin_dataSG_seedling_surv_trt_merge, family = binomial)
qqPlot(resid(int_surv_model))
hist(resid(int_surv_model))

Anova(int_surv_model, type=3)
emmeans(int_surv_model, ~precip)
emmeans(int_surv_model, pairwise~soil_root)
emmeans(int_surv_model, pairwise~precip|soil_root)

fin_fin_dataSG_seedling_surv_trt_soil_root_precip_g=fin_dataSG_seedling_surv_trt_merge %>% group_by(soil_root,precip)
int_surv_soil_root_precip=summarise_at(fin_fin_dataSG_seedling_surv_trt_soil_root_precip_g, "int_surv", mean)

ggplot(int_surv_soil_root_precip, aes(precip,int_surv))+geom_bar(position="dodge",stat="identity", color="black",aes(fill=soil_root))+
  scale_fill_manual(values = c("gray", "white", "black"),labels=c("Bulk","Sterile","Rhizo"))+
  scale_x_discrete(labels=c("Ambient","Drought"))+ylab("First week Transplant Survival")


summary(fin_dataSG_seedling_surv_trt)
summary(data_SG_biomass)
data_SG_biomass_seedling_surv_trt=merge(fin_dataSG_seedling_surv_trt, data_SG_biomass, by="Plant_Number", all.x = T)
data_SG_biomass_seedling_surv_trt$soil_root=with(data_SG_biomass_seedling_surv_trt, interaction(soil_status,root_association))
summary(data_SG_biomass_seedling_surv_trt)

#now lets look at biomass with dead as zero

data_SG_biomass_seedling_surv_trt$total_biomass[data_SG_biomass_seedling_surv_trt$total_biomass==0]=NA
data_SG_biomass_seedling_surv_trt$total_biomass_0=data_SG_biomass_seedling_surv_trt$total_biomass
data_SG_biomass_seedling_surv_trt$total_biomass_0[is.na(data_SG_biomass_seedling_surv_trt$total_biomass_0)]=0

data_SG_biomass_seedling_surv_trt$root_weight_g_0=data_SG_biomass_seedling_surv_trt$root_weight_g
data_SG_biomass_seedling_surv_trt$root_weight_g_0[is.na(data_SG_biomass_seedling_surv_trt$root_weight_g_0)]=0

data_SG_biomass_seedling_surv_trt$shoot_weight_g_0=data_SG_biomass_seedling_surv_trt$shoot_weight_g
data_SG_biomass_seedling_surv_trt$shoot_weight_g_0[is.na(data_SG_biomass_seedling_surv_trt$shoot_weight_g_0)]=0

summary(data_SG_biomass_seedling_surv_trt)


#now lets look at total biomass



seedling_biomas_model= lm(log(total_biomass_0+.001)~precip*soil_root+as.factor(block), data= data_SG_biomass_seedling_surv_trt)
qqPlot(resid(seedling_biomas_model))
hist(resid(seedling_biomas_model))

Anova(seedling_biomas_model, type=3)
#precip            31.95  1  13.4201 0.0004457 ***
#soil_root         23.30  2   4.8939 0.0098833 ** 
#precip:soil_root  15.88  2   3.3353 0.0406192 *  

emmeans(seedling_biomas_model, pairwise~soil_root)
emmeans(seedling_biomas_model, pairwise~soil_root|precip)
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



data_SG_biomass_seedling_surv_trt_soil_root_precip_g=data_SG_biomass_seedling_surv_trt %>% group_by(soil_root,precip)
seedling_total_bio_soil_root_precip=summarise_at(data_SG_biomass_seedling_surv_trt_soil_root_precip_g, 
                                        "total_biomass_0", funs(n(),mean,sd,se=sd(.)/sqrt(n())))

ggplot(seedling_total_bio_soil_root_precip, aes(x=precip,y=mean, ymin = mean-se, ymax= mean+se,fill=soil_root))+
  geom_bar(position="dodge",stat="identity", color="black",aes(fill=soil_root))+
  geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+
  scale_fill_manual(values = c("gray", "white", "black"),labels=c("Bulk","Sterile","Rhizo"))+
  scale_x_discrete(labels=c("Ambient","Drought"))+ylab("Total biomass transplants")+
  geom_text(aes(y=mean+se+0.02, label=n),position=position_dodge(width=0.9))+theme_bw()



#Roots alone
summary(data_SG_biomass_seedling_surv_trt)
seedling_root_model= lm(log(root_weight_g_0+.001)~precip*soil_root+as.factor(block), data= data_SG_biomass_seedling_surv_trt)
qqPlot(resid(seedling_root_model))
hist(resid(seedling_root_model))

Anova(seedling_root_model, type=3)
#precip             24.14  1   8.6958 0.00418 ** 
#soil_root          20.25  2   3.6482 0.03046 * 
#precip:soil_root   16.45  2   2.9639 0.05731 .  


emmeans(seedling_root_model, pairwise~soil_root)
emmeans(seedling_root_model, pairwise~soil_root|precip)
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



data_SG_biomass_seedling_surv_trt_soil_root_precip_g=data_SG_biomass_seedling_surv_trt %>% group_by(soil_root,precip)
seedling_root_soil_root_precip=summarise_at(data_SG_biomass_seedling_surv_trt_soil_root_precip_g, 
                                   "root_weight_g_0", funs(n(),mean,sd,se=sd(.)/sqrt(n())))

ggplot(seedling_root_soil_root_precip, aes(x=precip,y=mean, ymin = mean-se, ymax= mean+se,fill=soil_root))+
  geom_bar(position="dodge",stat="identity", color="black",aes(fill=soil_root))+
  geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+
  scale_fill_manual(values = c("gray", "white", "black"),labels=c("Bulk","Sterile","Rhizo"))+
  scale_x_discrete(labels=c("Ambient","Drought"))+ylab("Root biomass transplants")+
  geom_text(aes(y=mean+se+0.02, label=n),position=position_dodge(width=0.9))+theme_bw()

#with sandy roots removed
colnames(data_SG_biomass_seedling_surv_trt)
data_SG_biomass_seedling_surv_trt_N_sandy=subset(data_SG_biomass_seedling_surv_trt, sandy_root=="N")

#now lets look at total biomass



seedling_biomas_model_N_sandy= lm(log(total_biomass_0+0.001)~precip*soil_root+as.factor(block), data= data_SG_biomass_seedling_surv_trt_N_sandy)
qqPlot(resid(seedling_biomas_model_N_sandy))
hist(resid(seedling_biomas_model_N_sandy))

Anova(seedling_biomas_model_N_sandy, type=3)
#precip            19.52  1   9.0147  0.003623 ** 
#soil_root         26.57  2   6.1367  0.003384 **  

emmeans(seedling_biomas_model_N_sandy, pairwise~soil_root)
emmeans(seedling_biomas_model_N_sandy, pairwise~soil_root|precip)
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



data_SG_biomass_seedling_surv_trt_N_sandy_soil_root_precip_g=data_SG_biomass_seedling_surv_trt_N_sandy %>% group_by(soil_root,precip)
seedling_total_bio_N_sandy_soil_root_precip=summarise_at(data_SG_biomass_seedling_surv_trt_N_sandy_soil_root_precip_g, 
                                                 "total_biomass_0", funs(n(),mean,sd,se=sd(.)/sqrt(n())))

ggplot(seedling_total_bio_N_sandy_soil_root_precip, aes(x=precip,y=mean, ymin = mean-se, ymax= mean+se,fill=soil_root))+
  geom_bar(position="dodge",stat="identity", color="black",aes(fill=soil_root))+
  geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+
  scale_fill_manual(values = c("gray", "white", "black"),labels=c("Bulk","Sterile","Rhizo"))+
  scale_x_discrete(labels=c("Ambient","Drought"))+ylab("Total biomass transplants (with sandy samples removed)")+
  geom_text(aes(y=mean+se+0.02, label=n),position=position_dodge(width=0.9))+theme_bw()



#Roots alone
summary(data_SG_biomass_seedling_surv_trt_N_sandy)
seedling_root_model_N_sandy= lm(log(root_weight_g_0+.001)~precip*soil_root+as.factor(block), data= data_SG_biomass_seedling_surv_trt_N_sandy)
qqPlot(resid(seedling_root_model_N_sandy))
hist(resid(seedling_root_model_N_sandy))

Anova(seedling_root_model_N_sandy, type=3)
#precip             14.10  1   5.3558 0.02336 *  
#soil_root          22.21  2   4.2191 0.01830 *  


emmeans(seedling_root_model_N_sandy, pairwise~soil_root)
emmeans(seedling_root_model_N_sandy, pairwise~soil_root|precip)
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



data_SG_biomass_seedling_surv_trt_N_sandy_soil_root_precip_g=data_SG_biomass_seedling_surv_trt_N_sandy %>% group_by(soil_root,precip)
seedling_root_N_sandy_soil_root_precip=summarise_at(data_SG_biomass_seedling_surv_trt_N_sandy_soil_root_precip_g, 
                                            "root_weight_g_0", funs(n(),mean,sd,se=sd(.)/sqrt(n())))

ggplot(seedling_root_N_sandy_soil_root_precip, aes(x=precip,y=mean, ymin = mean-se, ymax= mean+se,fill=soil_root))+
  geom_bar(position="dodge",stat="identity", color="black",aes(fill=soil_root))+
  geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+
  scale_fill_manual(values = c("gray", "white", "black"),labels=c("Bulk","Sterile","Rhizo"))+
  scale_x_discrete(labels=c("Ambient","Drought"))+ylab("Root biomass transplants (with sandy samples removed)")+
  geom_text(aes(y=mean+se+0.02, label=n),position=position_dodge(width=0.9))+theme_bw()




#Shoots alone
summary(data_SG_biomass_seedling_surv_trt)
seedling_shoot_model= lm(sqrt(shoot_weight_g_0+.001)~precip*soil_root+as.factor(block), data= data_SG_biomass_seedling_surv_trt)
qqPlot(resid(seedling_shoot_model))
hist(resid(seedling_shoot_model))

Anova(seedling_shoot_model, type=3)
#precip           0.08441  1   67.5994 2.970e-12 ***
#soil_root        0.15209  2   60.8972 < 2.2e-16 ***
#precip:soil_root 0.05837  2   23.3708 1.016e-08 ***

emmeans(seedling_shoot_model, ~precip)
emmeans(seedling_shoot_model, pairwise~soil_root)
emmeans(seedling_shoot_model, pairwise~soil_root|precip)
#precip = A:
# L.B - S.B  0.14946329 0.2537568 34   0.589  0.8269
#L.B - L.R -0.25671123 0.1956448 34  -1.312  0.3983
#S.B - L.R -0.40617452 0.2432545 34  -1.670  0.2314

#precip = D:
#  contrast     estimate        SE df t.ratio p.value
#L.B - S.B -1.41470888 0.3800696 34  -3.722  0.0020
#L.B - L.R  0.08175893 0.2424079 34   0.337  0.9393
#S.B - L.R  1.49646782 0.4097436 34   3.652  0.0024



data_SG_biomass_seedling_surv_trt_soil_root_precip_g=data_SG_biomass_seedling_surv_trt %>% group_by(soil_root,precip)
seedling_shoot_soil_root_precip=summarise_at(data_SG_biomass_seedling_surv_trt_soil_root_precip_g, 
                                    "shoot_weight_g_0", funs(n(),mean,sd,se=sd(.)/sqrt(n())))

ggplot(seedling_shoot_soil_root_precip, aes(x=precip,y=mean, ymin = mean-se, ymax= mean+se,fill=soil_root))+
  geom_bar(position="dodge",stat="identity", color="black",aes(fill=soil_root))+
  geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+
  scale_fill_manual(values = c("gray", "white", "black"),labels=c("Bulk","Sterile","Rhizo"))+
  scale_x_discrete(labels=c("Ambient","Drought"))+ylab("Shoot biomass transplants")+
  geom_text(aes(y=mean+se+0.005, label=n),position=position_dodge(width=0.9))+theme_bw()


#filter out the dead seedlings

data_SG_biomass_seedling_surv_trt_na=subset(data_SG_biomass_seedling_surv_trt,total_biomass!=0)
summary(data_SG_biomass_seedling_surv_trt_na)

#now lets look at total biomass



seedling_biomas_model_na= lm(log(total_biomass_0)~precip*soil_root+as.factor(block), data= data_SG_biomass_seedling_surv_trt_na)
qqPlot(resid(seedling_biomas_model_na))
hist(resid(seedling_biomas_model_na))

Anova(seedling_biomas_model_na, type=3)
#soil_root         27.89  2  11.2867 5.642e-05 ***
#precip:soil_root   8.41  2   3.4050    0.0388 *  

emmeans(seedling_biomas_model_na, pairwise~soil_root)
emmeans(seedling_biomas_model_na, pairwise~soil_root|precip)
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



data_SG_biomass_seedling_surv_trt_na_soil_root_precip_g=data_SG_biomass_seedling_surv_trt_na %>% group_by(soil_root,precip)
seedling_total_bio_soil_root_precip_na=summarise_at(data_SG_biomass_seedling_surv_trt_na_soil_root_precip_g, 
                                                 "total_biomass_0", funs(n(),mean,sd,se=sd(.)/sqrt(n())))

ggplot(seedling_total_bio_soil_root_precip_na, aes(x=precip,y=mean, ymin = mean-se, ymax= mean+se,fill=soil_root))+
  geom_bar(position="dodge",stat="identity", color="black",aes(fill=soil_root))+
  geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+
  scale_fill_manual(values = c("gray", "white", "black"),labels=c("Bulk","Sterile","Rhizo"))+
  scale_x_discrete(labels=c("Ambient","Drought"))+ylab("Total biomass transplants (na removed)")+
  geom_text(aes(y=mean+se+0.02, label=n),position=position_dodge(width=0.9))+theme_bw()



#Roots alone
summary(data_SG_biomass_seedling_surv_trt_na)
seedling_root_model_na= lm(sqrt(root_weight_g_0)~precip*soil_root+as.factor(block), data= data_SG_biomass_seedling_surv_trt_na)
qqPlot(resid(seedling_root_model_na))
hist(resid(seedling_root_model_na))

Anova(seedling_root_model_na, type=3)
#soil_root        0.5631  2   8.1806 0.0006416 ***
#precip:soil_root 0.2938  2   4.2685 0.0178163 *  


emmeans(seedling_root_model_na, pairwise~soil_root)
emmeans(seedling_root_model_na, pairwise~soil_root|precip)
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



data_SG_biomass_seedling_surv_trt_na_soil_root_precip_g=data_SG_biomass_seedling_surv_trt_na %>% group_by(soil_root,precip)
seedling_root_soil_root_precip_na=summarise_at(data_SG_biomass_seedling_surv_trt_na_soil_root_precip_g, 
                                            "root_weight_g_0", funs(n(),mean,sd,se=sd(.)/sqrt(n())))

ggplot(seedling_root_soil_root_precip_na, aes(x=precip,y=mean, ymin = mean-se, ymax= mean+se,fill=soil_root))+
  geom_bar(position="dodge",stat="identity", color="black",aes(fill=soil_root))+
  geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+
  scale_fill_manual(values = c("gray", "white", "black"),labels=c("Bulk","Sterile","Rhizo"))+
  scale_x_discrete(labels=c("Ambient","Drought"))+ylab("Root biomass transplants (na removed)")+
  geom_text(aes(y=mean+se+0.02, label=n),position=position_dodge(width=0.9))+theme_bw()



#Shoots alone
summary(data_SG_biomass_seedling_surv_trt_na)
seedling_shoot_model_na= lm(sqrt(shoot_weight_g_0)~precip*soil_root+as.factor(block), data= data_SG_biomass_seedling_surv_trt_na)
qqPlot(resid(seedling_shoot_model_na))
hist(resid(seedling_shoot_model_na))

Anova(seedling_shoot_model_na, type=3)
#precip           0.03748  1   46.7229 2.489e-09 ***
#soil_root        0.16183  2  100.8740 < 2.2e-16 ***
#precip:soil_root 0.03762  2   23.4478 1.605e-08 ***

emmeans(seedling_shoot_model_na, ~precip)
emmeans(seedling_shoot_model_na, pairwise~soil_root)
emmeans(seedling_shoot_model_na, pairwise~soil_root|precip)
#precip = A:
# L.B - S.B  0.14946329 0.2537568 34   0.589  0.8269
#L.B - L.R -0.25671123 0.1956448 34  -1.312  0.3983
#S.B - L.R -0.40617452 0.2432545 34  -1.670  0.2314

#precip = D:
#  contrast     estimate        SE df t.ratio p.value
#L.B - S.B -1.41470888 0.3800696 34  -3.722  0.0020
#L.B - L.R  0.08175893 0.2424079 34   0.337  0.9393
#S.B - L.R  1.49646782 0.4097436 34   3.652  0.0024



data_SG_biomass_seedling_surv_trt_na_soil_root_precip_g=data_SG_biomass_seedling_surv_trt_na %>% group_by(soil_root,precip)
seedling_shoot_soil_root_precip_na=summarise_at(data_SG_biomass_seedling_surv_trt_na_soil_root_precip_g, 
                                             "shoot_weight_g_0", funs(n(),mean,sd,se=sd(.)/sqrt(n())))

ggplot(seedling_shoot_soil_root_precip_na, aes(x=precip,y=mean, ymin = mean-se, ymax= mean+se,fill=soil_root))+
  geom_bar(position="dodge",stat="identity", color="black",aes(fill=soil_root))+
  geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+
  scale_fill_manual(values = c("gray", "white", "black"),labels=c("Bulk","Sterile","Rhizo"))+
  scale_x_discrete(labels=c("Ambient","Drought"))+ylab("Shoot biomass transplants (na removed)")+
  geom_text(aes(y=mean+se+0.005, label=n),position=position_dodge(width=0.9))+theme_bw()


data_SG_biomass_seedling_surv_trt_na$root_shoot=data_SG_biomass_seedling_surv_trt_na$root_weight_g/data_SG_biomass_seedling_surv_trt_na$shoot_weight_g


summary(data_SG_biomass_seedling_surv_trt_na)

root_shoot_model_na= lm(log(root_shoot)~precip*soil_root+as.factor(block), data= data_SG_biomass_seedling_surv_trt_na)
qqPlot(resid(root_shoot_model_na))
hist(resid(root_shoot_model_na))

Anova(root_shoot_model_na, type=3)
#precip           0.03748  1   46.7229 2.489e-09 ***
#soil_root        0.16183  2  100.8740 < 2.2e-16 ***
#precip:soil_root 0.03762  2   23.4478 1.605e-08 ***

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



data_SG_biomass_seedling_surv_trt_na_soil_root_precip_g=data_SG_biomass_seedling_surv_trt_na %>% group_by(soil_root,precip)
seedling_root_shoot_soil_root_precip_na=summarise_at(data_SG_biomass_seedling_surv_trt_na_soil_root_precip_g, 
                                                "root_shoot", funs(n(),mean,sd,se=sd(.)/sqrt(n())))

ggplot(seedling_root_shoot_soil_root_precip_na, aes(x=precip,y=mean, ymin = mean-se, ymax= mean+se,fill=soil_root))+
  geom_bar(position="dodge",stat="identity", color="black",aes(fill=soil_root))+
  geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+
  scale_fill_manual(values = c("gray", "white", "black"),labels=c("Bulk","Sterile","Rhizo"))+
  scale_x_discrete(labels=c("Ambient","Drought"))+ylab("Root:shoot (na removed)")+
  geom_text(aes(y=mean+se+1, label=n),position=position_dodge(width=0.9))+theme_bw()


#remove the outlier
data_SG_biomass_seedling_surv_trt_na_out=subset(data_SG_biomass_seedling_surv_trt_na, root_shoot<50)
summary(data_SG_biomass_seedling_surv_trt_na_out)

root_shoot_model_na_out= lm(sqrt(root_shoot)~precip*soil_root+as.factor(block), data= data_SG_biomass_seedling_surv_trt_na_out)
qqPlot(resid(root_shoot_model_na_out))
hist(resid(root_shoot_model_na_out))

Anova(root_shoot_model_na_out, type=3)
#precip             4.33  1   3.0869 0.08336 . 

emmeans(root_shoot_model_na_out, ~precip)
emmeans(root_shoot_model_na_out, pairwise~soil_root)
emmeans(root_shoot_model_na_out, pairwise~soil_root|precip)




data_SG_biomass_seedling_surv_trt_na_soil_root_precip_out_g=data_SG_biomass_seedling_surv_trt_na_out %>% group_by(soil_root,precip)
seedling_root_shoot_soil_root_precip_na_out=summarise_at(data_SG_biomass_seedling_surv_trt_na_soil_root_precip_out_g, 
                                                     "root_shoot", funs(n(),mean,sd,se=sd(.)/sqrt(n())))

ggplot(seedling_root_shoot_soil_root_precip_na_out, aes(x=precip,y=mean, ymin = mean-se, ymax= mean+se,fill=soil_root))+
  geom_bar(position="dodge",stat="identity", color="black",aes(fill=soil_root))+
  geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+
  scale_fill_manual(values = c("gray", "white", "black"),labels=c("Bulk","Sterile","Rhizo"))+
  scale_x_discrete(labels=c("Ambient","Drought"))+ylab("Root:shoot (na removed)")+
  geom_text(aes(y=mean+se+1, label=n),position=position_dodge(width=0.9))+theme_bw()


#Final survival in the drought treatment only since all pots had live seedlings in the ambient treatments

summary(data_SG_biomass_seedling_surv_trt)


data_SG_biomass_seedling_surv_trt$surv=data_SG_biomass_seedling_surv_trt$total_biomass
data_SG_biomass_seedling_surv_trt$surv[is.na(data_SG_biomass_seedling_surv_trt$surv)]=0
data_SG_biomass_seedling_surv_trt$surv[data_SG_biomass_seedling_surv_trt$surv!=0]=1

data_SG_biomass_seedling_surv_drought=subset(data_SG_biomass_seedling_surv_trt, precip=="D")

#survival of the transplant

drought_surv_model= glm(surv~soil_root+as.factor(block), data= data_SG_biomass_seedling_surv_drought, family = binomial)
qqPlot(resid(drought_surv_model))
hist(resid(drought_surv_model))

Anova(drought_surv_model, type=3)


emmeans(drought_surv_model, pairwise~soil_root)


#Bars by soil_root

data_SG_biomass_seedling_surv_drought_g=group_by(data_SG_biomass_seedling_surv_drought, soil_root)
drought_surv_soil_root=summarise_at(data_SG_biomass_seedling_surv_drought_g, "surv", mean)

ggplot(drought_surv_soil_root, aes(soil_root,surv))+geom_bar(stat="identity", color="black",aes(fill=soil_root))+
  scale_fill_manual(values = c("gray", "white", "black"),labels=c("Bulk","Sterile","Rhizo"))+
  scale_x_discrete(limits=c("S.B", "L.B", "L.R"),labels=c("Sterile","Bulk","Rhizo"))+ylab("Survival under drought treatment")


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

delta_height_trans_model= lm((delta_height)~precip*soil_root+as.factor(block), data= transplant_SG_height_plant_harvest_trt)
qqPlot(resid(delta_height_trans_model))
hist(resid(delta_height_trans_model))

Anova(delta_height_trans_model, type=3)
#precip            770.3  1  64.4430 7.140e-12 ***
#soil_root        1074.7  2  44.9547 8.220e-14 ***
#as.factor(block)  151.0  4   3.1587    0.0183 *  
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

delta_height_trans_model_0= lm((delta_height)~precip*soil_root+as.factor(block), data= transplant_SG_height_plant_harvest_trt_0)
qqPlot(resid(delta_height_trans_model_0))
hist(resid(delta_height_trans_model_0))

Anova(delta_height_trans_model_0, type=3)
#precip            339.9  1  47.7230 2.211e-09 ***
#soil_root        1210.5  2  84.9733 < 2.2e-16 ***
#as.factor(block)   66.5  4   2.3352   0.06437 .  
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

hist(stdres(height_biomass))
qqPlot(stdres(height_biomass))

summary(height_biomass)

plot(data_SG_germinates$height_cm, data_SG_germinates$total_biomass)
abline(height_biomass)

height_leafnum = data_SG_germinates$height_cm * data_SG_germinates$leaf_num
leafnum_height = lm(total_biomass~height_leafnum+0, data= data_SG_germinates)

hist(stdres(leafnum_height))
qqPlot(stdres(leafnum_height))

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

biomass_time_model = lmer(est_biomass~date*soiltrt_rootass*precip+as.factor(block) + (1|Plant_Number), data = data_SG_estbiomass_germs)
Anova(biomass_time_model, type=3)

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

potweight_time_model = lmer(potweight_g~date*precip*+as.factor(block) + (1|Plant_Number), data = dataSG_potweight_trt)
Anova(potweight_time_model, type=3)



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
hist(stdres(One_anova_Height))
qqPlot(stdres(One_anova_Height))
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
hist(stdres(Two_anova_Height))
qqPlot(stdres(Two_anova_Height))
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

anova_biomass<-aov(Root_Shoot~soil_status*precip*root_association+as.factor(block), data = data_SG_germinates)
#Look at the residuals to make sure they are normal
hist(stdres(anova_biomass))
qqPlot(stdres(anova_biomass))
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
hist(stdres(reg_Height))
qqPlot(stdres(reg_Height))
#It mostly falls within the confidence intervals so I think it is normal!

#you can also use Box Cox to determine what transformation you should use
boxCox(reg_Height)
#Since the confidence interval overlaps with 1, there is no need to transform

summary(reg_Height) #this will give you the statistics for the regression
#Run the next 2 lines together to plot the results
plot(dataFP_Tri_PW$potweight_g,dataFP_Tri_PW$height_cm) 
abline(reg_Height)
