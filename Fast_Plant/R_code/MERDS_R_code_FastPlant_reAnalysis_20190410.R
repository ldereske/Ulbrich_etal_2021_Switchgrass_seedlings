#DO THESE FIRST!!!
#Import dataset into R from your directory
#Here is my directory

data_FP_biomass <- read.csv("D:/MERDS_2018/merds/Fast_Plant/R_data/Fast_plant_harvest_weight_R.csv")
FP_trt <- read.csv("D:/MERDS_2018/merds/Fast_Plant/R_data/FastPlant_trt.csv")
#now we need to merge the data file with the trt file


summary(data_FP_biomass)
summary(FP_trt)




data_FP_biomass_trt=merge(data_FP_biomass, FP_trt, by="Plant_Number", all.y = T)



summary(data_FP_biomass_trt)






library(ggplot2)
library(car)
library(MASS)
library(dplyr)
library(emmeans)
library(lme4)
options(contrasts=c("contr.sum", "contr.poly"))


#now lets look at total biomass biomass

#make a total biomass column
colnames(data_FP_biomass_trt)

data_FP_biomass_trt$tot_bio=data_FP_biomass_trt$root_weight+data_FP_biomass_trt$shoot_weight

FP_tot_bio_model= lm(sqrt(tot_bio+0.0001)~precip*soil_status*inocul+as.factor(block), data= data_FP_biomass_trt)
qqPlot(resid(FP_tot_bio_model))
hist(resid(FP_tot_bio_model))

Anova(FP_tot_bio_model, type=3,singular.ok =T)
#as.factor(block)          0.18679   4   3.5696 0.007555 **
#precip:inocul             0.09764   3   2.4878 0.061185 . 


emmeans(FP_tot_bio_model, pairwise~precip|inocul)
#$contrasts
#inocul = BI:
#  contrast estimate     SE  df t.ratio p.value
#A - D     -0.0675 0.0303 235 -2.226  0.0270 

#inocul = GI:
#  contrast estimate     SE  df t.ratio p.value
#A - D     -0.0488 0.0298 235 -1.638  0.1028 

#inocul = PM:
#  contrast estimate     SE  df t.ratio p.value
#A - D      0.0306 0.0310 235  0.986  0.3253 

#inocul = S:
#  contrast estimate     SE  df t.ratio p.value
#A - D      nonEst     NA  NA     NA      NA 

#inocul = SDS:
#  contrast estimate     SE  df t.ratio p.value
#A - D      0.0441 0.0298 235  1.480  0.1403 

#Results are averaged over the levels of: soil_status, block 

FP_tot_bio_model.2= lm(sqrt(tot_bio+0.0001)~precip*inocul+as.factor(block), data= data_FP_biomass_trt)
qqPlot(resid(FP_tot_bio_model.2))
hist(resid(FP_tot_bio_model.2))

Anova(FP_tot_bio_model.2, type=3)
#precip           0.0672   1  4.7349 0.030521 *  
#inocul           0.2440   4  4.2996 0.002225 ** 
#as.factor(block) 0.1933   4  3.4073 0.009836 ** 
#precip:inocul    0.1455   4  2.5645 0.038962 *
Anova(FP_tot_bio_model.2, type=2,singular.ok =T)

emmeans(FP_tot_bio_model.2, pairwise~precip|inocul)
# contrast estimate     SE  df t.ratio p.value
#A - D     -0.0687 0.0316 243 -2.176  0.0305 


data_FP_biomass_trt_NA=data_FP_biomass_trt[complete.cases(data_FP_biomass_trt[,"tot_bio"]),]
data_FP_biomass_trt_NA_g=data_FP_biomass_trt_NA %>% group_by(precip,inocul)
tot_bio_data_FP_biomass_trt_sum=summarise_at(data_FP_biomass_trt_NA_g, 
                                       "tot_bio", funs(n(),mean,sd,se=sd(.)/sqrt(n())))
FP_Isolates_order=c("BI","GI","PM","SDS","S")
ggplot(tot_bio_data_FP_biomass_trt_sum, aes(x=inocul,y=mean, ymin = mean-se, ymax= mean+se,fill=precip))+
  geom_bar(position="dodge",stat="identity", color="black",aes(fill=precip))+
  geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+
  scale_x_discrete(labels=c("Bad Isolates","Good Isolates", "Product Microbes","SDS Product","Sterile"),limits=FP_Isolates_order)+
  ylab("Total biomass")+scale_fill_manual(values = c("navy", "red4"),labels=c("Ambient","Drought"))+
  geom_text(aes(y=mean+se+0.02, label=n),position=position_dodge(width=0.9))+theme_bw()

#  
#


#I want to see if the inculation interacts with the soil status
#need to remove the sterile treatments since it is not a part of the factorial
data_FP_biomass_trt_live=subset(data_FP_biomass_trt, inocul!="S")
FP_tot_bio_model.3= lm(sqrt(tot_bio+0.0001)~precip*soil_status*inocul+as.factor(block), data= data_FP_biomass_trt_live)
qqPlot(resid(FP_tot_bio_model.3))
hist(resid(FP_tot_bio_model.3))

Anova(FP_tot_bio_model.3, type=3)
#soil_status                0.2774   1   23.6161 2.299e-06 ***
#as.factor(block)           0.1544   4    3.2868   0.01223 *  
#precip:inocul              0.1344   3    3.8154   0.01083 * 


emmeans(FP_tot_bio_model.3, ~soil_status)


data_FP_biomass_trt_NA_g2=data_FP_biomass_trt_NA %>% group_by(precip,inocul,soil_status)
tot_bio_data_FP_biomass_trt_sum2=summarise_at(data_FP_biomass_trt_NA_g2, 
                                             "tot_bio", funs(n(),mean,sd,se=sd(.)/sqrt(n())))
FP_Isolates_order=c("BI","GI","PM","SDS","S")
ggplot(tot_bio_data_FP_biomass_trt_sum2, aes(x=inocul,y=mean, ymin = mean-se, ymax= mean+se,fill=interaction(precip,soil_status)))+
  geom_bar(position="dodge",stat="identity", color="black",aes(fill=interaction(precip,soil_status)))+
  geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+
  scale_x_discrete(labels=c("Bad Isolates","Good Isolates", "Product Microbes","SDS Product","Sterile"),limits=FP_Isolates_order)+
  ylab("Total biomass")+scale_fill_manual(values = c("navy", "red4","blue","red"),labels=c("Ambient-Live","Drought-Live","Ambient-Sterile","Drought-Sterile"))+
  geom_text(aes(y=mean+se+0.02, label=n),position=position_dodge(width=0.9))+theme_bw()



#Roots alone

FP_root_bio_model= lm(sqrt(root_weight+0.0001)~precip*soil_status*inocul+as.factor(block), data= data_FP_biomass_trt)
qqPlot(resid(FP_root_bio_model))
hist(resid(FP_root_bio_model))

Anova(FP_root_bio_model, type=3,singular.ok =T)
#as.factor(block)          0.1652   4   2.9449 0.021097 * 
#precip:inocul             0.1751   3   4.1629 0.006748 **


emmeans(FP_root_bio_model, pairwise~precip|inocul)


FP_root_bio_model.2= lm(sqrt(root_weight+0.0001)~precip*inocul+as.factor(block), data= data_FP_biomass_trt)
qqPlot(resid(FP_root_bio_model.2))
hist(resid(FP_root_bio_model.2))

Anova(FP_root_bio_model.2, type=3)
#precip            0.1016   1   6.9570  0.008888 ** 
#inocul            0.4931   4   8.4438 2.148e-06 ***
#as.factor(block)  0.1729   4   2.9615  0.020467 *  
#precip:inocul     0.1879   4   3.2181  0.013439 *  


emmeans(FP_root_bio_model.2, pairwise~precip|inocul)
# A - D    -0.09924 0.0320 243 -3.098  0.0022 

#inocul = GI:
#  contrast estimate     SE  df t.ratio p.value
#A - D    -0.07376 0.0315 243 -2.344  0.0199 


data_FP_biomass_trt_NA=data_FP_biomass_trt[complete.cases(data_FP_biomass_trt[,"tot_bio"]),]
data_FP_biomass_trt_NA_g=data_FP_biomass_trt_NA %>% group_by(precip,inocul)
root_bio_data_FP_biomass_trt_sum=summarise_at(data_FP_biomass_trt_NA_g, 
                                             "root_weight", funs(n(),mean,sd,se=sd(.)/sqrt(n())))
FP_Isolates_order=c("BI","GI","PM","SDS","S")
ggplot(root_bio_data_FP_biomass_trt_sum, aes(x=inocul,y=mean, ymin = mean-se, ymax= mean+se,fill=precip))+
  geom_bar(position="dodge",stat="identity", color="black",aes(fill=precip))+
  geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+
  scale_x_discrete(labels=c("Bad Isolates","Good Isolates", "Product Microbes","SDS Product","Sterile"),limits=FP_Isolates_order)+
  ylab("Root biomass")+scale_fill_manual(values = c("navy", "red4"),labels=c("Ambient","Drought"))+
  geom_text(aes(y=mean+se+0.02, label=n),position=position_dodge(width=0.9))+theme_bw()


#I want to see if the inculation interacts with the soil status
#need to remove the sterile treatments since it is not a part of the factorial
data_FP_biomass_trt_live=subset(data_FP_biomass_trt, inocul!="S")
FP_root_bio_model.3= lm(sqrt(root_weight+0.0001)~precip*soil_status*inocul+as.factor(block), data= data_FP_biomass_trt_live)
qqPlot(resid(FP_root_bio_model.3))
hist(resid(FP_root_bio_model.3))

Anova(FP_root_bio_model.3, type=3)
#precip                     0.0605   1   4.9586  0.027023 *  
#soil_status                0.1194   1   9.7865  0.002008 ** 
#as.factor(block)           0.1385   4   2.8375  0.025395 *  
#precip:inocul              0.1749   3   4.7802  0.003037 ** 

emmeans(FP_root_bio_model.3, ~soil_status)


data_FP_biomass_trt_NA_g2=data_FP_biomass_trt_NA %>% group_by(precip,inocul,soil_status)
root_bio_data_FP_biomass_trt_sum2=summarise_at(data_FP_biomass_trt_NA_g2, 
                                              "root_weight", funs(n(),mean,sd,se=sd(.)/sqrt(n())))
FP_Isolates_order=c("BI","GI","PM","SDS","S")
ggplot(root_bio_data_FP_biomass_trt_sum2, aes(x=inocul,y=mean, ymin = mean-se, ymax= mean+se,fill=interaction(precip,soil_status)))+
  geom_bar(position="dodge",stat="identity", color="black",aes(fill=interaction(precip,soil_status)))+
  geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+
  scale_x_discrete(labels=c("Bad Isolates","Good Isolates", "Product Microbes","SDS Product","Sterile"),limits=FP_Isolates_order)+
  ylab("Root biomass")+scale_fill_manual(values = c("navy", "red4","blue","red"),labels=c("Ambient-Live","Drought-Live","Ambient-Sterile","Drought-Sterile"))+
  geom_text(aes(y=mean+se+0.02, label=n),position=position_dodge(width=0.9))+theme_bw()



#Shoots alone
FP_shoot_bio_model= lm(sqrt(shoot_weight+0.0001)~precip*soil_status*inocul+as.factor(block), data= data_FP_biomass_trt)
qqPlot(resid(FP_shoot_bio_model))
hist(resid(FP_shoot_bio_model))

Anova(FP_shoot_bio_model, type=3,singular.ok =T)
#as.factor(block)          0.02562   4   3.2316 0.01319 *


emmeans(FP_shoot_bio_model, pairwise~precip|inocul)


FP_shoot_bio_model.2= lm(sqrt(shoot_weight+0.0001)~precip*inocul+as.factor(block), data= data_FP_biomass_trt)
qqPlot(resid(FP_shoot_bio_model.2))
hist(resid(FP_shoot_bio_model.2))

Anova(FP_shoot_bio_model.2, type=3)
#precip           0.0636   1   22.7127 3.233e-06 ***
#inocul           0.0352   4    3.1417   0.01523 *  
#as.factor(block) 0.0262   4    2.3395   0.05582 .  


emmeans(FP_shoot_bio_model.2, pairwise~precip|inocul)
#$contrasts
#inocul = BI:
#  contrast estimate     SE  df t.ratio p.value
#A - D      0.0324 0.0140 244 2.313   0.0216 

#inocul = GI:
#  contrast estimate     SE  df t.ratio p.value
#A - D      0.0277 0.0138 244 2.010   0.0455 


data_FP_biomass_trt_NA=data_FP_biomass_trt[complete.cases(data_FP_biomass_trt[,"tot_bio"]),]
data_FP_biomass_trt_NA_g=data_FP_biomass_trt_NA %>% group_by(precip,inocul)
shoot_bio_data_FP_biomass_trt_sum=summarise_at(data_FP_biomass_trt_NA_g, 
                                              "shoot_weight", funs(n(),mean,sd,se=sd(.)/sqrt(n())))
FP_Isolates_order=c("BI","GI","PM","SDS","S")
ggplot(shoot_bio_data_FP_biomass_trt_sum, aes(x=inocul,y=mean, ymin = mean-se, ymax= mean+se,fill=precip))+
  geom_bar(position="dodge",stat="identity", color="black",aes(fill=precip))+
  geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+
  scale_x_discrete(labels=c("Bad Isolates","Good Isolates", "Product Microbes","SDS Product","Sterile"),limits=FP_Isolates_order)+
  ylab("Shoot biomass")+scale_fill_manual(values = c("navy", "red4"),labels=c("Ambient","Drought"))+
  geom_text(aes(y=mean+se+0.002, label=n),position=position_dodge(width=0.9))+theme_bw()


#I want to see if the inculation interacts with the soil status
#need to remove the sterile treatments since it is not a part of the factorial
data_FP_biomass_trt_live=subset(data_FP_biomass_trt, inocul!="S")
FP_shoot_bio_model.3= lm(sqrt(shoot_weight+0.0001)~precip*soil_status*inocul+as.factor(block), data= data_FP_biomass_trt_live)
qqPlot(resid(FP_shoot_bio_model.3))
hist(resid(FP_shoot_bio_model.3))

Anova(FP_shoot_bio_model.3, type=3)
#precip                    0.0676   1   33.5768 2.469e-08 ***
#soil_status               0.1919   1   95.2530 < 2.2e-16 ***
#as.factor(block)          0.0215   4    2.6703   0.03323 *  
#precip:soil_status        0.0119   1    5.9212   0.01579 *  
emmeans(FP_root_bio_model.3, pairwise~precip|soil_status)


data_FP_biomass_trt_NA_g2=data_FP_biomass_trt_NA %>% group_by(precip,inocul,soil_status)
shoot_bio_data_FP_biomass_trt_sum2=summarise_at(data_FP_biomass_trt_NA_g2, 
                                               "shoot_weight", funs(n(),mean,sd,se=sd(.)/sqrt(n())))
FP_Isolates_order=c("BI","GI","PM","SDS","S")
ggplot(shoot_bio_data_FP_biomass_trt_sum2, aes(x=inocul,y=mean, ymin = mean-se, ymax= mean+se,fill=interaction(precip,soil_status)))+
  geom_bar(position="dodge",stat="identity", color="black",aes(fill=interaction(precip,soil_status)))+
  geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+
  scale_x_discrete(labels=c("Bad Isolates","Good Isolates", "Product Microbes","SDS Product","Sterile"),limits=FP_Isolates_order)+
  ylab("Shoot biomass")+scale_fill_manual(values = c("navy", "red4","blue","red"),labels=c("Ambient-Live","Drought-Live","Ambient-Sterile","Drought-Sterile"))+
  geom_text(aes(y=mean+se+0.002, label=n),position=position_dodge(width=0.9))+theme_bw()



#let's look at root biomass but with all of the samples flagged for sand removed
#
colnames(data_FP_biomass_trt)
data_FP_biomass_trt_no_sand=subset(data_FP_biomass_trt, sandy_root=="N")

FP_root_NS_bio_model= lm(sqrt(root_weight+0.0001)~precip*soil_status*inocul+as.factor(block), data= data_FP_biomass_trt_no_sand)
qqPlot(resid(FP_root_NS_bio_model))
hist(resid(FP_root_NS_bio_model))

Anova(FP_root_NS_bio_model, type=3,singular.ok =T)
#inocul                    0.11555   3   3.6266 0.01390 *
#as.factor(block)          0.11583   4   2.7264 0.03039 *


emmeans(FP_root_NS_bio_model, pairwise~precip|inocul)


FP_root_NS_bio_model.2= lm(sqrt(root_weight+0.0001)~precip*inocul+as.factor(block), data= data_FP_biomass_trt_no_sand)
qqPlot(resid(FP_root_NS_bio_model.2))
hist(resid(FP_root_NS_bio_model.2))

Anova(FP_root_NS_bio_model.2, type=3)
#precip           0.0749   1   6.6177 0.01076 *  
#inocul           0.1369   4   3.0236 0.01871 *  
#as.factor(block) 0.1298   4   2.8675 0.02412 *  


emmeans(FP_root_NS_bio_model.2, pairwise~precip|inocul)
#$contrasts
#inocul = BI:
#  contrast estimate     SE  df t.ratio p.value
#A - D    -0.05398 0.0299 217 -1.807  0.0722 

#inocul = GI:
#  contrast estimate     SE  df t.ratio p.value
#A - D    -0.07404 0.0277 217 -2.672  0.0081 


data_FP_biomass_trt_no_sand_NA=data_FP_biomass_trt_no_sand[complete.cases(data_FP_biomass_trt_no_sand[,"tot_bio"]),]
data_FP_biomass_trt_no_sand_NA_g=data_FP_biomass_trt_no_sand_NA %>% group_by(precip,inocul)
root_data_FP_biomass_trt_no_sand_sum=summarise_at(data_FP_biomass_trt_no_sand_NA_g, 
                                              "root_weight", funs(n(),mean,sd,se=sd(.)/sqrt(n())))
FP_Isolates_order=c("BI","GI","PM","SDS","S")
ggplot(root_data_FP_biomass_trt_no_sand_sum, aes(x=inocul,y=mean, ymin = mean-se, ymax= mean+se,fill=precip))+
  geom_bar(position="dodge",stat="identity", color="black",aes(fill=precip))+
  geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+
  scale_x_discrete(labels=c("Bad Isolates","Good Isolates", "Product Microbes","SDS Product","Sterile"),limits=FP_Isolates_order)+
  ylab("Root biomass")+scale_fill_manual(values = c("navy", "red4"),labels=c("Ambient","Drought"))+
  geom_text(aes(y=mean+se+0.02, label=n),position=position_dodge(width=0.9))+theme_bw()


#I want to see if the inculation interacts with the soil status
#need to remove the sterile treatments since it is not a part of the factorial
data_FP_biomass_trt_no_sand_live=subset(data_FP_biomass_trt_no_sand, inocul!="S")
FP_root_NS_bio_model.3= lm(sqrt(root_weight+0.0001)~precip*soil_status*inocul+as.factor(block), data= data_FP_biomass_trt_no_sand_live)
qqPlot(resid(FP_root_NS_bio_model.3))
hist(resid(FP_root_NS_bio_model.3))

Anova(FP_root_NS_bio_model.3, type=3)
#precip                    0.0522   1   5.1655 0.0241248 *  
#soil_status               0.1374   1  13.6025 0.0002925 ***
#inocul                    0.1157   3   3.8201 0.0108617 *  
#as.factor(block)          0.1019   4   2.5222 0.0423638 *  
#precip:soil_status:inocul 0.0659   3   2.1752 0.0922345 . 

emmeans(FP_root_NS_bio_model.3, ~soil_status)

data_FP_biomass_trt_no_sand_NA=data_FP_biomass_trt_no_sand[complete.cases(data_FP_biomass_trt_no_sand[,"tot_bio"]),]
data_FP_biomass_trt_no_sand_NA_g2=data_FP_biomass_trt_no_sand_NA %>% group_by(precip,inocul,soil_status)
root_data_FP_biomass_trt_no_sand_sum2=summarise_at(data_FP_biomass_trt_no_sand_NA_g2, 
                                                  "root_weight", funs(n(),mean,sd,se=sd(.)/sqrt(n())))
FP_Isolates_order=c("BI","GI","PM","SDS","S")
ggplot(root_data_FP_biomass_trt_no_sand_sum2, aes(x=inocul,y=mean, ymin = mean-se, ymax= mean+se,fill=interaction(precip,soil_status)))+
  geom_bar(position="dodge",stat="identity", color="black",aes(fill=interaction(precip,soil_status)))+
  geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+
  scale_x_discrete(labels=c("Bad Isolates","Good Isolates", "Product Microbes","SDS Product","Sterile"),limits=FP_Isolates_order)+
  ylab("Root biomass")+scale_fill_manual(values = c("navy", "red4","blue","red"),labels=c("Ambient-Live","Drought-Live","Ambient-Sterile","Drought-Sterile"))+
  geom_text(aes(y=mean+se+0.01, label=n),position=position_dodge(width=0.9))+theme_bw()



#let look at root to shoot ratio

#first we need to calculate the ratio 
data_FP_biomass_trt$root_shoot=data_FP_biomass_trt$root_weight/data_FP_biomass_trt$shoot_weight

#Shoots alone
FP_root_shoot_model= lm(log(root_shoot)~precip*soil_status*inocul+as.factor(block), data= data_FP_biomass_trt)
qqPlot(resid(FP_root_shoot_model))
hist(resid(FP_root_shoot_model))

Anova(FP_root_shoot_model, type=3,singular.ok =T)
#precip:inocul              11.384   3   3.6638 0.013139 * 
#precip:soil_status:inocul  14.346   3   4.6172 0.003726 **


emmeans(FP_shoot_bio_model, pairwise~precip|inocul|soil_status)
"$contrasts
inocul = BI, soil_status = L:
  contrast estimate     SE  df t.ratio p.value
A - D     0.02067 0.0168 236 1.228   0.2208 

inocul = GI, soil_status = L:
  contrast estimate     SE  df t.ratio p.value
A - D     0.02797 0.0163 236 1.721   0.0866 

inocul = PM, soil_status = L:
  contrast estimate     SE  df t.ratio p.value
A - D     0.02971 0.0169 236 1.761   0.0796 

inocul = S, soil_status = L:
  contrast estimate     SE  df t.ratio p.value
A - D      nonEst     NA  NA    NA       NA 

inocul = SDS, soil_status = L:
  contrast estimate     SE  df t.ratio p.value
A - D     0.00129 0.0163 236 0.079   0.9370 

inocul = BI, soil_status = S:
  contrast estimate     SE  df t.ratio p.value
A - D     0.04562 0.0165 236 2.757   0.0063 

inocul = GI, soil_status = S:
  contrast estimate     SE  df t.ratio p.value
A - D     0.02510 0.0165 236 1.517   0.1306 

inocul = PM, soil_status = S:
  contrast estimate     SE  df t.ratio p.value
A - D     0.07723 0.0169 236 4.576   <.0001 

inocul = S, soil_status = S:
  contrast estimate     SE  df t.ratio p.value
A - D     0.02693 0.0173 236 1.561   0.1199 

inocul = SDS, soil_status = S:
  contrast estimate     SE  df t.ratio p.value
A - D     0.04686 0.0165 236 2.832   0.0050 "

FP_root_shoot_model.2= lm(log(root_shoot)~precip*inocul+as.factor(block), data= data_FP_biomass_trt)
qqPlot(resid(FP_root_shoot_model.2))
hist(resid(FP_root_shoot_model.2))

Anova(FP_root_shoot_model.2, type=3)
#precip            35.762   1  31.3463 6.1e-08 ***
#precip:inocul     10.541   4   2.3097  0.0587 .


emmeans(FP_root_shoot_model.2, pairwise~precip|inocul)
#$contrasts
#inocul = BI:
#  contrast estimate    SE  df t.ratio p.value
#A - D      -1.274 0.297 231 -4.290  <.0001 

#inocul = GI:
#  contrast estimate    SE  df t.ratio p.value
#A - D      -1.113 0.288 231 -3.862  0.0001 

#inocul = PM:
#  contrast estimate    SE  df t.ratio p.value
#A - D      -0.607 0.289 231 -2.097  0.0371 

#inocul = S:
#  contrast estimate    SE  df t.ratio p.value
#A - D      -0.838 0.414 231 -2.025  0.0440 
#inocul = SDS:
#contrast estimate    SE  df t.ratio p.value
#A - D      -0.157 0.286 231 -0.550  0.5826 


data_FP_biomass_trt_NA=data_FP_biomass_trt[complete.cases(data_FP_biomass_trt[,"root_shoot"]),]
data_FP_biomass_trt_NA_g=data_FP_biomass_trt_NA %>% group_by(precip,inocul)
root_shoot_data_FP_biomass_trt_sum=summarise_at(data_FP_biomass_trt_NA_g, 
                                               "root_shoot", funs(n(),mean,sd,se=sd(.)/sqrt(n())))
FP_Isolates_order=c("BI","GI","PM","SDS","S")
ggplot(root_shoot_data_FP_biomass_trt_sum, aes(x=inocul,y=mean, ymin = mean-se, ymax= mean+se,fill=precip))+
  geom_bar(position="dodge",stat="identity", color="black",aes(fill=precip))+
  geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+
  scale_x_discrete(labels=c("Bad Isolates","Good Isolates", "Product Microbes","SDS Product","Sterile"),limits=FP_Isolates_order)+
  ylab("Root:Shoot ratio")+scale_fill_manual(values = c("navy", "red4"),labels=c("Ambient","Drought"))+
  geom_text(aes(y=mean+se+0.5, label=n),position=position_dodge(width=0.9))+theme_bw()


#I want to see if the inculation interacts with the soil status
#need to remove the sterile treatments since it is not a part of the factorial
data_FP_biomass_trt_live=subset(data_FP_biomass_trt, inocul!="S")
FP_root_shoot_bio_model.3= lm(log(root_shoot)~precip*soil_status*inocul+as.factor(block), data= data_FP_biomass_trt_live)
qqPlot(resid(FP_root_shoot_bio_model.3))
hist(resid(FP_root_shoot_bio_model.3))

Anova(FP_root_shoot_bio_model.3, type=3)
#precip                     33.624   1  34.4030 1.854e-08 ***
#soil_status                15.588   1  15.9487 9.169e-05 *** 
#precip:inocul              11.394   3   3.8858  0.009950 ** 
#precip:soil_status:inocul  14.233   3   4.8542  0.002789 ** 
emmeans(FP_root_shoot_bio_model.3, pairwise~precip|inocul|soil_status)
"$contrasts
inocul = BI, soil_status = L:
 contrast estimate    SE  df t.ratio p.value
 A - D      -1.004 0.405 198 -2.478  0.0141 

inocul = GI, soil_status = L:
 contrast estimate    SE  df t.ratio p.value
 A - D      -2.000 0.388 198 -5.154  <.0001 

inocul = PM, soil_status = L:
 contrast estimate    SE  df t.ratio p.value
 A - D      -0.514 0.383 198 -1.341  0.1813 

inocul = SDS, soil_status = L:
 contrast estimate    SE  df t.ratio p.value
 A - D       0.372 0.368 198  1.012  0.3129 

inocul = BI, soil_status = S:
 contrast estimate    SE  df t.ratio p.value
 A - D      -1.536 0.375 198 -4.099  0.0001 

inocul = GI, soil_status = S:
 contrast estimate    SE  df t.ratio p.value
 A - D      -0.306 0.368 198 -0.831  0.4067 

inocul = PM, soil_status = S:
 contrast estimate    SE  df t.ratio p.value
 A - D      -0.706 0.375 198 -1.883  0.0611 

inocul = SDS, soil_status = S:
 contrast estimate    SE  df t.ratio p.value
 A - D      -0.616 0.383 198 -1.608  0.1093 "

data_FP_biomass_trt_NA_g2=data_FP_biomass_trt_NA %>% group_by(precip,inocul,soil_status)
root_shoot_data_FP_biomass_trt_sum2=summarise_at(data_FP_biomass_trt_NA_g2, 
                                                "root_shoot", funs(n(),mean,sd,se=sd(.)/sqrt(n())))
FP_Isolates_order=c("BI","GI","PM","SDS","S")
ggplot(root_shoot_data_FP_biomass_trt_sum2, aes(x=inocul,y=mean, ymin = mean-se, ymax= mean+se,fill=interaction(precip,soil_status)))+
  geom_bar(position="dodge",stat="identity", color="black",aes(fill=interaction(precip,soil_status)))+
  geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+
  scale_x_discrete(labels=c("Bad Isolates","Good Isolates", "Product Microbes","SDS Product","Sterile"),limits=FP_Isolates_order)+
  ylab("Root:Shoot ratio")+scale_fill_manual(values = c("navy", "red4","blue","red"),labels=c("Ambient-Live","Drought-Live","Ambient-Sterile","Drought-Sterile"))+
  geom_text(aes(y=mean+se+0.7, label=n),position=position_dodge(width=0.9))+theme_bw()
