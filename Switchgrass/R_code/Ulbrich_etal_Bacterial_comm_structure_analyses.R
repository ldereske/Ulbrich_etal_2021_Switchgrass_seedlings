######Here is the Code for analyzing the bacterial community composition of MERDS Switchgrass Experiment#######
#Lukas (belldere@msu.edu) processed the sequences and subset the larger run to only have samples from MERDS (MERDs 56 of 317 samples in the run)

library(phyloseq)
library(RVAideMemoire) # for tests after PERMANOVA
library(ggplot2)
library(ggrepel)
library(dplyr)
library(gridExtra)
library(car)
library(emmeans)
library(vegan)
library(ggpubr)
library(otuSummary)
source(system.file("extdata/lm_phyloseq.R", package = "microbiome"))
options(contrasts=c("contr.sum", "contr.poly"))

AICc.PERMANOVA <- function(adonis.model) {
  
  # check to see if object is an adonis model...
  
  if (!(adonis.model$aov.tab[1,1] >= 1))
    stop("object not output of adonis {vegan} ")
  
  # Ok, now extract appropriate terms from the adonis model
  # Calculating AICc using residual sum of squares (RSS) since I don't think that adonis returns something I can use as a liklihood function...
  
  RSS <- adonis.model$aov.tab[rownames(adonis.model$aov.tab) == "Residuals", "SumsOfSqs"]
  MSE <- adonis.model$aov.tab[rownames(adonis.model$aov.tab) == "Residuals", "MeanSqs"]
  
  k <- ncol(adonis.model$model.matrix)# + 1 # add one for error variance
  
  nn <- nrow(adonis.model$model.matrix)
  
  # AIC : 2*k + n*ln(RSS)
  # AICc: AIC + [2k(k+1)]/(n-k-1)
  
  # based on https://en.wikipedia.org/wiki/Akaike_information_criterion;
  # https://www.researchgate.net/post/What_is_the_AIC_formula;
  # http://avesbiodiv.mncn.csic.es/estadistica/ejemploaic.pdf
  
  # AIC.g is generalized version of AIC = 2k + n [Ln( 2(pi) RSS/n ) + 1]
  # AIC.pi = k + n [Ln( 2(pi) RSS/(n-k) ) +1],
  
  AIC <- 2*k + nn*log(RSS)
  AIC.g <- 2*k + nn * (1 + log( 2 * pi * RSS / nn))
  AIC.MSE <- 2*k + nn * log(MSE)
  AIC.pi <- k + nn*(1 + log( 2*pi*RSS/(nn-k) )   )
  AICc <- AIC + (2*k*(k + 1))/(nn - k - 1)
  AICc.MSE <- AIC.MSE + (2*k*(k + 1))/(nn - k - 1)
  AICc.pi <- AIC.pi + (2*k*(k + 1))/(nn - k - 1)
  
  output <- list("AIC" = AIC, "AIC.g" = AIC.g, "AICc" = AICc,
                 "AIC.MSE" = AIC.MSE, "AICc.MSE" = AICc.MSE,
                 "AIC.pi" = AIC.pi, "AICc.pi" = AICc.pi, "k" = k)
  
  return(output)   
  
}




#####Begin OTU based community####
#Tree 
run_20190617_16S.tree = read_tree("D:/MMPRNT_16S_016-018/USEARCH/tree_OTU_20190617_16S-V4_PE250_NWK.NWK")
#Let's load in the SILVA 123 version of the phyloseq object
load("D:/MERDS_2018/merds/Switchgrass/R_data/phyl_obj_SILVA_MERDS.RData")
#root the tree with a random node
head(taxa_names(phyl_SILVA_MERDS))
phyl_SILVA_MERDS=phyloseq(tax_table(phyl_SILVA_MERDS),otu_table(phyl_SILVA_MERDS),sample_data(phyl_SILVA_MERDS), phy_tree(run_20190617_16S.tree))
phy_tree(phyl_SILVA_MERDS)<-ape::root(phy_tree(phyl_SILVA_MERDS), "OTU3039", resolve.root=TRUE)
ntaxa(phyl_SILVA_MERDS)
#11892
sum(taxa_sums(phyl_SILVA_MERDS))
#980588
mean(sample_sums(phyl_SILVA_MERDS))
#17510.5
min(sample_sums(phyl_SILVA_MERDS))
#242
max(sample_sums(phyl_SILVA_MERDS))
#30633
sort(sample_sums(phyl_SILVA_MERDS))
#  MERDSSG42  MERDSSG85
#       242      10703


get_taxa_unique(SILVA_MERDS_data, taxonomic.rank = "Domain")
SILVA_MERDS_data<-subset_taxa(SILVA_MERDS_data,Domain!="")
ntaxa(SILVA_MERDS_data)
#11887
sum(taxa_sums(SILVA_MERDS_data))
#980577
get_taxa_unique(SILVA_MERDS_data, taxonomic.rank = "Domain")
SILVA_MERDS_data<-subset_taxa(SILVA_MERDS_data,Class!="c:Chloroplast")
ntaxa(SILVA_MERDS_data)
#11797
sum(otu_table(SILVA_MERDS_data))
#978323

SILVA_MERDS_data<-subset_taxa(SILVA_MERDS_data,Family!="f:Mitochondria")
ntaxa(SILVA_MERDS_data)
#11774
sum(otu_table(SILVA_MERDS_data))
#977652
sort(sample_sums(SILVA_MERDS_data))
# MERDSSG42  MERDSSG85 MERDSSG103
#       242      10694      11490
df=data.frame(otu_table(SILVA_MERDS_data))
df$OTU=row.names(df)

subset(df, OTU=="OTU9208")
#remove sample with 239 read
SILVA_MERDS_trunc=prune_samples(sample_sums(SILVA_MERDS_data) > 10000, SILVA_MERDS_data)
nrow(sample_data(SILVA_MERDS_trunc))
#55

SILVA_MERDS_trunc3=prune_taxa(taxa_sums(SILVA_MERDS_trunc) > 2, SILVA_MERDS_trunc)
ntaxa(SILVA_MERDS_trunc3)
#8864
sum(otu_table(SILVA_MERDS_trunc3))
#972894

######END Processing#### 

######Rarefied community Analyses####

#Let's rarefy the data
#SILVA_MERDS_rar=rarefy_even_depth(SILVA_MERDS_trunc3, sample.size= 10000, rngseed = T)
#498OTUs were removed because they are no longer 
#present in any sample after random subsampling

#save(SILVA_MERDS_rar, file = "D:/MERDS_2018/merds/Switchgrass/R_data/SILVA_MERDS_rar_phylo_obj.RData")
load("D:/MERDS_2018/merds/Switchgrass/R_data/SILVA_MERDS_rar_phylo_obj.RData")
#SILVA_MERDS_rar=phyloseq(otu_table(SILVA_MERDS_rar),tax_table(SILVA_MERDS_rar),sample_data(merds_map_trt_bio_nit),phy_tree(run_20190617_16S.tree))
#phy_tree(SILVA_MERDS_rar)<-ape::root(phy_tree(SILVA_MERDS_rar), "OTU3039", resolve.root=TRUE)

#####Change from intial####

SILVA_MERDS_rar_map=sample_data(SILVA_MERDS_rar)
head(SILVA_MERDS_rar_map)
SILVA_MERDS_rar_WU_dis=distance(SILVA_MERDS_rar,method = "wunifrac")
SILVA_MERDS_rar_WU_dis_M <- matrixConvert(SILVA_MERDS_rar_WU_dis, 
                                          colname = c("sample1", "sample2", "w_unifrac"))

head(SILVA_MERDS_rar_WU_dis_M)
nrow(SILVA_MERDS_rar_WU_dis_M)
#1485
SILVA_MERDS_rar_WU_dis_M$s1_s2=with(SILVA_MERDS_rar_WU_dis_M, interaction(sample1,sample2))

SILVA_MERDS_rar_WU_dis_trt0=merge(SILVA_MERDS_rar_WU_dis_M,SILVA_MERDS_rar_map[,c("soil_status","root_association","block",
                                                                                  "precip","life_stage","soil_root")], 
                                  by.x = "sample1",by.y = "row.names")
nrow(SILVA_MERDS_rar_WU_dis_trt0)
#1485
colnames(SILVA_MERDS_rar_WU_dis_trt0)[5:10]=c("s1_soil_status","s1_root_association","s1_block",
                                              "s1_precip","s1_life_stage","s1_soil_root")



SILVA_MERDS_rar_WU_dis_trt=merge(SILVA_MERDS_rar_WU_dis_trt0,SILVA_MERDS_rar_map[,c("soil_status","root_association","block",
                                                                                    "precip","life_stage","soil_root")], 
                                 by.x = "sample2",by.y = "row.names")


nrow(SILVA_MERDS_rar_WU_dis_trt)
#1485
colnames(SILVA_MERDS_rar_WU_dis_trt)[11:16]=c("s2_soil_status","s2_root_association","s2_block",
                                              "s2_precip","s2_life_stage","s2_soil_root")

SILVA_MERDS_rar_WU_dis_trt$s1_s2_life_stage=with(SILVA_MERDS_rar_WU_dis_trt, interaction(s1_life_stage,s2_life_stage))
unique(SILVA_MERDS_rar_WU_dis_trt$s1_s2_life_stage)

SILVA_MERDS_rar_WU_dis_trt$s1_s2_root_association=with(SILVA_MERDS_rar_WU_dis_trt, interaction(s1_root_association,s2_root_association))
unique(SILVA_MERDS_rar_WU_dis_trt$s1_s2_root_association)

SILVA_MERDS_rar_WU_dis_trt$s1_s2_block=with(SILVA_MERDS_rar_WU_dis_trt, interaction(s1_block,s2_block))
unique(SILVA_MERDS_rar_WU_dis_trt$s1_s2_block)

SILVA_MERDS_rar_WU_dis_trt$time=ifelse(SILVA_MERDS_rar_WU_dis_trt$s1_s2_life_stage=="Start.Start","Start","End")

unique(SILVA_MERDS_rar_WU_dis_trt$s2_life_stage)
SILVA_MERDS_rar_WU_dis_start=subset(SILVA_MERDS_rar_WU_dis_trt,s1_s2_life_stage=="S.Start"|
                                      s1_s2_life_stage=="G.Start"|s1_s2_life_stage=="Start.Start")
nrow(SILVA_MERDS_rar_WU_dis_start)
#404
unique(SILVA_MERDS_rar_WU_dis_start$s1_life_stage)
unique(SILVA_MERDS_rar_WU_dis_start$s2_life_stage)

SILVA_MERDS_rar_WU_dis_start %>% group_by(time,s1_soil_root,s1_life_stage,s1_root_association,
                                          s2_root_association) %>% summarise_at(vars(w_unifrac),c(~mean(.),~n()))

SILVA_MERDS_rar_WU_dis_start_wi=subset(SILVA_MERDS_rar_WU_dis_start,s1_s2_root_association=="B.B"|
                                         s1_s2_root_association=="R.R")

nrow(SILVA_MERDS_rar_WU_dis_start_wi)
#200
unique(SILVA_MERDS_rar_WU_dis_start$s1_block)

SILVA_MERDS_rar_WU_dis_start_wi_bl=subset(SILVA_MERDS_rar_WU_dis_start_wi,time=="Start"|s1_s2_block=="2.2"|
                                            s1_s2_block=="3.3"|s1_s2_block=="4.4"|s1_s2_block=="5.5")
nrow(SILVA_MERDS_rar_WU_dis_start_wi_bl)
#59


head(SILVA_MERDS_rar_WU_dis_start_wi_bl)
unique(SILVA_MERDS_rar_WU_dis_start_wi_bl$s1_s2_life_stage)




#####Transplant Experimental Distance from start community seedling####

SILVA_MERDS_rar_WU_dis_start_wi_bl_trans_exp=subset(SILVA_MERDS_rar_WU_dis_start_wi_bl,s1_s2_life_stage!="S.Start")
nrow(SILVA_MERDS_rar_WU_dis_start_wi_bl_trans_exp)
#36
unique(SILVA_MERDS_rar_WU_dis_start_wi_bl_trans_exp$s1_life_stage)

trans_WU_intercept_exp=SILVA_MERDS_rar_WU_dis_start_wi_bl_trans_exp %>% group_by(time,s1_soil_root,s1_precip,s1_soil_root) %>% 
  summarise_at(vars(w_unifrac),c(~mean(.),~n(),se=~sd(.)/sqrt(n())))
unique(SILVA_MERDS_rar_WU_dis_start_wi_bl_trans_exp$s1_soil_root)

ggplot(SILVA_MERDS_rar_WU_dis_start_wi_bl_trans_exp)+
  geom_point(aes(x=factor(time, levels = c("Start","End")),y=w_unifrac,fill=s1_precip, shape=s1_soil_root), size=3)+
  geom_segment(aes(x=1,xend=2,y=subset(trans_WU_intercept_exp,time=="Start"&s1_soil_root=="L.B")$mean,
                   yend=subset(trans_WU_intercept_exp,s1_soil_root=="L.B"&s1_precip=="A")$mean), linetype="solid", color="turquoise",size=1.5)+
  geom_segment(aes(x=1,xend=2,y=subset(trans_WU_intercept_exp,time=="Start"&s1_soil_root=="L.B")$mean,
                   yend=subset(trans_WU_intercept_exp,s1_soil_root=="L.B"&s1_precip=="D")$mean), linetype="solid", color="firebrick1",size=1.5)+
  geom_segment(aes(x=1,xend=2,y=subset(trans_WU_intercept_exp,time=="Start"&s1_soil_root=="L.R")$mean,
                   yend=subset(trans_WU_intercept_exp,s1_soil_root=="L.R"&s1_precip=="A")$mean), linetype="dashed", color="turquoise",size=1.5)+
  geom_segment(aes(x=1,xend=2,y=subset(trans_WU_intercept_exp,time=="Start"&s1_soil_root=="L.R")$mean,
                   yend=subset(trans_WU_intercept_exp,s1_soil_root=="L.R"&s1_precip=="D")$mean), linetype="dashed", color="firebrick1",size=1.5)+
  geom_segment(aes(x=1,xend=2,y=subset(trans_WU_intercept_exp,time=="Start"&s1_soil_root=="L.R")$mean,
                   yend=subset(trans_WU_intercept_exp,s1_soil_root=="S.B"&s1_precip=="A")$mean), linetype="dotted", color="turquoise",size=1.5)+
  geom_segment(aes(x=1,xend=2,y=subset(trans_WU_intercept_exp,time=="Start"&s1_soil_root=="L.R")$mean,
                   yend=subset(trans_WU_intercept_exp,s1_soil_root=="S.B"&s1_precip=="D")$mean), linetype="dotted", color="firebrick1",size=1.5)+
  geom_point(data = trans_WU_intercept_exp,aes(x=factor(time, levels = c("Start","End")),y=mean,fill=s1_precip, shape=s1_soil_root), size=7)+
  geom_errorbar(data = trans_WU_intercept_exp,aes(x=factor(time, levels = c("Start","End")),y=mean,ymin=mean-se,ymax=mean+se,color=s1_precip),
                width=.1,size=1)+
  scale_y_continuous(name = "Weighted Unifrac distance")+scale_color_manual(values = c("turquoise","firebrick1","black"))+
  scale_shape_manual(values= c(21,23,24),name=NULL)+scale_fill_manual(values = c("turquoise","firebrick1","black"))+
  theme_bw()+theme(axis.text = element_text(size = 18),axis.title.x = element_blank(), axis.title.y = element_text(size = 22),
                   legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#comm_turnover_line_seedling
#900x750

#####Presence Distance from start community seedling####

SILVA_MERDS_rar_WU_dis_wi_bl_trans_pres=subset(SILVA_MERDS_rar_WU_dis_start_wi_bl_trans_exp,time!="Start"&s1_soil_root!="L.R")
nrow(SILVA_MERDS_rar_WU_dis_wi_bl_trans_pres)
#16

WU_pair_trans_mod_pres= lm(w_unifrac~s1_soil_root*s1_precip, data= SILVA_MERDS_rar_WU_dis_wi_bl_trans_pres)
qqPlot(resid(WU_pair_trans_mod_pres))
hist(resid(WU_pair_trans_mod_pres))
shapiro.test(resid(WU_pair_trans_mod_pres))
#0.2615

Anova(WU_pair_trans_mod_pres, type=3)
#s1_soil_root           0.19065  1  159.2094 2.759e-08 ***

emmeans(WU_pair_trans_mod_pres, pairwise~s1_precip*s1_soil_root,adjust="none")




#####Origin Distance from start community seedling####
SILVA_MERDS_rar_WU_dis_start_wi_bl_trans=subset(SILVA_MERDS_rar_WU_dis_start_wi_bl,s1_soil_root!="S.B"&s1_s2_life_stage!="S.Start")
nrow(SILVA_MERDS_rar_WU_dis_start_wi_bl_trans)
#28
trans_WU_intercept=SILVA_MERDS_rar_WU_dis_start_wi_bl_trans %>% group_by(time,s1_soil_root,s1_precip,s1_soil_root) %>% 
  summarise_at(vars(w_unifrac),c(~mean(.),~n(),se=~sd(.)/sqrt(n())))

ggplot(SILVA_MERDS_rar_WU_dis_start_wi_bl_trans)+
  geom_point(aes(x=factor(time, levels = c("Start","End")),y=w_unifrac,fill=s1_precip, shape=s1_soil_root), size=3)+
  geom_segment(aes(x=1,xend=2,y=subset(trans_WU_intercept,time=="Start"&s1_soil_root=="L.B")$mean,
                   yend=subset(trans_WU_intercept,s1_soil_root=="L.B"&s1_precip=="A")$mean), linetype="solid", color="turquoise",size=1.5)+
  geom_segment(aes(x=1,xend=2,y=subset(trans_WU_intercept,time=="Start"&s1_soil_root=="L.B")$mean,
                   yend=subset(trans_WU_intercept,s1_soil_root=="L.B"&s1_precip=="D")$mean), linetype="solid", color="firebrick1",size=1.5)+
  geom_segment(aes(x=1,xend=2,y=subset(trans_WU_intercept,time=="Start"&s1_soil_root=="L.R")$mean,
                   yend=subset(trans_WU_intercept,s1_soil_root=="L.R"&s1_precip=="A")$mean), linetype="dashed", color="turquoise",size=1.5)+
  geom_segment(aes(x=1,xend=2,y=subset(trans_WU_intercept,time=="Start"&s1_soil_root=="L.R")$mean,
                   yend=subset(trans_WU_intercept,s1_soil_root=="L.R"&s1_precip=="D")$mean), linetype="dashed", color="firebrick1",size=1.5)+
  geom_point(data = trans_WU_intercept,aes(x=factor(time, levels = c("Start","End")),y=mean,fill=s1_precip, shape=s1_soil_root), size=7)+
  geom_errorbar(data = trans_WU_intercept,aes(x=factor(time, levels = c("Start","End")),y=mean,ymin=mean-se,ymax=mean+se,color=s1_precip),
                width=.1,size=1)+
  scale_y_continuous(name = "Weighted Unifrac distance",limits = c(0.075,0.35),breaks = c(0.10,0.15,0.20,0.25,0.30,0.35))+scale_color_manual(values = c("turquoise","firebrick1","black"))+
  scale_shape_manual(values= c(21,23,24),name=NULL)+scale_fill_manual(values = c("turquoise","firebrick1","black"))+
  theme_bw()+theme(axis.text = element_text(size = 18),axis.title.x = element_blank(), axis.title.y = element_text(size = 22),
                   legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#comm_turnover_line_seedling_soil_origin
#900x750
SILVA_MERDS_rar_WU_dis_wi_bl_trans=subset(SILVA_MERDS_rar_WU_dis_start_wi_bl_trans,time!="Start")

WU_pair_trans_mod= lm(w_unifrac~s1_soil_root*s1_precip, data= SILVA_MERDS_rar_WU_dis_wi_bl_trans)
qqPlot(resid(WU_pair_trans_mod))
hist(resid(WU_pair_trans_mod))
shapiro.test(resid(WU_pair_trans_mod))
#0.6766

Anova(WU_pair_trans_mod, type=3)
#s1_soil_root:s1_precip 0.01454  1   5.1174   0.04304 * 

emmeans(WU_pair_trans_mod, pairwise~s1_precip*s1_soil_root,adjust="none")




#####Seed Experimental Distance from start community germination####
SILVA_MERDS_rar_WU_dis_start_wi_bl_seed_exp=subset(SILVA_MERDS_rar_WU_dis_start_wi_bl,s1_soil_root!="G.B"&s1_s2_life_stage!="G.Start")
nrow(SILVA_MERDS_rar_WU_dis_start_wi_bl_seed_exp)
#35

unique(SILVA_MERDS_rar_WU_dis_start_wi_bl_seed_exp$s1_life_stage)

seed_WU_intercept_exp=SILVA_MERDS_rar_WU_dis_start_wi_bl_seed_exp %>% group_by(time,s1_soil_root,s1_precip,s1_soil_root) %>%
  summarise_at(vars(w_unifrac),c(~mean(.),~n(),se=~sd(.)/sqrt(n())))
unique(SILVA_MERDS_rar_WU_dis_start_wi_bl_seed_exp$s1_soil_root)

ggplot(SILVA_MERDS_rar_WU_dis_start_wi_bl_seed_exp)+
  geom_point(aes(x=factor(time, levels = c("Start","End")),y=w_unifrac,fill=s1_precip, shape=s1_soil_root), size=3)+
  geom_segment(aes(x=1,xend=2,y=subset(seed_WU_intercept_exp,time=="Start"&s1_soil_root=="L.B")$mean,
                   yend=subset(seed_WU_intercept_exp,s1_soil_root=="L.B"&s1_precip=="A")$mean), linetype="solid", color="turquoise",size=1.5)+
  geom_segment(aes(x=1,xend=2,y=subset(seed_WU_intercept_exp,time=="Start"&s1_soil_root=="L.B")$mean,
                   yend=subset(seed_WU_intercept_exp,s1_soil_root=="L.B"&s1_precip=="D")$mean), linetype="solid", color="firebrick1",size=1.5)+
  geom_segment(aes(x=1,xend=2,y=subset(seed_WU_intercept_exp,time=="Start"&s1_soil_root=="L.R")$mean,
                   yend=subset(seed_WU_intercept_exp,s1_soil_root=="L.R"&s1_precip=="A")$mean), linetype="dashed", color="turquoise",size=1.5)+
  geom_segment(aes(x=1,xend=2,y=subset(seed_WU_intercept_exp,time=="Start"&s1_soil_root=="L.R")$mean,
                   yend=subset(seed_WU_intercept_exp,s1_soil_root=="L.R"&s1_precip=="D")$mean), linetype="dashed", color="firebrick1",size=1.5)+
  geom_segment(aes(x=1,xend=2,y=subset(seed_WU_intercept_exp,time=="Start"&s1_soil_root=="L.R")$mean,
                   yend=subset(seed_WU_intercept_exp,s1_soil_root=="S.B"&s1_precip=="A")$mean), linetype="dotted", color="turquoise",size=1.5)+
  geom_segment(aes(x=1,xend=2,y=subset(seed_WU_intercept_exp,time=="Start"&s1_soil_root=="L.R")$mean,
                   yend=subset(seed_WU_intercept_exp,s1_soil_root=="S.B"&s1_precip=="D")$mean), linetype="dotted", color="firebrick1",size=1.5)+
  geom_point(data = seed_WU_intercept_exp,aes(x=factor(time, levels = c("Start","End")),y=mean,fill=s1_precip, shape=s1_soil_root), size=7)+
  geom_errorbar(data = seed_WU_intercept_exp,aes(x=factor(time, levels = c("Start","End")),y=mean,ymin=mean-se,ymax=mean+se,color=s1_precip),
                width=.1,size=1)+
  scale_y_continuous(name = "Weighted Unifrac distance")+scale_color_manual(values = c("turquoise","firebrick1","black"))+
  scale_shape_manual(values= c(21,23,24),name=NULL)+scale_fill_manual(values = c("turquoise","firebrick1","black"))+
  theme_bw()+theme(axis.text = element_text(size = 18),axis.title.x = element_blank(), axis.title.y = element_text(size = 22),
                   legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#comm_turnover_line_germination

#####Presence Distance from start community germination####

SILVA_MERDS_rar_WU_dis_wi_bl_seed_pres=subset(SILVA_MERDS_rar_WU_dis_start_wi_bl_seed_exp,time!="Start"&s1_soil_root!="L.R")
nrow(SILVA_MERDS_rar_WU_dis_wi_bl_seed_pres)
#15

WU_pair_seed_mod_pres= lm(w_unifrac~s1_soil_root*s1_precip, data= SILVA_MERDS_rar_WU_dis_wi_bl_seed_pres)
qqPlot(resid(WU_pair_seed_mod_pres))
hist(resid(WU_pair_seed_mod_pres))
shapiro.test(resid(WU_pair_seed_mod_pres))
#0.8458

Anova(WU_pair_seed_mod_pres, type=3)
#s1_soil_root           0.06906  1  36.5262 8.384e-05 ***

emmeans(WU_pair_seed_mod_pres, pairwise~s1_precip*s1_soil_root,adjust="none")


#####origin Distance from start community germination####

SILVA_MERDS_rar_WU_dis_start_wi_bl_seed=subset(SILVA_MERDS_rar_WU_dis_start_wi_bl,s1_soil_root!="G.B"&s1_s2_life_stage!="G.Start"&s1_soil_root!="S.B")
nrow(SILVA_MERDS_rar_WU_dis_start_wi_bl_seed)
#28
seed_WU_intercept=SILVA_MERDS_rar_WU_dis_start_wi_bl_seed %>% group_by(time,s1_soil_root,s1_precip,s1_soil_root) %>% 
  summarise_at(vars(w_unifrac),c(~mean(.),~n(),se=~sd(.)/sqrt(n())))

ggplot(SILVA_MERDS_rar_WU_dis_start_wi_bl_seed)+
  geom_point(aes(x=factor(time, levels = c("Start","End")),y=w_unifrac,fill=s1_precip, shape=s1_soil_root), size=3)+
  geom_segment(aes(x=1,xend=2,y=subset(seed_WU_intercept,time=="Start"&s1_soil_root=="L.B")$mean,
                   yend=subset(seed_WU_intercept,s1_soil_root=="L.B"&s1_precip=="A")$mean), linetype="solid", color="turquoise",size=1.5)+
  geom_segment(aes(x=1,xend=2,y=subset(seed_WU_intercept,time=="Start"&s1_soil_root=="L.B")$mean,
                   yend=subset(seed_WU_intercept,s1_soil_root=="L.B"&s1_precip=="D")$mean), linetype="solid", color="firebrick1",size=1.5)+
  geom_segment(aes(x=1,xend=2,y=subset(seed_WU_intercept,time=="Start"&s1_soil_root=="L.R")$mean,
                   yend=subset(seed_WU_intercept,s1_soil_root=="L.R"&s1_precip=="A")$mean), linetype="dashed", color="turquoise",size=1.5)+
  geom_segment(aes(x=1,xend=2,y=subset(seed_WU_intercept,time=="Start"&s1_soil_root=="L.R")$mean,
                   yend=subset(seed_WU_intercept,s1_soil_root=="L.R"&s1_precip=="D")$mean), linetype="dashed", color="firebrick1",size=1.5)+
  geom_point(data = seed_WU_intercept,aes(x=factor(time, levels = c("Start","End")),y=mean,fill=s1_precip, shape=s1_soil_root), size=7)+
  geom_errorbar(data = seed_WU_intercept,aes(x=factor(time, levels = c("Start","End")),y=mean,ymin=mean-se,ymax=mean+se,color=s1_precip),
                width=.1,size=1)+
  scale_y_continuous(name = "Weighted Unifrac distance",limits = c(0.075,0.35),breaks = c(0.10,0.15,0.20,0.25,0.30,0.35))+scale_color_manual(values = c("turquoise","firebrick1","black"))+
  scale_shape_manual(values= c(21,23,24),name=NULL)+scale_fill_manual(values = c("turquoise","firebrick1","black"))+
  theme_bw()+theme(axis.text = element_text(size = 18),axis.title.x = element_blank(), axis.title.y = element_text(size = 22), 
                   legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#comm_turnover_line_germination_soil_origin

SILVA_MERDS_rar_WU_dis_wi_bl_seed=subset(SILVA_MERDS_rar_WU_dis_start_wi_bl_seed,time!="Start")

WU_pair_seed_mod= lm(w_unifrac~s1_soil_root*s1_precip, data= SILVA_MERDS_rar_WU_dis_wi_bl_seed)
qqPlot(resid(WU_pair_seed_mod))
hist(resid(WU_pair_seed_mod))
shapiro.test(resid(WU_pair_seed_mod))
#0.8664

Anova(WU_pair_seed_mod, type=3)
#s1_soil_root           0.03289  1  17.5538  0.001255 ** 

emmeans(WU_pair_seed_mod, pairwise~s1_precip*s1_soil_root,adjust="none")



#####No Start####
#live
SILVA_MERDS_rar_live=subset_samples(SILVA_MERDS_rar, soil_status!="S"|life_stage=="Start")
nsamples(SILVA_MERDS_rar_live)
#40
SILVA_MERDS_rar_live_map=sample_data(SILVA_MERDS_rar_live)
#####Weighted Unifrac####
SILVA_MERDS_rar_live_NS=subset_samples(SILVA_MERDS_rar_live, life_stage!="Start")
nsamples(SILVA_MERDS_rar_live_NS)
#32
SILVA_MERDS_rar_live_NS_WU_ord=ordinate(SILVA_MERDS_rar_live_NS, method = "NMDS", distance = "wunifrac")
#*** Solution reached
#0.08729953  
plot_ordination(SILVA_MERDS_rar_live_NS,SILVA_MERDS_rar_live_NS_WU_ord, color="root_association",shape="life_stage")+geom_point(size=3)+
  theme_bw()



SILVA_MERDS_rar_live_NS_map=sample_data(SILVA_MERDS_rar_live_NS)
#SILVA_MERDS_rar_map$soil_root_stage=with(SILVA_MERDS_rar_map, interaction(soil_status,root_association,life_stage))
SILVA_MERDS_rar_live_NS_WU_dis=distance(SILVA_MERDS_rar_live_NS,method = "wunifrac")

(wu_full_P_mod=adonis(SILVA_MERDS_rar_live_NS_WU_dis~SILVA_MERDS_rar_live_NS_map$root_association*
                        SILVA_MERDS_rar_live_NS_map$life_stage*SILVA_MERDS_rar_live_NS_map$precip+
                        as.factor(SILVA_MERDS_rar_live_NS_map$block), permutations = 9999))
#SILVA_MERDS_rar_live_NS_map$root_association                                                                            1   0.03482 0.034816  3.3083 0.07598 0.0173 * 
#SILVA_MERDS_rar_live_NS_map$precip                                                                                      1   0.05128 0.051282  4.8729 0.11191 0.0026 **
#SILVA_MERDS_rar_live_NS_map$root_association:SILVA_MERDS_rar_live_NS_map$life_stage                                     1   0.04608 0.046077  4.3783 0.10055 0.0048 **
#SILVA_MERDS_rar_live_NS_map$root_association:SILVA_MERDS_rar_live_NS_map$life_stage:SILVA_MERDS_rar_live_NS_map$precip  1   0.02174 0.021736  2.0654 0.04743 0.0759 . 

AICc.PERMANOVA(wu_full_P_mod)
#$AIC
#[1] -26.3061


(wu_full_P_mod_no_block=adonis(SILVA_MERDS_rar_live_NS_WU_dis~SILVA_MERDS_rar_live_NS_map$root_association*
                                 SILVA_MERDS_rar_live_NS_map$life_stage*SILVA_MERDS_rar_live_NS_map$precip, permutations = 9999))
#SILVA_MERDS_rar_live_NS_map$root_association                                                                            1   0.03482 0.034816  3.1675 0.07598 0.0190 * 
#SILVA_MERDS_rar_live_NS_map$precip                                                                                      1   0.05128 0.051282  4.6656 0.11191 0.0036 **
#SILVA_MERDS_rar_live_NS_map$root_association:SILVA_MERDS_rar_live_NS_map$life_stage                                     1   0.04978 0.049783  4.5291 0.10864 0.0044 **
#SILVA_MERDS_rar_live_NS_map$root_association:SILVA_MERDS_rar_live_NS_map$life_stage:SILVA_MERDS_rar_live_NS_map$precip  1   0.02363 0.023628  2.1496 0.05156 0.0658 . 

AICc.PERMANOVA(wu_full_P_mod_no_block)
#$AIC
#-26.6422

#pairwise.perm.manova(SILVA_MERDS_rar_WU_dis, SILVA_MERDS_rar_map$soil_root_stage, nperm=2000)

#Transplant


SILVA_MERDS_rar_trans=subset_samples(SILVA_MERDS_rar, life_stage=="G"|life_stage=="Start")
SILVA_MERDS_rar_trans_map=sample_data(SILVA_MERDS_rar_trans)
summary(sample_data(SILVA_MERDS_rar_trans_map))
SILVA_MERDS_rar_trans_WU_ord=ordinate(SILVA_MERDS_rar_trans, method = "NMDS",distance = "wunifrac")
#*** Solution reached
#0.04807855 
SILVA_MERDS_rar_trans_ord_WU_points=merge(SILVA_MERDS_rar_trans_WU_ord$points,SILVA_MERDS_rar_trans_map,by="row.names")
nrow(SILVA_MERDS_rar_trans_ord_WU_points)
summary(SILVA_MERDS_rar_trans_ord_WU_points)
(trans_ord_p=ggplot(data=SILVA_MERDS_rar_trans_ord_WU_points, aes(x=MDS1,y=MDS2))+
    geom_point(size=5, aes(shape=soil_root, fill=precip))+
    theme_bw()+
    scale_shape_manual(values = c(21,23,24))+scale_fill_manual(values = c("turquoise","firebrick1","black"))+
    theme(axis.title = element_text(size = 20),axis.text = element_text(size = 18),legend.position = "none",
          panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

#900x750
#Live transplant
SILVA_MERDS_rar_live_trans=subset_samples(SILVA_MERDS_rar, soil_status!="S"&life_stage=="G"|life_stage=="Start")
SILVA_MERDS_rar_live_trans_map=sample_data(SILVA_MERDS_rar_live_trans)
summary(sample_data(SILVA_MERDS_rar_live_trans))
SILVA_MERDS_rar_live_trans_WU_ord=ordinate(SILVA_MERDS_rar_live_trans, method = "NMDS",distance = "wunifrac")
#*** Solution reached
#0.05617122  
SILVA_MERDS_rar_live_trans_ord_WU_points=merge(SILVA_MERDS_rar_live_trans_WU_ord$points,SILVA_MERDS_rar_live_trans_map,by="row.names")
nrow(SILVA_MERDS_rar_live_trans_ord_WU_points)
summary(SILVA_MERDS_rar_live_trans_ord_WU_points)
(live_trans_ord_p=ggplot(data=SILVA_MERDS_rar_live_trans_ord_WU_points, aes(x=MDS1,y=MDS2))+
    geom_point(size=5, aes(shape=root_association, fill=precip))+
    theme_bw()+
    scale_shape_manual(values = c(21,24))+scale_fill_manual(values = c("turquoise","firebrick1","black"))+
    theme(axis.title = element_text(size = 20),axis.text = element_text(size = 18),legend.position = "none",
          panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

#900x750
#Experiemental transplant
SILVA_MERDS_rar_exp_trans=subset_samples(SILVA_MERDS_rar, life_stage=="G")
SILVA_MERDS_rar_exp_trans_map=sample_data(SILVA_MERDS_rar_exp_trans)
summary(SILVA_MERDS_rar_exp_trans_map)
SILVA_MERDS_rar_exp_trans_WU_ord=ordinate(SILVA_MERDS_rar_exp_trans, method = "NMDS",distance = "wunifrac")
#*** Solution reached
#0.0533245  
SILVA_MERDS_rar_exp_trans_WU_ord_points=merge(SILVA_MERDS_rar_exp_trans_WU_ord$points,SILVA_MERDS_rar_exp_trans_map,by="row.names")
nrow(SILVA_MERDS_rar_exp_trans_WU_ord_points)
summary(SILVA_MERDS_rar_exp_trans_WU_ord_points)
(live_trans_ord_p=ggplot(data=SILVA_MERDS_rar_exp_trans_WU_ord_points, aes(x=MDS1,y=MDS2))+
    geom_point(size=5, aes(shape=soil_root, fill=precip))+
    theme_bw()+
    scale_shape_manual(values = c(21,23,24))+scale_fill_manual(values = c("turquoise","firebrick1"))+
    theme(axis.title = element_text(size = 20),axis.text = element_text(size = 18),legend.position = "none",
          panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

#Transplanted 
SILVA_MERDS_rar_trt=subset_samples(SILVA_MERDS_rar, life_stage!="Start")
SILVA_MERDS_rar_trt_live=subset_samples(SILVA_MERDS_rar_trt, soil_status!="S")
SILVA_MERDS_rar_trt_live_trans=subset_samples(SILVA_MERDS_rar_trt_live, life_stage=="G")

SILVA_MERDS_rar_trt_live_trans_ord=ordinate(SILVA_MERDS_rar_trt_live_trans, method = "NMDS",distance = "bray")
#*** Solution reached
#0.1607911
plot_ordination(SILVA_MERDS_rar_trt_live_trans,SILVA_MERDS_rar_trt_live_trans_ord, color="precip",shape="root_association")+geom_point(size=4)+
  geom_label_repel(size=3,aes(label = block))+theme_bw()


SILVA_MERDS_rar_trt_live_trans_map=sample_data(SILVA_MERDS_rar_trt_live_trans)
SILVA_MERDS_rar_trt_live_trans_map$soil_root=with(SILVA_MERDS_rar_trt_live_trans_map, interaction(soil_status,root_association))
SILVA_MERDS_rar_trt_live_transe_dis=distance(SILVA_MERDS_rar_trt_live_trans,method = "bray")

trans_perm=adonis(SILVA_MERDS_rar_trt_live_transe_dis~SILVA_MERDS_rar_trt_live_trans_map$soil_root*SILVA_MERDS_rar_trt_live_trans_map$precip
                  +as.factor(SILVA_MERDS_rar_trt_live_trans_map$block), permutations = 9999)
#SILVA_MERDS_rar_trt_live_trans_map$precip                                               1   0.26218 0.26218 2.13979 0.12422 0.0045 **
#as.factor(SILVA_MERDS_rar_trt_live_trans_map$block)                                     3   0.47584 0.15861 1.29453 0.22545 0.0808 . 
AICc.PERMANOVA(trans_perm)

trans_perm_no_block=adonis(SILVA_MERDS_rar_trt_live_transe_dis~SILVA_MERDS_rar_trt_live_trans_map$soil_root*SILVA_MERDS_rar_trt_live_trans_map$precip, 
                           permutations = 9999)
#SILVA_MERDS_rar_trt_live_trans_map$precip                                               1   0.26218 0.26218 1.99304 0.12422 0.0095 **
AICc.PERMANOVA(trans_perm_no_block)


#####Weighted Unifrac Live and Sterile Trans Comm Analyses####

SILVA_MERDS_rar_trt_trans_WU_ord=ordinate(SILVA_MERDS_rar_trt_trans, method = "NMDS",distance = "Wunifrac")
#*** Solution reached

#0.0533245
plot_ordination(SILVA_MERDS_rar_trt_trans,SILVA_MERDS_rar_trt_trans_WU_ord, color="precip",shape="soil_root")+geom_point(size=4)+
  geom_label_repel(size=3,aes(label = block))+theme_bw()


SILVA_MERDS_rar_trt_trans_map=sample_data(SILVA_MERDS_rar_trt_trans)
SILVA_MERDS_rar_trt_trans_map$soil_root=with(SILVA_MERDS_rar_trt_trans_map, interaction(soil_status,root_association))
SILVA_MERDS_rar_trt_trans_WU_dis=phyloseq::distance(SILVA_MERDS_rar_trt_trans,method = "Wunifrac")

adonis(SILVA_MERDS_rar_trt_trans_WU_dis~SILVA_MERDS_rar_trt_trans_map$soil_root*SILVA_MERDS_rar_trt_trans_map$precip
       +as.factor(SILVA_MERDS_rar_trt_trans_map$block), permutations = 9999)
#SILVA_MERDS_rar_trt_trans_map$soil_root                                       2   0.45684 0.228419 19.0401 0.58834 0.0001 ***
#SILVA_MERDS_rar_trt_trans_map$precip                                          1   0.04556 0.045556  3.7974 0.05867 0.0252 *  
#SILVA_MERDS_rar_trt_trans_map$soil_root:SILVA_MERDS_rar_trt_trans_map$precip  2   0.05295 0.026476  2.2069 0.06819 0.0733 .  

pairwise.perm.manova(SILVA_MERDS_rar_trt_trans_WU_dis,interaction(SILVA_MERDS_rar_trt_trans_map$soil_root,SILVA_MERDS_rar_trt_trans_map$precip),nperm = 9999)


#####Presence Weighted Unifrac Transplant####
head(sample_data(SILVA_MERDS_rar))

SILVA_MERDS_rar_pres_trans=subset_samples(SILVA_MERDS_rar, life_stage=="G"&life_stage!="Start"&
                                            root_association=="B")
nrow(sample_data(SILVA_MERDS_rar_pres_trans))
#16
unique(sample_data(SILVA_MERDS_rar_pres_trans)$soil_root)

SILVA_MERDS_rar_pres_trans_map=sample_data(SILVA_MERDS_rar_pres_trans)
SILVA_MERDS_rar_pres_trans_map$soil_root=with(SILVA_MERDS_rar_pres_trans_map, interaction(soil_status,root_association))
SILVA_MERDS_rar_pres_trans_WU_dis=distance(SILVA_MERDS_rar_pres_trans,method = "Wunifrac")

(pres_tran_WU=adonis(SILVA_MERDS_rar_pres_trans_WU_dis~SILVA_MERDS_rar_pres_trans_map$soil_root*SILVA_MERDS_rar_pres_trans_map$precip
                     +as.factor(SILVA_MERDS_rar_pres_trans_map$block), permutations = 9999))


AICc.PERMANOVA(pres_tran_WU)
#$AIC
#[1] -24.18841


pairwise.perm.manova(SILVA_MERDS_rar_pres_trans_WU_dis, SILVA_MERDS_rar_map$soil_root_stage, nperm=2000)


(pres_tran_WU_no_block=adonis(SILVA_MERDS_rar_pres_trans_WU_dis~SILVA_MERDS_rar_pres_trans_map$soil_root*
                                SILVA_MERDS_rar_pres_trans_map$precip, permutations = 9999))


AICc.PERMANOVA(pres_tran_WU_no_block)
#$AIC
#[1] -23.65435


#####Origin Weighted Unifrac Transplant####
head(sample_data(SILVA_MERDS_rar))

SILVA_MERDS_rar_orig_trans=subset_samples(SILVA_MERDS_rar, life_stage=="G"&life_stage!="Start"&
                                            soil_status =="L")
nrow(sample_data(SILVA_MERDS_rar_orig_trans))
#16
unique(sample_data(SILVA_MERDS_rar_orig_trans)$soil_root)

SILVA_MERDS_rar_orig_trans_map=sample_data(SILVA_MERDS_rar_orig_trans)
SILVA_MERDS_rar_orig_trans_map$soil_root=with(SILVA_MERDS_rar_orig_trans_map, interaction(soil_status,root_association))
SILVA_MERDS_rar_orig_trans_WU_dis=distance(SILVA_MERDS_rar_orig_trans,method = "Wunifrac")

(orig_tran_WU=adonis(SILVA_MERDS_rar_orig_trans_WU_dis~SILVA_MERDS_rar_orig_trans_map$soil_root*SILVA_MERDS_rar_orig_trans_map$precip
                     +as.factor(SILVA_MERDS_rar_orig_trans_map$block), permutations = 9999))


AICc.PERMANOVA(orig_tran_WU)
#$AIC
#[1] -21.24361


pairwise.perm.manova(SILVA_MERDS_rar_orig_trans_WU_dis, SILVA_MERDS_rar_map$soil_root_stage, nperm=2000)


(orig_tran_WU_no_block=adonis(SILVA_MERDS_rar_orig_trans_WU_dis~SILVA_MERDS_rar_orig_trans_map$soil_root*
                                SILVA_MERDS_rar_orig_trans_map$precip, permutations = 9999))


AICc.PERMANOVA(orig_tran_WU_no_block)
#$AIC
#[1] -23.33317

#SIMPER Analyses

#need the OTU table

SILVA_MERDS_rar_trt_live_trans=prune_taxa(taxa_sums(SILVA_MERDS_rar_trt_live_trans) > 0, SILVA_MERDS_rar_trt_live_trans)
ntaxa(SILVA_MERDS_rar_trt_live_trans)
#5499

SILVA_MERDS_rar_trt_live_trans_OTU=t(otu_table(SILVA_MERDS_rar_trt_live_trans))
colnames(SILVA_MERDS_rar_trt_live_trans_OTU)
row.names(SILVA_MERDS_rar_trt_live_trans_OTU)
length(colnames(SILVA_MERDS_rar_trt_live_trans_OTU))
#5499
length(row.names(SILVA_MERDS_rar_trt_live_trans_OTU))
#16

nrow(SILVA_MERDS_rar_trt_live_trans_map)
#16


#life stage
SILVA_MERDS_rar_trt_live_trans_precip.simp <- with(SILVA_MERDS_rar_trt_live_trans_map, simper(SILVA_MERDS_rar_trt_live_trans_OTU, precip,permutations=999))
summary(SILVA_MERDS_rar_trt_live_trans_precip.simp,ordered = T)
SILVA_MERDS_rar_trt_live_trans_precip.simp_mat_num=as.data.frame(cbind(as.numeric(SILVA_MERDS_rar_trt_live_trans_precip.simp$D_A$average),as.numeric(SILVA_MERDS_rar_trt_live_trans_precip.simp$D_A$ava),
                                                                       as.numeric(SILVA_MERDS_rar_trt_live_trans_precip.simp$D_A$avb),as.numeric(SILVA_MERDS_rar_trt_live_trans_precip.simp$D_A$p)))

head(SILVA_MERDS_rar_trt_live_trans_precip.simp_mat_num)
row.names(SILVA_MERDS_rar_trt_live_trans_precip.simp_mat_num)=SILVA_MERDS_rar_trt_live_trans_precip.simp$D_A$species
summary(SILVA_MERDS_rar_trt_live_trans_precip.simp_mat_num)
colnames(SILVA_MERDS_rar_trt_live_trans_precip.simp_mat_num)[c(1:4)]=c("average","av_Drought","av_Ambient","pval")


#add in the taxonomy 
SILVA_MERDS_rar_trt_live_trans_tax=tax_table(SILVA_MERDS_rar_trt_live_trans)

SILVA_MERDS_rar_trt_live_trans_precip.simp_mat=merge(SILVA_MERDS_rar_trt_live_trans_precip.simp_mat_num,SILVA_MERDS_rar_trt_live_trans_tax, by="row.names", all.x = T)
head(SILVA_MERDS_rar_trt_live_trans_precip.simp_mat)
SILVA_MERDS_rar_trt_live_trans_precip.simp_mat$FDR_adj=p.adjust(SILVA_MERDS_rar_trt_live_trans_precip.simp_mat$pval,method = "fdr")
SILVA_MERDS_rar_trt_live_trans_precip.simp_mat_sig=subset(SILVA_MERDS_rar_trt_live_trans_precip.simp_mat, pval<0.05)
nrow(SILVA_MERDS_rar_trt_live_trans_precip.simp_mat_sig)
#133

colnames(SILVA_MERDS_rar_trt_live_trans_precip.simp_mat_sig)[1]="OTU"
write.csv(SILVA_MERDS_rar_trt_live_trans_precip.simp_mat_sig, "D:/MERDS_2018/merds/Switchgrass/R_data/SILVA_MERDS_rar_trt_live_trans_precip.simp_mat_sig.csv")


#How abundanqt are these taxa in the community
SILVA_MERDS_rar_trt_trans=subset_samples(SILVA_MERDS_rar, life_stage=="G"&life_stage!="Start")
#No sterile
SILVA_MERDS_rar_trt_live_trans=subset_samples(SILVA_MERDS_rar_trt_trans, soil_root!="S.B")
nrow(sample_data(SILVA_MERDS_rar_trt_live_trans))
#16

SILVA_MERDS_rar_trt_live_trans_precip.simp_mat_sig=read.csv("D:/MERDS_2018/merds/Switchgrass/R_data/SILVA_MERDS_rar_trt_live_trans_precip.simp_mat_sig.csv")
head(SILVA_MERDS_rar_trt_live_trans_precip.simp_mat_sig)
summary(SILVA_MERDS_rar_trt_live_trans_precip.simp_mat_sig)

SILVA_MERDS_rar_trt_live_trans_simp=prune_taxa(as.character(SILVA_MERDS_rar_trt_live_trans_precip.simp_mat_sig$OTU),SILVA_MERDS_rar_trt_live_trans)
ntaxa(SILVA_MERDS_rar_trt_live_trans_simp)


SILVA_MERDS_rar_trt_live_trans_simp_df=data.frame(taxa_sums(SILVA_MERDS_rar_trt_live_trans_simp))
colnames(SILVA_MERDS_rar_trt_live_trans_simp_df)="taxa_sum"

SILVA_MERDS_rar_trt_live_trans_precip.simp_mat_sig_sum=merge(SILVA_MERDS_rar_trt_live_trans_precip.simp_mat_sig,SILVA_MERDS_rar_trt_live_trans_simp_df, by.x="OTU", by.y="row.names")
head(SILVA_MERDS_rar_trt_live_trans_precip.simp_mat_sig_sum)

SILVA_MERDS_rar_trt_live_trans_precip.simp_mat_sig_sum %>% group_by(Phylum)%>%summarise_at("taxa_sum",~sum(.))
sum(taxa_sums(SILVA_MERDS_rar_trt_live_trans))

#Actino
2592/160000


#Germination


SILVA_MERDS_rar_seed=subset_samples(SILVA_MERDS_rar, life_stage=="S"|life_stage=="Start")
SILVA_MERDS_rar_seed_map=sample_data(SILVA_MERDS_rar_seed)
summary(sample_data(SILVA_MERDS_rar_seed))
SILVA_MERDS_rar_seed_WU_ord=ordinate(SILVA_MERDS_rar_seed, method = "NMDS",distance = "wunifrac")
#*** Solution reached
#0.05767356  
SILVA_MERDS_rar_seed_WU_ord_points=merge(SILVA_MERDS_rar_seed_WU_ord$points,SILVA_MERDS_rar_seed_map,by="row.names")
nrow(SILVA_MERDS_rar_seed_WU_ord_points)
summary(SILVA_MERDS_rar_seed_WU_ord_points)
(seed_ord_p=ggplot(data=SILVA_MERDS_rar_seed_WU_ord_points, aes(x=MDS1,y=MDS2))+
    geom_point(size=5, aes(shape=soil_root, fill=precip))+
    theme_bw()+
    scale_shape_manual(values = c(21,23,24))+scale_fill_manual(values = c("turquoise","firebrick1","black"))+
    theme(axis.title = element_text(size = 20),axis.text = element_text(size = 18),legend.position = "none",
          panel.grid.major = element_blank(), panel.grid.minor = element_blank()))
#nmds_germination_w_sterile_start

#Live germination
SILVA_MERDS_rar_live_seed=subset_samples(SILVA_MERDS_rar, soil_status!="S"&life_stage=="S"|life_stage=="Start")
SILVA_MERDS_rar_live_seed_map=sample_data(SILVA_MERDS_rar_live_seed)
summary(sample_data(SILVA_MERDS_rar_live_seed))
SILVA_MERDS_rar_live_seed_WU_ord=ordinate(SILVA_MERDS_rar_live_seed, method = "NMDS",distance = "wunifrac")
#*** Solution reached
#0.05023149   
SILVA_MERDS_rar_live_seed_ord_WU_points=merge(SILVA_MERDS_rar_live_seed_WU_ord$points,SILVA_MERDS_rar_live_seed_map,by="row.names")
nrow(SILVA_MERDS_rar_live_seed_ord_WU_points)
summary(SILVA_MERDS_rar_live_seed_ord_WU_points)
(live_seed_ord_p=ggplot(data=SILVA_MERDS_rar_live_seed_ord_WU_points, aes(x=MDS1,y=MDS2))+
    geom_point(size=5, aes(shape=root_association, fill=precip))+
    theme_bw()+
    scale_shape_manual(values = c(21,24))+scale_fill_manual(values = c("turquoise","firebrick1","black"))+
    theme(axis.title = element_text(size = 20),axis.text = element_text(size = 18),legend.position = "none",
          panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

#nmds_germination_w_start
#Experiemental germination
SILVA_MERDS_rar_exp_seed=subset_samples(SILVA_MERDS_rar, life_stage=="S")
SILVA_MERDS_rar_exp_seed_map=sample_data(SILVA_MERDS_rar_exp_seed)
summary(SILVA_MERDS_rar_exp_seed_map)
SILVA_MERDS_rar_exp_seed_WU_ord=ordinate(SILVA_MERDS_rar_exp_seed, method = "NMDS",distance = "wunifrac")
#*** Solution reached
#0.06719555  
SILVA_MERDS_rar_exp_seed_WU_ord_points=merge(SILVA_MERDS_rar_exp_seed_WU_ord$points,SILVA_MERDS_rar_exp_seed_map,by="row.names")
nrow(SILVA_MERDS_rar_exp_seed_WU_ord_points)
summary(SILVA_MERDS_rar_exp_seed_WU_ord_points)
(live_seed_ord_p=ggplot(data=SILVA_MERDS_rar_exp_seed_WU_ord_points, aes(x=MDS1,y=MDS2))+
    geom_point(size=5, aes(shape=soil_root, fill=precip))+
    theme_bw()+
    scale_shape_manual(values = c(21,23,24))+scale_fill_manual(values = c("turquoise","firebrick1"))+
    theme(axis.title = element_text(size = 20),axis.text = element_text(size = 18),legend.position = "none",
          panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

#seed
SILVA_MERDS_rar_trt=subset_samples(SILVA_MERDS_rar, life_stage!="Start")
SILVA_MERDS_rar_trt_live=subset_samples(SILVA_MERDS_rar_trt, soil_status!="S")
SILVA_MERDS_rar_trt_live_seed=subset_samples(SILVA_MERDS_rar_trt_live, life_stage=="S")

SILVA_MERDS_rar_trt_live_seed_ord=ordinate(SILVA_MERDS_rar_trt_live_seed, method = "NMDS",distance = "bray")
#*** Solution reached
#0.1359863 
plot_ordination(SILVA_MERDS_rar_trt_live_seed,SILVA_MERDS_rar_trt_live_seed_ord, color="precip",shape="root_association")+geom_point(size=4)+
  geom_label_repel(size=3,aes(label = block))+theme_bw()


SILVA_MERDS_rar_trt_live_seed_map=sample_data(SILVA_MERDS_rar_trt_live_seed)
SILVA_MERDS_rar_trt_live_seed_map$soil_root=with(SILVA_MERDS_rar_trt_live_seed_map, interaction(soil_status,root_association))
SILVA_MERDS_rar_trt_live_seed_dis=distance(SILVA_MERDS_rar_trt_live_seed,method = "bray")

seed_perm=adonis(SILVA_MERDS_rar_trt_live_seed_dis~SILVA_MERDS_rar_trt_live_seed_map$soil_root*SILVA_MERDS_rar_trt_live_seed_map$precip
                 +as.factor(SILVA_MERDS_rar_trt_live_seed_map$block), permutations = 9999)
#SILVA_MERDS_rar_trt_live_seed_map$soil_root                                           1   0.38611 0.38611  3.3937 0.16782 0.0001 ***
#SILVA_MERDS_rar_trt_live_seed_map$precip                                              1   0.33423 0.33423  2.9378 0.14527 0.0005 ***
AICc.PERMANOVA(seed_perm)

seed_perm_no_block=adonis(SILVA_MERDS_rar_trt_live_seed_dis~SILVA_MERDS_rar_trt_live_seed_map$soil_root*SILVA_MERDS_rar_trt_live_seed_map$precip,
                          permutations = 9999)
#SILVA_MERDS_rar_trt_live_seed_map$soil_root                                           1   0.38611 0.38611  3.2748 0.16782 0.0001 ***
#SILVA_MERDS_rar_trt_live_seed_map$precip                                              1   0.33423 0.33423  2.8348 0.14527 0.0006 ***

AICc.PERMANOVA(seed_perm_no_block)

#####Seeds####

SILVA_MERDS_rar_trt_seed=subset_samples(SILVA_MERDS_rar, life_stage=="S"&life_stage!="Start")
nrow(sample_data(SILVA_MERDS_rar_trt_seed))
#23


#####Weighted Unifrac Live and Sterile seed Comm Analyses####

SILVA_MERDS_rar_trt_seed_WU_ord=ordinate(SILVA_MERDS_rar_trt_seed, method = "NMDS",distance = "Wunifrac")
#*** Solution reached

#0.06719555
plot_ordination(SILVA_MERDS_rar_trt_seed,SILVA_MERDS_rar_trt_seed_WU_ord, color="precip",shape="soil_root")+geom_point(size=4)+
  geom_label_repel(size=3,aes(label = block))+theme_bw()


SILVA_MERDS_rar_trt_seed_map=sample_data(SILVA_MERDS_rar_trt_seed)
SILVA_MERDS_rar_trt_seed_map$soil_root=with(SILVA_MERDS_rar_trt_seed_map, interaction(soil_status,root_association))
SILVA_MERDS_rar_trt_seed_WU_dis=distance(SILVA_MERDS_rar_trt_seed,method = "Wunifrac")

adonis(SILVA_MERDS_rar_trt_seed_WU_dis~SILVA_MERDS_rar_trt_seed_map$soil_root*SILVA_MERDS_rar_trt_seed_map$precip
       +as.factor(SILVA_MERDS_rar_trt_seed_map$block), permutations = 9999)
#SILVA_MERDS_rar_trt_seed_map$soil_root                                      2   0.41316 0.206579 15.9172 0.58087 0.0001 ***
#SILVA_MERDS_rar_trt_seed_map$precip                                         1   0.04542 0.045418  3.4995 0.06385 0.0243 *  



#####Presence Weighted Unifrac Seed####
head(sample_data(SILVA_MERDS_rar))

SILVA_MERDS_rar_pres_seed=subset_samples(SILVA_MERDS_rar, life_stage=="S"&life_stage!="Start"&
                                           root_association=="B")
nrow(sample_data(SILVA_MERDS_rar_pres_seed))
#15
unique(sample_data(SILVA_MERDS_rar_pres_seed)$soil_root)

SILVA_MERDS_rar_pres_seed_map=sample_data(SILVA_MERDS_rar_pres_seed)
SILVA_MERDS_rar_pres_seed_map$soil_root=with(SILVA_MERDS_rar_pres_seed_map, interaction(soil_status,root_association))
SILVA_MERDS_rar_pres_seed_WU_dis=distance(SILVA_MERDS_rar_pres_seed,method = "Wunifrac")

(pres_seed_WU=adonis(SILVA_MERDS_rar_pres_seed_WU_dis~SILVA_MERDS_rar_pres_seed_map$soil_root*SILVA_MERDS_rar_pres_seed_map$precip
                     +as.factor(SILVA_MERDS_rar_pres_seed_map$block), permutations = 9999))


AICc.PERMANOVA(pres_seed_WU)
#$AIC
#[1] -18.45658


pairwise.perm.manova(SILVA_MERDS_rar_pres_seed_WU_dis, SILVA_MERDS_rar_map$soil_root_stage, nperm=2000)


(pres_seed_WU_no_block=adonis(SILVA_MERDS_rar_pres_seed_WU_dis~SILVA_MERDS_rar_pres_seed_map$soil_root*
                                SILVA_MERDS_rar_pres_seed_map$precip, permutations = 9999))


AICc.PERMANOVA(pres_seed_WU_no_block)
#$AIC
#[1] -19.09256


#####Origin Weighted Unifrac Seed####
head(sample_data(SILVA_MERDS_rar))

SILVA_MERDS_rar_orig_seed=subset_samples(SILVA_MERDS_rar, life_stage=="S"&life_stage!="Start"&
                                           soil_status =="L")
nrow(sample_data(SILVA_MERDS_rar_orig_seed))
#16
unique(sample_data(SILVA_MERDS_rar_orig_seed)$soil_root)

SILVA_MERDS_rar_orig_seed_map=sample_data(SILVA_MERDS_rar_orig_seed)
SILVA_MERDS_rar_orig_seed_map$soil_root=with(SILVA_MERDS_rar_orig_seed_map, interaction(soil_status,root_association))
SILVA_MERDS_rar_orig_seed_WU_dis=distance(SILVA_MERDS_rar_orig_seed,method = "Wunifrac")

(orig_seed_WU=adonis(SILVA_MERDS_rar_orig_seed_WU_dis~SILVA_MERDS_rar_orig_seed_map$soil_root*SILVA_MERDS_rar_orig_seed_map$precip
                     +as.factor(SILVA_MERDS_rar_orig_seed_map$block), permutations = 9999))


AICc.PERMANOVA(orig_seed_WU)
#$AIC
#[1] -24.17012


pairwise.perm.manova(SILVA_MERDS_rar_orig_seed_WU_dis, SILVA_MERDS_rar_map$soil_root_stage, nperm=2000)


(orig_seed_WU_no_block=adonis(SILVA_MERDS_rar_orig_seed_WU_dis~SILVA_MERDS_rar_orig_seed_map$soil_root*
                                SILVA_MERDS_rar_orig_seed_map$precip, permutations = 9999))


AICc.PERMANOVA(orig_seed_WU_no_block)
#$AIC
#[1] -25.56769


#SIMPER Analyses

#need the OTU table

SILVA_MERDS_rar_trt_live_seed=prune_taxa(taxa_sums(SILVA_MERDS_rar_trt_live_seed) > 0, SILVA_MERDS_rar_trt_live_seed)
ntaxa(SILVA_MERDS_rar_trt_live_seed)
#5229

SILVA_MERDS_rar_trt_live_seed_OTU=t(otu_table(SILVA_MERDS_rar_trt_live_seed))
colnames(SILVA_MERDS_rar_trt_live_seed_OTU)
row.names(SILVA_MERDS_rar_trt_live_seed_OTU)
length(colnames(SILVA_MERDS_rar_trt_live_seed_OTU))
#5229
length(row.names(SILVA_MERDS_rar_trt_live_seed_OTU))
#16

nrow(SILVA_MERDS_rar_trt_live_seed_map)
#16


#Soil association
SILVA_MERDS_rar_trt_live_seed_soil.simp <- with(SILVA_MERDS_rar_trt_live_seed_map, simper(SILVA_MERDS_rar_trt_live_seed_OTU, soil_root ,permutations=999))
summary(SILVA_MERDS_rar_trt_live_seed_soil.simp,ordered = T)
SILVA_MERDS_rar_trt_live_seed_soil.simp_mat_num=as.data.frame(cbind(as.numeric(SILVA_MERDS_rar_trt_live_seed_soil.simp$L.R_L.B$average),as.numeric(SILVA_MERDS_rar_trt_live_seed_soil.simp$L.R_L.B$ava),
                                                                    as.numeric(SILVA_MERDS_rar_trt_live_seed_soil.simp$L.R_L.B$avb),as.numeric(SILVA_MERDS_rar_trt_live_seed_soil.simp$L.R_L.B$p)))

head(SILVA_MERDS_rar_trt_live_seed_soil.simp_mat_num)
row.names(SILVA_MERDS_rar_trt_live_seed_soil.simp_mat_num)=SILVA_MERDS_rar_trt_live_seed_soil.simp$L.R_L.B$species
summary(SILVA_MERDS_rar_trt_live_seed_soil.simp_mat_num)
colnames(SILVA_MERDS_rar_trt_live_seed_soil.simp_mat_num)[c(1:4)]=c("average","av_Rhizosphere","av_Bulk","pval")
SILVA_MERDS_rar_trt_live_seed_soil.simp_mat_num$tax=row.names(SILVA_MERDS_rar_trt_live_seed_soil.simp_mat_num)

#add in the taxonomy 
SILVA_MERDS_rar_trt_live_seed_tax=as.data.frame(tax_table(SILVA_MERDS_rar_trt_live_seed))
nrow(SILVA_MERDS_rar_trt_live_seed_tax)
#5229
SILVA_MERDS_rar_trt_live_seed_tax$tax=row.names(SILVA_MERDS_rar_trt_live_seed_tax)

attributes(SILVA_MERDS_rar_trt_live_seed_tax)
SILVA_MERDS_rar_trt_live_seed_soil.simp_mat=join(SILVA_MERDS_rar_trt_live_seed_soil.simp_mat_num,SILVA_MERDS_rar_trt_live_seed_tax, by="tax", type = "left")

head(SILVA_MERDS_rar_trt_live_seed_soil.simp_mat)
SILVA_MERDS_rar_trt_live_seed_soil.simp_mat$FDR_adj=p.adjust(SILVA_MERDS_rar_trt_live_seed_soil.simp_mat$pval,method = "fdr")
SILVA_MERDS_rar_trt_live_seed_soil.simp_mat_sig=subset(SILVA_MERDS_rar_trt_live_seed_soil.simp_mat, pval<0.05)
nrow(SILVA_MERDS_rar_trt_live_seed_soil.simp_mat_sig)
#128

colnames(SILVA_MERDS_rar_trt_live_seed_soil.simp_mat_sig)[5]="OTU"
write.csv(SILVA_MERDS_rar_trt_live_seed_soil.simp_mat_sig, "D:/MERDS_2018/merds/Switchgrass/R_data/SILVA_MERDS_rar_trt_live_seed_soil_root.simp_mat_sig.csv")

#precip
SILVA_MERDS_rar_trt_live_seed_precip.simp <- with(SILVA_MERDS_rar_trt_live_seed_map, simper(SILVA_MERDS_rar_trt_live_seed_OTU, precip,permutations=999))
summary(SILVA_MERDS_rar_trt_live_seed_precip.simp,ordered = T)
SILVA_MERDS_rar_trt_live_seed_precip.simp_mat_num=as.data.frame(cbind(as.numeric(SILVA_MERDS_rar_trt_live_seed_precip.simp$D_A$average),as.numeric(SILVA_MERDS_rar_trt_live_seed_precip.simp$D_A$ava),
                                                                      as.numeric(SILVA_MERDS_rar_trt_live_seed_precip.simp$D_A$avb),as.numeric(SILVA_MERDS_rar_trt_live_seed_precip.simp$D_A$p)))

head(SILVA_MERDS_rar_trt_live_seed_precip.simp_mat_num)
row.names(SILVA_MERDS_rar_trt_live_seed_precip.simp_mat_num)=SILVA_MERDS_rar_trt_live_seed_precip.simp$D_A$species
summary(SILVA_MERDS_rar_trt_live_seed_precip.simp_mat_num)
colnames(SILVA_MERDS_rar_trt_live_seed_precip.simp_mat_num)[c(1:4)]=c("average","av_Drought","av_Ambient","pval")


#add in the taxonomy 
SILVA_MERDS_rar_trt_live_seed_tax=tax_table(SILVA_MERDS_rar_trt_live_seed)

SILVA_MERDS_rar_trt_live_seed_precip.simp_mat=merge(SILVA_MERDS_rar_trt_live_seed_precip.simp_mat_num,SILVA_MERDS_rar_trt_live_seed_tax, by="row.names", all.x = T)
head(SILVA_MERDS_rar_trt_live_seed_precip.simp_mat)
SILVA_MERDS_rar_trt_live_seed_precip.simp_mat$FDR_adj=p.adjust(SILVA_MERDS_rar_trt_live_seed_precip.simp_mat$pval,method = "fdr")
SILVA_MERDS_rar_trt_live_seed_precip.simp_mat_sig=subset(SILVA_MERDS_rar_trt_live_seed_precip.simp_mat, pval<0.05)
nrow(SILVA_MERDS_rar_trt_live_seed_precip.simp_mat_sig)
#93

colnames(SILVA_MERDS_rar_trt_live_seed_precip.simp_mat_sig)[1]="OTU"
write.csv(SILVA_MERDS_rar_trt_live_seed_precip.simp_mat_sig, "D:/MERDS_2018/merds/Switchgrass/R_data/SILVA_MERDS_rar_trt_live_seed_precip.simp_mat_sig.csv")

#How abundanqt are these taxa in the community
SILVA_MERDS_rar_trt_seed=subset_samples(SILVA_MERDS_rar, life_stage=="S"&life_stage!="Start")
#No sterile
SILVA_MERDS_rar_trt_live_seed=subset_samples(SILVA_MERDS_rar_trt_seed, soil_root!="S.B")
nrow(sample_data(SILVA_MERDS_rar_trt_live_seed))
#16

SILVA_MERDS_rar_trt_live_seed_precip.simp_mat_sig=read.csv("D:/MERDS_2018/merds/Switchgrass/R_data/SILVA_MERDS_rar_trt_live_seed_precip.simp_mat_sig.csv")
head(SILVA_MERDS_rar_trt_live_seed_precip.simp_mat_sig)
summary(SILVA_MERDS_rar_trt_live_seed_precip.simp_mat_sig)

SILVA_MERDS_rar_trt_live_seed_simp=prune_taxa(as.character(SILVA_MERDS_rar_trt_live_seed_precip.simp_mat_sig$OTU),SILVA_MERDS_rar_trt_live_seed)
ntaxa(SILVA_MERDS_rar_trt_live_seed_simp)


SILVA_MERDS_rar_trt_live_seed_simp_df=data.frame(taxa_sums(SILVA_MERDS_rar_trt_live_seed_simp))
colnames(SILVA_MERDS_rar_trt_live_seed_simp_df)="taxa_sum"

SILVA_MERDS_rar_trt_live_seed_precip.simp_mat_sig_sum=merge(SILVA_MERDS_rar_trt_live_seed_precip.simp_mat_sig,SILVA_MERDS_rar_trt_live_seed_simp_df, by.x="OTU", by.y="row.names")
head(SILVA_MERDS_rar_trt_live_seed_precip.simp_mat_sig_sum)

SILVA_MERDS_rar_trt_live_seed_precip.simp_mat_sig_sum %>% group_by(Phylum)%>%summarise_at("taxa_sum",~sum(.))
sum(taxa_sums(SILVA_MERDS_rar_trt_live_seed))

#Actino
2768/160000
#0.0173

#####Sterile Analyses ####
#Sterile 
SILVA_MERDS_rar_sterile=subset_samples(SILVA_MERDS_rar, soil_status=="S")
ntaxa(SILVA_MERDS_rar_sterile)
#8366
SILVA_MERDS_rar_sterile=prune_taxa(taxa_sums(SILVA_MERDS_rar_sterile) > 0, SILVA_MERDS_rar_sterile)
ntaxa(SILVA_MERDS_rar_sterile)
#930
SILVA_MERDS_rar_sterile_map=sample_data(SILVA_MERDS_rar_sterile)
SILVA_MERDS_rar_sterile_ord=ordinate(SILVA_MERDS_rar_sterile, method = "NMDS",distance = "bray")
#*** Solution reached
#0.1223411
plot_ordination(SILVA_MERDS_rar_sterile,SILVA_MERDS_rar_sterile_ord, color="precip",shape="life_stage")+geom_point(size=4)+
  geom_label_repel(size=3,aes(label = block))+theme_bw()

SILVA_MERDS_rar_sterile_map=sample_data(SILVA_MERDS_rar_sterile)
SILVA_MERDS_rar_sterile_map$soil_root=with(SILVA_MERDS_rar_sterile_map, interaction(soil_status,root_association))
SILVA_MERDS_rar_sterile_dis=distance(SILVA_MERDS_rar_sterile,method = "bray")

sterile_perm=adonis(SILVA_MERDS_rar_sterile_dis~SILVA_MERDS_rar_sterile_map$precip*SILVA_MERDS_rar_sterile_map$life_stage
                    +as.factor(SILVA_MERDS_rar_sterile_map$block), permutations = 9999)

"Terms added sequentially (first to last)

                                                                          Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
SILVA_MERDS_rar_sterile_map$precip                                         1   0.38409 0.38409  2.8619 0.13372 0.0072
SILVA_MERDS_rar_sterile_map$life_stage                                     1   0.72564 0.72564  5.4069 0.25263 0.0001
as.factor(SILVA_MERDS_rar_sterile_map$block)                               3   0.43799 0.14600  1.0878 0.15249 0.3466
SILVA_MERDS_rar_sterile_map$precip:SILVA_MERDS_rar_sterile_map$life_stage  1   0.25094 0.25094  1.8698 0.08737 0.0581
Residuals                                                                  8   1.07367 0.13421         0.37380       
Total                                                                     14   2.87234                 1.00000         
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 "

AICc.PERMANOVA(sterile_perm)

sterile_perm_no_block=adonis(SILVA_MERDS_rar_sterile_dis~SILVA_MERDS_rar_sterile_map$precip*SILVA_MERDS_rar_sterile_map$life_stage, 
                             permutations = 9999)

"Terms added sequentially (first to last)

                                                                          Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
SILVA_MERDS_rar_sterile_map$precip                                         1   0.38409 0.38409  2.7180 0.13372 0.0086 ** 
SILVA_MERDS_rar_sterile_map$life_stage                                     1   0.72564 0.72564  5.1349 0.25263 0.0001 ***
SILVA_MERDS_rar_sterile_map$precip:SILVA_MERDS_rar_sterile_map$life_stage  1   0.20814 0.20814  1.4729 0.07246 0.1328    
Residuals                                                                 11   1.55446 0.14131         0.54118           
Total                                                                     14   2.87234                 1.00000  "

AICc.PERMANOVA(sterile_perm_no_block)



#####Sterile SIMPER Analyses####

#need the OTU table
SILVA_MERDS_rar_sterile_otu=t(otu_table(SILVA_MERDS_rar_sterile))

colnames(SILVA_MERDS_rar_sterile_otu)
row.names(SILVA_MERDS_rar_sterile_otu)

#life stage
SILVA_MERDS_rar_sterile.simp <- with(SILVA_MERDS_rar_sterile_map, simper(SILVA_MERDS_rar_sterile_otu, life_stage,permutations=999))
summary(SILVA_MERDS_rar_sterile.simp,ordered = T)
SILVA_MERDS_rar_sterile.simp_LS_mat_num=as.data.frame(cbind(as.numeric(SILVA_MERDS_rar_sterile.simp$G_S$average),as.numeric(SILVA_MERDS_rar_sterile.simp$G_S$ava),
                                                            as.numeric(SILVA_MERDS_rar_sterile.simp$G_S$avb),as.numeric(SILVA_MERDS_rar_sterile.simp$G_S$p)))

head(SILVA_MERDS_rar_sterile.simp_LS_mat_num)
row.names(SILVA_MERDS_rar_sterile.simp_LS_mat_num)=SILVA_MERDS_rar_sterile.simp$G_S$species
summary(SILVA_MERDS_rar_sterile.simp_LS_mat_num)
colnames(SILVA_MERDS_rar_sterile.simp_LS_mat_num)[c(1:4)]=c("average","av_Trans","av_Seed","pval")


#add in the taxonomy 
SILVA_MERDS_rar_sterile_tax=tax_table(SILVA_MERDS_rar_sterile)

SILVA_MERDS_rar_sterile.simp_LS_mat=merge(SILVA_MERDS_rar_sterile.simp_LS_mat_num,SILVA_MERDS_rar_sterile_tax, by="row.names", all.x = T)
head(SILVA_MERDS_rar_sterile.simp_LS_mat)
SILVA_MERDS_rar_sterile.simp_LS_mat$FDR_adj=p.adjust(SILVA_MERDS_rar_sterile.simp_LS_mat$pval,method = "fdr")
SILVA_MERDS_rar_sterile.simp_LS_mat_sig=subset(SILVA_MERDS_rar_sterile.simp_LS_mat, pval<0.05)
nrow(SILVA_MERDS_rar_sterile.simp_LS_mat_sig)
#59

colnames(SILVA_MERDS_rar_sterile.simp_LS_mat_sig)[1]="OTU"


write.csv(SILVA_MERDS_rar_sterile.simp_LS_mat_sig, "D:/MERDS_2018/merds/Switchgrass/R_data/SILVA_MERDS_rar_sterile.simp_life_stage_mat_sig.csv")



#Precip treatment
SILVA_MERDS_rar_sterile_precip.simp <- with(SILVA_MERDS_rar_sterile_map, simper(SILVA_MERDS_rar_sterile_otu, precip,permutations=999))
summary(SILVA_MERDS_rar_sterile_precip.simp,ordered = T)
SILVA_MERDS_rar_sterile_precip.simp_mat_num=as.data.frame(cbind(as.numeric(SILVA_MERDS_rar_sterile_precip.simp$D_A$average),as.numeric(SILVA_MERDS_rar_sterile_precip.simp$D_A$ava),
                                                                as.numeric(SILVA_MERDS_rar_sterile_precip.simp$D_A$avb),as.numeric(SILVA_MERDS_rar_sterile_precip.simp$D_A$p)))

head(SILVA_MERDS_rar_sterile_precip.simp_mat_num)
row.names(SILVA_MERDS_rar_sterile_precip.simp_mat_num)=SILVA_MERDS_rar_sterile_precip.simp$D_A$species
summary(SILVA_MERDS_rar_sterile_precip.simp_mat_num)
colnames(SILVA_MERDS_rar_sterile_precip.simp_mat_num)[c(1:4)]=c("average","av_Drought","av_Ambient","pval")


#add in the taxonomy 
SILVA_MERDS_rar_sterile_tax=tax_table(SILVA_MERDS_rar_sterile)

SILVA_MERDS_rar_sterile_precip.simp_mat=merge(SILVA_MERDS_rar_sterile_precip.simp_mat_num,SILVA_MERDS_rar_sterile_tax, by="row.names", all.x = T)
head(SILVA_MERDS_rar_sterile_precip.simp_mat)
SILVA_MERDS_rar_sterile_precip.simp_mat$FDR_adj=p.adjust(SILVA_MERDS_rar_sterile_precip.simp_mat$pval,method = "fdr")
SILVA_MERDS_rar_sterile_precip.simp_sig=subset(SILVA_MERDS_rar_sterile_precip.simp_mat, pval<0.05)
nrow(SILVA_MERDS_rar_sterile_precip.simp_sig)
#30

colnames(SILVA_MERDS_rar_sterile_precip.simp_sig)[1]="OTU"
write.csv(SILVA_MERDS_rar_sterile_precip.simp_sig, "D:/MERDS_2018/merds/Switchgrass/R_data/SILVA_MERDS_rar_sterile.simp_precip_mat_sig.csv")


#####Monoderm Analyses####

#####Actinobacteria Analyses######
get_taxa_unique(SILVA_MERDS_rar, taxonomic.rank="Phylum")
SILVA_MERDS_rar_map=sample_data(SILVA_MERDS_rar)
#p:Actinobacteria

Actino_SILVA_MERDS_rar <- subset_taxa(SILVA_MERDS_rar, Phylum=="p:Actinobacteria")
get_taxa_unique(Actino_SILVA_MERDS_rar, taxonomic.rank="Phylum")
Actino_all.reads=sample_sums(Actino_SILVA_MERDS_rar)
all.reads=sample_sums(SILVA_MERDS_rar)
Actino_soil_all.reads=cbind(Actino_all.reads,all.reads)
colnames(Actino_soil_all.reads)<-c("Actino.reads","total.reads")
Actino_soil_all.reads=merge(Actino_soil_all.reads, SILVA_MERDS_rar_map, by ="row.names")
head(Actino_soil_all.reads)
Actino_soil_all.reads=mutate(Actino_soil_all.reads, Actino.prop=Actino.reads/total.reads)

Actino_soil_all.reads_live_NS=subset(Actino_soil_all.reads,life_stage!="Start"&soil_root!="S.B")
nrow(Actino_soil_all.reads_live_NS)
#32
#Prop of Actinobacteria
Actin_model_live_NS= lm((Actino.prop)~soil_root*life_stage*precip+as.factor(block), data= Actino_soil_all.reads_live_NS)
qqPlot(resid(Actin_model_live_NS))
hist(resid(Actin_model_live_NS))
boxCox(Actin_model_live_NS)
shapiro.test(resid(Actin_model_live_NS))
#0.765
Anova(Actin_model_live_NS, type=3)
#soil_root        0.041116  2  37.4430 1.472e-06 ***
#precip           0.006325  1  11.5192  0.004007 ** 
#soil_root:precip 0.005279  2   4.8078  0.024357 *

emmeans(Actin_model_live_NS, pairwise~soil_root)


emmeans(Actin_model_live_NS, pairwise~soil_root|life_stage)


Actin_model_live_NS_no_block= lm((Actino.prop)~soil_root*life_stage*precip, data= Actino_soil_all.reads_live_NS)
qqPlot(resid(Actin_model_live_NS_no_block))
hist(resid(Actin_model_live_NS_no_block))
boxCox(Actin_model_live_NS_no_block)
shapiro.test(resid(Actin_model_live_NS_no_block))
#0.8128
Anova(Actin_model_live_NS_no_block, type=3)
#soil_root                   0.00499  1   6.2139 0.0199731 * 
#precip                      0.01637  1  20.3709 0.0001431 ***
#soil_root:life_stage        0.00527  1   6.5542 0.0171814 * 
#soil_root:life_stage:precip 0.00284  1   3.5310 0.0724307 .  

emmeans(Actin_model_live_NS_no_block, pairwise~soil_root)


emmeans(Actin_model_live_NS_no_block, pairwise~soil_root*life_stage,adjust="fdr")

AIC(Actin_model_live_NS,Actin_model_live_NS_no_block)

Actino_soil_trans.reads %>% group_by(soil_root) %>% summarise_at("Actino.prop", funs(n(),mean,sd,se=sd(.)/sqrt(n())))


Actin_MERDS_rar_trt_trans_precip_soil_g=Actino_soil_trans.reads %>% group_by(soil_root,precip)
Actino_prop_precip_trans=summarise_at(Actin_MERDS_rar_trt_trans_precip_soil_g, 
                                      "Actino.prop", funs(n(),mean,sd,se=sd(.)/sqrt(n())))
treatment_order=c("S.B","L.B","L.R")
(trans_obs_Actin_p=ggplot(Actino_prop_precip_trans, aes(x=precip,y=mean, ymin = mean-se, ymax= mean+se,fill=factor(soil_root,levels = treatment_order)))+
    geom_bar(position="dodge",stat="identity", color="black",aes(fill=factor(soil_root,levels = treatment_order)))+
    geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+
    scale_fill_manual(values = c( "white","lightgray", "darkgrey"),labels=c("Sterile","Bulk","Rhizo"), limits=treatment_order, name=NULL)+
    scale_x_discrete(labels=c("Ambient","Drought"), name=NULL)+ylab("Actinobacteria reads (proportion)")+
    geom_text(aes(y=mean-se-.01, label=n),size=10,position=position_dodge(width=0.9))+theme_bw()+
    theme(axis.title = element_text(size = 23), axis.text = element_text(size = 23),
          legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

#actino_seedling_prop


#####Firmicutes Analyses######
get_taxa_unique(SILVA_MERDS_rar, taxonomic.rank="Phylum")
SILVA_MERDS_rar_map=sample_data(SILVA_MERDS_rar)
#p:Actinobacteria

Firm_SILVA_MERDS_rar <- subset_taxa(SILVA_MERDS_rar, Phylum=="p:Firmicutes")
get_taxa_unique(Firm_SILVA_MERDS_rar, taxonomic.rank="Phylum")
Firm_all.reads=sample_sums(Firm_SILVA_MERDS_rar)
all.reads=sample_sums(SILVA_MERDS_rar)
Firm_soil_all.reads=cbind(Firm_all.reads,all.reads)
colnames(Firm_soil_all.reads)<-c("Firm.reads","total.reads")
Firm_soil_all.reads=merge(Firm_soil_all.reads, SILVA_MERDS_rar_map, by ="row.names")
head(Firm_soil_all.reads)
Firm_soil_all.reads=mutate(Firm_soil_all.reads, Firm.prop=Firm.reads/total.reads)

Firm_soil_all.reads_live_NS=subset(Firm_soil_all.reads,life_stage!="Start"&soil_root!="S.B")
nrow(Firm_soil_all.reads_live_NS)
#32
#Prop of Firmicutes
Firm_model_live_NS= lm(logit(Firm.prop)~soil_root*life_stage*precip+as.factor(block), data= Firm_soil_all.reads_live_NS)
qqPlot(resid(Firm_model_live_NS))
hist(resid(Firm_model_live_NS))
boxCox(Firm_model_live_NS)
shapiro.test(resid(Firm_model_live_NS))
#0.2612
Anova(Firm_model_live_NS, type=3)
#soil_root                     0.73  1    3.4654  0.076722 .  
#life_stage                    0.65  1    3.0876  0.093464 .  
#precip                        2.09  1    9.8758  0.004915 ** 

emmeans(Firm_model_live_NS, pairwise~soil_root)


emmeans(Firm_model_live_NS, pairwise~soil_root|life_stage)


Firm_model_live_NS_no_block= lm(log(Firm.prop)~soil_root*life_stage*precip, data= Firm_soil_all.reads_live_NS)
qqPlot(resid(Firm_model_live_NS_no_block))
hist(resid(Firm_model_live_NS_no_block))
boxCox(Firm_model_live_NS_no_block)
shapiro.test(resid(Firm_model_live_NS_no_block))
#0.03282
Anova(Firm_model_live_NS_no_block, type=3)
#soil_root                     0.84  1    3.8643 0.06099 .  
#life_stage                    0.75  1    3.4581 0.07525 .  
#precip                        2.24  1   10.3387 0.00370 **   

emmeans(Firm_model_live_NS_no_block, pairwise~soil_root)


emmeans(Firm_model_live_NS_no_block, pairwise~soil_root*life_stage,adjust="fdr")

AIC(Firm_model_live_NS,Firm_model_live_NS_no_block)


#####Actinobacteria Analyses######
get_taxa_unique(SILVA_MERDS_rar, taxonomic.rank="Phylum")
SILVA_MERDS_rar_map=sample_data(SILVA_MERDS_rar)
#p:Actinobacteria

Actino_SILVA_MERDS_rar <- subset_taxa(SILVA_MERDS_rar, Phylum=="p:Actinobacteria")
get_taxa_unique(Actino_SILVA_MERDS_rar, taxonomic.rank="Phylum")
Actino_all.reads=sample_sums(Actino_SILVA_MERDS_rar)
all.reads=sample_sums(SILVA_MERDS_rar)
Actino_soil_all.reads=cbind(Actino_all.reads,all.reads)
colnames(Actino_soil_all.reads)<-c("Actino.reads","total.reads")
Actino_soil_all.reads=merge(Actino_soil_all.reads, SILVA_MERDS_rar_map, by ="row.names")
head(Actino_soil_all.reads)
Actino_soil_all.reads=mutate(Actino_soil_all.reads, Actino.prop=Actino.reads/total.reads)

Actino_soil_all.reads_live_NS=subset(Actino_soil_all.reads,life_stage!="Start"&soil_root!="S.B")
nrow(Actino_soil_all.reads_live_NS)
#32
#Prop of Actinobacteria
Actin_model_live_NS= lm((Actino.prop)~soil_root*life_stage*precip+as.factor(block), data= Actino_soil_all.reads_live_NS)
qqPlot(resid(Actin_model_live_NS))
hist(resid(Actin_model_live_NS))
boxCox(Actin_model_live_NS)
shapiro.test(resid(Actin_model_live_NS))
#0.765
Anova(Actin_model_live_NS, type=3)
#soil_root        0.041116  2  37.4430 1.472e-06 ***
#precip           0.006325  1  11.5192  0.004007 ** 
#soil_root:precip 0.005279  2   4.8078  0.024357 *

emmeans(Actin_model_live_NS, pairwise~soil_root)


emmeans(Actin_model_live_NS, pairwise~soil_root|life_stage)


Actin_model_live_NS_no_block= lm((Actino.prop)~soil_root*life_stage*precip, data= Actino_soil_all.reads_live_NS)
qqPlot(resid(Actin_model_live_NS_no_block))
hist(resid(Actin_model_live_NS_no_block))
boxCox(Actin_model_live_NS_no_block)
shapiro.test(resid(Actin_model_live_NS_no_block))
#0.8128
Anova(Actin_model_live_NS_no_block, type=3)
#soil_root                   0.00499  1   6.2139 0.0199731 * 
#precip                      0.01637  1  20.3709 0.0001431 ***
#soil_root:life_stage        0.00527  1   6.5542 0.0171814 * 
#soil_root:life_stage:precip 0.00284  1   3.5310 0.0724307 .  

emmeans(Actin_model_live_NS_no_block, pairwise~soil_root)


emmeans(Actin_model_live_NS_no_block, pairwise~soil_root*life_stage,adjust="fdr")

AIC(Actin_model_live_NS,Actin_model_live_NS_no_block)

Actino_soil_trans.reads %>% group_by(soil_root) %>% summarise_at("Actino.prop", funs(n(),mean,sd,se=sd(.)/sqrt(n())))


Actin_MERDS_rar_trt_trans_precip_soil_g=Actino_soil_trans.reads %>% group_by(soil_root,precip)
Actino_prop_precip_trans=summarise_at(Actin_MERDS_rar_trt_trans_precip_soil_g, 
                                      "Actino.prop", funs(n(),mean,sd,se=sd(.)/sqrt(n())))
treatment_order=c("S.B","L.B","L.R")
(trans_obs_Actin_p=ggplot(Actino_prop_precip_trans, aes(x=precip,y=mean, ymin = mean-se, ymax= mean+se,fill=factor(soil_root,levels = treatment_order)))+
    geom_bar(position="dodge",stat="identity", color="black",aes(fill=factor(soil_root,levels = treatment_order)))+
    geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+
    scale_fill_manual(values = c( "white","lightgray", "darkgrey"),labels=c("Sterile","Bulk","Rhizo"), limits=treatment_order, name=NULL)+
    scale_x_discrete(labels=c("Ambient","Drought"), name=NULL)+ylab("Actinobacteria reads (proportion)")+
    geom_text(aes(y=mean-se-.01, label=n),size=10,position=position_dodge(width=0.9))+theme_bw()+
    theme(axis.title = element_text(size = 23), axis.text = element_text(size = 23),
          legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

#actino_seedling_prop

#####Firmicutes Analyses######
get_taxa_unique(SILVA_MERDS_rar, taxonomic.rank="Phylum")
SILVA_MERDS_rar_map=sample_data(SILVA_MERDS_rar)
#p:Firmicutes

Firm_SILVA_MERDS_rar <- subset_taxa(SILVA_MERDS_rar, Phylum=="p:Firmicutes")
get_taxa_unique(Firm_SILVA_MERDS_rar, taxonomic.rank="Phylum")
Firm_all.reads=sample_sums(Firm_SILVA_MERDS_rar)
all.reads=sample_sums(SILVA_MERDS_rar)
Firm_soil_all.reads=cbind(Firm_all.reads,all.reads)
colnames(Firm_soil_all.reads)<-c("Firm.reads","total.reads")
Firm_soil_all.reads=merge(Firm_soil_all.reads, SILVA_MERDS_rar_map, by ="row.names")
head(Firm_soil_all.reads)
Firm_soil_all.reads=mutate(Firm_soil_all.reads, Firm.prop=Firm.reads/total.reads)

Firm_soil_all.reads_live_NS=subset(Firm_soil_all.reads,life_stage!="Start"&soil_root!="S.B")
nrow(Firm_soil_all.reads_live_NS)
#32
#Prop of Firmicutes
Firm_model_live_NS= lm((Firm.prop)~soil_root*life_stage*precip+as.factor(block), data= Firm_soil_all.reads_live_NS)
qqPlot(resid(Firm_model_live_NS))
hist(resid(Firm_model_live_NS))
boxCox(Firm_model_live_NS)
shapiro.test(resid(Firm_model_live_NS))
#0.2612
Anova(Firm_model_live_NS, type=3)
#soil_root                     0.73  1    3.4654  0.076722 .  
#life_stage                    0.65  1    3.0876  0.093464 .  
#precip                        2.09  1    9.8758  0.004915 ** 

emmeans(Firm_model_live_NS, pairwise~soil_root)


emmeans(Firm_model_live_NS, pairwise~soil_root|life_stage)


Firm_model_live_NS_no_block= lm(logit(Firm.prop)~soil_root*life_stage*precip, data= Firm_soil_all.reads_live_NS)
qqPlot(resid(Firm_model_live_NS_no_block))
hist(resid(Firm_model_live_NS_no_block))
boxCox(Firm_model_live_NS_no_block)
shapiro.test(resid(Firm_model_live_NS_no_block))
#0.02957
Anova(Firm_model_live_NS_no_block, type=3)
#soil_root                     0.87  1    3.8779 0.06058 .  
#life_stage                    0.78  1    3.4875 0.07409 .  
#precip                        2.33  1   10.3530 0.00368 **   

emmeans(Firm_model_live_NS_no_block, pairwise~soil_root)


emmeans(Firm_model_live_NS_no_block, pairwise~soil_root*life_stage,adjust="fdr")

AIC(Firm_model_live_NS,Firm_model_live_NS_no_block)




#####Chloroflexi Analyses######
get_taxa_unique(SILVA_MERDS_rar, taxonomic.rank="Phylum")
SILVA_MERDS_rar_map=sample_data(SILVA_MERDS_rar)
#p:Chloroflexi

Chloro_SILVA_MERDS_rar <- subset_taxa(SILVA_MERDS_rar, Phylum=="p:Chloroflexi")
get_taxa_unique(Chloro_SILVA_MERDS_rar, taxonomic.rank="Phylum")
Chloro_all.reads=sample_sums(Chloro_SILVA_MERDS_rar)
all.reads=sample_sums(SILVA_MERDS_rar)
Chloro_soil_all.reads=cbind(Chloro_all.reads,all.reads)
colnames(Chloro_soil_all.reads)<-c("Chloro.reads","total.reads")
Chloro_soil_all.reads=merge(Chloro_soil_all.reads, SILVA_MERDS_rar_map, by ="row.names")
head(Chloro_soil_all.reads)
Chloro_soil_all.reads=mutate(Chloro_soil_all.reads, Chloro.prop=Chloro.reads/total.reads)

Chloro_soil_all.reads_live_NS=subset(Chloro_soil_all.reads,life_stage!="Start"&soil_root!="S.B")
nrow(Chloro_soil_all.reads_live_NS)
#32
#Prop of Chloroflexi
Chloro_model_live_NS= lm((Chloro.prop)~soil_root*life_stage*precip+as.factor(block), data= Chloro_soil_all.reads_live_NS)
qqPlot(resid(Chloro_model_live_NS))
hist(resid(Chloro_model_live_NS))
boxCox(Chloro_model_live_NS)
shapiro.test(resid(Chloro_model_live_NS))
#0.5024
Anova(Chloro_model_live_NS, type=3)
#soil_root:life_stage:precip 0.0004782  1   5.1282   0.03425 *  

emmeans(Chloro_model_live_NS, pairwise~soil_root)


emmeans(Chloro_model_live_NS, pairwise~soil_root*life_stage*precip)


Chloro_model_live_NS_no_block= lm((Chloro.prop)~soil_root*life_stage*precip, data= Chloro_soil_all.reads_live_NS)
qqPlot(resid(Chloro_model_live_NS_no_block))
hist(resid(Chloro_model_live_NS_no_block))
boxCox(Chloro_model_live_NS_no_block)
shapiro.test(resid(Chloro_model_live_NS_no_block))
#0.3977
Anova(Chloro_model_live_NS_no_block, type=3)
#soil_root                   0.0004285  1   4.8615   0.03729 *  
#soil_root:life_stage:precip 0.0005160  1   5.8542   0.02348 *  

emmeans(Chloro_model_live_NS_no_block, pairwise~soil_root)


emmeans(Chloro_model_live_NS_no_block, pairwise~soil_root*life_stage*precip,adjust="fdr")

AIC(Chloro_model_live_NS,Chloro_model_live_NS_no_block)



#####Transplant Actinobacteria Analyses######
get_taxa_unique(SILVA_MERDS_rar, taxonomic.rank="Phylum")
SILVA_MERDS_rar_map=sample_data(SILVA_MERDS_rar)
#p:Actinobacteria

Actino_SILVA_MERDS_rar <- subset_taxa(SILVA_MERDS_rar, Phylum=="p:Actinobacteria")
get_taxa_unique(Actino_SILVA_MERDS_rar, taxonomic.rank="Phylum")
Actino_all.reads=sample_sums(Actino_SILVA_MERDS_rar)
all.reads=sample_sums(SILVA_MERDS_rar)
Actino_soil_all.reads=cbind(Actino_all.reads,all.reads)
colnames(Actino_soil_all.reads)<-c("Actino.reads","total.reads")
Actino_soil_all.reads=merge(Actino_soil_all.reads, SILVA_MERDS_rar_map, by ="row.names")
head(Actino_soil_all.reads)
Actino_soil_all.reads=mutate(Actino_soil_all.reads, Actino.prop=Actino.reads/total.reads)

Actino_soil_all.reads_trt_trans=subset(Actino_soil_all.reads, life_stage=="G"&life_stage!="Start")
nrow(sample_data(Actino_soil_all.reads_trt_trans))
#24

#Prop of Actinobacteria
Actin_model_trans= lm((Actino.prop)~soil_root*precip+as.factor(block), data= Actino_soil_all.reads_trt_trans)
qqPlot(resid(Actin_model_trans))
hist(resid(Actin_model_trans))
boxCox(Actin_model_trans)
shapiro.test(resid(Actin_model_trans))
#0.7867
Anova(Actin_model_trans, type=3)
#soil_root        0.041116  2  37.4430 1.472e-06 ***
#precip           0.006325  1  11.5192  0.004007 ** 
#soil_root:precip 0.005279  2   4.8078  0.024357 *

emmeans(Actin_model_trans, pairwise~soil_root)


emmeans(Actin_model_trans, pairwise~soil_root|precip)

#####Presence Actinobacteria#####
Actino_soil_all.reads_pres_trans=subset(Actino_soil_all.reads_trt_trans, root_association =="B")
nrow(Chloro_soil_all.reads_pres_trans)

Actino_model_trans_pres= lm((Actino.prop)~soil_root*precip+as.factor(block), data= Actino_soil_all.reads_pres_trans)
qqPlot(resid(Actino_model_trans_pres))
hist(resid(Actino_model_trans_pres))
shapiro.test(resid(Actino_model_trans_pres))
#0.06574
Anova(Actino_model_trans_pres, type=3)
#soil_root        0.031073  1 110.4824 2.360e-06 ***
#precip           0.007102  1  25.2528  0.000714 ***
#as.factor(block) 0.002543  3   3.0142  0.086858 .  
#soil_root:precip 0.004157  1  14.7807  0.003939 **   


emmeans(Actino_model_trans_pres, pairwise~soil_root|precip)
emmeans(Actino_model_trans_pres, pairwise~soil_root)
#
emmeans(Actino_model_trans_pres, pairwise~soil_root*precip,adjust="none")


Actino_model_trans_pres_no_block= lm(Actino.prop~soil_root*precip, data= Actino_soil_all.reads_pres_trans)
qqPlot(resid(Actino_model_trans_pres_no_block))
hist(resid(Actino_model_trans_pres_no_block))
shapiro.test(resid(Actino_model_trans_pres_no_block))
#0.1999

Anova(Actino_model_trans_pres_no_block, type=3)
#soil_root        0.031073  1  73.4812 1.841e-06 ***
#precip           0.007102  1  16.7955  0.001477 ** 
#soil_root:precip 0.004157  1   9.8305  0.008606 ** 

emmeans(Actino_model_trans_pres_no_block, pairwise~soil_root*precip,adjust="none")


AIC(Actino_model_trans_pres,Actino_model_trans_pres_no_block)


#####Origin Actinobacteria#####
Actino_soil_all.reads_orig_trans=subset(Actino_soil_all.reads_trt_trans, soil_status =="L")
nrow(Actino_soil_all.reads_orig_trans)

Actino_model_trans_orig= lm((Actino.prop)~soil_root*precip+as.factor(block), data= Actino_soil_all.reads_orig_trans)
qqPlot(resid(Actino_model_trans_orig))
hist(resid(Actino_model_trans_orig))
shapiro.test(resid(Actino_model_trans_orig))
#0.4982
Anova(Actino_model_trans_orig, type=3)
#precip           0.007656  1   9.1173   0.01449 *  
#soil_root:precip 0.003752  1   4.4675   0.06369 . 


emmeans(Actino_model_trans_orig, pairwise~soil_root|precip)
emmeans(Actino_model_trans_orig, pairwise~soil_root)
#



Actino_model_trans_orig_no_block= lm((Actino.prop)~soil_root*precip, data= Actino_soil_all.reads_orig_trans)
qqPlot(resid(Actino_model_trans_orig_no_block))
hist(resid(Actino_model_trans_orig_no_block))
shapiro.test(resid(Actino_model_trans_orig_no_block))
#0.503

Anova(Actino_model_trans_orig_no_block, type=3)
#precip           0.007656  1   8.5921   0.01258 *  
#soil_root:precip 0.003752  1   4.2101   0.06267 .

emmeans(Actino_model_trans_orig_no_block, pairwise~soil_root*precip,adjust="none")


AIC(Actino_model_trans_orig,Actino_model_trans_orig_no_block)







Actino_soil_all.reads_trt_trans %>% group_by(soil_root,precip) %>% summarise_at("Actino.prop", funs(n(),mean,sd,se=sd(.)/sqrt(n())))


Actin_MERDS_rar_trt_trans_precip_soil_g=Actino_soil_all.reads_trt_trans %>% group_by(soil_root,precip)
Actino_prop_precip_trans=summarise_at(Actin_MERDS_rar_trt_trans_precip_soil_g, 
                                      "Actino.prop", funs(n(),mean,sd,se=sd(.)/sqrt(n())))
treatment_order=c("S.B","L.B","L.R")
(trans_obs_Actin_p=ggplot(Actino_prop_precip_trans, aes(precip,mean,fill=factor(soil_root,levels = treatment_order)))+
    geom_bar(position="dodge",stat="identity", color="black",aes(fill=factor(soil_root,levels = treatment_order)))+
    geom_errorbar(aes(ymin = mean-se, ymax= mean+se),position=position_dodge(width=0.9), width=0.2, size=1)+
    scale_fill_manual(values = c( "white","lightgray", "darkgrey"),labels=c("Sterile","Bulk","Rhizo"), limits=treatment_order, name=NULL)+
    scale_x_discrete(labels=c("Ambient","Drought"), name=NULL)+ylab("Actinobacteria reads (proportion)")+
    theme_bw()+
    theme(axis.title = element_text(size = 23), axis.text = element_text(size = 23),
          legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

#actino_seedling_prop

#####Transplant Firmicutes Analyses######
get_taxa_unique(SILVA_MERDS_rar, taxonomic.rank="Phylum")
SILVA_MERDS_rar_map=sample_data(SILVA_MERDS_rar)
#p:Firmicutes

Firm_SILVA_MERDS_rar <- subset_taxa(SILVA_MERDS_rar, Phylum=="p:Firmicutes")
get_taxa_unique(Firm_SILVA_MERDS_rar, taxonomic.rank="Phylum")
Firm_all.reads=sample_sums(Firm_SILVA_MERDS_rar)
all.reads=sample_sums(SILVA_MERDS_rar)
Firm_soil_all.reads=cbind(Firm_all.reads,all.reads)
colnames(Firm_soil_all.reads)<-c("Firm.reads","total.reads")
Firm_soil_all.reads=merge(Firm_soil_all.reads, SILVA_MERDS_rar_map, by ="row.names")
head(Firm_soil_all.reads)
Firm_soil_all.reads=mutate(Firm_soil_all.reads, Firm.prop=Firm.reads/total.reads)

Firm_soil_all.reads_trt_trans=subset(Firm_soil_all.reads, life_stage=="G"&life_stage!="Start")
nrow(sample_data(Firm_soil_all.reads_trt_trans))
#24

#Prop of Firmicutes
Firm_model_trans= lm((Firm.prop)~soil_root*precip+as.factor(block), data= Firm_soil_all.reads_trt_trans)
qqPlot(resid(Firm_model_trans))
hist(resid(Firm_model_trans))
boxCox(Firm_model_trans)
shapiro.test(resid(Firm_model_trans))
#0.8717
Anova(Firm_model_trans, type=3)
#soil_root         87546  2  4.4342   0.03069 *  

emmeans(Firm_model_trans, pairwise~soil_root)


emmeans(Firm_model_trans, pairwise~soil_root|precip)


Firm_model_trans_no_block= lm((Firm.prop)~soil_root*precip, data= Firm_soil_all.reads_trt_trans)
qqPlot(resid(Firm_model_trans_no_block))
hist(resid(Firm_model_trans_no_block))
boxCox(Firm_model_trans_no_block)
shapiro.test(resid(Firm_model_trans_no_block))
#0.04238
Anova(Firm_model_trans_no_block, type=3)
#soil_root         87546  2  3.9268   0.03844 *  

emmeans(Firm_model_trans_no_block, pairwise~soil_root)


emmeans(Firm_model_trans_no_block, pairwise~soil_root|precip)



#####Presence Firmicutes#####
Firm_soil_all.reads_pres_trans=subset(Firm_soil_all.reads_trt_trans, root_association =="B")
nrow(Firm_soil_all.reads_pres_trans)

Firm_model_trans_pres= lm((Firm.prop)~soil_root*precip+as.factor(block), data= Firm_soil_all.reads_pres_trans)
qqPlot(resid(Firm_model_trans_pres))
hist(resid(Firm_model_trans_pres))
shapiro.test(resid(Firm_model_trans_pres))
#0.3029
Anova(Firm_model_trans_pres, type=3)
#soil_root         81083  1  8.3748   0.01777 *  


emmeans(Firm_model_trans_pres, pairwise~soil_root|precip)
emmeans(Firm_model_trans_pres, pairwise~soil_root)
#



Firm_model_trans_pres_no_block= lm(Firm.prop~soil_root*precip, data= Firm_soil_all.reads_pres_trans)
qqPlot(resid(Firm_model_trans_pres_no_block))
hist(resid(Firm_model_trans_pres_no_block))
shapiro.test(resid(Firm_model_trans_pres_no_block))
#0.09909

Anova(Firm_model_trans_pres_no_block, type=3)
#soil_root        0.0008108  1  7.4154    0.0185 * 

emmeans(Firm_model_trans_pres_no_block, pairwise~soil_root*precip,adjust="none")


AIC(Firm_model_trans_pres,Firm_model_trans_pres_no_block)


#####Origin Firmicutes#####
Firm_soil_all.reads_orig_trans=subset(Firm_soil_all.reads_trt_trans, soil_status =="L")
nrow(Firm_soil_all.reads_orig_trans)

Firm_model_trans_orig= lm((Firm.prop)~soil_root*precip+as.factor(block), data= Firm_soil_all.reads_orig_trans)
qqPlot(resid(Firm_model_trans_orig))
hist(resid(Firm_model_trans_orig))
shapiro.test(resid(Firm_model_trans_orig))
#0.362
Anova(Firm_model_trans_orig, type=3)
#nada 


emmeans(Firm_model_trans_orig, pairwise~soil_root|precip)
emmeans(Firm_model_trans_orig, pairwise~soil_root)
#



Firm_model_trans_orig_no_block= lm((Firm.prop)~soil_root*precip, data= Firm_soil_all.reads_orig_trans)
qqPlot(resid(Firm_model_trans_orig_no_block))
hist(resid(Firm_model_trans_orig_no_block))
shapiro.test(resid(Firm_model_trans_orig_no_block))
#0.03271

Anova(Firm_model_trans_orig_no_block, type=3)
#nada   

emmeans(Firm_model_trans_orig_no_block, pairwise~soil_root*precip,adjust="none")


AIC(Firm_model_trans_orig,Firm_model_trans_orig_no_block)



Firm_soil_all.reads_trt_trans %>% group_by(soil_root,precip) %>% summarise_at("Firm.prop", funs(n(),mean,sd,se=sd(.)/sqrt(n())))


Firm_MERDS_rar_trt_trans_precip_soil_g=Firm_soil_all.reads_trt_trans %>% group_by(soil_root,precip)
Firm_prop_precip_trans=summarise_at(Firm_MERDS_rar_trt_trans_precip_soil_g, 
                                    "Firm.prop", funs(n(),mean,sd,se=sd(.)/sqrt(n())))
treatment_order=c("S.B","L.B","L.R")
(trans_obs_Firm_p=ggplot(Firm_prop_precip_trans, aes(x=precip,y=mean, ymin = mean-se, ymax= mean+se,fill=factor(soil_root,levels = treatment_order)))+
    geom_bar(position="dodge",stat="identity", color="black",aes(fill=factor(soil_root,levels = treatment_order)))+
    geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+
    scale_fill_manual(values = c( "white","lightgray", "darkgrey"),labels=c("Sterile","Bulk","Rhizo"), limits=treatment_order, name=NULL)+
    scale_x_discrete(labels=c("Ambient","Drought"), name=NULL)+ylab("Firmicutes reads (proportion)")+
    theme_bw()+
    theme(axis.title = element_text(size = 23), axis.text = element_text(size = 23),
          legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank()))
#9x8




#####Transplant Chloroflexi Analyses######
get_taxa_unique(SILVA_MERDS_rar, taxonomic.rank="Phylum")
SILVA_MERDS_rar_map=sample_data(SILVA_MERDS_rar)
#p:Chloroflexi

Chloro_SILVA_MERDS_rar <- subset_taxa(SILVA_MERDS_rar, Phylum=="p:Chloroflexi")
get_taxa_unique(Chloro_SILVA_MERDS_rar, taxonomic.rank="Phylum")
Chloro_all.reads=sample_sums(Chloro_SILVA_MERDS_rar)
all.reads=sample_sums(SILVA_MERDS_rar)
Chloro_soil_all.reads=cbind(Chloro_all.reads,all.reads)
colnames(Chloro_soil_all.reads)<-c("Chloro.reads","total.reads")
Chloro_soil_all.reads=merge(Chloro_soil_all.reads, SILVA_MERDS_rar_map, by ="row.names")
head(Chloro_soil_all.reads)
Chloro_soil_all.reads=mutate(Chloro_soil_all.reads, Chloro.prop=Chloro.reads/total.reads)

Chloro_soil_all.reads_trt_trans=subset(Chloro_soil_all.reads, life_stage=="G"&life_stage!="Start")
nrow(sample_data(Chloro_soil_all.reads_trt_trans))
#24

#Prop of Chloroflexi
Chloro_model_trans= lm((Chloro.prop)~soil_root*precip+as.factor(block), data= Chloro_soil_all.reads_trt_trans)
qqPlot(resid(Chloro_model_trans))
hist(resid(Chloro_model_trans))
boxCox(Chloro_model_trans)
shapiro.test(resid(Chloro_model_trans))
#0.2371
Anova(Chloro_model_trans, type=3)
#soil_root        0.0047352  2  25.5809 1.466e-05 ***
#soil_root:precip 0.0006476  2   3.4986   0.05661 . 

emmeans(Chloro_model_trans, pairwise~soil_root)


emmeans(Chloro_model_trans, pairwise~soil_root|precip)


Chloro_model_trans_no_block= lm((Chloro.prop)~soil_root*precip, data= Chloro_soil_all.reads_trt_trans)
qqPlot(resid(Chloro_model_trans_no_block))
hist(resid(Chloro_model_trans_no_block))
boxCox(Chloro_model_trans_no_block)
shapiro.test(resid(Chloro_model_trans_no_block))
#0.02399
Anova(Chloro_model_trans_no_block, type=3)
#soil_root        0.0047352  2  28.6065 2.575e-06 ***
#soil_root:precip 0.0006476  2   3.9124   0.03883 *  

emmeans(Chloro_model_trans_no_block, pairwise~soil_root)


emmeans(Chloro_model_trans_no_block, pairwise~soil_root|precip)



#####Presence Chloroflexi#####
Chloro_soil_all.reads_pres_trans=subset(Chloro_soil_all.reads_trt_trans, root_association =="B")
nrow(Chloro_soil_all.reads_pres_trans)

Chloro_model_trans_pres= lm((Chloro.prop)~soil_root*precip+as.factor(block), data= Chloro_soil_all.reads_pres_trans)
qqPlot(resid(Chloro_model_trans_pres))
hist(resid(Chloro_model_trans_pres))
shapiro.test(resid(Chloro_model_trans_pres))
#0.2041
Anova(Chloro_model_trans_pres, type=3)
#soil_root        0.0032433  1 163.8675 4.432e-07 ***
#precip           0.0001010  1   5.1031   0.05026 .  
#soil_root:precip 0.0000941  1   4.7539   0.05714 .  


emmeans(Chloro_model_trans_pres, pairwise~soil_root|precip)
emmeans(Chloro_model_trans_pres, pairwise~soil_root)
#



Chloro_model_trans_pres_no_block= lm((Chloro.prop)~soil_root*precip, data= Chloro_soil_all.reads_pres_trans)
qqPlot(resid(Chloro_model_trans_pres_no_block))
hist(resid(Chloro_model_trans_pres_no_block))
shapiro.test(resid(Chloro_model_trans_pres_no_block))
#0.06719

Anova(Chloro_model_trans_pres_no_block, type=3)
#soil_root        0.0032433  1 159.3924 2.741e-08 ***
#precip           0.0001010  1   4.9638   0.04578 *  
#soil_root:precip 0.0000941  1   4.6241   0.05261 .  

emmeans(Chloro_model_trans_pres_no_block, pairwise~soil_root*precip,adjust="none")


AIC(Chloro_model_trans_pres,Chloro_model_trans_pres_no_block)


#####Origin Chloroflexi#####
Chloro_soil_all.reads_orig_trans=subset(Chloro_soil_all.reads_trt_trans, soil_status =="L")
nrow(Chloro_soil_all.reads_orig_trans)

Chloro_model_trans_orig= lm((Chloro.prop)~soil_root*precip+as.factor(block), data= Chloro_soil_all.reads_orig_trans)
qqPlot(resid(Chloro_model_trans_orig))
hist(resid(Chloro_model_trans_orig))
shapiro.test(resid(Chloro_model_trans_orig))
#0.5008
Anova(Chloro_model_trans_orig, type=3)
#soil_root:precip 0.0006363  1  4.2793   0.06851 .   


emmeans(Chloro_model_trans_orig, pairwise~soil_root|precip)
emmeans(Chloro_model_trans_orig, pairwise~soil_root)
#



Chloro_model_trans_orig_no_block= lm((Chloro.prop)~soil_root*precip, data= Chloro_soil_all.reads_orig_trans)
qqPlot(resid(Chloro_model_trans_orig_no_block))
hist(resid(Chloro_model_trans_orig_no_block))
shapiro.test(resid(Chloro_model_trans_orig_no_block))
#0.3459

Anova(Chloro_model_trans_orig_no_block, type=3)
#soil_root:precip 0.0006363  1   5.1282   0.04286 *   

emmeans(Chloro_model_trans_orig_no_block, pairwise~soil_root*precip,adjust="none")


AIC(Chloro_model_trans_orig,Chloro_model_trans_orig_no_block)



Chloro_soil_all.reads_trt_trans %>% group_by(soil_root,precip) %>% summarise_at("Chloro.prop", funs(n(),mean,sd,se=sd(.)/sqrt(n())))


Chloro_MERDS_rar_trt_trans_precip_soil_g=Chloro_soil_all.reads_trt_trans %>% group_by(soil_root,precip)
Chloro_prop_precip_trans=summarise_at(Chloro_MERDS_rar_trt_trans_precip_soil_g, 
                                      "Chloro.prop", funs(n(),mean,sd,se=sd(.)/sqrt(n())))
treatment_order=c("S.B","L.B","L.R")
(trans_obs_Chloro_p=ggplot(Chloro_prop_precip_trans, aes(x=precip,y=mean, ymin = mean-se, ymax= mean+se,fill=factor(soil_root,levels = treatment_order)))+
    geom_bar(position="dodge",stat="identity", color="black",aes(fill=factor(soil_root,levels = treatment_order)))+
    geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+
    scale_fill_manual(values = c( "white","lightgray", "darkgrey"),labels=c("Sterile","Bulk","Rhizo"), limits=treatment_order, name=NULL)+
    scale_x_discrete(labels=c("Ambient","Drought"), name=NULL)+ylab("Chloroflexi reads (proportion)")+
    theme_bw()+
    theme(axis.title = element_text(size = 23), axis.text = element_text(size = 23),
          legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank()))
#9x8


#####Seed Actinobacteria Analyses######
get_taxa_unique(SILVA_MERDS_rar, taxonomic.rank="Phylum")
SILVA_MERDS_rar_map=sample_data(SILVA_MERDS_rar)
#p:Actinobacteria

Actino_SILVA_MERDS_rar <- subset_taxa(SILVA_MERDS_rar, Phylum=="p:Actinobacteria")
get_taxa_unique(Actino_SILVA_MERDS_rar, taxonomic.rank="Phylum")
Actino_all.reads=sample_sums(Actino_SILVA_MERDS_rar)
all.reads=sample_sums(SILVA_MERDS_rar)
Actino_soil_all.reads=cbind(Actino_all.reads,all.reads)
colnames(Actino_soil_all.reads)<-c("Actino.reads","total.reads")
Actino_soil_all.reads=merge(Actino_soil_all.reads, SILVA_MERDS_rar_map, by ="row.names")
head(Actino_soil_all.reads)
Actino_soil_all.reads=mutate(Actino_soil_all.reads, Actino.prop=Actino.reads/total.reads)

Actino_soil_all.reads_trt_seed=subset(Actino_soil_all.reads, life_stage=="S"&life_stage!="Start")
nrow(sample_data(Actino_soil_all.reads_trt_seed))
#24

#Prop of Actinobacteria
Actin_model_seed= lm((Actino.prop)~soil_root*precip+as.factor(block), data= Actino_soil_all.reads_trt_seed)
qqPlot(resid(Actin_model_seed))
hist(resid(Actin_model_seed))
boxCox(Actin_model_seed)
shapiro.test(resid(Actin_model_seed))
#0.4285
Anova(Actin_model_seed, type=3)
#soil_root        0.017457  2  11.4078 0.0011499 ** 
#precip           0.015790  1  20.6358 0.0004601 ***

emmeans(Actin_model_seed, pairwise~soil_root)


emmeans(Actin_model_seed, pairwise~soil_root|precip)

#####Presence Actinobacteria#####
Actino_soil_all.reads_pres_seed=subset(Actino_soil_all.reads_trt_seed, root_association =="B")
nrow(Actino_soil_all.reads_pres_seed)

Actino_model_seed_pres= lm((Actino.prop)~soil_root*precip+as.factor(block), data= Actino_soil_all.reads_pres_seed)
qqPlot(resid(Actino_model_seed_pres))
hist(resid(Actino_model_seed_pres))
shapiro.test(resid(Actino_model_seed_pres))
#0.8234
Anova(Actino_model_seed_pres, type=3)
#precip           0.010103  1  12.8508  0.007138 ** 


emmeans(Actino_model_seed_pres, pairwise~soil_root|precip)
emmeans(Actino_model_seed_pres, pairwise~soil_root)
#
emmeans(Actino_model_seed_pres, pairwise~soil_root*precip,adjust="none")


Actino_model_seed_pres_no_block= lm(Actino.prop~soil_root*precip, data= Actino_soil_all.reads_pres_seed)
qqPlot(resid(Actino_model_seed_pres_no_block))
hist(resid(Actino_model_seed_pres_no_block))
shapiro.test(resid(Actino_model_seed_pres_no_block))
#0.3869

Anova(Actino_model_seed_pres_no_block, type=3)
#soil_root        0.002289  1   3.2460  0.099041 .  
#precip           0.009860  1  13.9801  0.003272 ** 

emmeans(Actino_model_seed_pres_no_block, pairwise~soil_root*precip,adjust="none")


AIC(Actino_model_seed_pres,Actino_model_seed_pres_no_block)


#####Origin Actinobacteria#####
Actino_soil_all.reads_orig_seed=subset(Actino_soil_all.reads_trt_seed, soil_status =="L")
nrow(Actino_soil_all.reads_orig_seed)

Actino_model_seed_orig= lm((Actino.prop)~soil_root*precip+as.factor(block), data= Actino_soil_all.reads_orig_seed)
qqPlot(resid(Actino_model_seed_orig))
hist(resid(Actino_model_seed_orig))
shapiro.test(resid(Actino_model_seed_orig))
# 0.6791
Anova(Actino_model_seed_orig, type=3)
#soil_root        0.008624  1  11.1390  0.008695 ** 
#precip           0.007252  1   9.3665  0.013561 *  


emmeans(Actino_model_seed_orig, pairwise~soil_root|precip)
emmeans(Actino_model_seed_orig, pairwise~soil_root)
#



Actino_model_seed_orig_no_block= lm((Actino.prop)~soil_root*precip, data= Actino_soil_all.reads_orig_seed)
qqPlot(resid(Actino_model_seed_orig_no_block))
hist(resid(Actino_model_seed_orig_no_block))
shapiro.test(resid(Actino_model_seed_orig_no_block))
#0.8191

Anova(Actino_model_seed_orig_no_block, type=3)
#soil_root        0.010257  1  14.3287  0.002598 ** 
#precip           0.008728  1  12.1935  0.004448 **

emmeans(Actino_model_seed_orig_no_block, pairwise~soil_root*precip,adjust="none")


AIC(Actino_model_seed_orig,Actino_model_seed_orig_no_block)







Actino_soil_all.reads_trt_seed %>% group_by(soil_root,precip) %>% summarise_at("Actino.prop", funs(n(),mean,sd,se=sd(.)/sqrt(n())))


Actin_MERDS_rar_trt_seed_precip_soil_g=Actino_soil_all.reads_trt_seed %>% group_by(soil_root,precip)
Actino_prop_precip_seed=summarise_at(Actin_MERDS_rar_trt_seed_precip_soil_g, 
                                     "Actino.prop", funs(n(),mean,sd,se=sd(.)/sqrt(n())))
treatment_order=c("S.B","L.B","L.R")
(seed_obs_Actin_p=ggplot(Actino_prop_precip_seed, aes(x=precip,y=mean, ymin = mean-se, ymax= mean+se,fill=factor(soil_root,levels = treatment_order)))+
    geom_bar(position="dodge",stat="identity", color="black",aes(fill=factor(soil_root,levels = treatment_order)))+
    geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+
    scale_fill_manual(values = c( "white","lightgray", "darkgrey"),labels=c("Sterile","Bulk","Rhizo"), limits=treatment_order, name=NULL)+
    scale_x_discrete(labels=c("Ambient","Drought"), name=NULL)+ylab("Actinobacteria reads (proportion)")+
    theme_bw()+
    theme(axis.title = element_text(size = 23), axis.text = element_text(size = 23),
          legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

#actino_seedling_prop

#####Seed Firmicutes Analyses######
get_taxa_unique(SILVA_MERDS_rar, taxonomic.rank="Phylum")
SILVA_MERDS_rar_map=sample_data(SILVA_MERDS_rar)
#p:Firmicutes

Firm_SILVA_MERDS_rar <- subset_taxa(SILVA_MERDS_rar, Phylum=="p:Firmicutes")
get_taxa_unique(Firm_SILVA_MERDS_rar, taxonomic.rank="Phylum")
Firm_all.reads=sample_sums(Firm_SILVA_MERDS_rar)
all.reads=sample_sums(SILVA_MERDS_rar)
Firm_soil_all.reads=cbind(Firm_all.reads,all.reads)
colnames(Firm_soil_all.reads)<-c("Firm.reads","total.reads")
Firm_soil_all.reads=merge(Firm_soil_all.reads, SILVA_MERDS_rar_map, by ="row.names")
head(Firm_soil_all.reads)
Firm_soil_all.reads=mutate(Firm_soil_all.reads, Firm.prop=Firm.reads/total.reads)

Firm_soil_all.reads_trt_seed=subset(Firm_soil_all.reads, life_stage=="S"&life_stage!="Start")
nrow(sample_data(Firm_soil_all.reads_trt_seed))
#23

#Prop of Firmicutes
Firm_model_seed= lm(logit(Firm.prop)~soil_root*precip+as.factor(block), data= Firm_soil_all.reads_trt_seed)
qqPlot(resid(Firm_model_seed))
hist(resid(Firm_model_seed))
boxCox(Firm_model_seed)
shapiro.test(resid(Firm_model_seed))
#0.6711
Anova(Firm_model_seed, type=3)
#soil_root          8.859  2   16.0143 0.0002408 ***
#precip             3.569  1   12.9042 0.0029435 ** 

emmeans(Firm_model_seed, pairwise~soil_root)


emmeans(Firm_model_seed, pairwise~soil_root|precip)


Firm_model_seed_no_block= lm(logit(Firm.prop)~soil_root*precip, data= Firm_soil_all.reads_trt_seed)
qqPlot(resid(Firm_model_seed_no_block))
hist(resid(Firm_model_seed_no_block))
boxCox(Firm_model_seed_no_block)
shapiro.test(resid(Firm_model_seed_no_block))
#0.01997
Anova(Firm_model_seed_no_block, type=3)
#soil_root          9.421  2   18.0351 6.275e-05 ***
#precip             3.674  1   14.0690  0.001592 **  

emmeans(Firm_model_seed_no_block, pairwise~soil_root)


emmeans(Firm_model_seed_no_block, pairwise~soil_root|precip)



#####Presence Firmicutes#####
Firm_soil_all.reads_pres_seed=subset(Firm_soil_all.reads_trt_seed, root_association =="B")
nrow(Firm_soil_all.reads_pres_seed)

Firm_model_seed_pres= lm((Firm.prop)~soil_root*precip+as.factor(block), data= Firm_soil_all.reads_pres_seed)
qqPlot(resid(Firm_model_seed_pres))
hist(resid(Firm_model_seed_pres))
shapiro.test(resid(Firm_model_seed_pres))
#0.321
Anova(Firm_model_seed_pres, type=3)
#soil_root        0.0112213  1 10.3570 0.01227 *


emmeans(Firm_model_seed_pres, pairwise~soil_root|precip)
emmeans(Firm_model_seed_pres, pairwise~soil_root)
#



Firm_model_seed_pres_no_block= lm(Firm.prop~soil_root*precip, data= Firm_soil_all.reads_pres_seed)
qqPlot(resid(Firm_model_seed_pres_no_block))
hist(resid(Firm_model_seed_pres_no_block))
shapiro.test(resid(Firm_model_seed_pres_no_block))
#0.1103

Anova(Firm_model_seed_pres_no_block, type=3)
#soil_root        0.0127406  1 13.6009 0.0035760 ** 
#precip           0.0037143  1  3.9651 0.0718735 .  

emmeans(Firm_model_seed_pres_no_block, pairwise~soil_root*precip,adjust="none")


AIC(Firm_model_seed_pres,Firm_model_seed_pres_no_block)


#####Origin Firmicutes#####
Firm_soil_all.reads_orig_seed=subset(Firm_soil_all.reads_trt_seed, soil_status =="L")
nrow(Firm_soil_all.reads_orig_seed)

Firm_model_seed_orig= lm(logit(Firm.prop)~soil_root*precip+as.factor(block), data= Firm_soil_all.reads_orig_seed)
qqPlot(resid(Firm_model_seed_orig))
hist(resid(Firm_model_seed_orig))
shapiro.test(resid(Firm_model_seed_orig))
#0.9278
Anova(Firm_model_seed_orig, type=3)
#soil_root          0.683  1    3.8970  0.079817 .  
#precip             1.902  1   10.8548  0.009309 **    


emmeans(Firm_model_seed_orig, pairwise~soil_root|precip)
emmeans(Firm_model_seed_orig, pairwise~soil_root)
#
emmeans(Firm_model_seed_orig, pairwise~soil_root*precip,adjust="none")


Firm_model_seed_orig_no_block= lm(logit(Firm.prop)~soil_root*precip, data= Firm_soil_all.reads_orig_seed)
qqPlot(resid(Firm_model_seed_orig_no_block))
hist(resid(Firm_model_seed_orig_no_block))
shapiro.test(resid(Firm_model_seed_orig_no_block))
#0.01279

Anova(Firm_model_seed_orig_no_block, type=3)
#soil_root          0.992  1    4.4686  0.056140 .  
#precip             2.429  1   10.9491  0.006237 ** 

emmeans(Firm_model_seed_orig_no_block, pairwise~soil_root*precip,adjust="none")


AIC(Firm_model_seed_orig,Firm_model_seed_orig_no_block)



Firm_soil_all.reads_trt_seed %>% group_by(soil_root,precip) %>% summarise_at("Firm.prop", funs(n(),mean,sd,se=sd(.)/sqrt(n())))


Firm_MERDS_rar_trt_seed_precip_soil_g=Firm_soil_all.reads_trt_seed %>% group_by(soil_root,precip)
Firm_prop_precip_seed=summarise_at(Firm_MERDS_rar_trt_seed_precip_soil_g, 
                                   "Firm.prop", funs(n(),mean,sd,se=sd(.)/sqrt(n())))
treatment_order=c("S.B","L.B","L.R")
(seed_obs_Firm_p=ggplot(Firm_prop_precip_seed, aes(x=precip,y=mean, ymin = mean-se, ymax= mean+se,fill=factor(soil_root,levels = treatment_order)))+
    geom_bar(position="dodge",stat="identity", color="black",aes(fill=factor(soil_root,levels = treatment_order)))+
    geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+
    scale_fill_manual(values = c( "white","lightgray", "darkgrey"),labels=c("Sterile","Bulk","Rhizo"), limits=treatment_order, name=NULL)+
    scale_x_discrete(labels=c("Ambient","Drought"), name=NULL)+ylab("Firmicutes reads (proportion)")+
    theme_bw()+
    theme(axis.title = element_text(size = 23), axis.text = element_text(size = 23),
          legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank()))
#9x8




#####Seed Chloroflexi Analyses######
get_taxa_unique(SILVA_MERDS_rar, taxonomic.rank="Phylum")
SILVA_MERDS_rar_map=sample_data(SILVA_MERDS_rar)
#p:Chloroflexi

Chloro_SILVA_MERDS_rar <- subset_taxa(SILVA_MERDS_rar, Phylum=="p:Chloroflexi")
get_taxa_unique(Chloro_SILVA_MERDS_rar, taxonomic.rank="Phylum")
Chloro_all.reads=sample_sums(Chloro_SILVA_MERDS_rar)
all.reads=sample_sums(SILVA_MERDS_rar)
Chloro_soil_all.reads=cbind(Chloro_all.reads,all.reads)
colnames(Chloro_soil_all.reads)<-c("Chloro.reads","total.reads")
Chloro_soil_all.reads=merge(Chloro_soil_all.reads, SILVA_MERDS_rar_map, by ="row.names")
head(Chloro_soil_all.reads)
Chloro_soil_all.reads=mutate(Chloro_soil_all.reads, Chloro.prop=Chloro.reads/total.reads)

Chloro_soil_all.reads_trt_seed=subset(Chloro_soil_all.reads, life_stage=="S"&life_stage!="Start")
nrow(sample_data(Chloro_soil_all.reads_trt_seed))
#23

#Prop of Chloroflexi
Chloro_model_seed= lm((Chloro.prop)~soil_root*precip+as.factor(block), data= Chloro_soil_all.reads_trt_seed)
qqPlot(resid(Chloro_model_seed))
hist(resid(Chloro_model_seed))
boxCox(Chloro_model_seed)
shapiro.test(resid(Chloro_model_seed))
#0.1894
Anova(Chloro_model_seed, type=3)
#soil_root        0.0031325  2  37.2816 2.467e-06 ***
#precip           0.0001444  1   3.4382   0.08489 .  

emmeans(Chloro_model_seed, pairwise~soil_root)


emmeans(Chloro_model_seed, pairwise~soil_root|precip)


Chloro_model_seed_no_block= lm((Chloro.prop)~soil_root*precip, data= Chloro_soil_all.reads_trt_seed)
qqPlot(resid(Chloro_model_seed_no_block))
hist(resid(Chloro_model_seed_no_block))
boxCox(Chloro_model_seed_no_block)
shapiro.test(resid(Chloro_model_seed_no_block))
#0.1189
Anova(Chloro_model_seed_no_block, type=3)
#soil_root        0.0037148  2  50.3757 7.171e-08 ***
#precip           0.0001522  1   4.1287   0.05809 .  

emmeans(Chloro_model_seed_no_block, pairwise~soil_root)


emmeans(Chloro_model_seed_no_block, pairwise~soil_root|precip)



#####Presence Chloroflexi#####
Chloro_soil_all.reads_pres_seed=subset(Chloro_soil_all.reads_trt_seed, root_association =="B")
nrow(Chloro_soil_all.reads_pres_seed)

Chloro_model_seed_pres= lm((Chloro.prop)~soil_root*precip+as.factor(block), data= Chloro_soil_all.reads_pres_seed)
qqPlot(resid(Chloro_model_seed_pres))
hist(resid(Chloro_model_seed_pres))
shapiro.test(resid(Chloro_model_seed_pres))
#0.1476
Anova(Chloro_model_seed_pres, type=3)
#soil_root        0.00128615  1 71.9519 2.858e-05 *** 
#as.factor(block) 0.00018786  3  3.5032   0.06936 . 

emmeans(Chloro_model_seed_pres, pairwise~soil_root|precip)
emmeans(Chloro_model_seed_pres, pairwise~soil_root)
#
emmeans(Chloro_model_seed_pres, pairwise~soil_root*precip,adjust="none")


Chloro_model_seed_pres_no_block= lm(Chloro.prop~soil_root*precip, data= Chloro_soil_all.reads_pres_seed)
qqPlot(resid(Chloro_model_seed_pres_no_block))
hist(resid(Chloro_model_seed_pres_no_block))
shapiro.test(resid(Chloro_model_seed_pres_no_block))
#0.05492

Anova(Chloro_model_seed_pres_no_block, type=3)
#soil_root        0.00139501  1 46.3793 2.914e-05 ***

emmeans(Chloro_model_seed_pres_no_block, pairwise~soil_root*precip,adjust="none")


AIC(Chloro_model_seed_pres,Chloro_model_seed_pres_no_block)


#####Origin Chloroflexi#####
Chloro_soil_all.reads_orig_seed=subset(Chloro_soil_all.reads_trt_seed, soil_status =="L")
nrow(Chloro_soil_all.reads_orig_seed)

Chloro_model_seed_orig= lm((Chloro.prop)~soil_root*precip+as.factor(block), data= Chloro_soil_all.reads_orig_seed)
qqPlot(resid(Chloro_model_seed_orig))
hist(resid(Chloro_model_seed_orig))
shapiro.test(resid(Chloro_model_seed_orig))
#0.5532
Anova(Chloro_model_seed_orig, type=3)
#soil_root        0.0005386  1   8.3936   0.01768 *  
#precip           0.0002297  1   3.5795   0.09105 .   


emmeans(Chloro_model_seed_orig, pairwise~soil_root|precip)
emmeans(Chloro_model_seed_orig, pairwise~soil_root)
#



Chloro_model_seed_orig_no_block= lm((Chloro.prop)~soil_root*precip, data= Chloro_soil_all.reads_orig_seed)
qqPlot(resid(Chloro_model_seed_orig_no_block))
hist(resid(Chloro_model_seed_orig_no_block))
shapiro.test(resid(Chloro_model_seed_orig_no_block))
#0.6406


Anova(Chloro_model_seed_orig_no_block, type=3)
#soil_root        0.0005905  1  11.3101  0.005642 ** 
#precip           0.0002592  1   4.9648  0.045759 *  

emmeans(Chloro_model_seed_orig_no_block, pairwise~soil_root*precip,adjust="none")


AIC(Chloro_model_seed_orig,Chloro_model_seed_orig_no_block)



Chloro_soil_all.reads_trt_seed %>% group_by(soil_root,precip) %>% summarise_at("Chloro.prop", funs(n(),mean,sd,se=sd(.)/sqrt(n())))


Chloro_MERDS_rar_trt_seed_precip_soil_g=Chloro_soil_all.reads_trt_seed %>% group_by(soil_root,precip)
Chloro_prop_precip_seed=summarise_at(Chloro_MERDS_rar_trt_seed_precip_soil_g, 
                                     "Chloro.prop", funs(n(),mean,sd,se=sd(.)/sqrt(n())))
treatment_order=c("S.B","L.B","L.R")
(seed_obs_Chloro_p=ggplot(Chloro_prop_precip_seed, aes(x=precip,y=mean, ymin = mean-se, ymax= mean+se,fill=factor(soil_root,levels = treatment_order)))+
    geom_bar(position="dodge",stat="identity", color="black",aes(fill=factor(soil_root,levels = treatment_order)))+
    geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+
    scale_fill_manual(values = c( "white","lightgray", "darkgrey"),labels=c("Sterile","Bulk","Rhizo"), limits=treatment_order, name=NULL)+
    scale_x_discrete(labels=c("Ambient","Drought"), name=NULL)+ylab("Chloroflexi reads (proportion)")+
    theme_bw()+
    theme(axis.title = element_text(size = 23), axis.text = element_text(size = 23),
          legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank()))
#9x8



#####Transplant Stack bar graphs####

SILVA_MERDS_rar_trt_trans=subset_samples(SILVA_MERDS_rar, life_stage=="G"&life_stage!="Start")
nrow(sample_data(SILVA_MERDS_rar_trt_trans))
#24
SILVA_MERDS_rar_trt_trans_map=sample_data(SILVA_MERDS_rar_trt_trans)
SILVA_MERDS_rar_trt_trans_map$soil_root_precip=with(SILVA_MERDS_rar_trt_trans_map, interaction(soil_root, precip))

#Make a new Phyloseq obj with all of the associated data
SILVA_MERDS_rar_trt_trans=phyloseq(otu_table(SILVA_MERDS_rar_trt_trans),tax_table(SILVA_MERDS_rar_trt_trans),
                                   SILVA_MERDS_rar_trt_trans_map)
ntaxa(SILVA_MERDS_rar_trt_trans)
#merge OTUs by the soil and precipitation treatment
SILVA_MERDS_rar_trt_trans_fact=merge_samples(SILVA_MERDS_rar_trt_trans, "soil_root_precip")
sample_names(SILVA_MERDS_rar_trt_trans_fact)     

#combine the reads at Phylum level
get_taxa_unique(SILVA_MERDS_rar_trt_trans, taxonomic.rank="Phylum")
#44
(SILVA_MERDS_rar_trt_trans_fact.phylum<-tax_glom(SILVA_MERDS_rar_trt_trans_fact, taxrank="Phylum"))



#subset so there is only the top ten most abundant phyla
SILVA_TopPHYL_trans = names(sort(taxa_sums(SILVA_MERDS_rar_trt_trans_fact.phylum), TRUE)[1:10])
SILVA_MERDS_rar_trt_trans.T10 = prune_taxa(SILVA_TopPHYL_trans, SILVA_MERDS_rar_trt_trans_fact.phylum)
SILVA_trans_name_T10=get_taxa_unique(SILVA_MERDS_rar_trt_trans.T10, taxonomic.rank="Phylum")
SILVA_trans_name_T10_sep <- colsplit(SILVA_trans_name_T10, ":", c("letter", "Phyl_name"))


#transform the read counts to prop of total reads

SILVA_MERDS_rar_trt_trans_fact.phylum.prop=transform_sample_counts(SILVA_MERDS_rar_trt_trans_fact.phylum, function(x)x/sum(x))

taxon_positions=c("S.B.A","L.B.A", "L.R.A", "S.B.D","L.B.D" ,"L.R.D")
SILVA_MERDS_rar_trt_trans_10_prop = prune_taxa(SILVA_TopPHYL_trans, SILVA_MERDS_rar_trt_trans_fact.phylum.prop)

SILVA_MERDS_rar_trt_trans_10_prop_otu=as.data.frame(t(otu_table(SILVA_MERDS_rar_trt_trans_10_prop)))

#create an other taxa category
SILVA_taxon_sums_trans=c(sum(SILVA_MERDS_rar_trt_trans_10_prop_otu$L.B.A),sum(SILVA_MERDS_rar_trt_trans_10_prop_otu$S.B.A),
                         sum(SILVA_MERDS_rar_trt_trans_10_prop_otu$L.R.A),sum(SILVA_MERDS_rar_trt_trans_10_prop_otu$L.B.D),
                         sum(SILVA_MERDS_rar_trt_trans_10_prop_otu$S.B.D),sum(SILVA_MERDS_rar_trt_trans_10_prop_otu$L.R.D))
SILVA_tran_other_spp=c(as.numeric(1-SILVA_taxon_sums_trans))
SILVA_MERDS_rar_trt_trans_10_prop_OTU=rbind(SILVA_MERDS_rar_trt_trans_10_prop_otu,SILVA_tran_other_spp)
summary(SILVA_MERDS_rar_trt_trans_10_prop_OTU)
SILVA_MERDS_rar_trt_trans_10_prop_OTU[,"Phylum"]=c(as.character(SILVA_trans_name_T10_sep$Phyl_name),"Other Phyla")
summary(SILVA_MERDS_rar_trt_trans_10_prop_OTU)
SILVA_MERDS_rar_trt_trans_10_prop_OTU_M=melt(SILVA_MERDS_rar_trt_trans_10_prop_OTU,id="Phylum")
phyl_order=c(sort(as.character(SILVA_trans_name_T10_sep$Phyl_name)),"Other Phyla")
summary(SILVA_MERDS_rar_trt_trans_10_prop_OTU_M)



(p_bact_trans_color=ggplot(SILVA_MERDS_rar_trt_trans_10_prop_OTU_M,aes(x=variable,y=value,fill=factor(Phylum, levels=phyl_order)))+
    geom_bar(aes( fill=factor(Phylum, levels=phyl_order)), stat="identity", position="stack",color="black")+theme_bw()+
    theme(axis.text.y=element_text(size=18),axis.text.x=element_text(size=18),
          axis.title=element_text(size=20),panel.grid.major=element_blank(),legend.text = element_text(size=16), legend.title = element_text(size=20),
          panel.grid.minor=element_blank())+xlab(NULL)+ylab("Proportion")+
    scale_x_discrete(limits = taxon_positions,labels=c("Sterile\nAmbient","Bulk\nAmbient", 
                                                       "Rhizo\nAmbient",
                                                       "Sterile\nDrought",
                                                       "Bulk\nDrought", 
                                                       "Rhizo\nDrought"))+scale_fill_brewer(palette="Paired")+
    guides(fill=guide_legend(title="Phyla")))
#1000x700
#Stack bar graphs
SILVA_MERDS_rar_trt_seed_map=sample_data(SILVA_MERDS_rar_trt_seed)
SILVA_MERDS_rar_trt_seed_map$soil_root_precip=with(SILVA_MERDS_rar_trt_seed_map, interaction(soil_root, precip))

#Make a new Phyloseq obj with all of the associated data
SILVA_MERDS_rar_trt_seed=phyloseq(otu_table(SILVA_MERDS_rar_trt_seed),tax_table(SILVA_MERDS_rar_trt_seed),
                                  SILVA_MERDS_rar_trt_seed_map)
ntaxa(SILVA_MERDS_rar_trt_seed)
#8366
#merge OTUs by the soil and precipitation treatment
SILVA_MERDS_rar_trt_seed_fact=merge_samples(SILVA_MERDS_rar_trt_seed, "soil_root_precip")
sample_names(SILVA_MERDS_rar_trt_seed_fact)     

#combine the reads at Phylum level
get_taxa_unique(SILVA_MERDS_rar_trt_seed, taxonomic.rank="Phylum")
#44
(SILVA_MERDS_rar_trt_seed_fact.phylum<-tax_glom(SILVA_MERDS_rar_trt_seed_fact, taxrank="Phylum"))



#subset so there is only the top ten most abundant phyla
SILVA_TopPHYL_seed = names(sort(taxa_sums(SILVA_MERDS_rar_trt_seed_fact.phylum), TRUE)[1:10])
SILVA_MERDS_rar_trt_seed.T10 = prune_taxa(SILVA_TopPHYL_seed, SILVA_MERDS_rar_trt_seed_fact.phylum)
SILVA_seed_name_T10=get_taxa_unique(SILVA_MERDS_rar_trt_seed.T10, taxonomic.rank="Phylum")
SILVA_seed_name_T10_sep <- colsplit(SILVA_seed_name_T10, ":", c("letter", "Phyl_name"))


#seedform the read counts to prop of total reads

SILVA_MERDS_rar_trt_seed_fact.phylum.prop=transform_sample_counts(SILVA_MERDS_rar_trt_seed_fact.phylum, function(x)x/sum(x))

taxon_positions=c("S.B.A","L.B.A", "L.R.A", "S.B.D","L.B.D" ,"L.R.D")
SILVA_MERDS_rar_trt_seed_10_prop = prune_taxa(SILVA_TopPHYL_seed, SILVA_MERDS_rar_trt_seed_fact.phylum.prop)

SILVA_MERDS_rar_trt_seed_10_prop_otu=as.data.frame(t(otu_table(SILVA_MERDS_rar_trt_seed_10_prop)))

#create an other taxa category
SILVA_taxon_sums_seed=c(sum(SILVA_MERDS_rar_trt_seed_10_prop_otu$L.B.A),sum(SILVA_MERDS_rar_trt_seed_10_prop_otu$S.B.A),
                        sum(SILVA_MERDS_rar_trt_seed_10_prop_otu$L.R.A),sum(SILVA_MERDS_rar_trt_seed_10_prop_otu$L.B.D),
                        sum(SILVA_MERDS_rar_trt_seed_10_prop_otu$S.B.D),sum(SILVA_MERDS_rar_trt_seed_10_prop_otu$L.R.D))
SILVA_seed_other_spp=c(as.numeric(1-SILVA_taxon_sums_seed))
SILVA_MERDS_rar_trt_seed_10_prop_OTU=rbind(SILVA_MERDS_rar_trt_seed_10_prop_otu,SILVA_seed_other_spp)
summary(SILVA_MERDS_rar_trt_seed_10_prop_OTU)
SILVA_MERDS_rar_trt_seed_10_prop_OTU[,"Phylum"]=c(as.character(SILVA_seed_name_T10_sep$Phyl_name),"Other Phyla")
summary(SILVA_MERDS_rar_trt_seed_10_prop_OTU)
SILVA_MERDS_rar_trt_seed_10_prop_OTU_M=melt(SILVA_MERDS_rar_trt_seed_10_prop_OTU,id="Phylum")
phyl_order_seed=c(sort(as.character(SILVA_seed_name_T10_sep$Phyl_name)),"Other Phyla")
summary(SILVA_MERDS_rar_trt_seed_10_prop_OTU_M)



(p_bact_seed_color=ggplot(SILVA_MERDS_rar_trt_seed_10_prop_OTU_M,aes(x=variable,y=value,fill=factor(Phylum, levels=phyl_order_seed)))+
    geom_bar(aes( fill=factor(Phylum, levels=phyl_order_seed)), stat="identity", position="stack",color="black")+theme_bw()+
    theme(axis.text.y=element_text(size=18),axis.text.x=element_text(size=18),
          axis.title=element_text(size=20),panel.grid.major=element_blank(),legend.text = element_text(size=16), legend.title = element_text(size=20),
          panel.grid.minor=element_blank())+xlab(NULL)+ylab("Proportion")+
    scale_x_discrete(limits = taxon_positions,labels=c("Sterile\nAmbient","Bulk\nAmbient", 
                                                       "Rhizo\nAmbient",
                                                       "Sterile\nDrought",
                                                       "Bulk\nDrought", 
                                                       "Rhizo\nDrought"))+scale_fill_brewer(palette="Paired")+
    guides(fill=guide_legend(title="Phyla")))

#1000x700

#####Combined taxon stacked bargraph#####
ggarrange(p_bact_trans_color,p_bact_seed_color,ncol = 2, common.legend=F)
#20x7.38



#####Diversity####
alpha_meas = c("Observed", "Chao1", "Shannon", "InvSimpson")
SILVA_MERDS_rar_map=sample_data(SILVA_MERDS_rar)
SILVA_MERDS_rar.divfil=estimate_richness(SILVA_MERDS_rar,measures=alpha_meas)

SILVA_MERDS_rar.divfil=merge(SILVA_MERDS_rar.divfil, SILVA_MERDS_rar_map, by ="row.names")
#bact.soilE.t.divfil=mutate(bact.soilE.t.divfil, pielou=Shannon*(1/log(Observed)))
head(SILVA_MERDS_rar.divfil)
row.names(SILVA_MERDS_rar.divfil)=SILVA_MERDS_rar.divfil$Row.names
SILVA_MERDS_rar.divfil$Row.names=NULL

ggplot(SILVA_MERDS_rar.divfil, aes(x=soil_status, y=Observed))+geom_boxplot(aes(color=interaction(root_association,precip,life_stage)))+
  scale_y_continuous(name="Richness")+theme_bw()
SILVA_MERDS_rar.divfil$soil_root=with(SILVA_MERDS_rar.divfil, interaction(soil_status,root_association))
SILVA_MERDS_rar.divfil$soil_root_stage=with(SILVA_MERDS_rar.divfil, interaction(soil_status,root_association,life_stage))


#####Transplant Diversity #####
SILVA_MERDS_rar.divfil_trt=subset(SILVA_MERDS_rar.divfil, life_stage!="Start")
#
#let's just look at the treatments
SILVA_MERDS_rar.divfil_trt_trans=subset(SILVA_MERDS_rar.divfil_trt, life_stage=="G")
nrow(SILVA_MERDS_rar.divfil_trt_trans)
#24
bact_obs_rich_model_trans= lm(log(Observed)~soil_root*precip+as.factor(block), data= SILVA_MERDS_rar.divfil_trt_trans)
qqPlot(resid(bact_obs_rich_model_trans))
hist(resid(bact_obs_rich_model_trans))
shapiro.test(resid(bact_obs_rich_model_trans))
#0.9577
Anova(bact_obs_rich_model_trans, type=3)
#soil_root          17.05  2   177.4953 3.625e-11 ***
#precip              0.28  1     5.8152   0.02916 *   
SILVA_MERDS_rar.divfil_trt_trans %>% group_by(precip) %>% summarise_at("Observed", funs(n(),mean,sd,se=sd(.)/sqrt(n())))

emmeans(bact_obs_rich_model_trans, pairwise~soil_root|precip)
emmeans(bact_obs_rich_model_trans, pairwise~soil_root)
#$contrasts
#precip = A:
#  contrast  estimate    SE df t.ratio p.value
#L.B - S.B   1.7243 0.155 15  11.126 <.0001 
#L.B - L.R  -0.0981 0.155 15  -0.633 0.8044 
#S.B - L.R  -1.8224 0.155 15 -11.759 <.0001 

#precip = D:
#  contrast  estimate    SE df t.ratio p.value
#L.B - S.B   1.8919 0.155 15  12.207 <.0001 
#L.B - L.R   0.1792 0.155 15   1.156 0.4961 
#S.B - L.R  -1.7127 0.155 15 -11.051 <.0001 

bact_obs_rich_model_trans_no_block= lm(log(Observed)~soil_root*precip, data= SILVA_MERDS_rar.divfil_trt_trans)
qqPlot(resid(bact_obs_rich_model_trans_no_block))
hist(resid(bact_obs_rich_model_trans_no_block))
shapiro.test(resid(bact_obs_rich_model_trans_no_block))
#0.8875

Anova(bact_obs_rich_model_trans_no_block, type=3)
#soil_root          17.05  2   195.2992 6.249e-13 ***
#precip              0.28  1     6.3985   0.02098 * 

emmeans(bact_obs_rich_model_trans_no_block, pairwise~soil_root|precip)


AIC(bact_obs_rich_model_trans,bact_obs_rich_model_trans_no_block)

SILVA_MERDS_rar.divfil_trt_trans %>% group_by(soil_root,precip) %>% summarise_at("Observed", funs(n(),mean,sd,se=sd(.)/sqrt(n())))


SILVA_MERDS_rar.divfil_trt_trans_precip_soil_g=SILVA_MERDS_rar.divfil_trt_trans %>% group_by(soil_root,precip)
obs_rich_precip_trans=summarise_at(SILVA_MERDS_rar.divfil_trt_trans_precip_soil_g, 
                                   "Observed", funs(n(),mean,sd,se=sd(.)/sqrt(n())))

ggplot(obs_rich_precip_trans, aes(x=precip,y=mean,ymin = mean-se, ymax= mean+se,fill=soil_root))+
  geom_bar(position="dodge",stat="identity", color="black",aes(fill=soil_root))+
  geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+
  scale_fill_manual(values = c("gray", "white", "black"),labels=c("Bulk","Sterile","Rhizo"))+
  scale_x_discrete(labels=c("Ambient","Drought"))+ylab("Bacterial richness")+
  geom_text(aes(y=mean+se+50, label=n),position=position_dodge(width=0.9))+theme_bw()

ggplot(SILVA_MERDS_rar.divfil_trt, aes(x=precip, y=Observed))+geom_boxplot(aes(color=precip))+theme_bw()


#####Trans Observed Richness graph####
treatment_order=c("S.B","L.B","L.R")
(trans_obs_richness_p=ggplot(obs_rich_precip_trans, aes(x=precip,y=mean, ymin = mean-se, ymax= mean+se,fill=factor(soil_root,levels = treatment_order)))+
    geom_bar(position="dodge",stat="identity", color="black",aes(fill=factor(soil_root,levels = treatment_order)))+
    geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+ylim(c(0,1850))+
    scale_fill_manual(values = c( "white","lightgray", "darkgrey"),labels=c("Sterile","Bulk","Rhizo"), limits=treatment_order, name=NULL)+
    scale_x_discrete(labels=c("Ambient","Drought"), name=NULL)+ylab("Bacterial richness")+
    geom_text(aes(y=80, label=n),size=10,position=position_dodge(width=0.9))+theme_bw()+
    theme(axis.title = element_text(size = 23), axis.text = element_text(size = 23),
          legend.position = c(0.85,.9), legend.text=element_text(size=20),
          legend.background = element_rect(size=0.5,linetype="solid",colour ="black")))


(trans_obs_richness_p2=ggplot(obs_rich_precip_trans,aes(x=precip,fill=factor(soil_root,levels = treatment_order)))+
    geom_bar(aes(y=mean),position="dodge",stat="identity", color="black")+
    geom_errorbar(aes(y=mean, ymin = mean-se, ymax= mean+se),position=position_dodge(width=0.9), width=0.2, size=1)+ylim(c(0,1950))+
    scale_fill_manual(values = c( "white","lightgray", "darkgrey"),labels=c("Sterile","Bulk","Rhizo"), limits=treatment_order, name=NULL)+
    scale_x_discrete(labels=c("Ambient","Drought"), name=NULL)+ylab("Bacterial richness")+
    geom_text(aes(y=80, label=n),size=10,position=position_dodge(width=0.9))+theme_bw()+
    theme(axis.title = element_text(size = 23), axis.text = element_text(size = 23),
          legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

#800x750
#####Presence Richness#####
SILVA_MERDS_rar.divfil_pres_trans=subset(SILVA_MERDS_rar.divfil_trt_trans, root_association =="B")
nrow(SILVA_MERDS_rar.divfil_pres_trans)

bact_obs_rich_model_trans_pres= lm(log(Observed)~soil_root*precip+as.factor(block), data= SILVA_MERDS_rar.divfil_pres_trans)
qqPlot(resid(bact_obs_rich_model_trans_pres))
hist(resid(bact_obs_rich_model_trans_pres))
shapiro.test(resid(bact_obs_rich_model_trans_pres))
#0.09522
Anova(bact_obs_rich_model_trans_pres, type=3)
#soil_root          17.05  2   177.4953 3.625e-11 ***
#precip              0.28  1     5.8152   0.02916 *   


emmeans(bact_obs_rich_model_trans_pres, pairwise~soil_root|precip)
emmeans(bact_obs_rich_model_trans_pres, pairwise~soil_root)
#



bact_obs_rich_model_trans_pres_no_block= lm(log(Observed)~soil_root*precip, data= SILVA_MERDS_rar.divfil_pres_trans)
qqPlot(resid(bact_obs_rich_model_trans_pres_no_block))
hist(resid(bact_obs_rich_model_trans_pres_no_block))
shapiro.test(resid(bact_obs_rich_model_trans_pres_no_block))
#0.5744

Anova(bact_obs_rich_model_trans_pres_no_block, type=3)
#soil_root         13.08  1   300.6577 7.332e-10 ***

emmeans(bact_obs_rich_model_trans_pres_no_block, pairwise~soil_root*precip,adjust="none")


AIC(bact_obs_rich_model_trans_pres,bact_obs_rich_model_trans_pres_no_block)


#####Origin Richness#####
SILVA_MERDS_rar.divfil_orig_trans=subset(SILVA_MERDS_rar.divfil_trt_trans, soil_status =="L")
nrow(SILVA_MERDS_rar.divfil_orig_trans)

bact_obs_rich_model_trans_orig= lm((Observed)~soil_root*precip+as.factor(block), data= SILVA_MERDS_rar.divfil_orig_trans)
qqPlot(resid(bact_obs_rich_model_trans_orig))
hist(resid(bact_obs_rich_model_trans_orig))
shapiro.test(resid(bact_obs_rich_model_trans_orig))
#0.09522
Anova(bact_obs_rich_model_trans_orig, type=3)
#precip             302500  1   3.5969   0.09039 .  


emmeans(bact_obs_rich_model_trans_orig, pairwise~soil_root|precip)
emmeans(bact_obs_rich_model_trans_orig, pairwise~soil_root)
#



bact_obs_rich_model_trans_orig_no_block= lm((Observed)~soil_root*precip, data= SILVA_MERDS_rar.divfil_orig_trans)
qqPlot(resid(bact_obs_rich_model_trans_orig_no_block))
hist(resid(bact_obs_rich_model_trans_orig_no_block))
shapiro.test(resid(bact_obs_rich_model_trans_orig_no_block))
#0.834

Anova(bact_obs_rich_model_trans_orig_no_block, type=3)
#precip             302500  1   4.3880   0.05808 .  

emmeans(bact_obs_rich_model_trans_orig_no_block, pairwise~soil_root*precip,adjust="none")


AIC(bact_obs_rich_model_trans_orig,bact_obs_rich_model_trans_pres_no_block)


#Simpson
bact_inv_simp_model_trans= lm((InvSimpson)^(-.5)~soil_root*precip+as.factor(block), data= SILVA_MERDS_rar.divfil_trt_trans)
qqPlot(resid(bact_inv_simp_model_trans))
hist(resid(bact_inv_simp_model_trans))
boxCox(bact_inv_simp_model_trans)
shapiro.test(resid(bact_inv_simp_model_trans))
#0.5401

Anova(bact_inv_simp_model_trans, type=3)
#soil_root        0.06758  2   7.9341   0.00446 ** 
#precip           0.01918  1   4.5038   0.05088 .  
#as.factor(block) 0.00777  3   0.6081   0.61999    

emmeans(bact_inv_simp_model_trans, pairwise~soil_root)
emmeans(bact_inv_simp_model_trans, pairwise~soil_root|precip)
#$contrasts
#precip = A:
#contrast  estimate     SE df t.ratio p.value
#L.B - S.B  -0.0779 0.0461 15 -1.688  0.2418 
#L.B - L.R   0.0324 0.0461 15  0.703  0.7654 
#S.B - L.R   0.1103 0.0461 15  2.391  0.0734 

#precip = D:
#  contrast  estimate     SE df t.ratio p.value
#L.B - S.B  -0.1640 0.0461 15 -3.553  0.0077 
#L.B - L.R  -0.0708 0.0461 15 -1.535  0.3035 
#S.B - L.R   0.0931 0.0461 15  2.018  0.1420 

bact_inv_simp_model_trans_no_block= lm((InvSimpson)^(-.5)~soil_root*precip, data= SILVA_MERDS_rar.divfil_trt_trans)
qqPlot(resid(bact_inv_simp_model_trans_no_block))
hist(resid(bact_inv_simp_model_trans_no_block))
boxCox(bact_inv_simp_model_trans_no_block)
shapiro.test(resid(bact_inv_simp_model_trans_no_block))
#0.3567

Anova(bact_inv_simp_model_trans_no_block, type=3)
#soil_root        0.06758  2   8.4886  0.002532 ** 
#precip           0.01918  1   4.8185  0.041510 * 

emmeans(bact_inv_simp_model_trans_no_block, pairwise~soil_root|precip)

AIC(bact_inv_simp_model_trans,bact_inv_simp_model_trans_no_block)

SILVA_MERDS_rar.divfil_trt_trans %>% group_by(precip,soil_root) %>% summarise_at("InvSimpson", funs(n(),mean,sd,se=sd(.)/sqrt(n())))

SILVA_MERDS_rar.divfil_trt_trans_precip_soil_g=SILVA_MERDS_rar.divfil_trt_trans %>% group_by(soil_root,precip)
inv_simp_precip_trans=summarise_at(SILVA_MERDS_rar.divfil_trt_trans_precip_soil_g, 
                                   "InvSimpson", funs(n(),mean,sd,se=sd(.)/sqrt(n())))

ggplot(inv_simp_precip_trans, aes(x=precip,y=mean,ymin = mean-se, ymax= mean+se,fill=soil_root))+
  geom_bar(position="dodge",stat="identity", color="black",aes(fill=soil_root))+
  geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+
  scale_fill_manual(values = c("gray", "white", "black"),labels=c("Bulk","Sterile","Rhizo"))+
  scale_x_discrete(labels=c("Ambient","Drought"))+ylab("Inv Simpson bacteria")+
  geom_text(aes(y=mean+se+5, label=n),position=position_dodge(width=0.9))+theme_bw()


#####Trans inv Simpson graph####

(trans_inv_simpson_p=ggplot(inv_simp_precip_trans, aes(x=precip,y=mean, ymin = mean-se, ymax= mean+se,fill=factor(soil_root,levels = treatment_order)))+
   geom_bar(position="dodge",stat="identity", color="black",aes(fill=factor(soil_root,levels = treatment_order)))+
   geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+ylim(c(-10,200))+
   scale_fill_manual(values = c( "white","lightgray", "darkgrey"),labels=c("Sterile","Bulk","Rhizo"), limits=treatment_order, name=NULL)+
   scale_x_discrete(labels=c("Ambient","Drought"), name=NULL)+ylab("inv Simpson bacteria")+
   geom_text(aes(y=-7, label=n),size=10,position=position_dodge(width=0.9))+theme_bw()+
   theme(axis.title = element_text(size = 23), axis.text = element_text(size = 23),
         legend.position = c(0.85,.9), legend.text=element_text(size=20),
         legend.background = element_rect(size=0.5,linetype="solid",colour ="black")))

(trans_inv_simpson_p2=ggplot(inv_simp_precip_trans, aes(x=precip,fill=factor(soil_root,levels = treatment_order)))+
    geom_bar(aes(y=mean),position="dodge",stat="identity", color="black")+
    geom_errorbar(aes(y=mean, ymin = mean-se, ymax= mean+se),
                  position=position_dodge(width=0.9), width=0.2, size=1)+ylim(c(-10,220))+
    scale_fill_manual(values = c( "white","lightgray", "darkgrey"),labels=c("Sterile","Bulk","Rhizo"), limits=treatment_order, name=NULL)+
    scale_x_discrete(labels=c("Ambient","Drought"), name=NULL)+ylab("inv Simpson bacteria")+
    geom_text(aes(y=-7, label=n),size=10,position=position_dodge(width=0.9))+theme_bw()+
    theme(axis.title = element_text(size = 23), axis.text = element_text(size = 23),
          legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank()))


#800x750

#####Presence Simpson#####
SILVA_MERDS_rar.divfil_pres_trans=subset(SILVA_MERDS_rar.divfil_trt_trans, root_association =="B")
nrow(SILVA_MERDS_rar.divfil_pres_trans)

#Simpson
bact_inv_simp_model_trans_pres= lm((InvSimpson)^(-.5)~soil_root*precip+as.factor(block), data= SILVA_MERDS_rar.divfil_pres_trans)
qqPlot(resid(bact_inv_simp_model_trans_pres))
hist(resid(bact_inv_simp_model_trans_pres))
boxCox(bact_inv_simp_model_trans_pres)
shapiro.test(resid(bact_inv_simp_model_trans_pres))
#0.9784

Anova(bact_inv_simp_model_trans_pres, type=3)
#soil_root        0.05849  1  25.1118 0.0007278 ***

emmeans(bact_inv_simp_model_trans_pres, pairwise~soil_root)
emmeans(bact_inv_simp_model_trans_pres, pairwise~soil_root|precip)
#

bact_inv_simp_model_trans_pres_no_block= lm((InvSimpson)^(-.5)~soil_root*precip, data= SILVA_MERDS_rar.divfil_pres_trans)
qqPlot(resid(bact_inv_simp_model_trans_pres_no_block))
hist(resid(bact_inv_simp_model_trans_pres_no_block))
boxCox(bact_inv_simp_model_trans_pres_no_block)
shapiro.test(resid(bact_inv_simp_model_trans_pres_no_block))
# 0.4765

Anova(bact_inv_simp_model_trans_pres_no_block, type=3)
#soil_root        0.05849  1  22.8916 0.0004451 ***

emmeans(bact_inv_simp_model_trans_pres_no_block, pairwise~soil_root*precip, adjust="none")

AIC(bact_inv_simp_model_trans_pres,bact_inv_simp_model_trans_pres_no_block)


#####Origin Simpson#####
SILVA_MERDS_rar.divfil_orig_trans=subset(SILVA_MERDS_rar.divfil_trt_trans, soil_status =="L")
nrow(SILVA_MERDS_rar.divfil_orig_trans)

#Simpson
bact_inv_simp_model_trans_orig= lm((InvSimpson)~soil_root*precip+as.factor(block), data= SILVA_MERDS_rar.divfil_orig_trans)
qqPlot(resid(bact_inv_simp_model_trans_orig))
hist(resid(bact_inv_simp_model_trans_orig))
boxCox(bact_inv_simp_model_trans_orig)
shapiro.test(resid(bact_inv_simp_model_trans_orig))
#0.9785

Anova(bact_inv_simp_model_trans_orig, type=3)
#nada

emmeans(bact_inv_simp_model_trans_orig, pairwise~soil_root)
emmeans(bact_inv_simp_model_trans_orig, pairwise~soil_root|precip)
#

bact_inv_simp_model_trans_orig_no_block= lm((InvSimpson)~soil_root*precip, data= SILVA_MERDS_rar.divfil_orig_trans)
qqPlot(resid(bact_inv_simp_model_trans_orig_no_block))
hist(resid(bact_inv_simp_model_trans_orig_no_block))
boxCox(bact_inv_simp_model_trans_orig_no_block)
shapiro.test(resid(bact_inv_simp_model_trans_orig_no_block))
# 0.4765

Anova(bact_inv_simp_model_trans_orig_no_block, type=3)
#soil_root:precip  13801  1  3.6577 0.0799711 . 

emmeans(bact_inv_simp_model_trans_orig_no_block, pairwise~soil_root*precip, adjust="none")

AIC(bact_inv_simp_model_trans_orig,bact_inv_simp_model_trans_orig_no_block)



#Seeds

#let's just look at the treatments
SILVA_MERDS_rar.divfil_trt_seed=subset(SILVA_MERDS_rar.divfil_trt, life_stage=="S")
nrow(SILVA_MERDS_rar.divfil_trt_seed)
#23
bact_obs_rich_model_seed= lm(log(Observed)~soil_root*precip+as.factor(block), data= SILVA_MERDS_rar.divfil_trt_seed)
qqPlot(resid(bact_obs_rich_model_seed))
hist(resid(bact_obs_rich_model_seed))
boxCox(bact_obs_rich_model_seed)
shapiro.test(resid(bact_obs_rich_model_seed))
#0.5636
Anova(bact_obs_rich_model_seed, type=3)
#soil_root         16.28  2   318.5927 2.123e-12 ***


emmeans(bact_obs_rich_model_seed, pairwise~soil_root)


emmeans(bact_obs_rich_model_seed, pairwise~soil_root|precip)
#$contrasts
#precip = A:
#  contrast  estimate    SE df t.ratio p.value
#L.B - S.B    1.724 0.113 14  15.253 <.0001 
#L.B - L.R   -0.301 0.113 14  -2.664 0.0458 
#S.B - L.R   -2.025 0.113 14 -17.917 <.0001 

#precip = D:
#  contrast  estimate    SE df t.ratio p.value
#L.B - S.B    1.883 0.129 14  14.602 <.0001 
#L.B - L.R   -0.290 0.116 14  -2.499 0.0622 
#S.B - L.R   -2.173 0.138 14 -15.697 <.0001 

bact_obs_rich_model_seed_no_block= lm(log(Observed)~soil_root*precip, data= SILVA_MERDS_rar.divfil_trt_seed)
qqPlot(resid(bact_obs_rich_model_seed_no_block))
hist(resid(bact_obs_rich_model_seed_no_block))
boxCox(bact_obs_rich_model_seed_no_block)
shapiro.test(resid(bact_obs_rich_model_seed_no_block))
#0.6233
Anova(bact_obs_rich_model_seed_no_block, type=3)
#soil_root         18.95  2   344.8049 1.741e-14 ***

emmeans(bact_obs_rich_model_seed_no_block, pairwise~soil_root|precip)
"$contrasts
precip = A:
 contrast  estimate    SE df t.ratio p.value
 L.B - S.B    1.724 0.117 17  14.706 <.0001 
 L.B - L.R   -0.301 0.117 17  -2.569 0.0497 
 S.B - L.R   -2.025 0.117 17 -17.275 <.0001 

precip = D:
 contrast  estimate    SE df t.ratio p.value
 L.B - S.B    1.897 0.127 17  14.985 <.0001 
 L.B - L.R   -0.323 0.117 17  -2.752 0.0345 
 S.B - L.R   -2.220 0.127 17 -17.533 <.0001 
"

AIC(bact_obs_rich_model_seed,bact_obs_rich_model_seed_no_block)

SILVA_MERDS_rar.divfil_trt_seed %>% group_by(precip,soil_root) %>% summarise_at("Observed", funs(n(),mean,sd,se=sd(.)/sqrt(n())))


SILVA_MERDS_rar.divfil_trt_seed_precip_soil_g=SILVA_MERDS_rar.divfil_trt_seed %>% group_by(soil_root,precip)
obs_rich_precip_seed=summarise_at(SILVA_MERDS_rar.divfil_trt_seed_precip_soil_g, 
                                  "Observed", funs(n(),mean,sd,se=sd(.)/sqrt(n())))

ggplot(obs_rich_precip_seed, aes(x=precip,y=mean,ymin = mean-se, ymax= mean+se,fill=soil_root))+
  geom_bar(position="dodge",stat="identity", color="black",aes(fill=soil_root))+
  geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+
  scale_fill_manual(values = c("gray", "white", "black"),labels=c("Bulk","Sterile","Rhizo"))+
  scale_x_discrete(labels=c("Ambient","Drought"))+ylab("Bacterial richness")+
  geom_text(aes(y=mean+se+50, label=n),position=position_dodge(width=0.9))+theme_bw()

ggplot(SILVA_MERDS_rar.divfil_trt, aes(x=precip, y=Observed))+geom_boxplot(aes(color=precip))+theme_bw()


#####Seed Oberserved Richness graph####

(seed_obs_richness_p=ggplot(obs_rich_precip_seed, aes(x=precip,y=mean, ymin = mean-se, ymax= mean+se,fill=factor(soil_root,levels = treatment_order)))+
   geom_bar(position="dodge",stat="identity", color="black",aes(fill=factor(soil_root,levels = treatment_order)))+
   geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+ylim(c(0,1850))+
   scale_fill_manual(values = c( "white","lightgray", "darkgrey"),labels=c("Sterile","Bulk","Rhizo"), limits=treatment_order, name=NULL)+
   scale_x_discrete(labels=c("Ambient","Drought"), name=NULL)+ylab("Bacterial richness")+
   geom_text(aes(y=80, label=n),size=10,position=position_dodge(width=0.9))+theme_bw()+
   theme(axis.title.x = element_text(size = 23), axis.text.x = element_text(size = 23),
         axis.title.y = element_blank(), axis.text.y = element_blank(),
         legend.position = c(0.85,.9), legend.text=element_text(size=20),
         legend.background = element_rect(size=0.5,linetype="solid",colour ="black")))


(seed_obs_richness_p2=ggplot(obs_rich_precip_seed, aes(x=precip,fill=factor(soil_root,levels = treatment_order)))+
    geom_bar(aes(y=mean),position="dodge",stat="identity", color="black")+
    geom_errorbar(aes(y=mean, ymin = mean-se, ymax= mean+se),position=position_dodge(width=0.9), width=0.2, size=1)+ylim(c(0,1950))+
    scale_fill_manual(values = c( "white","lightgray", "darkgrey"),labels=c("Sterile","Bulk","Rhizo"), limits=treatment_order, name=NULL)+
    scale_x_discrete(labels=c("Ambient","Drought"), name=NULL)+ylab("Bacterial richness")+
    geom_text(aes(y=80, label=n),size=10,position=position_dodge(width=0.9))+theme_bw()+
    theme(axis.title.x = element_text(size = 23), axis.text.x = element_text(size = 23),
          axis.title.y = element_text(size = 23), axis.text.y = element_text(size = 23),
          legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank()))
#bact_richness_germination

#800x750

#####Combined obs rich graph#####
ggarrange(trans_obs_richness_p,seed_obs_richness_p,ncol = 2,  legend = "none",widths = c(1,.90))
#15x7.38

#####Presence Richness#####
SILVA_MERDS_rar.divfil_pres_seed=subset(SILVA_MERDS_rar.divfil_trt_seed, root_association =="B")
nrow(SILVA_MERDS_rar.divfil_pres_seed)

bact_obs_rich_model_seed_pres= lm((Observed)~soil_root*precip+as.factor(block), data= SILVA_MERDS_rar.divfil_pres_seed)
qqPlot(resid(bact_obs_rich_model_seed_pres))
hist(resid(bact_obs_rich_model_seed_pres))
shapiro.test(resid(bact_obs_rich_model_seed_pres))
#0.9939
Anova(bact_obs_rich_model_seed_pres, type=3)
#soil_root        2508022  1 173.9216 1.041e-06 ***
#as.factor(block)  204800  3   4.7340   0.03497 *  

emmeans(bact_obs_rich_model_seed_pres, pairwise~soil_root|precip)
emmeans(bact_obs_rich_model_seed_pres, pairwise~soil_root)
#

emmeans(bact_obs_rich_model_seed_pres, pairwise~soil_root*precip,adjust="none")


bact_obs_rich_model_seed_pres_no_block= lm((Observed)~soil_root*precip, data= SILVA_MERDS_rar.divfil_pres_seed)
qqPlot(resid(bact_obs_rich_model_seed_pres_no_block))
hist(resid(bact_obs_rich_model_seed_pres_no_block))
shapiro.test(resid(bact_obs_rich_model_seed_pres_no_block))
#0.5744

Anova(bact_obs_rich_model_seed_pres_no_block, type=3)
#soil_root        2742486  1  94.2249 9.944e-07 ***

emmeans(bact_obs_rich_model_seed_pres_no_block, pairwise~soil_root*precip,adjust="none")


AIC(bact_obs_rich_model_seed_pres,bact_obs_rich_model_seed_pres_no_block)


#####Origin Richness#####
SILVA_MERDS_rar.divfil_orig_seed=subset(SILVA_MERDS_rar.divfil_trt_seed, soil_status =="L")
nrow(SILVA_MERDS_rar.divfil_orig_seed)

bact_obs_rich_model_seed_orig= lm((Observed)~soil_root*precip+as.factor(block), data= SILVA_MERDS_rar.divfil_orig_seed)
qqPlot(resid(bact_obs_rich_model_seed_orig))
hist(resid(bact_obs_rich_model_seed_orig))
shapiro.test(resid(bact_obs_rich_model_seed_orig))
#0.5203
Anova(bact_obs_rich_model_seed_orig, type=3)
#soil_root          444042  1  10.9369 0.009126 ** 


emmeans(bact_obs_rich_model_seed_orig, pairwise~soil_root|precip)
emmeans(bact_obs_rich_model_seed_orig, pairwise~soil_root)
#



bact_obs_rich_model_seed_orig_no_block= lm((Observed)~soil_root*precip, data= SILVA_MERDS_rar.divfil_orig_seed)
qqPlot(resid(bact_obs_rich_model_seed_orig_no_block))
hist(resid(bact_obs_rich_model_seed_orig_no_block))
shapiro.test(resid(bact_obs_rich_model_seed_orig_no_block))
#0.7012

Anova(bact_obs_rich_model_seed_orig_no_block, type=3)
#soil_root          502327  1  13.8057  0.00295 **  

emmeans(bact_obs_rich_model_seed_orig_no_block, pairwise~soil_root*precip,adjust="none")


AIC(bact_obs_rich_model_seed_orig,bact_obs_rich_model_seed_orig_no_block)

#Simpson
bact_inv_simp_model_seed= lm(log(InvSimpson)~soil_root*precip+as.factor(block), data= SILVA_MERDS_rar.divfil_trt_seed)
qqPlot(resid(bact_inv_simp_model_seed))
hist(resid(bact_inv_simp_model_seed))
boxCox(bact_inv_simp_model_seed)
shapiro.test(resid(bact_inv_simp_model_seed))
#0.4042
Anova(bact_inv_simp_model_seed, type=3)
#soil_root         17.314  2  11.9446 0.0009404 ***


emmeans(bact_inv_simp_model_seed, pairwise~soil_root)
# contrast  estimate    SE df t.ratio p.value
#L.B - S.B    0.519 0.457 14  1.136  0.5087 
#L.B - L.R   -1.645 0.431 14 -3.815  0.0050 
#S.B - L.R   -2.164 0.476 14 -4.547  0.0012

emmeans(bact_inv_simp_model_seed, pairwise~soil_root|precip)

#$contrasts
#precip = A:
#contrast  estimate    SE df t.ratio p.value
#L.B - S.B    0.824 0.602 14  1.368  0.3832 
#L.B - L.R   -1.711 0.602 14 -2.842  0.0328 
#S.B - L.R   -2.535 0.602 14 -4.211  0.0024 

#precip = D:
#  contrast  estimate    SE df t.ratio p.value
#L.B - S.B    0.214 0.687 14  0.312  0.9480 
#L.B - L.R   -1.579 0.617 14 -2.557  0.0559 
#S.B - L.R   -1.793 0.737 14 -2.432  0.0702 


bact_inv_simp_model_seed_no_block= lm(log(InvSimpson)~soil_root*precip, data= SILVA_MERDS_rar.divfil_trt_seed)
qqPlot(resid(bact_inv_simp_model_seed_no_block))
hist(resid(bact_inv_simp_model_seed_no_block))
boxCox(bact_inv_simp_model_seed_no_block)
shapiro.test(resid(bact_inv_simp_model_seed_no_block))
#0.4735
Anova(bact_inv_simp_model_seed_no_block, type=3)
#soil_root         21.167  2  16.4729 0.0001051 ***


emmeans(bact_inv_simp_model_seed_no_block, pairwise~soil_root|precip)



AIC(bact_inv_simp_model_seed,bact_inv_simp_model_seed_no_block)

SILVA_MERDS_rar.divfil_trt_seed_precip_soil_g=SILVA_MERDS_rar.divfil_trt_seed %>% group_by(soil_root,precip)

SILVA_MERDS_rar.divfil_trt_seed %>% group_by(precip,soil_root) %>% summarise_at("InvSimpson", funs(n(),mean,sd,se=sd(.)/sqrt(n())))

inv_simp_precip_seed=summarise_at(SILVA_MERDS_rar.divfil_trt_seed_precip_soil_g, 
                                  "InvSimpson", funs(n(),mean,sd,se=sd(.)/sqrt(n())))

ggplot(inv_simp_precip_seed, aes(x=precip,y=mean,ymin = mean-se, ymax= mean+se,fill=soil_root))+
  geom_bar(position="dodge",stat="identity", color="black",aes(fill=soil_root))+
  geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+
  scale_fill_manual(values = c("gray", "white", "black"),labels=c("Bulk","Sterile","Rhizo"))+
  scale_x_discrete(labels=c("Ambient","Drought"))+ylab("Inv Simpson bacteria")+
  geom_text(aes(y=mean+se+5, label=n),position=position_dodge(width=0.9))+theme_bw()


#####Seed inv Simpson graph####

(seed_inv_simpson_p=ggplot(inv_simp_precip_seed, aes(x=precip,y=mean, ymin = mean-se, ymax= mean+se,fill=factor(soil_root,levels = treatment_order)))+
   geom_bar(position="dodge",stat="identity", color="black",aes(fill=factor(soil_root,levels = treatment_order)))+
   geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+ylim(c(-10,200))+
   scale_fill_manual(values = c( "white","lightgray", "darkgrey"),labels=c("Sterile","Bulk","Rhizo"), limits=treatment_order, name=NULL)+
   scale_x_discrete(labels=c("Ambient","Drought"), name=NULL)+ylab("inv Simpson bacteria")+
   geom_text(aes(y=-7, label=n),size=10,position=position_dodge(width=0.9))+theme_bw()+
   theme(axis.title.x = element_text(size = 23), axis.text.x = element_text(size = 23),
         axis.title.y = element_blank(), axis.text.y = element_blank(),
         legend.position = c(0.85,.9), legend.text=element_text(size=20),
         legend.background = element_rect(size=0.5,linetype="solid",colour ="black")))


(seed_inv_simpson_p=ggplot(inv_simp_precip_seed, aes(x=precip,y=mean, ymin = mean-se, ymax= mean+se,fill=factor(soil_root,levels = treatment_order)))+
    geom_bar(position="dodge",stat="identity", color="black",aes(fill=factor(soil_root,levels = treatment_order)))+
    geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+ylim(c(-10,200))+
    scale_fill_manual(values = c( "white","lightgray", "darkgrey"),labels=c("Sterile","Bulk","Rhizo"), limits=treatment_order, name=NULL)+
    scale_x_discrete(labels=c("Ambient","Drought"), name=NULL)+ylab("inv Simpson bacteria")+
    geom_text(aes(y=-7, label=n),size=10,position=position_dodge(width=0.9))+theme_bw()+
    theme(axis.title.x = element_text(size = 23), axis.text.x = element_text(size = 23),
          axis.title.y = element_text(size = 23), axis.text.y = element_text(size = 23),
          legend.position = "none"))

(seed_inv_simpson_p2=ggplot(inv_simp_precip_seed, aes(x=precip,fill=factor(soil_root,levels = treatment_order)))+
    geom_bar(aes(y=mean),position="dodge",stat="identity", color="black")+
    geom_errorbar(aes(y=mean, ymin = mean-se, ymax= mean+se),position=position_dodge(width=0.9), width=0.2, size=1)+ylim(c(-10,220))+
    scale_fill_manual(values = c( "white","lightgray", "darkgrey"),labels=c("Sterile","Bulk","Rhizo"), limits=treatment_order, name=NULL)+
    scale_x_discrete(labels=c("Ambient","Drought"), name=NULL)+ylab("inv Simpson bacteria")+
    geom_text(aes(y=-7, label=n),size=10,position=position_dodge(width=0.9))+theme_bw()+
    theme(axis.title.x = element_text(size = 23), axis.text.x = element_text(size = 23),
          axis.title.y = element_text(size = 23), axis.text.y = element_text(size = 23),
          legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank()))


#800x750

#####Combined inv simp graph#####
ggarrange(trans_inv_simpson_p,seed_inv_simpson_p,ncol = 2,  legend = "none",widths = c(1,.90))
#15x7.38


#####Presence Simpson#####
SILVA_MERDS_rar.divfil_pres_seed=subset(SILVA_MERDS_rar.divfil_trt_seed, root_association =="B")
nrow(SILVA_MERDS_rar.divfil_pres_seed)

bact_inv_simp_model_seed_pres= lm(log(InvSimpson)~soil_root*precip+as.factor(block), data= SILVA_MERDS_rar.divfil_pres_seed)
qqPlot(resid(bact_inv_simp_model_seed_pres))
hist(resid(bact_inv_simp_model_seed_pres))
boxCox(bact_inv_simp_model_seed_pres)
shapiro.test(resid(bact_inv_simp_model_seed_pres))
#0.08898
Anova(bact_inv_simp_model_seed_pres, type=3)
#soil_root         17.314  2  11.9446 0.0009404 ***


emmeans(bact_inv_simp_model_seed_pres, pairwise~soil_root)
#

emmeans(bact_inv_simp_model_seed_pres, pairwise~soil_root|precip)

#


bact_inv_simp_model_seed_pres_no_block= lm(log(InvSimpson)~soil_root*precip, data= SILVA_MERDS_rar.divfil_pres_seed)
qqPlot(resid(bact_inv_simp_model_seed_pres_no_block))
hist(resid(bact_inv_simp_model_seed_pres_no_block))
boxCox(bact_inv_simp_model_seed_pres_no_block)
shapiro.test(resid(bact_inv_simp_model_seed_pres_no_block))
#0.4478
Anova(bact_inv_simp_model_seed_pres_no_block, type=3)
#soil_root         21.167  2  16.4729 0.0001051 ***


emmeans(bact_inv_simp_model_seed_pres_no_block, pairwise~soil_root|precip)



AIC(bact_inv_simp_model_seed_pres,bact_inv_simp_model_seed_pres_no_block)


#####Origin Simpson#####
SILVA_MERDS_rar.divfil_orig_seed=subset(SILVA_MERDS_rar.divfil_trt_seed, soil_status =="L")
nrow(SILVA_MERDS_rar.divfil_orig_seed)

bact_inv_simp_model_seed_orig= lm((InvSimpson)~soil_root*precip+as.factor(block), data= SILVA_MERDS_rar.divfil_orig_seed)
qqPlot(resid(bact_inv_simp_model_seed_orig))
hist(resid(bact_inv_simp_model_seed_orig))
boxCox(bact_inv_simp_model_seed_orig)
shapiro.test(resid(bact_inv_simp_model_seed_orig))
#0.9253
Anova(bact_inv_simp_model_seed_orig, type=3)
#soil_root         25230  1  8.9340 0.015226 * 


emmeans(bact_inv_simp_model_seed_orig, pairwise~soil_root)
#

emmeans(bact_inv_simp_model_seed_orig, pairwise~soil_root|precip)

#


bact_inv_simp_model_seed_orig_no_block= lm((InvSimpson)~soil_root*precip, data= SILVA_MERDS_rar.divfil_orig_seed)
qqPlot(resid(bact_inv_simp_model_seed_orig_no_block))
hist(resid(bact_inv_simp_model_seed_orig_no_block))
boxCox(bact_inv_simp_model_seed_orig_no_block)
shapiro.test(resid(bact_inv_simp_model_seed_orig_no_block))
#0.8517
Anova(bact_inv_simp_model_seed_orig_no_block, type=3)
#soil_root         27358  1 12.5228 0.0040802 ** 


emmeans(bact_inv_simp_model_seed_orig_no_block, pairwise~soil_root|precip)

emmeans(bact_inv_simp_model_seed_orig_no_block, pairwise~soil_root*precip,adjust="none")

AIC(bact_inv_simp_model_seed_orig,bact_inv_simp_model_seed_orig_no_block)
