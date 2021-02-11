######Here is the Code for analyzing the bacterial community composition of MERDS Switchgrass Experiment#######
#Lukas (belldere@msu.edu) processed the sequences and subset the larger run to only have samples from MERDS (MERDs 56 of 317 samples in the run)

library(phyloseq)
library(RVAideMemoire) # for tests after PERMANOVA
library(ggplot2)
library(ggrepel)
library(microbiome)
library(dplyr)
library(limma)
library(gridExtra)
library(car)
library(emmeans)
library(vegan)
library(ggpubr)
library(DESeq2)
library(otuSummary)
source(system.file("extdata/lm_phyloseq.R", package = "microbiome"))
options(contrasts=c("contr.sum", "contr.poly"))
data_SG_biomass <- read.csv("D:/MERDS_2018/merds/Switchgrass/R_data/SG_TotalBiomass.csv")
SG_inorg_N<- read.csv("D:/MERDS_2018/merds/Switchgrass/R_data/SG_NO3NH4_ug_gdrysoil.csv", header = T)
summary(SG_inorg_N)
dataSG_seed_surv <- read.csv("D:/MERDS_2018/merds/Switchgrass/R_data/SG_surv_Seed_germ.csv")
summary(dataSG_seed_surv)
dataSG_seed_surv[,8:14]=NULL

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

#Tree 
run_20190617_16S.tree = read_tree("D:/MMPRNT_16S_016-018/USEARCH/tree_OTU_20190617_16S-V4_PE250_NWK.NWK")


#GTDB based taxonomy https://gtdb.ecogenomic.org/

taxa_raw_GTDBr89_MERDs= read.delim("D:/MERDS_2018/merds/Switchgrass/R_data/taxonomy.tsv",sep = c("\t"),header = T)
head(taxa_raw_GTDBr89_MERDs)
taxa_raw_GTDBr89_MERDs_sep=taxa_raw_GTDBr89_MERDs %>% separate("Taxon", c("Domain","Phylum","Class","Order","Family","Genus","Species"),sep = ";")
row.names(taxa_raw_GTDBr89_MERDs_sep)=taxa_raw_GTDBr89_MERDs_sep$Feature.ID
taxa_raw_GTDBr89_MERDs_sep$Feature.ID=NULL
taxa_raw_GTDBr89_MERDs_sep[is.na(taxa_raw_GTDBr89_MERDs_sep)] <- "UNKNOWN"
head(taxa_raw_GTDBr89_MERDs_sep)
min(taxa_raw_GTDBr89_MERDs_sep$Confidence)
unique(taxa_raw_GTDBr89_MERDs_sep$Domain)
unique(taxa_raw_GTDBr89_MERDs_sep$Phylum)
taxa_raw_GTDBr89_MERDs_sep$Phylum=gsub("a_A","a",taxa_raw_GTDBr89_MERDs_sep$Phylum)
taxa_raw_GTDBr89_MERDs_sep$Phylum=gsub("a_B","a",taxa_raw_GTDBr89_MERDs_sep$Phylum)

taxa_raw_GTDBr89_MERDs_sep$Phylum=gsub("s_A","s",taxa_raw_GTDBr89_MERDs_sep$Phylum)
taxa_raw_GTDBr89_MERDs_sep$Phylum=gsub("s_B","s",taxa_raw_GTDBr89_MERDs_sep$Phylum)
taxa_raw_GTDBr89_MERDs_sep$Phylum=gsub("s_C","s",taxa_raw_GTDBr89_MERDs_sep$Phylum)
taxa_raw_GTDBr89_MERDs_sep$Phylum=gsub("s_D","s",taxa_raw_GTDBr89_MERDs_sep$Phylum)
taxa_raw_GTDBr89_MERDs_sep$Phylum=gsub("s_E","s",taxa_raw_GTDBr89_MERDs_sep$Phylum)
taxa_raw_GTDBr89_MERDs_sep$Phylum=gsub("s_G","s",taxa_raw_GTDBr89_MERDs_sep$Phylum)
taxa_raw_GTDBr89_MERDs_sep$Phylum=gsub("s_I","s",taxa_raw_GTDBr89_MERDs_sep$Phylum)
taxa_raw_GTDBr89_MERDs_sep$Phylum=gsub("s_K","s",taxa_raw_GTDBr89_MERDs_sep$Phylum)


taxa_raw_GTDBr89_MERDs_sep_mat=as.matrix(taxa_raw_GTDBr89_MERDs_sep)
head(taxa_raw_GTDBr89_MERDs_sep_mat)
TAXA_GTDBr89_MERDs_confid=tax_table(taxa_raw_GTDBr89_MERDs_sep_mat)

#####Begin OTU based community####

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

#I updated the mapping file with the treatments from the experiment


merds_map_trt=read.csv("D:/MERDS_2018/merds/Switchgrass/R_data/map_MERDS_trt.csv", header = T, row.names = 1)
summary(merds_map_trt)
#Need to convert plant number to numeric
merds_map_trt$Plant_Number <- as.numeric(levels(merds_map_trt$Plant_Number)) [merds_map_trt$Plant_Number] # changing factor to numeric
summary(merds_map_trt)

#Let's also add in the biomass data and nitrogen data and time to germination 
merds_map_trt_bio=merge(merds_map_trt,data_SG_biomass, by="Plant_Number",all.x = T)
summary(merds_map_trt_bio)

#Biomass
merds_map_trt_bio_nit=merge(merds_map_trt_bio,SG_inorg_N, by="Plant_Number",all.x = T)
summary(merds_map_trt_bio_nit)
merds_map_trt_bio_nit$surv_germ=merds_map_trt_bio_nit$shoot_weight_g
merds_map_trt_bio_nit$surv_germ[merds_map_trt_bio_nit$surv_germ>0]=1
merds_map_trt_bio_nit$surv_germ[is.na(merds_map_trt_bio_nit$surv_germ)]=0
merds_map_trt_bio_nit$total_biomass[is.na(merds_map_trt_bio_nit$total_biomass)]=0
merds_map_trt_bio$shoot_weight_g[is.na(merds_map_trt_bio_nit$shoot_weight_g)]=0
merds_map_trt_bio_nit$root_weight_g[is.na(merds_map_trt_bio_nit$root_weight_g)]=0
row.names(merds_map_trt_bio_nit)=merds_map_trt_bio$SampleID
merds_map_trt_bio_nit$soil_root=with(merds_map_trt_bio_nit, interaction(soil_status,root_association))
summary(merds_map_trt_bio_nit)
#Germination
dataSG_seed_surv$pot_w_germ=dataSG_seed_surv$num_germinates



dataSG_seed_surv$pot_w_germ[dataSG_seed_surv$pot_w_germ>0]=1
dataSG_seed_surv_1=subset(dataSG_seed_surv, pot_w_germ==1)

dataSG_seed_surv_1_pot_g=group_by(dataSG_seed_surv_1, Plant_Number)

dataSG_seed_surv_first_germ=summarise_at(dataSG_seed_surv_1_pot_g, "exp_days", min)
summary(dataSG_seed_surv_first_germ)
colnames(dataSG_seed_surv_first_germ)[2]="days_to_germ"
merds_map_trt_bio_nit_germ=merge(merds_map_trt_bio_nit, dataSG_seed_surv_first_germ, by="Plant_Number", all.x = T)
#Let's add in the germination results



summary(dataSG_seed_surv)
dataSG_seed_surv_1=subset(dataSG_seed_surv, pot_w_germ==1)

dataSG_seed_surv_1_pot_g=group_by(dataSG_seed_surv_1, Plant_Number)

dataSG_seed_surv_1_pot_first_germ=summarise_at(dataSG_seed_surv_1_pot_g, "exp_days", min)


dataSG_seed_1st_germ_trt=merge(dataSG_seed_surv_1_pot_first_germ, merds_map_trt, by="Plant_Number", all.y = T)
dataSG_seed_1st_germ_trt$soil_root=with(dataSG_seed_1st_germ_trt, interaction(soil_status,root_association))
summary(dataSG_seed_1st_germ_trt)
row.names(dataSG_seed_1st_germ_trt)=dataSG_seed_1st_germ_trt$SampleID

#Quick analyses of Nitrogen diff
merds_map_trt_bio_nit %>% group_by(soil_root)  %>% summarise_at(c("NH4ppm_negto0","NO3ppm"), ~mean(., na.rm = TRUE))
#  soil_root NH4ppm_negto0 NO3ppm
#<fct>             <dbl>  <dbl>
#1 L.B              0.0188  0.807
#2 S.B              0.348   0.620
#3 L.R              0.0291  0.791


#what does the distribution of germination look like

merds_map_trt_bio_nit_germ=subset(merds_map_trt_bio_nit,life_stage=="S" )
nrow(merds_map_trt_bio_nit_germ)
#24

merds_map_trt_bio_nit_germ %>% group_by(soil_root,precip)  %>% summarise_at("surv_germ", ~mean(., na.rm = TRUE))
head(merds_map_trt_bio_nit_germ)

#Let's look at the number of germinants

#write.csv(fin_dataSG_seed_surv_trt_tot_seedling, "D:/MERDS_2018/merds/Switchgrass/R_data/fin_dataSG_seed_surv_trt_tot_seedling.csv")

fin_dataSG_seed_surv_trt_tot_seedling=read.csv("D:/MERDS_2018/merds/Switchgrass/R_data/fin_dataSG_seed_surv_trt_tot_seedling.csv")
head(fin_dataSG_seed_surv_trt_tot_seedling)

fin_dataSG_seed_surv_trt_tot_seedling_micro=merge(merds_map_trt_bio_nit_germ[,c("Plant_Number","surv_germ")],fin_dataSG_seed_surv_trt_tot_seedling,
                                                  by="Plant_Number",all.x=T)
nrow(fin_dataSG_seed_surv_trt_tot_seedling_micro)
#24

hist(fin_dataSG_seed_surv_trt_tot_seedling_micro$tot_num_germ)
hist(fin_dataSG_seed_surv_trt_tot_seedling_micro$t0_germ)
SG_total_germination_soil_root_precip=fin_dataSG_seed_surv_trt_tot_seedling_micro %>% group_by(soil_root,precip)  %>% summarise_at("tot_num_germ", c(~mean(.),se=~sd(.)/sqrt(n())))


treatment_order=c("S.B","L.B","L.R")
(total_germination_p2=ggplot(SG_total_germination_soil_root_precip, aes(x=precip,y=mean,fill=factor(soil_root,levels = treatment_order)))+
    geom_bar(position="dodge",stat="identity", color="black",aes(fill=factor(soil_root,levels = treatment_order)))+
    geom_errorbar(aes( ymin = mean-se, ymax= mean+se),position=position_dodge(width=0.9), width=0.2, size=1)+
    scale_fill_manual(values = c( "white","lightgray", "darkgrey"),labels=c("Sterile","Bulk","Rhizo"), limits=treatment_order, name=NULL)+
    scale_x_discrete(labels=c("Ambient","Drought"), name=NULL)+ylab("Total number of germinants")+
    theme_bw()+theme(axis.title = element_text(size = 23), axis.text = element_text(size = 23),
          legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank()))


#####Begin SILVA dataset v123####
######Begin Processing#### 
#now put the new mapping file into our phyloseq obj

SILVA_MERDS_data=phyloseq(otu_table(phyl_SILVA_MERDS),tax_table(phyl_SILVA_MERDS),sample_data(merds_map_trt_bio_nit),phy_tree(phyl_SILVA_MERDS))
nrow(sample_data(SILVA_MERDS_data))
#56
ntaxa(SILVA_MERDS_data)
#11892
sum(taxa_sums(SILVA_MERDS_data))
#980588
mean(sample_sums(SILVA_MERDS_data))
#17510.5
min(sample_sums(SILVA_MERDS_data))
#242
max(sample_sums(SILVA_MERDS_data))
#30633
sort(sample_sums(SILVA_MERDS_data))
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

#

######Begin NON-Rarefied community Analyses####
#let's look at raw ordination

SILVA_MERDS_ord=ordinate(SILVA_MERDS_trunc3, method = "NMDS",distance = "bray")
#*** Solution reached
#0.0275055
plot_ordination(SILVA_MERDS_trunc3,SILVA_MERDS_ord, color="root_association",shape="life_stage", label = "block")



alpha_meas = c("Observed", "Chao1", "Shannon", "InvSimpson")
SILVA_MERDS_trunc3_map=sample_data(SILVA_MERDS_trunc3)
SILVA_MERDS_trunc3.divfil=estimate_richness(SILVA_MERDS_trunc3,measures=alpha_meas)

SILVA_MERDS_trunc3.divfil=merge(SILVA_MERDS_trunc3.divfil, SILVA_MERDS_trunc3_map, by ="row.names")
#bact.soilE.t.divfil=mutate(bact.soilE.t.divfil, pielou=Shannon*(1/log(Observed)))
head(SILVA_MERDS_trunc3.divfil)
row.names(SILVA_MERDS_trunc3.divfil)=SILVA_MERDS_trunc3.divfil$Row.names
SILVA_MERDS_trunc3.divfil$Row.names=NULL
SILVA_MERDS_trunc3.divfil=merge(SILVA_MERDS_trunc3.divfil, sample_sums(SILVA_MERDS_trunc3), by ="row.names")
head(SILVA_MERDS_trunc3.divfil)
colnames(SILVA_MERDS_trunc3.divfil)[colnames(SILVA_MERDS_trunc3.divfil)=="y"]="read_abun"
row.names(SILVA_MERDS_trunc3.divfil)=SILVA_MERDS_trunc3.divfil$Row.names
SILVA_MERDS_trunc3.divfil$Row.names=NULL


ggplot(SILVA_MERDS_trunc3.divfil, aes(x=soil_status, y=read_abun))+geom_boxplot(aes(color=interaction(root_association,precip,life_stage)))+theme_bw()


SILVA_MERDS_trunc_sterile=subset_samples(SILVA_MERDS_trunc3, soil_status=="S")


SILVA_MERDS_sterile_ord=ordinate(SILVA_MERDS_trunc_sterile, method = "NMDS",distance = "bray")
#*** Solution reached
#0.09317059
plot_ordination(SILVA_MERDS_trunc_sterile,SILVA_MERDS_sterile_ord, color="precip",shape="life_stage", label = "block")

#see if there is a correlation between biomass and OTUs
SILVA_MERDS_trunc_sterile_core=core(SILVA_MERDS_trunc_sterile, detection = 0,prevalence = .75)
sample_sums(SILVA_MERDS_trunc_sterile_core)
SILVA_MERDS_trunc_sterile_core_map=sample_data(SILVA_MERDS_trunc_sterile_core)

biomas_mod_otus=lm_phyloseq(SILVA_MERDS_trunc_sterile_core, "total_biomass")
#Warning message:
#In transform(x, transformation) :
#  OTU table contains zeroes. Using log10(1 + x) transform.
"              logFC   AveExpr          t    P.Value adj.P.Val         B
OTU43     0.7509709 0.5854686  2.6450082 0.01668531 0.6877790 -4.394635
OTU78     0.7309627 0.7291038  2.4132150 0.02697173 0.6877790 -4.426389
OTU502    0.6450591 1.3997288  1.9770004 0.06394035 0.9490907 -4.485419
OTU1      0.2577273 2.4357622  1.3438871 0.19607232 0.9490907 -4.563083
OTU140    0.4008119 1.3817514  1.2711434 0.22023878 0.9490907 -4.570912
OTU2019   0.4418689 0.8933106  1.1892425 0.25015842 0.9490907 -4.579370
OTU148   -0.4028822 1.0668022 -1.1680958 0.25836483 0.9490907 -4.581489
OTU23918 -0.5095917 1.0461334 -1.1051741 0.28397941 0.9490907 -4.587629
OTU883   -0.3835969 2.3962545 -1.0213419 0.32094163 0.9490907 -4.595407
OTU752    0.1923123 1.6356522  0.9929961 0.33418281 0.9490907 -4.597928"
par(mfrow=c(2,2))
plot(get_sample(SILVA_MERDS_trunc_sterile_core,"OTU43"),SILVA_MERDS_trunc_sterile_core_map$total_biomass)
plot(get_sample(SILVA_MERDS_trunc_sterile_core,"OTU78"),SILVA_MERDS_trunc_sterile_core_map$total_biomass)
plot(get_sample(SILVA_MERDS_trunc_sterile_core,"OTU502"),SILVA_MERDS_trunc_sterile_core_map$total_biomass)
plot(get_sample(SILVA_MERDS_trunc_sterile_core,"OTU1"),SILVA_MERDS_trunc_sterile_core_map$total_biomass)
par(mfrow=c(2,2))
plot(get_sample(SILVA_MERDS_trunc_sterile_core,"OTU140"),SILVA_MERDS_trunc_sterile_core_map$total_biomass)
plot(get_sample(SILVA_MERDS_trunc_sterile_core,"OTU2019"),SILVA_MERDS_trunc_sterile_core_map$total_biomass)
plot(get_sample(SILVA_MERDS_trunc_sterile_core,"OTU148"),SILVA_MERDS_trunc_sterile_core_map$total_biomass)
plot(get_sample(SILVA_MERDS_trunc_sterile_core,"OTU23918"),SILVA_MERDS_trunc_sterile_core_map$total_biomass)
par(mfrow=c(2,1))
plot(get_sample(SILVA_MERDS_trunc_sterile_core,"OTU883"),SILVA_MERDS_trunc_sterile_core_map$total_biomass)
plot(get_sample(SILVA_MERDS_trunc_sterile_core,"OTU752"),SILVA_MERDS_trunc_sterile_core_map$total_biomass)


#How does this compare to the the OTU abundance versus biomass in live soil
SILVA_MERDS_trunc3_live_trt=subset_samples(SILVA_MERDS_trunc3, soil_status=="L"&life_stage!="Start")
nrow(sample_data(SILVA_MERDS_trunc3_live_trt))
#32
summary(sample_data(SILVA_MERDS_trunc3_live_trt))
SILVA_MERDS_trunc3_live_trt_map=sample_data(SILVA_MERDS_trunc3_live_trt)
par(mfrow=c(2,2))
plot(get_sample(SILVA_MERDS_trunc3_live_trt,"OTU43"),SILVA_MERDS_trunc3_live_trt_map$total_biomass)
plot(get_sample(SILVA_MERDS_trunc3_live_trt,"OTU78"),SILVA_MERDS_trunc3_live_trt_map$total_biomass)
plot(get_sample(SILVA_MERDS_trunc3_live_trt,"OTU502"),SILVA_MERDS_trunc3_live_trt_map$total_biomass)
plot(get_sample(SILVA_MERDS_trunc3_live_trt,"OTU1"),SILVA_MERDS_trunc3_live_trt_map$total_biomass)
par(mfrow=c(2,2))
plot(get_sample(SILVA_MERDS_trunc3_live_trt,"OTU2019"),SILVA_MERDS_trunc3_live_trt_map$total_biomass)
plot(get_sample(SILVA_MERDS_trunc3_live_trt,"OTU148"),SILVA_MERDS_trunc3_live_trt_map$total_biomass)
plot(get_sample(SILVA_MERDS_trunc3_live_trt,"OTU23918"),SILVA_MERDS_trunc3_live_trt_map$total_biomass)
plot(get_sample(SILVA_MERDS_trunc3_live_trt,"OTU883"),SILVA_MERDS_trunc3_live_trt_map$total_biomass)
par(mfrow=c(2,1))
plot(get_sample(SILVA_MERDS_trunc3_live_trt,"OTU752"),SILVA_MERDS_trunc3_live_trt_map$total_biomass)
plot(get_sample(SILVA_MERDS_trunc3_live_trt,"OTU24472"),SILVA_MERDS_trunc3_live_trt_map$total_biomass)


#how frequent are these OTUs in the rhizo and bulk starting soil

SILVA_MERDS_trunc3_start=subset_samples(SILVA_MERDS_trunc3, life_stage=="Start")
biomass_indic=c("OTU43","OTU78","OTU502","OTU1","OTU2019","OTU148","OTU23918","OTU883","OTU752","OTU24472")
SILVA_MERDS_trunc3_start_biomss_otu=t(otu_table(prune_taxa(biomass_indic,SILVA_MERDS_trunc3_start)))

SILVA_MERDS_trunc3_start_biomss_otu_trt=merge(SILVA_MERDS_trunc3_start_biomss_otu,sample_data(SILVA_MERDS_trunc3_start),by="row.names")

grid.arrange(ggplot(SILVA_MERDS_trunc3_start_biomss_otu_trt,aes(x=root_association,y=OTU43))+geom_boxplot(),
ggplot(SILVA_MERDS_trunc3_start_biomss_otu_trt,aes(x=root_association,y=OTU78))+geom_boxplot(),
ggplot(SILVA_MERDS_trunc3_start_biomss_otu_trt,aes(x=root_association,y=OTU502))+geom_boxplot(),
ggplot(SILVA_MERDS_trunc3_start_biomss_otu_trt,aes(x=root_association,y=OTU1))+geom_boxplot(), nrow=2,ncol=2)

grid.arrange(ggplot(SILVA_MERDS_trunc3_start_biomss_otu_trt,aes(x=root_association,y=OTU2019))+geom_boxplot(),
ggplot(SILVA_MERDS_trunc3_start_biomss_otu_trt,aes(x=root_association,y=OTU148))+geom_boxplot(),
ggplot(SILVA_MERDS_trunc3_start_biomss_otu_trt,aes(x=root_association,y=OTU23918))+geom_boxplot(),
ggplot(SILVA_MERDS_trunc3_start_biomss_otu_trt,aes(x=root_association,y=OTU883))+geom_boxplot(), nrow=2,ncol=2)


grid.arrange(ggplot(SILVA_MERDS_trunc3_start_biomss_otu_trt,aes(x=root_association,y=OTU752))+geom_boxplot(),
ggplot(SILVA_MERDS_trunc3_start_biomss_otu_trt,aes(x=root_association,y=OTU24472))+geom_boxplot(), nrow=2)

#there is a large outlier in both datasets that maybe skewing the results

#see if there is a correlation between biomass and OTUs
SILVA_MERDS_trunc_sterile_out=subset_samples(SILVA_MERDS_trunc_sterile, total_biomass<1)
SILVA_MERDS_trunc_sterile_out_core=core(SILVA_MERDS_trunc_sterile_out, detection = 0,prevalence = .75)
sample_sums(SILVA_MERDS_trunc_sterile_out_core)
SILVA_MERDS_trunc_sterile_out_core_map=sample_data(SILVA_MERDS_trunc_sterile_out_core)

biomas_mod_otus_otu=lm_phyloseq(SILVA_MERDS_trunc_sterile_out_core, "total_biomass")
#Warning message:
#In transform(x, transformation) :
#  OTU table contains zeroes. Using log10(1 + x) transform.
"             logFC   AveExpr         t    P.Value adj.P.Val         B
OTU148  -1.5066849 1.0589958 -2.689344 0.01563565 0.7661469 -4.533616
OTU2019  1.4267018 0.8926121  2.254355 0.03783274 0.7811411 -4.550973
OTU115  -2.1834756 1.6589504 -2.089652 0.05218602 0.7811411 -4.557522
OTU139  -1.6435907 1.2320660 -1.894277 0.07554546 0.7811411 -4.565162
OTU502   1.0398543 1.3559353  1.727591 0.10240739 0.7811411 -4.571497
OTU63   -1.2931896 2.1090492 -1.596441 0.12903756 0.7811411 -4.576314
OTU900   0.6975020 1.2261122  1.470375 0.15995043 0.7811411 -4.580770
OTU106   1.0074927 2.2038149  1.405929 0.17798032 0.7811411 -4.582970
OTU4     1.2710933 1.9247360  1.399929 0.17973994 0.7811411 -4.583172
OTU2839  0.8381175 0.9993585  1.325547 0.20274307 0.7811411 -4.585633"
par(mfrow=c(2,2))
plot(get_sample(SILVA_MERDS_trunc_sterile_out_core,"OTU148"),SILVA_MERDS_trunc_sterile_out_core_map$total_biomass)
plot(get_sample(SILVA_MERDS_trunc_sterile_out_core,"OTU2019"),SILVA_MERDS_trunc_sterile_out_core_map$total_biomass)
plot(get_sample(SILVA_MERDS_trunc_sterile_out_core,"OTU115"),SILVA_MERDS_trunc_sterile_out_core_map$total_biomass)
plot(get_sample(SILVA_MERDS_trunc_sterile_out_core,"OTU139"),SILVA_MERDS_trunc_sterile_out_core_map$total_biomass)
par(mfrow=c(2,2))
plot(get_sample(SILVA_MERDS_trunc_sterile_out_core,"OTU502"),SILVA_MERDS_trunc_sterile_out_core_map$total_biomass)
plot(get_sample(SILVA_MERDS_trunc_sterile_out_core,"OTU63"),SILVA_MERDS_trunc_sterile_out_core_map$total_biomass)
plot(get_sample(SILVA_MERDS_trunc_sterile_out_core,"OTU900"),SILVA_MERDS_trunc_sterile_out_core_map$total_biomass)
plot(get_sample(SILVA_MERDS_trunc_sterile_out_core,"OTU106"),SILVA_MERDS_trunc_sterile_out_core_map$total_biomass)
par(mfrow=c(2,1))
plot(get_sample(SILVA_MERDS_trunc_sterile_out_core,"OTU4"),SILVA_MERDS_trunc_sterile_out_core_map$total_biomass)
plot(get_sample(SILVA_MERDS_trunc_sterile_out_core,"OTU2839"),SILVA_MERDS_trunc_sterile_out_core_map$total_biomass)


#How does this compare to the the OTU abundance versus biomass in live soil
SILVA_MERDS_trunc3_live_trt=subset_samples(SILVA_MERDS_trunc3, soil_status=="L"&life_stage!="Start")
SILVA_MERDS_trunc3_live_trt_out=subset_samples(SILVA_MERDS_trunc3_live_trt, total_biomass<1)
nrow(sample_data(SILVA_MERDS_trunc3_live_trt_out))
#31
summary(sample_data(SILVA_MERDS_trunc3_live_trt_out))
SILVA_MERDS_trunc3_live_trt_out_map=sample_data(SILVA_MERDS_trunc3_live_trt_out)
par(mfrow=c(2,2))
plot(get_sample(SILVA_MERDS_trunc3_live_trt_out,"OTU148"),SILVA_MERDS_trunc3_live_trt_out_map$total_biomass)
plot(get_sample(SILVA_MERDS_trunc3_live_trt_out,"OTU2019"),SILVA_MERDS_trunc3_live_trt_out_map$total_biomass)
plot(get_sample(SILVA_MERDS_trunc3_live_trt_out,"OTU115"),SILVA_MERDS_trunc3_live_trt_out_map$total_biomass)
plot(get_sample(SILVA_MERDS_trunc3_live_trt_out,"OTU139"),SILVA_MERDS_trunc3_live_trt_out_map$total_biomass)
par(mfrow=c(2,2))
plot(get_sample(SILVA_MERDS_trunc3_live_trt_out,"OTU502"),SILVA_MERDS_trunc3_live_trt_out_map$total_biomass)
plot(get_sample(SILVA_MERDS_trunc3_live_trt_out,"OTU63"),SILVA_MERDS_trunc3_live_trt_out_map$total_biomass)
plot(get_sample(SILVA_MERDS_trunc3_live_trt_out,"OTU900"),SILVA_MERDS_trunc3_live_trt_out_map$total_biomass)
plot(get_sample(SILVA_MERDS_trunc3_live_trt_out,"OTU106"),SILVA_MERDS_trunc3_live_trt_out_map$total_biomass)
par(mfrow=c(2,1))
plot(get_sample(SILVA_MERDS_trunc3_live_trt_out,"OTU4"),SILVA_MERDS_trunc3_live_trt_out_map$total_biomass)
plot(get_sample(SILVA_MERDS_trunc3_live_trt_out,"OTU2839"),SILVA_MERDS_trunc3_live_trt_out_map$total_biomass)




#Let look at transplants and see if there is a correlation between biomass and OTUs
sample_data(SILVA_MERDS_trunc_sterile)
SILVA_MERDS_trunc_sterile_trans=subset_samples(SILVA_MERDS_trunc_sterile, life_stage=="G")
sample_data(SILVA_MERDS_trunc_sterile_trans)

#Only taxa in 75% of samples
SILVA_MERDS_trunc_sterile_trans_core=core(SILVA_MERDS_trunc_sterile_trans, detection = 0,prevalence = .75)
sample_sums(SILVA_MERDS_trunc_sterile_trans_core)


biomas_mod_otus=lm_phyloseq(SILVA_MERDS_trunc_sterile_trans_core, "total_biomass")
#Warning message:
#In transform(x, transformation) :
#  OTU table contains zeroes. Using log10(1 + x) transform.
"              logFC   AveExpr         t     P.Value adj.P.Val         B
OTU883   -1.0760738 2.3961867 -3.759898 0.002902171 0.1704741 -1.440618
OTU43     1.0106521 0.7150720  3.524517 0.004427898 0.1704741 -1.797326
OTU78     0.8413813 0.7579017  2.903998 0.013698904 0.3036855 -2.757026
OTU12410  0.7028079 0.6318266  2.717903 0.019242181 0.3036855 -3.045739
OTU57    -0.7719324 2.2764640 -2.632344 0.022486918 0.3036855 -3.177849
OTU193    1.4985374 1.5014569  2.604288 0.023663804 0.3036855 -3.221041
OTU552    0.8886797 1.3904565  2.336573 0.038380859 0.3681663 -3.628536
OTU377   -0.6054359 1.8774871 -2.255027 0.044403555 0.3681663 -3.750422
OTU2833  -0.6852767 1.2696164 -2.219543 0.047297452 0.3681663 -3.803045
OTU242   -0.5914900 0.6902468 -2.213430 0.047813799 0.3681663 -3.812083"

SILVA_MERDS_trunc_sterile_trans_core_map=sample_data(SILVA_MERDS_trunc_sterile_trans_core)
SILVA_MERDS_trunc_sterile_trans_core_map$total_biomass

plot(get_sample(SILVA_MERDS_rar_sterile_trans,"OTU883"),SILVA_MERDS_trunc_sterile_trans_core_map$total_biomass)
plot(get_sample(SILVA_MERDS_rar_sterile_trans,"OTU43"),SILVA_MERDS_trunc_sterile_trans_core_map$total_biomass)
plot(get_sample(SILVA_MERDS_rar_sterile_trans,"OTU78"),SILVA_MERDS_trunc_sterile_trans_core_map$total_biomass)
plot(get_sample(SILVA_MERDS_rar_sterile_trans,"OTU12410"),SILVA_MERDS_trunc_sterile_trans_core_map$total_biomass)
plot(get_sample(SILVA_MERDS_rar_sterile_trans,"OTU57"),SILVA_MERDS_trunc_sterile_trans_core_map$total_biomass)
plot(get_sample(SILVA_MERDS_rar_sterile_trans,"OTU193"),SILVA_MERDS_trunc_sterile_trans_core_map$total_biomass)
plot(get_sample(SILVA_MERDS_rar_sterile_trans,"OTU552"),SILVA_MERDS_trunc_sterile_trans_core_map$total_biomass)
plot(get_sample(SILVA_MERDS_rar_sterile_trans,"OTU377"),SILVA_MERDS_trunc_sterile_trans_core_map$total_biomass)
plot(get_sample(SILVA_MERDS_rar_sterile_trans,"OTU2833"),SILVA_MERDS_trunc_sterile_trans_core_map$total_biomass)
plot(get_sample(SILVA_MERDS_rar_sterile_trans,"OTU242"),SILVA_MERDS_trunc_sterile_trans_core_map$total_biomass)

#####Trans Origin DESeq2####

unique(sample_data(SILVA_MERDS_trunc3)$soil_root)
unrar_SILVA_MERDS_rar_orig_trans_bulk=subset_samples(SILVA_MERDS_trunc3, life_stage=="G"&life_stage!="Start"&
                                                 soil_root=="L.B")
unrar_SILVA_MERDS_rar_orig_trans_bulk=prune_taxa(taxa_sums(unrar_SILVA_MERDS_rar_orig_trans_bulk) > 0, unrar_SILVA_MERDS_rar_orig_trans_bulk)
nsamples(unrar_SILVA_MERDS_rar_orig_trans_bulk)
#8
unique(sample_data(SILVA_MERDS_trunc3)$precip)

#orig_trans_bulk_ds2 <- phyloseq_to_deseq2(unrar_SILVA_MERDS_rar_orig_trans_bulk,design =  ~ precip)
unrar_SILVA_MERDS_rar_orig_trans_bulk_OTU=data.frame(otu_table(unrar_SILVA_MERDS_rar_orig_trans_bulk))
unrar_SILVA_MERDS_rar_orig_trans_bulk_OTU[1:8,1:8]

unrar_SILVA_MERDS_rar_orig_trans_bulk_precip=data.frame(sample_data(unrar_SILVA_MERDS_rar_orig_trans_bulk)[,c("precip","soil_root")])
summary(unrar_SILVA_MERDS_rar_orig_trans_bulk_precip)
unrar_SILVA_MERDS_rar_orig_trans_bulk_precip$precip=factor(unrar_SILVA_MERDS_rar_orig_trans_bulk_precip$precip,levels = c("A","D"))

orig_trans_bulk_dds <- DESeqDataSetFromMatrix(countData = unrar_SILVA_MERDS_rar_orig_trans_bulk_OTU,
                              colData = unrar_SILVA_MERDS_rar_orig_trans_bulk_precip,
                              design = ~ precip)

orig_trans_bulk_dds <- DESeq(orig_trans_bulk_dds)
orig_trans_bulk_res <- results(orig_trans_bulk_dds)
orig_trans_bulk_df <- as.data.frame(orig_trans_bulk_res)
orig_trans_bulk_df$taxon <- rownames(orig_trans_bulk_df)
orig_trans_bulk_df <- orig_trans_bulk_df %>% arrange(log2FoldChange, padj)

library(knitr)
print(head(kable((orig_trans_bulk_df))))


#Code from Jennifer Jones
orig_trans_bulk_dds <- DESeqDataSetFromMatrix(countData = unrar_SILVA_MERDS_rar_orig_trans_bulk_OTU,
                                              colData = unrar_SILVA_MERDS_rar_orig_trans_bulk_precip,
                                              design = ~ precip)

# calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(orig_trans_bulk_dds), 1, gm_mean)
diagdds = estimateSizeFactors(orig_trans_bulk_dds, geoMeans = geoMeans)
diagdds = DESeq(diagdds, fitType="parametric")
#Note: The default multiple-inference correction is Benjamini-Hochberg, and occurs within the DESeq function.

#I should probably play around with differenttypes of fit and tests on this
#diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

diagdds

res = results(diagdds, cooksCutoff = FALSE)
res = res[order(res$padj), ]
res
alpha = 0.01 
sigtab = res[which(res$pvalue < alpha), ] # screening p value by alpha
dim(sigtab)
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(unrar_SILVA_MERDS_rar_orig_trans_bulk)[rownames(sigtab), ], "matrix"))
#if you get an error here, it's probably because there aren't any OTUs that remain after screening based on alpha. 
dim(sigtab)
head(sigtab)
names(sigtab)

sigtab$Phylum


SILVA_MERDS_rar_orig_trans_rhizo=subset_samples(SILVA_MERDS_rar, life_stage=="G"&life_stage!="Start"&
                                                  soil_root=="L.R")
SILVA_MERDS_rar_orig_trans_rhizo=prune_taxa(taxa_sums(SILVA_MERDS_rar_orig_trans_rhizo) > 0, SILVA_MERDS_rar_orig_trans_rhizo)
nsamples(SILVA_MERDS_rar_orig_trans_rhizo)
#8

orig_trans_rhizo_ds2 <- phyloseq_to_deseq2(SILVA_MERDS_rar_orig_trans_rhizo, ~ precip)
orig_trans_rhizo_dds <- DESeq(orig_trans_rhizo_ds2)
orig_trans_rhizo_res <- results(orig_trans_rhizo_dds)
orig_trans_rhizo_df <- as.data.frame(orig_trans_rhizo_res)
orig_trans_rhizo_df$taxon <- rownames(orig_trans_rhizo_df)
orig_trans_rhizo_df <- orig_trans_rhizo_df %>% arrange(log2FoldChange, padj)

library(knitr)
print(head(kable((orig_trans_rhizo_df))))

#Code from Jennifer Jones
orig_trans_bulk_dds <- DESeqDataSetFromMatrix(countData = unrar_SILVA_MERDS_rar_orig_trans_bulk_OTU,
                                              colData = unrar_SILVA_MERDS_rar_orig_trans_bulk_precip,
                                              design = ~ precip)

# calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(orig_trans_bulk_dds), 1, gm_mean)
diagdds = estimateSizeFactors(orig_trans_bulk_dds, geoMeans = geoMeans)
diagdds = DESeq(diagdds, fitType="parametric")
#Note: The default multiple-inference correction is Benjamini-Hochberg, and occurs within the DESeq function.

#I should probably play around with differenttypes of fit and tests on this
#diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

diagdds

res = results(diagdds, cooksCutoff = FALSE)
res = res[order(res$padj), ]
res
alpha = 0.01 
sigtab = res[which(res$padj < alpha), ] # screening p value by alpha
dim(sigtab)
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(unrar_SILVA_MERDS_rar_orig_trans_bulk)[rownames(sigtab), ], "matrix"))
#if you get an error here, it's probably because there aren't any OTUs that remain after screening based on alpha. 
dim(sigtab)
head(sigtab)
names(sigtab)


#####Trans Origin DESeq2####

unique(sample_data(SILVA_MERDS_trunc3)$soil_root)
unrar_SILVA_MERDS_rar_orig_seed_bulk=subset_samples(SILVA_MERDS_trunc3, life_stage=="S"&life_stage!="Start"&
                                                       soil_root=="L.B")
unrar_SILVA_MERDS_rar_orig_seed_bulk=prune_taxa(taxa_sums(unrar_SILVA_MERDS_rar_orig_seed_bulk) > 10, unrar_SILVA_MERDS_rar_orig_seed_bulk)
nsamples(unrar_SILVA_MERDS_rar_orig_seed_bulk)
#8
unique(sample_data(SILVA_MERDS_trunc3)$precip)


#####
#orig_trans_bulk_ds2 <- phyloseq_to_deseq2(unrar_SILVA_MERDS_rar_orig_trans_bulk,design =  ~ precip)
unrar_SILVA_MERDS_rar_orig_trans_bulk_OTU=data.frame(otu_table(unrar_SILVA_MERDS_rar_orig_trans_bulk))
unrar_SILVA_MERDS_rar_orig_trans_bulk_OTU[1:8,1:8]

unrar_SILVA_MERDS_rar_orig_trans_bulk_precip=data.frame(sample_data(unrar_SILVA_MERDS_rar_orig_trans_bulk)[,c("precip","soil_root")])
summary(unrar_SILVA_MERDS_rar_orig_trans_bulk_precip)
unrar_SILVA_MERDS_rar_orig_trans_bulk_precip$precip=factor(unrar_SILVA_MERDS_rar_orig_trans_bulk_precip$precip,levels = c("A","D"))

orig_trans_bulk_dds <- DESeqDataSetFromMatrix(countData = unrar_SILVA_MERDS_rar_orig_trans_bulk_OTU,
                                              colData = unrar_SILVA_MERDS_rar_orig_trans_bulk_precip,
                                              design = ~ precip)
#####

orig_seed_rhizo_ds2 <- phyloseq_to_deseq2(unrar_SILVA_MERDS_rar_orig_seed_bulk, ~ precip)
orig_seed_bulk_dds <- DESeq(orig_seed_rhizo_ds2)
orig_seed_bulk_res <- results(orig_seed_bulk_dds)
orig_seed_bulk_df <- as.data.frame(orig_seed_bulk_res)
orig_seed_bulk_df$taxon <- rownames(orig_seeds_bulk_df)
orig_seed_bulk_df <- orig_seeds_bulk_df %>% arrange(log2FoldChange, padj)

library(knitr)
print(head(kable((orig_trans_bulk_df))))


#Code from Jennifer Jones


# calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(orig_seed_rhizo_ds2), 1, gm_mean)
diagdds = estimateSizeFactors(orig_seed_rhizo_ds2, geoMeans = geoMeans)
diagdds = DESeq(diagdds, fitType="parametric")
#Note: The default multiple-inference correction is Benjamini-Hochberg, and occurs within the DESeq function.

#I should probably play around with differenttypes of fit and tests on this
#diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

diagdds

res = results(diagdds, cooksCutoff = FALSE)
res = res[order(res$padj), ]
res
alpha = 0.01 
sigtab = res[which(res$pvalue < alpha), ] # screening p value by alpha
dim(sigtab)
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(unrar_SILVA_MERDS_rar_orig_seed_bulk)[rownames(sigtab), ], "matrix"))
#if you get an error here, it's probably because there aren't any OTUs that remain after screening based on alpha. 
dim(sigtab)
head(sigtab)
names(sigtab)

sigtab$Phylum


SILVA_MERDS_rar_orig_trans_rhizo=subset_samples(SILVA_MERDS_rar, life_stage=="G"&life_stage!="Start"&
                                                  soil_root=="L.R")
SILVA_MERDS_rar_orig_trans_rhizo=prune_taxa(taxa_sums(SILVA_MERDS_rar_orig_trans_rhizo) > 0, SILVA_MERDS_rar_orig_trans_rhizo)
nsamples(SILVA_MERDS_rar_orig_trans_rhizo)
#8

orig_trans_rhizo_ds2 <- phyloseq_to_deseq2(SILVA_MERDS_rar_orig_trans_rhizo, ~ precip)
orig_trans_rhizo_dds <- DESeq(orig_trans_rhizo_ds2)
orig_trans_rhizo_res <- results(orig_trans_rhizo_dds)
orig_trans_rhizo_df <- as.data.frame(orig_trans_rhizo_res)
orig_trans_rhizo_df$taxon <- rownames(orig_trans_rhizo_df)
orig_trans_rhizo_df <- orig_trans_rhizo_df %>% arrange(log2FoldChange, padj)

library(knitr)
print(head(kable((orig_trans_rhizo_df))))

#Code from Jennifer Jones
orig_trans_bulk_dds <- DESeqDataSetFromMatrix(countData = unrar_SILVA_MERDS_rar_orig_trans_bulk_OTU,
                                              colData = unrar_SILVA_MERDS_rar_orig_trans_bulk_precip,
                                              design = ~ precip)

# calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(orig_trans_bulk_dds), 1, gm_mean)
diagdds = estimateSizeFactors(orig_trans_bulk_dds, geoMeans = geoMeans)
diagdds = DESeq(diagdds, fitType="parametric")
#Note: The default multiple-inference correction is Benjamini-Hochberg, and occurs within the DESeq function.

#I should probably play around with differenttypes of fit and tests on this
#diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

diagdds

res = results(diagdds, cooksCutoff = FALSE)
res = res[order(res$padj), ]
res
alpha = 0.01 
sigtab = res[which(res$padj < alpha), ] # screening p value by alpha
dim(sigtab)
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(unrar_SILVA_MERDS_rar_orig_trans_bulk)[rownames(sigtab), ], "matrix"))
#if you get an error here, it's probably because there aren't any OTUs that remain after screening based on alpha. 
dim(sigtab)
head(sigtab)
names(sigtab)





######END NON-Rarefied community Analyses####

######Rarefied community Analyses####

#Let's rarefy the data
#SILVA_MERDS_rar=rarefy_even_depth(SILVA_MERDS_trunc3, sample.size= 10000, rngseed = T)
#498OTUs were removed because they are no longer 
#present in any sample after random subsampling

#save(SILVA_MERDS_rar, file = "D:/MERDS_2018/merds/Switchgrass/R_data/SILVA_MERDS_rar_phylo_obj.RData")
load("D:/MERDS_2018/merds/Switchgrass/R_data/SILVA_MERDS_rar_phylo_obj.RData")
#SILVA_MERDS_rar=phyloseq(otu_table(SILVA_MERDS_rar),tax_table(SILVA_MERDS_rar),sample_data(merds_map_trt_bio_nit),phy_tree(run_20190617_16S.tree))
#phy_tree(SILVA_MERDS_rar)<-ape::root(phy_tree(SILVA_MERDS_rar), "OTU3039", resolve.root=TRUE)

SILVA_MERDS_rar_ord=ordinate(SILVA_MERDS_rar, method = "NMDS",distance = "bray")
#*** Solution reached
#Warning message:
#  In metaMDS(veganifyOTU(physeq), distance, ...) :
#  stress is (nearly) zero: you may have insufficient data
#7.6596e-05 
plot_ordination(SILVA_MERDS_rar,SILVA_MERDS_rar_ord, color="root_association",shape="life_stage")+geom_point(size=3)+
  theme_bw()



SILVA_MERDS_rar_map=sample_data(SILVA_MERDS_rar)
SILVA_MERDS_rar_map$soil_root_stage=with(SILVA_MERDS_rar_map, interaction(soil_status,root_association,life_stage))
SILVA_MERDS_rar_dis=distance(SILVA_MERDS_rar,method = "bray")

adonis(SILVA_MERDS_rar_dis~SILVA_MERDS_rar_map$soil_root_stage+as.factor(SILVA_MERDS_rar_map$block), permutations = 9999)
#                                     Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#SILVA_MERDS_rar_map$soil_root_stage   7    7.2229 1.03184  7.4654 0.51823 0.0001 ***
#as.factor(SILVA_MERDS_rar_map$block)  3    0.6331 0.21103  1.5268 0.04542 0.0639 .  
#Residuals                            44    6.0816 0.13822         0.43634           
#Total                                54   13.9376                 1.00000      

pairwise.perm.manova(SILVA_MERDS_rar_dis, SILVA_MERDS_rar_map$soil_root_stage, nperm=2000)



"	Pairwise comparisons using permutation MANOVAs on a distance matrix 

data:  SILVA_MERDS_rar_dis by SILVA_MERDS_rar_map$soil_root_stage
2000 permutations 

          L.B.G  S.B.G  L.R.G  L.B.S  S.B.S  L.R.S  L.B.Start
S.B.G     0.0020 -      -      -      -      -      -        
L.R.G     0.6613 0.0020 -      -      -      -      -        
L.B.S     0.0950 0.0020 0.1181 -      -      -      -        
S.B.S     0.0023 0.0020 0.0023 0.0020 -      -      -        
L.R.S     0.0614 0.0023 0.2298 0.0020 0.0023 -      -        
L.B.Start 0.0037 0.0037 0.0037 0.0023 0.0053 0.0037 -        
L.R.Start 0.0037 0.0037 0.0042 0.0020 0.0070 0.0037 0.9425   

P value adjustment method: fdr "


#####Jaccard####

SILVA_MERDS_rar_J_ord=ordinate(SILVA_MERDS_rar, method = "NMDS", distance = "jaccard", binary = TRUE)
#*** Solution reached
#Warning message:
#  In metaMDS(veganifyOTU(physeq), distance, ...) :
#  stress is (nearly) zero: you may have insufficient data
#7.099987e-05
plot_ordination(SILVA_MERDS_rar,SILVA_MERDS_rar_J_ord, color="root_association",shape="life_stage")+geom_point(size=3)+
  theme_bw()



SILVA_MERDS_rar_map=sample_data(SILVA_MERDS_rar)
SILVA_MERDS_rar_map$soil_root_stage=with(SILVA_MERDS_rar_map, interaction(soil_status,root_association,life_stage))
SILVA_MERDS_rar_JAC_dis=distance(SILVA_MERDS_rar,method = "jaccard", binary = TRUE)

adonis(SILVA_MERDS_rar_JAC_dis~SILVA_MERDS_rar_map$soil_root_stage+as.factor(SILVA_MERDS_rar_map$block), permutations = 9999)
#                                     Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#SILVA_MERDS_rar_map$soil_root_stage   7    6.5457 0.93509  3.8726 0.36094 0.0001 ***
#as.factor(SILVA_MERDS_rar_map$block)  3    0.9648 0.32160  1.3319 0.05320 0.0605 .  
#Residuals                            44   10.6244 0.24146         0.58586           
#Total                                54   18.1349                 1.00000       

pairwise.perm.manova(SILVA_MERDS_rar_JAC_dis, SILVA_MERDS_rar_map$soil_root_stage, nperm=2000)

"data:  SILVA_MERDS_rar_JAC_dis by SILVA_MERDS_rar_map$soil_root_stage
2000 permutations 

          L.B.G  S.B.G  L.R.G  L.B.S  S.B.S  L.R.S  L.B.Start
S.B.G     0.0023 -      -      -      -      -      -        
L.R.G     0.8313 0.0023 -      -      -      -      -        
L.B.S     0.1697 0.0030 0.3750 -      -      -      -        
S.B.S     0.0023 0.0023 0.0030 0.0023 -      -      -        
L.R.S     0.1697 0.0023 0.5124 0.0172 0.0030 -      -        
L.B.Start 0.0030 0.0030 0.0077 0.0030 0.0059 0.0049 -        
L.R.Start 0.0044 0.0030 0.0100 0.0030 0.0054 0.0037 0.8486   

P value adjustment method: fdr "


#####Weighted Unifrac####

SILVA_MERDS_rar_WU_ord=ordinate(SILVA_MERDS_rar, method = "NMDS", distance = "wunifrac")
#*** Solution reached
#0.0769077 
plot_ordination(SILVA_MERDS_rar,SILVA_MERDS_rar_WU_ord, color="root_association",shape="life_stage")+geom_point(size=3)+
  theme_bw()



SILVA_MERDS_rar_map=sample_data(SILVA_MERDS_rar)
SILVA_MERDS_rar_map$soil_root_stage=with(SILVA_MERDS_rar_map, interaction(soil_status,root_association,life_stage))
SILVA_MERDS_rar_WU_dis=distance(SILVA_MERDS_rar,method = "wunifrac")

adonis(SILVA_MERDS_rar_WU_dis~SILVA_MERDS_rar_map$soil_root_stage+as.factor(SILVA_MERDS_rar_map$block), permutations = 9999)
#                                     Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)    
#SILVA_MERDS_rar_map$soil_root_stage   7   1.30035 0.185765  13.367 0.66753 0.0001 ***
#as.factor(SILVA_MERDS_rar_map$block)  3   0.03619 0.012062   0.868 0.01858 0.5280    
#Residuals                            44   0.61146 0.013897         0.31389           
#Total                                54   1.94800                  1.00000        

pairwise.perm.manova(SILVA_MERDS_rar_WU_dis, SILVA_MERDS_rar_map$soil_root_stage, nperm=2000)

"data:  SILVA_MERDS_rar_WU_dis by SILVA_MERDS_rar_map$soil_root_stage
2000 permutations 

          L.B.G  S.B.G  L.R.G  L.B.S  S.B.S  L.R.S  L.B.Start
S.B.G     0.0016 -      -      -      -      -      -        
L.R.G     0.7206 0.0016 -      -      -      -      -        
L.B.S     0.0468 0.0016 0.1954 -      -      -      -        
S.B.S     0.0016 0.0054 0.0016 0.0016 -      -      -        
L.R.S     0.0952 0.0016 0.0952 0.0016 0.0016 -      -        
L.B.Start 0.0049 0.0049 0.0087 0.0047 0.0063 0.0028 -        
L.R.Start 0.0047 0.0049 0.0134 0.0049 0.0063 0.0049 0.6250   

P value adjustment method: fdr "



#####Betadisp####
SILVA_MERDS_rar_map=sample_data(SILVA_MERDS_rar)
head(SILVA_MERDS_rar_map)

SILVA_MERDS_rar_map$soil_root_stage=with(SILVA_MERDS_rar_map, interaction(soil_status,root_association,life_stage))
unique(SILVA_MERDS_rar_map$soil_root_stage)
unique(SILVA_MERDS_rar_map$precip)
SILVA_MERDS_rar_WU_dis=distance(SILVA_MERDS_rar,method = "wunifrac")
SILVA_MERDS_rar_WU_dis_betamod=betadisper(SILVA_MERDS_rar_WU_dis, 
                                                               with(SILVA_MERDS_rar_map,interaction(soil_root_stage,precip)))


mean(SILVA_MERDS_rar_WU_dis_betamod$distances)
#0.08777776
sd(SILVA_MERDS_rar_WU_dis_betamod$distances)
#0.0317882
sd(SILVA_MERDS_rar_WU_dis_betamod$distances)/sqrt(length(SILVA_MERDS_rar_WU_dis_betamod$distances))
#0.00428632


head(SILVA_MERDS_rar_WU_dis_betamod)
SILVA_MERDS_rar_WU_dis_betamod_grp=data.frame(SILVA_MERDS_rar_WU_dis_betamod$distances,as.factor(SILVA_MERDS_rar_WU_dis_betamod$group))

colnames(SILVA_MERDS_rar_WU_dis_betamod_grp)=c("W_unifrac_betadisp","betagroup")
nrow(SILVA_MERDS_rar_WU_dis_betamod_grp)
#55

SILVA_MERDS_rar_WU_dis_betamod_map=merge(SILVA_MERDS_rar_WU_dis_betamod_grp,SILVA_MERDS_rar_map,by="row.names")
nrow(SILVA_MERDS_rar_WU_dis_betamod_map)
#55
head(SILVA_MERDS_rar_WU_dis_betamod_map)
SILVA_MERDS_rar_WU_dis_betamod_map$time=ifelse(SILVA_MERDS_rar_WU_dis_betamod_map$precip=="Start","Start","End")

#Transplant

#####Experimental Betadisp Transplant####

SILVA_MERDS_rar_WU_dis_betamod_map_trans=subset(SILVA_MERDS_rar_WU_dis_betamod_map,life_stage == "G"|life_stage=="Start")

unique(SILVA_MERDS_rar_WU_dis_betamod_map_trans$betagroup)
SILVA_MERDS_rar_WU_dis_betamod_map_trans$W_unifrac_betadisp

wu_beta_inter_exp=SILVA_MERDS_rar_WU_dis_betamod_map_trans %>% group_by(betagroup) %>% summarise_at(vars(W_unifrac_betadisp),c(~mean(.),~n()))

ggplot(SILVA_MERDS_rar_WU_dis_betamod_map_trans,aes(x=factor(time, levels = c("Start","End")),y=W_unifrac_betadisp))+
  geom_point(aes(fill=precip, shape=soil_root), size=5)+scale_y_continuous(name = "Betadispersion \nWeighted Unifrac distance")+
  scale_shape_manual(values = c(21,23,24), name=NULL)+
  scale_fill_manual(values = c("blue","red","black"))+
  geom_segment(aes(x=1,xend=2,y=subset(wu_beta_inter_exp,betagroup=="L.B.Start.Start")$mean,
                   yend=subset(wu_beta_inter_exp,betagroup=="L.B.G.A")$mean), linetype="solid", color="blue",size=1.5)+
  geom_segment(aes(x=1,xend=2,y=subset(wu_beta_inter_exp,betagroup=="L.B.Start.Start")$mean,
                   yend=subset(wu_beta_inter_exp,betagroup=="L.B.G.D")$mean), linetype="solid", color="red",size=1.5)+
  geom_segment(aes(x=1,xend=2,y=subset(wu_beta_inter_exp,betagroup=="L.R.Start.Start")$mean,
                   yend=subset(wu_beta_inter_exp,betagroup=="L.R.G.A")$mean), linetype="dashed", color="blue",size=1.5)+
  geom_segment(aes(x=1,xend=2,y=subset(wu_beta_inter_exp,betagroup=="L.R.Start.Start")$mean,
                   yend=subset(wu_beta_inter_exp,betagroup=="L.R.G.D")$mean), linetype="dashed", color="red",size=1.5)+
  geom_segment(aes(x=1,xend=2,y=subset(wu_beta_inter_exp,betagroup=="L.B.Start.Start")$mean,
                   yend=subset(wu_beta_inter_exp,betagroup=="S.B.G.A")$mean), linetype="dotted", color="blue",size=1.5)+
  geom_segment(aes(x=1,xend=2,y=subset(wu_beta_inter_exp,betagroup=="L.R.Start.Start")$mean,
                   yend=subset(wu_beta_inter_exp,betagroup=="S.B.G.D")$mean), linetype="dotted", color="red",size=1.5)+
  theme_bw()+theme(axis.text = element_text(size = 18),axis.title.x = element_blank(), axis.title.y = element_text(size = 22),
                   legend.position = "none")
#9.5x7

ggplot(SILVA_MERDS_rar_WU_dis_betamod_map_trans,aes(x=factor(time, levels = c("Start","End")),y=W_unifrac_betadisp))+
  geom_point(aes(fill=precip, shape=soil_root), size=5)+scale_y_continuous(name = "Betadispersion \nWeighted Unifrac distance")+
  scale_shape_manual(values = c(21,23,24), name=NULL)+
  scale_fill_manual(values = c("blue","red","black"))+
  geom_segment(aes(x=1,xend=2,y=subset(wu_beta_inter_exp,betagroup=="L.B.Start.Start")$mean,
                   yend=subset(wu_beta_inter_exp,betagroup=="L.B.G.A")$mean), linetype="solid", color="blue",size=1.5)+
  geom_segment(aes(x=1,xend=2,y=subset(wu_beta_inter_exp,betagroup=="L.B.Start.Start")$mean,
                   yend=subset(wu_beta_inter_exp,betagroup=="L.B.G.D")$mean), linetype="solid", color="red",size=1.5)+
  geom_segment(aes(x=1,xend=2,y=subset(wu_beta_inter_exp,betagroup=="L.R.Start.Start")$mean,
                   yend=subset(wu_beta_inter_exp,betagroup=="L.R.G.A")$mean), linetype="dashed", color="blue",size=1.5)+
  geom_segment(aes(x=1,xend=2,y=subset(wu_beta_inter_exp,betagroup=="L.R.Start.Start")$mean,
                   yend=subset(wu_beta_inter_exp,betagroup=="L.R.G.D")$mean), linetype="dashed", color="red",size=1.5)+
  geom_segment(aes(x=1,xend=2,y=subset(wu_beta_inter_exp,betagroup=="L.B.Start.Start")$mean,
                   yend=subset(wu_beta_inter_exp,betagroup=="S.B.G.A")$mean), linetype="dotted", color="blue",size=1.5)+
  geom_segment(aes(x=1,xend=2,y=subset(wu_beta_inter_exp,betagroup=="L.R.Start.Start")$mean,
                   yend=subset(wu_beta_inter_exp,betagroup=="S.B.G.D")$mean), linetype="dotted", color="red",size=1.5)+
  theme_bw()+theme(axis.text = element_text(size = 18),axis.title.x = element_blank(), axis.title.y = element_text(size = 22),
                   legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#900x750


#####Presence Betadisp Transplant####
SILVA_MERDS_rar_WU_dis_betamod_map_trans=subset(SILVA_MERDS_rar_WU_dis_betamod_map,life_stage == "G"|life_stage=="Start")
SILVA_MERDS_rar_WU_dis_betamod_map_trans_B=subset(SILVA_MERDS_rar_WU_dis_betamod_map_trans,root_association=="B")

 
SILVA_MERDS_rar_WU_dis_betamod_map_trans_pres=subset(SILVA_MERDS_rar_WU_dis_betamod_map_trans_B,time!="Start")
nrow(SILVA_MERDS_rar_WU_dis_betamod_map_trans_pres)
#16

WU_beta_trans_mod_pres= lm(W_unifrac_betadisp~soil_root*precip, data= SILVA_MERDS_rar_WU_dis_betamod_map_trans_pres)
qqPlot(resid(WU_beta_trans_mod_pres))
hist(resid(WU_beta_trans_mod_pres))
shapiro.test(resid(WU_beta_trans_mod_pres))
#0.9686

Anova(WU_beta_trans_mod_pres, type=3)
#Nada

emmeans(WU_beta_trans_mod_pres, pairwise~s1_precip|s1_soil_root)


#####Origin Betadisp Transplant####


SILVA_MERDS_rar_WU_dis_betamod_map_trans=subset(SILVA_MERDS_rar_WU_dis_betamod_map,life_stage == "G"|life_stage=="Start")
SILVA_MERDS_rar_WU_dis_betamod_map_trans_L=subset(SILVA_MERDS_rar_WU_dis_betamod_map_trans,soil_status=="L")
unique(SILVA_MERDS_rar_WU_dis_betamod_map_trans_L$betagroup)
SILVA_MERDS_rar_WU_dis_betamod_map_trans_L$W_unifrac_betadisp

wu_beta_inter=SILVA_MERDS_rar_WU_dis_betamod_map_trans_L %>% group_by(betagroup) %>% summarise_at(vars(W_unifrac_betadisp),c(~mean(.),~n()))

ggplot(SILVA_MERDS_rar_WU_dis_betamod_map_trans_L,aes(x=factor(time, levels = c("Start","End")),y=W_unifrac_betadisp))+
  geom_point(aes(fill=precip, shape=soil_root), size=5)+scale_y_continuous(name = "Betadispersion \nWeighted Unifrac distance")+
  scale_shape_manual(values = c(21,24), name=NULL,labels=c("Ambient","Drought","start"))+
  scale_fill_manual(values = c("blue","red","black"),labels=c("Bulk","Rhizo"),)+
  geom_segment(aes(x=1,xend=2,y=subset(wu_beta_inter,betagroup=="L.B.Start.Start")$mean,
                   yend=subset(wu_beta_inter,betagroup=="L.B.G.A")$mean), linetype="solid", color="blue",size=1.5)+
  geom_segment(aes(x=1,xend=2,y=subset(wu_beta_inter,betagroup=="L.B.Start.Start")$mean,
                   yend=subset(wu_beta_inter,betagroup=="L.B.G.D")$mean), linetype="solid", color="red",size=1.5)+
  geom_segment(aes(x=1,xend=2,y=subset(wu_beta_inter,betagroup=="L.R.Start.Start")$mean,
                   yend=subset(wu_beta_inter,betagroup=="L.R.G.A")$mean), linetype="dashed", color="blue",size=1.5)+
  geom_segment(aes(x=1,xend=2,y=subset(wu_beta_inter,betagroup=="L.R.Start.Start")$mean,
                   yend=subset(wu_beta_inter,betagroup=="L.R.G.D")$mean), linetype="dashed", color="red",size=1.5)+
  theme_bw()+theme(axis.text = element_text(size = 18),axis.title.x = element_blank(), axis.title.y = element_text(size = 22),
                   legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#9.5x7
#+
#beta_disp_line_seedling_soil_origin
#geom_text(aes(x=2,y=W_unifrac_betadisp, label=Row.names))
subset(SILVA_MERDS_rar_WU_dis_betamod_map_trans_L,Row.names=="MERDSSG145")


SILVA_MERDS_rar_WU_dis_betamod_map_trans_ST=subset(SILVA_MERDS_rar_WU_dis_betamod_map_trans_L,time!="Start")
nrow(SILVA_MERDS_rar_WU_dis_betamod_map_trans_ST)
#16

WU_beta_trans_mod= lm(W_unifrac_betadisp~soil_root*precip, data= SILVA_MERDS_rar_WU_dis_betamod_map_trans_ST)
qqPlot(resid(WU_beta_trans_mod))
hist(resid(WU_beta_trans_mod))
shapiro.test(resid(WU_beta_trans_mod))
#0.1002

Anova(WU_beta_trans_mod, type=3)
#Nada

emmeans(WU_beta_trans_mod, pairwise~s1_precip|s1_soil_root)


#Let's remove the outlier community MERDSSG145

SILVA_MERDS_rar_WU_dis_betamod_map_trans_L_sub=subset(SILVA_MERDS_rar_WU_dis_betamod_map_trans_L,Row.names!="MERDSSG145")
nrow(SILVA_MERDS_rar_WU_dis_betamod_map_trans_L_sub)
#23

wu_beta_inter_sub=SILVA_MERDS_rar_WU_dis_betamod_map_trans_L_sub %>% group_by(betagroup) %>% summarise_at(vars(W_unifrac_betadisp),c(~mean(.),~n()))

ggplot(SILVA_MERDS_rar_WU_dis_betamod_map_trans_L_sub,aes(x=factor(time, levels = c("Start","End")),y=W_unifrac_betadisp))+
  geom_point(aes(color=soil_root, shape=precip), size=3)+scale_y_continuous(name = "Betadispersion \nWeighted Unifrac distance")+
  scale_shape_discrete(name=NULL,labels=c("Ambient","Drought","start"))+scale_colour_manual(values = c("grey","black"),labels=c("Bulk","Rhizo"),)+
  geom_segment(aes(x=1,xend=2,y=subset(wu_beta_inter_sub,betagroup=="L.B.Start.Start")$mean,
                   yend=subset(wu_beta_inter_sub,betagroup=="L.B.G.A")$mean), linetype="solid", color="grey",size=1.5)+
  geom_segment(aes(x=1,xend=2,y=subset(wu_beta_inter_sub,betagroup=="L.B.Start.Start")$mean,
                   yend=subset(wu_beta_inter_sub,betagroup=="L.B.G.D")$mean), linetype="dashed", color="grey",size=1.5)+
  geom_segment(aes(x=1,xend=2,y=subset(wu_beta_inter_sub,betagroup=="L.R.Start.Start")$mean,
                   yend=subset(wu_beta_inter_sub,betagroup=="L.R.G.A")$mean), linetype="solid", color="black",size=1.5)+
  geom_segment(aes(x=1,xend=2,y=subset(wu_beta_inter_sub,betagroup=="L.R.Start.Start")$mean,
                   yend=subset(wu_beta_inter_sub,betagroup=="L.R.G.D")$mean), linetype="dashed", color="black",size=1.5)+
  theme_bw()+theme(axis.text = element_text(size = 18),axis.title.x = element_blank(), axis.title.y = element_text(size = 22))
#+
#geom_text(aes(x=2,y=W_unifrac_betadisp, label=Row.names))



SILVA_MERDS_rar_WU_dis_betamod_map_trans_ST_s=subset(SILVA_MERDS_rar_WU_dis_betamod_map_trans_L_sub,time!="Start")
nrow(SILVA_MERDS_rar_WU_dis_betamod_map_trans_ST_s)
#15

WU_beta_trans_s_mod= lm(W_unifrac_betadisp~soil_root*precip, data= SILVA_MERDS_rar_WU_dis_betamod_map_trans_ST_s)
qqPlot(resid(WU_beta_trans_s_mod))
hist(resid(WU_beta_trans_s_mod))
shapiro.test(resid(WU_beta_trans_s_mod))
#0.6797

Anova(WU_beta_trans_s_mod, type=3)
#Nada

emmeans(WU_beta_trans_s_mod, pairwise~s1_precip|s1_soil_root)







#Seed

#####Experimental Betadisp Seed####



SILVA_MERDS_rar_WU_dis_betamod_map_seed=subset(SILVA_MERDS_rar_WU_dis_betamod_map,life_stage == "S"|life_stage=="Start")



unique(SILVA_MERDS_rar_WU_dis_betamod_map_seed$betagroup)
SILVA_MERDS_rar_WU_dis_betamod_map_seed$W_unifrac_betadisp

wu_beta_inter_seed_exp=SILVA_MERDS_rar_WU_dis_betamod_map_seed %>% group_by(betagroup) %>% summarise_at(vars(W_unifrac_betadisp),c(~mean(.),~n()))

ggplot(SILVA_MERDS_rar_WU_dis_betamod_map_seed,aes(x=factor(time, levels = c("Start","End")),y=W_unifrac_betadisp))+
  geom_point(aes(fill=precip, shape=soil_root), size=5)+scale_y_continuous(name = "Betadispersion \nWeighted Unifrac distance")+
  scale_shape_manual(values = c(21,23,24), name=NULL)+
  scale_fill_manual(values = c("blue","red","black"))+
  geom_segment(aes(x=1,xend=2,y=subset(wu_beta_inter_seed_exp,betagroup=="L.B.Start.Start")$mean,
                   yend=subset(wu_beta_inter_seed_exp,betagroup=="L.B.S.A")$mean), linetype="solid", color="blue",size=1.5)+
  geom_segment(aes(x=1,xend=2,y=subset(wu_beta_inter_seed_exp,betagroup=="L.B.Start.Start")$mean,
                   yend=subset(wu_beta_inter_seed_exp,betagroup=="L.B.S.D")$mean), linetype="solid", color="red",size=1.5)+
  geom_segment(aes(x=1,xend=2,y=subset(wu_beta_inter_seed_exp,betagroup=="L.R.Start.Start")$mean,
                   yend=subset(wu_beta_inter_seed_exp,betagroup=="L.R.S.A")$mean), linetype="dashed", color="blue",size=1.5)+
  geom_segment(aes(x=1,xend=2,y=subset(wu_beta_inter_seed_exp,betagroup=="L.R.Start.Start")$mean,
                   yend=subset(wu_beta_inter_seed_exp,betagroup=="L.R.S.D")$mean), linetype="dashed", color="red",size=1.5)+
  geom_segment(aes(x=1,xend=2,y=subset(wu_beta_inter_seed_exp,betagroup=="L.B.Start.Start")$mean,
                   yend=subset(wu_beta_inter_seed_exp,betagroup=="S.B.S.A")$mean), linetype="dotted", color="blue",size=1.5)+
  geom_segment(aes(x=1,xend=2,y=subset(wu_beta_inter_seed_exp,betagroup=="L.R.Start.Start")$mean,
                   yend=subset(wu_beta_inter_seed_exp,betagroup=="S.B.S.D")$mean), linetype="dotted", color="red",size=1.5)+
  theme_bw()+theme(axis.text = element_text(size = 18),axis.title.x = element_blank(), axis.title.y = element_text(size = 22),
                   legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#9.5x7
#beta_disp_line_germination
#####Presence Betadisp seed####
SILVA_MERDS_rar_WU_dis_betamod_map_seed_exp=subset(SILVA_MERDS_rar_WU_dis_betamod_map,life_stage == "S")
SILVA_MERDS_rar_WU_dis_betamod_map_seed_pres=subset(SILVA_MERDS_rar_WU_dis_betamod_map_seed_exp,root_association=="B")


unique(SILVA_MERDS_rar_WU_dis_betamod_map_seed_pres$betagroup)
nrow(SILVA_MERDS_rar_WU_dis_betamod_map_seed_pres)
#12

WU_beta_seed_mod_pres= lm(W_unifrac_betadisp~soil_root*precip, data= SILVA_MERDS_rar_WU_dis_betamod_map_seed_pres)
qqPlot(resid(WU_beta_seed_mod_pres))
hist(resid(WU_beta_seed_mod_pres))
shapiro.test(resid(WU_beta_seed_mod_pres))
#0.9925

Anova(WU_beta_seed_mod_pres, type=3)
#Nada

emmeans(WU_beta_seed_mod_pres, pairwise~s1_precip|s1_soil_root)



#####Origin Betadisp seed####
SILVA_MERDS_rar_WU_dis_betamod_map_seed=subset(SILVA_MERDS_rar_WU_dis_betamod_map,life_stage == "S"|life_stage=="Start")
SILVA_MERDS_rar_WU_dis_betamod_map_seed_L=subset(SILVA_MERDS_rar_WU_dis_betamod_map_seed,soil_status=="L")
unique(SILVA_MERDS_rar_WU_dis_betamod_map_seed_L$betagroup)


wu_beta_seed_inter=SILVA_MERDS_rar_WU_dis_betamod_map_seed_L %>% group_by(betagroup) %>% summarise_at(vars(W_unifrac_betadisp),c(~mean(.),~n()))

ggplot(SILVA_MERDS_rar_WU_dis_betamod_map_seed_L,aes(x=factor(time, levels = c("Start","End")),y=W_unifrac_betadisp))+
  geom_point(aes(fill=precip, shape=soil_root), size=5)+scale_y_continuous(name = "Betadispersion \nWeighted Unifrac distance")+
  scale_shape_manual(name=NULL, values= c(21,24))+scale_fill_manual(values = c("blue","red","black"))+
  geom_segment(aes(x=1,xend=2,y=subset(wu_beta_seed_inter,betagroup=="L.B.Start.Start")$mean,
                   yend=subset(wu_beta_seed_inter,betagroup=="L.B.S.A")$mean), linetype="solid", color="blue",size=1.5)+
  geom_segment(aes(x=1,xend=2,y=subset(wu_beta_seed_inter,betagroup=="L.B.Start.Start")$mean,
                   yend=subset(wu_beta_seed_inter,betagroup=="L.B.S.D")$mean), linetype="solid", color="red",size=1.5)+
  geom_segment(aes(x=1,xend=2,y=subset(wu_beta_seed_inter,betagroup=="L.R.Start.Start")$mean,
                   yend=subset(wu_beta_seed_inter,betagroup=="L.R.S.A")$mean), linetype="dashed", color="blue",size=1.5)+
  geom_segment(aes(x=1,xend=2,y=subset(wu_beta_seed_inter,betagroup=="L.R.Start.Start")$mean,
                   yend=subset(wu_beta_seed_inter,betagroup=="L.R.S.D")$mean), linetype="dashed", color="red",size=1.5)+
  theme_bw()+theme(axis.text = element_text(size = 18),axis.title.x = element_blank(), axis.title.y = element_text(size = 22),
                   legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#+
#  geom_text(aes(x=2,y=W_unifrac_betadisp, label=Row.names))
subset(SILVA_MERDS_rar_WU_dis_betamod_map_seed_L,Row.names=="MERDSSG64")


SILVA_MERDS_rar_WU_dis_betamod_map_seed_ST=subset(SILVA_MERDS_rar_WU_dis_betamod_map_seed_L,time!="Start")
nrow(SILVA_MERDS_rar_WU_dis_betamod_map_seed_ST)
#16

WU_beta_seed_mod= lm(W_unifrac_betadisp~soil_root*precip, data= SILVA_MERDS_rar_WU_dis_betamod_map_seed_ST)
qqPlot(resid(WU_beta_seed_mod))
hist(resid(WU_beta_seed_mod))
shapiro.test(resid(WU_beta_seed_mod))
#0.5671

Anova(WU_beta_seed_mod, type=3)
#Nada

emmeans(WU_beta_seed_mod, pairwise~s1_precip|s1_soil_root)

#Remove outlier

SILVA_MERDS_rar_WU_dis_betamod_map_seed_L_sub=subset(SILVA_MERDS_rar_WU_dis_betamod_map_seed_L,Row.names!="MERDSSG64")
nrow(SILVA_MERDS_rar_WU_dis_betamod_map_seed_L_sub)
#23

wu_beta_seed_S_inter=SILVA_MERDS_rar_WU_dis_betamod_map_seed_L_sub %>% group_by(betagroup) %>% summarise_at(vars(W_unifrac_betadisp),c(~mean(.),~n()))

ggplot(SILVA_MERDS_rar_WU_dis_betamod_map_seed_L_sub,aes(x=factor(time, levels = c("Start","End")),y=W_unifrac_betadisp))+
  geom_point(aes(color=soil_root, shape=precip), size=3)+scale_y_continuous(name = "Betadispersion \nWeighted Unifrac distance")+
  scale_shape_discrete(name=NULL,labels=c("Ambient","Drought","start"))+scale_colour_manual(values = c("grey","black"),labels=c("Bulk","Rhizo"),)+
  geom_segment(aes(x=1,xend=2,y=subset(wu_beta_seed_S_inter,betagroup=="L.B.Start.Start")$mean,
                   yend=subset(wu_beta_seed_S_inter,betagroup=="L.B.S.A")$mean), linetype="solid", color="grey",size=1.5)+
  geom_segment(aes(x=1,xend=2,y=subset(wu_beta_seed_S_inter,betagroup=="L.B.Start.Start")$mean,
                   yend=subset(wu_beta_seed_S_inter,betagroup=="L.B.S.D")$mean), linetype="dashed", color="grey",size=1.5)+
  geom_segment(aes(x=1,xend=2,y=subset(wu_beta_seed_S_inter,betagroup=="L.R.Start.Start")$mean,
                   yend=subset(wu_beta_seed_S_inter,betagroup=="L.R.S.A")$mean), linetype="solid", color="black",size=1.5)+
  geom_segment(aes(x=1,xend=2,y=subset(wu_beta_seed_inter,betagroup=="L.R.Start.Start")$mean,
                   yend=subset(wu_beta_seed_S_inter,betagroup=="L.R.S.D")$mean), linetype="dashed", color="black",size=1.5)+
  theme_bw()+theme(axis.text = element_text(size = 18),axis.title.x = element_blank(), axis.title.y = element_text(size = 22))
#+
#  geom_text(aes(x=2,y=W_unifrac_betadisp, label=Row.names))



SILVA_MERDS_rar_WU_dis_betamod_map_seed_ST_s=subset(SILVA_MERDS_rar_WU_dis_betamod_map_seed_ST,Row.names!="MERDSSG64")
nrow(SILVA_MERDS_rar_WU_dis_betamod_map_seed_ST_s)
#15

WU_beta_seed_s_mod= lm(W_unifrac_betadisp~soil_root*precip, data= SILVA_MERDS_rar_WU_dis_betamod_map_seed_ST_s)
qqPlot(resid(WU_beta_seed_s_mod))
hist(resid(WU_beta_seed_s_mod))
shapiro.test(resid(WU_beta_seed_s_mod))
#0.854

Anova(WU_beta_seed_s_mod, type=3)
#Nada

emmeans(WU_beta_seed_s_mod, pairwise~s1_precip|s1_soil_root)

#####Change from intial####

SILVA_MERDS_rar_map=sample_data(SILVA_MERDS_rar)
head(SILVA_MERDS_rar_map)
SILVA_MERDS_rar_WU_dis_M <- matrixConvert(SILVA_MERDS_rar_WU_dis, 
                                             colname = c("sample1", "sample2", "w_unifrac"))

head(SILVA_MERDS_rar_WU_dis_M)
nrow(SILVA_MERDS_rar_WU_dis_M)
SILVA_MERDS_rar_WU_dis_M$s1_s2=with(SILVA_MERDS_rar_WU_dis_M, interaction(sample1,sample2))

SILVA_MERDS_rar_WU_dis_trt0=merge(SILVA_MERDS_rar_WU_dis_M,SILVA_MERDS_rar_map[,c("soil_status","root_association","block",
                                                                                  "precip","life_stage","soil_root")], 
                                  by.x = "sample1",by.y = "row.names")
nrow(SILVA_MERDS_rar_WU_dis_trt0)

colnames(SILVA_MERDS_rar_WU_dis_trt0)[5:10]=c("s1_soil_status","s1_root_association","s1_block",
                                               "s1_precip","s1_life_stage","s1_soil_root")



SILVA_MERDS_rar_WU_dis_trt=merge(SILVA_MERDS_rar_WU_dis_trt0,SILVA_MERDS_rar_map[,c("soil_status","root_association","block",
                                                                                  "precip","life_stage","soil_root")], 
                                  by.x = "sample2",by.y = "row.names")


nrow(SILVA_MERDS_rar_WU_dis_trt)

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

WU_intercept=SILVA_MERDS_rar_WU_dis_start_wi_bl %>% group_by(time,s1_soil_root,s1_life_stage,s1_precip) %>% summarise_at(vars(w_unifrac),c(~mean(.),~n()))
subset(WU_intercept,s1_soil_root=="L.B"&s1_life_stage=="G"&s1_precip=="A")$mean

ggplot(SILVA_MERDS_rar_WU_dis_start_wi_bl,aes(x=factor(time, levels = c("Start","End")),y=w_unifrac))+
  geom_point(aes(color=s1_soil_root, shape=s1_life_stage), size=2)+
  geom_segment(aes(x=1,xend=2,y=subset(WU_intercept,time=="Start"&s1_soil_root=="L.B")$mean,
                   yend=subset(WU_intercept,s1_soil_root=="L.B"&s1_life_stage=="G"&s1_precip=="A")$mean), linetype="solid", color="red")+
  geom_text(aes(x=2.1, y=subset(WU_intercept,s1_soil_root=="L.B"&s1_life_stage=="G"&s1_precip=="A")$mean, label="Bulk Seedling Ambient"))+
  geom_segment(aes(x=1,xend=2,y=subset(WU_intercept,time=="Start"&s1_soil_root=="L.B")$mean,
                   yend=subset(WU_intercept,s1_soil_root=="L.B"&s1_life_stage=="G"&s1_precip=="D")$mean), linetype="solid", color="red")+
  geom_text(aes(x=2.1, y=subset(WU_intercept,s1_soil_root=="L.B"&s1_life_stage=="G"&s1_precip=="D")$mean, label="Bulk Seedling Drought"))+
  geom_segment(aes(x=1,xend=2,y=subset(WU_intercept,time=="Start"&s1_soil_root=="L.B")$mean,
                   yend=subset(WU_intercept,s1_soil_root=="L.B"&s1_life_stage=="S"&s1_precip=="A")$mean), linetype="longdash", color="red")+
  geom_text(aes(x=2.1, y=subset(WU_intercept,s1_soil_root=="L.B"&s1_life_stage=="S"&s1_precip=="A")$mean, label="Bulk Germination Ambient"))+
  geom_segment(aes(x=1,xend=2,y=subset(WU_intercept,time=="Start"&s1_soil_root=="L.B")$mean,
                   yend=subset(WU_intercept,s1_soil_root=="L.B"&s1_life_stage=="S"&s1_precip=="D")$mean), linetype="longdash", color="red")+
  geom_text(aes(x=2.1, y=subset(WU_intercept,s1_soil_root=="L.B"&s1_life_stage=="S"&s1_precip=="D")$mean+.01, label="Bulk Germination Drought"))+
  geom_segment(aes(x=1,xend=2,y=subset(WU_intercept,time=="Start"&s1_soil_root=="L.B")$mean,
                   yend=subset(WU_intercept,s1_soil_root=="S.B"&s1_life_stage=="G"&s1_precip=="A")$mean), linetype="solid", color="green")+
  geom_segment(aes(x=1,xend=2,y=subset(WU_intercept,time=="Start"&s1_soil_root=="L.B")$mean,
                   yend=subset(WU_intercept,s1_soil_root=="S.B"&s1_life_stage=="G"&s1_precip=="D")$mean), linetype="solid", color="green")+
  geom_segment(aes(x=1,xend=2,y=subset(WU_intercept,time=="Start"&s1_soil_root=="L.B")$mean,
                   yend=subset(WU_intercept,s1_soil_root=="S.B"&s1_life_stage=="S"&s1_precip=="A")$mean), linetype="longdash", color="green")+
  geom_segment(aes(x=1,xend=2,y=subset(WU_intercept,time=="Start"&s1_soil_root=="L.B")$mean,
                   yend=subset(WU_intercept,s1_soil_root=="S.B"&s1_life_stage=="S"&s1_precip=="D")$mean), linetype="longdash", color="green")+
  geom_segment(aes(x=1,xend=2,y=subset(WU_intercept,time=="Start"&s1_soil_root=="L.R")$mean,
                   yend=subset(WU_intercept,s1_soil_root=="L.R"&s1_life_stage=="G"&s1_precip=="A")$mean), linetype="solid", color="blue")+
  geom_text(aes(x=2.1, y=subset(WU_intercept,s1_soil_root=="L.R"&s1_life_stage=="G"&s1_precip=="A")$mean-.01, label="Rhizo Seedling Ambient"))+
  geom_segment(aes(x=1,xend=2,y=subset(WU_intercept,time=="Start"&s1_soil_root=="L.R")$mean,
                   yend=subset(WU_intercept,s1_soil_root=="L.R"&s1_life_stage=="G"&s1_precip=="D")$mean), linetype="solid", color="blue")+
  geom_text(aes(x=2.1, y=subset(WU_intercept,s1_soil_root=="L.R"&s1_life_stage=="G"&s1_precip=="D")$mean-.01, label="Rhizo Seedling Drought"))+
  geom_segment(aes(x=1,xend=2,y=subset(WU_intercept,time=="Start"&s1_soil_root=="L.R")$mean,
                   yend=subset(WU_intercept,s1_soil_root=="L.R"&s1_life_stage=="S"&s1_precip=="A")$mean), linetype="longdash", color="blue")+
  geom_text(aes(x=2.1, y=subset(WU_intercept,s1_soil_root=="L.R"&s1_life_stage=="S"&s1_precip=="A")$mean-.01, label="Rhizo Germination Ambient"))+
  geom_segment(aes(x=1,xend=2,y=subset(WU_intercept,time=="Start"&s1_soil_root=="L.R")$mean,
                   yend=subset(WU_intercept,s1_soil_root=="L.R"&s1_life_stage=="S"&s1_precip=="D")$mean), linetype="longdash", color="blue")+
  geom_text(aes(x=2.1, y=subset(WU_intercept,s1_soil_root=="L.R"&s1_life_stage=="S"&s1_precip=="D")$mean, label="Rhizo Germination Drought"))


#Live

head(SILVA_MERDS_rar_WU_dis_start_wi_bl)
unique(SILVA_MERDS_rar_WU_dis_start_wi_bl$s1_s2_life_stage)

#Transplant


#####Experimental Distance from start community seedling####

SILVA_MERDS_rar_WU_dis_start_wi_bl_trans_exp=subset(SILVA_MERDS_rar_WU_dis_start_wi_bl,s1_s2_life_stage!="S.Start")
nrow(SILVA_MERDS_rar_WU_dis_start_wi_bl_trans_exp)
#36
unique(SILVA_MERDS_rar_WU_dis_start_wi_bl_trans_exp$s1_life_stage)

trans_WU_intercept_exp=SILVA_MERDS_rar_WU_dis_start_wi_bl_trans_exp %>% group_by(time,s1_soil_root,s1_precip,s1_soil_root) %>% summarise_at(vars(w_unifrac),c(~mean(.),~n()))
unique(SILVA_MERDS_rar_WU_dis_start_wi_bl_trans_exp$s1_soil_root)

ggplot(SILVA_MERDS_rar_WU_dis_start_wi_bl_trans_exp,aes(x=factor(time, levels = c("Start","End")),y=w_unifrac))+
  geom_point(aes(fill=s1_precip, shape=s1_soil_root), size=5)+scale_y_continuous(name = "Weighted Unifrac distance")+
  scale_shape_manual(values= c(21,23,24),name=NULL)+scale_fill_manual(values = c("blue","red","black"))+
  geom_segment(aes(x=1,xend=2,y=subset(trans_WU_intercept_exp,time=="Start"&s1_soil_root=="L.B")$mean,
                   yend=subset(trans_WU_intercept_exp,s1_soil_root=="L.B"&s1_precip=="A")$mean), linetype="solid", color="blue",size=1.5)+
  geom_segment(aes(x=1,xend=2,y=subset(trans_WU_intercept_exp,time=="Start"&s1_soil_root=="L.B")$mean,
                   yend=subset(trans_WU_intercept_exp,s1_soil_root=="L.B"&s1_precip=="D")$mean), linetype="solid", color="red",size=1.5)+
  geom_segment(aes(x=1,xend=2,y=subset(trans_WU_intercept_exp,time=="Start"&s1_soil_root=="L.R")$mean,
                   yend=subset(trans_WU_intercept_exp,s1_soil_root=="L.R"&s1_precip=="A")$mean), linetype="dashed", color="blue",size=1.5)+
  geom_segment(aes(x=1,xend=2,y=subset(trans_WU_intercept_exp,time=="Start"&s1_soil_root=="L.R")$mean,
                   yend=subset(trans_WU_intercept_exp,s1_soil_root=="L.R"&s1_precip=="D")$mean), linetype="dashed", color="red",size=1.5)+
  geom_segment(aes(x=1,xend=2,y=subset(trans_WU_intercept_exp,time=="Start"&s1_soil_root=="L.R")$mean,
                   yend=subset(trans_WU_intercept_exp,s1_soil_root=="S.B"&s1_precip=="A")$mean), linetype="dotted", color="blue",size=1.5)+
  geom_segment(aes(x=1,xend=2,y=subset(trans_WU_intercept_exp,time=="Start"&s1_soil_root=="L.R")$mean,
                   yend=subset(trans_WU_intercept_exp,s1_soil_root=="S.B"&s1_precip=="D")$mean), linetype="dotted", color="red",size=1.5)+
  theme_bw()+theme(axis.text = element_text(size = 18),axis.title.x = element_blank(), axis.title.y = element_text(size = 22),
                   legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#comm_turnover_line_seedling

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
trans_WU_intercept=SILVA_MERDS_rar_WU_dis_start_wi_bl_trans %>% group_by(time,s1_soil_root,s1_precip,s1_soil_root) %>% summarise_at(vars(w_unifrac),c(~mean(.),~n()))

ggplot(SILVA_MERDS_rar_WU_dis_start_wi_bl_trans,aes(x=factor(time, levels = c("Start","End")),y=w_unifrac))+
  geom_point(aes(fill=s1_precip, shape=s1_soil_root), size=5)+scale_y_continuous(name = "Weighted Unifrac distance")+
  scale_shape_manual(values= c(21,24),name=NULL)+scale_fill_manual(values = c("blue","red","black"))+
  geom_segment(aes(x=1,xend=2,y=subset(trans_WU_intercept,time=="Start"&s1_soil_root=="L.B")$mean,
                   yend=subset(trans_WU_intercept,s1_soil_root=="L.B"&s1_precip=="A")$mean), linetype="solid", color="blue",size=1.5)+
  geom_segment(aes(x=1,xend=2,y=subset(trans_WU_intercept,time=="Start"&s1_soil_root=="L.B")$mean,
                   yend=subset(trans_WU_intercept,s1_soil_root=="L.B"&s1_precip=="D")$mean), linetype="solid", color="red",size=1.5)+
  geom_segment(aes(x=1,xend=2,y=subset(trans_WU_intercept,time=="Start"&s1_soil_root=="L.R")$mean,
                   yend=subset(trans_WU_intercept,s1_soil_root=="L.R"&s1_precip=="A")$mean), linetype="dashed", color="blue",size=1.5)+
  geom_segment(aes(x=1,xend=2,y=subset(trans_WU_intercept,time=="Start"&s1_soil_root=="L.R")$mean,
                   yend=subset(trans_WU_intercept,s1_soil_root=="L.R"&s1_precip=="D")$mean), linetype="dashed", color="red",size=1.5)+
  theme_bw()+theme(axis.text = element_text(size = 18),axis.title.x = element_blank(), axis.title.y = element_text(size = 22),
                   legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#comm_turnover_line_seedling_soil_origin
SILVA_MERDS_rar_WU_dis_wi_bl_trans=subset(SILVA_MERDS_rar_WU_dis_start_wi_bl_trans,time!="Start")

WU_pair_trans_mod= lm(w_unifrac~s1_soil_root*s1_precip, data= SILVA_MERDS_rar_WU_dis_wi_bl_trans)
qqPlot(resid(WU_pair_trans_mod))
hist(resid(WU_pair_trans_mod))
shapiro.test(resid(WU_pair_trans_mod))
#0.6766

Anova(WU_pair_trans_mod, type=3)
#s1_soil_root:s1_precip 0.01454  1   5.1174   0.04304 * 

emmeans(WU_pair_trans_mod, pairwise~s1_precip*s1_soil_root,adjust="none")

#Remove the betadisp outlier MERDSSG145
head(SILVA_MERDS_rar_WU_dis_start_wi_bl_trans)
SILVA_MERDS_rar_WU_dis_start_wi_bl_trans_sub=subset(SILVA_MERDS_rar_WU_dis_start_wi_bl_trans,sample1!="MERDSSG145")

nrow(SILVA_MERDS_rar_WU_dis_start_wi_bl_trans_sub)
#27
trans_WU_S_intercept=SILVA_MERDS_rar_WU_dis_start_wi_bl_trans_sub %>% group_by(time,s1_soil_root,s1_precip,s1_soil_root) %>% summarise_at(vars(w_unifrac),c(~mean(.),~n()))

ggplot(SILVA_MERDS_rar_WU_dis_start_wi_bl_trans_sub,aes(x=factor(time, levels = c("Start","End")),y=w_unifrac))+
  geom_point(aes(color=s1_soil_root, shape=s1_precip), size=3)+scale_y_continuous(name = "Weighted Unifrac distance")+
  scale_shape_discrete(name=NULL,labels=c("Ambient","Drought","start"))+scale_colour_manual(values = c("grey","black"),labels=c("Bulk","Rhizo"),)+
  geom_segment(aes(x=1,xend=2,y=subset(trans_WU_S_intercept,time=="Start"&s1_soil_root=="L.B")$mean,
                   yend=subset(trans_WU_S_intercept,s1_soil_root=="L.B"&s1_precip=="A")$mean), linetype="solid", color="grey",size=1.5)+
  geom_segment(aes(x=1,xend=2,y=subset(trans_WU_S_intercept,time=="Start"&s1_soil_root=="L.B")$mean,
                   yend=subset(trans_WU_S_intercept,s1_soil_root=="L.B"&s1_precip=="D")$mean), linetype="dashed", color="grey",size=1.5)+
  geom_segment(aes(x=1,xend=2,y=subset(trans_WU_S_intercept,time=="Start"&s1_soil_root=="L.R")$mean,
                   yend=subset(trans_WU_S_intercept,s1_soil_root=="L.R"&s1_precip=="A")$mean), linetype="solid", color="black",size=1.5)+
  geom_segment(aes(x=1,xend=2,y=subset(trans_WU_S_intercept,time=="Start"&s1_soil_root=="L.R")$mean,
                   yend=subset(trans_WU_S_intercept,s1_soil_root=="L.R"&s1_precip=="D")$mean), linetype="dashed", color="black",size=1.5)+
  theme_bw()+theme(axis.text = element_text(size = 18),axis.title.x = element_blank(), axis.title.y = element_text(size = 22))

SILVA_MERDS_rar_WU_dis_wi_bl_trans_sub=subset(SILVA_MERDS_rar_WU_dis_start_wi_bl_trans_sub,time!="Start")
nrow(SILVA_MERDS_rar_WU_dis_wi_bl_trans_sub)
#15

WU_pair_trans_s_mod= lm(w_unifrac~s1_soil_root*s1_precip, data= SILVA_MERDS_rar_WU_dis_wi_bl_trans_sub)
qqPlot(resid(WU_pair_trans_s_mod))
hist(resid(WU_pair_trans_s_mod))
shapiro.test(resid(WU_pair_trans_s_mod))
#0.3983

Anova(WU_pair_trans_s_mod, type=3)
#s1_precip              0.00954  1   5.1011  0.045205 *  
#s1_soil_root:s1_precip 0.02192  1  11.7186  0.005691 ** 

emmeans(WU_pair_trans_s_mod, pairwise~s1_precip*s1_soil_root)



#Seed


#####Experimental Distance from start community germination####
SILVA_MERDS_rar_WU_dis_start_wi_bl_seed_exp=subset(SILVA_MERDS_rar_WU_dis_start_wi_bl,s1_soil_root!="G.B"&s1_s2_life_stage!="G.Start")
nrow(SILVA_MERDS_rar_WU_dis_start_wi_bl_seed_exp)
#35

unique(SILVA_MERDS_rar_WU_dis_start_wi_bl_seed_exp$s1_life_stage)

seed_WU_intercept_exp=SILVA_MERDS_rar_WU_dis_start_wi_bl_seed_exp %>% group_by(time,s1_soil_root,s1_precip,s1_soil_root) %>% summarise_at(vars(w_unifrac),c(~mean(.),~n()))
unique(SILVA_MERDS_rar_WU_dis_start_wi_bl_seed_exp$s1_soil_root)

ggplot(SILVA_MERDS_rar_WU_dis_start_wi_bl_seed_exp,aes(x=factor(time, levels = c("Start","End")),y=w_unifrac))+
  geom_point(aes(fill=s1_precip, shape=s1_soil_root), size=5)+scale_y_continuous(name = "Weighted Unifrac distance")+
  scale_shape_manual(values= c(21,23,24),name=NULL)+scale_fill_manual(values = c("blue","red","black"))+
  geom_segment(aes(x=1,xend=2,y=subset(seed_WU_intercept_exp,time=="Start"&s1_soil_root=="L.B")$mean,
                   yend=subset(seed_WU_intercept_exp,s1_soil_root=="L.B"&s1_precip=="A")$mean), linetype="solid", color="blue",size=1.5)+
  geom_segment(aes(x=1,xend=2,y=subset(seed_WU_intercept_exp,time=="Start"&s1_soil_root=="L.B")$mean,
                   yend=subset(seed_WU_intercept_exp,s1_soil_root=="L.B"&s1_precip=="D")$mean), linetype="solid", color="red",size=1.5)+
  geom_segment(aes(x=1,xend=2,y=subset(seed_WU_intercept_exp,time=="Start"&s1_soil_root=="L.R")$mean,
                   yend=subset(seed_WU_intercept_exp,s1_soil_root=="L.R"&s1_precip=="A")$mean), linetype="dashed", color="blue",size=1.5)+
  geom_segment(aes(x=1,xend=2,y=subset(seed_WU_intercept_exp,time=="Start"&s1_soil_root=="L.R")$mean,
                   yend=subset(seed_WU_intercept_exp,s1_soil_root=="L.R"&s1_precip=="D")$mean), linetype="dashed", color="red",size=1.5)+
  geom_segment(aes(x=1,xend=2,y=subset(seed_WU_intercept_exp,time=="Start"&s1_soil_root=="L.R")$mean,
                   yend=subset(seed_WU_intercept_exp,s1_soil_root=="S.B"&s1_precip=="A")$mean), linetype="dotted", color="blue",size=1.5)+
  geom_segment(aes(x=1,xend=2,y=subset(seed_WU_intercept_exp,time=="Start"&s1_soil_root=="L.R")$mean,
                   yend=subset(seed_WU_intercept_exp,s1_soil_root=="S.B"&s1_precip=="D")$mean), linetype="dotted", color="red",size=1.5)+
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
seed_WU_intercept=SILVA_MERDS_rar_WU_dis_start_wi_bl_seed %>% group_by(time,s1_soil_root,s1_precip,s1_soil_root) %>% summarise_at(vars(w_unifrac),c(~mean(.),~n()))

ggplot(SILVA_MERDS_rar_WU_dis_start_wi_bl_seed,aes(x=factor(time, levels = c("Start","End")),y=w_unifrac))+
  geom_point(aes(fill=s1_precip, shape=s1_soil_root), size=5)+scale_y_continuous(name = "Weighted Unifrac distance")+
  scale_shape_manual(values= c(21,24),name=NULL)+scale_fill_manual(values = c("blue","red", "black"))+
  geom_segment(aes(x=1,xend=2,y=subset(seed_WU_intercept,time=="Start"&s1_soil_root=="L.B")$mean,
                   yend=subset(seed_WU_intercept,s1_soil_root=="L.B"&s1_precip=="A")$mean), linetype="solid", color="blue",size=1.5)+
  geom_segment(aes(x=1,xend=2,y=subset(seed_WU_intercept,time=="Start"&s1_soil_root=="L.B")$mean,
                   yend=subset(seed_WU_intercept,s1_soil_root=="L.B"&s1_precip=="D")$mean), linetype="solid", color="red",size=1.5)+
  geom_segment(aes(x=1,xend=2,y=subset(seed_WU_intercept,time=="Start"&s1_soil_root=="L.R")$mean,
                   yend=subset(seed_WU_intercept,s1_soil_root=="L.R"&s1_precip=="A")$mean), linetype="dashed", color="blue",size=1.5)+
  geom_segment(aes(x=1,xend=2,y=subset(seed_WU_intercept,time=="Start"&s1_soil_root=="L.R")$mean,
                   yend=subset(seed_WU_intercept,s1_soil_root=="L.R"&s1_precip=="D")$mean), linetype="dashed", color="red",size=1.5)+
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

#Remove the outlier MERDSSG64
SILVA_MERDS_rar_WU_dis_start_wi_bl_seed_sub=subset(SILVA_MERDS_rar_WU_dis_start_wi_bl_seed,sample1!="MERDSSG64")

nrow(SILVA_MERDS_rar_WU_dis_start_wi_bl_seed_sub)
#27
seed_WU_s_intercept=SILVA_MERDS_rar_WU_dis_start_wi_bl_seed_sub %>% group_by(time,s1_soil_root,s1_precip,s1_soil_root) %>% summarise_at(vars(w_unifrac),c(~mean(.),~n()))

ggplot(SILVA_MERDS_rar_WU_dis_start_wi_bl_seed_sub,aes(x=factor(time, levels = c("Start","End")),y=w_unifrac))+
  geom_point(aes(color=s1_soil_root, shape=s1_precip), size=3)+scale_y_continuous(name = "Weighted Unifrac distance")+
  scale_shape_discrete(name=NULL,labels=c("Ambient","Drought","start"))+scale_colour_manual(values = c("grey","black"),labels=c("Bulk","Rhizo"),)+
  geom_segment(aes(x=1,xend=2,y=subset(seed_WU_s_intercept,time=="Start"&s1_soil_root=="L.B")$mean,
                   yend=subset(seed_WU_s_intercept,s1_soil_root=="L.B"&s1_precip=="A")$mean), linetype="solid", color="grey",size=1.5)+
  geom_segment(aes(x=1,xend=2,y=subset(seed_WU_s_intercept,time=="Start"&s1_soil_root=="L.B")$mean,
                   yend=subset(seed_WU_s_intercept,s1_soil_root=="L.B"&s1_precip=="D")$mean), linetype="dashed", color="grey",size=1.5)+
  geom_segment(aes(x=1,xend=2,y=subset(seed_WU_s_intercept,time=="Start"&s1_soil_root=="L.R")$mean,
                   yend=subset(seed_WU_s_intercept,s1_soil_root=="L.R"&s1_precip=="A")$mean), linetype="solid", color="black",size=1.5)+
  geom_segment(aes(x=1,xend=2,y=subset(seed_WU_s_intercept,time=="Start"&s1_soil_root=="L.R")$mean,
                   yend=subset(seed_WU_s_intercept,s1_soil_root=="L.R"&s1_precip=="D")$mean), linetype="dashed", color="black",size=1.5)+
  theme_bw()+theme(axis.text = element_text(size = 18),axis.title.x = element_blank(), axis.title.y = element_text(size = 22))

SILVA_MERDS_rar_WU_dis_wi_bl_seed_sub=subset(SILVA_MERDS_rar_WU_dis_start_wi_bl_seed_sub,time!="Start")

WU_pair_seed_s_mod= lm(w_unifrac~s1_soil_root*s1_precip, data= SILVA_MERDS_rar_WU_dis_wi_bl_seed_sub)
qqPlot(resid(WU_pair_seed_s_mod))
hist(resid(WU_pair_seed_s_mod))
shapiro.test(resid(WU_pair_seed_s_mod))
#0.8738

Anova(WU_pair_seed_s_mod, type=3)
#s1_soil_root           0.03987  1  31.1425 0.000165 ***

emmeans(WU_pair_seed_s_mod, pairwise~s1_precip|s1_soil_root)



#####UNWeighted Unifrac####

SILVA_MERDS_rar_unWU_ord=ordinate(SILVA_MERDS_rar, method = "NMDS", distance = "unifrac")
#*** Solution reached
#0.02830039 
plot_ordination(SILVA_MERDS_rar,SILVA_MERDS_rar_unWU_ord, color="root_association",shape="life_stage")+geom_point(size=3)+
  theme_bw()



SILVA_MERDS_rar_map=sample_data(SILVA_MERDS_rar)
SILVA_MERDS_rar_map$soil_root_stage=with(SILVA_MERDS_rar_map, interaction(soil_status,root_association,life_stage))
SILVA_MERDS_rar_unWU_dis=distance(SILVA_MERDS_rar,method = "unifrac")

adonis(SILVA_MERDS_rar_unWU_dis~SILVA_MERDS_rar_map$soil_root_stage+as.factor(SILVA_MERDS_rar_map$block), permutations = 9999)
#                                     Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#SILVA_MERDS_rar_map$soil_root_stage   7    6.2849 0.89784  5.8406 0.46299 0.0001 ***
#as.factor(SILVA_MERDS_rar_map$block)  3    0.5260 0.17533  1.1405 0.03875 0.2469    
#Residuals                            44    6.7638 0.15372         0.49827           
#Total                                54   13.5747                 1.00000       

pairwise.perm.manova(SILVA_MERDS_rar_unWU_dis, SILVA_MERDS_rar_map$soil_root_stage, nperm=2000)

"data:  SILVA_MERDS_rar_unWU_dis by SILVA_MERDS_rar_map$soil_root_stage
2000 permutations 

          L.B.G  S.B.G  L.R.G  L.B.S  S.B.S  L.R.S  L.B.Start
S.B.G     0.0028 -      -      -      -      -      -        
L.R.G     0.6787 0.0028 -      -      -      -      -        
L.B.S     0.0627 0.0028 0.2647 -      -      -      -        
S.B.S     0.0028 0.0038 0.0028 0.0028 -      -      -        
L.R.S     0.1685 0.0028 0.4575 0.0178 0.0028 -      -        
L.B.Start 0.0028 0.0044 0.0044 0.0043 0.0052 0.0028 -        
L.R.Start 0.0049 0.0044 0.0084 0.0043 0.0107 0.0052 0.5872   

P value adjustment method: fdr "



#####Weighted Unifrac####

SILVA_MERDS_rar_WU_ord=ordinate(SILVA_MERDS_rar, method = "NMDS", distance = "wunifrac")
#*** Solution reached
#0.0769077 
plot_ordination(SILVA_MERDS_rar,SILVA_MERDS_rar_WU_ord, color="root_association",shape="life_stage")+geom_point(size=3)+
  theme_bw()



SILVA_MERDS_rar_map=sample_data(SILVA_MERDS_rar)
SILVA_MERDS_rar_map$soil_root_stage=with(SILVA_MERDS_rar_map, interaction(soil_status,root_association,life_stage))
SILVA_MERDS_rar_WU_dis=distance(SILVA_MERDS_rar,method = "wunifrac")

adonis(SILVA_MERDS_rar_WU_dis~SILVA_MERDS_rar_map$soil_root_stage+as.factor(SILVA_MERDS_rar_map$block), permutations = 9999)
#                                     Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)    
#SILVA_MERDS_rar_map$soil_root_stage   7   1.30035 0.185765  13.367 0.66753 0.0001 ***
#as.factor(SILVA_MERDS_rar_map$block)  3   0.03619 0.012062   0.868 0.01858 0.5280    
#Residuals                            44   0.61146 0.013897         0.31389           
#Total                                54   1.94800                  1.00000        

pairwise.perm.manova(SILVA_MERDS_rar_WU_dis, SILVA_MERDS_rar_map$soil_root_stage, nperm=2000)


#####NMDS Compound Graph####
#Let's Graph this

#first what does the combined transplant and seed-start NMDS look like

#live
SILVA_MERDS_rar_live=subset_samples(SILVA_MERDS_rar, soil_status!="S"|life_stage=="Start")
nsamples(SILVA_MERDS_rar_live)
#40
SILVA_MERDS_rar_live_map=sample_data(SILVA_MERDS_rar_live)
SILVA_MERDS_rar_live_ord=ordinate(SILVA_MERDS_rar_live, method = "NMDS",distance = "bray")
#*** Solution reached
#0.1523099 
SILVA_MERDS_rar_live_ord_points=merge(SILVA_MERDS_rar_live_ord$points,SILVA_MERDS_rar_live_map,by="row.names")

(live_trans_ord_p=ggplot(data=SILVA_MERDS_rar_live_ord_points, aes(x=MDS1,y=MDS2))+
    geom_point(size=5, aes(shape=root_association, fill=precip))+stat_ellipse(level=.5, type = "norm",aes(color=interaction(root_association,precip)))+
    theme_bw()+scale_fill_manual(values = c("white","dark grey","black"))+scale_shape_manual(values = c(21,24))+
    ggtitle(label = "Live Soils")+theme(legend.position = "none",plot.title = element_text(hjust = 0.5))+geom_label_repel(size=3,aes(label = life_stage)))

(live_trans_ord_p=ggplot(data=SILVA_MERDS_rar_live_ord_points, aes(x=MDS1,y=MDS2))+
    geom_point(size=5, aes(shape=root_association, fill=precip))+
    theme_bw()+scale_fill_manual(values = c("white","dark grey","black"))+scale_shape_manual(values = c(21,24))+
    ggtitle(label = "Live Soils")+theme(legend.position = "none",plot.title = element_text(hjust = 0.5))+geom_label_repel(size=3,aes(label = life_stage)))

#####No Start####

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


#####UnWeighted Unifrac####
SILVA_MERDS_rar_live_NS=subset_samples(SILVA_MERDS_rar_live, life_stage!="Start")
nsamples(SILVA_MERDS_rar_live_NS)
#32
SILVA_MERDS_rar_live_NS_unWU_ord=ordinate(SILVA_MERDS_rar_live_NS, method = "NMDS", distance = "unifrac")
#*** No convergence -- monoMDS stopping criteria:
#1: no. of iterations >= maxit
#19: stress ratio > sratmax
#0.212476   
plot_ordination(SILVA_MERDS_rar_live_NS,SILVA_MERDS_rar_live_NS_unWU_ord, color="root_association",shape="life_stage")+geom_point(size=3)+
  theme_bw()



SILVA_MERDS_rar_live_NS_map=sample_data(SILVA_MERDS_rar_live_NS)
#SILVA_MERDS_rar_map$soil_root_stage=with(SILVA_MERDS_rar_map, interaction(soil_status,root_association,life_stage))
SILVA_MERDS_rar_live_NS_unWU_dis=distance(SILVA_MERDS_rar_live_NS,method = "unifrac")

(UNwu_full_P_mod=adonis(SILVA_MERDS_rar_live_NS_unWU_dis~SILVA_MERDS_rar_live_NS_map$root_association*
                        SILVA_MERDS_rar_live_NS_map$life_stage*SILVA_MERDS_rar_live_NS_map$precip+
                          as.factor(SILVA_MERDS_rar_live_NS_map$block), permutations = 9999))
#SILVA_MERDS_rar_live_NS_map$root_association                                                                            1    0.1712 0.17117 1.19316 0.03505 0.0754 .  
#SILVA_MERDS_rar_live_NS_map$life_stage                                                                                  1    0.1714 0.17139 1.19469 0.03509 0.0752 .  
#SILVA_MERDS_rar_live_NS_map$precip                                                                                      1    0.3293 0.32931 2.29557 0.06743 0.0001 ***
#as.factor(SILVA_MERDS_rar_live_NS_map$block)                                                                            3    0.5943 0.19811 1.38099 0.12169 0.0002 ***
#SILVA_MERDS_rar_live_NS_map$root_association:SILVA_MERDS_rar_live_NS_map$life_stage                                     1    0.1714 0.17140 1.19480 0.03509 0.0735 .  
AICc.PERMANOVA(UNwu_full_P_mod)
#$AIC
#[1] 57.28953


(UNwu_full_P_mod_no_block=adonis(SILVA_MERDS_rar_live_NS_unWU_dis~SILVA_MERDS_rar_live_NS_map$root_association*
                                 SILVA_MERDS_rar_live_NS_map$life_stage*SILVA_MERDS_rar_live_NS_map$precip, permutations = 9999))

#SILVA_MERDS_rar_live_NS_map$precip                                                                                      1    0.3293 0.32931 2.19547 0.06743 0.0001 ***
#SILVA_MERDS_rar_live_NS_map$root_association:SILVA_MERDS_rar_live_NS_map$life_stage                                     1    0.1753 0.17526 1.16839 0.03588 0.0991 .  

AICc.PERMANOVA(UNwu_full_P_mod_no_block)
#$AIC
#[1] 56.98928


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
    scale_shape_manual(values = c(21,23,24))+scale_fill_manual(values = c("blue","red","black"))+
    theme(axis.title = element_text(size = 20),axis.text = element_text(size = 18),legend.position = "none",
  panel.grid.major = element_blank(), panel.grid.minor = element_blank()))


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
    scale_shape_manual(values = c(21,24))+scale_fill_manual(values = c("blue","red","black"))+
    theme(axis.title = element_text(size = 20),axis.text = element_text(size = 18),legend.position = "none",
          panel.grid.major = element_blank(), panel.grid.minor = element_blank()))


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
    scale_shape_manual(values = c(21,23,24))+scale_fill_manual(values = c("blue","red"))+
    theme(axis.title = element_text(size = 20),axis.text = element_text(size = 18),legend.position = "none",
          panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

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
    scale_shape_manual(values = c(21,24))+scale_fill_manual(values = c("blue","red","black"))+
    theme(axis.title = element_text(size = 20),axis.text = element_text(size = 18),legend.position = "none",
          panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

#Live transplant
SILVA_MERDS_rar_live_trans=subset_samples(SILVA_MERDS_rar, soil_status!="S"&life_stage=="G"|life_stage=="Start")
SILVA_MERDS_rar_live_trans_map=sample_data(SILVA_MERDS_rar_live_trans)
summary(sample_data(SILVA_MERDS_rar_live_trans))
SILVA_MERDS_rar_live_trans_ord=ordinate(SILVA_MERDS_rar_live_trans, method = "NMDS",distance = "bray")
#*** Solution reached
#0.1242186 
SILVA_MERDS_rar_live_trans_ord_points=merge(SILVA_MERDS_rar_live_trans_ord$points,SILVA_MERDS_rar_live_trans_map,by="row.names")
nrow(SILVA_MERDS_rar_live_trans_ord_points)
summary(SILVA_MERDS_rar_live_trans_ord_points)
(live_trans_ord_p=ggplot(data=SILVA_MERDS_rar_live_trans_ord_points, aes(x=MDS1,y=MDS2))+
    geom_point(size=5, aes(shape=root_association, fill=precip))+
    theme_bw()+
    scale_shape_manual(values = c(21,24))+scale_fill_manual(values = c("white","darkgrey","black"))+
    ggtitle(label = "Transplant")+theme(plot.title = element_text(hjust = 0.5)))
(live_trans_ord_p_pot_num=ggplot(data=SILVA_MERDS_rar_live_trans_ord_points, aes(x=MDS1,y=MDS2))+
    geom_point(size=5, aes(shape=root_association, fill=precip))+
  geom_label_repel(size=3,aes(label = Plant_Number))+theme_bw()+scale_fill_manual(values = c("white","dark grey","black"))+scale_shape_manual(values = c(21,24))+
  ggtitle(label = "Microbial Soils")+theme(plot.title = element_text(hjust = 0.5)))

#Let's look at the env fit vectors

(trans_meta_vec <-envfit(SILVA_MERDS_rar_live_trans_ord, data.frame(SILVA_MERDS_rar_live_trans_map[,c("total_biomass","ug_N_NO3_g_dry_soil",
                                                                                           "ug_N_NH4_g_dry_soil_negto0",
                                                                                           "percent_soil_moisture_dry_weight")]), perm=9999, na.rm=TRUE)) #I think i added the environmental data
trans_meta_vec.df <-as.data.frame(trans_meta_vec$vectors$arrows*sqrt(trans_meta_vec$vectors$r))
trans_meta_vec.df$huh <-rownames(trans_meta_vec.df)
scores(trans_meta_vec, "vectors")
summary(trans_meta_vec)

(live_trans_ord_p=ggplot(data=SILVA_MERDS_rar_live_trans_ord_points, aes(x=MDS1,y=MDS2))+
    geom_point(size=5, aes(shape=root_association, fill=precip))+geom_label_repel(size=3,aes(label = ug_N_NO3_g_dry_soil))+
    theme_bw()+geom_segment(data=trans_meta_vec.df, aes(x=0,xend=NMDS1,y=0,yend=NMDS2),
                            arrow = arrow(length=unit(0.5, "cm")), colour="grey")+
    geom_text_repel (data=trans_meta_vec.df, aes(x=NMDS1, y=NMDS2, label=huh), size=5)+
    scale_shape_manual(values = c(21,24))+scale_fill_manual(values = c("white","dark grey","black"))+
    ggtitle(label = "Transplant")+theme(plot.title = element_text(hjust = 0.5)))

#Overlay the total biomass

(live_trans_ord_p=ggplot(data=SILVA_MERDS_rar_live_trans_ord_points, aes(x=MDS1,y=MDS2))+
    geom_point(size=5, aes(shape=root_association, fill=precip))+
    theme_bw()+geom_label_repel(size=3,aes(label = total_biomass))+
    scale_shape_manual(values = c(21,24))+scale_fill_manual(values = c("white","dark grey","black"))+
    ggtitle(label = "Transplant")+theme(plot.title = element_text(hjust = 0.5)))

#Overlay the soil moisture

(live_trans_ord_p=ggplot(data=SILVA_MERDS_rar_live_trans_ord_points, aes(x=MDS1,y=MDS2))+
    geom_point(size=5, aes(shape=root_association, fill=precip))+
    theme_bw()+geom_label_repel(size=3,aes(label = percent_soil_moisture_dry_weight))+
    scale_shape_manual(values = c(21,24))+scale_fill_manual(values = c("white","dark grey","black"))+
    ggtitle(label = "Transplant")+theme(plot.title = element_text(hjust = 0.5)))


#ellipse

(live_trans_ord_p=ggplot(data=SILVA_MERDS_rar_live_trans_ord_points, aes(x=MDS1,y=MDS2))+
    stat_ellipse(level=.5,geom="polygon",alpha=.2,aes(fill=interaction(precip,root_association)))+
    geom_point(size=5, aes(shape=root_association, fill=interaction(precip,root_association)))+
    theme_bw()+scale_fill_brewer(palette = "Paired")+
    scale_shape_manual(values = c(21,24))+
    ggtitle(label = "Transplant")+theme(plot.title = element_text(hjust = 0.5)))

#Hull to group the points

grp.start_bulk <- SILVA_MERDS_rar_live_trans_ord_points[SILVA_MERDS_rar_live_trans_ord_points$precip == "Start"&
                                                           SILVA_MERDS_rar_live_trans_ord_points$root_association == "B", 
                                                         ][chull(SILVA_MERDS_rar_live_trans_ord_points[SILVA_MERDS_rar_live_trans_ord_points$life_stage == 
                                                                                                         "Start"&
                                                                                                         SILVA_MERDS_rar_live_trans_ord_points$root_association == "B", c("MDS1", "MDS2")]), ] 
nrow(grp.start_bulk)
#6
grp.start_rhizo <- SILVA_MERDS_rar_live_trans_ord_points[SILVA_MERDS_rar_live_trans_ord_points$precip == "Start"&
                                                          SILVA_MERDS_rar_live_trans_ord_points$root_association == "R", 
                                                        ][chull(SILVA_MERDS_rar_live_trans_ord_points[SILVA_MERDS_rar_live_trans_ord_points$life_stage == 
                                                                                                        "Start"&
                                                                                                        SILVA_MERDS_rar_live_trans_ord_points$root_association == "R", c("MDS1", "MDS2")]), ]  # hull values 
nrow(grp.start_rhizo)
#3

grp.Amb_Bulk=SILVA_MERDS_rar_live_trans_ord_points[SILVA_MERDS_rar_live_trans_ord_points$precip == "A"&
                                                     SILVA_MERDS_rar_live_trans_ord_points$root_association == "B", 
                                                   ][chull(SILVA_MERDS_rar_live_trans_ord_points[SILVA_MERDS_rar_live_trans_ord_points$precip == "A"&
                                                                                                   SILVA_MERDS_rar_live_trans_ord_points$root_association == "B",c("MDS1", "MDS2")]), ]  # hull values for grp A
nrow(grp.Amb_Bulk)

grp.Amb_Rhizo=SILVA_MERDS_rar_live_trans_ord_points[SILVA_MERDS_rar_live_trans_ord_points$precip == "A"&
                                                     SILVA_MERDS_rar_live_trans_ord_points$root_association == "R", 
                                                   ][chull(SILVA_MERDS_rar_live_trans_ord_points[SILVA_MERDS_rar_live_trans_ord_points$precip == "A"&
                                                                                                   SILVA_MERDS_rar_live_trans_ord_points$root_association == "R",c("MDS1", "MDS2")]), ]  # hull values for grp A
nrow(grp.Amb_Rhizo)


grp.Dro_Bulk=SILVA_MERDS_rar_live_trans_ord_points[SILVA_MERDS_rar_live_trans_ord_points$precip == "D"&
                                                     SILVA_MERDS_rar_live_trans_ord_points$root_association == "B", 
                                                   ][chull(SILVA_MERDS_rar_live_trans_ord_points[SILVA_MERDS_rar_live_trans_ord_points$precip == "D"&
                                                                                                   SILVA_MERDS_rar_live_trans_ord_points$root_association == "B",c("MDS1", "MDS2")]), ]  # hull values for grp A
nrow(grp.Dro_Bulk)

grp.Dro_Rhizo=SILVA_MERDS_rar_live_trans_ord_points[SILVA_MERDS_rar_live_trans_ord_points$precip == "D"&
                                                      SILVA_MERDS_rar_live_trans_ord_points$root_association == "R", 
                                                    ][chull(SILVA_MERDS_rar_live_trans_ord_points[SILVA_MERDS_rar_live_trans_ord_points$precip == "D"&
                                                                                                    SILVA_MERDS_rar_live_trans_ord_points$root_association == "R",c("MDS1", "MDS2")]), ]  # hull values for grp A
nrow(grp.Dro_Rhizo)

hull_trans_live=rbind(grp.start_bulk,grp.start_rhizo,grp.Amb_Bulk,grp.Amb_Rhizo,grp.Dro_Bulk,grp.Dro_Rhizo)
nrow(hull_trans_live)
#23


(live_trans_ord_hull_p=ggplot(data=SILVA_MERDS_rar_live_trans_ord_points, aes(x=MDS1,y=MDS2))+
    geom_polygon(data=hull_trans_live,aes(alpha=precip,x=MDS1,y=MDS2,fill=interaction(precip,root_association),
                                    group=interaction(precip,root_association)))+
    scale_alpha_manual(values = c(0.9,0.5,0.6))+
    geom_point(size=5, aes(shape=root_association, fill=interaction(precip,root_association)))+
    theme_bw()+scale_fill_manual(values = c("gray95","dark grey","black","gray95", "dark grey","black"))+
    scale_shape_manual(values = c(21,24)))

#+
#  ggtitle(label = "Transplant")+theme(plot.title = element_text(hjust = 0.5)))
SILVA_MERDS_rar_live_trans_ord_points$precip_soil_fact=with(SILVA_MERDS_rar_live_trans_ord_points,interaction(precip,root_association))
levels(SILVA_MERDS_rar_live_trans_ord_points$precip_soil_fact) <- c(levels(SILVA_MERDS_rar_live_trans_ord_points$precip_soil_fact), "Start")
SILVA_MERDS_rar_live_trans_ord_points$precip_soil_fact[SILVA_MERDS_rar_live_trans_ord_points$precip_soil_fact=="Start.R"] <- "Start"
SILVA_MERDS_rar_live_trans_ord_points$precip_soil_fact[SILVA_MERDS_rar_live_trans_ord_points$precip_soil_fact=="Start.B"] <- "Start"

grp.start <- SILVA_MERDS_rar_live_trans_ord_points[SILVA_MERDS_rar_live_trans_ord_points$precip_soil_fact == "Start", 
                                                        ][chull(SILVA_MERDS_rar_live_trans_ord_points[SILVA_MERDS_rar_live_trans_ord_points$precip_soil_fact == 
                                                                                                        "Start", c("MDS1", "MDS2")]), ] 
nrow(grp.start)
#8


grp.Amb_Bulk=SILVA_MERDS_rar_live_trans_ord_points[SILVA_MERDS_rar_live_trans_ord_points$precip_soil_fact == "A.B", 
                                                   ][chull(SILVA_MERDS_rar_live_trans_ord_points[SILVA_MERDS_rar_live_trans_ord_points$precip_soil_fact == 
                                                                                                   "A.B", c("MDS1", "MDS2")]), ]
nrow(grp.Amb_Bulk)

grp.Amb_Rhizo=SILVA_MERDS_rar_live_trans_ord_points[SILVA_MERDS_rar_live_trans_ord_points$precip_soil_fact == "A.R", 
                                                    ][chull(SILVA_MERDS_rar_live_trans_ord_points[SILVA_MERDS_rar_live_trans_ord_points$precip_soil_fact == 
                                                                                                    "A.R", c("MDS1", "MDS2")]), ]  # hull values for grp A
nrow(grp.Amb_Rhizo)


grp.Dro_Bulk=SILVA_MERDS_rar_live_trans_ord_points[SILVA_MERDS_rar_live_trans_ord_points$precip_soil_fact == "D.B", 
                                                   ][chull(SILVA_MERDS_rar_live_trans_ord_points[SILVA_MERDS_rar_live_trans_ord_points$precip_soil_fact == 
                                                                                                   "D.B", c("MDS1", "MDS2")]), ]
nrow(grp.Dro_Bulk)

grp.Dro_Rhizo=SILVA_MERDS_rar_live_trans_ord_points[SILVA_MERDS_rar_live_trans_ord_points$precip_soil_fact == "D.R", 
                                                    ][chull(SILVA_MERDS_rar_live_trans_ord_points[SILVA_MERDS_rar_live_trans_ord_points$precip_soil_fact == 
                                                                                                    "D.R", c("MDS1", "MDS2")]), ]
nrow(grp.Dro_Rhizo)

hull_trans_live=rbind(grp.start,grp.Amb_Bulk,grp.Amb_Rhizo,grp.Dro_Bulk,grp.Dro_Rhizo)
nrow(hull_trans_live)
#22

(live_trans_ord_hull2_p=ggplot(data=SILVA_MERDS_rar_live_trans_ord_points, aes(x=MDS1,y=MDS2))+
    geom_polygon(data=hull_trans_live,aes(alpha=precip,x=MDS1,y=MDS2,fill=(precip),
                                          group=(precip_soil_fact)))+
    scale_alpha_manual(values = c(0.9,0.5,0.6))+
    geom_point(size=5, aes(shape=root_association, fill=(precip)))+
    theme_bw()+scale_fill_manual(values = c("gray95","dark grey","black"))+
    scale_shape_manual(values = c(21,24)))


#sterile transplants
SILVA_MERDS_rar_stl_trans=subset_samples(SILVA_MERDS_rar, soil_status=="S"&life_stage=="G")
summary(sample_data(SILVA_MERDS_rar_stl_trans))
SILVA_MERDS_rar_stl_ord=ordinate(SILVA_MERDS_rar_stl_trans, method = "NMDS",distance = "bray")
#*** Solution reached
#Warning message:
#In metaMDS(veganifyOTU(physeq), distance, ...) :
#  stress is (nearly) zero: you may have insufficient data
#0.1223408 
(stl_trans_ord_p=plot_ordination(SILVA_MERDS_rar_stl_trans,SILVA_MERDS_rar_stl_trans_ord,shape="precip")+geom_point(size=5, aes(shape=precip, fill=precip))+
    theme_bw()+scale_fill_manual(values = c("white","dark grey","black"))+scale_shape_manual(values = c(21,24,23)))
(stl_trans_ord_p_pot_num=plot_ordination(SILVA_MERDS_rar_stl_trans,SILVA_MERDS_rar_stl_trans_ord,shape="precip")+geom_point(size=5, aes(shape=precip, fill=precip))+
    geom_label_repel(size=3,aes(label = block))+theme_bw()+scale_fill_manual(values = c("white","dark grey","black"))+scale_shape_manual(values = c(21,24))+
    ggtitle(label = "Microbial Soils")+theme(plot.title = element_text(hjust = 0.5)))




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
    scale_shape_manual(values = c(21,23,24))+scale_fill_manual(values = c("blue","red","black"))+
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
    scale_shape_manual(values = c(21,24))+scale_fill_manual(values = c("blue","red","black"))+
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
    scale_shape_manual(values = c(21,23,24))+scale_fill_manual(values = c("blue","red"))+
    theme(axis.title = element_text(size = 20),axis.text = element_text(size = 18),legend.position = "none",
          panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

nmds_germination

#Live seed
SILVA_MERDS_rar_live_seed=subset_samples(SILVA_MERDS_rar, soil_status!="S"&life_stage=="S"|life_stage=="Start")
SILVA_MERDS_rar_live_seed_map=sample_data(SILVA_MERDS_rar_live_seed)
summary(sample_data(SILVA_MERDS_rar_live_seed))
SILVA_MERDS_rar_live_seed_ord=ordinate(SILVA_MERDS_rar_live_seed, method = "NMDS",distance = "bray")
#*** Solution reached
#0.1153055
SILVA_MERDS_rar_live_seed_ord_points=merge(SILVA_MERDS_rar_live_seed_ord$points,SILVA_MERDS_rar_live_seed_map,by="row.names")
nrow(SILVA_MERDS_rar_live_seed_ord_points)
summary(SILVA_MERDS_rar_live_seed_ord_points)
(live_seed_ord_p=ggplot(data=SILVA_MERDS_rar_live_seed_ord_points, aes(x=MDS1,y=MDS2))+geom_point(size=5, aes(shape=root_association, fill=precip))+
    theme_bw()+scale_fill_manual(values = c("white","dark grey","black"))+scale_shape_manual(values = c(21,24))+
    ggtitle(label = "Seedling")+theme(plot.title = element_text(hjust = 0.5), axis.title.y = element_blank()))
(live_seed_ord_p_pot_num=ggplot(data=SILVA_MERDS_rar_live_seed_ord_points, aes(x=MDS1,y=MDS2))+geom_point(size=5, aes(shape=root_association, fill=precip))+
    geom_label_repel(size=3,aes(label = Plant_Number))+theme_bw()+scale_fill_manual(values = c("white","dark grey","black"))+scale_shape_manual(values = c(21,24))+
    ggtitle(label = "Microbial Soils")+theme(plot.title = element_text(hjust = 0.5)))



#Let's look at the env fit vectors

(seed_meta_vec <-envfit(SILVA_MERDS_rar_live_seed_ord, data.frame(SILVA_MERDS_rar_live_seed_map[,c("total_biomass","ug_N_NO3_g_dry_soil",
                                                                                                      "ug_N_NH4_g_dry_soil_negto0",
                                                                                                      "percent_soil_moisture_dry_weight",
                                                                                                   "days_to_germ")]), perm=9999, na.rm=TRUE)) #I think i added the environmental data
"***VECTORS

NMDS1    NMDS2     r2 Pr(>r)
total_biomass                    -0.21596  0.97640 0.2486 0.4086
ug_N_NO3_g_dry_soil               0.95922 -0.28266 0.3180 0.2607
ug_N_NH4_g_dry_soil_negto0       -0.99795 -0.06395 0.0110 0.9667
percent_soil_moisture_dry_weight -0.97492 -0.22257 0.0398 0.8628
days_to_germ                     -0.68761  0.72608 0.0953 0.6966
Permutation: free
Number of permutations: 9999"

trans_meta_vec.df <-as.data.frame(trans_meta_vec$vectors$arrows*sqrt(trans_meta_vec$vectors$r))
trans_meta_vec.df$huh <-rownames(trans_meta_vec.df)
scores(trans_meta_vec, "vectors")
summary(trans_meta_vec)

#overlay time to germination
(live_seed_ord_p=ggplot(data=SILVA_MERDS_rar_live_seed_ord_points, aes(x=MDS1,y=MDS2))+geom_point(size=5, aes(shape=root_association, fill=precip))+
    theme_bw()+scale_fill_manual(values = c("white","dark grey","black"))+scale_shape_manual(values = c(21,24))+
    geom_label_repel(size=3,aes(label = days_to_germ))+
    ggtitle(label = "Seedling")+theme(plot.title = element_text(hjust = 0.5), axis.title.y = element_blank()))



#Hull to group the points

grp.seed_start <- SILVA_MERDS_rar_live_seed_ord_points[SILVA_MERDS_rar_live_seed_ord_points$precip == "Start", 
                                                   ][chull(SILVA_MERDS_rar_live_seed_ord_points[SILVA_MERDS_rar_live_seed_ord_points$life_stage == 
                                                                                                   "Start", c("MDS1", "MDS2")]), ]  # hull values 
nrow(grp.seed_start)
#5
grp.seed_start_rhizo <- SILVA_MERDS_rar_live_seed_ord_points[SILVA_MERDS_rar_live_seed_ord_points$precip == "Start"&
                                                                SILVA_MERDS_rar_live_seed_ord_points$root_association == "R", 
                                                         ][chull(SILVA_MERDS_rar_live_seed_ord_points[SILVA_MERDS_rar_live_seed_ord_points$life_stage == 
                                                                                                         "Start"&
                                                                                                        SILVA_MERDS_rar_live_seed_ord_points$root_association == "R", c("MDS1", "MDS2")]), ]  # hull values 
nrow(grp.seed_start_rhizo)
#3

grp.seed_Amb_Bulk=SILVA_MERDS_rar_live_seed_ord_points[SILVA_MERDS_rar_live_seed_ord_points$precip == "A"&
                                                         SILVA_MERDS_rar_live_seed_ord_points$root_association == "B", 
                                                   ][chull(SILVA_MERDS_rar_live_seed_ord_points[SILVA_MERDS_rar_live_seed_ord_points$precip == "A"&
                                                                                                  SILVA_MERDS_rar_live_seed_ord_points$root_association == "B",c("MDS1", "MDS2")]), ]  # hull values for grp A
nrow(grp.seed_Amb_Bulk)

grp.seed_Amb_Rhizo=SILVA_MERDS_rar_live_seed_ord_points[SILVA_MERDS_rar_live_seed_ord_points$precip == "A"&
                                                      SILVA_MERDS_rar_live_seed_ord_points$root_association == "R", 
                                                    ][chull(SILVA_MERDS_rar_live_seed_ord_points[SILVA_MERDS_rar_live_seed_ord_points$precip == "A"&
                                                                                                    SILVA_MERDS_rar_live_seed_ord_points$root_association == "R",c("MDS1", "MDS2")]), ]  # hull values for grp A
nrow(grp.seed_Amb_Rhizo)


grp.seed_Dro_Bulk=SILVA_MERDS_rar_live_seed_ord_points[SILVA_MERDS_rar_live_seed_ord_points$precip == "D"&
                                                         SILVA_MERDS_rar_live_seed_ord_points$root_association == "B", 
                                                   ][chull(SILVA_MERDS_rar_live_seed_ord_points[SILVA_MERDS_rar_live_seed_ord_points$precip == "D"&
                                                                                                  SILVA_MERDS_rar_live_seed_ord_points$root_association == "B",c("MDS1", "MDS2")]), ]  # hull values for grp A
nrow(grp.seed_Dro_Bulk)

grp.seed_Dro_Rhizo=SILVA_MERDS_rar_live_seed_ord_points[SILVA_MERDS_rar_live_seed_ord_points$precip == "D"&
                                                          SILVA_MERDS_rar_live_seed_ord_points$root_association == "R", 
                                                    ][chull(SILVA_MERDS_rar_live_seed_ord_points[SILVA_MERDS_rar_live_seed_ord_points$precip == "D"&
                                                                                                   SILVA_MERDS_rar_live_seed_ord_points$root_association == "R",c("MDS1", "MDS2")]), ]  # hull values for grp A
nrow(grp.seed_Dro_Rhizo)

hull_seed_live=rbind(grp.seed_start,grp.seed_Amb_Bulk,grp.seed_Amb_Rhizo,grp.seed_Dro_Bulk,grp.seed_Dro_Rhizo)
nrow(hull_seed_live)
#20


(live_seed_ord_hull_p=ggplot(data=SILVA_MERDS_rar_live_seed_ord_points, aes(x=MDS1,y=MDS2))+
    geom_polygon(data=hull_seed_live,aes(alpha=precip,x=MDS1,y=MDS2,fill=interaction(precip,root_association),
                                          group=interaction(precip,root_association)))+
    scale_alpha_manual(values = c(0.9,0.5,0.5))+
    geom_point(size=5, aes(shape=root_association, fill=interaction(precip,root_association)))+
    theme_bw()+scale_fill_manual(values = c("gray95","dark grey","black","gray95", "dark grey","black"))+
    scale_shape_manual(values = c(21,24))+
    ggtitle(label = "Seed-Start")+theme(plot.title = element_text(hjust = 0.5)))

#new version of Hull
SILVA_MERDS_rar_live_seed_ord_points$precip_soil_fact=with(SILVA_MERDS_rar_live_seed_ord_points,interaction(precip,root_association))
levels(SILVA_MERDS_rar_live_seed_ord_points$precip_soil_fact) <- c(levels(SILVA_MERDS_rar_live_seed_ord_points$precip_soil_fact), "Start")
SILVA_MERDS_rar_live_seed_ord_points$precip_soil_fact[SILVA_MERDS_rar_live_seed_ord_points$precip_soil_fact=="Start.R"] <- "Start"
SILVA_MERDS_rar_live_seed_ord_points$precip_soil_fact[SILVA_MERDS_rar_live_seed_ord_points$precip_soil_fact=="Start.B"] <- "Start"

grp.seed_start <- SILVA_MERDS_rar_live_seed_ord_points[SILVA_MERDS_rar_live_seed_ord_points$precip_soil_fact == "Start", 
                                                   ][chull(SILVA_MERDS_rar_live_seed_ord_points[SILVA_MERDS_rar_live_seed_ord_points$precip_soil_fact == 
                                                                                                   "Start", c("MDS1", "MDS2")]), ] 
nrow(grp.seed_start)
#5


grp.seed_Amb_Bulk=SILVA_MERDS_rar_live_seed_ord_points[SILVA_MERDS_rar_live_seed_ord_points$precip_soil_fact == "A.B", 
                                                   ][chull(SILVA_MERDS_rar_live_seed_ord_points[SILVA_MERDS_rar_live_seed_ord_points$precip_soil_fact == 
                                                                                                   "A.B", c("MDS1", "MDS2")]), ]
nrow(grp.seed_Amb_Bulk)
#4
grp.seed_Amb_Rhizo=SILVA_MERDS_rar_live_seed_ord_points[SILVA_MERDS_rar_live_seed_ord_points$precip_soil_fact == "A.R", 
                                                    ][chull(SILVA_MERDS_rar_live_seed_ord_points[SILVA_MERDS_rar_live_seed_ord_points$precip_soil_fact == 
                                                                                                    "A.R", c("MDS1", "MDS2")]), ]  # hull values for grp A
nrow(grp.seed_Amb_Rhizo)
#4

grp.seed_Dro_Bulk=SILVA_MERDS_rar_live_seed_ord_points[SILVA_MERDS_rar_live_seed_ord_points$precip_soil_fact == "D.B", 
                                                   ][chull(SILVA_MERDS_rar_live_seed_ord_points[SILVA_MERDS_rar_live_seed_ord_points$precip_soil_fact == 
                                                                                                   "D.B", c("MDS1", "MDS2")]), ]
nrow(grp.seed_Dro_Bulk)

grp.seed_Dro_Rhizo=SILVA_MERDS_rar_live_seed_ord_points[SILVA_MERDS_rar_live_seed_ord_points$precip_soil_fact == "D.R", 
                                                    ][chull(SILVA_MERDS_rar_live_seed_ord_points[SILVA_MERDS_rar_live_seed_ord_points$precip_soil_fact == 
                                                                                                    "D.R", c("MDS1", "MDS2")]), ]
nrow(grp.seed_Dro_Rhizo)
#4
hull_seed_live=rbind(grp.seed_start,grp.seed_Amb_Bulk,grp.seed_Amb_Rhizo,grp.seed_Dro_Bulk,grp.seed_Dro_Rhizo)
nrow(hull_seed_live)
#20

(live_seed_ord_hull2_p=ggplot(data=SILVA_MERDS_rar_live_seed_ord_points, aes(x=MDS1,y=MDS2))+
    geom_polygon(data=hull_seed_live,aes(alpha=precip,x=MDS1,y=MDS2,fill=(precip),
                                          group=(precip_soil_fact)))+
    scale_alpha_manual(values = c(0.9,0.5,0.6))+
    geom_point(size=5, aes(shape=root_association, fill=(precip)))+
    theme_bw()+scale_fill_manual(values = c("gray95","dark grey","black"))+
    scale_shape_manual(values = c(21,24)))



#sterile seed
SILVA_MERDS_rar_stl_seed=subset_samples(SILVA_MERDS_rar, soil_status=="S"&life_stage=="S")
summary(sample_data(SILVA_MERDS_rar_stl_seed))
SILVA_MERDS_rar_stl_seed_ord=ordinate(SILVA_MERDS_rar_stl_seed, method = "NMDS",distance = "bray")
#*** Solution reached
#Warning message:
#In metaMDS(veganifyOTU(physeq), distance, ...) :
#  stress is (nearly) zero: you may have insufficient data
#8.3787e-05 
(stl_seed_ord_p=plot_ordination(SILVA_MERDS_rar_stl_seed,SILVA_MERDS_rar_stl_seed_ord,shape="precip")+geom_point(size=5, aes(shape=precip, fill=precip))+
    theme_bw()+scale_fill_manual(values = c("white","dark grey","black"))+scale_shape_manual(values = c(21,24,23)))
(stl_seed_ord_p_pot_num=plot_ordination(SILVA_MERDS_rar_stl_seed,SILVA_MERDS_rar_stl_seed_ord,shape="precip")+geom_point(size=5, aes(shape=precip, fill=precip))+
    geom_label_repel(size=3,aes(label = block))+theme_bw()+scale_fill_manual(values = c("white","dark grey","black"))+scale_shape_manual(values = c(21,24))+
    ggtitle(label = "Sterile")+theme(plot.title = element_text(hjust = 0.5)))

ggarrange(live_trans_ord_p,live_seed_ord_p,ncol = 2,  legend = "none")
#10x5.5

ggarrange(live_trans_ord_hull2_p,live_seed_ord_hull2_p,ncol = 2,  legend = "none")
#10x5.5
#Remove the sterile since it seems to be driving a lot of this

SILVA_MERDS_rar_live=subset_samples(SILVA_MERDS_rar, soil_status!="S")

SILVA_MERDS_rar_live_ord=ordinate(SILVA_MERDS_rar_live, method = "NMDS",distance = "bray")
#*** No convergence -- monoMDS stopping criteria:
#20: stress ratio > sratmax
#0.1523099
plot_ordination(SILVA_MERDS_rar_live,SILVA_MERDS_rar_live_ord, color="root_association",shape="precip")+geom_point(size=3)+
  geom_label_repel(size=3,aes(label = Plant_Number))+theme_bw()
plot_ordination(SILVA_MERDS_rar_live,SILVA_MERDS_rar_live_ord, color="root_association",shape="precip")+geom_point(size=3)+
  geom_label_repel(size=3,aes(label = block))+theme_bw()

SILVA_MERDS_rar_live_map=sample_data(SILVA_MERDS_rar_live)
SILVA_MERDS_rar_live_map$soil_root_stage_precip=with(SILVA_MERDS_rar_live_map, interaction(soil_status,root_association,life_stage,precip))
SILVA_MERDS_rar_live_map$soil_root=with(SILVA_MERDS_rar_live_map, interaction(soil_status,root_association))
summary(SILVA_MERDS_rar_live_map)
SILVA_MERDS_rar_live_dis=distance(SILVA_MERDS_rar_live,method = "bray")

adonis(SILVA_MERDS_rar_live_dis~SILVA_MERDS_rar_live_map$root_association+SILVA_MERDS_rar_live_map$precip+SILVA_MERDS_rar_live_map$life_stage
       +as.factor(SILVA_MERDS_rar_live_map$block), permutations = 9999)
#SILVA_MERDS_rar_live_map$root_association  1    0.2394 0.23940  2.0762 0.03660 0.0245 *  
#SILVA_MERDS_rar_live_map$precip            2    1.7520 0.87599  7.5969 0.26785 0.0001 ***
#as.factor(SILVA_MERDS_rar_live_map$block)  3    0.7091 0.23636  2.0498 0.10841 0.0020 ** 

adonis(SILVA_MERDS_rar_live_dis~SILVA_MERDS_rar_live_map$soil_root_stage_precip
       +as.factor(SILVA_MERDS_rar_live_map$block), permutations = 9999)
#SILVA_MERDS_rar_live_map$soil_root_stage_precip  9    2.9071 0.32301  2.9639 0.44444 0.0001 ***
#as.factor(SILVA_MERDS_rar_live_map$block)        3    0.6914 0.23045  2.1146 0.10570 0.0014 ** 
pairwise.perm.manova(SILVA_MERDS_rar_live_dis, SILVA_MERDS_rar_live_map$soil_root_stage_precip, nperm=2000)


"> pairwise.perm.manova(SILVA_MERDS_rar_live_dis, SILVA_MERDS_rar_live_map$soil_root_stage_precip, nperm=2000)

	Pairwise comparisons using permutation MANOVAs on a distance matrix 

data:  SILVA_MERDS_rar_live_dis by SILVA_MERDS_rar_live_map$soil_root_stage_precip
2000 permutations 

                L.B.G.A L.R.G.A L.B.S.A L.R.S.A L.B.G.D L.R.G.D L.B.S.D L.R.S.D L.B.Start.Start
L.R.G.A         0.443   -       -       -       -       -       -       -       -              
L.B.S.A         0.148   0.124   -       -       -       -       -       -       -              
L.R.S.A         0.100   0.866   0.100   -       -       -       -       -       -              
L.B.G.D         0.066   0.520   0.066   0.174   -       -       -       -       -              
L.R.G.D         0.066   0.174   0.066   0.115   0.296   -       -       -       -              
L.B.S.D         0.066   0.100   0.066   0.066   0.142   0.628   -       -       -              
L.R.S.D         0.066   0.148   0.066   0.117   0.402   0.194   0.100   -       -              
L.B.Start.Start 0.066   0.066   0.066   0.066   0.066   0.066   0.066   0.066   -              
L.R.Start.Start 0.066   0.066   0.067   0.066   0.066   0.066   0.066   0.066   0.936          
"


#Only the starting community 
#Are they significantly different?
SILVA_MERDS_rar_start=subset_samples(SILVA_MERDS_rar, life_stage=="Start")

SILVA_MERDS_rar_start_ord=ordinate(SILVA_MERDS_rar_start, method = "NMDS",distance = "bray")
#*** Solution reached
#0.03181609
plot_ordination(SILVA_MERDS_rar_start,SILVA_MERDS_rar_start_ord, color="root_association",shape="as.factor(block)")+geom_point(size=3)+
  geom_label_repel(size=3,aes(label = block))+theme_bw()

SILVA_MERDS_rar_start_map=sample_data(SILVA_MERDS_rar_start)
SILVA_MERDS_rar_start_dis=distance(SILVA_MERDS_rar_start,method = "bray")
adonis(SILVA_MERDS_rar_start_dis~SILVA_MERDS_rar_start_map$root_association
       +as.factor(SILVA_MERDS_rar_start_map$block), permutations = 9999)

#SILVA_MERDS_rar_start_map$root_association  1   0.07476 0.074756 0.76287 0.10453 0.7764
#as.factor(SILVA_MERDS_rar_start_map$block)  3   0.34643 0.115477 1.17844 0.48441 0.2405




#Only treatments

SILVA_MERDS_rar_trt=subset_samples(SILVA_MERDS_rar, life_stage!="Start")




SILVA_MERDS_rar_trt_ord=ordinate(SILVA_MERDS_rar_trt, method = "NMDS",distance = "bray")
#*** Solution reached
#Warning message:
#  In metaMDS(veganifyOTU(physeq), distance, ...) :
#  stress is (nearly) zero: you may have insufficient data
#8.121468e-05
plot_ordination(SILVA_MERDS_rar_trt,SILVA_MERDS_rar_trt_ord, color="precip",shape="life_stage", label = "block")

SILVA_MERDS_rar_trt_map=sample_data(SILVA_MERDS_rar_trt)
SILVA_MERDS_rar_trt_map$soil_root=with(SILVA_MERDS_rar_trt_map, interaction(soil_status,root_association))
SILVA_MERDS_rar_trt_dis=distance(SILVA_MERDS_rar_trt,method = "bray")

adonis(SILVA_MERDS_rar_trt_dis~SILVA_MERDS_rar_trt_map$soil_root*SILVA_MERDS_rar_trt_map$precip*SILVA_MERDS_rar_trt_map$life_stage
       +as.factor(SILVA_MERDS_rar_trt_map$block), permutations = 9999)
"Terms added sequentially (first to last)

                                                                                                    Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
SILVA_MERDS_rar_trt_map$soil_root                                                                    2    4.3301 2.16507 17.4498 0.37696 0.0001 ***
SILVA_MERDS_rar_trt_map$precip                                                                       1    0.5459 0.54588  4.3997 0.04752 0.0017 ** 
SILVA_MERDS_rar_trt_map$life_stage                                                                   1    0.3409 0.34091  2.7477 0.02968 0.0160 *  
as.factor(SILVA_MERDS_rar_trt_map$block)                                                             3    0.5862 0.19539  1.5748 0.05103 0.0603 .  
SILVA_MERDS_rar_trt_map$soil_root:SILVA_MERDS_rar_trt_map$precip                                     2    0.4514 0.22568  1.8189 0.03929 0.0451 *  
SILVA_MERDS_rar_trt_map$soil_root:SILVA_MERDS_rar_trt_map$life_stage                                 2    0.7489 0.37447  3.0181 0.06520 0.0025 ** 
SILVA_MERDS_rar_trt_map$precip:SILVA_MERDS_rar_trt_map$life_stage                                    1    0.1749 0.17487  1.4094 0.01522 0.1648    
SILVA_MERDS_rar_trt_map$soil_root:SILVA_MERDS_rar_trt_map$precip:SILVA_MERDS_rar_trt_map$life_stage  2    0.3382 0.16910  1.3629 0.02944 0.1523    
Residuals                                                                                           32    3.9704 0.12407         0.34564           
Total                                                                                               46   11.4869                 1.00000   "


#Remove the sterile since it seems to be driving a lot of this

SILVA_MERDS_rar_trt_live=subset_samples(SILVA_MERDS_rar_trt, soil_status!="S")

SILVA_MERDS_rar_trt_live_ord=ordinate(SILVA_MERDS_rar_trt_live, method = "NMDS",distance = "bray")
#*** Solution reached
#0.1815625
plot_ordination(SILVA_MERDS_rar_trt_live,SILVA_MERDS_rar_trt_live_ord, color="precip",shape="life_stage", label = "block")


SILVA_MERDS_rar_trt_live_map=sample_data(SILVA_MERDS_rar_trt_live)
SILVA_MERDS_rar_trt_live_map$soil_root=with(SILVA_MERDS_rar_trt_live_map, interaction(soil_status,root_association))
SILVA_MERDS_rar_trt_live_dis=distance(SILVA_MERDS_rar_trt_live,method = "bray")

adonis(SILVA_MERDS_rar_trt_live_dis~SILVA_MERDS_rar_trt_live_map$soil_root*SILVA_MERDS_rar_trt_live_map$precip*SILVA_MERDS_rar_trt_live_map$life_stage
       +as.factor(SILVA_MERDS_rar_trt_live_map$block), permutations = 9999)


"Terms added sequentially (first to last)

                                                                                                                    Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
SILVA_MERDS_rar_trt_live_map$soil_root                                                                              1    0.2775 0.27753  2.4486 0.06084 0.0024 ** 
SILVA_MERDS_rar_trt_live_map$precip                                                                                 1    0.4882 0.48819  4.3072 0.10701 0.0001 ***
SILVA_MERDS_rar_trt_live_map$life_stage                                                                             1    0.1505 0.15055  1.3282 0.03300 0.1419    
as.factor(SILVA_MERDS_rar_trt_live_map$block)                                                                       3    0.6329 0.21098  1.8614 0.13874 0.0016 ** 
SILVA_MERDS_rar_trt_live_map$soil_root:SILVA_MERDS_rar_trt_live_map$precip                                          1    0.1334 0.13340  1.1770 0.02924 0.2427    
SILVA_MERDS_rar_trt_live_map$soil_root:SILVA_MERDS_rar_trt_live_map$life_stage                                      1    0.2150 0.21495  1.8965 0.04712 0.0187 *  
SILVA_MERDS_rar_trt_live_map$precip:SILVA_MERDS_rar_trt_live_map$life_stage                                         1    0.1049 0.10486  0.9251 0.02299 0.5141    
SILVA_MERDS_rar_trt_live_map$soil_root:SILVA_MERDS_rar_trt_live_map$precip:SILVA_MERDS_rar_trt_live_map$life_stage  1    0.1793 0.17929  1.5818 0.03930 0.0605 .  
Residuals                                                                                                          21    2.3802 0.11334         0.52176           
Total                                                                                                              31    4.5619                 1.00000           
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
> "

#there is a marginally significant three way interaction with life stage 
#since our other analyses are split by life stage I am going to split them here 

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

#beta dispersion
trans_perm.mod=betadisper(SILVA_MERDS_rar_trt_live_transe_dis, with(SILVA_MERDS_rar_trt_live_trans_map, interaction(soil_root,precip)))
plot(trans_perm.mod)
boxplot(trans_perm.mod)
anova(trans_perm.mod)
(mod.HSDtrans<- TukeyHSD(trans_perm.mod))
plot(mod.HSDtrans)
trans_perm.mod_perm <- permutest(trans_perm.mod, permutations = 9999, pairwise = TRUE)
trans_perm.mod_stat <- permustats(trans_perm.mod_perm)
summary(trans_perm.mod_stat)


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

#beta dispersion
seed_perm.mod=betadisper(SILVA_MERDS_rar_trt_live_seed_dis, with(SILVA_MERDS_rar_trt_live_seed_map, interaction(soil_root,precip)))
plot(seed_perm.mod)
boxplot(seed_perm.mod)
anova(seed_perm.mod)
(mod.HSDseed <- TukeyHSD(seed_perm.mod))
plot(mod.HSD015)
seed_perm.mod_perm <- permutest(seed_perm.mod, permutations = 9999, pairwise = TRUE)
seed_perm.mod_stat <- permustats(seed_perm.mod_perm)
summary(seed_perm.mod_stat)



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

(sterile_ord_p=plot_ordination(SILVA_MERDS_rar_sterile,SILVA_MERDS_rar_sterile_ord,shape="life_stage")+geom_point(size=5, aes(shape=life_stage, fill=precip))+
    theme_bw()+scale_fill_manual(values = c("white","dark grey","black"))+scale_shape_manual(values = c(22,23))+
    ggtitle(label = "Sterile")+theme(plot.title = element_text(hjust = 0.5)))

(sterile_ord_p=plot_ordination(SILVA_MERDS_rar_sterile,SILVA_MERDS_rar_sterile_ord,shape="life_stage")+geom_point(size=5, aes(shape=life_stage, fill=precip))+
    theme_bw()+scale_fill_manual(values = c("white","dark grey","black"))+scale_shape_manual(values = c(22,23))+
    ggtitle(label = "Sterile")+theme(plot.title = element_text(hjust = 0.5)))

#Let's look at the env fit vectors

(sterile_meta_vec <-envfit(SILVA_MERDS_rar_sterile_ord, data.frame(SILVA_MERDS_rar_sterile_map[,c("total_biomass","ug_N_NO3_g_dry_soil",
                                                                                                      "ug_N_NH4_g_dry_soil_negto0",
                                                                                                      "percent_soil_moisture_dry_weight")]), perm=9999, na.rm=TRUE)) #I think i added the environmental data
#percent_soil_moisture_dry_weight -0.60170  0.79872 0.5765 0.0231 *
sterile_meta_vec.df <-as.data.frame(sterile_meta_vec$vectors$arrows*sqrt(sterile_meta_vec$vectors$r))
sterile_meta_vec.df$huh <-rownames(sterile_meta_vec.df)
scores(sterile_meta_vec, "vectors")
summary(sterile_meta_vec.df)

(sterile_ord_p=plot_ordination(SILVA_MERDS_rar_sterile,SILVA_MERDS_rar_sterile_ord)+geom_point(size=5, aes(shape=life_stage, fill=precip))+
    theme_bw()+geom_text_repel(size=5, aes(label=percent_soil_moisture_dry_weight))+
    scale_fill_manual(values = c("white","dark grey","black"))+scale_shape_manual(values = c(22,23))+
    ggtitle(label = "Sterile")+geom_segment(data=sterile_meta_vec.df, aes(x=0,xend=NMDS1,y=0,yend=NMDS2),
                                            arrow = arrow(length=unit(0.5, "cm")), colour="grey")+
    geom_text_repel (data=sterile_meta_vec.df, aes(x=NMDS1, y=NMDS2, label=huh), size=5)+theme(plot.title = element_text(hjust = 0.5)))





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

#SIMPER Analyses

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



#Let look at transplants and see if there is a correlation between biomass and OTUs

SILVA_MERDS_rar_sterile_trans=subset_samples(SILVA_MERDS_rar_sterile, life_stage=="G")
sample_data(SILVA_MERDS_rar_sterile_trans)
source(system.file("extdata/lm_phyloseq.R", package = "microbiome"))
biomas_mod_otus=lm_phyloseq(SILVA_MERDS_rar_sterile_trans, "total_biomass")
#Warning messages:
#1: In transform(x, transformation) :
#  OTU table contains zeroes. Using log10(1 + x) transform.
#2: In fitFDist(var, df1 = df, covariate = covariate) :
#  More than half of residual variances are exactly zero: eBayes unreliable
"logFC    AveExpr        t      P.Value   adj.P.Val        B
OTU259   1.4906662 0.27848584 9.925403 2.800264e-05 0.003587517 2.036929
OTU2653  0.9784821 0.18279975 9.925295 2.800460e-05 0.003587517 2.036849
OTU1215  0.5654503 0.10563726 9.924914 2.801145e-05 0.003587517 2.036572
OTU13103 0.5654503 0.10563726 9.924914 2.801145e-05 0.003587517 2.036572
OTU1866  0.5654503 0.10563726 9.924914 2.801145e-05 0.003587517 2.036572
OTU28945 0.5654503 0.10563726 9.924914 2.801145e-05 0.003587517 2.036572
OTU18008 0.5654503 0.10563726 9.924914 2.801145e-05 0.003587517 2.036572
OTU10541 0.5206565 0.09726891 9.924812 2.801330e-05 0.003587517 2.036497
OTU8891  0.5206565 0.09726891 9.924812 2.801330e-05 0.003587517 2.036497
OTU625   0.4676768 0.08737125 9.924651 2.801620e-05 0.003587517 2.036380"


get_sample(SILVA_MERDS_rar_sterile_trans,"OTU259")
get_sample(SILVA_MERDS_rar_sterile_trans,"OTU2653")
get_sample(SILVA_MERDS_rar_sterile_trans,"OTU1215")
get_sample(SILVA_MERDS_rar_sterile_trans,"OTU13103")
get_sample(SILVA_MERDS_rar_sterile_trans,"OTU1866")
get_sample(SILVA_MERDS_rar_sterile_trans,"OTU28945")
get_sample(SILVA_MERDS_rar_sterile_trans,"OTU18008")
get_sample(SILVA_MERDS_rar_sterile_trans,"OTU10541")
get_sample(SILVA_MERDS_rar_sterile_trans,"OTU8891")
get_sample(SILVA_MERDS_rar_sterile_trans,"OTU625")


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

#richness
bact_obs_rich_full_model= lm(log(Observed)~soil_root_stage+as.factor(block), data= SILVA_MERDS_rar.divfil)
qqPlot(resid(bact_obs_rich_full_model))
hist(resid(bact_obs_rich_full_model))

Anova(bact_obs_rich_full_model, type=3)
#soil_root_stage    43.36  7   174.6772 <2e-16 ***

emmeans(bact_obs_rich_full_model, pairwise~soil_root_stage, adjust="fdr")
"> emmeans(bact_obs_rich_full_model, pairwise~soil_root_stage)
$emmeans
soil_root_stage emmean     SE df lower.CL upper.CL
L.B.G             7.18 0.0666 44     7.04     7.31
S.B.G             5.37 0.0666 44     5.24     5.50
L.R.G             7.14 0.0666 44     7.00     7.27
L.B.S             6.92 0.0666 44     6.78     7.05
S.B.S             5.13 0.0722 44     4.99     5.28
L.R.S             7.21 0.0672 44     7.08     7.35
L.B.Start         7.63 0.0942 44     7.44     7.82
L.R.Start         7.64 0.0942 44     7.45     7.83

Results are averaged over the levels of: block 
Results are given on the log (not the response) scale. 
Confidence level used: 0.95 

$contrasts
contrast              estimate     SE df t.ratio p.value
L.B.G - S.B.G           1.8081 0.0942 44  19.203 <.0001 
L.B.G - L.R.G           0.0405 0.0942 44   0.430 0.7205 
L.B.G - L.B.S           0.2618 0.0942 44   2.780 0.0101 
L.B.G - S.B.S           2.0439 0.0982 44  20.807 <.0001 
L.B.G - L.R.S          -0.0343 0.0946 44  -0.362 0.7456 
L.B.G - L.B.Start      -0.4541 0.1153 44  -3.938 0.0005 
L.B.G - L.R.Start      -0.4584 0.1153 44  -3.975 0.0004 
S.B.G - L.R.G          -1.7675 0.0942 44 -18.772 <.0001 
S.B.G - L.B.S          -1.5463 0.0942 44 -16.422 <.0001 
S.B.G - S.B.S           0.2358 0.0982 44   2.401 0.0252 
S.B.G - L.R.S          -1.8423 0.0946 44 -19.476 <.0001 
S.B.G - L.B.Start      -2.2622 0.1153 44 -19.617 <.0001 
S.B.G - L.R.Start      -2.2665 0.1153 44 -19.654 <.0001 
L.R.G - L.B.S           0.2212 0.0942 44   2.350 0.0272 
L.R.G - S.B.S           2.0034 0.0982 44  20.394 <.0001 
L.R.G - L.R.S          -0.0748 0.0946 44  -0.791 0.4855 
L.R.G - L.B.Start      -0.4946 0.1153 44  -4.289 0.0002 
L.R.G - L.R.Start      -0.4989 0.1153 44  -4.327 0.0002 
L.B.S - S.B.S           1.7821 0.0982 44  18.142 <.0001 
L.B.S - L.R.S          -0.2960 0.0946 44  -3.129 0.0041 
L.B.S - L.B.Start      -0.7159 0.1153 44  -6.208 <.0001 
L.B.S - L.R.Start      -0.7202 0.1153 44  -6.245 <.0001 
S.B.S - L.R.S          -2.0781 0.0996 44 -20.863 <.0001 
S.B.S - L.B.Start      -2.4980 0.1187 44 -21.050 <.0001 
S.B.S - L.R.Start      -2.5023 0.1187 44 -21.086 <.0001 
L.R.S - L.B.Start      -0.4199 0.1157 44  -3.630 0.0010 
L.R.S - L.R.Start      -0.4242 0.1157 44  -3.667 0.0010 
L.B.Start - L.R.Start  -0.0043 0.1332 44  -0.032 0.9744 

Results are averaged over the levels of: block 
Results are given on the log (not the response) scale. 
P value adjustment: fdr method for 28 tests 
"


ggplot(SILVA_MERDS_rar.divfil, aes(x=soil_status, y=InvSimpson))+geom_boxplot(aes(color=interaction(root_association,precip,life_stage)))+
  scale_y_continuous(name="Inv Simpson")+theme_bw()




#simpson

bact_simpson_full_model= lm(log(InvSimpson)~soil_root_stage+as.factor(block), data= SILVA_MERDS_rar.divfil)
qqPlot(resid(bact_simpson_full_model))
hist(resid(bact_simpson_full_model))

Anova(bact_simpson_full_model, type=3)
#soil_root_stage   62.85  7   15.2104 6.502e-10 ***

emmeans(bact_simpson_full_model, pairwise~soil_root_stage, adjust ="fdr")
"$contrasts
 contrast              estimate    SE df t.ratio p.value
 L.B.G - S.B.G          1.30099 0.384 44  3.387  0.0027 
 L.B.G - L.R.G          0.09023 0.384 44  0.235  0.8456 
 L.B.G - L.B.S          1.17698 0.384 44  3.064  0.0058 
 L.B.G - S.B.S          1.75439 0.401 44  4.377  0.0002 
 L.B.G - L.R.S         -0.48820 0.386 44 -1.265  0.2480 
 L.B.G - L.B.Start     -1.59752 0.470 44 -3.395  0.0027 
 L.B.G - L.R.Start     -1.58906 0.470 44 -3.377  0.0027 
 S.B.G - L.R.G         -1.21076 0.384 44 -3.152  0.0048 
 S.B.G - L.B.S         -0.12401 0.384 44 -0.323  0.8059 
 S.B.G - S.B.S          0.45340 0.401 44  1.131  0.2958 
 S.B.G - L.R.S         -1.78918 0.386 44 -4.636  0.0001 
 S.B.G - L.B.Start     -2.89850 0.470 44 -6.160  <.0001 
 S.B.G - L.R.Start     -2.89005 0.470 44 -6.143  <.0001 
 L.R.G - L.B.S          1.08674 0.384 44  2.829  0.0103 
 L.R.G - S.B.S          1.66416 0.401 44  4.152  0.0004 
 L.R.G - L.R.S         -0.57843 0.386 44 -1.499  0.1796 
 L.R.G - L.B.Start     -1.68775 0.470 44 -3.587  0.0019 
 L.R.G - L.R.Start     -1.67930 0.470 44 -3.569  0.0019 
 L.B.S - S.B.S          0.57741 0.401 44  1.441  0.1908 
 L.B.S - L.R.S         -1.66517 0.386 44 -4.315  0.0002 
 L.B.S - L.B.Start     -2.77449 0.470 44 -5.897  <.0001 
 L.B.S - L.R.Start     -2.76604 0.470 44 -5.879  <.0001 
 S.B.S - L.R.S         -2.24258 0.406 44 -5.518  <.0001 
 S.B.S - L.B.Start     -3.35190 0.484 44 -6.923  <.0001 
 S.B.S - L.R.Start     -3.34345 0.484 44 -6.906  <.0001 
 L.R.S - L.B.Start     -1.10932 0.472 44 -2.350  0.0324 
 L.R.S - L.R.Start     -1.10087 0.472 44 -2.333  0.0324 
 L.B.Start - L.R.Start  0.00845 0.543 44  0.016  0.9877 

Results are averaged over the levels of: block 
Results are given on the log (not the response) scale. 
P value adjustment: fdr method for 28 tests "


ggplot(SILVA_MERDS_rar.divfil, aes(x=soil_status, y=InvSimpson))+geom_boxplot(aes(color=interaction(root_association,precip,life_stage)))+
  scale_y_continuous(name="Inv Simpson")+theme_bw()

#####Trans initial community####
#Let's just look at the initial community

SILVA_MERDS_rar.divfil_start=subset(SILVA_MERDS_rar.divfil, life_stage=="Start")

SILVA_MERDS_rar.divfil_start %>% group_by(soil_root) %>% summarise_at("Observed",funs(n(),mean,sd,se=sd(.)/sqrt(n())))

SILVA_MERDS_rar.divfil_start %>% group_by(soil_root) %>% summarise_at("InvSimpson",funs(n(),mean,sd,se=sd(.)/sqrt(n())))


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
Actino_soil_all.reads_start=subset(Actino_soil_all.reads,life_stage=="Start")
nrow(Actino_soil_all.reads_start)
#8
#Prop of Actinobacteria
Actino_model_start= lm(logit(Actino.prop)~soil_root+as.factor(block), data= Actino_soil_all.reads_start)
qqPlot(resid(Actino_model_start))
hist(resid(Actino_model_start))
boxCox(Actino_model_start)
shapiro.test(resid(Actino_model_start))
#0.0519
Anova(Actino_model_start, type=3)
#nada

emmeans(Actino_model_start, pairwise~soil_root)


emmeans(Actino_model_start, pairwise~soil_root*life_stage*precip)

Actino_model_start_no_block= lm(logit(Actino.prop)~soil_root, data= Actino_soil_all.reads_start)
qqPlot(resid(Actino_model_start_no_block))
hist(resid(Actino_model_start_no_block))
boxCox(Actino_model_start_no_block)
shapiro.test(resid(Actino_model_start_no_block))
#0.7783
Anova(Actino_model_start_no_block, type=3)
#nada

AIC(Actino_model_start,Actino_model_start_no_block)

Actino_soil_all.reads_start %>% group_by(soil_root) %>% summarise_at("Actino.prop",funs(n(),mean,sd,se=sd(.)/sqrt(n())))


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

Firm_soil_all.reads_start=subset(Firm_soil_all.reads,life_stage=="Start")
nrow(Firm_soil_all.reads_start)
#8
#Prop of Firmicutes
Firm_model_start= lm(logit(Firm.prop)~soil_root+as.factor(block), data= Firm_soil_all.reads_start)
qqPlot(resid(Firm_model_start))
hist(resid(Firm_model_start))
boxCox(Firm_model_start)
shapiro.test(resid(Firm_model_start))
#0.03837
Anova(Firm_model_start, type=3)
#nada

emmeans(Firm_model_start, pairwise~soil_root)


emmeans(Firm_model_start, pairwise~soil_root*life_stage*precip)

Firm_model_start_no_block= lm(logit(Firm.prop)~soil_root, data= Firm_soil_all.reads_start)
qqPlot(resid(Firm_model_start_no_block))
hist(resid(Firm_model_start_no_block))
boxCox(Firm_model_start_no_block)
shapiro.test(resid(Firm_model_start_no_block))
#0.354
Anova(Firm_model_start_no_block, type=3)
#soil_root     0.360  1    6.6726   0.04159 * 

AIC(Firm_model_start,Firm_model_start_no_block)

Firm_soil_all.reads_start %>% group_by(soil_root) %>% summarise_at("Firm.prop",funs(n(),mean,sd,se=sd(.)/sqrt(n())))

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

Chloro_soil_all.reads_start=subset(Chloro_soil_all.reads,life_stage=="Start")
nrow(Chloro_soil_all.reads_start)
#8
#Prop of Chloroflexi
Chloro_model_start= lm((Chloro.prop)~soil_root+as.factor(block), data= Chloro_soil_all.reads_start)
qqPlot(resid(Chloro_model_start))
hist(resid(Chloro_model_start))
boxCox(Chloro_model_start)
shapiro.test(resid(Chloro_model_start))
#0.6608
Anova(Chloro_model_start, type=3)
#nada

emmeans(Chloro_model_start, pairwise~soil_root)


emmeans(Chloro_model_start, pairwise~soil_root*life_stage*precip)

Chloro_model_start_no_block= lm((Chloro.prop)~soil_root, data= Chloro_soil_all.reads_start)
qqPlot(resid(Chloro_model_start_no_block))
hist(resid(Chloro_model_start_no_block))
boxCox(Chloro_model_start_no_block)
shapiro.test(resid(Chloro_model_start_no_block))
#0.09412
Anova(Chloro_model_start_no_block, type=3)
#nada

AIC(Chloro_model_start,Chloro_model_start_no_block)

Chloro_soil_all.reads_start %>% group_by(soil_root) %>% summarise_at("Chloro.prop",funs(n(),mean,sd,se=sd(.)/sqrt(n())))

#let's just look at the treatments
SILVA_MERDS_rar.divfil_trt=subset(SILVA_MERDS_rar.divfil, life_stage!="Start")
nrow(SILVA_MERDS_rar.divfil_trt)
#47
bact_obs_rich_model= lm(log(Observed)~life_stage*soil_root*precip+as.factor(block), data= SILVA_MERDS_rar.divfil_trt)
qqPlot(resid(bact_obs_rich_model))
hist(resid(bact_obs_rich_model))

Anova(bact_obs_rich_model, type=3)
#life_stage                     0.25  1     7.0966 0.01199 *  
#soil_root                     34.17  2   492.6814 < 2e-16 ***
#precip                         0.26  1     7.4928 0.01003 *  
#life_stage:soil_root           0.28  2     3.9931 0.02830 *  


emmeans(bact_obs_rich_model, pairwise~soil_root|life_stage)
#$contrasts
#life_stage = G:
#  contrast  estimate     SE df t.ratio p.value
#L.B - S.B   1.8081 0.0931 32  19.419 <.0001 
#L.B - L.R   0.0405 0.0931 32   0.435 0.9012 
#S.B - L.R  -1.7675 0.0931 32 -18.983 <.0001 

#life_stage = S:
#  contrast  estimate     SE df t.ratio p.value
#L.B - S.B   1.7934 0.0982 32  18.262 <.0001 
#L.B - L.R  -0.2940 0.0936 32  -3.139 0.0099 
#S.B - L.R  -2.0874 0.1001 32 -20.856 <.0001 




SILVA_MERDS_rar.divfil_trt_stage_soil_g=SILVA_MERDS_rar.divfil_trt %>% group_by(soil_root,life_stage)
obs_rich_stage=summarise_at(SILVA_MERDS_rar.divfil_trt_stage_soil_g, 
                            "Observed", funs(n(),mean,sd,se=sd(.)/sqrt(n())))

ggplot(obs_rich_stage, aes(x=life_stage,y=mean,ymin = mean-se, ymax= mean+se,fill=soil_root))+
  geom_bar(position="dodge",stat="identity", color="black",aes(fill=soil_root))+
  geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+
  scale_fill_manual(values = c("gray", "white", "black"),labels=c("Bulk","Sterile","Rhizo"))+
  scale_x_discrete(labels=c("Transplant","Seed"))+ylab("Bacterial richness")+
  geom_text(aes(y=mean+se+20, label=n),position=position_dodge(width=0.9))+theme_bw()

ggplot(SILVA_MERDS_rar.divfil_trt, aes(x=precip, y=Observed))+geom_boxplot(aes(color=precip))+theme_bw()




#Simpson
bact_inv_simp_model= lm(sqrt(InvSimpson)~life_stage*soil_root*precip+as.factor(block), data= SILVA_MERDS_rar.divfil_trt)
qqPlot(resid(bact_inv_simp_model))
hist(resid(bact_inv_simp_model))

Anova(bact_inv_simp_model, type=3,singular.ok =T)
#soil_root                    210.09  2  13.9318 4.444e-05 ***
#life_stage:soil_root          52.10  2   3.4549    0.0438 *   

emmeans(bact_inv_simp_model, pairwise~soil_root|life_stage)
#$contrasts
#life_stage = G:
#  contrast  estimate   SE df t.ratio p.value
#L.B - S.B    3.889 1.37 32  2.833  0.0210 
#L.B - L.R   -0.285 1.37 32 -0.207  0.9766 
#S.B - L.R   -4.174 1.37 32 -3.040  0.0127 

#life_stage = S:
#  contrast  estimate   SE df t.ratio p.value
#L.B - S.B    1.052 1.45 32  0.726  0.7498 
#L.B - L.R   -5.389 1.38 32 -3.903  0.0013 
#S.B - L.R   -6.441 1.48 32 -4.364  0.0004

SILVA_MERDS_rar.divfil_trt_stage_soil_g=SILVA_MERDS_rar.divfil_trt %>% group_by(soil_root,life_stage)
inv_simp_stage=summarise_at(SILVA_MERDS_rar.divfil_trt_stage_soil_g, 
                                   "InvSimpson", funs(n(),mean,sd,se=sd(.)/sqrt(n())))

ggplot(inv_simp_stage, aes(x=life_stage,y=mean,ymin = mean-se, ymax= mean+se,fill=soil_root))+
  geom_bar(position="dodge",stat="identity", color="black",aes(fill=soil_root))+
  geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+
  scale_fill_manual(values = c("gray", "white", "black"),labels=c("Bulk","Sterile","Rhizo"))+
  scale_x_discrete(labels=c("Transplant","Seed"))+ylab("Inv Simpson bacteria")+
  geom_text(aes(y=mean+se+5, label=n),position=position_dodge(width=0.9))+theme_bw()


#####Transplants####

SILVA_MERDS_rar_trt_trans=subset_samples(SILVA_MERDS_rar, life_stage=="G"&life_stage!="Start")
nrow(sample_data(SILVA_MERDS_rar_trt_trans))
#24


#####Bray Live and Sterile Trans Comm Analyses####
SILVA_MERDS_rar_trt_trans_ord=ordinate(SILVA_MERDS_rar_trt_trans, method = "NMDS",distance = "bray")
#*** Solution reached
#Warning message:
#In metaMDS(veganifyOTU(physeq), distance, ...) :
#stress is (nearly) zero: you may have insufficient data
#7.561319e-05
plot_ordination(SILVA_MERDS_rar_trt_trans,SILVA_MERDS_rar_trt_trans_ord, color="precip",shape="soil_root")+geom_point(size=4)+
  geom_label_repel(size=3,aes(label = block))+theme_bw()


SILVA_MERDS_rar_trt_trans_map=sample_data(SILVA_MERDS_rar_trt_trans)
SILVA_MERDS_rar_trt_trans_map$soil_root=with(SILVA_MERDS_rar_trt_trans_map, interaction(soil_status,root_association))
SILVA_MERDS_rar_trt_trans_dis=distance(SILVA_MERDS_rar_trt_trans,method = "bray")

adonis(SILVA_MERDS_rar_trt_trans_dis~SILVA_MERDS_rar_trt_trans_map$soil_root*SILVA_MERDS_rar_trt_trans_map$precip
       +as.factor(SILVA_MERDS_rar_trt_trans_map$block), permutations = 9999)
#SILVA_MERDS_rar_trt_trans_map$soil_root                                       2    2.6421 1.32104 10.9459 0.47836 0.0001 ***
#SILVA_MERDS_rar_trt_trans_map$precip                                          1    0.2382 0.23816  1.9734 0.04312 0.0874 .  





#####Jaccard Live and Sterile Trans Comm Analyses####

SILVA_MERDS_rar_trt_trans_J_ord=ordinate(SILVA_MERDS_rar_trt_trans, method = "NMDS", distance = "jaccard", binary = TRUE)
#Warning message:
#In metaMDS(veganifyOTU(physeq), distance, ...) :
#  stress is (nearly) zero: you may have insufficient data
#7.661321e-05 
plot_ordination(SILVA_MERDS_rar_trt_trans,SILVA_MERDS_rar_trt_trans_J_ord, color="precip",shape="soil_root")+geom_point(size=4)+
  geom_label_repel(size=3,aes(label = block))+theme_bw()


SILVA_MERDS_rar_trt_trans_map=sample_data(SILVA_MERDS_rar_trt_trans)
SILVA_MERDS_rar_trt_trans_map$soil_root=with(SILVA_MERDS_rar_trt_trans_map, interaction(soil_status,root_association))
SILVA_MERDS_rar_trt_trans_J_dis=distance(SILVA_MERDS_rar_trt_trans,method = "jaccard", binary = TRUE)

adonis(SILVA_MERDS_rar_trt_trans_J_dis~SILVA_MERDS_rar_trt_trans_map$soil_root*SILVA_MERDS_rar_trt_trans_map$precip
       +as.factor(SILVA_MERDS_rar_trt_trans_map$block), permutations = 9999)
#SILVA_MERDS_rar_trt_trans_map$soil_root                                       2    2.5660 1.28299  5.6229 0.33357 0.0001 ***
#SILVA_MERDS_rar_trt_trans_map$precip                                          1    0.3559 0.35589  1.5598 0.04626 0.0968 . 

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




#####UnWeighted Unifrac Live and Sterile Trans Comm Analyses####

SILVA_MERDS_rar_trt_trans_unWU_ord=ordinate(SILVA_MERDS_rar_trt_trans, method = "NMDS",distance = "unifrac")
#*** Solution reached
#Warning message:
#In metaMDS(ps.dist) :
#  stress is (nearly) zero: you may have insufficient data
#8.884679e-05
plot_ordination(SILVA_MERDS_rar_trt_trans,SILVA_MERDS_rar_trt_trans_unWU_ord, color="precip",shape="soil_root")+geom_point(size=4)+
  geom_label_repel(size=3,aes(label = block))+theme_bw()


SILVA_MERDS_rar_trt_trans_map=sample_data(SILVA_MERDS_rar_trt_trans)
SILVA_MERDS_rar_trt_trans_map$soil_root=with(SILVA_MERDS_rar_trt_trans_map, interaction(soil_status,root_association))
SILVA_MERDS_rar_trt_trans_unWU_dis=distance(SILVA_MERDS_rar_trt_trans,method = "unifrac")

adonis(SILVA_MERDS_rar_trt_trans_unWU_dis~SILVA_MERDS_rar_trt_trans_map$soil_root*SILVA_MERDS_rar_trt_trans_map$precip
       +as.factor(SILVA_MERDS_rar_trt_trans_map$block), permutations = 9999)
#SILVA_MERDS_rar_trt_trans_map$soil_root                                       2    2.4418 1.22089  8.6757 0.42595 0.0001 ***
#SILVA_MERDS_rar_trt_trans_map$precip                                          1    0.2791 0.27909  1.9833 0.04869 0.0582 .  


#####Presence UnWeighted Unifrac Transplant####
head(sample_data(SILVA_MERDS_rar))

SILVA_MERDS_rar_pres_trans=subset_samples(SILVA_MERDS_rar, life_stage=="G"&life_stage!="Start"&
                                            root_association=="B")
nrow(sample_data(SILVA_MERDS_rar_pres_trans))
#16
unique(sample_data(SILVA_MERDS_rar_pres_trans)$soil_root)

SILVA_MERDS_rar_pres_trans_map=sample_data(SILVA_MERDS_rar_pres_trans)
SILVA_MERDS_rar_pres_trans_map$soil_root=with(SILVA_MERDS_rar_pres_trans_map, interaction(soil_status,root_association))
SILVA_MERDS_rar_pres_trans_UnWU_dis=distance(SILVA_MERDS_rar_pres_trans,method = "unifrac")

(pres_tran_UnWU=adonis(SILVA_MERDS_rar_pres_trans_UnWU_dis~SILVA_MERDS_rar_pres_trans_map$soil_root*SILVA_MERDS_rar_pres_trans_map$precip
                     +as.factor(SILVA_MERDS_rar_pres_trans_map$block), permutations = 9999))


AICc.PERMANOVA(pres_tran_UnWU)
#$AIC
#[1] 17.33034


pairwise.perm.manova(SILVA_MERDS_rar_pres_trans_UnWU_dis, SILVA_MERDS_rar_map$soil_root_stage, nperm=2000)


(pres_tran_UnWU_no_block=adonis(SILVA_MERDS_rar_pres_trans_UnWU_dis~SILVA_MERDS_rar_pres_trans_map$soil_root*
                                SILVA_MERDS_rar_pres_trans_map$precip, permutations = 9999))


AICc.PERMANOVA(pres_tran_UnWU_no_block)
#$AIC
#[1] 16.20397


#####Origin UnWeighted Unifrac Transplant####
head(sample_data(SILVA_MERDS_rar))

SILVA_MERDS_rar_orig_trans=subset_samples(SILVA_MERDS_rar, life_stage=="G"&life_stage!="Start"&
                                            soil_status =="L")
nrow(sample_data(SILVA_MERDS_rar_orig_trans))
#16
unique(sample_data(SILVA_MERDS_rar_orig_trans)$soil_root)

SILVA_MERDS_rar_orig_trans_map=sample_data(SILVA_MERDS_rar_orig_trans)
SILVA_MERDS_rar_orig_trans_map$soil_root=with(SILVA_MERDS_rar_orig_trans_map, interaction(soil_status,root_association))
SILVA_MERDS_rar_orig_trans_UnWU_dis=distance(SILVA_MERDS_rar_orig_trans,method = "unifrac")

(orig_tran_UnWU=adonis(SILVA_MERDS_rar_orig_trans_UnWU_dis~SILVA_MERDS_rar_orig_trans_map$soil_root*SILVA_MERDS_rar_orig_trans_map$precip
                     +as.factor(SILVA_MERDS_rar_orig_trans_map$block), permutations = 9999))


AICc.PERMANOVA(orig_tran_UnWU)
#$AIC
#[1] 17.73641


pairwise.perm.manova(SILVA_MERDS_rar_orig_trans_UnWU_dis, SILVA_MERDS_rar_map$soil_root_stage, nperm=2000)


(orig_tran_UnWU_no_block=adonis(SILVA_MERDS_rar_orig_trans_UnWU_dis~SILVA_MERDS_rar_orig_trans_map$soil_root*
                                SILVA_MERDS_rar_orig_trans_map$precip, permutations = 9999))


AICc.PERMANOVA(orig_tran_UnWU_no_block)
#$AIC
#[1] 17.10015




#####Origin DESeq2####

unique(sample_data(SILVA_MERDS_rar)$soil_root)
SILVA_MERDS_rar_orig_trans_bulk=subset_samples(SILVA_MERDS_rar, life_stage=="G"&life_stage!="Start"&
                                            soil_root=="L.B")
SILVA_MERDS_rar_orig_trans_bulk=prune_taxa(taxa_sums(SILVA_MERDS_rar_orig_trans_bulk) > 0, SILVA_MERDS_rar_orig_trans_bulk)
nsamples(SILVA_MERDS_rar_orig_trans_bulk)
#8

orig_trans_bulk_ds2 <- phyloseq_to_deseq2(SILVA_MERDS_rar_orig_trans_bulk, ~ precip)
orig_trans_bulk_dds <- DESeq(orig_trans_bulk_ds2)
orig_trans_bulk_res <- results(orig_trans_bulk_dds)
orig_trans_bulk_df <- as.data.frame(orig_trans_bulk_res)
orig_trans_bulk_df$taxon <- rownames(orig_trans_bulk_df)
orig_trans_bulk_df <- orig_trans_bulk_df %>% arrange(log2FoldChange, padj)

library(knitr)
print(head(kable((orig_trans_bulk_df))))


SILVA_MERDS_rar_orig_trans_rhizo=subset_samples(SILVA_MERDS_rar, life_stage=="G"&life_stage!="Start"&
                                                 soil_root=="L.R")
SILVA_MERDS_rar_orig_trans_rhizo=prune_taxa(taxa_sums(SILVA_MERDS_rar_orig_trans_rhizo) > 0, SILVA_MERDS_rar_orig_trans_rhizo)
nsamples(SILVA_MERDS_rar_orig_trans_rhizo)
#8

orig_trans_rhizo_ds2 <- phyloseq_to_deseq2(SILVA_MERDS_rar_orig_trans_rhizo, ~ precip)
orig_trans_rhizo_dds <- DESeq(orig_trans_rhizo_ds2)
orig_trans_rhizo_res <- results(orig_trans_rhizo_dds)
orig_trans_rhizo_df <- as.data.frame(orig_trans_rhizo_res)
orig_trans_rhizo_df$taxon <- rownames(orig_trans_rhizo_df)
orig_trans_rhizo_df <- orig_trans_rhizo_df %>% arrange(log2FoldChange, padj)

library(knitr)
print(head(kable((orig_trans_rhizo_df))))




#Stack bar graphs
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

#####Using GTDB for classification####
#Make a new Phyloseq obj with all of the associated data
GTDB_MERDS_rar_trt_trans=phyloseq(otu_table(SILVA_MERDS_rar_trt_trans),TAXA_GTDBr89_MERDs_confid,
                                   SILVA_MERDS_rar_trt_trans_map)
ntaxa(GTDB_MERDS_rar_trt_trans)
#8366
#merge OTUs by the soil and precipitation treatment
GTDB_MERDS_rar_trt_trans_fact=merge_samples(GTDB_MERDS_rar_trt_trans, "soil_root_precip")
sample_names(GTDB_MERDS_rar_trt_trans_fact)     

#combine the reads at Phylum level
get_taxa_unique(GTDB_MERDS_rar_trt_trans, taxonomic.rank="Phylum")
#42
(GTDB_MERDS_rar_trt_trans_fact.phylum<-tax_glom(GTDB_MERDS_rar_trt_trans_fact, taxrank="Phylum"))



#subset so there is only the top ten most abundant phyla
GTDB_TopPHYL_trans = names(sort(taxa_sums(GTDB_MERDS_rar_trt_trans_fact.phylum), TRUE)[1:10])
GTDB_MERDS_rar_trt_trans.T10 = prune_taxa(GTDB_TopPHYL_trans, GTDB_MERDS_rar_trt_trans_fact.phylum)
GTDB_trans_name_T10=get_taxa_unique(GTDB_MERDS_rar_trt_trans.T10, taxonomic.rank="Phylum")
GTDB_trans_name_T10_sep <- colsplit(GTDB_trans_name_T10, "__", c("letter", "Phyl_name"))


#transform the read counts to prop of total reads

GTDB_MERDS_rar_trt_trans_fact.phylum.prop=transform_sample_counts(GTDB_MERDS_rar_trt_trans_fact.phylum, function(x)x/sum(x))

taxon_positions=c("S.B.A","L.B.A", "L.R.A", "S.B.D","L.B.D" ,"L.R.D")
GTDB_MERDS_rar_trt_trans_10_prop = prune_taxa(GTDB_TopPHYL_trans, GTDB_MERDS_rar_trt_trans_fact.phylum.prop)

GTDB_MERDS_rar_trt_trans_10_prop_otu=as.data.frame(t(otu_table(GTDB_MERDS_rar_trt_trans_10_prop)))

#create an other taxa category
GTDB_taxon_sums_trans=c(sum(GTDB_MERDS_rar_trt_trans_10_prop_otu$L.B.A),sum(GTDB_MERDS_rar_trt_trans_10_prop_otu$S.B.A),
                         sum(GTDB_MERDS_rar_trt_trans_10_prop_otu$L.R.A),sum(GTDB_MERDS_rar_trt_trans_10_prop_otu$L.B.D),
                         sum(GTDB_MERDS_rar_trt_trans_10_prop_otu$S.B.D),sum(GTDB_MERDS_rar_trt_trans_10_prop_otu$L.R.D))
GTDB_tran_other_spp=c(as.numeric(1-GTDB_taxon_sums_trans))
GTDB_MERDS_rar_trt_trans_10_prop_OTU=rbind(GTDB_MERDS_rar_trt_trans_10_prop_otu,GTDB_tran_other_spp)
summary(GTDB_MERDS_rar_trt_trans_10_prop_OTU)
GTDB_MERDS_rar_trt_trans_10_prop_OTU[,"Phylum"]=c(as.character(GTDB_trans_name_T10_sep$Phyl_name),"Other Phyla")
summary(GTDB_MERDS_rar_trt_trans_10_prop_OTU)
GTDB_MERDS_rar_trt_trans_10_prop_OTU_M=melt(GTDB_MERDS_rar_trt_trans_10_prop_OTU,id="Phylum")
phyl_order=c(sort(as.character(GTDB_trans_name_T10_sep$Phyl_name)),"Other Phyla")
summary(GTDB_MERDS_rar_trt_trans_10_prop_OTU_M)



(p_bact_T10_v2.1_color=ggplot(GTDB_MERDS_rar_trt_trans_10_prop_OTU_M,aes(x=variable,y=value,fill=factor(Phylum, levels=phyl_order)))+
    geom_bar(aes( fill=factor(Phylum, levels=phyl_order)), stat="identity", position="stack",color="black")+theme_bw()+
    theme(axis.text.y=element_text(size=18),axis.text.x=element_text(size=18),
          axis.title=element_text(size=20),panel.grid.major=element_blank(),legend.text = element_text(size=16), legend.title = element_text(size=20),
          panel.grid.minor=element_blank())+xlab(NULL)+ylab("Proportion")+
    scale_x_discrete(limits = taxon_positions,labels=c("Sterile\nAmbient","Bulk\nAmbient", 
                                                       "Rhizo\nAmbient",
                                                       "Sterile\nDrought",
                                                       "Bulk\nDrought", 
                                                       "Rhizo\nDrought"))+scale_fill_brewer(palette="Paired")+
    guides(fill=guide_legend(title="Phyla"))+ ggtitle(label = "Genome Taxonomy DataBase"))
#####


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


#####Using RDP for classification####
#Make a new Phyloseq obj with all of the associated data
RDP_MERDS_rar_trt_trans=phyloseq(otu_table(SILVA_MERDS_rar_trt_trans),tax_table(phyl_RDP_MERDS),
                                  SILVA_MERDS_rar_trt_trans_map)
ntaxa(RDP_MERDS_rar_trt_trans)
#8366
#merge OTUs by the soil and precipitation treatment
RDP_MERDS_rar_trt_trans_fact=merge_samples(RDP_MERDS_rar_trt_trans, "soil_root_precip")
sample_names(RDP_MERDS_rar_trt_trans_fact)     

#combine the reads at Phylum level
get_taxa_unique(RDP_MERDS_rar_trt_trans, taxonomic.rank="Phylum")
#30
(RDP_MERDS_rar_trt_trans_fact.phylum<-tax_glom(RDP_MERDS_rar_trt_trans_fact, taxrank="Phylum"))



#subset so there is only the top ten most abundant phyla
RDP_TopPHYL_trans = names(sort(taxa_sums(RDP_MERDS_rar_trt_trans_fact.phylum), TRUE)[1:10])
RDP_MERDS_rar_trt_trans.T10 = prune_taxa(RDP_TopPHYL_trans, RDP_MERDS_rar_trt_trans_fact.phylum)
RDP_trans_name_T10=get_taxa_unique(RDP_MERDS_rar_trt_trans.T10, taxonomic.rank="Phylum")
RDP_trans_name_T10_sep <- colsplit(RDP_trans_name_T10, ":", c("letter", "Phyl_name"))


#transform the read counts to prop of total reads

RDP_MERDS_rar_trt_trans_fact.phylum.prop=transform_sample_counts(RDP_MERDS_rar_trt_trans_fact.phylum, function(x)x/sum(x))

taxon_positions=c("S.B.A","L.B.A", "L.R.A", "S.B.D","L.B.D" ,"L.R.D")
RDP_MERDS_rar_trt_trans_10_prop = prune_taxa(RDP_TopPHYL_trans, RDP_MERDS_rar_trt_trans_fact.phylum.prop)

RDP_MERDS_rar_trt_trans_10_prop_otu=as.data.frame(t(otu_table(RDP_MERDS_rar_trt_trans_10_prop)))

#create an other taxa category
RDP_taxon_sums_trans=c(sum(RDP_MERDS_rar_trt_trans_10_prop_otu$L.B.A),sum(RDP_MERDS_rar_trt_trans_10_prop_otu$S.B.A),
                        sum(RDP_MERDS_rar_trt_trans_10_prop_otu$L.R.A),sum(RDP_MERDS_rar_trt_trans_10_prop_otu$L.B.D),
                        sum(RDP_MERDS_rar_trt_trans_10_prop_otu$S.B.D),sum(RDP_MERDS_rar_trt_trans_10_prop_otu$L.R.D))
RDP_tran_other_spp=c(as.numeric(1-RDP_taxon_sums_trans))
RDP_MERDS_rar_trt_trans_10_prop_OTU=rbind(RDP_MERDS_rar_trt_trans_10_prop_otu,RDP_tran_other_spp)
summary(RDP_MERDS_rar_trt_trans_10_prop_OTU)
RDP_MERDS_rar_trt_trans_10_prop_OTU[,"Phylum"]=c(as.character(RDP_trans_name_T10_sep$Phyl_name),"Other Phyla")
summary(RDP_MERDS_rar_trt_trans_10_prop_OTU)
RDP_MERDS_rar_trt_trans_10_prop_OTU_M=melt(RDP_MERDS_rar_trt_trans_10_prop_OTU,id="Phylum")
phyl_order=c(sort(as.character(RDP_trans_name_T10_sep$Phyl_name)),"Other Phyla")
summary(RDP_MERDS_rar_trt_trans_10_prop_OTU_M)



(p_bact_T10_v2.1_color=ggplot(RDP_MERDS_rar_trt_trans_10_prop_OTU_M,aes(x=variable,y=value,fill=factor(Phylum, levels=phyl_order)))+
    geom_bar(aes( fill=factor(Phylum, levels=phyl_order)), stat="identity", position="stack",color="black")+theme_bw()+
    theme(axis.text.y=element_text(size=18),axis.text.x=element_text(size=18),
          axis.title=element_text(size=20),panel.grid.major=element_blank(),legend.text = element_text(size=16), legend.title = element_text(size=20),
          panel.grid.minor=element_blank())+xlab(NULL)+ylab("Proportion")+
    scale_x_discrete(limits = taxon_positions,labels=c("Sterile\nAmbient","Bulk\nAmbient", 
                                                       "Rhizo\nAmbient",
                                                       "Sterile\nDrought",
                                                       "Bulk\nDrought", 
                                                       "Rhizo\nDrought"))+scale_fill_brewer(palette="Paired")+
    guides(fill=guide_legend(title="Phyla"))+ ggtitle(label = "RDP DataBase"))
#####


#####Transplant live Community analyses####
SILVA_MERDS_rar_trt_trans=subset_samples(SILVA_MERDS_rar, life_stage=="G"&life_stage!="Start")
#No sterile
SILVA_MERDS_rar_trt_live_trans=subset_samples(SILVA_MERDS_rar_trt_trans, soil_root!="S.B")
nrow(sample_data(SILVA_MERDS_rar_trt_live_trans))
#16

#####Bray Live Trans Comm Analyses####

SILVA_MERDS_rar_trt_live_trans_ord=ordinate(SILVA_MERDS_rar_trt_live_trans, method = "NMDS",distance = "bray")
#*** No convergence -- monoMDS stopping criteria:
#20: stress ratio > sratmax
#0.1607911
plot_ordination(SILVA_MERDS_rar_trt_live_trans,SILVA_MERDS_rar_trt_live_trans_ord, color="precip",shape="root_association")+geom_point(size=4)+
  geom_label_repel(size=3,aes(label = block))+theme_bw()


SILVA_MERDS_rar_trt_live_trans_map=sample_data(SILVA_MERDS_rar_trt_live_trans)
SILVA_MERDS_rar_trt_live_trans_map$soil_root=with(SILVA_MERDS_rar_trt_live_trans_map, interaction(soil_status,root_association))
SILVA_MERDS_rar_trt_live_transe_dis=distance(SILVA_MERDS_rar_trt_live_trans,method = "bray")

adonis(SILVA_MERDS_rar_trt_live_transe_dis~SILVA_MERDS_rar_trt_live_trans_map$soil_root*SILVA_MERDS_rar_trt_live_trans_map$precip
       +as.factor(SILVA_MERDS_rar_trt_live_trans_map$block), permutations = 9999)
#SILVA_MERDS_rar_trt_live_trans_map$precip                                               1   0.26218 0.26218 2.13979 0.12422 0.0045 **
#as.factor(SILVA_MERDS_rar_trt_live_trans_map$block)                                     3   0.47584 0.15861 1.29453 0.22545 0.0808 . 

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


#####Jaccard Live Trans Comm Analyses####

SILVA_MERDS_rar_trt_live_trans_J_ord=ordinate(SILVA_MERDS_rar_trt_live_trans, method = "NMDS", distance = "jaccard", binary = TRUE)
#*** Solution reached
#0.1554819
plot_ordination(SILVA_MERDS_rar_trt_live_trans,SILVA_MERDS_rar_trt_live_trans_J_ord, color="precip",shape="root_association")+geom_point(size=4)+
  geom_label_repel(size=3,aes(label = block))+theme_bw()


SILVA_MERDS_rar_trt_live_trans_map=sample_data(SILVA_MERDS_rar_trt_live_trans)
SILVA_MERDS_rar_trt_live_trans_map$soil_root=with(SILVA_MERDS_rar_trt_live_trans_map, interaction(soil_status,root_association))
SILVA_MERDS_rar_trt_live_trans_J__dis=distance(SILVA_MERDS_rar_trt_live_trans,method = "jaccard", binary = TRUE)

adonis(SILVA_MERDS_rar_trt_live_trans_J__dis~SILVA_MERDS_rar_trt_live_trans_map$soil_root*SILVA_MERDS_rar_trt_live_trans_map$precip
       +as.factor(SILVA_MERDS_rar_trt_live_trans_map$block), permutations = 9999)
#SILVA_MERDS_rar_trt_live_trans_map$precip                                               1    0.3355 0.33546  1.4496 0.08996 0.0018 **
#as.factor(SILVA_MERDS_rar_trt_live_trans_map$block)                                     3    0.8458 0.28193  1.2183 0.22683 0.0036 **

#####Weighted Unifrac Live Trans Comm Analyses####

SILVA_MERDS_rar_trt_live_trans_WU_ord=ordinate(SILVA_MERDS_rar_trt_live_trans, method = "NMDS", distance = "wunifrac")
#*** Solution reached
#0.06950258
plot_ordination(SILVA_MERDS_rar_trt_live_trans,SILVA_MERDS_rar_trt_live_trans_WU_ord, color="precip",shape="root_association")+geom_point(size=4)+
  geom_label_repel(size=3,aes(label = block))+theme_bw()

SILVA_MERDS_rar_trt_live_trans_WU_ord$
SILVA_MERDS_rar_trt_live_trans_map=sample_data(SILVA_MERDS_rar_trt_live_trans)
SILVA_MERDS_rar_trt_live_trans_map$soil_root=with(SILVA_MERDS_rar_trt_live_trans_map, interaction(soil_status,root_association))
SILVA_MERDS_rar_trt_live_trans_WU_dis=distance(SILVA_MERDS_rar_trt_live_trans,method = "wunifrac")

adonis(SILVA_MERDS_rar_trt_live_trans_WU_dis~SILVA_MERDS_rar_trt_live_trans_map$soil_root*SILVA_MERDS_rar_trt_live_trans_map$precip
       +as.factor(SILVA_MERDS_rar_trt_live_trans_map$block), permutations = 9999)
#SILVA_MERDS_rar_trt_live_trans_map$precip                                               1  0.031179 0.031179 2.53940 0.15160 0.0444 *
#SILVA_MERDS_rar_trt_live_trans_map$soil_root:SILVA_MERDS_rar_trt_live_trans_map$precip  1  0.024895 0.024895 2.02760 0.12105 0.0864 .


#####UnWeighted Unifrac Live Trans Comm Analyses####

SILVA_MERDS_rar_trt_live_trans_unWU_ord=ordinate(SILVA_MERDS_rar_trt_live_trans, method = "NMDS", distance = "unifrac")
#*** Solution reached
#0.1692574
plot_ordination(SILVA_MERDS_rar_trt_live_trans,SILVA_MERDS_rar_trt_live_trans_unWU_ord, color="precip",shape="root_association")+geom_point(size=4)+
  geom_label_repel(size=3,aes(label = block))+theme_bw()


SILVA_MERDS_rar_trt_live_trans_map=sample_data(SILVA_MERDS_rar_trt_live_trans)
SILVA_MERDS_rar_trt_live_trans_map$soil_root=with(SILVA_MERDS_rar_trt_live_trans_map, interaction(soil_status,root_association))
SILVA_MERDS_rar_trt_live_trans_unWU_dis=distance(SILVA_MERDS_rar_trt_live_trans,method = "unifrac")

adonis(SILVA_MERDS_rar_trt_live_trans_unWU_dis~SILVA_MERDS_rar_trt_live_trans_map$soil_root*SILVA_MERDS_rar_trt_live_trans_map$precip
       +as.factor(SILVA_MERDS_rar_trt_live_trans_map$block), permutations = 9999)
#SILVA_MERDS_rar_trt_live_trans_map$precip                                               1   0.27333 0.27333  1.9477 0.11746 0.0001 ***
#as.factor(SILVA_MERDS_rar_trt_live_trans_map$block)                                     3   0.50303 0.16768  1.1948 0.21616 0.0349 *  


#Stack bar graphs
SILVA_MERDS_rar_trt_live_trans_map=sample_data(SILVA_MERDS_rar_trt_live_trans)
SILVA_MERDS_rar_trt_live_trans_map$soil_root_precip=with(SILVA_MERDS_rar_trt_live_trans_map, interaction(soil_root, precip))

#Make a new Phyloseq obj with all of the associated data
SILVA_MERDS_rar_trt_live_trans=phyloseq(otu_table(SILVA_MERDS_rar_trt_live_trans),tax_table(SILVA_MERDS_rar_trt_live_trans),
                                 SILVA_MERDS_rar_trt_live_trans_map)
#merge OTUs by the soil and precipitation treatment
SILVA_MERDS_rar_trt_live_trans_fact=merge_samples(SILVA_MERDS_rar_trt_live_trans, "soil_root_precip")
sample_names(SILVA_MERDS_rar_trt_live_trans_fact)     

#combine the reads at Phylum level
get_taxa_unique(SILVA_MERDS_rar_trt_live_trans, taxonomic.rank="Phylum")
#38
(SILVA_MERDS_rar_trt_live_trans_fact.phylum<-tax_glom(SILVA_MERDS_rar_trt_live_trans_fact, taxrank="Phylum"))



#subset so there is only the top ten most abundant phyla
SILVA_TopPHYL_trans_live = names(sort(taxa_sums(SILVA_MERDS_rar_trt_live_trans_fact.phylum), TRUE)[1:10])

SILVA_trans_name_T10=get_taxa_unique(SILVA_trans_T10, taxonomic.rank="Phylum")
SILVA_trans_name_T10_sep <- colsplit(SILVA_trans_name_T10, ":", c("letter", "Phyl_name"))


#transform the read counts to prop of total reads

SILVA_MERDS_rar_trt_live_trans_fact.phylum.prop=transform_sample_counts(SILVA_MERDS_rar_trt_live_trans_fact.phylum, function(x)x/sum(x))

taxon_positions_live=c("L.B.A", "L.R.A", "L.B.D" ,"L.R.D")
SILVA_MERDS_rar_trt_live_trans_10_prop = prune_taxa(SILVA_TopPHYL_trans_live, SILVA_MERDS_rar_trt_live_trans_fact.phylum.prop)

SILVA_MERDS_rar_trt_live_trans_10_prop_otu=as.data.frame(t(otu_table(SILVA_MERDS_rar_trt_live_trans_10_prop)))

#create an other taxa category
SILVA_taxon_sums_L_trans=c(sum(SILVA_MERDS_rar_trt_live_trans_10_prop_otu$L.B.A),sum(SILVA_MERDS_rar_trt_live_trans_10_prop_otu$L.R.A),
             sum(SILVA_MERDS_rar_trt_live_trans_10_prop_otu$L.B.D),sum(SILVA_MERDS_rar_trt_live_trans_10_prop_otu$L.R.D))
SILVA_L_tran_other_spp=c(as.numeric(1-SILVA_taxon_sums_L_trans))
SILVA_MERDS_rar_trt_live_trans_10_prop_OTU=rbind(SILVA_MERDS_rar_trt_live_trans_10_prop_otu,SILVA_L_tran_other_spp)
summary(SILVA_MERDS_rar_trt_live_trans_10_prop_OTU)
SILVA_MERDS_rar_trt_live_trans_10_prop_OTU[,"Phylum"]=c(as.character(SILVA_trans_name_T10_sep$Phyl_name),"Other Phyla")
summary(SILVA_MERDS_rar_trt_live_trans_10_prop_OTU)
SILVA_MERDS_rar_trt_live_trans_10_prop_OTU_M=melt(SILVA_MERDS_rar_trt_live_trans_10_prop_OTU,id="Phylum")
phyl_order=c(sort(as.character(SILVA_trans_name_T10_sep$Phyl_name)),"Other Phyla")
summary(SILVA_MERDS_rar_trt_live_trans_10_prop_OTU_M)



(p_bact_T10_v2.1_color=ggplot(SILVA_MERDS_rar_trt_live_trans_10_prop_OTU_M,aes(x=variable,y=value,fill=factor(Phylum, levels=phyl_order)))+
    geom_bar(aes( fill=factor(Phylum, levels=phyl_order)), stat="identity", position="stack",color="black")+theme_bw()+
    theme(axis.text.y=element_text(size=18),axis.text.x=element_blank(),
          axis.title=element_text(size=20),panel.grid.major=element_blank(),legend.text = element_text(size=16), legend.title = element_text(size=20),
          panel.grid.minor=element_blank())+xlab(NULL)+ylab("Proportion")+
    scale_x_discrete(limits = taxon_positions_live,labels=c("Bulk\nAmbient", 
                                                       "Rhizosphere\nAmbient",
                                                       "Bulk\nDrought", 
                                                       "Rhizosphere\nDrought"))+scale_fill_brewer(palette="Paired")+
    guides(fill=guide_legend(title="Phyla")))






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


(trans_obs_richness_p2=ggplot(obs_rich_precip_trans, aes(x=precip,y=mean, ymin = mean-se, ymax= mean+se,fill=factor(soil_root,levels = treatment_order)))+
    geom_bar(position="dodge",stat="identity", color="black",aes(fill=factor(soil_root,levels = treatment_order)))+
    geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+ylim(c(0,1950))+
    scale_fill_manual(values = c( "white","lightgray", "darkgrey"),labels=c("Sterile","Bulk","Rhizo"), limits=treatment_order, name=NULL)+
    scale_x_discrete(labels=c("Ambient","Drought"), name=NULL)+ylab("Bacterial richness")+
    geom_text(aes(y=80, label=n),size=10,position=position_dodge(width=0.9))+theme_bw()+
    theme(axis.title = element_text(size = 23), axis.text = element_text(size = 23),
          legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank()))


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

(trans_inv_simpson_p2=ggplot(inv_simp_precip_trans, aes(x=precip,y=mean, ymin = mean-se, ymax= mean+se,fill=factor(soil_root,levels = treatment_order)))+
    geom_bar(position="dodge",stat="identity", color="black",aes(fill=factor(soil_root,levels = treatment_order)))+
    geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+ylim(c(-10,220))+
    scale_fill_manual(values = c( "white","lightgray", "darkgrey"),labels=c("Sterile","Bulk","Rhizo"), limits=treatment_order, name=NULL)+
    scale_x_discrete(labels=c("Ambient","Drought"), name=NULL)+ylab("inv Simpson bacteria")+
    geom_text(aes(y=-7, label=n),size=10,position=position_dodge(width=0.9))+theme_bw()+
    theme(axis.title = element_text(size = 23), axis.text = element_text(size = 23),
          legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank()))



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

#####Are shared OTUs dominant####

#Live only
SILVA_MERDS_rar_trt_live_trans_shared=core(SILVA_MERDS_rar_trt_live_trans, detection = 0,prevalence = 0.9999999999)
SILVA_MERDS_rar_trt_live_trans_shared_sampl_sum=sample_sums(SILVA_MERDS_rar_trt_live_trans_shared)
ntaxa(SILVA_MERDS_rar_trt_live_trans_shared)
#128
SILVA_MERDS_rar_trt_live_trans_map=sample_data(SILVA_MERDS_rar_trt_live_trans)
SILVA_MERDS_rar_trt_live_trans_map$soil_root=with(SILVA_MERDS_rar_trt_live_trans_map, interaction(soil_status,root_association))

shared_SILVA_MERDS_rar_trt_trans_sampl_sum=merge(SILVA_MERDS_rar_trt_live_trans_map,SILVA_MERDS_rar_trt_live_trans_shared_sampl_sum, by="row.names")


shared_SILVA_MERDS_rar_trt_trans_sampl_sum_precip_soil_g=shared_SILVA_MERDS_rar_trt_trans_sampl_sum %>% group_by(soil_root,precip)
shared_read_abund_precip_trans=summarise_at(shared_SILVA_MERDS_rar_trt_trans_sampl_sum_precip_soil_g, 
                                           "y", list(~n(),~mean,~sd,se=~sd(.)/sqrt(n())))

ggplot(shared_read_abund_precip_trans, aes(x=precip,y=mean,ymin = mean-se, ymax= mean+se,fill=soil_root))+
  geom_bar(position="dodge",stat="identity", color="black",aes(fill=soil_root))+
  geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+
  scale_fill_manual(values = c("light gray", "dark grey"),labels=c("Bulk","Rhizo"))+
  scale_x_discrete(labels=c("Ambient","Drought"))+ylab("Shared Taxa Read Abun")+
  geom_text(aes(y=mean+se+200, label=n),position=position_dodge(width=0.9))+theme_bw()


#Diversity of the Shared taxa
alpha_meas = c("Observed", "Shannon", "InvSimpson")
SILVA_MERDS_rar_trt_live_trans_map=sample_data(SILVA_MERDS_rar_trt_live_trans)
SILVA_MERDS_rar_trt_live_trans_map$soil_root=with(SILVA_MERDS_rar_trt_live_trans_map, interaction(soil_status,root_association))

SILVA_MERDS_rar_trt_live_trans_shared.divfil=estimate_richness(SILVA_MERDS_rar_trt_live_trans_shared,measures=alpha_meas)

SILVA_MERDS_rar_trt_live_trans_shared.divfil=merge(SILVA_MERDS_rar_trt_live_trans_shared.divfil, SILVA_MERDS_rar_trt_live_trans_map, by ="row.names")
#bact.soilE.t.divfil=mutate(bact.soilE.t.divfil, pielou=Shannon*(1/log(Observed)))
head(SILVA_MERDS_rar_trt_live_trans_shared.divfil)
row.names(SILVA_MERDS_rar_trt_live_trans_shared.divfil)=SILVA_MERDS_rar_trt_live_trans_shared.divfil$Row.names
SILVA_MERDS_rar_trt_live_trans_shared.divfil$Row.names=NULL



SILVA_MERDS_rar_trt_live_trans_shared.divfil_g=SILVA_MERDS_rar_trt_live_trans_shared.divfil %>% group_by(soil_root,precip)
shared_inv_simp_precip_trans=summarise_at(SILVA_MERDS_rar_trt_live_trans_shared.divfil_g, 
                                         "InvSimpson", list(~n(),~mean,~sd,se=~sd(.)/sqrt(n())))

ggplot(shared_inv_simp_precip_trans, aes(x=precip,y=mean,ymin = mean-se, ymax= mean+se,fill=soil_root))+
  geom_bar(position="dodge",stat="identity", color="black",aes(fill=soil_root))+
  geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+
  scale_fill_manual(values = c("light gray", "dark grey"),labels=c("Bulk","Rhizo"))+
  scale_x_discrete(labels=c("Ambient","Drought"))+ylab("Shared inv Simpson")+
  geom_text(aes(y=mean+se+1, label=n),position=position_dodge(width=0.9))+theme_bw()


#Found in at least one sample per treatment
SILVA_MERDS_rar_trt_live_trans_map=sample_data(SILVA_MERDS_rar_trt_live_trans)
SILVA_MERDS_rar_trt_live_trans_map$soil_root=with(SILVA_MERDS_rar_trt_live_trans_map, interaction(soil_status,root_association))
summary(SILVA_MERDS_rar_trt_live_trans_map)

SILVA_MERDS_rar_trt_live_trans_b_a=subset_samples(SILVA_MERDS_rar_trt_live_trans, root_association=="B"&precip=="A")
SILVA_MERDS_rar_trt_live_trans_b_a=prune_taxa(taxa_sums(SILVA_MERDS_rar_trt_live_trans_b_a) > 0, SILVA_MERDS_rar_trt_live_trans_b_a)
ntaxa(SILVA_MERDS_rar_trt_live_trans_b_a)
#3108
tran_bulk.amb.n<-taxa_names(SILVA_MERDS_rar_trt_live_trans_b_a)
length(tran_bulk.amb.n)
#3108

SILVA_MERDS_rar_trt_live_trans_r_a=subset_samples(SILVA_MERDS_rar_trt_live_trans, root_association=="R"&precip=="A")
SILVA_MERDS_rar_trt_live_trans_r_a=prune_taxa(taxa_sums(SILVA_MERDS_rar_trt_live_trans_r_a) > 0, SILVA_MERDS_rar_trt_live_trans_r_a)
ntaxa(SILVA_MERDS_rar_trt_live_trans_r_a)
#3463
trans_rhizo.amb.n<-taxa_names(SILVA_MERDS_rar_trt_live_trans_r_a)
length(trans_rhizo.amb.n)
#3463

trans.amb_b_in_r<-tran_bulk.amb.n[tran_bulk.amb.n %in% trans_rhizo.amb.n]
length(trans.amb_b_in_r) ## How many OTUs
#2032


SILVA_MERDS_rar_trt_live_trans_b_d=subset_samples(SILVA_MERDS_rar_trt_live_trans, root_association=="B"&precip=="D")
SILVA_MERDS_rar_trt_live_trans_b_d=prune_taxa(taxa_sums(SILVA_MERDS_rar_trt_live_trans_b_d) > 0, SILVA_MERDS_rar_trt_live_trans_b_d)
ntaxa(SILVA_MERDS_rar_trt_live_trans_b_d)
#2780
trans_bulk.drought.n<-taxa_names(SILVA_MERDS_rar_trt_live_trans_b_d)
length(trans_bulk.drought.n)
#2780

SILVA_MERDS_rar_trt_live_trans_r_d=subset_samples(SILVA_MERDS_rar_trt_live_trans, root_association=="R"&precip=="D")
SILVA_MERDS_rar_trt_live_trans_r_d=prune_taxa(taxa_sums(SILVA_MERDS_rar_trt_live_trans_r_d) > 0, SILVA_MERDS_rar_trt_live_trans_r_d)
ntaxa(SILVA_MERDS_rar_trt_live_trans_r_d)
#2405
trans_rhizo.drought.n<-taxa_names(SILVA_MERDS_rar_trt_live_trans_r_d)
length(trans_rhizo.drought.n)
#2405


trans.drought_b_in_r<-trans_bulk.drought.n[trans_bulk.drought.n %in% trans_rhizo.drought.n]
length(trans.drought_b_in_r) ## How many OTUs
#1593

trans.drought_in_amb<-trans.amb_b_in_r[trans.amb_b_in_r %in% trans.drought_b_in_r]
length(trans.drought_in_amb) ## How many OTUs
#1163


SILVA_MERDS_rar_trt_live_trans_sh<-prune_taxa(trans.drought_in_amb,SILVA_MERDS_rar_trt_live_trans)
SILVA_MERDS_rar_trt_live_trans_sh_sampl_sum=sample_sums(SILVA_MERDS_rar_trt_live_trans_sh)
ntaxa(SILVA_MERDS_rar_trt_live_trans_sh)
#1163
SILVA_MERDS_rar_trt_live_trans_map=sample_data(SILVA_MERDS_rar_trt_live_trans)
SILVA_MERDS_rar_trt_live_trans_map$soil_root=with(SILVA_MERDS_rar_trt_live_trans_map, interaction(soil_status,root_association))

sh_SILVA_MERDS_rar_trt_trans_sampl_sum=merge(SILVA_MERDS_rar_trt_live_trans_map,SILVA_MERDS_rar_trt_live_trans_sh_sampl_sum, by="row.names")


sh_SILVA_MERDS_rar_trt_trans_sampl_sum_precip_soil_g=sh_SILVA_MERDS_rar_trt_trans_sampl_sum %>% group_by(soil_root,precip)
sh_read_abund_precip_trans=summarise_at(sh_SILVA_MERDS_rar_trt_trans_sampl_sum_precip_soil_g, 
                                       "y", list(~n(),~mean,~sd,se=~sd(.)/sqrt(n())))

ggplot(sh_read_abund_precip_trans, aes(x=precip,y=mean,ymin = mean-se, ymax= mean+se,fill=soil_root))+
  geom_bar(position="dodge",stat="identity", color="black",aes(fill=soil_root))+
  geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+
  scale_fill_manual(values = c("light gray", "dark grey"),labels=c("Bulk","Rhizo"))+
  scale_x_discrete(labels=c("Ambient","Drought"))+ylab("Shared Taxa Read Abun")+
  geom_text(aes(y=mean+se+500, label=n),position=position_dodge(width=0.9))+theme_bw()


#Diversity of the Shared taxa
alpha_meas = c("Observed", "Shannon", "InvSimpson")
SILVA_MERDS_rar_trt_live_trans_map=sample_data(SILVA_MERDS_rar_trt_live_trans)
SILVA_MERDS_rar_trt_live_trans_map$soil_root=with(SILVA_MERDS_rar_trt_live_trans_map, interaction(soil_status,root_association))

SILVA_MERDS_rar_trt_live_trans_sh.divfil=estimate_richness(SILVA_MERDS_rar_trt_live_trans_sh,measures=alpha_meas)

SILVA_MERDS_rar_trt_live_trans_sh.divfil=merge(SILVA_MERDS_rar_trt_live_trans_sh.divfil, SILVA_MERDS_rar_trt_live_trans_map, by ="row.names")
#bact.soilE.t.divfil=mutate(bact.soilE.t.divfil, pielou=Shannon*(1/log(Observed)))
head(SILVA_MERDS_rar_trt_live_trans_sh.divfil)
row.names(SILVA_MERDS_rar_trt_live_trans_sh.divfil)=SILVA_MERDS_rar_trt_live_trans_sh.divfil$Row.names
SILVA_MERDS_rar_trt_live_trans_sh.divfil$Row.names=NULL



SILVA_MERDS_rar_trt_live_trans_sh.divfil_g=SILVA_MERDS_rar_trt_live_trans_sh.divfil %>% group_by(soil_root,precip)
sh_inv_simp_precip_trans=summarise_at(SILVA_MERDS_rar_trt_live_trans_sh.divfil_g, 
                                     "InvSimpson", list(~n(),~mean,~sd,se=~sd(.)/sqrt(n())))

ggplot(sh_inv_simp_precip_trans, aes(x=precip,y=mean,ymin = mean-se, ymax= mean+se,fill=soil_root))+
  geom_bar(position="dodge",stat="identity", color="black",aes(fill=soil_root))+
  geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+
  scale_fill_manual(values = c("light gray", "dark grey"),labels=c("Bulk","Rhizo"))+
  scale_x_discrete(labels=c("Ambient","Drought"))+ylab("Shared inv Simpson")+
  geom_text(aes(y=mean+se+5, label=n),position=position_dodge(width=0.9))+theme_bw()



######Trans DOMINANT TAXA Sterile and live####

SILVA_MERDS_rar_trans=subset_samples(SILVA_MERDS_rar, life_stage=="G")
summary(sample_data(SILVA_MERDS_rar_trans))
SILVA_MERDS_rar_trans_shared=core(SILVA_MERDS_rar_trans, detection = 0,prevalence = 0.9999999999)

SILVA_MERDS_rar_trans_shared_sampl_sum=sample_sums(SILVA_MERDS_rar_trans_shared)
ntaxa(SILVA_MERDS_rar_trans_shared)
#10
SILVA_MERDS_rar_trans_shared_map=sample_data(SILVA_MERDS_rar_trans_shared)
SILVA_MERDS_rar_trans_shared_map$soil_root=with(SILVA_MERDS_rar_trans_shared_map, interaction(soil_status,root_association))

SILVA_MERDS_rar_trans_shared_sampl_sum=merge(SILVA_MERDS_rar_trans_shared_map,SILVA_MERDS_rar_trans_shared_sampl_sum, by="row.names")




treatment_order=c("S.B","L.B","L.R")
SILVA_MERDS_rar_trans_shared_sampl_sum_precip_soil_g=SILVA_MERDS_rar_trans_shared_sampl_sum %>% group_by(soil_root,precip)
shared_read_abund_precip_trans_W_stl=summarise_at(SILVA_MERDS_rar_trans_shared_sampl_sum_precip_soil_g, 
                                                  "y", list(~n(),~mean(.),~sd(.),se=~sd(.)/sqrt(n())))

(core_trans_taxa=ggplot(shared_read_abund_precip_trans_W_stl, aes(x=precip,y=mean, ymin = mean-se, ymax= mean+se,fill=factor(soil_root,levels = treatment_order)))+
  geom_bar(position="dodge",stat="identity", color="black",aes(fill=factor(soil_root,levels = treatment_order)))+
  geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+
  scale_fill_manual(values = c( "white","light gray", "dark grey"),labels=c("Sterile","Bulk","Rhizo"), limits=treatment_order, name=NULL)+
  scale_x_discrete(labels=c("Ambient","Drought"), name=NULL)+ylab("Core taxa read # (in all samples)")+
  geom_text(aes(y=500, label=n),size=10,position=position_dodge(width=0.9))+theme_bw()+theme(axis.title = element_text(size = 23), axis.text = element_text(size = 23),
                                                                                      legend.position = c(0.85,.9), legend.text=element_text(size=20),
                                                                                      legend.background = element_rect(size=0.5,linetype="solid",colour ="black")))


#Diversity of the Shared taxa
alpha_meas = c("Observed", "Shannon", "InvSimpson")
SILVA_MERDS_rar_trans_shared_map=sample_data(SILVA_MERDS_rar_trans_shared)
SILVA_MERDS_rar_trans_shared_map$soil_root=with(SILVA_MERDS_rar_trans_shared_map, interaction(soil_status,root_association))

SILVA_MERDS_rar_trans_shared.divfil=estimate_richness(SILVA_MERDS_rar_trans_shared,measures=alpha_meas)

SILVA_MERDS_rar_trans_shared.divfil=merge(SILVA_MERDS_rar_trans_shared.divfil, SILVA_MERDS_rar_trans_shared_map, by ="row.names")
#bact.soilE.t.divfil=mutate(bact.soilE.t.divfil, pielou=Shannon*(1/log(Observed)))
head(SILVA_MERDS_rar_trans_shared.divfil)
row.names(SILVA_MERDS_rar_trans_shared.divfil)=SILVA_MERDS_rar_trans_shared.divfil$Row.names
SILVA_MERDS_rar_trans_shared.divfil$Row.names=NULL



SILVA_MERDS_rar_trans_shared.divfil_g=SILVA_MERDS_rar_trans_shared.divfil %>% group_by(soil_root,precip)
shared_inv_simp_precip_trans_w_stl=summarise_at(SILVA_MERDS_rar_trans_shared.divfil_g, 
                                          "InvSimpson", list(~n(),~mean(.),~sd(.),se=~sd(.)/sqrt(n())))

ggplot(shared_inv_simp_precip_trans_w_stl, aes(x=precip,y=mean,ymin = mean-se, ymax= mean+se,fill=factor(soil_root,levels = treatment_order)))+
  geom_bar(position="dodge",stat="identity", color="black",aes(fill=factor(soil_root,levels = treatment_order)))+
  geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+
  scale_fill_manual(values = c( "white","light gray", "dark grey"),labels=c("Sterile","Bulk","Rhizo"), limits=treatment_order, name=NULL)+
  scale_x_discrete(labels=c("Ambient","Drought"), name=NULL)+ylab("Core diversity invSimpson \n(in all samples)")+
  geom_text(aes(y=1, label=n),size=10,position=position_dodge(width=0.9))+theme_bw()+
  theme(axis.title = element_text(size = 23), axis.text = element_text(size = 23),
        legend.position = c(0.85,.9), legend.text=element_text(size=20),
        legend.background = element_rect(size=0.5,linetype="solid",colour ="black"))+
  ylim(c(0,8.5))


#Found in at least one sample per treatment
SILVA_MERDS_rar_trt_live_trans_map=sample_data(SILVA_MERDS_rar_trt_live_trans)
SILVA_MERDS_rar_trt_live_trans_map$soil_root=with(SILVA_MERDS_rar_trt_live_trans_map, interaction(soil_status,root_association))
summary(SILVA_MERDS_rar_trt_live_trans_map)

SILVA_MERDS_rar_trt_live_trans_b_a=subset_samples(SILVA_MERDS_rar_trt_live_trans, root_association=="B"&precip=="A")
SILVA_MERDS_rar_trt_live_trans_b_a=prune_taxa(taxa_sums(SILVA_MERDS_rar_trt_live_trans_b_a) > 0, SILVA_MERDS_rar_trt_live_trans_b_a)
ntaxa(SILVA_MERDS_rar_trt_live_trans_b_a)
#3108
tran_bulk.amb.n<-taxa_names(SILVA_MERDS_rar_trt_live_trans_b_a)
length(tran_bulk.amb.n)
#3108

SILVA_MERDS_rar_trt_live_trans_r_a=subset_samples(SILVA_MERDS_rar_trt_live_trans, root_association=="R"&precip=="A")
SILVA_MERDS_rar_trt_live_trans_r_a=prune_taxa(taxa_sums(SILVA_MERDS_rar_trt_live_trans_r_a) > 0, SILVA_MERDS_rar_trt_live_trans_r_a)
ntaxa(SILVA_MERDS_rar_trt_live_trans_r_a)
#3463
trans_rhizo.amb.n<-taxa_names(SILVA_MERDS_rar_trt_live_trans_r_a)
length(trans_rhizo.amb.n)
#3463

trans.amb_b_in_r<-tran_bulk.amb.n[tran_bulk.amb.n %in% trans_rhizo.amb.n]
length(trans.amb_b_in_r) ## How many OTUs
#2032


SILVA_MERDS_rar_trt_live_trans_b_d=subset_samples(SILVA_MERDS_rar_trt_live_trans, root_association=="B"&precip=="D")
SILVA_MERDS_rar_trt_live_trans_b_d=prune_taxa(taxa_sums(SILVA_MERDS_rar_trt_live_trans_b_d) > 0, SILVA_MERDS_rar_trt_live_trans_b_d)
ntaxa(SILVA_MERDS_rar_trt_live_trans_b_d)
#2780
trans_bulk.drought.n<-taxa_names(SILVA_MERDS_rar_trt_live_trans_b_d)
length(trans_bulk.drought.n)
#2780

SILVA_MERDS_rar_trt_live_trans_r_d=subset_samples(SILVA_MERDS_rar_trt_live_trans, root_association=="R"&precip=="D")
SILVA_MERDS_rar_trt_live_trans_r_d=prune_taxa(taxa_sums(SILVA_MERDS_rar_trt_live_trans_r_d) > 0, SILVA_MERDS_rar_trt_live_trans_r_d)
ntaxa(SILVA_MERDS_rar_trt_live_trans_r_d)
#2405
trans_rhizo.drought.n<-taxa_names(SILVA_MERDS_rar_trt_live_trans_r_d)
length(trans_rhizo.drought.n)
#2405


trans.drought_b_in_r<-trans_bulk.drought.n[trans_bulk.drought.n %in% trans_rhizo.drought.n]
length(trans.drought_b_in_r) ## How many OTUs
#1593

trans.drought_in_amb<-trans.amb_b_in_r[trans.amb_b_in_r %in% trans.drought_b_in_r]
length(trans.drought_in_amb) ## How many OTUs
#1163


SILVA_MERDS_rar_trt_live_trans_sh<-prune_taxa(trans.drought_in_amb,SILVA_MERDS_rar_trt_live_trans)
SILVA_MERDS_rar_trt_live_trans_sh_sampl_sum=sample_sums(SILVA_MERDS_rar_trt_live_trans_sh)
ntaxa(SILVA_MERDS_rar_trt_live_trans_sh)
#1163
SILVA_MERDS_rar_trt_live_trans_map=sample_data(SILVA_MERDS_rar_trt_live_trans)
SILVA_MERDS_rar_trt_live_trans_map$soil_root=with(SILVA_MERDS_rar_trt_live_trans_map, interaction(soil_status,root_association))

sh_SILVA_MERDS_rar_trt_trans_sampl_sum=merge(SILVA_MERDS_rar_trt_live_trans_map,SILVA_MERDS_rar_trt_live_trans_sh_sampl_sum, by="row.names")


sh_SILVA_MERDS_rar_trt_trans_sampl_sum_precip_soil_g=sh_SILVA_MERDS_rar_trt_trans_sampl_sum %>% group_by(soil_root,precip)
sh_read_abund_precip_trans=summarise_at(sh_SILVA_MERDS_rar_trt_trans_sampl_sum_precip_soil_g, 
                                        "y", list(~n(),~mean,~sd,se=~sd(.)/sqrt(n())))

ggplot(sh_read_abund_precip_trans, aes(x=precip,y=mean,ymin = mean-se, ymax= mean+se,fill=soil_root))+
  geom_bar(position="dodge",stat="identity", color="black",aes(fill=soil_root))+
  geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+
  scale_fill_manual(values = c("light gray", "dark grey"),labels=c("Bulk","Rhizo"))+
  scale_x_discrete(labels=c("Ambient","Drought"))+ylab("Shared Taxa Read Abun")+
  geom_text(aes(y=mean+se+500, label=n),position=position_dodge(width=0.9))+theme_bw()


#Diversity of the Shared taxa
alpha_meas = c("Observed", "Shannon", "InvSimpson")
SILVA_MERDS_rar_trt_live_trans_map=sample_data(SILVA_MERDS_rar_trt_live_trans)
SILVA_MERDS_rar_trt_live_trans_map$soil_root=with(SILVA_MERDS_rar_trt_live_trans_map, interaction(soil_status,root_association))

SILVA_MERDS_rar_trt_live_trans_sh.divfil=estimate_richness(SILVA_MERDS_rar_trt_live_trans_sh,measures=alpha_meas)

SILVA_MERDS_rar_trt_live_trans_sh.divfil=merge(SILVA_MERDS_rar_trt_live_trans_sh.divfil, SILVA_MERDS_rar_trt_live_trans_map, by ="row.names")
#bact.soilE.t.divfil=mutate(bact.soilE.t.divfil, pielou=Shannon*(1/log(Observed)))
head(SILVA_MERDS_rar_trt_live_trans_sh.divfil)
row.names(SILVA_MERDS_rar_trt_live_trans_sh.divfil)=SILVA_MERDS_rar_trt_live_trans_sh.divfil$Row.names
SILVA_MERDS_rar_trt_live_trans_sh.divfil$Row.names=NULL



SILVA_MERDS_rar_trt_live_trans_sh.divfil_g=SILVA_MERDS_rar_trt_live_trans_sh.divfil %>% group_by(soil_root,precip)
sh_inv_simp_precip_trans=summarise_at(SILVA_MERDS_rar_trt_live_trans_sh.divfil_g, 
                                      "InvSimpson", list(~n(),~mean,~sd,se=~sd(.)/sqrt(n())))

ggplot(sh_inv_simp_precip_trans, aes(x=precip,y=mean,ymin = mean-se, ymax= mean+se,fill=soil_root))+
  geom_bar(position="dodge",stat="identity", color="black",aes(fill=soil_root))+
  geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+
  scale_fill_manual(values = c("light gray", "dark grey"),labels=c("Bulk","Rhizo"))+
  scale_x_discrete(labels=c("Ambient","Drought"))+ylab("Shared inv Simpson")+
  geom_text(aes(y=mean+se+5, label=n),position=position_dodge(width=0.9))+theme_bw()




#####Seeds####

SILVA_MERDS_rar_trt_seed=subset_samples(SILVA_MERDS_rar, life_stage=="S"&life_stage!="Start")
nrow(sample_data(SILVA_MERDS_rar_trt_seed))
#23

#####Bray Live and Sterile Trans Comm Analyses####

SILVA_MERDS_rar_trt_seed_ord=ordinate(SILVA_MERDS_rar_trt_seed, method = "NMDS",distance = "bray")
#*** Solution reached
#Warning message:
#In metaMDS(veganifyOTU(physeq), distance, ...) :
#stress is (nearly) zero: you may have insufficient data
#8.513795e-05
plot_ordination(SILVA_MERDS_rar_trt_seed,SILVA_MERDS_rar_trt_seed_ord, color="precip",shape="soil_root")+geom_point(size=4)+
  geom_label_repel(size=3,aes(label = block))+theme_bw()


SILVA_MERDS_rar_trt_seed_map=sample_data(SILVA_MERDS_rar_trt_seed)
SILVA_MERDS_rar_trt_seed_map$soil_root=with(SILVA_MERDS_rar_trt_seed_map, interaction(soil_status,root_association))
SILVA_MERDS_rar_trt_seed_dis=distance(SILVA_MERDS_rar_trt_seed,method = "bray")

adonis(SILVA_MERDS_rar_trt_seed_dis~SILVA_MERDS_rar_trt_seed_map$soil_root*SILVA_MERDS_rar_trt_seed_map$precip
       +as.factor(SILVA_MERDS_rar_trt_seed_map$block), permutations = 9999)
#SILVA_MERDS_rar_trt_seed_map$soil_root                                      2    2.4268 1.21338  9.2823 0.43181 0.0001 ***
#SILVA_MERDS_rar_trt_seed_map$precip                                         1    0.4933 0.49326  3.7734 0.08777 0.0046 **  
#SILVA_MERDS_rar_trt_seed_map$soil_root:SILVA_MERDS_rar_trt_seed_map$precip  2    0.4361 0.21807  1.6682 0.07761 0.0792 .  


#####Jaccard Live and Sterile Trans Comm Analyses####

SILVA_MERDS_rar_trt_seed_J_ord=ordinate(SILVA_MERDS_rar_trt_seed, method = "NMDS", distance = "jaccard", binary = TRUE)
#Warning message:
#In metaMDS(veganifyOTU(physeq), distance, ...) :
#  stress is (nearly) zero: you may have insufficient data
#7.751021e-05
plot_ordination(SILVA_MERDS_rar_trt_seed,SILVA_MERDS_rar_trt_seed_J_ord, color="precip",shape="soil_root")+geom_point(size=4)+
  geom_label_repel(size=3,aes(label = block))+theme_bw()


SILVA_MERDS_rar_trt_seed_map=sample_data(SILVA_MERDS_rar_trt_seed)
SILVA_MERDS_rar_trt_seed_map$soil_root=with(SILVA_MERDS_rar_trt_seed_map, interaction(soil_status,root_association))
SILVA_MERDS_rar_trt_seed_J_dis=distance(SILVA_MERDS_rar_trt_seed,method = "jaccard", binary = TRUE)

adonis(SILVA_MERDS_rar_trt_seed_J_dis~SILVA_MERDS_rar_trt_seed_map$soil_root*SILVA_MERDS_rar_trt_seed_map$precip
       +as.factor(SILVA_MERDS_rar_trt_seed_map$block), permutations = 9999)
#SILVA_MERDS_rar_trt_seed_map$soil_root                                      2    2.3676 1.18381  4.8965 0.31628 0.0001 ***
#SILVA_MERDS_rar_trt_seed_map$precip                                         1    0.3449 0.34488  1.4265 0.04607 0.0913 .  

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










#####UnWeighted Unifrac Live and Sterile seed Comm Analyses####

SILVA_MERDS_rar_trt_seed_unWU_ord=ordinate(SILVA_MERDS_rar_trt_seed, method = "NMDS",distance = "unifrac")
#*** Solution reached
#Warning message:
#In metaMDS(ps.dist) :
#  stress is (nearly) zero: you may have insufficient data
#7.563441e-05
plot_ordination(SILVA_MERDS_rar_trt_seed,SILVA_MERDS_rar_trt_seed_unWU_ord, color="precip",shape="soil_root")+geom_point(size=4)+
  geom_label_repel(size=3,aes(label = block))+theme_bw()


SILVA_MERDS_rar_trt_seed_map=sample_data(SILVA_MERDS_rar_trt_seed)
SILVA_MERDS_rar_trt_seed_map$soil_root=with(SILVA_MERDS_rar_trt_seed_map, interaction(soil_status,root_association))
SILVA_MERDS_rar_trt_seed_unWU_dis=distance(SILVA_MERDS_rar_trt_seed,method = "unifrac")

adonis(SILVA_MERDS_rar_trt_seed_unWU_dis~SILVA_MERDS_rar_trt_seed_map$soil_root*SILVA_MERDS_rar_trt_seed_map$precip
       +as.factor(SILVA_MERDS_rar_trt_seed_map$block), permutations = 9999)
#SILVA_MERDS_rar_trt_seed_map$soil_root                                      2    2.4126 1.20629  8.1897 0.42986 0.0001 ***
#SILVA_MERDS_rar_trt_seed_map$precip                                         1    0.2466 0.24663  1.6744 0.04394 0.0885 .  


#####Presence UnWeighted Unifrac Seed####
head(sample_data(SILVA_MERDS_rar))

SILVA_MERDS_rar_pres_seed=subset_samples(SILVA_MERDS_rar, life_stage=="S"&life_stage!="Start"&
                                           root_association=="B")
nrow(sample_data(SILVA_MERDS_rar_pres_seed))
#15
unique(sample_data(SILVA_MERDS_rar_pres_seed)$soil_root)

SILVA_MERDS_rar_pres_seed_map=sample_data(SILVA_MERDS_rar_pres_seed)
SILVA_MERDS_rar_pres_seed_map$soil_root=with(SILVA_MERDS_rar_pres_seed_map, interaction(soil_status,root_association))
SILVA_MERDS_rar_pres_seed_UnWU_dis=distance(SILVA_MERDS_rar_pres_seed,method = "unifrac")

(pres_seed_UnWU=adonis(SILVA_MERDS_rar_pres_seed_UnWU_dis~SILVA_MERDS_rar_pres_seed_map$soil_root*SILVA_MERDS_rar_pres_seed_map$precip
                     +as.factor(SILVA_MERDS_rar_pres_seed_map$block), permutations = 9999))


AICc.PERMANOVA(pres_seed_UnWU)
#$AIC
#[1] 16.00398


pairwise.perm.manova(SILVA_MERDS_rar_pres_seed_UnWU_dis, SILVA_MERDS_rar_map$soil_root_stage, nperm=2000)


(pres_seed_UnWU_no_block=adonis(SILVA_MERDS_rar_pres_seed_UnWU_dis~SILVA_MERDS_rar_pres_seed_map$soil_root*
                                SILVA_MERDS_rar_pres_seed_map$precip, permutations = 9999))


AICc.PERMANOVA(pres_seed_UnWU_no_block)
#$AIC
#[1] 15.3049


#####Origin UnWeighted Unifrac Seed####
head(sample_data(SILVA_MERDS_rar))

SILVA_MERDS_rar_orig_seed=subset_samples(SILVA_MERDS_rar, life_stage=="S"&life_stage!="Start"&
                                           soil_status =="L")
nrow(sample_data(SILVA_MERDS_rar_orig_seed))
#16
unique(sample_data(SILVA_MERDS_rar_orig_seed)$soil_root)

SILVA_MERDS_rar_orig_seed_map=sample_data(SILVA_MERDS_rar_orig_seed)
SILVA_MERDS_rar_orig_seed_map$soil_root=with(SILVA_MERDS_rar_orig_seed_map, interaction(soil_status,root_association))
SILVA_MERDS_rar_orig_seed_UnWU_dis=distance(SILVA_MERDS_rar_orig_seed,method = "unifrac")

(orig_seed_UnWU=adonis(SILVA_MERDS_rar_orig_seed_UnWU_dis~SILVA_MERDS_rar_orig_seed_map$soil_root*SILVA_MERDS_rar_orig_seed_map$precip
                     +as.factor(SILVA_MERDS_rar_orig_seed_map$block), permutations = 9999))


AICc.PERMANOVA(orig_seed_UnWU)
#$AIC
#[1] 18.48733


pairwise.perm.manova(SILVA_MERDS_rar_orig_seed_UnWU_dis, SILVA_MERDS_rar_map$soil_root_stage, nperm=2000)


(orig_seed_UnWU_no_block=adonis(SILVA_MERDS_rar_orig_seed_UnWU_dis~SILVA_MERDS_rar_orig_seed_map$soil_root*
                                SILVA_MERDS_rar_orig_seed_map$precip, permutations = 9999))


AICc.PERMANOVA(orig_seed_UnWU_no_block)
#$AIC
#[1] 17.70276






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



#####Seed live Community analyses####

#####Bray Live Seed Comm Analyses####

#No sterile
SILVA_MERDS_rar_trt_live_seed=subset_samples(SILVA_MERDS_rar_trt_seed, soil_root!="S.B")

SILVA_MERDS_rar_trt_live_seed_ord=ordinate(SILVA_MERDS_rar_trt_live_seed, method = "NMDS",distance = "bray")
#*** No convergence -- monoMDS stopping criteria:
#20: stress ratio > sratmax
#0.1359863 
plot_ordination(SILVA_MERDS_rar_trt_live_seed,SILVA_MERDS_rar_trt_live_seed_ord, color="precip",shape="root_association")+geom_point(size=4)+
  geom_label_repel(size=3,aes(label = block))+theme_bw()


SILVA_MERDS_rar_trt_live_seed_map=sample_data(SILVA_MERDS_rar_trt_live_seed)
SILVA_MERDS_rar_trt_live_seed_map$soil_root=with(SILVA_MERDS_rar_trt_live_seed_map, interaction(soil_status,root_association))
SILVA_MERDS_rar_trt_live_seed_dis=distance(SILVA_MERDS_rar_trt_live_seed,method = "bray")

adonis(SILVA_MERDS_rar_trt_live_seed_dis~SILVA_MERDS_rar_trt_live_seed_map$soil_root*SILVA_MERDS_rar_trt_live_seed_map$precip
       +as.factor(SILVA_MERDS_rar_trt_live_seed_map$block), permutations = 9999)
#SILVA_MERDS_rar_trt_live_seed_map$soil_root                                           1   0.38611 0.38611  3.3937 0.16782 0.0001 ***
#SILVA_MERDS_rar_trt_live_seed_map$precip                                              1   0.33423 0.33423  2.9378 0.14527 0.0004 ***


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
#####Jaccard Live Seed Comm Analyses####

SILVA_MERDS_rar_trt_live_seed_J_ord=ordinate(SILVA_MERDS_rar_trt_live_seed, method = "NMDS", distance = "jaccard", binary = TRUE)
#*** Solution reached
#0.1515835
plot_ordination(SILVA_MERDS_rar_trt_live_seed,SILVA_MERDS_rar_trt_live_seed_J_ord, color="precip",shape="root_association")+geom_point(size=4)+
  geom_label_repel(size=3,aes(label = block))+theme_bw()


SILVA_MERDS_rar_trt_live_seed_map=sample_data(SILVA_MERDS_rar_trt_live_seed)
SILVA_MERDS_rar_trt_live_seed_map$soil_root=with(SILVA_MERDS_rar_trt_live_seed_map, interaction(soil_status,root_association))
SILVA_MERDS_rar_trt_live_seed_J__dis=distance(SILVA_MERDS_rar_trt_live_seed,method = "jaccard", binary = TRUE)

adonis(SILVA_MERDS_rar_trt_live_seed_J__dis~SILVA_MERDS_rar_trt_live_seed_map$soil_root*SILVA_MERDS_rar_trt_live_seed_map$precip
       +as.factor(SILVA_MERDS_rar_trt_live_seed_map$block), permutations = 9999)
#SILVA_MERDS_rar_trt_live_seed_map$soil_root                                           1    0.3042 0.30420 1.28868 0.07996 0.0040 **
#SILVA_MERDS_rar_trt_live_seed_map$precip                                              1    0.3102 0.31015 1.31390 0.08153 0.0032 **
#as.factor(SILVA_MERDS_rar_trt_live_seed_map$block)                                    3    0.8328 0.27759 1.17596 0.21890 0.0032 **


#####Weighted Unifrac Live seed Comm Analyses####

SILVA_MERDS_rar_trt_live_seed_WU_ord=ordinate(SILVA_MERDS_rar_trt_live_seed, method = "NMDS", distance = "wunifrac")
#*** Solution reached
#0.05507781
plot_ordination(SILVA_MERDS_rar_trt_live_seed,SILVA_MERDS_rar_trt_live_seed_WU_ord, color="precip",shape="root_association")+geom_point(size=4)+
  geom_label_repel(size=3,aes(label = block))+theme_bw()


SILVA_MERDS_rar_trt_live_seed_map=sample_data(SILVA_MERDS_rar_trt_live_seed)
SILVA_MERDS_rar_trt_live_seed_map$soil_root=with(SILVA_MERDS_rar_trt_live_seed_map, interaction(soil_status,root_association))
SILVA_MERDS_rar_trt_live_seed_WU_dis=distance(SILVA_MERDS_rar_trt_live_seed,method = "wunifrac")

adonis(SILVA_MERDS_rar_trt_live_seed_WU_dis~SILVA_MERDS_rar_trt_live_seed_map$soil_root*SILVA_MERDS_rar_trt_live_seed_map$precip
       +as.factor(SILVA_MERDS_rar_trt_live_seed_map$block), permutations = 9999)
#SILVA_MERDS_rar_trt_live_seed_map$soil_root                                           1  0.079780 0.079780  7.4430 0.31261 0.0002 ***
#SILVA_MERDS_rar_trt_live_seed_map$precip                                              1  0.033106 0.033106  3.0886 0.12972 0.0286 *  


#####UnWeighted Unifrac Live seed Comm Analyses####

SILVA_MERDS_rar_trt_live_seed_unWU_ord=ordinate(SILVA_MERDS_rar_trt_live_seed, method = "NMDS", distance = "unifrac")
#*** Solution reached
#0.1689273
plot_ordination(SILVA_MERDS_rar_trt_live_seed,SILVA_MERDS_rar_trt_live_seed_unWU_ord, color="precip",shape="root_association")+geom_point(size=4)+
  geom_label_repel(size=3,aes(label = block))+theme_bw()


SILVA_MERDS_rar_trt_live_seed_map=sample_data(SILVA_MERDS_rar_trt_live_seed)
SILVA_MERDS_rar_trt_live_seed_map$soil_root=with(SILVA_MERDS_rar_trt_live_seed_map, interaction(soil_status,root_association))
SILVA_MERDS_rar_trt_live_seed_unWU_dis=distance(SILVA_MERDS_rar_trt_live_seed,method = "unifrac")

adonis(SILVA_MERDS_rar_trt_live_seed_unWU_dis~SILVA_MERDS_rar_trt_live_seed_map$soil_root*SILVA_MERDS_rar_trt_live_seed_map$precip
       +as.factor(SILVA_MERDS_rar_trt_live_seed_map$block), permutations = 9999)
#SILVA_MERDS_rar_trt_live_seed_map$soil_root                                           1   0.20284 0.20284 1.37881 0.08501 0.0054 ** 
#SILVA_MERDS_rar_trt_live_seed_map$precip                                              1   0.21860 0.21860 1.48592 0.09161 0.0002 ***
#as.factor(SILVA_MERDS_rar_trt_live_seed_map$block)                                    3   0.51058 0.17019 1.15688 0.21398 0.0219 *  


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


(seed_obs_richness_p2=ggplot(obs_rich_precip_seed, aes(x=precip,y=mean, ymin = mean-se, ymax= mean+se,fill=factor(soil_root,levels = treatment_order)))+
    geom_bar(position="dodge",stat="identity", color="black",aes(fill=factor(soil_root,levels = treatment_order)))+
    geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+ylim(c(0,1850))+
    scale_fill_manual(values = c( "white","lightgray", "darkgrey"),labels=c("Sterile","Bulk","Rhizo"), limits=treatment_order, name=NULL)+
    scale_x_discrete(labels=c("Ambient","Drought"), name=NULL)+ylab("Bacterial richness")+
    geom_text(aes(y=80, label=n),size=10,position=position_dodge(width=0.9))+theme_bw()+
    theme(axis.title.x = element_text(size = 23), axis.text.x = element_text(size = 23),
          axis.title.y = element_text(size = 23), axis.text.y = element_text(size = 23),
          legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank()))
#bact_richness_germination


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

(seed_inv_simpson_p2=ggplot(inv_simp_precip_seed, aes(x=precip,y=mean, ymin = mean-se, ymax= mean+se,fill=factor(soil_root,levels = treatment_order)))+
    geom_bar(position="dodge",stat="identity", color="black",aes(fill=factor(soil_root,levels = treatment_order)))+
    geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+ylim(c(-10,200))+
    scale_fill_manual(values = c( "white","lightgray", "darkgrey"),labels=c("Sterile","Bulk","Rhizo"), limits=treatment_order, name=NULL)+
    scale_x_discrete(labels=c("Ambient","Drought"), name=NULL)+ylab("inv Simpson bacteria")+
    geom_text(aes(y=-7, label=n),size=10,position=position_dodge(width=0.9))+theme_bw()+
    theme(axis.title.x = element_text(size = 23), axis.text.x = element_text(size = 23),
          axis.title.y = element_text(size = 23), axis.text.y = element_text(size = 23),
          legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank()))


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


#####Are shared OTUs dominant####
SILVA_MERDS_rar_trt_live_seed_shared=core(SILVA_MERDS_rar_trt_live_seed, detection = 0,prevalence = 0.9999999999)
SILVA_MERDS_rar_trt_live_seed_shared_sampl_sum=sample_sums(SILVA_MERDS_rar_trt_live_seed_shared)
ntaxa(SILVA_MERDS_rar_trt_live_seed_shared)
#107
SILVA_MERDS_rar_trt_live_seed_map=sample_data(SILVA_MERDS_rar_trt_live_seed)
SILVA_MERDS_rar_trt_live_seed_map$soil_root=with(SILVA_MERDS_rar_trt_live_seed_map, interaction(soil_status,root_association))

shared_SILVA_MERDS_rar_trt_seed_sampl_sum=merge(SILVA_MERDS_rar_trt_live_seed_map,SILVA_MERDS_rar_trt_live_seed_shared_sampl_sum, by="row.names")


shared_SILVA_MERDS_rar_trt_seed_sampl_sum_precip_soil_g=shared_SILVA_MERDS_rar_trt_seed_sampl_sum %>% group_by(soil_root,precip)
shared_read_abund_precip_seed=summarise_at(shared_SILVA_MERDS_rar_trt_seed_sampl_sum_precip_soil_g, 
                                  "y", list(~n(),~mean,~sd,se=~sd(.)/sqrt(n())))

ggplot(shared_read_abund_precip_seed, aes(x=precip,y=mean,ymin = mean-se, ymax= mean+se,fill=soil_root))+
  geom_bar(position="dodge",stat="identity", color="black",aes(fill=soil_root))+
  geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+
  scale_fill_manual(values = c("light gray", "dark grey"),labels=c("Bulk","Rhizo"))+
  scale_x_discrete(labels=c("Ambient","Drought"))+ylab("Shared Taxa Read Abun")+
  geom_text(aes(y=mean+se+100, label=n),position=position_dodge(width=0.9))+theme_bw()



source(system.file("extdata/lm_phyloseq.R", package = "microbiome"))
days_to_germ_mod_core_otus=lm_phyloseq(SILVA_MERDS_rar_trt_live_seed_shared, "days_to_germ")

"              logFC   AveExpr         t     P.Value adj.P.Val         B
OTU203    0.3816643 0.8371876  3.084710 0.005364603 0.5189346 -3.143445
OTU565    0.3102816 0.6502313  2.828936 0.009699713 0.5189346 -3.364949
OTU496   -0.1735765 0.6949704 -2.213914 0.037372806 0.8378279 -3.879686
OTU41     0.2347644 0.8932648  2.152420 0.042456020 0.8378279 -3.928634
OTU9      0.1848904 1.6610634  2.021643 0.055391871 0.8378279 -4.030560
OTU669   -0.2458820 1.0946091 -1.951879 0.063640635 0.8378279 -4.083599
OTU861    0.1908916 1.1128776  1.898917 0.070608557 0.8378279 -4.123189
OTU29761  0.2230778 0.6913682  1.863623 0.075615246 0.8378279 -4.149231
OTU403    0.1911195 0.7680598  1.777272 0.089186455 0.8378279 -4.211731
OTU523    0.1753177 1.1708752  1.733702 0.096797100 0.8378279 -4.242580"

par(mfrow=c(2,2))
plot(get_sample(SILVA_MERDS_rar_trt_live_seed_shared,"OTU203"),SILVA_MERDS_rar_trt_live_seed_map$days_to_germ)
plot(get_sample(SILVA_MERDS_rar_trt_live_seed_shared,"OTU565"),SILVA_MERDS_rar_trt_live_seed_map$days_to_germ)
plot(get_sample(SILVA_MERDS_rar_trt_live_seed_shared,"OTU496"),SILVA_MERDS_rar_trt_live_seed_map$days_to_germ)
plot(get_sample(SILVA_MERDS_rar_trt_live_seed_shared,"OTU41"),SILVA_MERDS_rar_trt_live_seed_map$days_to_germ)
par(mfrow=c(2,2))
plot(get_sample(SILVA_MERDS_rar_trt_live_seed_shared,"OTU9"),SILVA_MERDS_rar_trt_live_seed_map$days_to_germ)
plot(get_sample(SILVA_MERDS_rar_trt_live_seed_shared,"OTU669"),SILVA_MERDS_rar_trt_live_seed_map$days_to_germ)
plot(get_sample(SILVA_MERDS_rar_trt_live_seed_shared,"OTU861"),SILVA_MERDS_rar_trt_live_seed_map$days_to_germ)
plot(get_sample(SILVA_MERDS_rar_trt_live_seed_shared,"OTU29761"),SILVA_MERDS_rar_trt_live_seed_map$days_to_germ)
par(mfrow=c(2,1))
plot(get_sample(SILVA_MERDS_rar_trt_live_seed_shared,"OTU403"),SILVA_MERDS_rar_trt_live_seed_map$days_to_germ)
plot(get_sample(SILVA_MERDS_rar_trt_live_seed_shared,"OTU523"),SILVA_MERDS_rar_trt_live_seed_map$days_to_germ)



#Diversity of the Shared taxa
alpha_meas = c("Observed", "Shannon", "InvSimpson")
SILVA_MERDS_rar_trt_live_seed_map=sample_data(SILVA_MERDS_rar_trt_live_seed)
SILVA_MERDS_rar_trt_live_seed_map$soil_root=with(SILVA_MERDS_rar_trt_live_seed_map, interaction(soil_status,root_association))

SILVA_MERDS_rar_trt_live_seed_shared.divfil=estimate_richness(SILVA_MERDS_rar_trt_live_seed_shared,measures=alpha_meas)

SILVA_MERDS_rar_trt_live_seed_shared.divfil=merge(SILVA_MERDS_rar_trt_live_seed_shared.divfil, SILVA_MERDS_rar_trt_live_seed_map, by ="row.names")
#bact.soilE.t.divfil=mutate(bact.soilE.t.divfil, pielou=Shannon*(1/log(Observed)))
head(SILVA_MERDS_rar_trt_live_seed_shared.divfil)
row.names(SILVA_MERDS_rar_trt_live_seed_shared.divfil)=SILVA_MERDS_rar_trt_live_seed_shared.divfil$Row.names
SILVA_MERDS_rar_trt_live_seed_shared.divfil$Row.names=NULL



SILVA_MERDS_rar_trt_live_seed_shared.divfil_g=SILVA_MERDS_rar_trt_live_seed_shared.divfil %>% group_by(soil_root,precip)
shared_inv_simp_precip_seed=summarise_at(SILVA_MERDS_rar_trt_live_seed_shared.divfil_g, 
                                           "InvSimpson", list(~n(),~mean,~sd,se=~sd(.)/sqrt(n())))

ggplot(shared_inv_simp_precip_seed, aes(x=precip,y=mean,ymin = mean-se, ymax= mean+se,fill=soil_root))+
  geom_bar(position="dodge",stat="identity", color="black",aes(fill=soil_root))+
  geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+
  scale_fill_manual(values = c("light gray", "dark grey"),labels=c("Bulk","Rhizo"))+
  scale_x_discrete(labels=c("Ambient","Drought"))+ylab("Shared inv Simpson")+
  geom_text(aes(y=mean+se+1, label=n),position=position_dodge(width=0.9))+theme_bw()

#Found in at least one sample per treatment
SILVA_MERDS_rar_trt_live_seed_map=sample_data(SILVA_MERDS_rar_trt_live_seed)
SILVA_MERDS_rar_trt_live_seed_map$soil_root=with(SILVA_MERDS_rar_trt_live_seed_map, interaction(soil_status,root_association))
summary(SILVA_MERDS_rar_trt_live_seed_map)

SILVA_MERDS_rar_trt_live_seed_b_a=subset_samples(SILVA_MERDS_rar_trt_live_seed, root_association=="B"&precip=="A")
SILVA_MERDS_rar_trt_live_seed_b_a=prune_taxa(taxa_sums(SILVA_MERDS_rar_trt_live_seed_b_a) > 0, SILVA_MERDS_rar_trt_live_seed_b_a)
ntaxa(SILVA_MERDS_rar_trt_live_seed_b_a)
#2459
bulk.amb.n<-taxa_names(SILVA_MERDS_rar_trt_live_seed_b_a)
length(bulk.amb.n)
#2459

SILVA_MERDS_rar_trt_live_seed_r_a=subset_samples(SILVA_MERDS_rar_trt_live_seed, root_association=="R"&precip=="A")
SILVA_MERDS_rar_trt_live_seed_r_a=prune_taxa(taxa_sums(SILVA_MERDS_rar_trt_live_seed_r_a) > 0, SILVA_MERDS_rar_trt_live_seed_r_a)
ntaxa(SILVA_MERDS_rar_trt_live_seed_r_a)
#3143
rhizo.amb.n<-taxa_names(SILVA_MERDS_rar_trt_live_seed_r_a)
length(rhizo.amb.n)
#3143

bact.amb_b_in_r<-bulk.amb.n[bulk.amb.n %in% rhizo.amb.n]
length(bact.amb_b_in_r) ## How many OTUs
#1653


SILVA_MERDS_rar_trt_live_seed_b_d=subset_samples(SILVA_MERDS_rar_trt_live_seed, root_association=="B"&precip=="D")
SILVA_MERDS_rar_trt_live_seed_b_d=prune_taxa(taxa_sums(SILVA_MERDS_rar_trt_live_seed_b_d) > 0, SILVA_MERDS_rar_trt_live_seed_b_d)
ntaxa(SILVA_MERDS_rar_trt_live_seed_b_d)
#2294
bulk.drought.n<-taxa_names(SILVA_MERDS_rar_trt_live_seed_b_d)
length(bulk.drought.n)
#2294

SILVA_MERDS_rar_trt_live_seed_r_d=subset_samples(SILVA_MERDS_rar_trt_live_seed, root_association=="R"&precip=="D")
SILVA_MERDS_rar_trt_live_seed_r_d=prune_taxa(taxa_sums(SILVA_MERDS_rar_trt_live_seed_r_d) > 0, SILVA_MERDS_rar_trt_live_seed_r_d)
ntaxa(SILVA_MERDS_rar_trt_live_seed_r_d)
#2992
rhizo.drought.n<-taxa_names(SILVA_MERDS_rar_trt_live_seed_r_d)
length(rhizo.drought.n)
#2992


bact.drought_b_in_r<-bulk.drought.n[bulk.drought.n %in% rhizo.drought.n]
length(bact.drought_b_in_r) ## How many OTUs
#1560

seed.drought_in_amb<-bact.amb_b_in_r[bact.amb_b_in_r %in% bact.drought_b_in_r]
length(seed.drought_in_amb) ## How many OTUs
#1025


SILVA_MERDS_rar_trt_live_seed_sh<-prune_taxa(seed.drought_in_amb,SILVA_MERDS_rar_trt_live_seed)

source(system.file("extdata/lm_phyloseq.R", package = "microbiome"))
days_to_germ_mod_shared_otus=lm_phyloseq(SILVA_MERDS_rar_trt_live_seed_sh, "days_to_germ")

"              logFC   AveExpr         t     P.Value adj.P.Val         B
OTU3210  -0.4579204 0.5392178 -3.605444 0.001121532 0.6235779 -2.374037
OTU1324   0.4308551 0.5710941  3.574931 0.001216737 0.6235779 -2.407061
OTU737    0.3725357 0.4327240  3.207952 0.003186411 0.7800252 -2.802048
OTU506   -0.3736884 0.3950206 -3.182830 0.003399086 0.7800252 -2.828846
OTU223   -0.3510161 0.5899634 -3.084685 0.004367294 0.7800252 -2.933097
OTU243   -0.3731732 0.4269269 -3.061416 0.004632686 0.7800252 -2.957701
OTU516    0.3553084 0.6002113  2.992041 0.005517915 0.7800252 -3.030767
OTU15587 -0.2989538 0.3539695 -2.952757 0.006088002 0.7800252 -3.071934
OTU600    0.3397441 0.4367629  2.870570 0.007465734 0.8502642 -3.157527
OTU203    0.3017443 0.9295849  2.741694 0.010230400 0.8692509 -3.290101"

par(mfrow=c(2,2))
plot(get_sample(SILVA_MERDS_rar_trt_live_seed_sh,"OTU3210"),SILVA_MERDS_rar_trt_live_seed_map$days_to_germ)
plot(get_sample(SILVA_MERDS_rar_trt_live_seed_sh,"OTU1324"),SILVA_MERDS_rar_trt_live_seed_map$days_to_germ)
plot(get_sample(SILVA_MERDS_rar_trt_live_seed_sh,"OTU737"),SILVA_MERDS_rar_trt_live_seed_map$days_to_germ)
plot(get_sample(SILVA_MERDS_rar_trt_live_seed_sh,"OTU506"),SILVA_MERDS_rar_trt_live_seed_map$days_to_germ)
par(mfrow=c(2,2))
plot(get_sample(SILVA_MERDS_rar_trt_live_seed_sh,"OTU223"),SILVA_MERDS_rar_trt_live_seed_map$days_to_germ)
plot(get_sample(SILVA_MERDS_rar_trt_live_seed_sh,"OTU243"),SILVA_MERDS_rar_trt_live_seed_map$days_to_germ)
plot(get_sample(SILVA_MERDS_rar_trt_live_seed_sh,"OTU516"),SILVA_MERDS_rar_trt_live_seed_map$days_to_germ)
plot(get_sample(SILVA_MERDS_rar_trt_live_seed_sh,"OTU15587"),SILVA_MERDS_rar_trt_live_seed_map$days_to_germ)
par(mfrow=c(2,1))
plot(get_sample(SILVA_MERDS_rar_trt_live_seed_sh,"OTU600"),SILVA_MERDS_rar_trt_live_seed_map$days_to_germ)
plot(get_sample(SILVA_MERDS_rar_trt_live_seed_sh,"OTU203"),SILVA_MERDS_rar_trt_live_seed_map$days_to_germ)


SILVA_MERDS_rar_trt_live_seed_sh_sampl_sum=sample_sums(SILVA_MERDS_rar_trt_live_seed_sh)
ntaxa(SILVA_MERDS_rar_trt_live_seed_sh)
#1025
SILVA_MERDS_rar_trt_live_seed_map=sample_data(SILVA_MERDS_rar_trt_live_seed)
SILVA_MERDS_rar_trt_live_seed_map$soil_root=with(SILVA_MERDS_rar_trt_live_seed_map, interaction(soil_status,root_association))

sh_SILVA_MERDS_rar_trt_seed_sampl_sum=merge(SILVA_MERDS_rar_trt_live_seed_map,SILVA_MERDS_rar_trt_live_seed_sh_sampl_sum, by="row.names")


sh_SILVA_MERDS_rar_trt_seed_sampl_sum_precip_soil_g=sh_SILVA_MERDS_rar_trt_seed_sampl_sum %>% group_by(soil_root,precip)
sh_read_abund_precip_seed=summarise_at(sh_SILVA_MERDS_rar_trt_seed_sampl_sum_precip_soil_g, 
                                           "y", list(~n(),~mean,~sd,se=~sd(.)/sqrt(n())))

ggplot(sh_read_abund_precip_seed, aes(x=precip,y=mean,ymin = mean-se, ymax= mean+se,fill=soil_root))+
  geom_bar(position="dodge",stat="identity", color="black",aes(fill=soil_root))+
  geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+
  scale_fill_manual(values = c("light gray", "dark grey"),labels=c("Bulk","Rhizo"))+
  scale_x_discrete(labels=c("Ambient","Drought"))+ylab("Shared Taxa Read Abun")+
  geom_text(aes(y=mean+se+300, label=n),position=position_dodge(width=0.9))+theme_bw()


#Diversity of the Shared taxa
alpha_meas = c("Observed", "Shannon", "InvSimpson")
SILVA_MERDS_rar_trt_live_seed_map=sample_data(SILVA_MERDS_rar_trt_live_seed)
SILVA_MERDS_rar_trt_live_seed_map$soil_root=with(SILVA_MERDS_rar_trt_live_seed_map, interaction(soil_status,root_association))

SILVA_MERDS_rar_trt_live_seed_sh.divfil=estimate_richness(SILVA_MERDS_rar_trt_live_seed_sh,measures=alpha_meas)

SILVA_MERDS_rar_trt_live_seed_sh.divfil=merge(SILVA_MERDS_rar_trt_live_seed_sh.divfil, SILVA_MERDS_rar_trt_live_seed_map, by ="row.names")
#bact.soilE.t.divfil=mutate(bact.soilE.t.divfil, pielou=Shannon*(1/log(Observed)))
head(SILVA_MERDS_rar_trt_live_seed_sh.divfil)
row.names(SILVA_MERDS_rar_trt_live_seed_sh.divfil)=SILVA_MERDS_rar_trt_live_seed_sh.divfil$Row.names
SILVA_MERDS_rar_trt_live_seed_sh.divfil$Row.names=NULL



SILVA_MERDS_rar_trt_live_seed_sh.divfil_g=SILVA_MERDS_rar_trt_live_seed_sh.divfil %>% group_by(soil_root,precip)
sh_inv_simp_precip_seed=summarise_at(SILVA_MERDS_rar_trt_live_seed_sh.divfil_g, 
                                         "InvSimpson", list(~n(),~mean,~sd,se=~sd(.)/sqrt(n())))

ggplot(sh_inv_simp_precip_seed, aes(x=precip,y=mean,ymin = mean-se, ymax= mean+se,fill=soil_root))+
  geom_bar(position="dodge",stat="identity", color="black",aes(fill=soil_root))+
  geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+
  scale_fill_manual(values = c("light gray", "dark grey"),labels=c("Bulk","Rhizo"))+
  scale_x_discrete(labels=c("Ambient","Drought"))+ylab("Shared inv Simpson")+
  geom_text(aes(y=mean+se+5, label=n),position=position_dodge(width=0.9))+theme_bw()


######Seed DOMINANT TAXA Sterile and live####

SILVA_MERDS_rar_seed=subset_samples(SILVA_MERDS_rar, life_stage=="S")
summary(sample_data(SILVA_MERDS_rar_seed))
SILVA_MERDS_rar_seed_shared=core(SILVA_MERDS_rar_seed, detection = 0,prevalence = 0.9999999999)

SILVA_MERDS_rar_seed_shared_sampl_sum=sample_sums(SILVA_MERDS_rar_seed_shared)
ntaxa(SILVA_MERDS_rar_seed_shared)
#10
SILVA_MERDS_rar_seed_shared_map=sample_data(SILVA_MERDS_rar_seed_shared)
SILVA_MERDS_rar_seed_shared_map$soil_root=with(SILVA_MERDS_rar_seed_shared_map, interaction(soil_status,root_association))

shared_SILVA_MERDS_rar_seed_sampl_sum=merge(SILVA_MERDS_rar_seed_shared_map,SILVA_MERDS_rar_seed_shared_sampl_sum, by="row.names")




treatment_order=c("S.B","L.B","L.R")
shared_SILVA_MERDS_rar_seed_sampl_sum_precip_soil_g=shared_SILVA_MERDS_rar_seed_sampl_sum %>% group_by(soil_root,precip)
shared_read_abund_precip_seed_W_stl=summarise_at(shared_SILVA_MERDS_rar_seed_sampl_sum_precip_soil_g, 
                                                  "y", list(~n(),~mean(.),~sd(.),se=~sd(.)/sqrt(n())))

(core_trans_taxa=ggplot(shared_read_abund_precip_seed_W_stl, aes(x=precip,y=mean, ymin = mean-se, ymax= mean+se,fill=factor(soil_root,levels = treatment_order)))+
    geom_bar(position="dodge",stat="identity", color="black",aes(fill=factor(soil_root,levels = treatment_order)))+
    geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+
    scale_fill_manual(values = c( "white","light gray", "dark grey"),labels=c("Sterile","Bulk","Rhizo"), limits=treatment_order, name=NULL)+
    scale_x_discrete(labels=c("Ambient","Drought"), name=NULL)+ylab("Core taxa read # (in all samples)")+
    geom_text(aes(y=500, label=n),size=10,position=position_dodge(width=0.9))+theme_bw()+theme(axis.title = element_text(size = 23), axis.text = element_text(size = 23),
                                                                                               legend.position = c(0.85,.9), legend.text=element_text(size=20),
                                                                                               legend.background = element_rect(size=0.5,linetype="solid",colour ="black")))


#Diversity of the Shared taxa
alpha_meas = c("Observed", "Shannon", "InvSimpson")
SILVA_MERDS_rar_seed_shared_map=sample_data(SILVA_MERDS_rar_seed_shared)
SILVA_MERDS_rar_seed_shared_map$soil_root=with(SILVA_MERDS_rar_seed_shared_map, interaction(soil_status,root_association))

SILVA_MERDS_rar_seed_shared.divfil=estimate_richness(SILVA_MERDS_rar_seed_shared,measures=alpha_meas)

SILVA_MERDS_rar_seed_shared.divfil=merge(SILVA_MERDS_rar_seed_shared.divfil, SILVA_MERDS_rar_seed_shared_map, by ="row.names")
#bact.soilE.t.divfil=mutate(bact.soilE.t.divfil, pielou=Shannon*(1/log(Observed)))
head(SILVA_MERDS_rar_seed_shared.divfil)
row.names(SILVA_MERDS_rar_seed_shared.divfil)=SILVA_MERDS_rar_seed_shared.divfil$Row.names
SILVA_MERDS_rar_seed_shared.divfil$Row.names=NULL



SILVA_MERDS_rar_seed_shared.divfil_g=SILVA_MERDS_rar_seed_shared.divfil %>% group_by(soil_root,precip)
shared_inv_simp_precip_seed_w_stl=summarise_at(SILVA_MERDS_rar_seed_shared.divfil_g, 
                                                "InvSimpson", list(~n(),~mean(.),~sd(.),se=~sd(.)/sqrt(n())))

ggplot(shared_inv_simp_precip_seed_w_stl, aes(x=precip,y=mean,ymin = mean-se, ymax= mean+se,fill=factor(soil_root,levels = treatment_order)))+
  geom_bar(position="dodge",stat="identity", color="black",aes(fill=factor(soil_root,levels = treatment_order)))+
  geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+
  scale_fill_manual(values = c( "white","light gray", "dark grey"),labels=c("Sterile","Bulk","Rhizo"), limits=treatment_order, name=NULL)+
  scale_x_discrete(labels=c("Ambient","Drought"), name=NULL)+ylab("Core diversity invSimpson \n(in all samples)")+
  geom_text(aes(y=1, label=n),size=10,position=position_dodge(width=0.9))+theme_bw()+
  theme(axis.title = element_text(size = 23), axis.text = element_text(size = 23),
        legend.position = c(0.85,.9), legend.text=element_text(size=20),
        legend.background = element_rect(size=0.5,linetype="solid",colour ="black"))+
  ylim(c(0,8.5))





#####End SILVA dataset v123####

#####Does using RDP change the numbers of species that are Classified as non bacterial####


load("D:/MERDS_2018/merds/Switchgrass/R_data/phyl_obj_RDP_MERDS.RData")

ntaxa(phyl_RDP_MERDS)
#11892
sum(taxa_sums(phyl_RDP_MERDS))
#980588
mean(sample_sums(phyl_RDP_MERDS))
#17510.5
min(sample_sums(phyl_RDP_MERDS))
#242
max(sample_sums(phyl_RDP_MERDS))
#30633
sort(sample_sums(phyl_RDP_MERDS))
#  MERDSSG42  MERDSSG85
#       242      10703

RDP_MERDS_data<-subset_taxa(phyl_RDP_MERDS,Class!="c:Chloroplast")
ntaxa(RDP_MERDS_data)
#11851 RDP
#11797 silva
sum(otu_table(RDP_MERDS_data))
#979020 RDP
#980577 silva

#THE RDP DATABASE DOES NOT HAVE MITOCHONDRIA. BLAH!!!!!!!

RDP_MERDS_data<-subset_taxa(RDP_MERDS_data,Family!="f:Mitochondria")
ntaxa(RDP_MERDS_data)
#11851 RDP
#11774 silva
sum(otu_table(RDP_MERDS_data))
#979020 RDP
#977652 silva

sort(sample_sums(RDP_MERDS_data))
# MERDSSG42  MERDSSG85 MERDSSG103
#       242      10697      11519

#THERE ARE not that different.
#####End OTU based community####


#####Begin ZOTU based community####

#Let's load in the SILVA 123 version of the phyloseq object
load("D:/MERDS_2018/merds/Switchgrass/R_data/Zotu_phyl_obj_SILVA_MERDS.RData")

ntaxa(Zotu_phyl_SILVA_MERDS)
#25557
sum(taxa_sums(Zotu_phyl_SILVA_MERDS))
#993801
mean(sample_sums(Zotu_phyl_SILVA_MERDS))
#17746.45
min(sample_sums(Zotu_phyl_SILVA_MERDS))
#244
max(sample_sums(Zotu_phyl_SILVA_MERDS))
#31305
sort(sample_sums(Zotu_phyl_SILVA_MERDS))
# MERDSSG42  MERDSSG85 MERDSSG103
#       244      10818      11682

#I updated the mapping file with the treatments from the experiment

merds_map_trt=read.csv("D:/MERDS_2018/merds/Switchgrass/R_data/map_MERDS_trt.csv", header = T, row.names = 1)
summary(merds_map_trt)
#Need to convert plant number to numeric
merds_map_trt$Plant_Number <- as.numeric(levels(merds_map_trt$Plant_Number)) [merds_map_trt$Plant_Number] # changing factor to numeric
summary(merds_map_trt)

#Let's also add in the biomass data
merds_map_trt_bio=merge(merds_map_trt,data_SG_biomass, by="Plant_Number",all.x = T)
summary(merds_map_trt_bio)
merds_map_trt_bio$surv_germ=merds_map_trt_bio$shoot_weight_g
merds_map_trt_bio$surv_germ[merds_map_trt_bio$surv_germ>0]=1
merds_map_trt_bio$surv_germ[is.na(merds_map_trt_bio$surv_germ)]=0
merds_map_trt_bio$total_biomass[is.na(merds_map_trt_bio$total_biomass)]=0
merds_map_trt_bio$shoot_weight_g[is.na(merds_map_trt_bio$shoot_weight_g)]=0
merds_map_trt_bio$root_weight_g[is.na(merds_map_trt_bio$root_weight_g)]=0
row.names(merds_map_trt_bio)=merds_map_trt_bio$SampleID




#####Begin SILVA dataset v123####
#now put the new mapping file into our phyloseq obj

Zotu_SILVA_MERDS_data=phyloseq(otu_table(Zotu_phyl_SILVA_MERDS),tax_table(Zotu_phyl_SILVA_MERDS),sample_data(merds_map_trt_bio))
nrow(sample_data(Zotu_SILVA_MERDS_data))
#56
ntaxa(Zotu_SILVA_MERDS_data)
#25557
sum(taxa_sums(Zotu_SILVA_MERDS_data))
#993801
mean(sample_sums(Zotu_SILVA_MERDS_data))
#17746.45
min(sample_sums(Zotu_SILVA_MERDS_data))
#244
max(sample_sums(Zotu_SILVA_MERDS_data))
#31305
sort(sample_sums(Zotu_SILVA_MERDS_data))
# MERDSSG42  MERDSSG85 MERDSSG103
#       244      10818      11682

get_taxa_unique(Zotu_SILVA_MERDS_data, taxonomic.rank = "Domain")
Zotu_SILVA_MERDS_data<-subset_taxa(Zotu_SILVA_MERDS_data,Domain!="")
ntaxa(Zotu_SILVA_MERDS_data)
#25550
sum(taxa_sums(Zotu_SILVA_MERDS_data))
#993701

Zotu_SILVA_MERDS_data<-subset_taxa(Zotu_SILVA_MERDS_data,Class!="c:Chloroplast")
ntaxa(Zotu_SILVA_MERDS_data)
#25415
sum(otu_table(Zotu_SILVA_MERDS_data))
#991458

Zotu_SILVA_MERDS_data<-subset_taxa(Zotu_SILVA_MERDS_data,Family!="f:Mitochondria")
ntaxa(Zotu_SILVA_MERDS_data)
#25384
sum(otu_table(Zotu_SILVA_MERDS_data))
#990782
sort(sample_sums(Zotu_SILVA_MERDS_data))
# MERDSSG42  MERDSSG85 MERDSSG103
#       244      10808      11642

#remove sample with 239 read
Zotu_SILVA_MERDS_trunc=prune_samples(sample_sums(Zotu_SILVA_MERDS_data) > 10000, Zotu_SILVA_MERDS_data)
nrow(sample_data(Zotu_SILVA_MERDS_trunc))
#55

Zotu_SILVA_MERDS_trunc3=prune_taxa(taxa_sums(Zotu_SILVA_MERDS_trunc) > 2, Zotu_SILVA_MERDS_trunc)
ntaxa(Zotu_SILVA_MERDS_trunc3)
#18683
sum(otu_table(Zotu_SILVA_MERDS_trunc3))
#980751

#let's look at raw ordination

Zotu_SILVA_MERDS_ord=ordinate(Zotu_SILVA_MERDS_trunc3, method = "NMDS",distance = "bray")
#*** Solution reached
#Warning message:
#  In metaMDS(veganifyOTU(physeq), distance, ...) :
#  stress is (nearly) zero: you may have insufficient data
#8.084044e-05
plot_ordination(Zotu_SILVA_MERDS_trunc3,Zotu_SILVA_MERDS_ord, color="root_association",shape="life_stage", label = "block")



alpha_meas = c("Observed", "Chao1", "Shannon", "InvSimpson")
Zotu_SILVA_MERDS_trunc3_map=sample_data(Zotu_SILVA_MERDS_trunc3)
Zotu_SILVA_MERDS_trunc3.divfil=estimate_richness(Zotu_SILVA_MERDS_trunc3,measures=alpha_meas)

Zotu_SILVA_MERDS_trunc3.divfil=merge(Zotu_SILVA_MERDS_trunc3.divfil, Zotu_SILVA_MERDS_trunc3_map, by ="row.names")
#bact.soilE.t.divfil=mutate(bact.soilE.t.divfil, pielou=Shannon*(1/log(Observed)))
head(Zotu_SILVA_MERDS_trunc3.divfil)
row.names(Zotu_SILVA_MERDS_trunc3.divfil)=Zotu_SILVA_MERDS_trunc3.divfil$Row.names
Zotu_SILVA_MERDS_trunc3.divfil$Row.names=NULL
Zotu_SILVA_MERDS_trunc3.divfil=merge(Zotu_SILVA_MERDS_trunc3.divfil, sample_sums(Zotu_SILVA_MERDS_trunc3), by ="row.names")
head(Zotu_SILVA_MERDS_trunc3.divfil)
colnames(Zotu_SILVA_MERDS_trunc3.divfil)[colnames(Zotu_SILVA_MERDS_trunc3.divfil)=="y"]="read_abun"
row.names(Zotu_SILVA_MERDS_trunc3.divfil)=Zotu_SILVA_MERDS_trunc3.divfil$Row.names
Zotu_SILVA_MERDS_trunc3.divfil$Row.names=NULL


ggplot(Zotu_SILVA_MERDS_trunc3.divfil, aes(x=soil_status, y=read_abun))+geom_boxplot(aes(color=interaction(root_association,precip,life_stage)))+theme_bw()


Zotu_SILVA_MERDS_trunc_sterile=subset_samples(Zotu_SILVA_MERDS_trunc3, soil_status=="S")


Zotu_SILVA_MERDS_sterile_ord=ordinate(Zotu_SILVA_MERDS_trunc_sterile, method = "NMDS",distance = "bray")
#*** Solution reached
#0.08564865 
plot_ordination(Zotu_SILVA_MERDS_trunc_sterile,Zotu_SILVA_MERDS_sterile_ord, color="precip",shape="life_stage", label = "block")

#see if there is a correlation between biomass and OTUs
Zotu_SILVA_MERDS_trunc_sterile_core=core(Zotu_SILVA_MERDS_trunc_sterile, detection = 0,prevalence = .75)
sample_sums(Zotu_SILVA_MERDS_trunc_sterile_core)
Zotu_SILVA_MERDS_trunc_sterile_core_map=sample_data(Zotu_SILVA_MERDS_trunc_sterile_core)

Zotu_biomas_mod_otus=lm_phyloseq(Zotu_SILVA_MERDS_trunc_sterile_core, "total_biomass")
#Warning message:
#In transform(x, transformation) :
#  OTU table contains zeroes. Using log10(1 + x) transform.
"              logFC   AveExpr         t    P.Value adj.P.Val         B
Zotu265   0.9450491 1.0565828  2.723638 0.01409950 0.6344773 -4.380538
Zotu1167  0.6450591 1.3997288  1.981166 0.06336253 0.7785255 -4.483191
Zotu1564  0.5367404 0.7717864  1.979673 0.06354324 0.7785255 -4.483392
Zotu68   -0.8068861 2.2968545 -1.915769 0.07172824 0.7785255 -4.491944
Zotu463  -0.5428312 0.6463718 -1.815268 0.08650284 0.7785255 -4.505194
Zotu13   -0.5679574 2.7558109 -1.642421 0.11818834 0.8864126 -4.527279
Zotu243   0.4655628 1.1270367  1.456382 0.16283627 0.9313585 -4.549793
Zotu1     0.2526384 2.4046019  1.328676 0.20086860 0.9313585 -4.564321
Zotu689   0.4823702 0.9289082  1.273984 0.21918020 0.9313585 -4.570278
Zotu6904 -0.3797663 0.6289936 -1.246745 0.22877384 0.9313585 -4.573182"
par(mfrow=c(2,2))
plot(get_sample(Zotu_SILVA_MERDS_trunc_sterile_core,"Zotu265"),Zotu_SILVA_MERDS_trunc_sterile_core_map$total_biomass)
plot(get_sample(Zotu_SILVA_MERDS_trunc_sterile_core,"Zotu1167"),Zotu_SILVA_MERDS_trunc_sterile_core_map$total_biomass)
plot(get_sample(Zotu_SILVA_MERDS_trunc_sterile_core,"Zotu1564"),Zotu_SILVA_MERDS_trunc_sterile_core_map$total_biomass)
plot(get_sample(Zotu_SILVA_MERDS_trunc_sterile_core,"Zotu68"),Zotu_SILVA_MERDS_trunc_sterile_core_map$total_biomass)
par(mfrow=c(2,2))
plot(get_sample(Zotu_SILVA_MERDS_trunc_sterile_core,"Zotu463"),Zotu_SILVA_MERDS_trunc_sterile_core_map$total_biomass)
plot(get_sample(Zotu_SILVA_MERDS_trunc_sterile_core,"Zotu13"),Zotu_SILVA_MERDS_trunc_sterile_core_map$total_biomass)
plot(get_sample(Zotu_SILVA_MERDS_trunc_sterile_core,"Zotu243"),Zotu_SILVA_MERDS_trunc_sterile_core_map$total_biomass)
plot(get_sample(Zotu_SILVA_MERDS_trunc_sterile_core,"Zotu1"),Zotu_SILVA_MERDS_trunc_sterile_core_map$total_biomass)
par(mfrow=c(2,1))
plot(get_sample(Zotu_SILVA_MERDS_trunc_sterile_core,"Zotu689"),Zotu_SILVA_MERDS_trunc_sterile_core_map$total_biomass)
plot(get_sample(Zotu_SILVA_MERDS_trunc_sterile_core,"Zotu6904"),Zotu_SILVA_MERDS_trunc_sterile_core_map$total_biomass)


#How does this compare to the the OTU abundance versus biomass in live soil
Zotu_SILVA_MERDS_trunc3_live_trt=subset_samples(Zotu_SILVA_MERDS_trunc3, soil_status=="L"&life_stage!="Start")
nrow(sample_data(Zotu_SILVA_MERDS_trunc3_live_trt))
#32
summary(sample_data(Zotu_SILVA_MERDS_trunc3_live_trt))
Zotu_SILVA_MERDS_trunc3_live_trt_map=sample_data(Zotu_SILVA_MERDS_trunc3_live_trt)
par(mfrow=c(2,2))
plot(get_sample(Zotu_SILVA_MERDS_trunc3_live_trt,"Zotu265"),Zotu_SILVA_MERDS_trunc3_live_trt_map$total_biomass)
plot(get_sample(Zotu_SILVA_MERDS_trunc3_live_trt,"Zotu1167"),Zotu_SILVA_MERDS_trunc3_live_trt_map$total_biomass)
plot(get_sample(Zotu_SILVA_MERDS_trunc3_live_trt,"Zotu1564"),Zotu_SILVA_MERDS_trunc3_live_trt_map$total_biomass)
plot(get_sample(Zotu_SILVA_MERDS_trunc3_live_trt,"Zotu68"),Zotu_SILVA_MERDS_trunc3_live_trt_map$total_biomass)
par(mfrow=c(2,2))
plot(get_sample(Zotu_SILVA_MERDS_trunc3_live_trt,"Zotu463"),Zotu_SILVA_MERDS_trunc3_live_trt_map$total_biomass)
plot(get_sample(Zotu_SILVA_MERDS_trunc3_live_trt,"Zotu13"),Zotu_SILVA_MERDS_trunc3_live_trt_map$total_biomass)
plot(get_sample(Zotu_SILVA_MERDS_trunc3_live_trt,"Zotu243"),Zotu_SILVA_MERDS_trunc3_live_trt_map$total_biomass)
plot(get_sample(Zotu_SILVA_MERDS_trunc3_live_trt,"Zotu1"),Zotu_SILVA_MERDS_trunc3_live_trt_map$total_biomass)
par(mfrow=c(2,1))
plot(get_sample(Zotu_SILVA_MERDS_trunc3_live_trt,"Zotu689"),Zotu_SILVA_MERDS_trunc3_live_trt_map$total_biomass)
plot(get_sample(Zotu_SILVA_MERDS_trunc3_live_trt,"Zotu6904"),Zotu_SILVA_MERDS_trunc3_live_trt_map$total_biomass)


#how frequent are these OTUs in the rhizo and bulk starting soil

Zotu_SILVA_MERDS_trunc3_start=subset_samples(Zotu_SILVA_MERDS_trunc3, life_stage=="Start")
Zotu_biomass_indic=c("Zotu265","Zotu1167","Zotu1564","Zotu68","Zotu463","Zotu13","Zotu243","Zotu1","Zotu689","Zotu6904")
Zotu_SILVA_MERDS_trunc3_start_biomss_otu=t(otu_table(prune_taxa(Zotu_biomass_indic,Zotu_SILVA_MERDS_trunc3_start)))

Zotu_SILVA_MERDS_trunc3_start_biomss_otu_trt=merge(Zotu_SILVA_MERDS_trunc3_start_biomss_otu,sample_data(Zotu_SILVA_MERDS_trunc3_start),by="row.names")

grid.arrange(ggplot(Zotu_SILVA_MERDS_trunc3_start_biomss_otu_trt,aes(x=root_association,y=Zotu265))+geom_boxplot(),
             ggplot(Zotu_SILVA_MERDS_trunc3_start_biomss_otu_trt,aes(x=root_association,y=Zotu1167))+geom_boxplot(),
             ggplot(Zotu_SILVA_MERDS_trunc3_start_biomss_otu_trt,aes(x=root_association,y=Zotu1564))+geom_boxplot(),
             ggplot(Zotu_SILVA_MERDS_trunc3_start_biomss_otu_trt,aes(x=root_association,y=Zotu68))+geom_boxplot(), nrow=2,ncol=2)

grid.arrange(ggplot(Zotu_SILVA_MERDS_trunc3_start_biomss_otu_trt,aes(x=root_association,y=Zotu463))+geom_boxplot(),
             ggplot(Zotu_SILVA_MERDS_trunc3_start_biomss_otu_trt,aes(x=root_association,y=Zotu13))+geom_boxplot(),
             ggplot(Zotu_SILVA_MERDS_trunc3_start_biomss_otu_trt,aes(x=root_association,y=Zotu243))+geom_boxplot(),
             ggplot(Zotu_SILVA_MERDS_trunc3_start_biomss_otu_trt,aes(x=root_association,y=Zotu1))+geom_boxplot(), nrow=2,ncol=2)


grid.arrange(ggplot(Zotu_SILVA_MERDS_trunc3_start_biomss_otu_trt,aes(x=root_association,y=Zotu689))+geom_boxplot(),
             ggplot(Zotu_SILVA_MERDS_trunc3_start_biomss_otu_trt,aes(x=root_association,y=Zotu6904))+geom_boxplot(), nrow=2)

#there is a large outlier in both datasets that maybe skewing the results

#see if there is a correlation between biomass and OTUs
SILVA_MERDS_trunc_sterile_out=subset_samples(SILVA_MERDS_trunc_sterile, total_biomass<1)
SILVA_MERDS_trunc_sterile_out_core=core(SILVA_MERDS_trunc_sterile_out, detection = 0,prevalence = .75)
sample_sums(SILVA_MERDS_trunc_sterile_out_core)
SILVA_MERDS_trunc_sterile_out_core_map=sample_data(SILVA_MERDS_trunc_sterile_out_core)

biomas_mod_otus_otu=lm_phyloseq(SILVA_MERDS_trunc_sterile_out_core, "total_biomass")
#Warning message:
#In transform(x, transformation) :
#  OTU table contains zeroes. Using log10(1 + x) transform.
"             logFC   AveExpr         t    P.Value adj.P.Val         B
OTU148  -1.5066849 1.0589958 -2.689344 0.01563565 0.7661469 -4.533616
OTU2019  1.4267018 0.8926121  2.254355 0.03783274 0.7811411 -4.550973
OTU115  -2.1834756 1.6589504 -2.089652 0.05218602 0.7811411 -4.557522
OTU139  -1.6435907 1.2320660 -1.894277 0.07554546 0.7811411 -4.565162
OTU502   1.0398543 1.3559353  1.727591 0.10240739 0.7811411 -4.571497
OTU63   -1.2931896 2.1090492 -1.596441 0.12903756 0.7811411 -4.576314
OTU900   0.6975020 1.2261122  1.470375 0.15995043 0.7811411 -4.580770
OTU106   1.0074927 2.2038149  1.405929 0.17798032 0.7811411 -4.582970
OTU4     1.2710933 1.9247360  1.399929 0.17973994 0.7811411 -4.583172
OTU2839  0.8381175 0.9993585  1.325547 0.20274307 0.7811411 -4.585633"
par(mfrow=c(2,2))
plot(get_sample(SILVA_MERDS_trunc_sterile_out_core,"OTU148"),SILVA_MERDS_trunc_sterile_out_core_map$total_biomass)
plot(get_sample(SILVA_MERDS_trunc_sterile_out_core,"OTU2019"),SILVA_MERDS_trunc_sterile_out_core_map$total_biomass)
plot(get_sample(SILVA_MERDS_trunc_sterile_out_core,"OTU115"),SILVA_MERDS_trunc_sterile_out_core_map$total_biomass)
plot(get_sample(SILVA_MERDS_trunc_sterile_out_core,"OTU139"),SILVA_MERDS_trunc_sterile_out_core_map$total_biomass)
par(mfrow=c(2,2))
plot(get_sample(SILVA_MERDS_trunc_sterile_out_core,"OTU502"),SILVA_MERDS_trunc_sterile_out_core_map$total_biomass)
plot(get_sample(SILVA_MERDS_trunc_sterile_out_core,"OTU63"),SILVA_MERDS_trunc_sterile_out_core_map$total_biomass)
plot(get_sample(SILVA_MERDS_trunc_sterile_out_core,"OTU900"),SILVA_MERDS_trunc_sterile_out_core_map$total_biomass)
plot(get_sample(SILVA_MERDS_trunc_sterile_out_core,"OTU106"),SILVA_MERDS_trunc_sterile_out_core_map$total_biomass)
par(mfrow=c(2,1))
plot(get_sample(SILVA_MERDS_trunc_sterile_out_core,"OTU4"),SILVA_MERDS_trunc_sterile_out_core_map$total_biomass)
plot(get_sample(SILVA_MERDS_trunc_sterile_out_core,"OTU2839"),SILVA_MERDS_trunc_sterile_out_core_map$total_biomass)


#How does this compare to the the OTU abundance versus biomass in live soil
SILVA_MERDS_trunc3_live_trt=subset_samples(SILVA_MERDS_trunc3, soil_status=="L"&life_stage!="Start")
SILVA_MERDS_trunc3_live_trt_out=subset_samples(SILVA_MERDS_trunc3_live_trt, total_biomass<1)
nrow(sample_data(SILVA_MERDS_trunc3_live_trt_out))
#31
summary(sample_data(SILVA_MERDS_trunc3_live_trt_out))
SILVA_MERDS_trunc3_live_trt_out_map=sample_data(SILVA_MERDS_trunc3_live_trt_out)
par(mfrow=c(2,2))
plot(get_sample(SILVA_MERDS_trunc3_live_trt_out,"OTU148"),SILVA_MERDS_trunc3_live_trt_out_map$total_biomass)
plot(get_sample(SILVA_MERDS_trunc3_live_trt_out,"OTU2019"),SILVA_MERDS_trunc3_live_trt_out_map$total_biomass)
plot(get_sample(SILVA_MERDS_trunc3_live_trt_out,"OTU115"),SILVA_MERDS_trunc3_live_trt_out_map$total_biomass)
plot(get_sample(SILVA_MERDS_trunc3_live_trt_out,"OTU139"),SILVA_MERDS_trunc3_live_trt_out_map$total_biomass)
par(mfrow=c(2,2))
plot(get_sample(SILVA_MERDS_trunc3_live_trt_out,"OTU502"),SILVA_MERDS_trunc3_live_trt_out_map$total_biomass)
plot(get_sample(SILVA_MERDS_trunc3_live_trt_out,"OTU63"),SILVA_MERDS_trunc3_live_trt_out_map$total_biomass)
plot(get_sample(SILVA_MERDS_trunc3_live_trt_out,"OTU900"),SILVA_MERDS_trunc3_live_trt_out_map$total_biomass)
plot(get_sample(SILVA_MERDS_trunc3_live_trt_out,"OTU106"),SILVA_MERDS_trunc3_live_trt_out_map$total_biomass)
par(mfrow=c(2,1))
plot(get_sample(SILVA_MERDS_trunc3_live_trt_out,"OTU4"),SILVA_MERDS_trunc3_live_trt_out_map$total_biomass)
plot(get_sample(SILVA_MERDS_trunc3_live_trt_out,"OTU2839"),SILVA_MERDS_trunc3_live_trt_out_map$total_biomass)




#Let look at transplants and see if there is a correlation between biomass and OTUs
sample_data(SILVA_MERDS_trunc_sterile)
SILVA_MERDS_trunc_sterile_trans=subset_samples(SILVA_MERDS_trunc_sterile, life_stage=="G")
sample_data(SILVA_MERDS_trunc_sterile_trans)

#Only taxa in 75% of samples
SILVA_MERDS_trunc_sterile_trans_core=core(SILVA_MERDS_trunc_sterile_trans, detection = 0,prevalence = .75)
sample_sums(SILVA_MERDS_trunc_sterile_trans_core)


biomas_mod_otus=lm_phyloseq(SILVA_MERDS_trunc_sterile_trans_core, "total_biomass")
#Warning message:
#In transform(x, transformation) :
#  OTU table contains zeroes. Using log10(1 + x) transform.
"              logFC   AveExpr         t     P.Value adj.P.Val         B
OTU883   -1.0760738 2.3961867 -3.759898 0.002902171 0.1704741 -1.440618
OTU43     1.0106521 0.7150720  3.524517 0.004427898 0.1704741 -1.797326
OTU78     0.8413813 0.7579017  2.903998 0.013698904 0.3036855 -2.757026
OTU12410  0.7028079 0.6318266  2.717903 0.019242181 0.3036855 -3.045739
OTU57    -0.7719324 2.2764640 -2.632344 0.022486918 0.3036855 -3.177849
OTU193    1.4985374 1.5014569  2.604288 0.023663804 0.3036855 -3.221041
OTU552    0.8886797 1.3904565  2.336573 0.038380859 0.3681663 -3.628536
OTU377   -0.6054359 1.8774871 -2.255027 0.044403555 0.3681663 -3.750422
OTU2833  -0.6852767 1.2696164 -2.219543 0.047297452 0.3681663 -3.803045
OTU242   -0.5914900 0.6902468 -2.213430 0.047813799 0.3681663 -3.812083"

SILVA_MERDS_trunc_sterile_trans_core_map=sample_data(SILVA_MERDS_trunc_sterile_trans_core)
SILVA_MERDS_trunc_sterile_trans_core_map$total_biomass

plot(get_sample(SILVA_MERDS_rar_sterile_trans,"OTU883"),SILVA_MERDS_trunc_sterile_trans_core_map$total_biomass)
plot(get_sample(SILVA_MERDS_rar_sterile_trans,"OTU43"),SILVA_MERDS_trunc_sterile_trans_core_map$total_biomass)
plot(get_sample(SILVA_MERDS_rar_sterile_trans,"OTU78"),SILVA_MERDS_trunc_sterile_trans_core_map$total_biomass)
plot(get_sample(SILVA_MERDS_rar_sterile_trans,"OTU12410"),SILVA_MERDS_trunc_sterile_trans_core_map$total_biomass)
plot(get_sample(SILVA_MERDS_rar_sterile_trans,"OTU57"),SILVA_MERDS_trunc_sterile_trans_core_map$total_biomass)
plot(get_sample(SILVA_MERDS_rar_sterile_trans,"OTU193"),SILVA_MERDS_trunc_sterile_trans_core_map$total_biomass)
plot(get_sample(SILVA_MERDS_rar_sterile_trans,"OTU552"),SILVA_MERDS_trunc_sterile_trans_core_map$total_biomass)
plot(get_sample(SILVA_MERDS_rar_sterile_trans,"OTU377"),SILVA_MERDS_trunc_sterile_trans_core_map$total_biomass)
plot(get_sample(SILVA_MERDS_rar_sterile_trans,"OTU2833"),SILVA_MERDS_trunc_sterile_trans_core_map$total_biomass)
plot(get_sample(SILVA_MERDS_rar_sterile_trans,"OTU242"),SILVA_MERDS_trunc_sterile_trans_core_map$total_biomass)



#Let's rarefy the data
Zotu_SILVA_MERDS_rar=rarefy_even_depth(Zotu_SILVA_MERDS_trunc3, sample.size= 10000, rngseed = T)
#1360OTUs were removed because they are no longer 
#present in any sample after random subsampling

Zotu_SILVA_MERDS_rar_ord=ordinate(Zotu_SILVA_MERDS_rar, method = "NMDS",distance = "bray")
#*** Solution reached
#Warning message:
#  In metaMDS(veganifyOTU(physeq), distance, ...) :
#  stress is (nearly) zero: you may have insufficient data
#7.894471e-05
plot_ordination(Zotu_SILVA_MERDS_rar,Zotu_SILVA_MERDS_rar_ord, color="root_association",shape="life_stage")+geom_point(size=3)+
  theme_bw()



Zotu_SILVA_MERDS_rar_map=sample_data(Zotu_SILVA_MERDS_rar)
Zotu_SILVA_MERDS_rar_map$soil_root_stage=with(Zotu_SILVA_MERDS_rar_map, interaction(soil_status,root_association,life_stage))
Zotu_SILVA_MERDS_rar_dis=distance(Zotu_SILVA_MERDS_rar,method = "bray")

adonis(Zotu_SILVA_MERDS_rar_dis~Zotu_SILVA_MERDS_rar_map$soil_root_stage+as.factor(Zotu_SILVA_MERDS_rar_map$block), permutations = 9999)
#                                     Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#Zotu_SILVA_MERDS_rar_map$soil_root_stage   7    7.3762 1.05375  5.6911 0.44845 0.0001 ***
#as.factor(Zotu_SILVA_MERDS_rar_map$block)  3    0.9251 0.30837  1.6654 0.05624 0.0195 *  
#Residuals                                 44    8.1469 0.18516         0.49531           
#Total                                     54   16.4483                 1.00000  

pairwise.perm.manova(Zotu_SILVA_MERDS_rar_dis, Zotu_SILVA_MERDS_rar_map$soil_root_stage, nperm=2000)



"	Pairwise comparisons using permutation MANOVAs on a distance matrix 

data:  Zotu_SILVA_MERDS_rar_dis by Zotu_SILVA_MERDS_rar_map$soil_root_stage
2000 permutations 

L.B.G  S.B.G  L.R.G  L.B.S  S.B.S  L.R.S  L.B.Start
S.B.G     0.0020 -      -      -      -      -      -        
L.R.G     0.8266 0.0020 -      -      -      -      -        
L.B.S     0.1797 0.0020 0.1108 -      -      -      -        
S.B.S     0.0023 0.0020 0.0023 0.0020 -      -      -        
L.R.S     0.0779 0.0023 0.3261 0.0020 0.0023 -      -        
L.B.Start 0.0037 0.0037 0.0037 0.0023 0.0053 0.0037 -        
L.R.Start 0.0037 0.0037 0.0042 0.0020 0.0070 0.0037 0.9220     

P value adjustment method: fdr "


#Remove the sterile since it seems to be driving a lot of this

Zotu_SILVA_MERDS_rar_live=subset_samples(Zotu_SILVA_MERDS_rar, soil_status!="S")

Zotu_SILVA_MERDS_rar_live_ord=ordinate(Zotu_SILVA_MERDS_rar_live, method = "NMDS",distance = "bray")
#*** Solution reached
#0.1439837 
plot_ordination(Zotu_SILVA_MERDS_rar_live,Zotu_SILVA_MERDS_rar_live_ord, color="root_association",shape="precip")+geom_point(size=3)+
  geom_label_repel(size=3,aes(label = block))+theme_bw()


Zotu_SILVA_MERDS_rar_live_map=sample_data(Zotu_SILVA_MERDS_rar_live)
Zotu_SILVA_MERDS_rar_live_map$soil_root=with(Zotu_SILVA_MERDS_rar_live_map, interaction(soil_status,root_association))
summary(Zotu_SILVA_MERDS_rar_live_map)
Zotu_SILVA_MERDS_rar_live_dis=distance(Zotu_SILVA_MERDS_rar_live,method = "bray")

adonis(Zotu_SILVA_MERDS_rar_live_dis~Zotu_SILVA_MERDS_rar_live_map$root_association+Zotu_SILVA_MERDS_rar_live_map$precip+Zotu_SILVA_MERDS_rar_live_map$life_stage
       +as.factor(Zotu_SILVA_MERDS_rar_live_map$block), permutations = 9999)
#Zotu_SILVA_MERDS_rar_live_map$root_association  1    0.3069 0.30688  1.7880 0.03448 0.0230 *  
#Zotu_SILVA_MERDS_rar_live_map$precip            2    1.8234 0.91171  5.3119 0.20486 0.0001 ***
#as.factor(Zotu_SILVA_MERDS_rar_live_map$block)  3    1.0696 0.35653  2.0772 0.12016 0.0003 ***



#Only treatments

SILVA_MERDS_rar_trt=subset_samples(SILVA_MERDS_rar, life_stage!="Start")




SILVA_MERDS_rar_trt_ord=ordinate(SILVA_MERDS_rar_trt, method = "NMDS",distance = "bray")
#*** Solution reached
#Warning message:
#  In metaMDS(veganifyOTU(physeq), distance, ...) :
#  stress is (nearly) zero: you may have insufficient data
#8.121468e-05
plot_ordination(SILVA_MERDS_rar_trt,SILVA_MERDS_rar_trt_ord, color="precip",shape="life_stage", label = "block")

SILVA_MERDS_rar_trt_map=sample_data(SILVA_MERDS_rar_trt)
SILVA_MERDS_rar_trt_map$soil_root=with(SILVA_MERDS_rar_trt_map, interaction(soil_status,root_association))
SILVA_MERDS_rar_trt_dis=distance(SILVA_MERDS_rar_trt,method = "bray")

adonis(SILVA_MERDS_rar_trt_dis~SILVA_MERDS_rar_trt_map$soil_root*SILVA_MERDS_rar_trt_map$precip*SILVA_MERDS_rar_trt_map$life_stage
       +as.factor(SILVA_MERDS_rar_trt_map$block), permutations = 9999)
"Terms added sequentially (first to last)

Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
SILVA_MERDS_rar_trt_map$soil_root                                                                    2    4.3301 2.16507 17.4498 0.37696 0.0001 ***
SILVA_MERDS_rar_trt_map$precip                                                                       1    0.5459 0.54588  4.3997 0.04752 0.0017 ** 
SILVA_MERDS_rar_trt_map$life_stage                                                                   1    0.3409 0.34091  2.7477 0.02968 0.0160 *  
as.factor(SILVA_MERDS_rar_trt_map$block)                                                             3    0.5862 0.19539  1.5748 0.05103 0.0603 .  
SILVA_MERDS_rar_trt_map$soil_root:SILVA_MERDS_rar_trt_map$precip                                     2    0.4514 0.22568  1.8189 0.03929 0.0451 *  
SILVA_MERDS_rar_trt_map$soil_root:SILVA_MERDS_rar_trt_map$life_stage                                 2    0.7489 0.37447  3.0181 0.06520 0.0025 ** 
SILVA_MERDS_rar_trt_map$precip:SILVA_MERDS_rar_trt_map$life_stage                                    1    0.1749 0.17487  1.4094 0.01522 0.1648    
SILVA_MERDS_rar_trt_map$soil_root:SILVA_MERDS_rar_trt_map$precip:SILVA_MERDS_rar_trt_map$life_stage  2    0.3382 0.16910  1.3629 0.02944 0.1523    
Residuals                                                                                           32    3.9704 0.12407         0.34564           
Total                                                                                               46   11.4869                 1.00000   "


#Remove the sterile since it seems to be driving a lot of this

SILVA_MERDS_rar_trt_live=subset_samples(SILVA_MERDS_rar_trt, soil_status!="S")

SILVA_MERDS_rar_trt_live_ord=ordinate(SILVA_MERDS_rar_trt_live, method = "NMDS",distance = "bray")
#*** Solution reached
#0.1815625
plot_ordination(SILVA_MERDS_rar_trt_live,SILVA_MERDS_rar_trt_live_ord, color="precip",shape="life_stage", label = "block")


SILVA_MERDS_rar_trt_live_map=sample_data(SILVA_MERDS_rar_trt_live)
SILVA_MERDS_rar_trt_live_map$soil_root=with(SILVA_MERDS_rar_trt_live_map, interaction(soil_status,root_association))
SILVA_MERDS_rar_trt_live_dis=distance(SILVA_MERDS_rar_trt_live,method = "bray")

adonis(SILVA_MERDS_rar_trt_live_dis~SILVA_MERDS_rar_trt_live_map$soil_root*SILVA_MERDS_rar_trt_live_map$precip*SILVA_MERDS_rar_trt_live_map$life_stage
       +as.factor(SILVA_MERDS_rar_trt_live_map$block), permutations = 9999)


"Terms added sequentially (first to last)

Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
SILVA_MERDS_rar_trt_live_map$soil_root                                                                              1    0.2775 0.27753  2.4486 0.06084 0.0024 ** 
SILVA_MERDS_rar_trt_live_map$precip                                                                                 1    0.4882 0.48819  4.3072 0.10701 0.0001 ***
SILVA_MERDS_rar_trt_live_map$life_stage                                                                             1    0.1505 0.15055  1.3282 0.03300 0.1419    
as.factor(SILVA_MERDS_rar_trt_live_map$block)                                                                       3    0.6329 0.21098  1.8614 0.13874 0.0016 ** 
SILVA_MERDS_rar_trt_live_map$soil_root:SILVA_MERDS_rar_trt_live_map$precip                                          1    0.1334 0.13340  1.1770 0.02924 0.2427    
SILVA_MERDS_rar_trt_live_map$soil_root:SILVA_MERDS_rar_trt_live_map$life_stage                                      1    0.2150 0.21495  1.8965 0.04712 0.0187 *  
SILVA_MERDS_rar_trt_live_map$precip:SILVA_MERDS_rar_trt_live_map$life_stage                                         1    0.1049 0.10486  0.9251 0.02299 0.5141    
SILVA_MERDS_rar_trt_live_map$soil_root:SILVA_MERDS_rar_trt_live_map$precip:SILVA_MERDS_rar_trt_live_map$life_stage  1    0.1793 0.17929  1.5818 0.03930 0.0605 .  
Residuals                                                                                                          21    2.3802 0.11334         0.52176           
Total                                                                                                              31    4.5619                 1.00000           
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
> "

#there is a marginally significant three way interaction with life stage 
#since our other analyses are split by life stage I am going to split them here 

#Transplanted 

SILVA_MERDS_rar_trt_live_trans=subset_samples(SILVA_MERDS_rar_trt_live, life_stage=="G")

SILVA_MERDS_rar_trt_live_trans_ord=ordinate(SILVA_MERDS_rar_trt_live_trans, method = "NMDS",distance = "bray")
#*** No convergence -- monoMDS stopping criteria:
#20: stress ratio > sratmax
#0.1607911
plot_ordination(SILVA_MERDS_rar_trt_live_trans,SILVA_MERDS_rar_trt_live_trans_ord, color="precip",shape="root_association")+geom_point(size=4)+
  geom_label_repel(size=3,aes(label = block))+theme_bw()


SILVA_MERDS_rar_trt_live_trans_map=sample_data(SILVA_MERDS_rar_trt_live_trans)
SILVA_MERDS_rar_trt_live_trans_map$soil_root=with(SILVA_MERDS_rar_trt_live_trans_map, interaction(soil_status,root_association))
SILVA_MERDS_rar_trt_live_transe_dis=distance(SILVA_MERDS_rar_trt_live_trans,method = "bray")

adonis(SILVA_MERDS_rar_trt_live_transe_dis~SILVA_MERDS_rar_trt_live_trans_map$soil_root*SILVA_MERDS_rar_trt_live_trans_map$precip
       +as.factor(SILVA_MERDS_rar_trt_live_trans_map$block), permutations = 9999)
#SILVA_MERDS_rar_trt_live_trans_map$precip                                               1   0.26218 0.26218 2.13979 0.12422 0.0045 **
#as.factor(SILVA_MERDS_rar_trt_live_trans_map$block)                                     3   0.47584 0.15861 1.29453 0.22545 0.0808 . 


#seed

SILVA_MERDS_rar_trt_live_seed=subset_samples(SILVA_MERDS_rar_trt_live, life_stage=="S")

SILVA_MERDS_rar_trt_live_seed_ord=ordinate(SILVA_MERDS_rar_trt_live_seed, method = "NMDS",distance = "bray")
#*** Solution reached
#0.1359863 
plot_ordination(SILVA_MERDS_rar_trt_live_seed,SILVA_MERDS_rar_trt_live_seed_ord, color="precip",shape="root_association")+geom_point(size=4)+
  geom_label_repel(size=3,aes(label = block))+theme_bw()


SILVA_MERDS_rar_trt_live_seed_map=sample_data(SILVA_MERDS_rar_trt_live_seed)
SILVA_MERDS_rar_trt_live_seed_map$soil_root=with(SILVA_MERDS_rar_trt_live_seed_map, interaction(soil_status,root_association))
SILVA_MERDS_rar_trt_live_seed_dis=distance(SILVA_MERDS_rar_trt_live_seed,method = "bray")

adonis(SILVA_MERDS_rar_trt_live_seed_dis~SILVA_MERDS_rar_trt_live_seed_map$soil_root*SILVA_MERDS_rar_trt_live_seed_map$precip
       +as.factor(SILVA_MERDS_rar_trt_live_seed_map$block), permutations = 9999)
#SILVA_MERDS_rar_trt_live_seed_map$soil_root                                           1   0.38611 0.38611  3.3937 0.16782 0.0001 ***
#SILVA_MERDS_rar_trt_live_seed_map$precip                                              1   0.33423 0.33423  2.9378 0.14527 0.0005 ***

#Sterile 
SILVA_MERDS_rar_sterile=subset_samples(SILVA_MERDS_rar, soil_status=="S")


SILVA_MERDS_rar_sterile_ord=ordinate(SILVA_MERDS_rar_sterile, method = "NMDS",distance = "bray")
#*** Solution reached
#0.1132872
plot_ordination(SILVA_MERDS_rar_sterile,SILVA_MERDS_rar_sterile_ord, color="precip",shape="life_stage")+geom_point(size=4)+
  geom_label_repel(size=3,aes(label = block))+theme_bw()




SILVA_MERDS_rar_sterile_map=sample_data(SILVA_MERDS_rar_sterile)
SILVA_MERDS_rar_sterile_map$soil_root=with(SILVA_MERDS_rar_sterile_map, interaction(soil_status,root_association))
SILVA_MERDS_rar_sterile_dis=distance(SILVA_MERDS_rar_sterile,method = "bray")

adonis(SILVA_MERDS_rar_sterile_dis~SILVA_MERDS_rar_sterile_map$precip*SILVA_MERDS_rar_sterile_map$life_stage
       +as.factor(SILVA_MERDS_rar_sterile_map$block), permutations = 9999)

"Terms added sequentially (first to last)

Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
SILVA_MERDS_rar_sterile_map$precip                                         1   0.38409 0.38409  2.8619 0.13372 0.0063 ** 
SILVA_MERDS_rar_sterile_map$life_stage                                     1   0.72564 0.72564  5.4069 0.25263 0.0001 ***
as.factor(SILVA_MERDS_rar_sterile_map$block)                               3   0.43799 0.14600  1.0878 0.15249 0.3478    
SILVA_MERDS_rar_sterile_map$precip:SILVA_MERDS_rar_sterile_map$life_stage  1   0.25094 0.25094  1.8698 0.08737 0.0601 .  
Residuals                                                                  8   1.07367 0.13421         0.37380           
Total                                                                     14   2.87234                 1.00000           
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 "

#Let look at transplants and see if there is a correlation between biomass and OTUs

SILVA_MERDS_rar_sterile_trans=subset_samples(SILVA_MERDS_rar_sterile, life_stage=="G")
sample_data(SILVA_MERDS_rar_sterile_trans)
source(system.file("extdata/lm_phyloseq.R", package = "microbiome"))
biomas_mod_otus=lm_phyloseq(SILVA_MERDS_rar_sterile_trans, "total_biomass")
#Warning messages:
#1: In transform(x, transformation) :
#  OTU table contains zeroes. Using log10(1 + x) transform.
#2: In fitFDist(var, df1 = df, covariate = covariate) :
#  More than half of residual variances are exactly zero: eBayes unreliable
"logFC    AveExpr        t      P.Value   adj.P.Val        B
OTU259   1.4906662 0.27848584 9.925403 2.800264e-05 0.003587517 2.036929
OTU2653  0.9784821 0.18279975 9.925295 2.800460e-05 0.003587517 2.036849
OTU1215  0.5654503 0.10563726 9.924914 2.801145e-05 0.003587517 2.036572
OTU13103 0.5654503 0.10563726 9.924914 2.801145e-05 0.003587517 2.036572
OTU1866  0.5654503 0.10563726 9.924914 2.801145e-05 0.003587517 2.036572
OTU28945 0.5654503 0.10563726 9.924914 2.801145e-05 0.003587517 2.036572
OTU18008 0.5654503 0.10563726 9.924914 2.801145e-05 0.003587517 2.036572
OTU10541 0.5206565 0.09726891 9.924812 2.801330e-05 0.003587517 2.036497
OTU8891  0.5206565 0.09726891 9.924812 2.801330e-05 0.003587517 2.036497
OTU625   0.4676768 0.08737125 9.924651 2.801620e-05 0.003587517 2.036380"


get_sample(SILVA_MERDS_rar_sterile_trans,"OTU259")
get_sample(SILVA_MERDS_rar_sterile_trans,"OTU2653")
get_sample(SILVA_MERDS_rar_sterile_trans,"OTU1215")
get_sample(SILVA_MERDS_rar_sterile_trans,"OTU13103")
get_sample(SILVA_MERDS_rar_sterile_trans,"OTU1866")
get_sample(SILVA_MERDS_rar_sterile_trans,"OTU28945")
get_sample(SILVA_MERDS_rar_sterile_trans,"OTU18008")
get_sample(SILVA_MERDS_rar_sterile_trans,"OTU10541")
get_sample(SILVA_MERDS_rar_sterile_trans,"OTU8891")
get_sample(SILVA_MERDS_rar_sterile_trans,"OTU625")



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


bact_obs_rich_full_model= lm(log(Observed)~soil_root_stage+as.factor(block), data= SILVA_MERDS_rar.divfil)
qqPlot(resid(bact_obs_rich_full_model))
hist(resid(bact_obs_rich_full_model))

Anova(bact_obs_rich_full_model, type=3)
#soil_root_stage    43.36  7   174.6772 <2e-16 ***

emmeans(bact_obs_rich_full_model, pairwise~soil_root_stage)
"> emmeans(bact_obs_rich_full_model, pairwise~soil_root_stage)
$emmeans
soil_root_stage emmean     SE df lower.CL upper.CL
L.B.G             7.18 0.0666 44     7.04     7.31
S.B.G             5.37 0.0666 44     5.24     5.50
L.R.G             7.14 0.0666 44     7.00     7.27
L.B.S             6.92 0.0666 44     6.78     7.05
S.B.S             5.13 0.0722 44     4.99     5.28
L.R.S             7.21 0.0672 44     7.08     7.35
L.B.Start         7.63 0.0942 44     7.44     7.82
L.R.Start         7.64 0.0942 44     7.45     7.83

Results are averaged over the levels of: block 
Results are given on the log (not the response) scale. 
Confidence level used: 0.95 

$contrasts
contrast              estimate     SE df t.ratio p.value
L.B.G - S.B.G           1.8081 0.0942 44  19.203 <.0001 
L.B.G - L.R.G           0.0405 0.0942 44   0.430 0.9999 
L.B.G - L.B.S           0.2618 0.0942 44   2.780 0.1265 
L.B.G - S.B.S           2.0439 0.0982 44  20.807 <.0001 
L.B.G - L.R.S          -0.0343 0.0946 44  -0.362 1.0000 
L.B.G - L.B.Start      -0.4541 0.1153 44  -3.938 0.0065 
L.B.G - L.R.Start      -0.4584 0.1153 44  -3.975 0.0058 
S.B.G - L.R.G          -1.7675 0.0942 44 -18.772 <.0001 
S.B.G - L.B.S          -1.5463 0.0942 44 -16.422 <.0001 
S.B.G - S.B.S           0.2358 0.0982 44   2.401 0.2665 
S.B.G - L.R.S          -1.8423 0.0946 44 -19.476 <.0001 
S.B.G - L.B.Start      -2.2622 0.1153 44 -19.617 <.0001 
S.B.G - L.R.Start      -2.2665 0.1153 44 -19.654 <.0001 
L.R.G - L.B.S           0.2212 0.0942 44   2.350 0.2911 
L.R.G - S.B.S           2.0034 0.0982 44  20.394 <.0001 
L.R.G - L.R.S          -0.0748 0.0946 44  -0.791 0.9928 
L.R.G - L.B.Start      -0.4946 0.1153 44  -4.289 0.0023 
L.R.G - L.R.Start      -0.4989 0.1153 44  -4.327 0.0020 
L.B.S - S.B.S           1.7821 0.0982 44  18.142 <.0001 
L.B.S - L.R.S          -0.2960 0.0946 44  -3.129 0.0567 
L.B.S - L.B.Start      -0.7159 0.1153 44  -6.208 <.0001 
L.B.S - L.R.Start      -0.7202 0.1153 44  -6.245 <.0001 
S.B.S - L.R.S          -2.0781 0.0996 44 -20.863 <.0001 
S.B.S - L.B.Start      -2.4980 0.1187 44 -21.050 <.0001 
S.B.S - L.R.Start      -2.5023 0.1187 44 -21.086 <.0001 
L.R.S - L.B.Start      -0.4199 0.1157 44  -3.630 0.0155 
L.R.S - L.R.Start      -0.4242 0.1157 44  -3.667 0.0139 
L.B.Start - L.R.Start  -0.0043 0.1332 44  -0.032 1.0000 

Results are averaged over the levels of: block 
Results are given on the log (not the response) scale. 
P value adjustment: tukey method for comparing a family of 8 estimates "


ggplot(SILVA_MERDS_rar.divfil, aes(x=soil_status, y=InvSimpson))+geom_boxplot(aes(color=interaction(root_association,precip,life_stage)))+
  scale_y_continuous(name="Inv Simpson")+theme_bw()


#let's just look at the treatments
SILVA_MERDS_rar.divfil_trt=subset(SILVA_MERDS_rar.divfil, life_stage!="Start")
nrow(SILVA_MERDS_rar.divfil_trt)
#47
bact_obs_rich_model= lm(log(Observed)~life_stage*soil_root*precip+as.factor(block), data= SILVA_MERDS_rar.divfil_trt)
qqPlot(resid(bact_obs_rich_model))
hist(resid(bact_obs_rich_model))

Anova(bact_obs_rich_model, type=3)
#life_stage                     0.25  1     7.0966 0.01199 *  
#soil_root                     34.17  2   492.6814 < 2e-16 ***
#precip                         0.26  1     7.4928 0.01003 *  
#life_stage:soil_root           0.28  2     3.9931 0.02830 *  


emmeans(bact_obs_rich_model, pairwise~soil_root|life_stage)
#$contrasts
#life_stage = G:
#  contrast  estimate     SE df t.ratio p.value
#L.B - S.B   1.8081 0.0931 32  19.419 <.0001 
#L.B - L.R   0.0405 0.0931 32   0.435 0.9012 
#S.B - L.R  -1.7675 0.0931 32 -18.983 <.0001 

#life_stage = S:
#  contrast  estimate     SE df t.ratio p.value
#L.B - S.B   1.7934 0.0982 32  18.262 <.0001 
#L.B - L.R  -0.2940 0.0936 32  -3.139 0.0099 
#S.B - L.R  -2.0874 0.1001 32 -20.856 <.0001 




SILVA_MERDS_rar.divfil_trt_stage_soil_g=SILVA_MERDS_rar.divfil_trt %>% group_by(soil_root,life_stage)
obs_rich_stage=summarise_at(SILVA_MERDS_rar.divfil_trt_stage_soil_g, 
                            "Observed", funs(n(),mean,sd,se=sd(.)/sqrt(n())))

ggplot(obs_rich_stage, aes(x=life_stage,y=mean,ymin = mean-se, ymax= mean+se,fill=soil_root))+
  geom_bar(position="dodge",stat="identity", color="black",aes(fill=soil_root))+
  geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+
  scale_fill_manual(values = c("gray", "white", "black"),labels=c("Bulk","Sterile","Rhizo"))+
  scale_x_discrete(labels=c("Transplant","Seed"))+ylab("Bacterial richness")+
  geom_text(aes(y=mean+se+20, label=n),position=position_dodge(width=0.9))+theme_bw()

ggplot(SILVA_MERDS_rar.divfil_trt, aes(x=precip, y=Observed))+geom_boxplot(aes(color=precip))+theme_bw()




#Simpson
bact_inv_simp_model= lm(sqrt(InvSimpson)~life_stage*soil_root*precip+as.factor(block), data= SILVA_MERDS_rar.divfil_trt)
qqPlot(resid(bact_inv_simp_model))
hist(resid(bact_inv_simp_model))

Anova(bact_inv_simp_model, type=3,singular.ok =T)
#soil_root                    210.09  2  13.9318 4.444e-05 ***
#life_stage:soil_root          52.10  2   3.4549    0.0438 *   

emmeans(bact_inv_simp_model, pairwise~soil_root|life_stage)
#$contrasts
#life_stage = G:
#  contrast  estimate   SE df t.ratio p.value
#L.B - S.B    3.889 1.37 32  2.833  0.0210 
#L.B - L.R   -0.285 1.37 32 -0.207  0.9766 
#S.B - L.R   -4.174 1.37 32 -3.040  0.0127 

#life_stage = S:
#  contrast  estimate   SE df t.ratio p.value
#L.B - S.B    1.052 1.45 32  0.726  0.7498 
#L.B - L.R   -5.389 1.38 32 -3.903  0.0013 
#S.B - L.R   -6.441 1.48 32 -4.364  0.0004

SILVA_MERDS_rar.divfil_trt_stage_soil_g=SILVA_MERDS_rar.divfil_trt %>% group_by(soil_root,life_stage)
inv_simp_stage=summarise_at(SILVA_MERDS_rar.divfil_trt_stage_soil_g, 
                            "InvSimpson", funs(n(),mean,sd,se=sd(.)/sqrt(n())))

ggplot(inv_simp_stage, aes(x=life_stage,y=mean,ymin = mean-se, ymax= mean+se,fill=soil_root))+
  geom_bar(position="dodge",stat="identity", color="black",aes(fill=soil_root))+
  geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+
  scale_fill_manual(values = c("gray", "white", "black"),labels=c("Bulk","Sterile","Rhizo"))+
  scale_x_discrete(labels=c("Transplant","Seed"))+ylab("Inv Simpson bacteria")+
  geom_text(aes(y=mean+se+5, label=n),position=position_dodge(width=0.9))+theme_bw()


#Transplants

#let's just look at the treatments
SILVA_MERDS_rar.divfil_trt_trans=subset(SILVA_MERDS_rar.divfil_trt, life_stage=="G")
nrow(SILVA_MERDS_rar.divfil_trt_trans)
#24
bact_obs_rich_model_trans= lm(log(Observed)~soil_root*precip+as.factor(block), data= SILVA_MERDS_rar.divfil_trt_trans)
qqPlot(resid(bact_obs_rich_model_trans))
hist(resid(bact_obs_rich_model_trans))

Anova(bact_obs_rich_model_trans, type=3)
#soil_root          17.05  2   177.4953 3.625e-11 ***
#precip              0.28  1     5.8152   0.02916 *   


emmeans(bact_obs_rich_model_trans, pairwise~soil_root|precip)
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




#Simpson
bact_inv_simp_model_trans= lm((InvSimpson)^(-.5)~soil_root*precip+as.factor(block), data= SILVA_MERDS_rar.divfil_trt_trans)
qqPlot(resid(bact_inv_simp_model_trans))
hist(resid(bact_inv_simp_model_trans))
boxCox(bact_inv_simp_model_trans)
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
   geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+ylim(c(0,115))+
   scale_fill_manual(values = c( "white","light gray", "dark grey"),labels=c("Sterile","Bulk","Rhizo"), limits=treatment_order, name=NULL)+
   scale_x_discrete(labels=c("Ambient","Drought"), name=NULL)+ylab("inv Simpson bacteria")+
   geom_text(aes(y=5, label=n),size=10,position=position_dodge(width=0.9))+theme_bw()+
   theme(axis.title = element_text(size = 23), axis.text = element_text(size = 23),
         legend.position = c(0.85,.9), legend.text=element_text(size=20),
         legend.background = element_rect(size=0.5,linetype="solid",colour ="black")))



#####Are shared OTUs dominant####
SILVA_MERDS_rar_trt_live_trans_shared=core(SILVA_MERDS_rar_trt_live_trans, detection = 0,prevalence = 0.9999999999)
SILVA_MERDS_rar_trt_live_trans_shared_sampl_sum=sample_sums(SILVA_MERDS_rar_trt_live_trans_shared)
ntaxa(SILVA_MERDS_rar_trt_live_trans_shared)
#128
SILVA_MERDS_rar_trt_live_trans_map=sample_data(SILVA_MERDS_rar_trt_live_trans)
SILVA_MERDS_rar_trt_live_trans_map$soil_root=with(SILVA_MERDS_rar_trt_live_trans_map, interaction(soil_status,root_association))

shared_SILVA_MERDS_rar_trt_trans_sampl_sum=merge(SILVA_MERDS_rar_trt_live_trans_map,SILVA_MERDS_rar_trt_live_trans_shared_sampl_sum, by="row.names")


shared_SILVA_MERDS_rar_trt_trans_sampl_sum_precip_soil_g=shared_SILVA_MERDS_rar_trt_trans_sampl_sum %>% group_by(soil_root,precip)
shared_read_abund_precip_trans=summarise_at(shared_SILVA_MERDS_rar_trt_trans_sampl_sum_precip_soil_g, 
                                            "y", list(~n(),~mean,~sd,se=~sd(.)/sqrt(n())))

ggplot(shared_read_abund_precip_trans, aes(x=precip,y=mean,ymin = mean-se, ymax= mean+se,fill=soil_root))+
  geom_bar(position="dodge",stat="identity", color="black",aes(fill=soil_root))+
  geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+
  scale_fill_manual(values = c("light gray", "dark grey"),labels=c("Bulk","Rhizo"))+
  scale_x_discrete(labels=c("Ambient","Drought"))+ylab("Shared Taxa Read Abun")+
  geom_text(aes(y=mean+se+200, label=n),position=position_dodge(width=0.9))+theme_bw()


#Diversity of the Shared taxa
alpha_meas = c("Observed", "Shannon", "InvSimpson")
SILVA_MERDS_rar_trt_live_trans_map=sample_data(SILVA_MERDS_rar_trt_live_trans)
SILVA_MERDS_rar_trt_live_trans_map$soil_root=with(SILVA_MERDS_rar_trt_live_trans_map, interaction(soil_status,root_association))

SILVA_MERDS_rar_trt_live_trans_shared.divfil=estimate_richness(SILVA_MERDS_rar_trt_live_trans_shared,measures=alpha_meas)

SILVA_MERDS_rar_trt_live_trans_shared.divfil=merge(SILVA_MERDS_rar_trt_live_trans_shared.divfil, SILVA_MERDS_rar_trt_live_trans_map, by ="row.names")
#bact.soilE.t.divfil=mutate(bact.soilE.t.divfil, pielou=Shannon*(1/log(Observed)))
head(SILVA_MERDS_rar_trt_live_trans_shared.divfil)
row.names(SILVA_MERDS_rar_trt_live_trans_shared.divfil)=SILVA_MERDS_rar_trt_live_trans_shared.divfil$Row.names
SILVA_MERDS_rar_trt_live_trans_shared.divfil$Row.names=NULL



SILVA_MERDS_rar_trt_live_trans_shared.divfil_g=SILVA_MERDS_rar_trt_live_trans_shared.divfil %>% group_by(soil_root,precip)
shared_inv_simp_precip_trans=summarise_at(SILVA_MERDS_rar_trt_live_trans_shared.divfil_g, 
                                          "InvSimpson", list(~n(),~mean,~sd,se=~sd(.)/sqrt(n())))

ggplot(shared_inv_simp_precip_trans, aes(x=precip,y=mean,ymin = mean-se, ymax= mean+se,fill=soil_root))+
  geom_bar(position="dodge",stat="identity", color="black",aes(fill=soil_root))+
  geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+
  scale_fill_manual(values = c("light gray", "dark grey"),labels=c("Bulk","Rhizo"))+
  scale_x_discrete(labels=c("Ambient","Drought"))+ylab("Shared inv Simpson")+
  geom_text(aes(y=mean+se+1, label=n),position=position_dodge(width=0.9))+theme_bw()


#Found in at least one sample per treatment
SILVA_MERDS_rar_trt_live_trans_map=sample_data(SILVA_MERDS_rar_trt_live_trans)
SILVA_MERDS_rar_trt_live_trans_map$soil_root=with(SILVA_MERDS_rar_trt_live_trans_map, interaction(soil_status,root_association))
summary(SILVA_MERDS_rar_trt_live_trans_map)

SILVA_MERDS_rar_trt_live_trans_b_a=subset_samples(SILVA_MERDS_rar_trt_live_trans, root_association=="B"&precip=="A")
SILVA_MERDS_rar_trt_live_trans_b_a=prune_taxa(taxa_sums(SILVA_MERDS_rar_trt_live_trans_b_a) > 0, SILVA_MERDS_rar_trt_live_trans_b_a)
ntaxa(SILVA_MERDS_rar_trt_live_trans_b_a)
#3108
tran_bulk.amb.n<-taxa_names(SILVA_MERDS_rar_trt_live_trans_b_a)
length(tran_bulk.amb.n)
#3108

SILVA_MERDS_rar_trt_live_trans_r_a=subset_samples(SILVA_MERDS_rar_trt_live_trans, root_association=="R"&precip=="A")
SILVA_MERDS_rar_trt_live_trans_r_a=prune_taxa(taxa_sums(SILVA_MERDS_rar_trt_live_trans_r_a) > 0, SILVA_MERDS_rar_trt_live_trans_r_a)
ntaxa(SILVA_MERDS_rar_trt_live_trans_r_a)
#3463
trans_rhizo.amb.n<-taxa_names(SILVA_MERDS_rar_trt_live_trans_r_a)
length(trans_rhizo.amb.n)
#3463

trans.amb_b_in_r<-tran_bulk.amb.n[tran_bulk.amb.n %in% trans_rhizo.amb.n]
length(trans.amb_b_in_r) ## How many OTUs
#2032


SILVA_MERDS_rar_trt_live_trans_b_d=subset_samples(SILVA_MERDS_rar_trt_live_trans, root_association=="B"&precip=="D")
SILVA_MERDS_rar_trt_live_trans_b_d=prune_taxa(taxa_sums(SILVA_MERDS_rar_trt_live_trans_b_d) > 0, SILVA_MERDS_rar_trt_live_trans_b_d)
ntaxa(SILVA_MERDS_rar_trt_live_trans_b_d)
#2780
trans_bulk.drought.n<-taxa_names(SILVA_MERDS_rar_trt_live_trans_b_d)
length(trans_bulk.drought.n)
#2780

SILVA_MERDS_rar_trt_live_trans_r_d=subset_samples(SILVA_MERDS_rar_trt_live_trans, root_association=="R"&precip=="D")
SILVA_MERDS_rar_trt_live_trans_r_d=prune_taxa(taxa_sums(SILVA_MERDS_rar_trt_live_trans_r_d) > 0, SILVA_MERDS_rar_trt_live_trans_r_d)
ntaxa(SILVA_MERDS_rar_trt_live_trans_r_d)
#2405
trans_rhizo.drought.n<-taxa_names(SILVA_MERDS_rar_trt_live_trans_r_d)
length(trans_rhizo.drought.n)
#2405


trans.drought_b_in_r<-trans_bulk.drought.n[trans_bulk.drought.n %in% trans_rhizo.drought.n]
length(trans.drought_b_in_r) ## How many OTUs
#1593

trans.drought_in_amb<-trans.amb_b_in_r[trans.amb_b_in_r %in% trans.drought_b_in_r]
length(trans.drought_in_amb) ## How many OTUs
#1163


SILVA_MERDS_rar_trt_live_trans_sh<-prune_taxa(trans.drought_in_amb,SILVA_MERDS_rar_trt_live_trans)
SILVA_MERDS_rar_trt_live_trans_sh_sampl_sum=sample_sums(SILVA_MERDS_rar_trt_live_trans_sh)
ntaxa(SILVA_MERDS_rar_trt_live_trans_sh)
#1163
SILVA_MERDS_rar_trt_live_trans_map=sample_data(SILVA_MERDS_rar_trt_live_trans)
SILVA_MERDS_rar_trt_live_trans_map$soil_root=with(SILVA_MERDS_rar_trt_live_trans_map, interaction(soil_status,root_association))

sh_SILVA_MERDS_rar_trt_trans_sampl_sum=merge(SILVA_MERDS_rar_trt_live_trans_map,SILVA_MERDS_rar_trt_live_trans_sh_sampl_sum, by="row.names")


sh_SILVA_MERDS_rar_trt_trans_sampl_sum_precip_soil_g=sh_SILVA_MERDS_rar_trt_trans_sampl_sum %>% group_by(soil_root,precip)
sh_read_abund_precip_trans=summarise_at(sh_SILVA_MERDS_rar_trt_trans_sampl_sum_precip_soil_g, 
                                        "y", list(~n(),~mean,~sd,se=~sd(.)/sqrt(n())))

ggplot(sh_read_abund_precip_trans, aes(x=precip,y=mean,ymin = mean-se, ymax= mean+se,fill=soil_root))+
  geom_bar(position="dodge",stat="identity", color="black",aes(fill=soil_root))+
  geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+
  scale_fill_manual(values = c("light gray", "dark grey"),labels=c("Bulk","Rhizo"))+
  scale_x_discrete(labels=c("Ambient","Drought"))+ylab("Shared Taxa Read Abun")+
  geom_text(aes(y=mean+se+500, label=n),position=position_dodge(width=0.9))+theme_bw()


#Diversity of the Shared taxa
alpha_meas = c("Observed", "Shannon", "InvSimpson")
SILVA_MERDS_rar_trt_live_trans_map=sample_data(SILVA_MERDS_rar_trt_live_trans)
SILVA_MERDS_rar_trt_live_trans_map$soil_root=with(SILVA_MERDS_rar_trt_live_trans_map, interaction(soil_status,root_association))

SILVA_MERDS_rar_trt_live_trans_sh.divfil=estimate_richness(SILVA_MERDS_rar_trt_live_trans_sh,measures=alpha_meas)

SILVA_MERDS_rar_trt_live_trans_sh.divfil=merge(SILVA_MERDS_rar_trt_live_trans_sh.divfil, SILVA_MERDS_rar_trt_live_trans_map, by ="row.names")
#bact.soilE.t.divfil=mutate(bact.soilE.t.divfil, pielou=Shannon*(1/log(Observed)))
head(SILVA_MERDS_rar_trt_live_trans_sh.divfil)
row.names(SILVA_MERDS_rar_trt_live_trans_sh.divfil)=SILVA_MERDS_rar_trt_live_trans_sh.divfil$Row.names
SILVA_MERDS_rar_trt_live_trans_sh.divfil$Row.names=NULL



SILVA_MERDS_rar_trt_live_trans_sh.divfil_g=SILVA_MERDS_rar_trt_live_trans_sh.divfil %>% group_by(soil_root,precip)
sh_inv_simp_precip_trans=summarise_at(SILVA_MERDS_rar_trt_live_trans_sh.divfil_g, 
                                      "InvSimpson", list(~n(),~mean,~sd,se=~sd(.)/sqrt(n())))

ggplot(sh_inv_simp_precip_trans, aes(x=precip,y=mean,ymin = mean-se, ymax= mean+se,fill=soil_root))+
  geom_bar(position="dodge",stat="identity", color="black",aes(fill=soil_root))+
  geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+
  scale_fill_manual(values = c("light gray", "dark grey"),labels=c("Bulk","Rhizo"))+
  scale_x_discrete(labels=c("Ambient","Drought"))+ylab("Shared inv Simpson")+
  geom_text(aes(y=mean+se+5, label=n),position=position_dodge(width=0.9))+theme_bw()


#Seeds

#let's just look at the treatments
SILVA_MERDS_rar.divfil_trt_seed=subset(SILVA_MERDS_rar.divfil_trt, life_stage=="S")
nrow(SILVA_MERDS_rar.divfil_trt_seed)
#23
bact_obs_rich_model_seed= lm(log(Observed)~soil_root*precip+as.factor(block), data= SILVA_MERDS_rar.divfil_trt_seed)
qqPlot(resid(bact_obs_rich_model_seed))
hist(resid(bact_obs_rich_model_seed))
boxCox(bact_obs_rich_model_seed)
Anova(bact_obs_rich_model_seed, type=3)
#soil_root         16.28  2   318.5927 2.123e-12 ***



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




#Simpson
bact_inv_simp_model_seed= lm(log(InvSimpson)~soil_root*precip+as.factor(block), data= SILVA_MERDS_rar.divfil_trt_seed)
qqPlot(resid(bact_inv_simp_model_seed))
hist(resid(bact_inv_simp_model_seed))
boxCox(bact_inv_simp_model_seed)
Anova(bact_inv_simp_model_seed, type=3)
#soil_root         17.314  2  11.9446 0.0009404 ***

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

SILVA_MERDS_rar.divfil_trt_seed_precip_soil_g=SILVA_MERDS_rar.divfil_trt_seed %>% group_by(soil_root,precip)
inv_simp_precip_seed=summarise_at(SILVA_MERDS_rar.divfil_trt_seed_precip_soil_g, 
                                  "InvSimpson", funs(n(),mean,sd,se=sd(.)/sqrt(n())))

ggplot(inv_simp_precip_seed, aes(x=precip,y=mean,ymin = mean-se, ymax= mean+se,fill=soil_root))+
  geom_bar(position="dodge",stat="identity", color="black",aes(fill=soil_root))+
  geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+
  scale_fill_manual(values = c("gray", "white", "black"),labels=c("Bulk","Sterile","Rhizo"))+
  scale_x_discrete(labels=c("Ambient","Drought"))+ylab("Inv Simpson bacteria")+
  geom_text(aes(y=mean+se+5, label=n),position=position_dodge(width=0.9))+theme_bw()


#####Are shared OTUs dominant####
SILVA_MERDS_rar_trt_live_seed_shared=core(SILVA_MERDS_rar_trt_live_seed, detection = 0,prevalence = 0.9999999999)
SILVA_MERDS_rar_trt_live_seed_shared_sampl_sum=sample_sums(SILVA_MERDS_rar_trt_live_seed_shared)
ntaxa(SILVA_MERDS_rar_trt_live_seed_shared)
#107
SILVA_MERDS_rar_trt_live_seed_map=sample_data(SILVA_MERDS_rar_trt_live_seed)
SILVA_MERDS_rar_trt_live_seed_map$soil_root=with(SILVA_MERDS_rar_trt_live_seed_map, interaction(soil_status,root_association))

shared_SILVA_MERDS_rar_trt_seed_sampl_sum=merge(SILVA_MERDS_rar_trt_live_seed_map,SILVA_MERDS_rar_trt_live_seed_shared_sampl_sum, by="row.names")


shared_SILVA_MERDS_rar_trt_seed_sampl_sum_precip_soil_g=shared_SILVA_MERDS_rar_trt_seed_sampl_sum %>% group_by(soil_root,precip)
shared_read_abund_precip_seed=summarise_at(shared_SILVA_MERDS_rar_trt_seed_sampl_sum_precip_soil_g, 
                                           "y", list(~n(),~mean,~sd,se=~sd(.)/sqrt(n())))

ggplot(shared_read_abund_precip_seed, aes(x=precip,y=mean,ymin = mean-se, ymax= mean+se,fill=soil_root))+
  geom_bar(position="dodge",stat="identity", color="black",aes(fill=soil_root))+
  geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+
  scale_fill_manual(values = c("light gray", "dark grey"),labels=c("Bulk","Rhizo"))+
  scale_x_discrete(labels=c("Ambient","Drought"))+ylab("Shared Taxa Read Abun")+
  geom_text(aes(y=mean+se+100, label=n),position=position_dodge(width=0.9))+theme_bw()


#Diversity of the Shared taxa
alpha_meas = c("Observed", "Shannon", "InvSimpson")
SILVA_MERDS_rar_trt_live_seed_map=sample_data(SILVA_MERDS_rar_trt_live_seed)
SILVA_MERDS_rar_trt_live_seed_map$soil_root=with(SILVA_MERDS_rar_trt_live_seed_map, interaction(soil_status,root_association))

SILVA_MERDS_rar_trt_live_seed_shared.divfil=estimate_richness(SILVA_MERDS_rar_trt_live_seed_shared,measures=alpha_meas)

SILVA_MERDS_rar_trt_live_seed_shared.divfil=merge(SILVA_MERDS_rar_trt_live_seed_shared.divfil, SILVA_MERDS_rar_trt_live_seed_map, by ="row.names")
#bact.soilE.t.divfil=mutate(bact.soilE.t.divfil, pielou=Shannon*(1/log(Observed)))
head(SILVA_MERDS_rar_trt_live_seed_shared.divfil)
row.names(SILVA_MERDS_rar_trt_live_seed_shared.divfil)=SILVA_MERDS_rar_trt_live_seed_shared.divfil$Row.names
SILVA_MERDS_rar_trt_live_seed_shared.divfil$Row.names=NULL



SILVA_MERDS_rar_trt_live_seed_shared.divfil_g=SILVA_MERDS_rar_trt_live_seed_shared.divfil %>% group_by(soil_root,precip)
shared_inv_simp_precip_seed=summarise_at(SILVA_MERDS_rar_trt_live_seed_shared.divfil_g, 
                                         "InvSimpson", list(~n(),~mean,~sd,se=~sd(.)/sqrt(n())))

ggplot(shared_inv_simp_precip_seed, aes(x=precip,y=mean,ymin = mean-se, ymax= mean+se,fill=soil_root))+
  geom_bar(position="dodge",stat="identity", color="black",aes(fill=soil_root))+
  geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+
  scale_fill_manual(values = c("light gray", "dark grey"),labels=c("Bulk","Rhizo"))+
  scale_x_discrete(labels=c("Ambient","Drought"))+ylab("Shared inv Simpson")+
  geom_text(aes(y=mean+se+1, label=n),position=position_dodge(width=0.9))+theme_bw()

#Found in at least one sample per treatment
SILVA_MERDS_rar_trt_live_seed_map=sample_data(SILVA_MERDS_rar_trt_live_seed)
SILVA_MERDS_rar_trt_live_seed_map$soil_root=with(SILVA_MERDS_rar_trt_live_seed_map, interaction(soil_status,root_association))
summary(SILVA_MERDS_rar_trt_live_seed_map)

SILVA_MERDS_rar_trt_live_seed_b_a=subset_samples(SILVA_MERDS_rar_trt_live_seed, root_association=="B"&precip=="A")
SILVA_MERDS_rar_trt_live_seed_b_a=prune_taxa(taxa_sums(SILVA_MERDS_rar_trt_live_seed_b_a) > 0, SILVA_MERDS_rar_trt_live_seed_b_a)
ntaxa(SILVA_MERDS_rar_trt_live_seed_b_a)
#2459
bulk.amb.n<-taxa_names(SILVA_MERDS_rar_trt_live_seed_b_a)
length(bulk.amb.n)
#2459

SILVA_MERDS_rar_trt_live_seed_r_a=subset_samples(SILVA_MERDS_rar_trt_live_seed, root_association=="R"&precip=="A")
SILVA_MERDS_rar_trt_live_seed_r_a=prune_taxa(taxa_sums(SILVA_MERDS_rar_trt_live_seed_r_a) > 0, SILVA_MERDS_rar_trt_live_seed_r_a)
ntaxa(SILVA_MERDS_rar_trt_live_seed_r_a)
#3143
rhizo.amb.n<-taxa_names(SILVA_MERDS_rar_trt_live_seed_r_a)
length(rhizo.amb.n)
#3143

bact.amb_b_in_r<-bulk.amb.n[bulk.amb.n %in% rhizo.amb.n]
length(bact.amb_b_in_r) ## How many OTUs
#1653


SILVA_MERDS_rar_trt_live_seed_b_d=subset_samples(SILVA_MERDS_rar_trt_live_seed, root_association=="B"&precip=="D")
SILVA_MERDS_rar_trt_live_seed_b_d=prune_taxa(taxa_sums(SILVA_MERDS_rar_trt_live_seed_b_d) > 0, SILVA_MERDS_rar_trt_live_seed_b_d)
ntaxa(SILVA_MERDS_rar_trt_live_seed_b_d)
#2294
bulk.drought.n<-taxa_names(SILVA_MERDS_rar_trt_live_seed_b_d)
length(bulk.drought.n)
#2294

SILVA_MERDS_rar_trt_live_seed_r_d=subset_samples(SILVA_MERDS_rar_trt_live_seed, root_association=="R"&precip=="D")
SILVA_MERDS_rar_trt_live_seed_r_d=prune_taxa(taxa_sums(SILVA_MERDS_rar_trt_live_seed_r_d) > 0, SILVA_MERDS_rar_trt_live_seed_r_d)
ntaxa(SILVA_MERDS_rar_trt_live_seed_r_d)
#2992
rhizo.drought.n<-taxa_names(SILVA_MERDS_rar_trt_live_seed_r_d)
length(rhizo.drought.n)
#2992


bact.drought_b_in_r<-bulk.drought.n[bulk.drought.n %in% rhizo.drought.n]
length(bact.drought_b_in_r) ## How many OTUs
#1560

seed.drought_in_amb<-bact.amb_b_in_r[bact.amb_b_in_r %in% bact.drought_b_in_r]
length(seed.drought_in_amb) ## How many OTUs
#1025


SILVA_MERDS_rar_trt_live_seed_sh<-prune_taxa(seed.drought_in_amb,SILVA_MERDS_rar_trt_live_seed)
SILVA_MERDS_rar_trt_live_seed_sh_sampl_sum=sample_sums(SILVA_MERDS_rar_trt_live_seed_sh)
ntaxa(SILVA_MERDS_rar_trt_live_seed_sh)
#1025
SILVA_MERDS_rar_trt_live_seed_map=sample_data(SILVA_MERDS_rar_trt_live_seed)
SILVA_MERDS_rar_trt_live_seed_map$soil_root=with(SILVA_MERDS_rar_trt_live_seed_map, interaction(soil_status,root_association))

sh_SILVA_MERDS_rar_trt_seed_sampl_sum=merge(SILVA_MERDS_rar_trt_live_seed_map,SILVA_MERDS_rar_trt_live_seed_sh_sampl_sum, by="row.names")


sh_SILVA_MERDS_rar_trt_seed_sampl_sum_precip_soil_g=sh_SILVA_MERDS_rar_trt_seed_sampl_sum %>% group_by(soil_root,precip)
sh_read_abund_precip_seed=summarise_at(sh_SILVA_MERDS_rar_trt_seed_sampl_sum_precip_soil_g, 
                                       "y", list(~n(),~mean,~sd,se=~sd(.)/sqrt(n())))

ggplot(sh_read_abund_precip_seed, aes(x=precip,y=mean,ymin = mean-se, ymax= mean+se,fill=soil_root))+
  geom_bar(position="dodge",stat="identity", color="black",aes(fill=soil_root))+
  geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+
  scale_fill_manual(values = c("light gray", "dark grey"),labels=c("Bulk","Rhizo"))+
  scale_x_discrete(labels=c("Ambient","Drought"))+ylab("Shared Taxa Read Abun")+
  geom_text(aes(y=mean+se+300, label=n),position=position_dodge(width=0.9))+theme_bw()


#Diversity of the Shared taxa
alpha_meas = c("Observed", "Shannon", "InvSimpson")
SILVA_MERDS_rar_trt_live_seed_map=sample_data(SILVA_MERDS_rar_trt_live_seed)
SILVA_MERDS_rar_trt_live_seed_map$soil_root=with(SILVA_MERDS_rar_trt_live_seed_map, interaction(soil_status,root_association))

SILVA_MERDS_rar_trt_live_seed_sh.divfil=estimate_richness(SILVA_MERDS_rar_trt_live_seed_sh,measures=alpha_meas)

SILVA_MERDS_rar_trt_live_seed_sh.divfil=merge(SILVA_MERDS_rar_trt_live_seed_sh.divfil, SILVA_MERDS_rar_trt_live_seed_map, by ="row.names")
#bact.soilE.t.divfil=mutate(bact.soilE.t.divfil, pielou=Shannon*(1/log(Observed)))
head(SILVA_MERDS_rar_trt_live_seed_sh.divfil)
row.names(SILVA_MERDS_rar_trt_live_seed_sh.divfil)=SILVA_MERDS_rar_trt_live_seed_sh.divfil$Row.names
SILVA_MERDS_rar_trt_live_seed_sh.divfil$Row.names=NULL



SILVA_MERDS_rar_trt_live_seed_sh.divfil_g=SILVA_MERDS_rar_trt_live_seed_sh.divfil %>% group_by(soil_root,precip)
sh_inv_simp_precip_seed=summarise_at(SILVA_MERDS_rar_trt_live_seed_sh.divfil_g, 
                                     "InvSimpson", list(~n(),~mean,~sd,se=~sd(.)/sqrt(n())))

ggplot(sh_inv_simp_precip_seed, aes(x=precip,y=mean,ymin = mean-se, ymax= mean+se,fill=soil_root))+
  geom_bar(position="dodge",stat="identity", color="black",aes(fill=soil_root))+
  geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+
  scale_fill_manual(values = c("light gray", "dark grey"),labels=c("Bulk","Rhizo"))+
  scale_x_discrete(labels=c("Ambient","Drought"))+ylab("Shared inv Simpson")+
  geom_text(aes(y=mean+se+5, label=n),position=position_dodge(width=0.9))+theme_bw()


#####End ZOTU based community####