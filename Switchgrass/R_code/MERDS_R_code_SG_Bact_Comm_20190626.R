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
source(system.file("extdata/lm_phyloseq.R", package = "microbiome"))
options(contrasts=c("contr.sum", "contr.poly"))
data_SG_biomass <- read.csv("D:/MERDS_2018/merds/Switchgrass/R_data/SG_TotalBiomass.csv")
SG_inorg_N<- read.csv("D:/MERDS_2018/merds/Switchgrass/R_data/SG_NO3NH4_ug_gdrysoil.csv", header = T)
summary(SG_inorg_N)
dataSG_seed_surv <- read.csv("D:/MERDS_2018/merds/Switchgrass/R_data/SG_surv_Seed_germ.csv")
summary(dataSG_seed_surv)
dataSG_seed_surv[,8:14]=NULL
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
dataSG_seed_surv_trt$pot_w_germ[dataSG_seed_surv$pot_w_germ>0]=1
dataSG_seed_surv_1=subset(dataSG_seed_surv, pot_w_germ==1)

dataSG_seed_surv_1_pot_g=group_by(dataSG_seed_surv_1, Plant_Number)

dataSG_seed_surv_first_germ=summarise_at(dataSG_seed_surv_1_pot_g, "exp_days", min)
summary(dataSG_seed_surv_first_germ)
colnames(dataSG_seed_surv_first_germ)[2]="days_to_germ"
merds_map_trt_bio_nit_germ=merge(merds_map_trt_bio_nit, dataSG_seed_surv_first_germ, by="Plant_Number", all.x = T)
#Let's add in the germination results


dataSG_seed_surv_trt$pot_w_germ=dataSG_seed_surv_trt$num_germinates
dataSG_seed_surv_trt$pot_w_germ[dataSG_seed_surv_trt$pot_w_germ>0]=1
summary(dataSG_seed_surv_trt)
dataSG_seed_surv_trt_1=subset(dataSG_seed_surv_trt, pot_w_germ==1)

dataSG_seed_surv_trt_1_pot_g=group_by(dataSG_seed_surv_trt_1, Plant_Number)

dataSG_seed_surv_trt_1_pot_first_germ=summarise_at(dataSG_seed_surv_trt_1_pot_g, "exp_days", min)


dataSG_seed_1st_germ_trt=merge(dataSG_seed_surv_trt_1_pot_first_germ, SG_trt, by="Plant_Number", all.y = T)
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
#####Begin SILVA dataset v123####
######Begin Processing#### 
#now put the new mapping file into our phyloseq obj

SILVA_MERDS_data=phyloseq(otu_table(phyl_SILVA_MERDS),tax_table(phyl_SILVA_MERDS),sample_data(merds_map_trt_bio_nit))
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

######END NON-Rarefied community Analyses####

######Rarefied community Analyses####

#Let's rarefy the data
#SILVA_MERDS_rar=rarefy_even_depth(SILVA_MERDS_trunc3, sample.size= 10000, rngseed = T)
#498OTUs were removed because they are no longer 
#present in any sample after random subsampling

#save(SILVA_MERDS_rar, file = "D:/MERDS_2018/merds/Switchgrass/R_data/SILVA_MERDS_rar_phylo_obj.RData")
load("D:/MERDS_2018/merds/Switchgrass/R_data/SILVA_MERDS_rar_phylo_obj.RData")
#SILVA_MERDS_rar=phyloseq(otu_table(SILVA_MERDS_rar),tax_table(SILVA_MERDS_rar),sample_data(dataSG_seed_1st_germ_trt))


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
#####NMDS Compound Graph####
#Let's Graph this

#Live transplant
SILVA_MERDS_rar_live_trans=subset_samples(SILVA_MERDS_rar, soil_status!="S"&life_stage=="G"|life_stage=="Start")
summary(sample_data(SILVA_MERDS_rar_live_trans))
SILVA_MERDS_rar_live_trans_ord=ordinate(SILVA_MERDS_rar_live_trans, method = "NMDS",distance = "bray")
#*** Solution reached
#0.1242186 
(live_trans_ord_p=plot_ordination(SILVA_MERDS_rar_live_trans,SILVA_MERDS_rar_live_trans_ord,shape="root_association")+geom_point(size=5, aes(shape=root_association, fill=precip))+
    theme_bw()+scale_fill_manual(values = c("white","dark grey","black"))+scale_shape_manual(values = c(21,24))+
    ggtitle(label = "Transplant")+theme(plot.title = element_text(hjust = 0.5)))
(live_trans_ord_p_pot_num=plot_ordination(SILVA_MERDS_rar_live_trans,SILVA_MERDS_rar_live_trans_ord,shape="root_association")+geom_point(size=5, aes(shape=root_association, fill=precip))+
  geom_label_repel(size=3,aes(label = Plant_Number))+theme_bw()+scale_fill_manual(values = c("white","dark grey","black"))+scale_shape_manual(values = c(21,24))+
  ggtitle(label = "Microbial Soils")+theme(plot.title = element_text(hjust = 0.5)))

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


#Live seed
SILVA_MERDS_rar_live_seed=subset_samples(SILVA_MERDS_rar, soil_status!="S"&life_stage=="S"|life_stage=="Start")
summary(sample_data(SILVA_MERDS_rar_live_seed))
SILVA_MERDS_rar_live_seed_ord=ordinate(SILVA_MERDS_rar_live_seed, method = "NMDS",distance = "bray")
#*** Solution reached
#0.1153055
(live_seed_ord_p=plot_ordination(SILVA_MERDS_rar_live_seed,SILVA_MERDS_rar_live_seed_ord,shape="root_association")+geom_point(size=5, aes(shape=root_association, fill=precip))+
    theme_bw()+scale_fill_manual(values = c("white","dark grey","black"))+scale_shape_manual(values = c(21,24))+
    ggtitle(label = "Seedling")+theme(plot.title = element_text(hjust = 0.5), axis.title.y = element_blank()))
(live_seed_ord_p_pot_num=plot_ordination(SILVA_MERDS_rar_live_seed,SILVA_MERDS_rar_live_seed_ord,shape="root_association")+geom_point(size=5, aes(shape=root_association, fill=precip))+
    geom_label_repel(size=3,aes(label = Plant_Number))+theme_bw()+scale_fill_manual(values = c("white","dark grey","black"))+scale_shape_manual(values = c(21,24))+
    ggtitle(label = "Microbial Soils")+theme(plot.title = element_text(hjust = 0.5)))

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


#####Sterile Analyses ####
#Sterile 
SILVA_MERDS_rar_sterile=subset_samples(SILVA_MERDS_rar, soil_status=="S")
ntaxa(SILVA_MERDS_rar_sterile)
#8366
SILVA_MERDS_rar_sterile=prune_taxa(taxa_sums(SILVA_MERDS_rar_sterile) > 0, SILVA_MERDS_rar_sterile)
ntaxa(SILVA_MERDS_rar_sterile)
#930
SILVA_MERDS_rar_sterile_ord=ordinate(SILVA_MERDS_rar_sterile, method = "NMDS",distance = "bray")
#*** Solution reached
#0.1223411
plot_ordination(SILVA_MERDS_rar_sterile,SILVA_MERDS_rar_sterile_ord, color="precip",shape="life_stage")+geom_point(size=4)+
  geom_label_repel(size=3,aes(label = block))+theme_bw()

(sterile_ord_p=plot_ordination(SILVA_MERDS_rar_sterile,SILVA_MERDS_rar_sterile_ord,shape="life_stage")+geom_point(size=5, aes(shape=life_stage, fill=precip))+
    theme_bw()+scale_fill_manual(values = c("white","dark grey","black"))+scale_shape_manual(values = c(22,23))+
    ggtitle(label = "Sterile")+theme(plot.title = element_text(hjust = 0.5)))


SILVA_MERDS_rar_sterile_map=sample_data(SILVA_MERDS_rar_sterile)
SILVA_MERDS_rar_sterile_map$soil_root=with(SILVA_MERDS_rar_sterile_map, interaction(soil_status,root_association))
SILVA_MERDS_rar_sterile_dis=distance(SILVA_MERDS_rar_sterile,method = "bray")

adonis(SILVA_MERDS_rar_sterile_dis~SILVA_MERDS_rar_sterile_map$precip*SILVA_MERDS_rar_sterile_map$life_stage
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
#No sterile
SILVA_MERDS_rar_trt_live_trans=subset_samples(SILVA_MERDS_rar_trt_trans, soil_root!="S.B")
nrow(sample_data(SILVA_MERDS_rar_trt_live_trans))
#16
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


#Diversity 
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
   scale_fill_manual(values = c( "white","light gray", "dark grey"),labels=c("Sterile","Bulk","Rhizo"), limits=treatment_order, name=NULL)+
   scale_x_discrete(labels=c("Ambient","Drought"), name=NULL)+ylab("Bacterial richness")+
   geom_text(aes(y=80, label=n),size=10,position=position_dodge(width=0.9))+theme_bw()+
   theme(axis.title = element_text(size = 23), axis.text = element_text(size = 23),
         legend.position = c(0.85,.9), legend.text=element_text(size=20),
         legend.background = element_rect(size=0.5,linetype="solid",colour ="black")))



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
   geom_errorbar(position=position_dodge(width=0.9), width=0.2, size=1)+ylim(c(-10,200))+
   scale_fill_manual(values = c( "white","light gray", "dark grey"),labels=c("Sterile","Bulk","Rhizo"), limits=treatment_order, name=NULL)+
   scale_x_discrete(labels=c("Ambient","Drought"), name=NULL)+ylab("inv Simpson bacteria")+
   geom_text(aes(y=-7, label=n),size=10,position=position_dodge(width=0.9))+theme_bw()+
   theme(axis.title = element_text(size = 23), axis.text = element_text(size = 23),
         legend.position = c(0.85,.9), legend.text=element_text(size=20),
         legend.background = element_rect(size=0.5,linetype="solid",colour ="black")))




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



#####Combined taxon stacked bargraph#####
ggarrange(p_bact_trans_color,p_bact_seed_color,ncol = 2, common.legend=F)
#20x7.38



#####Seed live Community analyses####
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
SILVA_MERDS_rar.divfil_trt_seed %>% group_by(soil_root) %>% summarise_at("Observed", funs(n(),mean,sd,se=sd(.)/sqrt(n())))


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
   scale_fill_manual(values = c( "white","light gray", "dark grey"),labels=c("Sterile","Bulk","Rhizo"), limits=treatment_order, name=NULL)+
   scale_x_discrete(labels=c("Ambient","Drought"), name=NULL)+ylab("Bacterial richness")+
   geom_text(aes(y=80, label=n),size=10,position=position_dodge(width=0.9))+theme_bw()+
   theme(axis.title.x = element_text(size = 23), axis.text.x = element_text(size = 23),
         axis.title.y = element_blank(), axis.text.y = element_blank(),
         legend.position = c(0.85,.9), legend.text=element_text(size=20),
         legend.background = element_rect(size=0.5,linetype="solid",colour ="black")))

#####Combined obs rich graph#####
ggarrange(trans_obs_richness_p,seed_obs_richness_p,ncol = 2,  legend = "none",widths = c(1,.90))
#15x7.38



#Simpson
bact_inv_simp_model_seed= lm(log(InvSimpson)~soil_root*precip+as.factor(block), data= SILVA_MERDS_rar.divfil_trt_seed)
qqPlot(resid(bact_inv_simp_model_seed))
hist(resid(bact_inv_simp_model_seed))
boxCox(bact_inv_simp_model_seed)
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

SILVA_MERDS_rar.divfil_trt_seed_precip_soil_g=SILVA_MERDS_rar.divfil_trt_seed %>% group_by(soil_root,precip)
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
   scale_fill_manual(values = c( "white","light gray", "dark grey"),labels=c("Sterile","Bulk","Rhizo"), limits=treatment_order, name=NULL)+
   scale_x_discrete(labels=c("Ambient","Drought"), name=NULL)+ylab("inv Simpson bacteria")+
   geom_text(aes(y=-7, label=n),size=10,position=position_dodge(width=0.9))+theme_bw()+
   theme(axis.title.x = element_text(size = 23), axis.text.x = element_text(size = 23),
         axis.title.y = element_blank(), axis.text.y = element_blank(),
         legend.position = c(0.85,.9), legend.text=element_text(size=20),
         legend.background = element_rect(size=0.5,linetype="solid",colour ="black")))


#####Combined inv simp graph#####
ggarrange(trans_inv_simpson_p,seed_inv_simpson_p,ncol = 2,  legend = "none",widths = c(1,.90))
#15x7.38


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