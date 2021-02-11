
biolog <- read.csv("D:/MERDS_2018/merds/Switchgrass/R_data/2018.11.01_Merds_SG_biolog_sampleday_absorbance_Matrix.csv", header = TRUE)
metadata <- read.csv("D:/MERDS_2018/merds/Switchgrass/R_data/2018.11.01_Merds_SG_biolog_sampleday_metadata.csv", header = TRUE)
colnames(metadata)
biolog_meta <- merge(metadata, biolog, by = "Sample_day")
dim(biolog_meta) #63 40 

library(RVAideMemoire) # for tests after PERMANOVA
library(usedist)
# subset for day 4  
biolog_meta.4 <- filter(biolog_meta, Day == 4)





#########################################
# PERMANOVA without control

biolog_meta.4_NoCntrl <- filter(biolog_meta.4, soil_root_association != "Control")


# Euclidean
set.seed(2) # this allows you to get the same dissimilarity matrix again and again
biolog_meta.4_NoCntrl_euc<- vegdist(biolog_meta.4_NoCntrl[,c(10:38)], method = "euclidian", binary = FALSE, diag = FALSE, upper = FALSE, na.rm = FALSE)

# Bray 
set.seed(2) # this allows you to get the same dissimilarity matrix again and again
biolog_meta.4_NoCntrl_bray<- vegdist(biolog_meta.4_NoCntrl[,c(10:38)], method = "bray", binary = FALSE, diag = FALSE, upper = FALSE, na.rm = FALSE)


# check betadisperson (differences in spead among the treatments)
betadisp <- betadisper(biolog_meta.4_NoCntrl_euc,group = biolog_meta.4_NoCntrl$soil_root_association,type = "centroid")
permutest(betadisp) 

#Response: Distances
#Df  Sum Sq  Mean Sq     F N.Perm Pr(>F)
#Groups     2 0.08828 0.044138 0.281    999  0.754
#Residuals 15 2.35605 0.157070

boxplot(betadisp)  

permanova <- adonis(biolog_meta.4_NoCntrl_euc~ biolog_meta.4_NoCntrl$soil_root_association*biolog_meta.4_NoCntrl$precip, permutations = 9999, method = "bray")
permanova 

### EUCLIDIAN 
#                                                  Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
#soil_root_association                           2    11.463  5.7314  #3.3595 0.31108  0.005 **
#precip                                              1     1.601  1.6015  0.9387 0.04346  0.477 #soil_root_association:biolog_meta.4_NoCntrl$precip  2     3.312  1.6560  0.9707 0.08988  0.464 
#Residuals                                                 12    20.472  1.7060         0.55558 #Total                                                          17    36.849           1.00000

#### BRAY 
#                                                Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)   
#soil_root_association                           2   0.25325 0.126623  3.5025 0.31900  0.003 **
#precip                                              1   0.04198 0.041976  1.1611 0.05288  0.297
#soil_root_association*$precip                       2   0.06482 0.032408  0.8964 0.08165  0.495 Residuals                                                12   0.43383 0.036152         0.54648 
#Total                                                17   0.79387                  1.00000          

permanova_block <- adonis(biolog_meta.4_NoCntrl_euc~ biolog_meta.4_NoCntrl$soil_root_association*biolog_meta.4_NoCntrl$precip+
                            as.factor(biolog_meta.4_NoCntrl$block), permutations = 9999, method = "bray")
permanova_block



biolog_meta.4_NoCntrl_euc_ord=metaMDS(biolog_meta.4_NoCntrl_euc)
#Run 20 stress 0.094445 
#*** Solution reached
#0.08842766

biolog_meta.4_NoCntrl_euc_ord_points=biolog_meta.4_NoCntrl_euc_ord$points

biolog_meta.4_NoCntrl_euc_ord_points_trt=
  data.frame(cbind(biolog_meta.4_NoCntrl_euc_ord_points,biolog_meta.4_NoCntrl[,c("soil_status","soil_root_association",
                                                            "root_association","block","trt")]))
ggplot(biolog_meta.4_NoCntrl_euc_ord_points_trt, aes(x=MDS1,y=MDS2))+geom_point(size=2,aes(color=soil_root_association))+
  theme_bw()
 
###################################3
# PAIRWISE PERMANOVA 

#What I have used in the past for pairwise.adonis
pairwise.perm.manova(biolog_meta.4_NoCntrl_euc, biolog_meta.4_NoCntrl$soil_root_association, nperm=2000)


  # Pairwise Adonis 
  #install.packages('devtools')
  #library(devtools)
  install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)

Permanova.matrix <- as.matrix(biolog_meta.4_NoCntrl_euc)
row.names(Permanova.matrix)
vector <- as.vector(biolog_meta.4_NoCntrl$soil_root_association)

pairwise.adonis(Permanova.matrix, vector, p.adjust.m = "fdr", perm = 2000)


#let's look at the distance between centroids

soil_dis_centroids=dist_groups(biolog_meta.4_NoCntrl_euc,biolog_meta.4_NoCntrl$soil_root_association)

soil_dis_centroids %>% group_by(Label) %>% summarise_at(vars(Distance), list(~mean(.),~sd(.)))
