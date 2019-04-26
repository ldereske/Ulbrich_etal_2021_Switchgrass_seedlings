#Here are the column names from the csv files

#SG_harvest_measurement.csv
`
Plant_Number= the pot number
height_cm= height in cm 
leaf_no = number of leaves
coverscale_1_5  = percent of leaf that is shriveled or necrotic (scale needed)
root_length= length of the roots in cm
notes_plantheight= notes from oberserver measuring plant height
notes_rootlength= notes from oberserver measuring root_length
observer= person who measure the given pot
`

#SG_height_time_combine.csv
`
Plant_Number= the pot number
height_cm= height in cm
leaf_num= number of leaves
coverscale_1_5 = percent of leaf that is shriveled or necrotic (scale needed)
notes=notes from oberserver
date= date of the given measurement m_dd
observer= person who measure the given pot or date
exp_period= period within the experiment (weekly_sur= weekly survey by a given oberserver, harvest= harvest of the 
experiment on 7/9/2018, planting= planting of the experiment on 6/12/2018, replacement= replacement of dead 
transplants in the first week)
date_raw = date in the mm/dd/yyy format
surv = observer determined dead plants 1 = alive 0 = dead, only for the harvest period
`

#SG_surv_Seed_germ.csv
`
Plant_Number= the pot number
date= date of the given measurement m_dd
Rdate= date in the mm/dd/yyy format
exp_days = days since beginning of experiment
num_germinates = number of germinants on given day
volunteer_dicot = number of dicots/volunteers removed from the pot
num_removed = number of switchgrass seedlings removed on a given day to return the pot to one seedling (we plant 10 
seeds)
`

#SG_surv_Seedling_replace.csv
`
Plant_Number= the pot number
date= date of the given measurement m_dd
Rdate= date in the mm/dd/yyy format
survival= survival of transplant 1= alive and 0=dead
Transplant_height_cm= height in cm of plant used to replace a dead seedling
Transplant_leaf= number of leaves of plant used to replace a dead seedling
volunteer_dicot = number of dicots/volunteers removed from the pot
volunteer_grass = number of grass/volunteers removed from the pot
`
#SG_TotalBiomass.csv
`
Plant_Number= the pot number
root_weight_g= weight of the harvested root in g
shoot_weight_g= weight of the harvested shoot in g
total_biomass= root_weight_g+shoot_weight_g
notes_roots=notes from oberserver weighing roots
notes_shoots=notes from oberserver weighing shoots
sandy_root= samples that the person weighing flagged as having a lot of sand attached to the roots
`

#SG_trt_sheet.csv
`
Plant_Number= the pot number
raw_label= the raw label with short code of the treatments
soil_status= status of soils was if they were L=alive (not autoclaved) S=sterile (autoclaved)
root_association= where the soil was originally sampled B=bulk (dug with a shovel in patchs without Switchgrass) 
R=rhizosphere soils (soils shook from live plant daug up in the plot)
block= pots were blocked by the plant that they recieved soil from and location in the greenhouse)
precip= A=ambient precipitation (3.4 mL per day 1997-2017 May to July average monthly precipitation = 93.15 mm) or 
D=drought precipitation (3.4 mL every three days 32.85 mm per month based of over the 2012 drought)
life_stage= pot recieved either S=seeds (10 per pot) or G=seedlings transplanted from sterile sand/vermiculite
trt= the raw label with short code of the treatments with pot number removed
`
