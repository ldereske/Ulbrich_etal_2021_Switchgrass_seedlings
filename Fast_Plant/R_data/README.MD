# Here are the column names and file info


### Fastplant_treatment_metadata and FastPlant_trt (formated for R)

+ Plant_Number = unique pot ID # from 1-270, links with all data
+ soil_status = L = live inocula added, S = sterile inocula added 
+ block = block 1-5, treatments grouped by block location in the greenhouse 
+ replicate = Each block/treatment/soil_status combination has 3 replicates, listed a 'a-c'
+ precip =  A=ambient precipitation (3.4 mL per day 1997-2017 May to July average monthly precipitation = 93.15 mm) or D=drought precipitation (3.4 mL every three days 32.85 mm per month based of over the 2012 drought)
+ inocul = 20mL of inocula that was added to the sterilized sand/vermiculte: SDS = commercial product spectrum DS (https://www.tainio.com/product/spectrum-ds-biological-soil-inoculant/) that aims to improve plant growth and resistance to drought; PM = 'product microbes' or  5 bacterial strains that we were able to culture and isolate from spectrum DS (important difference that there is no additional organic matter as there is with SDS); GI= 'good isolates' or 3 fungi and 3 bacteria that were shown to be drought resistant; BI= 'bad isolates' or 3 fungi and 3 bacteria that were found to be drought sensitive; S='sterile' autoclaved nanopure water; more details provided in inocula details textfile 
+ trt = Descriptive label used to describe all treatments (SoilStatus/block/replicate/precip_Inocul)


### Fastplant_biomass_data and Fast_plant_harvest_weight_R (formated for R)

+ Plant_Number = unique pot ID # from 1-270, links with treatment metadata
+ root_weight = dry root weight (g), in some cases there was no plant because it died/empty pot, in this case == 0, in some cases there was no data recorded (not sure why at this moment), then == NA
+ notes_root = notes for root_weight column
+ observer_root = individual that weighed roots & took notes, NA if data are missing
+ date_weighed_root = date that the roots were weighed (they were dried before weighing, but weighing took place over several months), NA if data are missing
+ entry_notes_root = additional notes from root weighing 
+ sandy_root = Y = presence of lots of sand/vermiculite on root (this could have increased weight), N = root is not sandy
+ shoot_weight = dry shoot weight (g), in some cases there was no plant because it died/empty pot, in this case == 0, in some cases there was no data recorded (not sure why at this moment), then == NA
+ notes_shoot = notes for shoot_weight column
+ observer_shoot = individual that weighed shoots & took notes, NA if data are missing
+ date_weighed_shoot = date that the shoots were weighed (they were dried before weighing, but weighing took place over several months), NA if data are missing
+ shoot_height_cm_harvest = shoot height (cm) measured day of harvest
+ coverscale_harvest = Qualitative estimate of % of leaves that were brown measured on harvest day, Coverscale was recorded on a scale of 0-5 with approximate percentage ranges of 0%, 1-5%, 6-25%, 26-75%, and 75-99%, and “dead”
+ harvest_notes = Additional notes from harvest, details plant health, etc.
+ survival_harvest = quantitiative index of survival at harvest, 1 = survived, 0 = dead, NA = no data
+ buds_714 = # of buds on plant on July 14th, NA if no data 
+ pods_714 = # of seed pods on plant on July 14th, NA if no data
+ shoot_height_cm_710 = Shoot height (cm) measured on July 10th, NA if no data
+ buds_710 = # of buds on plant on July 10th, NA if no data
+ pods_710 = # of seed pods on plant on July 10th, NA if no data
+ coverscale_710 = Qualitative estimate of % of leaves that were brown measured on July 10th, Coverscale was recorded on a scale of 0-5 with approximate percentage ranges of 0%, 1-5%, 6-25%, 26-75%, and 75-99%, and “dead”
+ notes_710 = Notes collected on July 10th, linked to other measurements taken on July 10th
+ shoot_height_cm_702 = shoot height (cm) measured on July 2nd, NA if no data
+ leaf_no_702 = # of leaves on July 2nd
+ leaf_length_702 = Leaf length (cm) on July 2nd  - I'm not sure why there are so many missing data here or how we measured leaf length 
+ coverscale_702 = Qualitative estimate of % of leaves that were brown measured on July 2nd, Coverscale was recorded on a scale of 0-5 with approximate percentage ranges of 0%, 1-5%, 6-25%, 26-75%, and 75-99%, and “dead”
+ notes_702 = Notes collected on July 2nd, linked to other measurements taken on July 2nd


### Fastplant_potweight

###### Weight of pots (including soils and plants) at beginning of expeirment and end of experiment. This can be used as an estimate for soil moisture, though will be confounded by plant biomass changes 
+ Plant_Number = unique pot ID # from 1-270, links with treatment metadata
+ Potweight_final_g = Final pot weight (g) on July 11th, 2018
+ Potweight_initial_g = Initial pot weight (g) on June 20th, 2018


### FastPlant_FastPlant_soilNassay_RawAbsorbanceData

###### Ammonium and Nitrate concentrations were determined using previously described 96-well plate assays with soil extracts from unfumigated .5M K2SO4 extracts;absorption values were read on a Synergy HTX plate reader (BioTek, Winooski, Vermont, USA)(ammonium - Sinsabaugh et al., 2000; nitrate - Patton & Kryskalla, 2011) Data in this file are the raw absorbance data and need to be converted to ppm using standard curve 
+ Sample = Plant_Number linked to metadata; STD# = standard curve reagent at ppm nitrate or ammonium; C# = control sample of just K2SO4; Empty - no sample; 'MMPRNT' - samples from another project, disregard this.
+ Well = Well in the 96well plate
+ 610nm = 610 nm absorbance for ammonium assay that should be converted to ppm NH4 
+ 540nm = 540 nm absorbance for nitrate assay that should be converted to ppm NO3 
+ 700nm = 700 nm absorbance for nitrate assay that should be converted to ppm NO3 
+ DilutionFactor = dilution factor of sample, 1 = no dilution of K2SO4 extracts for 
+ Plate = Plate # (there were 2 plates used for all samples for either assay; standard curves and controls should be unique to each plate) 
