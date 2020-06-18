
######################################
#### Clonality plot ####
######################################
#~/Rscript3.5.0 ~/res/scripts-orphan/get_ccf_subclonality_all_batches_expected_copies.R current_wes_somatic50.maf current_wes_somatic50.ccf.maf
######
##Clinical track
cli_info<-fread("~/bladder-paper/wes_clinical.csv") %>% filter(patient.id != "P-0005150")

patient.order<-fread("~/bladder-paper/patient_order_WES.csv") %>% filter(V1 != "P-0005150")

status_order<-c("Metastasis_Location","Type_Treatment","Time_Between")
cli_info$key=factor(cli_info$key,level=c(unique(cli_info$key)))
cli_info$value=factor(cli_info$value,level=c(unique(cli_info$value)))
cli_info$patient.id<-factor(cli_info$patient.id)

color_map<- c("#D9D9D9","#969696", "#FFF7EC","#FDD49E","#FC8D59","#D7301F","#7F0000","#DEEBF7" ,"#9ECAE1", "#4292C6", "#08519C")

track<- ggplot(cli_info, aes(x=patient.id,y=key))+ 
  geom_tile(aes(fill=value),width=0.9,height=0.9)+
  scale_y_discrete(limits=rev(status_order))+
  scale_x_discrete(limits=patient.order$V1) +
  scale_fill_manual(values=color_map,name = "Clinical Information",guide="legend")+ 
  guides(color = guide_legend(ncol = 1)) +
  theme(legend.key = element_rect(size = 2, color = "white"),
        legend.key.size = unit(1.5, "lines")) +
  theme(text = element_text(color = "gray20"),
        legend.title = element_blank(),
        legend.position = c("bottom"), # position the legend in the upper left
        #legend.justification = 0, # anchor point for legend.position.
        legend.text = element_text(size = 12, color = "gray10"),
        title = element_text(size = 12, face = "bold", color = "gray10"),
        axis.text.x = element_text(size =12, face = "bold",angle=90),
        axis.text.y = element_text(size =12, face = "bold"),
        #axis.text.y = element_blank(),
        # panel.grid.major.y = element_blank(),
        # panel.grid.major.x = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()
       # panel.grid.minor.y=element_line(colour='black', size=0),
       # panel.grid.major.y=element_line(colour='black', size=0),
       # panel.grid.minor.x=element_line(colour='black', size=0),
       # panel.grid.major.x=element_line(colour='black', size=0),
       # plot.margin = unit(c(1,1,1,1), 'lines')
       )
track
##############################
##############################
######################################
#### Clonality plot ####
######################################
#~/Rscript3.5.0 ~/res/scripts-orphan/get_ccf_subclonality_all_batches_expected_copies.R current_wes_somatic50.maf current_wes_somatic50.ccf.maf

#wesclonalitymaf = fread('/Users/chavans/juno/work/ccs/ccs_wes/Proj_07871_DFLOQ/Result/somatic/allzygositymafs.Private.clonality.maf')
wesclonalitymaf0 =  fread('~/bladder-paper/bladder_mut_somatic.ccf.freeze.050520.maf')
dim(wesclonalitymaf0)
length(Private(wesclonalitymaf0$Tumor_Sample_Barcode))
names(wesclonalitymaf0)

# wesclonalitymaf = wesclonalitymaf %>%
#   mutate(
#          clonal_call = ifelse(ccf_expected_copies_em > 0.8 | (ccf_expected_copies_em > 0.7 & ccf_expected_upper_em > 0.9), TRUE, FALSE),
#          clonal_call = ifelse(cf.em < (0.6*purity), "INDETERMINABLE", clonal_call)
#          )
# #table(wesclonalitymaf$clonal_call)
# wesclonalitymaf0 = wesclonalitymaf %>% select(-c(cf, tcn, lcn))
# setnames(wesclonalitymaf0, "cf.em","cf")
# setnames(wesclonalitymaf0, "lcn.em","lcn")
# setnames(wesclonalitymaf0, "tcn.em","tcn")
# wesclonalitymaf0 = wesclonalitymaf0 %>% select(-c(241:250))
# names(wesclonalitymaf0)
# wesclonalitymaf0 = wesclonalitymaf0 %>% rename_at(vars(ends_with("_em")), 
#                                                   funs(str_replace(., "_em", "")))
# names(wesclonalitymaf0)
# setnames(wesclonalitymaf0, "ccf_expected_upper","ccf_expected_copies_upper")
# setnames(wesclonalitymaf0, "ccf_expected_prob90","ccf_expected_copies_prob90")
# setnames(wesclonalitymaf0, "ccf_expected_lower","ccf_expected_copies_lower")
# setnames(wesclonalitymaf0, "ccf_expected_prob95","ccf_expected_copies_prob95")

# wesmaf2 = fread('/ifs/res/taylorlab/chavans/bladder_kdm6a/matched_pair_batch2/mut_somatic.maf') 
# 
# wesmaf0 = wesmaf2 %>%
#   mutate(
#   clonal_call = ifelse(ccf_expected_copies > 0.8 | (ccf_expected_copies > 0.7 & ccf_expected_copies_upper > 0.9), TRUE, FALSE),
#   clonal_call = ifelse(cf < (0.6*purity), "INDETERMINABLE", clonal_call)
#   )
# dim(wesmaf0)
# table(wesmaf0$clonal_call)
# length(Private(wesmaf0$Tumor_Sample_Barcode))
# 
# combined_maf = rbind.fill(wesclonalitymaf0, wesmaf0)
# length(Private(combined_maf$Tumor_Sample_Barcode))
# table(combined_maf$clonal_call)
# dim(combined_maf)

#write.table(combined_maf, '/ifs/res/taylorlab/chavans/bladder_kdm6a/bladder_mut_somatic.ccf.freeze.050520.maf', sep = "\t", row.names = FALSE, quote = FALSE, append = FALSE)

# 
# filter(wesclonalitymaf, Tumor_Sample_Barcode %in% NAsamples$Tumor_Sample_Barcode, tcn.em == 2, lcn.em == 1) %>%
#   group_by(Tumor_Sample_Barcode) %>%
#   dplyr::summarise(MutationBasedPurity = 2*median(t_var_freq))

# # A tibble: 5 x 2
# Tumor_Sample_Barcode MutationBasedPurity
# <chr>                              <dbl>
# 1 s_C_271D5P_M001_d                  0.116
# 2 s_C_271D5P_P001_d                  0.118
# 3 s_C_ADEV70_P001_d                  0.232
# 4 s_C_PDMVDR_P001_d                  0.118
# 5 s_C_W384MJ_M001_d                  0.129
# 
# filter(wesclonalitymaf, Tumor_Sample_Barcode %in% NAsamples$Tumor_Sample_Barcode, tcn.em == 2, lcn.em == 1) %>%
#   arrange(Tumor_Sample_Barcode) %>%
#   select(Tumor_Sample_Barcode, t_var_freq)


head(wesclonalitymaf0); dim(wesclonalitymaf0)


wesclonalitymaf = wesclonalitymaf0 %>% mutate(CLONALITY = ifelse(clonal_call == TRUE, TRUE, FALSE), 
                                             VAR_TAG = str_c(Chromosome,':',Start_Position,':',End_Position,':',Reference_Allele,':',Tumor_Seq_Allele2),
                                             PTN = substring(Tumor_Sample_Barcode,1,10),
                                             SAMPLE_TYPE = ifelse(substring(Tumor_Sample_Barcode,11,13)=="_M0","Metastasis","Primary"))

wesclonalitymaf %>% select(Tumor_Sample_Barcode, SAMPLE_TYPE) %>% table

#Fix the primary and met switch s_C_000796
wesclonalitymaf = wesclonalitymaf %>% mutate(SAMPLE_TYPE = ifelse(Tumor_Sample_Barcode %like% "s_C_000796_M0","Primary",SAMPLE_TYPE))
wesclonalitymaf = wesclonalitymaf %>% mutate(SAMPLE_TYPE = ifelse(Tumor_Sample_Barcode %like% "s_C_000796_P0","Metastasis",SAMPLE_TYPE))
wesclonalitymaf %>% select(Tumor_Sample_Barcode, SAMPLE_TYPE) %>% table

wes_vartag_table = wesclonalitymaf %>% distinct(.) %>% select(PTN, VAR_TAG) %>% group_by(PTN,  VAR_TAG) %>% dplyr::summarize(TAG_COUNT_PER_PTN = n())

wesclonalitymaf0 = inner_join(wesclonalitymaf, wes_vartag_table, by = c('PTN','VAR_TAG')) %>%
                   mutate(MUT_STATUS = ifelse(TAG_COUNT_PER_PTN >1, 'Shared', 'Private'))

sampleqc = fread("~/bladder-paper/SampleListPreAndPostQC-050520.txt") %>% filter(Exclude == 0)
dim(sampleqc) #51
wesclonalitymaf0_fil = filter(wesclonalitymaf0, PTN %in% sampleqc$Patient_ID, 
                                                Tumor_Sample_Barcode %in% sampleqc$TUMOR_ID) %>%
                       filter(PTN !="s_C_ADEV70") %>% #FACETS not available for the primary sample
                       filter(Tumor_Sample_Barcode != "s_C_FR895T_P003_d") %>% # Only 1 common mutation with P001
                       filter(Tumor_Sample_Barcode %ni% c('s_C_000796_M001_d', 's_C_002086_T001_d', 's_C_002088_T002_d', 's_C_FR895T_P002_d')) #Multiple primary or met per patient   
unique(wesclonalitymaf0_fil$Tumor_Sample_Barcode) #44

##Replace CMO ids with DMP ids
wesclonalitymaf0_fil = wesclonalitymaf0_fil %>% mutate(PTN = mapvalues(Tumor_Sample_Barcode, sampleqc$TUMOR_ID, sampleqc$DMP_PatID))
wesclonalitymaf0_fil %>% select(Tumor_Sample_Barcode,PTN) %>% distinct()

#filter(wesclonalitymaf0_fil, PTN =="s_C_LF6FEV") %>% select(Tumor_Sample_Barcode, SAMPLE_TYPE, CLONALITY) %>% table() #1 mut diff is okay.
#Rest are multiple samples per patient, ask Tim if we should pick a relevant one
# s_C_000796
# s_C_002086
# s_C_002088
# s_C_FR895T
# ##Use these:
# s_C_000796_M002_d
# s_C_000796_P001_d
# 
# s_C_002086_M001_d
# s_C_002086_T002_d
# 
# s_C_002088_M001_d
# s_C_002088_T002_d
# 
# s_C_FR895T_M001_d
# s_C_FR895T_P001_d
##Remove these:
#s_C_000796_M001_d, s_C_002086_T001_d, s_C_002088_T002_d, s_C_FR895T_P002_d
# wesclonalitymaf_shared_clonal = filter(wesclonalitymaf0_fil, ccf.call == "Clonal") %>% select(PTN, VAR_TAG) %>% group_by(PTN,  VAR_TAG) %>%
#   summarize(TAG_COUNT_PER_PTN = n()) %>% mutate(MUT_STATUS = ifelse(TAG_COUNT_PER_PTN >1, 'SHARED_CLONAL', 'UNIQUE_CLONAL')) %>% 
#   dcast(PTN ~ MUT_STATUS) %>% mutate(UNION_CLONAL = SHARED_CLONAL + UNIQUE_CLONAL)

###############
####### PLOT1 START1 #######
############
## Shared
#wesclonalitymaf_shared = wesclonalitymaf0_fil %>% select(PTN, MUT_STATUS) %>% distinct(.) %>%
#  dcast(PTN ~ MUT_STATUS) %>% mutate(UNION = Shared + Private, Fraction_Shared = Shared/UNION)
#group_by(PTN, Shared) %>% dplyr::dplyr::summarise(n()) %>% group_by(PTN) %>% dplyr::dplyr::summarise(sum(`n()`))
#select(Matched_Norm_Sample_Barcode, Tumor_Sample_Barcode, VAR_TAG)

wesclonalitymaf_shared0 = wesclonalitymaf0_fil %>% filter(MUT_STATUS == 'Shared') %>% 
  select(PTN, VAR_TAG) %>% distinct(.) %>% group_by(PTN) %>% dplyr::summarise(Shared = n())
wesclonalitymaf_unique = wesclonalitymaf0_fil %>% filter(MUT_STATUS != 'Shared') %>% 
  select(PTN, VAR_TAG) %>% distinct(.) %>% group_by(PTN) %>% dplyr::summarise(Private = n())
wesclonalitymaf_union = wesclonalitymaf0_fil %>%
  select(PTN, VAR_TAG) %>% distinct(.) %>% group_by(PTN) %>% dplyr::summarise(UNION = n())
wesclonalitymaf_shared = join_all(list(wesclonalitymaf_shared0, 
                                        wesclonalitymaf_unique, 
                                        wesclonalitymaf_union
), by = c('PTN'), type = "left")
  
## Private PRI MET
wesclonalitymaf_pri_unique = wesclonalitymaf0_fil %>% select(PTN, SAMPLE_TYPE, VAR_TAG, CLONALITY, MUT_STATUS) %>% 
  filter(SAMPLE_TYPE == "Primary", MUT_STATUS == "Private"
  ) %>% 
  group_by(PTN) %>%
  dplyr::summarize(Private_Primary = n()) 

wesclonalitymaf_met_unique = wesclonalitymaf0_fil %>% select(PTN, SAMPLE_TYPE, VAR_TAG, CLONALITY, MUT_STATUS) %>% 
  filter(SAMPLE_TYPE == "Metastasis", MUT_STATUS == "Private") %>% 
  group_by(PTN) %>%
  dplyr::summarize(Private_Metastasis = n()) 

## CALCULATE FRACTIONS
## JOIN ALL
wessubclonality = join_all(list(wesclonalitymaf_shared, 
                                wesclonalitymaf_pri_unique, 
                                wesclonalitymaf_met_unique
), by = c('PTN'), type = "left") 
names(wessubclonality); 
dim(wessubclonality)#22

## CALCULATE
wessubclonality_fracs = wessubclonality %>% mutate( 
  Shared = Shared/UNION,
  Private_Primary = Private_Primary/UNION, 
  Private_Metastasis = Private_Metastasis/UNION) %>%
  select(PTN,Shared,Private_Primary,Private_Metastasis)

wessubclonality_fracs

dim(wessubclonality_fracs)

wessubclonality_fracs_ = wessubclonality_fracs %>% select(names(wessubclonality_fracs)[1],
                                                          names(wessubclonality_fracs)[2:4]) %>% replace(., is.na(.), 0)
wessubclonality_fracs_.m = melt(wessubclonality_fracs_)
wessubclonality_fracs_.m$variable = factor(wessubclonality_fracs_.m$variable,
                                           levels=rev(c("Shared","Private_Primary","Private_Metastasis")))

## set the levels in order we want
wessubclonality_fracs_.m$PTN <- factor(wessubclonality_fracs_.m$PTN,levels = patient.order$V1)

p1 = ggplot(data=wessubclonality_fracs_.m) +
  geom_bar(aes(x=reorder(PTN,variable), y=value, fill=wessubclonality_fracs_.m$variable),stat='identity') +
  scale_fill_jco() + 
  #scale_x_discrete(limits=patient.order$V1) +
  #values = rev(c('tan4','tan','navy','lightblue','black','grey')) )+ #top to bottom
  ylab('Fraction of mutations') + xlab('') +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_blank(),
        axis.title.y = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0)),
    #axis.text.x = element_text(angle=90, hjust=1, size=NULL),
        legend.key.size = unit(.5, "cm"),
        legend.text = element_text(size = 12, face = "bold", color = "gray10"),
        #legend.text = element_text(size=9),
        legend.position='top',
        legend.title=element_blank(),
        panel.grid.minor.y=element_line(colour='black', size=0),
        panel.grid.major.y=element_line(colour='black', size=0),
        panel.grid.minor.x=element_line(colour='black', size=0),
        panel.grid.major.x=element_line(colour='black', size=0),
        plot.margin = unit(c(1,1,1,1), 'lines')) 
p1
#track

somePDFPath = "~/bladder-paper/Figure_WES_Shared_Unique_061820_part1.pdf"
pdf(file=somePDFPath)   
#grid.arrange(p1,p2,p3, ncol = 1, newpage = FALSE)
cowplot::plot_grid(p1,track, align = "v", ncol = 1, rel_heights = c(.6,.4))
dev.off()

###############
####### PLOT1 END #######
############

###############
####### PLOT2 START #######
############

## Shared CLONAL
wesclonalitymaf_shared_clonal_in_both = wesclonalitymaf0_fil %>% select(PTN, SAMPLE_TYPE, VAR_TAG, CLONALITY, MUT_STATUS) %>% 
  filter(MUT_STATUS == "Shared", CLONALITY == TRUE) %>% 
  group_by(PTN,VAR_TAG) %>%
  dplyr::summarise(TAG_COUNT_PER_PTN = n()) %>%
  mutate(BOTH = ifelse(TAG_COUNT_PER_PTN >=2, TRUE, FALSE)) %>% 
  filter(BOTH == TRUE) %>% distinct() %>%
  dplyr::summarize(CLONAL_BOTH = n())

## Shared SUBCLONAL
wesclonalitymaf_shared_sub_in_both = wesclonalitymaf0_fil %>% select(PTN, SAMPLE_TYPE, VAR_TAG, CLONALITY, MUT_STATUS) %>% 
  filter(MUT_STATUS == "Shared", CLONALITY == FALSE) %>% 
  group_by(PTN,VAR_TAG) %>%
  dplyr::summarise(TAG_COUNT_PER_PTN = n()) %>%
  mutate(BOTH = ifelse(TAG_COUNT_PER_PTN >=2, TRUE, FALSE)) %>%
  filter(BOTH == TRUE) %>% distinct() %>%
  dplyr::summarize(SUB_BOTH = n())

## Shared CLONAL IN ONE SUB IN OTHER
wesclonalitymaf_shared_clonal_in_pri = wesclonalitymaf0_fil %>% select(PTN, SAMPLE_TYPE, VAR_TAG, CLONALITY, MUT_STATUS) %>% 
  filter(MUT_STATUS == "Shared", CLONALITY == TRUE, SAMPLE_TYPE == "Primary") %>% 
  group_by(PTN,VAR_TAG) %>%
  dplyr::summarise(TAG_COUNT_PER_PTN = n()) %>%
  dplyr::summarize(CLONAL_PRI = n())

wesclonalitymaf_shared_clonal_in_met = wesclonalitymaf0_fil %>% select(PTN, SAMPLE_TYPE, VAR_TAG, CLONALITY, MUT_STATUS) %>% 
  filter(MUT_STATUS == "Shared", CLONALITY == TRUE, SAMPLE_TYPE == "Metastasis") %>% 
  group_by(PTN,VAR_TAG) %>%
  dplyr::summarise(TAG_COUNT_PER_PTN = n()) %>%
  dplyr::summarize(CLONAL_MET = n())

## CALCULATE FRACTIONS
## JOIN ALL
wessubclonality_allshared = join_all(list(wesclonalitymaf_shared, 
                                          wesclonalitymaf_shared_clonal_in_both,
                                          wesclonalitymaf_shared_sub_in_both, 
                                          wesclonalitymaf_shared_clonal_in_pri, 
                                          wesclonalitymaf_shared_clonal_in_met
), by = c('PTN'), type = "left") %>% mutate(
  CLONAL_PRI_ONLY = CLONAL_PRI - CLONAL_BOTH,
  CLONAL_MET_ONLY = CLONAL_MET - CLONAL_BOTH
) %>% select(PTN, Shared, Private, UNION, CLONAL_BOTH, SUB_BOTH, CLONAL_PRI_ONLY, CLONAL_MET_ONLY)
names(wessubclonality_allshared); 
dim(wessubclonality_allshared)

## CALCULATE
wessubclonality_fracs = wessubclonality_allshared %>% mutate( 
  Clonal_Both = CLONAL_BOTH/UNION,
  Subclonal_Both = SUB_BOTH/UNION, 
  Clonal_Primary = CLONAL_PRI_ONLY/UNION,
  Clonal_Metastasis = CLONAL_MET_ONLY/UNION) %>%
  select(PTN,Clonal_Both,Subclonal_Both,Clonal_Primary,Clonal_Metastasis)

wessubclonality_fracs

dim(wessubclonality_fracs)

wessubclonality_fracs_ = wessubclonality_fracs %>% select(names(wessubclonality_fracs)[1],
                                                          names(wessubclonality_fracs)[2:5]) %>% replace(., is.na(.), 0)
wessubclonality_fracs_.m = melt(wessubclonality_fracs_)
wessubclonality_fracs_.m$variable = factor(wessubclonality_fracs_.m$variable,
                                           levels=rev(c("Clonal_Both","Subclonal_Both","Clonal_Primary","Clonal_Metastasis")))

## set the levels in order we want
wessubclonality_fracs_.m$PTN <- factor(wessubclonality_fracs_.m$PTN,levels = patient.order$V1)
p2 = ggplot(data=wessubclonality_fracs_.m) +
  geom_bar(aes(x=PTN, y=value, fill=wessubclonality_fracs_.m$variable),stat='identity') +
  scale_fill_nejm() + 
  #scale_x_discrete(limits=patient.order$V1) +
  #values = rev(c('tan4','tan','navy','lightblue','black','grey')) )+ #top to bottom
  ylab('Fraction of mutations') + xlab('') +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_blank(),
        axis.title.y = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0)),
        #axis.text.x = element_text(angle=90, hjust=1, size=NULL),
        legend.key.size = unit(.5, "cm"),
        legend.text = element_text(size = 9, face = "bold", color = "gray10"),
        #legend.text = element_text(size=9),
        legend.position='top',
        legend.title=element_blank(),
        plot.margin = unit(c(1,1,1,1), 'lines')
        ) 
p2

somePDFPath = "~/bladder-paper/Figure_WES_Shared_Unique_061820_part2.pdf"
pdf(file=somePDFPath)   
#grid.arrange(p1,p2,p3, ncol = 1, newpage = FALSE)
cowplot::plot_grid(p2,track, align = "v", ncol = 1, rel_heights = c(0.6,.4))
dev.off()  
###############
####### PLOT2 END #######
############

###############
####### PLOT3 START #######
############
## Private CLONAL PRI MET
wesclonalitymaf_pri_unique = wesclonalitymaf0_fil %>% select(PTN, SAMPLE_TYPE, VAR_TAG, CLONALITY, MUT_STATUS) %>% 
  filter(SAMPLE_TYPE == "Primary", MUT_STATUS == "Private", CLONALITY == TRUE) %>% 
  group_by(PTN) %>%
  dplyr::summarize(UNIQUE_PRI_CLONAL = n())

wesclonalitymaf_met_unique = wesclonalitymaf0_fil %>% select(PTN, SAMPLE_TYPE, VAR_TAG, CLONALITY, MUT_STATUS) %>% 
  filter(SAMPLE_TYPE == "Metastasis", MUT_STATUS == "Private", CLONALITY == TRUE) %>% 
  group_by(PTN) %>%
  dplyr::summarize(UNIQUE_MET_CLONAL = n())

## Private SUBCLONAL PRI MET
wesclonalitymaf_pri_unique_sub = wesclonalitymaf0_fil %>% select(PTN, SAMPLE_TYPE, VAR_TAG, CLONALITY, MUT_STATUS) %>% 
  filter(SAMPLE_TYPE == "Primary", MUT_STATUS == "Private", CLONALITY == FALSE) %>% 
  group_by(PTN) %>%
  dplyr::summarize(UNIQUE_PRI_SUB = n())

wesclonalitymaf_met_unique_sub = wesclonalitymaf0_fil %>% select(PTN, SAMPLE_TYPE, VAR_TAG, CLONALITY, MUT_STATUS) %>% 
  filter(SAMPLE_TYPE == "Metastasis", MUT_STATUS == "Private", CLONALITY == FALSE) %>% 
  group_by(PTN) %>%
  dplyr::summarize(UNIQUE_MET_SUB = n()) 

## CALCULATE FRACTIONS
## JOIN ALL
wessubclonality_allunique = join_all(list(
                                wesclonalitymaf_shared,
                                wesclonalitymaf_pri_unique, 
                                wesclonalitymaf_met_unique,
                                wesclonalitymaf_pri_unique_sub,
                                wesclonalitymaf_met_unique_sub
), by = c('PTN'), type = "left") #%>% select(-c(BaitSet)) 
names(wessubclonality_allunique); 
dim(wessubclonality_allunique)
## CALCULATE
  wessubclonality_fracs = wessubclonality_allunique %>% mutate( 
                                Clonal_Primary = UNIQUE_PRI_CLONAL/UNION,
                                Clonal_Metastasis = UNIQUE_MET_CLONAL/UNION, 
                                Subclonal_Primary = UNIQUE_PRI_SUB/UNION,
                                Subclonal_Metastasis = UNIQUE_MET_SUB/UNION) %>%
    select(PTN,Clonal_Primary,Clonal_Metastasis,Subclonal_Primary,Subclonal_Metastasis)
  names(wessubclonality_fracs);
  dim(wessubclonality_fracs)
  
  wessubclonality_fracs_ = wessubclonality_fracs %>% select(names(wessubclonality_fracs)[1],names(wessubclonality_fracs)[2:5]) %>% replace(., is.na(.), 0)
  wessubclonality_fracs_.m = melt(wessubclonality_fracs_)
  
##PLOT IT AWAY  
wessubclonality_fracs_.m$variable = factor(wessubclonality_fracs_.m$variable,
                                             levels=rev(c("Clonal_Primary","Clonal_Metastasis",
                                                          "Subclonal_Primary","Subclonal_Metastasis")))
## set the levels in order we want
  wessubclonality_fracs_.m$PTN <- factor(wessubclonality_fracs_.m$PTN,levels = patient.order$V1)
  
p3 = ggplot(data=wessubclonality_fracs_.m) +
    geom_bar(aes(x=PTN, y=value, fill=wessubclonality_fracs_.m$variable),stat='identity') +
    scale_fill_aaas() + 
      #values = rev(c('tan4','tan','navy','lightblue','black','grey')) )+ #top to bottom
  ylab('Fraction of mutations') + xlab('') +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_blank(),
        axis.title.y = element_text(margin = margin(t = 0, r = 0.25, b = 0, l = 0)),
        #axis.text.x = element_text(angle=90, hjust=1, size=NULL),
        legend.key.size = unit(.5, "cm"),
        legend.text = element_text(size = 8, face = "bold", color = "gray10"),
        #legend.text = element_text(size=9),
        legend.position='top',
        legend.title=element_blank(),
        plot.margin = unit(c(1,1,1,1), 'lines')
  )
p3


somePDFPath = "~/bladder-paper/Figure_WES_Shared_Unique_061820_part3.pdf"
pdf(file=somePDFPath)   
#grid.arrange(p1,p2,p3, ncol = 1, newpage = FALSE)
cowplot::plot_grid(p3,track, align = "v", ncol = 1, rel_heights = c(0.6,.4))
#rel_heights = c(0.25,.25,.25,.25)
dev.off()  

###############
####### PLOT3 END #######
############


