
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
#wesclonalitymaf = fread('/Users/chavans/juno/work/ccs/ccs_wes/Proj_07871_DFLOQ/Result/somatic/allzygositymafs.unique.clonality.maf')
wesclonalitymaf0 =  fread('~/bladder-paper/bladder_mut_somatic.ccf.freeze.050520.maf')
dim(wesclonalitymaf0)
length(unique(wesclonalitymaf0$Tumor_Sample_Barcode))
names(wesclonalitymaf0)

head(wesclonalitymaf0); dim(wesclonalitymaf0)


wesclonalitymaf = wesclonalitymaf0 %>% mutate(CLONALITY = ifelse(clonal_call == TRUE, TRUE, FALSE), 
                                             VAR_TAG = str_c(Chromosome,':',Start_Position,':',End_Position,':',Reference_Allele,':',Tumor_Seq_Allele2),
                                             PTN = substring(Tumor_Sample_Barcode,1,10),
                                             SAMPLE_TYPE = ifelse(substring(Tumor_Sample_Barcode,11,13)=="_M0","Metastasis","Primary"))


#Fix the primary and met switch s_C_000796
wesclonalitymaf = wesclonalitymaf %>% mutate(SAMPLE_TYPE = ifelse(Tumor_Sample_Barcode %like% "s_C_000796_M0","Primary",SAMPLE_TYPE))
wesclonalitymaf = wesclonalitymaf %>% mutate(SAMPLE_TYPE = ifelse(Tumor_Sample_Barcode %like% "s_C_000796_P0","Metastasis",SAMPLE_TYPE))
wesclonalitymaf %>% select(Tumor_Sample_Barcode, SAMPLE_TYPE) %>% table

wes_vartag_table = wesclonalitymaf %>% distinct(.) %>% select(PTN, VAR_TAG) %>% group_by(PTN,  VAR_TAG) %>% dplyr::summarize(TAG_COUNT_PER_PTN = n())

wesclonalitymaf0 = inner_join(wesclonalitymaf, wes_vartag_table, by = c('PTN','VAR_TAG')) %>%
                   mutate(MUT_STATUS = ifelse(TAG_COUNT_PER_PTN >1, 'SHARED', 'UNIQUE'))

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
  Fraction_Shared = Shared/UNION,
  Fraction_Private_Primary = Private_Primary/UNION, 
  Fraction_Private_Metastasis = Private_Metastasis/UNION) %>%
  select(PTN,Fraction_Shared,Fraction_Private_Primary,Fraction_Private_Metastasis)

wessubclonality_fracs

dim(wessubclonality_fracs)

wessubclonality_fracs_ = wessubclonality_fracs %>% select(names(wessubclonality_fracs)[1],
                                                          names(wessubclonality_fracs)[2:4]) %>% replace(., is.na(.), 0)
wessubclonality_fracs_.m = melt(wessubclonality_fracs_)
wessubclonality_fracs_.m$variable = factor(wessubclonality_fracs_.m$variable,
                                           levels=rev(c("Fraction_Shared","Fraction_Private_Primary","Fraction_Private_Metastasis")))

## set the levels in order we want
wessubclonality_fracs_.m$PTN <- factor(wessubclonality_fracs_.m$PTN,levels = patient.order$V1)

p1 = ggplot(data=wessubclonality_fracs_.m) +
  geom_bar(aes(x=reorder(PTN,variable), y=value, fill=wessubclonality_fracs_.m$variable),stat='identity') +
  scale_fill_jco() + 
  #scale_x_discrete(limits=patient.order$V1) +
  #values = rev(c('tan4','tan','navy','lightblue','black','grey')) )+ #top to bottom
  ylab('Fraction of mutations') + xlab('') +
  theme_classic(base_size = 16) +
  theme(axis.text.x = element_blank(),
    #axis.text.x = element_text(angle=90, hjust=1, size=NULL),
        legend.key.size = unit(.5, "cm"),
        legend.text = element_text(size = 12, color = "gray10"),
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


