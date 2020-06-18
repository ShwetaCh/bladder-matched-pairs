
######################################
#### Clonality plot ####
######################################
#~/Rscript3.5.0 ~/res/scripts-orphan/get_ccf_subclonality_all_batches_expected_copies.R current_wes_somatic50.maf current_wes_somatic50.ccf.maf
######
##Clinical track
cli_info<-fread("~/bladder-paper/wes_clinical.csv")
status_order<-c("met.location","type.treatment","time.between")
cli_info$key=factor(cli_info$key,level=c(unique(cli_info$key)))
cli_info$value=factor(cli_info$value,level=c(unique(cli_info$value)))
cli_info$patient.id<-factor(cli_info$patient.id)
patient.order<-levels(unique(cli_info$patient.id)); patient.order

color_map<- c("#D9D9D9","#969696", "#FFF7EC","#FDD49E","#FC8D59","#D7301F","#7F0000","#DEEBF7" ,"#9ECAE1", "#4292C6", "#08519C")
track<- ggplot(cli_info, aes(x=patient.id,y=key))+ 
  geom_tile(aes(fill=value),width=0.9,height=0.9)+
  scale_y_discrete(limits=rev(status_order))+
  scale_x_discrete(limits=patient.order) +
  scale_fill_manual(values=color_map,name = "Clinical Information",guide="legend")+ 
  guides(color = guide_legend(ncol = 1)) +
  theme(legend.key = element_rect(size = 2, color = "white"),
        legend.key.size = unit(1.5, "lines"))+
  theme(text = element_text(color = "gray20"),
        legend.position = c("right"), # position the legend in the upper left
        legend.justification = 0, # anchor point for legend.position.
        legend.text = element_text(size = 12, color = "gray10"),
        title = element_text(size = 12, face = "bold", color = "gray10"),
        axis.text = element_text(size =12, face = "bold",angle=90),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.background = element_blank()
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
wesclonalitymaf0_fil = wesclonalitymaf0_fil %>% mutate(Tumor_Sample_Barcode = mapvalues(Tumor_Sample_Barcode, sampleqc$TUMOR_ID, sampleqc$DMP_PatID))
wesclonalitymaf0_fil %>% select(Tumor_Sample_Barcode) %>% distinct()



