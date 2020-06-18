######################################
#### TMB plot ####
######################################

sampleqc = fread("/Volumes/home/all_data/bladder-kdm6a-hannah/SampleListPreAndPostQC-050520.txt") %>% filter(Exclude == 0)
dim(sampleqc) #51
sampleqc_pass = filter(sampleqc, TUMOR_ID != "s_C_FR895T_P003_d") # Only 1 common mutation wuth P001
#filter(PTN !="s_C_ADEV70") %>% #FACETS not available for the primary sample
unique(sampleqc_pass$TUMOR_ID) #50

tmbdata = fread('/Users/chavans/juno/work/ccs/ccs_wes/Proj_07871_DFLOQ/Result/somatic/sample_data.txt') %>% 
  mutate(Tumor_Sample_Barcode = substring(sample,1,17),
         PTN = substring(Tumor_Sample_Barcode,1,10),
  SAMPLE_TYPE = ifelse(substring(Tumor_Sample_Barcode,11,13)=="_M0","Metastasis","Primary")) %>%
  filter(Tumor_Sample_Barcode != "s_C_FR895T_P003_d") %>% 
  select(SAMPLE_TYPE, TMB, PTN)
tmbdata$SAMPLE_TYPE <- factor(tmbdata$SAMPLE_TYPE, levels = c('Primary','Metastasis'),ordered = TRUE)
unique(tmbdata$Tumor_Sample_Barcode)

#wessamplelevel_pure_ = wessamplelevel_pure %>% select(Tumor_Sample_Barcode, SAMPLE_TYPE, TMB, LOW_PURITY, PTN)
# exclude_singletons = wessamplelevel_pure_ %>% group_by(PTN) %>% dplyr::dplyr::summarise(total = n()) %>% filter(total < 2) %>% select(PTN)
# wessamplelevel_pure__ = filter(wessamplelevel_pure_, PTN %nin% exclude_singletons$PTN) %>% select(SAMPLE_TYPE, TMB, PTN)
# length(unique(wessamplelevel_pure__$PTN)) #11

#intersect(wessamplelevel_pure_ %>% group_by(PTN) %>% dplyr::dplyr::summarise(total = n()) %>% filter(total < 2) %>% select(PTN),filter(wessamplelevel, LOW_PURITY == TRUE) %>% select(PTN))
pri = filter(tmbdata, SAMPLE_TYPE == "Primary") %>% arrange(TMB)
met = filter(tmbdata, SAMPLE_TYPE == "Metastasis") %>% arrange(TMB)
pri_met = left_join(pri, met, by = 'PTN')
setnames(pri_met, c('TMB.x', 'TMB.y'), c('Primary','Metastasis'))
somePDFPath = "/Volumes/home/all_data/bladder-kdm6a-hannah/wes-tmb-clean-updated_060220.pdf"
pdf(file=somePDFPath)  
#theme_classic(base_size = 16)
ggpaired(pri_met, cond1 = "Primary", cond2 = "Metastasis", color = "condition",
         line.color = "gray", line.size = 0.4) +
  scale_color_nejm() + ylab("Tumor Mutation Burden (Mutations Per MB)") + xlab("") + 
  scale_y_discrete(limits=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18")) +
  stat_compare_means(paired = TRUE, method = 'wilcox.test', method.args = list(alternative = "less")) + 
  theme(legend.position = 'bottom', legend.title = element_blank())  
dev.off()

BLS = fread('/Volumes/home/all_data/bladder-kdm6a-hannah/bladder_impact_set.txt')

imtmbdata = fread('/Volumes/home/all_data/bladder-kdm6a-hannah/impact_patient_tmb.txt') %>% 
  select(patient.id = `Patient ID`, SampleID = `Sample ID`, type = `Sample Type`, TMB.old = `Impact TMB Score`, Panel = `Gene Panel`,	MutCount = `Mutation Count`) %>% 
  filter(patient.id %in% BLS$PatientID1)
length(unique(imtmbdata$patient.id)) #80
dim(imtmbdata) #162


canonical_capture_area = c('v3'= 0.896665, 'v5' =1.016478, 'v6' = 1.139322)
im_tmb <- imtmbdata %>%
  group_by(SampleID) %>% 
  mutate(TMB = case_when(
    grepl("341", Panel) ~ MutCount/canonical_capture_area[['v3']], 
    grepl("410", Panel) ~ MutCount/canonical_capture_area[['v5']],
    grepl("468", Panel) ~ MutCount/canonical_capture_area[['v6']])
  ) %>%
  ungroup()

priim = im_tmb %>% filter(type == "Primary")
metim = im_tmb %>% filter(type == "Metastasis")
 
pri_met_im = left_join(priim, metim, by = 'patient.id')
dim(pri_met_im) #83

setnames(pri_met_im, c('TMB.x','TMB.y'),c('Primary','Metastasis'))
pri_met_im = filter(pri_met_im, !(is.na(Primary) | is.na(Metastasis)))
dim(pri_met_im) #81
#write.table(pri_met_im,'/Volumes/home/all_data/bladder-kdm6a-hannah/impact-tmb_061620.txt',row.names=F,quote=F,append=F,sep="\t")

somePDFPath = "/Volumes/home/all_data/bladder-kdm6a-hannah/impact-tmb-clean-updated_061620.pdf"
pdf(file=somePDFPath) 
ggpaired(pri_met_im, cond1 = "Primary", cond2 = "Metastasis", color = "condition",
         line.color = "gray", line.size = 0.4) +
  scale_color_nejm() + ylab("Tumor Mutation Burden (Mutations Per MB)") + xlab("") + 
  #scale_y_discrete(limits=c("5","10","20","30","40","50","60","70","80","90","100")) +
  stat_compare_means(paired = TRUE, method = 'wilcox.test', method.args = list(alternative = "less")) + 
  theme(legend.position = 'bottom', legend.title = element_blank())  
dev.off()

setdiff(BLS$PatientID,imtmbdata$patient.id)
