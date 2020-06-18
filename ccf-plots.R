#### CCF plots
#filter(sampleqc_pass, substring(DMP,1,9) %in% c('P-0012205', 'P-0048306', 'P-0046580')) %>% select(TUMOR_ID,DMP)
# s_C_002092_M002_d, KDM6A-12
# s_C_FXHF72_M001_d, s_C_FXHF72_P001_d
# s_C_LF6FEV_M001_d, s_C_LF6FEV_P001_d
# 
# 1 s_C_002092_M002_d P-0012205-T01-WES
# 2 s_C_FXHF72_M001_d P-0048306-T01-WES
# 3 s_C_FXHF72_P001_d P-0048306-T02-WES
# 4 s_C_LF6FEV_M001_d P-0046580-T02-WES
# 5 s_C_LF6FEV_P001_d P-0046580-T01-WES

wesclonalitymaf =  fread('/ifs/res/taylorlab/chavans/bladder_kdm6a/bladder_mut_somatic.ccf.freeze.050520.maf')
dim(wesclonalitymaf)
length(unique(wesclonalitymaf$Tumor_Sample_Barcode))
names(wesclonalitymaf)

head(wesclonalitymaf); dim(wesclonalitymaf)
wesclonalitymaf = wesclonalitymaf %>% mutate(CLONALITY = ifelse(clonal_call == TRUE, TRUE, FALSE), 
                                             VAR_TAG = str_c(Chromosome,':',Start_Position,':',End_Position,':',Reference_Allele,':',Tumor_Seq_Allele2),
                                             PTN = substring(Tumor_Sample_Barcode,1,10),
                                             SAMPLE_TYPE = ifelse(substring(Tumor_Sample_Barcode,11,13)=="_M0","Metastasis","Primary"))
wes_vartag_table = wesclonalitymaf %>% distinct(.) %>% select(PTN, VAR_TAG) %>% group_by(PTN,  VAR_TAG) %>% summarize(TAG_COUNT_PER_PTN = n())

wesclonalitymaf0 = inner_join(wesclonalitymaf, wes_vartag_table, by = c('PTN','VAR_TAG')) %>%
  mutate(MUT_STATUS = ifelse(TAG_COUNT_PER_PTN >1, 'SHARED', 'UNIQUE'))

sampleqc = fread("/Volumes/home/all_data/bladder-kdm6a-hannah/SampleListPreAndPostQC-050520.txt") %>% filter(Exclude == 0)
dim(sampleqc) #51
wesclonalitymaf0_fil = filter(wesclonalitymaf0, PTN %in% sampleqc$Patient_ID, 
                              Tumor_Sample_Barcode %in% sampleqc$TUMOR_ID) %>%
  filter(PTN !="s_C_ADEV70") %>% #FACETS not available for the primary sample
  filter(Tumor_Sample_Barcode != "s_C_FR895T_P003_d") # Only 1 common mutation wuth P001
unique(wesclonalitymaf0_fil$Tumor_Sample_Barcode) #48


####### PAT1
somePDFPath = "/Volumes/home/all_data/bladder-kdm6a-hannah/CCF_pat1_061020.pdf"
pdf(file=somePDFPath) 

x1 = filter(wesclonalitymaf0_fil, Tumor_Sample_Barcode == "s_C_FXHF72_P001_d") %>% 
  select(ccf_expected_copies, PTN, VAR_TAG)
#hist(x1$ccf_expected_copies)
y1 = filter(wesclonalitymaf0_fil, Tumor_Sample_Barcode == "s_C_FXHF72_M001_d") %>% 
  select(ccf_expected_copies, PTN, VAR_TAG)
#hist(y1$ccf_expected_copies)

pat1 = full_join(x1,y1,by = c("PTN","VAR_TAG")) %>% mutate(ccf_expected_copies.x = ifelse(is.na(ccf_expected_copies.x),0,ccf_expected_copies.x),
                                                           ccf_expected_copies.y = ifelse(is.na(ccf_expected_copies.y),0,ccf_expected_copies.y))
dim(pat1)
head(pat1)
setnames(pat1, c('ccf_expected_copies.x','ccf_expected_copies.y'),c('Primary','Metastasis'))
#plot(pat1$Primary, pat1$Metastasis)

ggplot(pat1, aes(x = Primary, y = Metastasis)) + 
  stat_density2d(aes(alpha = ..density.., fill = ..level..), geom = "raster", contour = FALSE, h = c(1,1)) +
  scale_fill_gradient(low = "mistyrose", high = "firebrick4") + xlim(-0.05,1.05) + ylim(-0.05,1.05) +
  geom_point(pch = 21, size = 2) + #fill = 'lightgray'  
  #scale_fill_manual(values = c('lightgreen','lightgray')) +
  #scale_color_manual(values = c('red','blue')) +
  #geom_abline(intercept = 0, slope = 1, lty = 2) +
  #geom_smooth(method=lm, se=FALSE, lty=1) + 
  #scale_y_continuous(trans = 'log10',  limits = c(1,500), breaks = c(0,1,2,5,10,20,50,100,200,500)) + 
  #scale_x_continuous(trans='log10', limits = c(1,500), breaks = c(0,1,2,5,10,20,50,100,200,500)) + 
  coord_fixed() + #xlab('TMBIMPACT + 1') + ylab('TMBWES + 1') + 
  theme_classic(base_size = 14) + ggtitle('P-0048306') + 
  theme(legend.title = element_blank(), legend.position = 'none', plot.margin = unit(c(1,1,1,1), 'lines'))
dev.off()
########################
####### PAT2
somePDFPath = "/Volumes/home/all_data/bladder-kdm6a-hannah/CCF_pat2_061020.pdf"
pdf(file=somePDFPath) 

x1 = filter(wesclonalitymaf0_fil, Tumor_Sample_Barcode == "s_C_LF6FEV_P001_d") %>% 
  select(ccf_expected_copies, PTN, VAR_TAG)
#hist(x1$ccf_expected_copies)
y1 = filter(wesclonalitymaf0_fil, Tumor_Sample_Barcode == "s_C_LF6FEV_M001_d") %>% 
  select(ccf_expected_copies, PTN, VAR_TAG)
#hist(y1$ccf_expected_copies)

pat1 = full_join(x1,y1,by = c("PTN","VAR_TAG")) %>% mutate(ccf_expected_copies.x = ifelse(is.na(ccf_expected_copies.x),0,ccf_expected_copies.x),
                                                           ccf_expected_copies.y = ifelse(is.na(ccf_expected_copies.y),0,ccf_expected_copies.y))
dim(pat1)
head(pat1)
setnames(pat1, c('ccf_expected_copies.x','ccf_expected_copies.y'),c('Primary','Metastasis'))
#plot(pat1$Primary, pat1$Metastasis)

ggplot(pat1, aes(x = Primary, y = Metastasis)) + 
  stat_density2d(aes(alpha = ..density.., fill = ..level..), geom = "raster", contour = FALSE, h = c(1,1)) +
  scale_fill_gradient(low = "mistyrose", high = "firebrick4") + xlim(-0.05,1.05) + ylim(-0.05,1.05) +
  geom_point(pch = 21, size = 2) + #fill = 'lightgray'  
  #scale_fill_manual(values = c('lightgreen','lightgray')) +
  #scale_color_manual(values = c('red','blue')) +
  #geom_abline(intercept = 0, slope = 1, lty = 2) +
  #geom_smooth(method=lm, se=FALSE, lty=1) + 
  #scale_y_continuous(trans = 'log10',  limits = c(1,500), breaks = c(0,1,2,5,10,20,50,100,200,500)) + 
  #scale_x_continuous(trans='log10', limits = c(1,500), breaks = c(0,1,2,5,10,20,50,100,200,500)) + 
  coord_fixed() + #xlab('TMBIMPACT + 1') + ylab('TMBWES + 1') + 
  theme_classic(base_size = 14) + ggtitle('P-0046580') + 
  theme(legend.title = element_blank(), legend.position = 'none', plot.margin = unit(c(1,1,1,1), 'lines'))
dev.off()
########################
####### PAT3
somePDFPath = "/Volumes/home/all_data/bladder-kdm6a-hannah/CCF_pat3_061020.pdf"
pdf(file=somePDFPath) 

x1 = filter(wesclonalitymaf0_fil, Tumor_Sample_Barcode == "s_C_002092_T001_d") %>% 
  select(ccf_expected_copies, PTN, VAR_TAG)
#hist(x1$ccf_expected_copies)
y1 = filter(wesclonalitymaf0_fil, Tumor_Sample_Barcode == "s_C_002092_M002_d") %>% 
  select(ccf_expected_copies, PTN, VAR_TAG)
#hist(y1$ccf_expected_copies)

pat1 = full_join(x1,y1,by = c("PTN","VAR_TAG")) %>% mutate(ccf_expected_copies.x = ifelse(is.na(ccf_expected_copies.x),0,ccf_expected_copies.x),
                                                           ccf_expected_copies.y = ifelse(is.na(ccf_expected_copies.y),0,ccf_expected_copies.y))
dim(pat1)
head(pat1)
setnames(pat1, c('ccf_expected_copies.x','ccf_expected_copies.y'),c('Primary','Metastasis'))
#plot(pat1$Primary, pat1$Metastasis)

ggplot(pat1, aes(x = Primary, y = Metastasis)) + 
  stat_density2d(aes(alpha = ..density.., fill = ..level..), geom = "raster", contour = FALSE, h = c(1,1)) +
  scale_fill_gradient(low = "mistyrose", high = "firebrick4") + xlim(-0.05,1.05) + ylim(-0.05,1.05) +
  geom_point(pch = 21, size = 2) + #fill = 'lightgray'  
  #scale_fill_manual(values = c('lightgreen','lightgray')) +
  #scale_color_manual(values = c('red','blue')) +
  #geom_abline(intercept = 0, slope = 1, lty = 2) +
  #geom_smooth(method=lm, se=FALSE, lty=1) + 
  #scale_y_continuous(trans = 'log10',  limits = c(1,500), breaks = c(0,1,2,5,10,20,50,100,200,500)) + 
  #scale_x_continuous(trans='log10', limits = c(1,500), breaks = c(0,1,2,5,10,20,50,100,200,500)) + 
  coord_fixed() + #xlab('TMBIMPACT + 1') + ylab('TMBWES + 1') + 
  theme_classic(base_size = 14) + ggtitle('P-0012205') + 
  theme(legend.title = element_blank(), legend.position = 'none', plot.margin = unit(c(1,1,1,1), 'lines'))
dev.off()
########################
