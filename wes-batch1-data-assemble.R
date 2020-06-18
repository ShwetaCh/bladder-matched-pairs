######################################
## WES TEMPO results 11/26/19 Figure-4
######################################

wessamplelevel = fread('/Users/chavans/juno/work/ccs/ccs_wes/Proj_07871_DFLOQ/Result/somatic/sample_data.txt')
wessamplelevel = wessamplelevel %>% mutate(Tumor_Sample_Barcode = str_split(sample,"__",simplify = TRUE)[,1],
                                           SAMPLE_TYPE = ifelse(grepl("M0",sample)=="TRUE","Metastasis","Primary"),
                                           PTN = substr(Tumor_Sample_Barcode,1,10))
dim(wessamplelevel); head(wessamplelevel) #sample

NAsamples = filter(wessamplelevel, is.na(purity)) %>% select(Tumor_Sample_Barcode)
NAsamples$Tumor_Sample_Barcode
[1] "s_C_271D5P_M001_d" "s_C_271D5P_P001_d" "s_C_ADEV70_P001_d"
[4] "s_C_PDMVDR_P001_d" "s_C_W384MJ_M001_d"

wessamplelevel = wessamplelevel %>% mutate(purity = ifelse(is.na(purity), -1, purity))
wessamplelevel = wessamplelevel %>% mutate(LOW_PURITY = ifelse(purity > 0.2 ,FALSE, TRUE)) 
filter(wessamplelevel, LOW_PURITY == TRUE) %>% select(PTN) %>% distinct(.)
length(Private(wessamplelevel$PTN))
table(wessamplelevel$LOW_PURITY); table(wessamplelevel$purity)

#14 are low purity

wesqc_align = fread('/Users/chavans/juno/work/ccs/ccs_wes/Proj_07871_DFLOQ/Result/qc/alignment_qc.txt')
dim(wesqc_align); head(wesqc_align) #Sample T
#filter(wesqc_align, MeanTargetCoverage <=50)
filter(wesqc_align, Sample %like% '_N0',  MedianTargetCoverage <=45) #Exclude this pair s_C_001390_N001_d	s_C_001390_M001_d, N coverage is 18.
filter(wesqc_align, !(Sample %like% '_N0'),  MedianTargetCoverage <=70) #Exclude s_C_001390_M001_d is 52.

wesqc_concor = fread('/Users/chavans/juno/work/ccs/ccs_wes/Proj_07871_DFLOQ/Result/qc/concordance_qc.txt')
dim(wesqc_concor); head(wesqc_concor) #Sample T_N
filter(wesqc_concor, Concordance <90) #s_C_001601_P001_d__s_C_001601_N002_d       92.03

wesqc_contam = fread('/Users/chavans/juno/work/ccs/ccs_wes/Proj_07871_DFLOQ/Result/qc/contamination_qc.txt')
dim(wesqc_contam); head(wesqc_contam) #Sample
filter(wesqc_contam, Contamination >= 5) #s_C_001601_P001_d__s_C_001601_N002_d           T s_C_001601_P001_d         5.518

wesmaf = fread('/Users/chavans/juno/work/ccs/ccs_wes/Proj_07871_DFLOQ/Result/somatic/mut_somatic.maf')
dim(wesmaf); head(wesmaf) #Tumor_Sample_Barcode
wesids = fread('/Users/chavans/juno/work/ccs/chavans/res/bladder_kdm6a/WES_subset_ids.txt')
dim(wesids); head(wesids) #SampleID

wessamplelevel_pure = filter(wessamplelevel, !(LOW_PURITY == TRUE), !(sample %like% 's_C_001601'), !(sample %like% 's_C_001390'))
Private(wessamplelevel_pure$PTN)

##Multiple samples
multi = wessamplelevel %>% group_by(PTN) %>% dplyr::dplyr::summarise(total = n()) %>% filter(total > 2) %>% select(PTN)
