samples = fread('/Volumes/home/all_data/sample_level.txt', header = TRUE)
names(samples)[1] = 'sampleid'
head(samples)
print(length(samples$sampleid))
project_home = '/ifs/res/taylorlab/chavans/bladder_kdm6a/'

impact_facets = fread('/ifs/res/taylorlab/impact/facets/facets_0.5.14/manifests/impact_facets_manifest_2019_09_07.txt') %>% 
  filter(tumor_sample %in% samples$sampleid) %>%
  filter(run_status == "complete") %>%
  select(tumor_sample, run_output_dir)
head(impact_facets)
dim(impact_facets) #1221
#write.table(impact_facets,'/ifs/res/taylorlab/chavans/bladder_kdm6a/facets_mapping_bl.txt', row.names = FALSE, quote = FALSE, sep = "\t", append = FALSE)

#load('/ifs/res/taylorlab/impact/facets/facets_0.5.14/manifests/impact_facets_qc_2019_09_07.Rdata')
#names(impact_facets_qc)
samples_with_atleast_one_good_fit = filter(impact_facets_qc, tumor_sample %in% 
                                             samples$sampleid, facets_suite_qc == TRUE, valid_purity_filter_pass ==TRUE) %>% 
                                             distinct(tumor_sample)
write.table(samples_with_atleast_one_good_fit,'/ifs/res/taylorlab/chavans/bladder_kdm6a/samples_with_atleast_one_good_facets_fit.txt', row.names = FALSE, quote = FALSE, sep = "\t", append = FALSE)

samples_with_no_good_fit = filter(impact_facets_qc, tumor_sample %in% 
                                             samples$sampleid, facets_suite_qc == FALSE) %>% 
                                             distinct(tumor_sample)

df = filter(impact_facets_qc, tumor_sample %in% samples$sampleid, facets_suite_qc == TRUE, valid_purity_filter_pass ==TRUE, fit_name == "default") %>% distinct(.)
dim(df); dim(df); length(Private(df$tumor_sample)) #800
impact_facets_good_df = df %>% filter(tumor_sample %in% Private(df$tumor_sample)) %>% select(tumor_sample, purity_run_prefix) %>% distinct(.)

df1 = filter(impact_facets_qc, tumor_sample %in% samples$sampleid, facets_suite_qc == TRUE, valid_purity_filter_pass ==TRUE, fit_name != "default") %>% distinct(.)
dim(df1); dim(df1); length(Private(df1$tumor_sample)) #1339
impact_facets_good_df1 = df1 %>% filter(tumor_sample %in% Private(df1$tumor_sample)) %>% select(tumor_sample, purity_run_prefix) %>% distinct(.)

df2 = setdiff(df1$tumor_sample, df$tumor_sample) #198
impact_facets_good_df2 = filter(impact_facets_good_df1, df1$tumor_sample %in% df2)

final_impact_facets = rbind(impact_facets_good_df,impact_facets_good_df2) #1157 multiple alt diplogr exists, deduplicated in excel, and updated the below file.
dim(final_impact_facets); head(final_impact_facets) 

#write.table(final_impact_facets,'/ifs/res/taylorlab/chavans/bladder_kdm6a/facets_mapping_good_bl.txt', row.names = FALSE, quote = FALSE, sep = "\t", append = FALSE)

impact_maf = fread('/ifs/res/taylorlab/data_repositories/dmp/mskimpact/data_mutations_extended.txt', skip = 1)
impact_maf_bl = impact_maf %>% filter(Tumor_Sample_Barcode %in% samples$sampleid)
dim(impact_maf_bl); head(impact_maf_bl); length(Private(impact_maf_bl$Tumor_Sample_Barcode)) #1215
write.table(impact_maf_bl,'/ifs/res/taylorlab/chavans/bladder_kdm6a/data_mutations_extended_bl.txt', row.names = FALSE, quote = FALSE, sep = "\t", append = FALSE)

#impact_clinical = fread('/ifs/res/taylorlab/data_repositories/dmp/mskimpact/data_clinical_sample.txt', skip = 4) %>%
# filter(SAMPLE_ID %in% samples$sampleid) %>%
#  select(SAMPLE_ID, SAMPLE_TYPE)
#head(impact_clinical); dim(impact_clinical)

ccf0 = fread('/ifs/res/taylorlab/chavans/bladder_kdm6a/data_mutations_extended_bl.mafanno.txt') %>% 
  mutate(VAF = as.numeric(t_alt_count/(t_alt_count+t_ref_count)))
ccf = left_join(samples, ccf0, by = c('sampleid' = 'Tumor_Sample_Barcode')) 
names(ccf); dim(ccf); length(Private(ccf$sampleid))
setnames(ccf, 'sampleid', 'Tumor_Sample_Barcode')
setnames(ccf,'SampleGrade', 'SAMPLE_TYPE')
head(ccf)
ccf_genes = ccf %>% 
  select(ccf_expected_copies_em, Hugo_Symbol, Tumor_Sample_Barcode, SAMPLE_TYPE) %>% 
  filter(Hugo_Symbol %in% c('ARID1A','KDM6A','TP53','RB1','FGFR3','CREBBP','KTM2C','PIK3CA','ERBB2'))


ccf_genes$SAMPLE_TYPE <- factor(ccf_genes$SAMPLE_TYPE, levels = c('Low','High','Metastatic'),ordered = TRUE)
p1 = ggplot(data = ccf_genes, aes(x = Hugo_Symbol, y = ccf_expected_copies_em, color = SAMPLE_TYPE)) + 
  geom_boxplot(width=0.5) + 
  theme_classic(base_size=14) + 
  scale_color_manual(values = c('steelblue','maroon','darkgrey')) +
  theme(legend.position = 'bottom', axis.title.x = element_blank(), axis.ticks.y = element_blank()) +
  ylab('Cancer_Clonal_Fraction') +
  ylim(c(0,1.15)) +
  stat_compare_means(aes(group = SAMPLE_TYPE, label = paste0("p = ", ..p.format..)))
p1
####
df2 <- data_summary(ccf_genes, varname="ccf_expected_copies_em", 
                    groupnames=c("SAMPLE_TYPE", "Hugo_Symbol"))
# Convert dose to a factor variable
df2$Hugo_Symbol=as.factor(df2$Hugo_Symbol)
head(df2)

df2$SAMPLE_TYPE <- factor(df2$SAMPLE_TYPE, levels = c('Primary','Metastasis'),ordered = TRUE)
p2 = ggplot(data = df2, aes(x = Hugo_Symbol, y = ccf_expected_copies_em, fill = SAMPLE_TYPE)) + 
  geom_bar(stat="identity", color="black", position=position_dodge(), width = 0.5) + 
  geom_errorbar(aes(ymin=ccf_expected_copies_em-sd, ymax=ccf_expected_copies_em+sd), width=.2, position=position_dodge(.4)) +
  theme_classic(base_size=14) + 
  scale_fill_manual(values = c('darkgray','black')) +
  theme(legend.position = 'bottom', axis.title.x = element_blank(), axis.ticks.y = element_blank()) +
  ylab('Cancer_Clonal_Fraction') 

grid.arrange(p1,p2, ncol = 1, heights = c(1,1), newpage = F)
ggsave('~/kdm6a_paper/Fig1B.pdf', plot = grid.arrange(p1,p2, ncol = 1, heights = c(1,1), newpage = F))
