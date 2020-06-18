######################################
#### Signatures plot ####
######################################

d = wessamplelevel %>% filter(Tumor_Sample_Barcode %in% Private(wesclonalitymaf0_fil$Tumor_Sample_Barcode)) %>%
  select(Tumor_Sample_Barcode, Number_of_Mutations, c(Signature_1:Signature_30)) %>% arrange(Tumor_Sample_Barcode)
dim(d); head(d)
min.perc=0.2; my.x.axis.text.size=10; leg.pos='right'
d = as.data.table(d)
setnames(d, names(d), c('sample_name','num_of_mutations',1:30))
setkey(d, sample_name)
e = melt(d[,c(1,3:32),with=F])
setnames(e, names(e), c('sample_name', 'signature', 'prop'))
e = as.data.table(e)
e = e[which(e$prop >= min.perc),]

other = e[,1-sum(prop), by=sample_name]
setnames(other, names(other),c('sample_name','prop'))
other = cbind(other, signature=31)
other = other[,c(1,3,2),with=F]
other$prop[which(other$prop < 0)] = 0

e = rbind(e,other,fill=T)

num.muts = d[,.(sample_name,num_of_mutations)]
setkey(num.muts)
setkey(e)

e = merge(e, num.muts, all.y=T)
e[is.na(signature),c('signature','prop'):=list(31,1)]
mmr_msi = colorRampPalette(c('#6BAED6', '#4292C6','#2171B5'))(5) # 6, 15, 20, 21, 26
aid_apobec = rev(brewer.pal(8, 'Oranges'))[2:4] # 2, 9, 13
t_c = rev(brewer.pal(8, 'Greens'))[2:4] # 5, 12, 16
c_t = rev(brewer.pal(8, 'Purples'))[4:5] # 11, 23
t_a = c("lightpink1","lightpink2") #rev(brewer.pal(8, 'Reds'))[4:5] # 25, 27
s = c(brewer.pal(12,'Set3')[c(1:3,7:12)], brewer.pal(8, 'Dark2')) # 1, 4, 7, 8, 14, 17, 18, 19, 22, 24, 28, 29, 30

my.cols = c(s[1], aid_apobec[1], s[2], 'wheat', t_c[1], mmr_msi[1], rev(brewer.pal(8, 'Reds'))[4], s[5], aid_apobec[2], s[4],  #10
            c_t[1], s[7], aid_apobec[3], s[8], mmr_msi[2], t_c[2], s[9], s[10], s[11], mmr_msi[3],         #20
            mmr_msi[4], s[12], c_t[2], s[13], t_a[1], mmr_msi[5], t_a[2], s[14], s[15], s[16], 'grey85')           #31

e$signature = sub('Signature.','',e$signature)
my.order = c(6,15,20,21,26,10,1,2,13,3,4,5,12,16,7,8,9,11,23,14,17,18,19,22,24,25,27,28,29,30,31)
my.order = rev(my.order)
my.order = my.order[which(my.order %in% e$signature)]

percent.aidapobec = e[signature %in% c(2,13),sum(prop),by=sample_name] 
sample.order = percent.aidapobec[order(percent.aidapobec$V1, decreasing=T,na.last =T),]$sample_name
no.aidapobec = setdiff(e$sample_name,sample.order)
sample.order = c(sample.order,no.aidapobec)

percent.mmrmsi = e[signature %in% c(6,15,20,21,26),sum(prop),by=sample_name]
sample.order = percent.mmrmsi[order(percent.mmrmsi$V1, decreasing=T,na.last =T),]$sample_name
no.mmrmsi = setdiff(e$sample_name,sample.order)

sample.order2 = e[sample_name %in% no.mmrmsi][order(prop,decreasing=T)]$sample_name
sample.order = c(sample.order,sample.order2)

e[,sample_name:=factor(sample_name, levels=Private(sample.order))]
e[,signature:=factor(signature,  levels=Private(my.order), ordered=T)]
setorder(e,sample_name,signature)

#print(e) 

num.muts$sample_name = factor(num.muts$sample_name, Private(sample.order))
setorder(num.muts,sample_name) 
#tips2$day <- factor(tips2$day,levels = c("Fri", "Sat", "Sun", "Thur"))
somePDFPath = "/Volumes/home/all_data/bladder-kdm6a-hannah/wes-signatures-clean-updated.pdf"
pdf(file=somePDFPath) 
ggplot(data=e) +
  geom_bar(aes(x=factor(sample_name, Private(sort(sample.order))), y=prop,fill=signature), stat='identity') +
  geom_text(aes(label='', x=sample_name, group=signature, y=prop, angle=90), size=3, position=position_fill(vjust=.5)) +
  geom_text(data=num.muts,aes(y=1,x=sample_name,label=num_of_mutations),size=2) +
  scale_fill_manual(values=my.cols[my.order], name='Signatures') +
  labs(x='', y='Proportion of mutations') +
  theme_classic() +
  theme(axis.text.x = element_text(angle=90, hjust=1, size=NULL),
        legend.key.size = unit(.5, "cm"),
        legend.text = element_text(size=6),
        legend.position=leg.pos,
        panel.grid.minor.y=element_line(colour='black', size=0),
        panel.grid.major.y=element_line(colour='black', size=0),
        panel.grid.minor.x=element_line(colour='black', size=0),
        panel.grid.major.x=element_line(colour='black', size=0),
        plot.margin = unit(c(1,1,1,1), 'lines')) +
  guides(fill = guide_legend(ncol = 1, reverse=T))
dev.off()
