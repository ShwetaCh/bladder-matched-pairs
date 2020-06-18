############################################################
## ccf Mut plot
############################################################
somePDFPath = "/Volumes/home/all_data/bladder-kdm6a-hannah/wes-mutccf-clean-updated.pdf"
pdf(file=somePDFPath)  
allplots = list() #*
patids = Private(wesclonalitymaf0_fil$PTN)
wesclonalitymaf0_fil = as.data.table(wesclonalitymaf0_fil)
for(i in 1:length(patids))
  local({ #*
    print(paste(patids[i], sample_purity))
    ptn_subset_maf = wesclonalitymaf0_fil[PTN==patids[i]]
    print(dim(ptn_subset_maf))
    p = do.plot(ptn_subset_maf) #*any()
    allplots[[i]] <<- p #*
    print(p) #*
  })
dev.off()

##############################################################
# ccf CNA plot
##############################################################
#Per patient plots

acn_list = list()
load("~/juno/work/ccs/ccs_wes/Proj_07871_DFLOQ/Result/somatic/facets/s_C_002092_T001_d__s_C_002092_N002_d/facets0.5.14c100pc500/s_C_002092_T001_d__s_C_002092_N002_d_purity.Rdata")
acn_list[[1]] = allele.copy.number(out = out, fit = fit)
load("~/juno/work/ccs/ccs_wes/Proj_07871_DFLOQ/Result/somatic/facets/s_C_002092_M002_d__s_C_002092_N001_d/facets0.5.14c100pc500/s_C_002092_M002_d__s_C_002092_N001_d_purity.Rdata")
acn_list[[2]] = allele.copy.number(out = out, fit = fit)
print(acn_list)
grid.arrange(acn_list[[2]], acn_list[[2]], nrow = 2)
grid.arrange(plotlist = acn_list, nrow = 2)

#cowplot::plot_grid(plotlist = acn_list,nrow = 2)
