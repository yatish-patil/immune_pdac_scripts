rm(list=ls())

library(pamr) # for doing PAM
### Running PAM analysis on immune specific genes nanostringPAT patient

# PAM ANalysis
  mainDir = "/Users/ypatil/Dropbox (ICR)/Yatish/Analysis/PDA/immune_pdac/new_again/"
  dataDir = "/Users/ypatil/Dropbox (ICR)/Yatish/Analysis/PDA/immune_pdac/new_again/data/from_NMF_fdr_to_PAM/"
  scrDir = "/Users/ypatil/Dropbox (ICR)/Yatish/Analysis/PDA/immune_pdac/new_again/scripts/immune_pdac_scripts/"
  
outDir <- paste0(mainDir,"/out_pam_analysis_2017_08_24/")
dir.create(outDir, showWarnings = FALSE)

dataFile = paste0(dataDir,"2017-08-24_7_subtypes_data_order_with_nanostring_PAT_samples_80_58_92_imm_sp_genes.txt")
imm_spec_genes_data  = read.delim(dataFile,header = TRUE,sep="\t",row.names = 1);dim(imm_spec_genes_data)

subtypeFile = paste0(dataDir,"2017-08-24_7_order_subtypes_file_with_nanostring_PAT_samples_80_58.txt")
subtypes = read.delim(subtypeFile,header = TRUE,sep = "\t");dim(subtypes)
subtypes = subtypes[order(subtypes$cophenetic_subtypes),]


t_imm_spec_genes_data = t(imm_spec_genes_data)[ order(match(rownames(t(imm_spec_genes_data)), subtypes$samples)), ]
with_subtypes_order_imm_spec_genes_data=t(t_imm_spec_genes_data);dim(with_subtypes_order_imm_spec_genes_data)

groupnanostringPAT = c(rep("1",length(which(subtypes$cophenetic_subtypes == 1))),
             rep("2",length(which(subtypes$cophenetic_subtypes == 2))),
             rep("3",length(which(subtypes$cophenetic_subtypes == 3))),
             rep("4",length(which(subtypes$cophenetic_subtypes == 4))),
             rep("5",length(which(subtypes$cophenetic_subtypes == 5))),
             rep("6",length(which(subtypes$cophenetic_subtypes == 6))),
             rep("7",length(which(subtypes$cophenetic_subtypes == 7)))
             )

designnanostringPAT = factor(groupnanostringPAT)

mydata <- list(x=with_subtypes_order_imm_spec_genes_data,y=designnanostringPAT,geneid=rownames(with_subtypes_order_imm_spec_genes_data),
               genenames=rownames(with_subtypes_order_imm_spec_genes_data))
set.seed(2)

mytrain <-pamr.train(mydata)    
mycv <- pamr.cv(mytrain,mydata, nfold=5)

pam_plot = paste0(outDir,Sys.Date(),"_7_subtypes_80_FDR_threshold_graph_IMM_PDAC_pam.pdf",sep="")
pdf(pam_plot,width=14,height=8)
pamr.plotcv(mycv)
dev.off()

conmatx<-pamr.confusion(mycv, threshold=0.38)
pamcen<-pcen<-pamr.listgenes(mytrain, mydata, threshold=0.38)
colnames(pamcen)<-c("genes",colnames(pamcen)[-1])
write.table(pamcen,paste0(outDir,Sys.Date(),"_IMM_PDAC_7_subtypes_80_PAM_centroids.txt"), sep="\t", quote=FALSE, row.names=FALSE)

