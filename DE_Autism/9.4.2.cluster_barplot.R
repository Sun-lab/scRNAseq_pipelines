setwd("/Volumes/SpecialData/fh_data/Data_PRJNA434002/")

cluster_tag_seq=1:17
file_tag_seq=c("PFC3k","PFC5k") #
dist_method_seq=c("klmean","jsd")
fit_method_seq=c("empirical","nbzinb")
F_method_seq=c("p","ps")
fit_tag_seq=c("","nb") #"","nb","zinb"
resid_flag_seq=c("", "logresid","adj")
covariate_flag_seq=c("", "readdepth")

pre_tag="dca" #c("dca","scvi")

perm_label_seq=0:10
ind_covariate_flag="ind"

perm_method=""

method_seq=c("DESeq","MAST","jsd_empirical","klmean_empirical","jsd_zinb","klmean_zinb","jsd_direct","klmean_direct")

library("RColorBrewer")
cur_col=brewer.pal(9,"Set1")
#get some data info
#get meta info
#get meta info
#input phenotype #note!!!Here we have to first deal with PFC
if(length(grep("PFC",file_tag))>0){
  if(is.na(unlist(strsplit(file_tag,"k"))[2])){
    tmeta=readRDS(paste0("../",dataset_folder,"/meta_PFC.rds"))
  }
  if(!is.na(unlist(strsplit(file_tag,"k"))[2])){
    tmeta=readRDS(paste0("../",dataset_folder,"/meta",unlist(strsplit(file_tag,"k"))[2],"_PFC.rds"))
  }
}
if(length(grep("PFC",file_tag))==0){
  if(is.na(unlist(strsplit(file_tag,"k"))[2])){
    tmeta=read.table(paste0("../",dataset_folder,"/meta.tsv"),header = TRUE, sep = "\t")
  }
  if(!is.na(unlist(strsplit(file_tag,"k"))[2])){
    tmeta=readRDS(paste0("../",dataset_folder,"/meta",unlist(strsplit(file_tag,"k"))[2],".rds"))
  }
}
cur_cluster=unique(tmeta$cluster)
#get gene name
sim_matrix_bulk_dca_sim_12_PFC5k=readRDS("7.Result/rawcount_matrix_bulk_12_PFC5k.rds")
gene_symbol_PFC5k=rownames(sim_matrix_bulk_dca_sim_12_PFC5k)


################### read in files and plot preparation ################

power_array=readRDS(paste0("./8.Result/final_power_array.rds"))
i_file=1 #PFC3k
i_F=1 #ps
i_fit=2 #zinb
#i_cluster=2 #2,9,12,14
i_resid=2 #logresid
i_cov=2 #cov
#[i_file,i_F,i_fit,i_cluster,i_resid,i_cov,i_perm_label,]
#i_cluster_seq=c(2,9,12,14)
#i_method_seq=c(1,2,7)

power_cur=power_array[i_file,i_F,i_fit,,i_resid,i_cov,1,]
rownames(power_cur)=cur_cluster
colnames(power_cur)=c("DESeq2","MAST","IDEAS")
#power_cur=power_cur[i_cluster_seq,]
power_cur=power_cur[order(power_cur[,3]),]

#png("~/Desktop/fh/1.Testing_scRNAseq/Power_of_different_cell_types.png",height = 6,width=8,units="in",res=300)
op=par(mar=c(5,8,3,2))
barplot(t(power_cur),beside = TRUE,col=rep(cur_col[1:8],17),horiz=TRUE,las=2,main="Power of Different Cell Types",xlab="Proportion of Genes with p-values <0.05" )
legend("bottomright",colnames(power_cur),pch=15,col=cur_col[1:8],bty="n")
par(op)
#dev.off()

