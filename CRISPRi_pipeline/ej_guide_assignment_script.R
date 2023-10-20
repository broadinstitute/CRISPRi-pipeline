
#in this version pooled all reps and lanes together
args = commandArgs(trailingOnly=TRUE)

#untar things upstream?

sample = args[1]
DATADIR = args[2]
min_thresh = args[3] #default should be 2, unless I gess for crop-seq, need to ask?
use_10x_calls = args[4] #yes/no --> if no, need a min_thresh
testing = args[5]
#testing = "no" #set auto to no later

#sample = "moi5_freshB"
#DATADIR = "/Users/ejagoda/Documents/HGRM/Crispri/230901_WTC11_ENCODE_DC_TAP_MOI5_FreshvsFrozen/Outputs_cellranger/Fresh-MOI5B/"
#min_thresh = 2
#use_10x_calls = "yes"
#testing = "yes" #only flag for testing to do first 100 guides
#source("/Users/ejagoda/Documents/HGRM/Crispri/Crispr_enhancer_analysis_tap_perturb_etc/v3/fromHelen_DC/seurat_functions.R") 

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(dplyr)
  library(ggplot2)
  library(reshape2)
  library(tidyr)
  library(stringr)
  library(gridExtra)
  library(grid)
})

source("/oak/stanford/groups/engreitz/Users/ejagoda/seurat_functions.R") 

output_path = paste0(DATADIR,"/outputs/")
if (file.exists(output_path) == F){
  dir.create(output_path)
}

#if tar on computer
#filtered_feature_bc_matrix_path = paste0(DATADIR,"/filtered_feature_bc_matrix/")
#crisprdir =  paste0(DATADIR,"/crispr_analysis/")

#if tar on server

filtered_feature_bc_matrix_path = paste0(DATADIR)
crisprdir =  paste0(DATADIR)


##guide assignment
if (exists("min_thresh") == F){
  min_thresh = 2
}

if (exists("use_10x_calls") == F & file.exists(crisprdir)){
  use_10x_calls = "yes"
}

if (exists("use_10x_calls") == F & file.exists(crisprdir) == F){
  use_10x_calls = "no"
}


cat(paste0("found crispr directory (", crisprdir, "): ", file.exists(crisprdir),"\nusing thresholds from 10x cell ranger file: ",use_10x_calls,"\nmin threshold: ",min_thresh,"\noutput dir: ",output_path,"\n"))

if (exists("use_10x_calls") == "yes" & file.exists(crisprdir) == F){
  cat ("error: can't use 10x calls because crisprdir not found")
  break
}


mtx <- Read10X(filtered_feature_bc_matrix_path)
if(class(mtx) == "list") {
  message("Loading Gene Expression Matrix and Guide Capture Matrix")
  gex.mtx <- mtx[["Gene Expression"]]
  guide.mtx <- mtx[["CRISPR Guide Capture"]]
} else {
  message("loading Gene Expression Matrix")
  gex.mtx <- mtx
}

if (testing == 'yes'){
  guide.mtx_cp = guide.mtx[1:100,]
} else{
  guide.mtx_cp = guide.mtx
}


if (use_10x_calls == "yes"){
thresh_file = read.table(paste0(crisprdir,"protospacer_umi_thresholds.csv"),header=T,sep = ',')

cells_per_guide = c()
thresh_uses = c()

for (i in 1:nrow((guide.mtx_cp))){
  #cat(paste0(i,"\t,"))
  guide = row.names(guide.mtx)[i]
  if (guide %in% thresh_file$Protospacer){
    thresh = thresh_file[thresh_file$Protospacer == guide,"UMI.threshold"]
    if (thresh < min_thresh){
      thresh_use = min_thresh
    }
    else{
      thresh_use = thresh
    }
  }else{
    thresh_use = min_thresh
  }
  thresh_uses = c(thresh_uses,thresh_use)
  binarize_calls = function(value){
    if (value < thresh_use){
      return(0)
    }
    if (value >= thresh_use){
      return(1)
    }
  }
  new_row = sapply(guide.mtx[i,],binarize_calls)
  guide.mtx_cp[i,] = new_row
  cells_per_guide = c(cells_per_guide,length(which(guide.mtx_cp[i,]== 1)))
}
}


if (use_10x_calls == "no"){
  thresh_use = min_thresh
  cells_per_guide = c()
  thresh_uses = c()
  for (i in 1:nrow((guide.mtx_cp))){
    #cat(paste0(i,"\t,"))
    guide = row.names(guide.mtx)[i]
    binarize_calls = function(value){
      if (value < thresh_use){
        return(0)
      }
      if (value >= thresh_use){
        return(1)
      }
    }
    new_row = sapply(guide.mtx[i,],binarize_calls)
    guide.mtx_cp[i,] = new_row
    cells_per_guide = c(cells_per_guide,length(which(guide.mtx_cp[i,]== 1)))
    thresh_uses = c(thresh_uses,thresh_use)
  }
}

pertub_status = data.frame(guide.mtx_cp)
pertub_status = data.frame(cbind(row.names(pertub_status),pertub_status))
colnames(pertub_status)[1] = "VECTOR"
row.names(pertub_status) = NULL
#colnames(pertub_status) = gsub(colnames(pertub_status),pattern = "\\.",replacement = "-") it's always gonna write it with a .
write.table(pertub_status, paste0(output_path,sample,"calls_based_on_10x_",use_10x_calls,"min_thresh_",min_thresh,"_perturb_status.txt"),quote = F,sep = '\t',row.names  = F)


per_guide_summary_table = data.frame(cbind(row.names(guide.mtx_cp),cells_per_guide,thresh_uses))
colnames(per_guide_summary_table) = c("guide","cells","UMI.threshold")
write.table(per_guide_summary_table,paste0(output_path,sample,"calls_based_on_10x_",use_10x_calls,"min_thresh_",min_thresh,"_per_guide_summary_table.txt"),quote = F,sep = '\t',row.names  = F)

png(paste0(output_path,sample,"calls_based_on_10x_",use_10x_calls,"min_thresh",min_thresh,"_Guide_umi_threshold_hist_log10.png"))
p = ggplot(per_guide_summary_table,aes(x = log10(as.numeric(paste0(UMI.threshold)))))+
  geom_histogram(bins = 50)+
  ylab("# guides")+
  xlab("log10 Guide Umi Threshold")+
  geom_vline(xintercept = log10(1))+
  labs(title = paste0(sample,"\nGuide Umi Threshold (log10)\nMean: ",round(mean(as.numeric(paste0(per_guide_summary_table$UMI.threshold))),2),
                      "\nMedian: ",round(median(as.numeric(paste0(per_guide_summary_table$UMI.threshold))),2),
                      "\nPerc 1: ",round(nrow(per_guide_summary_table[as.numeric(paste0(per_guide_summary_table$UMI.threshold)) == 1,])/nrow(per_guide_summary_table),3) * 100,"%",
                      "\nPerc min thresh (",min_thresh,"): ",round(nrow(per_guide_summary_table[as.numeric(paste0(per_guide_summary_table$UMI.threshold)) == min_thresh,])/nrow(per_guide_summary_table),3) * 100,"%",
                      "\nPerc <5: ",round(nrow(per_guide_summary_table[as.numeric(paste0(per_guide_summary_table$UMI.threshold)) < 5,])/nrow(per_guide_summary_table),3) * 100,"%" ))
print(p)
dev.off()

png(paste0(output_path,sample,"calls_based_on_10x_",use_10x_calls,"min_thresh",min_thresh,"_Cells_per_guide.png"))
p = ggplot(per_guide_summary_table,aes(x = as.numeric(paste0(cells))))+
  geom_histogram()+
  labs(title = paste0(sample,"\nCells per Guide\nMean = ",round(mean(as.numeric(paste0(per_guide_summary_table$cells))),3),
                      "\nMedian = ",median(as.numeric(paste0(per_guide_summary_table$cells))),
                      "\nPerc 0 = ", 100*(nrow(per_guide_summary_table[per_guide_summary_table$cells == 0,])/nrow(per_guide_summary_table)),"%"))+
  xlab("Cells")+
  ylab("Guides")
print(p)
dev.off()

col_sums = colSums(guide.mtx_cp)
guides_per_cell = data.frame(col_sums)
colnames(guides_per_cell) = "total_guides"
guides_per_cell$CBC = row.names(guides_per_cell)
row.names(guides_per_cell) = NULL

write.table(guides_per_cell,paste0(output_path,sample,"calls_based_on_10x_",use_10x_calls,"min_thresh",min_thresh,"_guides_per_cell.txt"),quote = F,sep = '\t')

png(paste0(output_path,sample,"_Guides_per_cell.png"))
p = ggplot(guides_per_cell,aes(x = total_guides))+
  geom_histogram()+
  labs(title = paste0(sample,"\nGuides per cell\nMean = ",round(mean(as.numeric(paste0(guides_per_cell$total_guides))),3),
         "\nMedian = ",median(as.numeric(paste0(guides_per_cell$total_guides))),
         "\nPerc 0 = ", 100*round(length(guides_per_cell[as.numeric(paste0(guides_per_cell$total_guides)) == 0,])/nrow(guides_per_cell),3),"%"))
print(p)
dev.off()



##end guide assignment
#testing trying to get umis per ambient vs assigned

get_assigned_guides = function(cell){
  guides = rownames(guide.mtx_cp)[which(guide.mtx_cp[,cell] == 1)]
  return(paste0(guides,collapse = "|"))
}
  

assigned_guides = sapply(colnames(guide.mtx_cp),get_assigned_guides)
guide_calls_tab = data.frame(cbind(colnames(guide.mtx_cp),assigned_guides))
row.names(guide_calls_tab) = NULL
colnames(guide_calls_tab) = c("CBC","assigned_guides")

protospacer_calls_per_cell_tab = merge(guide_calls_tab,guides_per_cell)
colnames(protospacer_calls_per_cell_tab) = c("cell_barcode","feature_call","num_features")
write.table(protospacer_calls_per_cell_tab[,c(1,3,2)],paste0(output_path,sample,"protospacer_calls_per_cell.csv"),row.names = F,quote = F,sep = ",")

get_umis_assigned = function(cell){
  guides = rownames(guide.mtx_cp)[which(guide.mtx_cp[,cell] == 1)]
  assigned_umis = sum(guide.mtx[guides,cell])
  return(assigned_umis)
}

get_umis_not_assigned = function(cell){
  guides = rownames(guide.mtx_cp)[which(guide.mtx_cp[,cell] == 0)]
  not_assigned_umis = sum(guide.mtx[guides,cell])
  return(not_assigned_umis)
}

assigned_umis_list = sapply(colnames(guide.mtx_cp),get_umis_assigned)
not_assigned_umis_list = sapply(colnames(guide.mtx_cp),get_umis_not_assigned)

umis_tab = data.frame(cbind(colnames(guide.mtx_cp),assigned_umis_list,not_assigned_umis_list))
row.names(umis_tab) = NULL
colnames(umis_tab)[1] = "CBC"
umis_tab$ratio = as.numeric(paste0(umis_tab$assigned_umis_list))/as.numeric(paste0(umis_tab$not_assigned_umis_list))

png(paste0(output_path,sample,"_assinged_vs_not_assinged_umis_cdf.png"))
p = ggplot(umis_tab,aes(x = log10(ratio)))+
  stat_ecdf()+
  xlab("Assigned guide umis/non-assigned guide umis")+
  ylab("Fraction of cells")+
  labs(title = paste0(sample," Assinged vs ambient guide umis(log10) \nFractions cells ratio < 1: ",round(nrow(umis_tab[umis_tab$ratio < 1,])/nrow(umis_tab),3),"\nmean (cells with at least 1 umi): ", mean(umis_tab[umis_tab$not_assigned_umis_list != "0" | umis_tab$assigned_umis_list != "0","ratio"]),
                      "\nmedian (cells with at least 1 umi) : ",round(median(umis_tab[umis_tab$not_assigned_umis_list != "0" | umis_tab$assigned_umis_list != "0","ratio"]),3),
                      "\ncells with no umis: ",nrow(umis_tab[umis_tab$not_assigned_umis_list == "0" & umis_tab$assigned_umis_list == "0",])))
print(p)
dev.off()

png(paste0(output_path,sample,"_assinged_vs_not_assinged_umis_boxplot.png"))
p = ggplot(umis_tab,aes(y = log10(ratio)))+
  geom_boxplot()+
  ylab("log 10 Assigned guide umis/non-assigned guide umis")+
  labs(title = paste0(sample," Assinged vs ambient guide umis(log10) \nFractions cells ratio < 1: ",round(nrow(umis_tab[umis_tab$ratio < 1,])/nrow(umis_tab),3),"\nmean (cells with at least 1 umi): ", mean(umis_tab[umis_tab$not_assigned_umis_list != "0" | umis_tab$assigned_umis_list != "0","ratio"]),
                      "\nmedian (cells with at least 1 umi) : ",round(median(umis_tab[umis_tab$not_assigned_umis_list != "0" | umis_tab$assigned_umis_list != "0","ratio"]),3),
                      "\ncells with no umis: ",nrow(umis_tab[umis_tab$not_assigned_umis_list == "0" & umis_tab$assigned_umis_list == "0",])))
print(p)
dev.off()
##make graphs of this etc this is really good good good
