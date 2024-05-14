
library("Seurat")
library('harmony') 
library(ggplot2)
library(pryr)
library(future) 

library(fgsea) ## not compatible when seurat is loaded!
library(ggplot2)
library(BiocParallel)
library(org.Hs.eg.db)


options(future.globals.maxSize = 30000 * 1024^2)
plan(multiprocess) 

###SET LIBRARY PATHS TO PERSONAL DIRECTORY AND WELL DIRECTORY 
.libPaths(c( "~/R/3.4.3-openblas-0.2.18-omp-gcc5.4.0", .libPaths()))
.libPaths(c("/users/immune-rep/kvi236/INSTALLED_PROGRAMS/R_MODULES", .libPaths()))
options(repos='http://cran.ma.imperial.ac.uk/')

concat = function(v) {
	res = ""
	for (i in 1:length(v)){res = paste0(res,v[i])}
	res
}
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
                     function(x) 
                       rgb(x[1], x[2], x[3], alpha=alpha)) }


batch = "RCC_doublets"
PLOTS = ""
out_dir = "~/Google_Drive/Projects/Single cell analysis/InfiltrImmune/PROJECT_CITE-seq CC/Doublet characterisation project/Doublet analysis/RCC_DGE/"
out_dir_raw = "~/Google_Drive/Projects/Single cell analysis/InfiltrImmune/PROJECT_CITE-seq CC/Doublet characterisation project/Doublet analysis/RCC_DGE/"
type = "B_cells"
analysis = "Doublet_DGE"

###################################### STEP 1 simulate doublets, perform differential gene expression, compare to randomised labels
#load the labelled single cell objects
files = c("~/Google_Drive/Projects/Single cell analysis/InfiltrImmune/PROJECT_CITE-seq CC/Doublet characterisation project/Doublet analysis/Tissue1_with_mito_percentage.RDS", 
          "~/Google_Drive/Projects/Single cell analysis/InfiltrImmune/PROJECT_CITE-seq CC/Doublet characterisation project/Doublet analysis/Tissue2_with_mito_percentage.RDS",
          "~/Google_Drive/Projects/Single cell analysis/InfiltrImmune/PROJECT_CITE-seq CC/Doublet characterisation project/Doublet analysis/Blood_with_mito_percentage.RDS")

merge = NULL
for(f in c(1:length(files))){
  pbmc = readRDS(file=files[f])
  pbmc = UpdateSeuratObject(object = pbmc)
  if(length(merge)==0){
    merge = pbmc
  }else{merge<- merge(merge,y = pbmc)}
  print (f)
}
  
head(merge@meta.data)
t_doublet_counts = table(merge@meta.data$orig.ident, merge@meta.data$Doublet_cell_type)

### iterate through cell combinations with enough doublets per patient 
doublet_types = table(merge@meta.data$Doublet_cell_type)
doublet_types = names(doublet_types[which(doublet_types>=10)])
doublet_types = doublet_types[which(doublet_types!="NA")]
cell_type = merge@meta.data$cell_refined_annotation
cell_ids = rownames(merge@meta.data)

repeat_samples = c("01_Kidney_met_panc1_biopsy_CD45p1", "02_Kidney_met_panc1_biopsy_CD45p2")
doublet_types_filtered = names(which(apply(t (t_doublet_counts [repeat_samples,doublet_types]), 1, min)>=5))
doublet_types = doublet_types_filtered

for(ind in c(1:length(doublet_types))){
  t1 = t_doublet_counts[,doublet_types[ind]]
  samples_run = names(which(t1>=5))
  for(s in c(1:length(samples_run))){
    ### get observed doublet gene counts
    doublets_sub = cell_ids[intersect(which(merge@meta.data$Doublet_cell_type == doublet_types[ind]), which(merge@meta.data$orig.ident==samples_run[s]))]
    cell_type_substituents = gsub("."," ",strsplit(doublet_types[ind],"..",fixed = T)[[1]], fixed = T)
    singlets1 = cell_ids[intersect(which(cell_type == cell_type_substituents[1]), which(merge@meta.data$orig.ident==samples_run[s]))]
    singlets2 = cell_ids[intersect(which(cell_type == cell_type_substituents[2]), which(merge@meta.data$orig.ident==samples_run[s]))]
    counts_doublets = merge @ assays$ RNA@ counts[,doublets_sub]
    
    #### get simulated doublet gene counts
    n_simulated_doublets = 500 # can change this
    sample1 = sample(singlets1, n_simulated_doublets, replace =T)
    sample2 = sample(singlets2, n_simulated_doublets, replace =T)
    counts_simulated = merge @ assays$ RNA@ counts[,sample1] + merge @ assays$ RNA@ counts[,sample2]
    
    #### merge and normalise/scale
    types = c(rep("observed_doublets", length(doublets_sub)), rep("simulated_doublets", length(sample1)))
    counts_all_doublets = cbind(counts_doublets, counts_simulated)
    data <- CreateSeuratObject(counts_all_doublets, project = "DGE")
    data@meta.data$types = types
    data <- NormalizeData(data)
    data <- ScaleData(data)
    
    #### perform DGE
    Idents(object = data) = data@meta.data$types
    lvels = levels(x = data)
    
    #if(length(unique(test@meta.data$pat_sample))>1){
    #    pbmc.markers <- FindMarkers(test, ident.1 = lvels[2], ident.2 = lvels[1],min.pct = 0.005, logfc.threshold = 0.25,test.use = "LR", latent.vars = 'pat_sample') 
    # }else{pbmc.markers <- FindMarkers(test, ident.1 = lvels[2], ident.2 = lvels[1],min.pct = 0.005, logfc.threshold = 0.25,test.use = "LR")}
    pbmc.markers1 <- FindMarkers(data, ident.1 = lvels[1], ident.2 = lvels[2],min.pct = 0.05, logfc.threshold = 0.1,test.use = "poisson")
    length(which(as.numeric(pbmc.markers1[,"p_val_adj"])<0.05))
    
    file_prefix = gsub(" ","_", concat(c(doublet_types[ind],"__", samples_run[s])), fixed = T)
    
    head(pbmc.markers1,50)
    out_file_table = concat(c(out_dir, "Doublet_DGE_overall_", file_prefix,"_", batch,".txt"))
    write.table(pbmc.markers1, file = out_file_table, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,col.names = T, qmethod = c("escape", "double"),fileEncoding = "")
    print (file_prefix)
    length(which(as.numeric(pbmc.markers1[,"p_val_adj"])<0.05))
    
    #### randomise
    runs = 5
    overall_randomised = NULL
    data@meta.data$types_random = sample(data@meta.data$types,length(data@meta.data$types), replace = F)
    for(r in c(1:runs)){
      data@meta.data$types_random = sample(data@meta.data$types_random,length(data@meta.data$types_random), replace = F)
      print(table(head(data@meta.data$types_random, 50)))
      table(data@meta.data$types_random, data@meta.data$types)
      Idents(object = data) = data@meta.data$types_random
      lvels = levels(x = data)
      pbmc.markers <- FindMarkers(data, ident.1 = lvels[1], ident.2 = lvels[2],min.pct = 0.05, logfc.threshold = 0.1,test.use = "poisson")
      pbmc.markers = cbind(pbmc.markers, r)
      if(length(overall_randomised)==0){
        overall_randomised = pbmc.markers
      }else{overall_randomised = rbind(overall_randomised, pbmc.markers)}
    }
    length(which(as.numeric(overall_randomised[,"p_val_adj"])<0.05))/runs
    
    file_prefix = gsub(" ","_", concat(c(doublet_types[ind],"__", samples_run[s])), fixed = T)
    
    head(pbmc.markers,50)
    out_file_table = concat(c(out_dir, "Doublet_DGE_randomised_overall_", file_prefix,"_", batch,".txt"))
    write.table(overall_randomised, file = out_file_table, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,col.names = T, qmethod = c("escape", "double"),fileEncoding = "")
  }
}  


#### perform DGE combining tissue 1 and 2
for(ind in c(1:length(doublet_types))){
  t1 = sum(t_doublet_counts[repeat_samples,doublet_types[ind]])
  if(sum(t1)>=10){
    ### get observed doublet gene counts
    doublets_sub = cell_ids[intersect(which(merge@meta.data$Doublet_cell_type == doublet_types[ind]), which(merge@meta.data$orig.ident %in% repeat_samples))]
    cell_type_substituents = gsub("."," ",strsplit(doublet_types[ind],"..",fixed = T)[[1]], fixed = T)
    
    counts_simulated =NULL
    source = merge@meta.data[doublets_sub,"orig.ident"]
    for(i in c(1:length(repeat_samples))){
      singlets1 = cell_ids[intersect(which(cell_type == cell_type_substituents[1]), which(merge@meta.data$orig.ident %in% repeat_samples[i]))]
      singlets2 = cell_ids[intersect(which(cell_type == cell_type_substituents[2]), which(merge@meta.data$orig.ident %in% repeat_samples[i]))]
      counts_doublets = merge @ assays$ RNA@ counts[,doublets_sub]
    
      #### get simulated doublet gene counts
      n_simulated_doublets = 250 # can change this, subsampled per sample
      sample1 = sample(singlets1, n_simulated_doublets, replace =T)
      sample2 = sample(singlets2, n_simulated_doublets, replace =T)
      source = c(source, rep(repeat_samples[i], length(sample1)))
      if(length(counts_simulated)==0){
        counts_simulated = merge @ assays$ RNA@ counts[,sample1] + merge @ assays$ RNA@ counts[,sample2]
      }else{counts_simulated = cbind(counts_simulated, (merge @ assays$ RNA@ counts[,sample1] + merge @ assays$ RNA@ counts[,sample2]))}
    }
      
    #### merge and normalise/scale
    
    colnames(counts_simulated) = apply(cbind(colnames(counts_simulated) , ceiling(runif(length(counts_simulated[1,]))*10000)), 1, paste, collapse = "")
    table(table( colnames(counts_simulated)))
    types = c(rep("observed_doublets", length(doublets_sub)), rep("simulated_doublets", length(counts_simulated[1,])))
    counts_all_doublets = cbind(counts_doublets, counts_simulated)
    data <- CreateSeuratObject(counts_all_doublets, project = "DGE")
    data@meta.data$types = types
    data@meta.data$source = source
    data <- NormalizeData(data)
    data <- ScaleData(data)
    
    #### perform DGE
    Idents(object = data) = data@meta.data$types
    lvels = levels(x = data)
    
    pbmc.markers1 <- FindMarkers(data, ident.1 = lvels[2], ident.2 = lvels[1],min.pct = 0.005, logfc.threshold = 0.25,test.use = "poisson", latent.vars = 'source') 
    #pbmc.markers1 <- FindMarkers(data, ident.1 = lvels[1], ident.2 = lvels[2],min.pct = 0.05, logfc.threshold = 0.1,test.use = "poisson")
    length(which(as.numeric(pbmc.markers1[,"p_val_adj"])<0.05))
    
    file_prefix = gsub(" ","_", concat(c(doublet_types[ind],"__combined_tissue")), fixed = T)
    
    head(pbmc.markers1,50)
    out_file_table = concat(c(out_dir, "Doublet_DGE_overall_", file_prefix,"_", batch,".txt"))
    write.table(pbmc.markers1, file = out_file_table, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,col.names = T, qmethod = c("escape", "double"),fileEncoding = "")
    print (file_prefix)
    length(which(as.numeric(pbmc.markers1[,"p_val_adj"])<0.05))
    
    #### randomise
    runs = 5
    overall_randomised = NULL
    data@meta.data$types_random = sample(data@meta.data$types,length(data@meta.data$types), replace = F)
    for(r in c(1:runs)){
      data@meta.data$types_random = sample(data@meta.data$types_random,length(data@meta.data$types_random), replace = F)
      print(table(head(data@meta.data$types_random, 50)))
      
      table(data@meta.data$types_random, data@meta.data$types)
      Idents(object = data) = data@meta.data$types_random
      lvels = levels(x = data)
      pbmc.markers <- FindMarkers(data, ident.1 = lvels[2], ident.2 = lvels[1],min.pct = 0.005, logfc.threshold = 0.25,test.use = "poisson", latent.vars = 'source') 
      print(head(pbmc.markers, 3))
      pbmc.markers = cbind(pbmc.markers, r)
      if(length(overall_randomised)==0){
        overall_randomised = pbmc.markers
      }else{overall_randomised = rbind(overall_randomised, pbmc.markers)}
    }
    length(which(as.numeric(overall_randomised[,"p_val_adj"])<0.05))/runs
    
    file_prefix = gsub(" ","_", concat(c(doublet_types[ind],"__combined_tissue")), fixed = T)
    
    head(pbmc.markers,50)
    out_file_table = concat(c(out_dir, "Doublet_DGE_randomised_overall_", file_prefix,"_", batch,".txt"))
    write.table(overall_randomised, file = out_file_table, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,col.names = T, qmethod = c("escape", "double"),fileEncoding = "")
  }
}  

###################################### STEP 2
######## parse output: 
## which outputs are robust between samples?
## which outputs are statistically greater than randomisation?
overall_markers_obs = NULL
overall_markers_rand = NULL
samples_run1 = c(sort(unique(merge@meta.data$orig.ident)), "combined_tissue")
for(ind in c(1:length(doublet_types_filtered))){
  t1 = t_doublet_counts[,doublet_types[ind]]
  samples_run = names(which(t1>=5))
  for(s in c(1:length(samples_run1))){
    file_prefix = gsub(" ","_", concat(c(doublet_types_filtered[ind],"__", samples_run1[s])), fixed = T)
    file_obs = concat(c(out_dir, "Doublet_DGE_overall_", file_prefix,"_", batch,".txt"))
    file_ran = concat(c(out_dir, "Doublet_DGE_randomised_overall_", file_prefix,"_", batch,".txt"))
    
    p <- as.matrix(read.csv(file_obs, head=T, sep="\t"))
    #w = which(as.numeric(p[,'p_val_adj'])<0.05)
    #if(length(w)<2){w = c(1:3)}
    p=cbind(samples_run1[s],doublet_types[ind], p[,])
    
    p1 <- as.matrix(read.csv(file_ran, head=T, sep="\t"))
    w = which(as.numeric(p1[,'p_val'])<0.05)
    if(length(w)<2){w = c(1:3)}
    p1=cbind(samples_run1[s], doublet_types[ind], p1[,])
    
    if(length(overall_markers_obs)==0){
      overall_markers_obs = p
      overall_markers_rand = p1
    }else{
      overall_markers_obs = rbind(overall_markers_obs, p)
      overall_markers_rand = rbind(overall_markers_rand, p1)}
  }
}
#### compare between tissue samples
# only assess doublet types that have 5 or more doublets in both repeats. 
t_count = table(overall_markers_obs[,2], overall_markers_obs[,1])
t(t_doublet_counts[repeat_samples, doublet_types_filtered])
t_count[doublet_types_filtered,repeat_samples]
t_count = t_count[which(apply(t_count[doublet_types_filtered,repeat_samples], 1, max)>0),]


fileout1=concat(c(out_dir,"Doublet_DGE_replication_", batch,".pdf"))
w=3
pdf(file=fileout1, height=w*1, width=w*3)
par(mfrow= c(1,3), mar = c(4.5,4.5,4,4))
l1= NULL
l2= NULL
l3= NULL
l4= NULL
l5= NULL
jaccard = NULL
plot_doublet_types_filtered = gsub(".", " ", gsub("..", " - ", doublet_types_filtered, fixed = T), fixed = T)
plot_doublet_types_filtered = gsub("memory", "mem. ", gsub("Conventional ", "Conv. ", plot_doublet_types_filtered, fixed = T), fixed = T)
for(ind in c(1:length(doublet_types_filtered))){
  w1 = intersect(which(overall_markers_obs[,2]==doublet_types_filtered[ind]), which(overall_markers_obs[,1]==repeat_samples[1]))
  w2 = intersect(which(overall_markers_obs[,2]==doublet_types_filtered[ind]), which(overall_markers_obs[,1]==repeat_samples[2]))
  w3 = intersect(which(overall_markers_obs[,2]==doublet_types_filtered[ind]), which(overall_markers_obs[,1]=="combined_tissue"))
  
  all_genes = sort(unique(rownames(overall_markers_obs)))
  x = rep(1,length(all_genes))
  names(x) = all_genes
  y = x
  z = x
  
  overall_markers_obs[w1[which(as.numeric(overall_markers_obs[w1,"avg_log2FC"]) > 2.5)],]
  
  x[rownames(overall_markers_obs[w1,])] = as.numeric(overall_markers_obs[w1,"avg_log2FC"])
  y[rownames(overall_markers_obs[w2,])] = as.numeric(overall_markers_obs[w2,"avg_log2FC"])
  z[rownames(overall_markers_obs[w3,])] = as.numeric(overall_markers_obs[w3,"avg_log2FC"])
  
  x[rownames(overall_markers_obs[w1,])] = as.numeric(overall_markers_obs[w1,"p_val_adj"])
  y[rownames(overall_markers_obs[w2,])] = as.numeric(overall_markers_obs[w2,"p_val_adj"])
  z[rownames(overall_markers_obs[w3,])] = as.numeric(overall_markers_obs[w3,"p_val_adj"])
  x[which(x==0)] = 1e-100
  y[which(y==0)] = 1e-100
  z[which(z==0)] = 1e-100
  
  
  cols = rep("grey", length(all_genes))
  cols[intersect(which(x<0.05),which(y<0.05))] = "darkgreen"
  cols = add.alpha(cols, alpha = 0.8)
  
  plot(x,y, log = 'xy', bg = cols, pch = 21, main = plot_doublet_types_filtered[ind], col = NA, xlab = "p-value (repeat 1)",ylab = "p-value (repeat 2)")
  
  names(x)[intersect(which(x<0.05),which(y<0.05))]
  
  l1 = c(l1, length(intersect(which(x<0.05),which(y<0.05))))
  l2 = c(l2, length(which(x<0.05)))
  l3 = c(l3, length(which(y<0.05)))
  l4 = c(l4, length(intersect(which(x<0.05),which(z<0.05))))
  l5 = c(l5, length(intersect(which(y<0.05),which(z<0.05))))
  # Function for computing Jaccard Similarity 
  jaccard_similarity <- function(A, B) { 
    intersection = length(intersect(A, B)) 
    union = length(A) + length(B) - intersection 
    return (intersection/union) 
  } 
  
  # Jaccard Similarity between sets, A and B 
  Jaccard_Similarity <- jaccard_similarity(which(x<0.05),which(y<0.05)) 
  
  jaccard = c(jaccard, Jaccard_Similarity)
  
  all_genes[intersect(which(x<0.05),which(y<0.05))]
  all_genes[intersect(which(x<0.05),which(y<0.05))]
  all_genes[which(z<0.05)]
}
dev.off()

l = cbind(l2,l3,l1,l4,l5)
colnames(l) = c("nDGE (tissue run 1)", "nDGE (tissue run 2)", "nDGE (overlapping run 1 & 2)", "nDGE (overlapping run 1 & combined)","nDGE (overlapping run 2 & combined)")
min_overlap = 100*l[,3]/apply(l[,c(1,2)], 1, min)
l = cbind(t(t_doublet_counts[repeat_samples, doublet_types_filtered]), l, min_overlap)
colnames(l)[c(1,2,length(l[1,]))] = c("n doublets (tissue run 1)", "n doublets (tissue run 2)", "% overlap of DGEs between run 1 & 2")
rownames(l) = plot_doublet_types_filtered

out_file_table = concat(c(out_dir, "Doublet_DGE_overlap_repeats_", batch,".txt"))
write.table(l, file = out_file_table, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,col.names = T, qmethod = c("escape", "double"),fileEncoding = "")

#############################
### compare between obs versus randomised

Boxplot_custom<-function(groups, main, width_plot, colsx, colsx1, range_y, range_x){
  factors = names(groups)
  max = max(c(unlist(groups), unlist(groups), range_y)*1.15)
  min = 0
  if(min(unlist(groups))<0){
    min = min(c(unlist(groups), unlist(groups))*1.15)}
  b = (max-min)*0.034
  ylab = ""
  draw_signif_lines = TRUE
  y = max(c(unlist(groups), unlist(groups))*1)+b
  max_width = range_x
  max_scale = max#min(c(max,100))
  range = max-min
  if(range>50){scale = c(-100:100)*20}
  if(range>200){scale = c(-100:100)*50}
  if(range<=50){scale = c(-100:100)*10}
  if(range<=30){scale = c(-100:100)*5}
  if(range <10){scale = c(-100:100)*2.5}
  if(range <5){scale = c(-100:100)*1}
  if(range <4){scale = c(-100:100)*0.5}
  if(range <1.5){scale = c(-100:1000)*0.2}
  if(range <0.5){scale = c(-100:100)*0.1}
  if(range <0.1){scale = c(-100:100)*0.01}
  if(range <0.01){scale = c(-100:100)*0.001}
  cex = 0.9
  Fun<-function(x){x}
  
  scale = scale[intersect(which(scale<= max_scale), which(scale>=min))]
  plot(c(0.5, max_width +0.5),c(min, max), pch=20, col="white",xlab="",ylab ="",cex=cex, cex.lab=cex+0.1,	cex.axis=cex,cex.main = cex, col.axis = "white",tck=0, mgp = c(2,0,0), main = main, axes = FALSE, ylim = c(min, max))
  mtext(side = 2, text = ylab, line = 2.8,cex= cex-0.1, las = 3, font = 1)
  mtext(side = 1, text = factors, line = 0.35,cex= cex-0.1,  at = c(1:length(factors)), las = 2, font = 1)
  segments(0.5,Fun(scale),length(groups)+0.5,Fun(scale),col = "grey",lwd = 1,lty = 3 )
  segments(0.5,Fun(scale),length(groups)+0.5,Fun(scale),col = "grey",lwd = 1,lty = 3 )
  mtext(side = 2, text = scale, line = 0.15,cex= cex-0.1,  at =Fun(scale), las = 2, font = 1)
  width = 0.38
  index = 1
  l = length(groups)
  l1 = length(groups[[1]])
  
  for(i in c(1:l)){
    points1=as.numeric(groups[[i]])
    box1<-c(as.numeric(quantile(points1))[3], as.numeric(quantile(points1, probs = c(0.1, 0.9))), as.numeric(quantile(points1))[c(2, 4)])	
    Draw_box_plot(box1,i,width,colsx[i],1, colsx1[i])
    points(rep(i, length(points1)),points1, pch =21, col=colsx[i],bg = colsx1[i], cex = 0.7)
  }
}

Draw_box_plot<-function(box,x,width,c,lwd,line_col){
  segments(x, box[2], x, box[3], col = line_col,lwd =lwd)
  segments(x-(width/2), box[2], x+(width/2), box[2], col = line_col,lwd =lwd)
  segments(x-(width/2), box[3], x+(width/2), box[3], col = line_col,lwd =lwd)
  rect(x-width, box[4], x+width, box[5], col = c,lwd =lwd, border = line_col)
  segments(x-width, box[1], x+width, box[1], col = line_col,lwd=2*lwd)}


fileout1=concat(c(out_dir,"Doublet_DGE_obs_versus_rand_", batch,".pdf"))
w=4.5
pdf(file=fileout1, height=w*1.3, width=w*4*0.55)
par(mfrow= c(1,4), mar = c(16,3,15,0.1))

for(s in c(1:length(samples_run1))){
  groups = NULL
  groups_obs = NULL
  z_score = NULL
  for(ind in c(1:length(doublet_types_filtered))){
    w1 = intersect(which(overall_markers_obs[,2]==doublet_types_filtered[ind]), which(overall_markers_obs[,1]==samples_run1[s]))
    w2 = intersect(which(overall_markers_rand[,2]==doublet_types_filtered[ind]), which(overall_markers_rand[,1]==samples_run1[s]))
    
    head(overall_markers_obs[w1,])
    head(overall_markers_rand[w2,])
    
    l1 = length(which(as.numeric(overall_markers_obs[w1,"p_val_adj"])<0.05))
    reps = sort(unique(as.numeric(overall_markers_rand[w2,"r"])))
    l2 = NULL
    for(r in reps){
      w3 =intersect(w2, which(as.numeric(overall_markers_rand[,"r"])==reps[r]))
      l2 = c(l2, length(which(as.numeric(overall_markers_rand[w3,"p_val_adj"])<0.05)))
    }
    groups = c(groups, list(l2))
    groups_obs = c(groups_obs, l1)
    z_score = c(z_score, (l1-mean(l2))/(sd(l2)))
  }
  names(groups) = doublet_types_filtered
  names(groups_obs) = doublet_types_filtered
  names(z_score) = doublet_types_filtered
  p_value_difference = 2*pnorm(q=z_score, lower.tail=FALSE)
  
  
  plot_doublet_types_filtered = gsub(".", " ", gsub("..", " - ", doublet_types_filtered, fixed = T), fixed = T)
  plot_doublet_types_filtered = gsub("memory", "mem. ", gsub("Conventional ", "Conv. ", plot_doublet_types_filtered, fixed = T), fixed = T)
  
  range_y = c(0, max(c(unlist(groups), unlist(groups_obs))))*0.9
  range_x = length(groups)
  cols= add.alpha(c(rep("red", length(groups))), alpha = 0.5)
  names(groups) = plot_doublet_types_filtered
  Boxplot_custom(groups, samples_run1[s], 0.45, cols,cols, range_y, range_x)
  cols1 = add.alpha(c(rep("black", length(groups))), alpha = 0.5)
  points(c(1:length(groups_obs)), groups_obs, pch = 24, col = cols1, bg =cols1,cex = 2)
  leg = apply(cbind("p-val=",signif(p_value_difference, digits = 2),"\nz-score=", signif(z_score, digits = 2)), 1, paste, collapse = "")
  mtext(leg, side= 3, line = 0.2, las = 2,at = c(1:length(groups)), cex = 0.7)
  
  }

dev.off()

#############################
### consider CD4 versus CD8 DGEs in combined

fileout1=concat(c(out_dir,"Doublet_DGE_CD4_vv_CD8_", batch,".pdf"))
w=4.25
pdf(file=fileout1, height=w*1, width=w*2)
par(mfrow= c(1,2), mar = c(5,5,6,6))

for(i in c(1)){#run
  use_filtered = c(1,3)
  s = 4 ### combined DGE
  w1 = intersect(which(overall_markers_obs[,2]==doublet_types_filtered[use_filtered[1]]), which(overall_markers_obs[,1]==samples_run1[s]))
  w2 = intersect(which(overall_markers_obs[,2]==doublet_types_filtered[use_filtered[2]]), which(overall_markers_obs[,1]==samples_run1[s]))
  
  head(overall_markers_obs[w1,])
  head(overall_markers_obs[w2,])
  
  p1= overall_markers_obs[w1[which(as.numeric(overall_markers_obs[w1,"p_val_adj"])<0.05)],]
  p2= overall_markers_obs[w2[which(as.numeric(overall_markers_obs[w2,"p_val_adj"])<0.05)],]
  
  all_genes = sort(unique(rownames(overall_markers_obs)))
  x = rep(1,length(all_genes))
  names(x) = all_genes
  y = x
  
  #x[rownames(overall_markers_obs[w1,])] = as.numeric(overall_markers_obs[w1,"avg_log2FC"])
  #y[rownames(overall_markers_obs[w2,])] = as.numeric(overall_markers_obs[w2,"avg_log2FC"])
  
  x[rownames(p1)] = as.numeric(p1[,"p_val_adj"])
  y[rownames(p2)] = as.numeric(p2[,"p_val_adj"])
  x[which(x==0)] = 1e-100
  y[which(y==0)] = 1e-100
  
  cols = rep("grey", length(all_genes))
  cols[intersect(which(x<0.05),which(y< 0.05))] = "darkgreen"
    cols[intersect(which(x<0.05),which(y> 0.05))] = "darkred"
      cols[intersect(which(x>0.05),which(y< 0.05))] = "darkblue"
        cols = add.alpha(cols, alpha = 0.5)
        plot(x,y, log = 'xy', bg = cols, pch = 21, main = "p-values", col = NA, xlab = concat(c("p-value (",plot_doublet_types_filtered[use_filtered[1]],")")),
             ylab = concat(c("p-value (",plot_doublet_types_filtered[use_filtered[2]],")")))
        
        
        ######
        all_genes = sort(unique(rownames(overall_markers_obs)))
        x1 = rep(0,length(all_genes))
        names(x1) = all_genes
        y1 = x1
        
        x1[rownames(p1)] = as.numeric(p1[,"avg_log2FC"])
        y1[rownames(p2)] = as.numeric(p2[,"avg_log2FC"])
        
        #cols = rep("grey", length(all_genes))
        #cols[intersect(which(x1< -0.5),which(y1< -0.5))] = "darkgreen"
        #cols[intersect(which(x1< -0.5),which(y1> 0.5))] = "darkred"
        #cols[intersect(which(x1> 0.5),which(y1< -0.5))] = "darkblue"
        
        #cols = add.alpha(cols, alpha = 0.5)
        plot(x1,y1, log = '', bg = cols, pch = 21, main = "Log2FC", col = NA, xlab = concat(c("Log2FC (",plot_doublet_types_filtered[use_filtered[1]],")")),
             ylab = concat(c("Log2FC (",plot_doublet_types_filtered[use_filtered[2]],")")))
        
        inter_p_value = names(x)[intersect(which(x<0.05),which(y<0.05))]
        q1_signif = names(x1)[intersect(which(x1<0), which(y1>0))]
        rangex =max(x1)- min(x1)
        shift = rangex*0.6
        sep = seq(from = min(x1)+shift, to = (min(x1)+rangex*0.3)+shift, length = length(q1_signif))
        
        #mtext(q1_signif, 3, line = 0.1, at = sep, cex = 0.9)
        #for(i in c(1:length(q1_signif))){
        #  segments(x1[q1_signif[i]],y1[q1_signif[i]], sep[i], max(y1)*1.3 ,col = "grey", lwd = 2, lty = 1)
        #}
        
        genes_of_interest1 = c("CD74", "IL7R","HLA-DRA")
        genes_of_interest1 = genes_of_interest1[order(y1[genes_of_interest1]) ]
        rangey =max(y1)- min(y1)
        shift = rangey*0.6
        sep = seq(from = min(y1)+shift+1, to = (min(y1)+rangex*0.3)+shift, length = length(genes_of_interest1))
        
        mtext(genes_of_interest1, 4, line = 0.1, at = sep, cex = 0.9, las = 1)
        for(i in c(1:length(genes_of_interest1))){
          segments(x1[genes_of_interest1[i]],y1[genes_of_interest1[i]], max(x1)*1.3 ,sep[i], col = "grey", lwd = 2, lty = 1)
        }
        #text(x1[genes_of_interest1],y1[genes_of_interest1], genes_of_interest1)
        
        
}
dev.off()


####### gene set enrichment between DGE groups
gene_groups = NULL
s = 4 # combined samples only
for(ind in c(1:length(doublet_types_filtered))){
  w1 = intersect(which(overall_markers_obs[,2]==doublet_types_filtered[ind]), which(overall_markers_obs[,1]==samples_run1[s]))
  head(overall_markers_obs[w1,])
  p1= overall_markers_obs[w1[which(as.numeric(overall_markers_obs[w1,"p_val_adj"])<0.05)],]
  p1 = p1[which(as.numeric(overall_markers_obs[w1,"p_val_adj"])<0.05),]
  p1 = p1[,c("avg_log2FC","p_val_adj")]
  gene_groups = c(gene_groups, list(p1))
}
use_filtered = c(1,3)
p1 = gene_groups[[use_filtered[1]]]
p2 = gene_groups[[use_filtered[2]]]
sig_up1 = p1[intersect(which(as.numeric(p1[,1])>0), which(as.numeric(p1[,2])<0.05)),]
sig_up2 = p2[intersect(which(as.numeric(p2[,1])>0), which(as.numeric(p2[,2])<0.05)),]
sig_down1 = p1[intersect(which(as.numeric(p1[,1])<0), which(as.numeric(p1[,2])<0.05)),]
sig_down2 = p2[intersect(which(as.numeric(p2[,1])<0), which(as.numeric(p2[,2])<0.05)),]
unique_up1 = sig_up1[setdiff(rownames(sig_up1), rownames(sig_up2)), ]
unique_up2 = sig_up2[setdiff(rownames(sig_up2), rownames(sig_up1)), ]

gene_groups = c(gene_groups, list(sig_up1), list(sig_up2),list(sig_down1),list(sig_down2),list(unique_up1), list(unique_up2))
names(gene_groups) = c(doublet_types_filtered, paste(doublet_types_filtered[use_filtered], "up"), paste(doublet_types_filtered[use_filtered], "down"),
                       paste(doublet_types_filtered[use_filtered], "up unique"))



Plot_cell_surface_genes<-function(gene_groups){
  file = "~/Google_Drive/Projects/Alchemab/Projects/Gene target screening/Version 1.0/STRINGDB/9606.protein.info.v12.0.txt"
  a <- as.matrix(read.csv(file, head=T, sep="\t"))
  
  file = "~/Google_Drive/Projects/Alchemab/Projects/Gene target screening/Version 1.0/STRINGDB/9606.protein.enrichment.terms.v12.0.txt"
  p <- as.matrix(read.csv(file, head=T, sep="\t"))
  
  indir = "~/Google_Drive/Projects/Alchemab/Projects/Gene target screening/Version 1.0/Human protein atlas/"
  file = concat(c(indir, "subcellular_location.tsv"))
  b <- as.matrix(read.csv(file, head=T, sep="\t"))
  
  input_directory = "~/Google_Drive/Projects/GITHUB/PDAC150K/DATA/RECEPTOR_LIGAND_DGE/"
  file =  concat(c(input_directory, "Receptor_ligand_fantom.gsc.riken.jp.txt"))
  p <- as.matrix(read.csv(file, head=T, sep="\t"))
  ligand = p[,"Ligand.ApprovedSymbol"]
  receptor = p[,"Receptor.ApprovedSymbol"]
  lr = cbind(ligand, receptor)
  
  
  all_potential_surface_genes = NULL
  for(gind in c(1:length(gene_groups))){
    all_genes = rownames(gene_groups[[gind]])
    ############### Get subcellular location of all hits and exclude TCR/IG genes
    exclude = all_genes[c(grep("^TRAV", all_genes), grep("^TRBV", all_genes), grep("^IGHV", all_genes), grep("^IGK", all_genes), grep("^IGL", all_genes), grep("^RPS", all_genes), grep("^RPL", all_genes))]
    all_genes[which(all_genes %in% a[,"preferred_name"]==F)]
    w = which(a[,"preferred_name"] %in% all_genes)
    id1 = a[w,"X.string_protein_id"]
    names(id1) = a[w,"preferred_name"]
    id2 = a[w,"preferred_name"]
    names(id2) = a[w,"X.string_protein_id"]
    
    p1 = p[which(p[,"X.string_protein_id"] %in% id1),]
    p2 = p1[which(p1[,"category"] %in% c("Cellular Component (Gene Ontology)")),]
    cell_surface_locations = "Plasma membrane"
    cell_surface = sort(unique(id2[unique(p2[which(p2[,"description"] %in% cell_surface_locations),"X.string_protein_id"])]))
    cell_surface_locations = c("Vesicle membrane","Tricellular tight junction",
                               "Plasma membrane","Focal adhesion sites","Cell Junctions",
                               "Spanning component of plasma membrane",
                               "Side of membrane",
                               "Secretory vesicle",
                               "Receptor complex",
                               "Protein complex involved in cell adhesion",
                               "Protein complex involved in cell-cell adhesion",
                               "Protein complex involved in cell-matrix adhesion",
                               "Plasma membrane",
                               "Plasma membrane bounded cell projection",
                               "Plasma membrane bounded cell projection cytoplasm", 
                               "Plasma membrane protein complex", 
                               "Plasma membrane proton-transporting V-type ATPase complex", 
                               "Plasma membrane raft",
                               "Plasma membrane region", 
                               "Plasma membrane signaling receptor complex",
                               "Plasma lipoprotein particle",
                               "MHC class I peptide loading complex",
                               "MHC class I protein complex", 
                               "MHC class Ib protein complex", 
                               "MHC class II protein complex", 
                               "MHC protein complex",
                               "Membrane attack complex", 
                               "Membrane coat", 
                               "Membrane protein complex",
                               "Junctional membrane complex", 
                               "Intrinsic component of plasma membrane",
                               "Intrinsic component of external side of plasma membrane",
                               "Integral component of plasma membrane",
                               "Extrinsic component of cytoplasmic side of plasma membrane",
                               "Extrinsic component of external side of plasma membrane",
                               "External side of apical plasma membrane", 
                               "External side of plasma membrane",
                               "Cell-cell contact zone", 
                               "Cell-cell junction", 
                               "Cell surface")
    
    cell_surface = sort(unique(id2[unique(p2[which(p2[,"description"] %in% cell_surface_locations),"X.string_protein_id"])]))
    
    
    ## check
    genes_not_included = all_genes[which(all_genes %in% b[which(b[,"Gene.name"] %in% all_genes),"Gene.name"]==F)]
    genes_not_included = setdiff(genes_not_included, exclude)
    
    subcellular_location=as.data.frame(b[which(b[,"Gene.name"] %in% all_genes),])
    locations = subcellular_location[,c("Enhanced","Supported","Approved")]
    #locations = sort(unique(c(locations[,1], locations[,2], locations[,3])))
    cell_surface_locations = c("Plasma membrane","Focal adhesion sites","Vesicles","Cell Junctions")
    potential_cell_interaction_gene = NULL
    for(i in c(1:length(subcellular_location[,1]))){
      for(j in c(1:length(cell_surface_locations))){
        w = grep(cell_surface_locations[j], locations[i,])
        if(length(w)!=0){
          potential_cell_interaction_gene = c(potential_cell_interaction_gene, i)
        }
      }
    }
    potential_cell_interaction_gene = unique(potential_cell_interaction_gene)
    potential_cell_interaction_gene = subcellular_location[potential_cell_interaction_gene,"Gene.name"]
    
    intersect(potential_cell_interaction_gene, cell_surface)
    setdiff(cell_surface, potential_cell_interaction_gene)
    setdiff(potential_cell_interaction_gene, cell_surface)
    
    all_potential_cell_interaction_gene = sort(unique(c(potential_cell_interaction_gene, cell_surface)))
    all_potential_cell_interaction_gene= setdiff(all_potential_cell_interaction_gene, exclude)
    
    all_lr_only =all_genes[which(all_genes %in% c(ligand, receptor))]
    all_potential_cell_interaction_gene = sort(unique(c(all_potential_cell_interaction_gene, all_lr_only)))
    
    all_potential_cell_interaction_gene1 = cbind(names(gene_groups)[gind],all_potential_cell_interaction_gene)
    if(length(all_potential_surface_genes)==0){ all_potential_surface_genes= all_potential_cell_interaction_gene1
    }else{all_potential_surface_genes= rbind(all_potential_surface_genes,all_potential_cell_interaction_gene1)}
  }
  
  out_file_table = concat(c(out_dir, "Doublet_DGE_poential_surface_proteins_", batch,".txt"))
  write.table(all_potential_surface_genes, file = out_file_table, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,col.names = T, qmethod = c("escape", "double"),fileEncoding = "")
  
  ### do volcano plots?
  library(ggplot2)
  library(ggrepel)
  for(ind in c(1:length(doublet_types_filtered))){
    t1 = t_doublet_counts[,doublet_types[ind]]
    for(s in c(1:length(samples_run1))){
      file_prefix = gsub(" ","_", concat(c(doublet_types_filtered[ind],"__", samples_run1[s])), fixed = T)
      file_obs = concat(c(out_dir, "Doublet_DGE_overall_", file_prefix,"_", batch,".txt"))
      file_ran = concat(c(out_dir, "Doublet_DGE_randomised_overall_", file_prefix,"_", batch,".txt"))
      
      p <- as.matrix(read.csv(file_obs, head=T, sep="\t"))
      df <- as.data.frame(p)
      
      de <- df[complete.cases(df), ]
      de$diffexpressed <- "NO"
      # if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
      de$diffexpressed[de$avg_log2FC > 0.4 & de$p_val_adj < 0.05] <- "UP"
      # if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
      de$diffexpressed[de$avg_log2FC < -0.4 & de$p_val_adj < 0.05] <- "DOWN"
      
      de$gene_symbol = rownames(de)
      de$delabel <- NA
      #de$delabel[de$diffexpressed != "NO"] <- de$gene_symbol[de$diffexpressed != "NO"]
      up = de$gene_symbol[which(de$diffexpressed =="UP")]
      
      w = intersect(which(de$gene_symbol %in% c(setdiff(up[grep("^CD",up)], up[grep("^CDK",up)]), unique(sort(unique(unlist(all_potential_surface_genes))))), which(de$diffexpressed != "NO")))
      de$delabel[w] <- de$gene_symbol[w]
      
      
      ggplot(data=de, aes(x=avg_log2FC, y=-log10(p_val_adj), col=diffexpressed, label=delabel)) +
        geom_point() + 
        theme_minimal() +
        geom_text_repel() +
        scale_color_manual(values=c("blue", "black", "red")) +
        geom_vline(xintercept=c(-0.4, 0.4), col="red") +
        geom_hline(yintercept=-log10(0.05), col="red")      
      
    }
}
}

Get_pathway_enrichment<-function(){
  library(fgsea) ## not compatible when seurat is loaded!
  pathways.hallmark <- gmtPathways("~/Google_Drive/Projects/Single cell analysis/InfiltrImmune/Validation_datasets/TGCA_RCC/msigdb.v6.2.symbols.gmt.txt")
  w = names(pathways.hallmark)[sort(unique(c(grep("KEGG_",names(pathways.hallmark)), grep("REACTOME_",names(pathways.hallmark)), grep("PID_",names(pathways.hallmark)))))]
  w = setdiff(w, grep("GO_",names(pathways.hallmark)))
  w = grep("GO_",names(pathways.hallmark))
  pathways.hallmark1 = pathways.hallmark[w]
  
  ### fgsea
  universe = unique(rownames(overall_markers_obs))
  pathways.hallmark_sub = NULL
  pathways.hallmark_sub_names = NULL
  for(path in c(1:length(pathways.hallmark1))){
    g = intersect(pathways.hallmark1[[path]], universe)
    if(length(g)>=3){
      pathways.hallmark_sub = c(pathways.hallmark_sub, list(g))
      pathways.hallmark_sub_names = c(pathways.hallmark_sub_names, names(pathways.hallmark1)[path])
    }}
  names(pathways.hallmark_sub) = pathways.hallmark_sub_names
  names(pathways.hallmark_sub)[  grep("HLA", pathways.hallmark_sub)]
  
  for(ind in c(1:length(gene_groups))){
    ord = gene_groups[[ind]]
    rank = as.numeric(ord[,"avg_log2FC"])
    names(rank) = rownames(ord)
    rank = rank[grep("RPL", names(rank), invert = T)]
    rank = rank[grep("RPS", names(rank), invert = T)]
    rank = sort(rank)
    fgseaRes <- fgsea(pathways.hallmark_sub, rank)#, minSize = 10)#, scoreType = "pos")
    range(fgseaRes$pval)
    fgseaResTidy = fgseaRes[which(fgseaRes$pval<0.05),]
    fgseaResTidy = fgseaResTidy[order(fgseaResTidy[,"pval"])]
    
    fgseaResTidy[grep("FOS", fgseaResTidy$leadingEdge),]
    fgseaResTidy$pathway[grep("LYMPH",fgseaResTidy$pathway)]
    
    
    if(ind %in% c(3,5)){
      pathways_to_plot = c("GO_RESPONSE_TO_CYTOKINE","GO_POSITIVE_REGULATION_OF_IMMUNE_SYSTEM_PROCESS","GO_REGULATION_OF_CELL_PROLIFERATION","GO_POSITIVE_REGULATION_OF_CELL_COMMUNICATION","GO_LYMPHOCYTE_MEDIATED_IMMUNITY")
      plotEnrichment(pathways.hallmark_sub[[pathways_to_plot[1]]],
                     rank) + labs(title=pathways_to_plot)
      
      fileout1=concat(c(out_dir,"Doublet_enrichment_", batch,"_",names(gene_groups)[ind],".pdf"))
      w=2.25
      pdf(file=fileout1, height=w*1, width=w*4.5)
      par(mfrow= c(1,1), mar = c(5,5,6,6))
      plotGseaTable(pathways.hallmark_sub[pathways_to_plot], rank, fgseaRes, 
                    gseaParam = 0.5,colwidths = c(8, 1, 0.8, 1.2, 1.2))
      dev.off()
      
    }
    
  
Get_pathway_enrichment_not_use<-function(){
  library(fgsea) ## not compatible when seurat is loaded!
  pathways.hallmark <- gmtPathways("~/Google_Drive/Projects/Single cell analysis/InfiltrImmune/Validation_datasets/TGCA_RCC/msigdb.v6.2.symbols.gmt.txt")
  w = names(pathways.hallmark)[sort(unique(c(grep("KEGG_",names(pathways.hallmark)), grep("REACTOME_",names(pathways.hallmark)), grep("PID_",names(pathways.hallmark)))))]
  w = setdiff(w, grep("GO_",names(pathways.hallmark)))
  pathways.hallmark1 = pathways.hallmark[w]
  
  
  library(ggplot2)
  library(ggrepel)
  for(ind in c(1:length(doublet_types_filtered))){
    t1 = t_doublet_counts[,doublet_types[ind]]
    for(s in c(4)){
      file_prefix = gsub(" ","_", concat(c(doublet_types_filtered[ind],"__", samples_run1[s])), fixed = T)
      file_obs = concat(c(out_dir, "Doublet_DGE_overall_", file_prefix,"_", batch,".txt"))
      file_ran = concat(c(out_dir, "Doublet_DGE_randomised_overall_", file_prefix,"_", batch,".txt"))
      p <- as.matrix(read.csv(file_obs, head=T, sep="\t"))
      
      ### fgsea
      universe = rownames(overall_markers_obs)
      pathways.hallmark_sub = NULL
      pathways.hallmark_sub_names = NULL
      for(path in c(1:length(pathways.hallmark1))){
        g = intersect(pathways.hallmark1[[path]], universe)
        if(length(g)>=4){
          pathways.hallmark_sub = c(pathways.hallmark_sub, list(g))
          pathways.hallmark_sub_names = c(pathways.hallmark_sub_names, names(pathways.hallmark1)[path])
        }}
      names(pathways.hallmark_sub) = pathways.hallmark_sub_names
      
      
      
      
      
      library(org.Hs.eg.db)
      universe_entrez = select(org.Hs.eg.db, 
                               keys = universe,
                               columns = c("ENTREZID", "SYMBOL"),
                               keytype = "SYMBOL")
      universe_entrez = universe_entrez[,"ENTREZID"]
      pathways= reactomePathways(universe_entrez)
      
      sig_threshold =0.05/length(LF_gene_mat[1,])
      list_pathways = NULL
      list_names = NULL
      for(f in c(1:length(LF_gene_mat[1,]))){
        c = 2
        if(c==1){
          ranks = sort(LF_gene_mat[,f])
          ranks = ranks[which(ranks!=0)]
          fgseaRes <- fgsea(pathways.hallmark_sub, ranks)
          range(fgseaRes$padj)
          fgseaResTidy = fgseaRes[which(fgseaRes$padj<0.05),]
          fgseaResTidy = fgseaResTidy[order(fgseaResTidy[,"pval"])]
          
          topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
          topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
          topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
          plotGseaTable(pathways.hallmark_sub[topPathways], ranks, fgseaRes, gseaParam=0.5)
          plotEnrichment(pathways.hallmark_sub[["GO_DEFENSE_RESPONSE"]],ranks) + labs(title="GO_DEFENSE_RESPONSE")
        }
        
      }
      
}    
  }
}

Somethingelse<-function(){

  
  
  
  
    
    ############### plot heatmap of cell surface genes
    colnames_use = colnames(mat_p_values)[grep("RANDOM", colnames(mat_p_values),invert = T)]
    colnames_use = colnames_use[grep("TOTAL", colnames_use,invert = T)]
    colnames_use1 = gsub("B CELL ","",colnames_use, fixed = T)
    colnames_use1 = gsub("T CELL ","",colnames_use1, fixed = T)
    colnames_use1 = gsub("B.cell.memory ","B.cell.memory [",colnames_use1, fixed = T)
    colnames_use1 = gsub(".."," - ",colnames_use1, fixed = T)
    colnames_use1 = gsub("."," ",colnames_use1, fixed = T)
    colnames_use1 = paste(colnames_use1,"]", sep = "")
    mat_p_values1 = mat_p_values[all_potential_cell_interaction_gene,colnames_use]
    mat_FC1 = mat_FC[all_potential_cell_interaction_gene,colnames_use]
    colnames(mat_p_values1 ) = colnames_use1
    colnames(mat_FC1 ) = colnames_use1
    
    a = apply(mat_p_values1, 1, function(x){length(which(x<0.05))})
    mat_p_values1 = mat_p_values1[which(a>=1),]
    mat_FC1 = mat_FC1[which(a>=1),]
    
    
    a = apply(mat_p_values1[,grep("T CELL", colnames_use)], 1, function(x){length(which(x<0.05))})
    b = apply(mat_p_values1[,grep("B CELL", colnames_use)], 1, function(x){length(which(x<0.05))})
    
    a[which(a>=2)]
    b[which(b>=2)]
    
    library(gplots)
    fileout1=concat(c(outdir, "/Plot_DGE_", batch,"_cell_surface.pdf"))
    w=6
    pdf(file=fileout1, height=w*2.9, width=w*0.5)
    par(mfrow= c(1,1), mar = c(12,4,4,4))
    heatmap.2(mat_FC1, col="bluered", scale="none", trace="none", cexCol = 0.8,cexRow = 0.8, margins=c(25,8), dendrogram = "row",Colv=FALSE)
    dev.off()
    
    
  
  a = apply(mat_FC1, 1, function(x){length(which(x>0))})
  mat_p_values1 = mat_p_values1[which(a>=1),]
  mat_FC1 = mat_FC1[which(a>=1),]
  
  library(gplots)
  fileout1=concat(c(outdir, "/Plot_DGE_", batch,"_cell_surface_up.pdf"))
  w=6
  pdf(file=fileout1, height=w*1.6, width=w*0.5)
  par(mfrow= c(1,1), mar = c(12,4,4,4))
  heatmap.2(mat_FC1, col="bluered", scale="none", trace="none", cexCol = 0.8,cexRow = 0.8, margins=c(25,8), dendrogram = "row",Colv=FALSE)
  dev.off()
  
  
  ############### TOTAL
  colnames_use = colnames(mat_p_values)[grep("RANDOM", colnames(mat_p_values),invert = T)]
  colnames_use = colnames_use[grep("TOTAL", colnames_use,invert = F)]
  colnames_use1 = gsub("B CELL ","",colnames_use, fixed = T)
  colnames_use1 = gsub("T CELL ","",colnames_use1, fixed = T)
  colnames_use1 = gsub("B.cell.memory ","B.cell.memory [",colnames_use1, fixed = T)
  colnames_use1 = gsub(".."," - ",colnames_use1, fixed = T)
  colnames_use1 = gsub("."," ",colnames_use1, fixed = T)
  colnames_use1 = paste(colnames_use1,"]", sep = "")
  colnames_use1[grep("B CELL",colnames_use)] = gsub("TOTAL","B cells", colnames_use1[grep("B CELL",colnames_use)])
  colnames_use1[grep("T CELL",colnames_use)] = gsub("TOTAL","T cells", colnames_use1[grep("T CELL",colnames_use)])
  mat_p_values1 = mat_p_values[all_potential_cell_interaction_gene,colnames_use]
  mat_FC1 = mat_FC[all_potential_cell_interaction_gene,colnames_use]
  colnames(mat_p_values1 ) = colnames_use1
  colnames(mat_FC1 ) = colnames_use1
  
  a = apply(mat_p_values1, 1, function(x){length(which(x<0.05))})
  mat_p_values1 = mat_p_values1[which(a>=1),]
  mat_FC1 = mat_FC1[which(a>=1),]
  
  
  a = apply(mat_p_values1[,grep("T CELL", colnames_use)], 1, function(x){length(which(x<0.05))})
  b = apply(mat_p_values1[,grep("B CELL", colnames_use)], 1, function(x){length(which(x<0.05))})
  
  a[which(a>=2)]
  b[which(b>=2)]
  
  library(gplots)
  fileout1=concat(c(outdir, "/Plot_DGE_", batch,"_cell_surface_TOTAL.pdf"))
  w=6
  pdf(file=fileout1, height=w*3.1, width=w*0.5)
  par(mfrow= c(1,1), mar = c(12,4,4,4))
  heatmap.2(mat_FC1, col="bluered", scale="none", trace="none", cexCol = 0.8,cexRow = 0.8, margins=c(25,8), dendrogram = "row",Colv=FALSE)
  dev.off()
  
  
  a = apply(mat_FC1, 1, function(x){length(which(x>0))})
  mat_p_values1 = mat_p_values1[which(a>=1),]
  mat_FC1 = mat_FC1[which(a>=1),]
  
  library(gplots)
  fileout1=concat(c(outdir, "/Plot_DGE_", batch,"_cell_surface_TOTAL_up.pdf"))
  w=6
  pdf(file=fileout1, height=w*1.6, width=w*0.5)
  par(mfrow= c(1,1), mar = c(12,4,4,4))
  heatmap.2(mat_FC1, col="bluered", scale="none", trace="none", cexCol = 0.8,cexRow = 0.8, margins=c(25,8), dendrogram = "row",Colv=FALSE)
  dev.off()
  
  ######## with receptor ligands too
  ########### get protein location of DGE genes
  Get_protein_locations<-function(){
    ### get receptors and ligands
    input_directory = "~/Google_Drive/Projects/GITHUB/PDAC150K/DATA/RECEPTOR_LIGAND_DGE/"
    file =  concat(c(input_directory, "Receptor_ligand_fantom.gsc.riken.jp.txt"))
    p <- as.matrix(read.csv(file, head=T, sep="\t"))
    ligand = p[,"Ligand.ApprovedSymbol"]
    receptor = p[,"Receptor.ApprovedSymbol"]
    cell_surface_receptor_ligands = sort(unique(c(ligand, receptor)))
    return(cbind(ligand, receptor))
  }
  
  ligand_receptor  = Get_protein_locations()
  all_potential_cell_interaction_gene = intersect(all_genes, sort(ligand_receptor))
  all_potential_cell_interaction_gene= sort(unique(c(all_potential_cell_interaction_gene, all_potential_cell_interaction_gene1)))
  
  ############### plot heatmap of cell surface genes
  colnames_use = colnames(mat_p_values)[grep("RANDOM", colnames(mat_p_values),invert = T)]
  colnames_use = colnames_use[grep("TOTAL", colnames_use,invert = T)]
  colnames_use1 = gsub("B CELL ","",colnames_use, fixed = T)
  colnames_use1 = gsub("T CELL ","",colnames_use1, fixed = T)
  colnames_use1 = gsub("B.cell.memory ","B.cell.memory [",colnames_use1, fixed = T)
  colnames_use1 = gsub(".."," - ",colnames_use1, fixed = T)
  colnames_use1 = gsub("."," ",colnames_use1, fixed = T)
  colnames_use1 = paste(colnames_use1,"]", sep = "")
  mat_p_values1 = mat_p_values[all_potential_cell_interaction_gene,colnames_use]
  mat_FC1 = mat_FC[all_potential_cell_interaction_gene,colnames_use]
  colnames(mat_p_values1 ) = colnames_use1
  colnames(mat_FC1 ) = colnames_use1
  
  a = apply(mat_p_values1, 1, function(x){length(which(x<0.05))})
  mat_p_values1 = mat_p_values1[which(a>=1),]
  mat_FC1 = mat_FC1[which(a>=1),]
  
  
  a = apply(mat_p_values1[,grep("T CELL", colnames_use)], 1, function(x){length(which(x<0.05))})
  b = apply(mat_p_values1[,grep("B CELL", colnames_use)], 1, function(x){length(which(x<0.05))})
  
  a[which(a>=2)]
  b[which(b>=2)]
  
  library(gplots)
  fileout1=concat(c(outdir, "/Plot_DGE_", batch,"_cell_surface_RL.pdf"))
  w=6
  pdf(file=fileout1, height=w*2.9, width=w*0.5)
  par(mfrow= c(1,1), mar = c(12,4,4,4))
  heatmap.2(mat_FC1, col="bluered", scale="none", trace="none", cexCol = 0.8,cexRow = 0.8, margins=c(25,8), dendrogram = "row",Colv=FALSE)
  dev.off()
  
  a = apply(mat_FC1, 1, function(x){length(which(x>0))})
  mat_p_values1 = mat_p_values1[which(a>=1),]
  mat_FC1 = mat_FC1[which(a>=1),]
  
  library(gplots)
  fileout1=concat(c(outdir, "/Plot_DGE_", batch,"_cell_surface_RL_up.pdf"))
  w=6
  pdf(file=fileout1, height=w*1.6, width=w*0.5)
  par(mfrow= c(1,1), mar = c(12,4,4,4))
  heatmap.2(mat_FC1, col="bluered", scale="none", trace="none", cexCol = 0.8,cexRow = 0.8, margins=c(25,8), dendrogram = "row",Colv=FALSE)
  dev.off()
  
  
}}

}

Plot_cell_receptor_ligand_genes<-function(all_genes){
  
  ########### get protein location of DGE genes
  Get_protein_locations<-function(){
    ### get receptors and ligands
    input_directory = "~/Google_Drive/Projects/GITHUB/PDAC150K/DATA/RECEPTOR_LIGAND_DGE/"
    file =  concat(c(input_directory, "Receptor_ligand_fantom.gsc.riken.jp.txt"))
    p <- as.matrix(read.csv(file, head=T, sep="\t"))
    ligand = p[,"Ligand.ApprovedSymbol"]
    receptor = p[,"Receptor.ApprovedSymbol"]
    cell_surface_receptor_ligands = sort(unique(c(ligand, receptor)))
    return(cbind(ligand, receptor))
  }
  
  ligand_receptor  = Get_protein_locations()
  all_potential_cell_interaction_gene = intersect(all_genes, sort(ligand_receptor))
  
  
  ############### plot heatmap of cell surface genes
  colnames_use = colnames(mat_p_values)[grep("RANDOM", colnames(mat_p_values),invert = T)]
  colnames_use = colnames_use[grep("TOTAL", colnames_use,invert = T)]
  colnames_use1 = gsub("B CELL ","",colnames_use, fixed = T)
  colnames_use1 = gsub("T CELL ","",colnames_use1, fixed = T)
  colnames_use1 = gsub("B.cell.memory ","B.cell.memory [",colnames_use1, fixed = T)
  colnames_use1 = gsub(".."," - ",colnames_use1, fixed = T)
  colnames_use1 = gsub("."," ",colnames_use1, fixed = T)
  colnames_use1 = paste(colnames_use1,"]", sep = "")
  mat_p_values1 = mat_p_values[all_potential_cell_interaction_gene,colnames_use]
  mat_FC1 = mat_FC[all_potential_cell_interaction_gene,colnames_use]
  colnames(mat_p_values1 ) = colnames_use1
  colnames(mat_FC1 ) = colnames_use1
  
  a = apply(mat_p_values1, 1, function(x){length(which(x<0.05))})
  mat_p_values1 = mat_p_values1[which(a>=1),]
  mat_FC1 = mat_FC1[which(a>=1),]
  
  
  a = apply(mat_p_values1[,grep("T CELL", colnames_use)], 1, function(x){length(which(x<0.05))})
  b = apply(mat_p_values1[,grep("B CELL", colnames_use)], 1, function(x){length(which(x<0.05))})
  
  a[which(a>=2)]
  b[which(b>=2)]
  
  library(gplots)
  fileout1=concat(c(outdir, "/Plot_DGE_", batch,"_receptor_ligand.pdf"))
  w=6
  pdf(file=fileout1, height=w*1.4, width=w*0.5)
  par(mfrow= c(1,1), mar = c(12,4,4,4))
  heatmap.2(mat_FC1, col="bluered", scale="none", trace="none", cexCol = 0.8,cexRow = 0.8, margins=c(25,8), dendrogram = "row",Colv=FALSE)
  dev.off()
}





