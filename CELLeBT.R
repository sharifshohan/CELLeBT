# #### Run on Oxford compute set up
# # srun -p short --cpus-per-task 10 --pty bash
module purge
module load R-bundle-Bioconductor/3.9-foss-2019a-R-3.6.0
module load HDF5/1.10.5-gompi-2019a
module load umap-learn/0.3.10-foss-2019a-Python-3.7.2
module load Seurat/3.1.2-foss-2019a-R-3.6.0
module load Harmony/1.0.0-foss-2019a-R-3.6.0
R

#change directory according to need
setwd("/well/immune-rep/users/npy727/B_T_doublet/percentage_b_t_cell/")
.libPaths(c( "~/R/3.4.3-openblas-0.2.18-omp-gcc5.4.0", .libPaths()))
.libPaths(c("/well/immune-rep/users/npy727/INSTALLED_PROGRAMS/R_MODULES", .libPaths()))
options(repos = 'http://cran.ma.imperial.ac.uk/')

library(ROCR)
library(viridisLite)
library(viridis)
library("Seurat")
library('harmony') 
library(ggplot2)
library(pryr)




library(future) 

#This is a concatenate function Rachael often uses
concat = function(v) {
  res = ""
  for (i in 1:length(v)){res = paste(res,v[i],sep="")}
  res
}

#Make a vector of colours for plotting
cols<- c(
  "green4", "#E31A1C", # red
          "dodgerblue2",
               "#6A3D9A", # purple
               "#FF7F00", # orange
               "black", "gold1",
               "skyblue2", "#FB9A99", # lt pink
               "palegreen2",
               "#CAB2D6", # lt purple
               "#FDBF6F", # lt orange
               "gray70", "khaki2",
               "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
               "darkturquoise", "green1", "yellow4", "yellow3",
               "darkorange4", "brown"
)
##Read the pbmc File
##RCC
pbmc<- readRDS("/well/immune-rep/shared/10X_GENOMICS/RCC_WORKING_DATA/FINAL//Seurat_full_annotation_RCCmetP_RCC.pbmc")

# #Make a column with the broad cell type taking only the T cell and B cells for later part
# broad_cell_type <- strsplit(pbmc@meta.data$cell_refined_annotation, " ", fixed = T)
# 
# for(i in c(1:length(broad_cell_type))){
#   broad_cell_type[i] = concat(c(broad_cell_type[[i]][1]," ",broad_cell_type[[i]][2]))
# }
# broad_cell_type = unlist(broad_cell_type)
# broad_cell_type <- ifelse(broad_cell_type =="- NA", "-", broad_cell_type)
# pbmc@meta.data$broad_cell_type <- broad_cell_type

#Changing the names for further analysis
#cell <- read.csv("/well/immune-rep/users/npy727/B_T_doublet/working_from_Rachaels_code/PDAC/Seurat_merad_annotations_PDAC150Ka_2.txt", sep = "\t") 
cell <- read.csv("/well/immune-rep/users/npy727/B_T_doublet/working_from_Rachaels_code/PDAC/cell_type_merging.csv")
cell$cell_types_merad <- as.character(cell$cell_types_merad)
cell$Level.1 <- as.character(cell$Level.1)
pbmc@meta.data$cell_type <- pbmc@meta.data$cell_refined_annotation
cell_type <- pbmc@meta.data$cell_type
#names(cell_type) <- rownames(pbmc@meta.data)
#head(cell_type)
cell_type_broad = cell_type
for(i in c(1:length(cell$Level.1))){
  print(i)
  cell_type_broad[which(cell_type_broad==cell$cell_types_merad[i])] = cell$Level.1[i]
}
pbmc@meta.data$cell_refined_annotation <- cell_type_broad


pbmc_real <- pbmc
#get the main sample directly
#get the main sample directly
pbmc<- subset(pbmc, subset = orig.ident == "04_Kidney_met_panc1_blood_CD45p1")
pbmc<- subset(pbmc, subset = orig.ident == "01_Kidney_met_panc1_biopsy_CD45p1")
pbmc<- subset(pbmc, subset = orig.ident == "02_Kidney_met_panc1_biopsy_CD45p2")
head(pbmc@meta.data)


table(pbmc@meta.data$cell_refined_annotation)
## filter pbmc for "B/T_cell" column and Get rid of the low number of cell and gamma/delta cells


Cells_to_keep <- function(pbmc, threshold=30){
  cells_to_keep <- table(pbmc@meta.data$cell_refined_annotation)
  cells_to_keep <- names(cells_to_keep[cells_to_keep>threshold])
  cells_to_keep <- cells_to_keep[cells_to_keep != "Gamma/delta" & cells_to_keep != "Activated B cell" & cells_to_keep != "CXCR4lo CD27lo memory B cell" &cells_to_keep != "T cell gdT"]
  for(i in 1:length(pbmc@meta.data$cell_refined_annotation)){
    for(j in 1:length(cells_to_keep)){
      if(pbmc@meta.data$cell_refined_annotation[i] == cells_to_keep[j]){
        pbmc@meta.data$cells_final[i] <- "Yes"
        break
      }else{
        pbmc@meta.data$cells_final[i] <- "No"
      }
    }
  }
  pbmc<- subset(pbmc, subset= cells_final == "Yes")
  return(pbmc)
}

pbmc<- Cells_to_keep(pbmc)

### get the user to define a column specifying which are B/T cells ("B/T_cell" or other)
broad_cell_type = pbmc@meta.data$broad_cell_type
cell_type_broad = broad_cell_type
## check that you have all the cell groups
cell_type_broad[c(grep("B cell", cell_type_broad), grep("T cell", cell_type_broad), grep("-", cell_type_broad))] = "B/T/D_cell"
#cell_type_broad[which(cell_type_broad!= "B/T_cell")]="-"
pbmc@meta.data$B_Tcell = cell_type_broad

pbmc<- subset(pbmc, subset = B_Tcell =="B/T/D_cell")


cell_ids = rownames(pbmc@meta.data)
#subset only the RCC met Pan Sample only from the pbmc Dataset
# sample = pbmc@meta.data$orig.ident
# sample[grep("panc1", sample)] <- "RCC"
# sample[which(sample != "RCC")] <- "-"
# pbmc@meta.data$sample <- sample
# pbmc<- subset(pbmc, subset = sample == "RCC")


## the get those overlapping with VDJ

#### get VDJ information

m_VDJ <- readRDS("/well/immune-rep/shared/10X_GENOMICS/RCC_WORKING_DATA//VDJ_information_RCCmetP_RCC.VDJ")
m_VDJ_BCR = m_VDJ$BCR
m_VDJ_TCR = m_VDJ$TCR
# #PDAC
# 
# m_VDJ = readRDS(file=concat(c("/well/immune-rep/shared/10X_GENOMICS/PDAK150K_WORKING_DATA/VDJ/VDJ_information_PDAC150Ka.VDJ")))
# m_VDJ_BCR = m_VDJ$BCR
# m_VDJ_TCR = m_VDJ$TCR

### match names (RCC met pan only)
cell_ids_filtered = cell_ids
#cell_ids_filtered = cell_ids[grep("panc1", cell_ids)] # GEX RCC met pan only IDs

#check: rownames(m_VDJ_BCR) == rownames(m_VDJ_TCR)
cell_ids_VDJ = rownames(m_VDJ_BCR)
cell_ids_VDJ_s = strsplit(cell_ids_VDJ, "||", fixed = T)

# for(i in c(1:length(cell_ids_VDJ_s))){
#   cell_ids_VDJ_s[i] = concat(c(cell_ids_VDJ_s[[i]][1]))
# }
# cell_ids_VDJ_s = unlist(cell_ids_VDJ_s)


for(i in c(1:length(cell_ids_VDJ_s))){
  cell_ids_VDJ_s[i] = concat(c(cell_ids_VDJ_s[[i]][2],"_",cell_ids_VDJ_s[[i]][2],"_",cell_ids_VDJ_s[[i]][1], "-1"))
}
cell_ids_VDJ_s = unlist(cell_ids_VDJ_s)



####For the RCC data
cell_ids_VDJ = rownames(m_VDJ_BCR)
cell_ids_VDJ_s = strsplit(cell_ids_VDJ, "||", fixed = T)
for(i in c(1:length(cell_ids_VDJ_s))){
  cell_ids_VDJ_s[i] = concat(c(cell_ids_VDJ_s[[i]][2],"_",cell_ids_VDJ_s[[i]][1]))
}
cell_ids_VDJ_s = unlist(cell_ids_VDJ_s)


# cell_ids_VDJ_s = strsplit(cell_ids_VDJ, "_", fixed = T)
# 
# for(i in c(1:length(cell_ids_VDJ_s))){
#   cell_ids_VDJ_s[i] = concat(c(cell_ids_VDJ_s[[i]][3]))
# }
# cell_ids_VDJ_s = unlist(cell_ids_VDJ_s)
# cell_ids_VDJ_s = strsplit(cell_ids_VDJ, "-", fixed = T)
# for(i in c(1:length(cell_ids_VDJ_s))){
#   cell_ids_VDJ_s[i] = concat(c(cell_ids_VDJ_s[[i]][1]))
# }
# 
# cell_ids_VDJ_s = unlist(cell_ids_VDJ_s)
# for(i in c(1:length(cell_ids_VDJ_s))){
#   cell_ids_VDJ_s[i] = concat(c(cell_ids_VDJ_s[[i]][2],"_",cell_ids_VDJ_s[[i]][2],"_",cell_ids_VDJ_s[[i]][1], "-1"))
# }
# cell_ids_VDJ_s = unlist(cell_ids_VDJ_s)
# length(intersect(cell_ids_VDJ_s[grep("panc1", cell_ids_VDJ_s)], cell_ids_filtered))
# length(cell_ids_filtered)
# length(cell_ids_VDJ_s[grep("panc1", cell_ids_VDJ_s)])

#relabel VDJ IDs to align with GEX
rownames(m_VDJ_BCR) = cell_ids_VDJ_s
rownames(m_VDJ_TCR) = cell_ids_VDJ_s

#include only  VDJ IDs that are present also with GEX data
inter= intersect(rownames(m_VDJ_TCR), cell_ids_filtered)
m_VDJ_BCR = m_VDJ_BCR[inter,]
m_VDJ_TCR = m_VDJ_TCR[inter,]

# ##Repeatation of same process
# ## only consider IDs from RCC met pan. 
# ids_BCR_filter = rownames(m_VDJ_BCR)[grep("Kidney_met_panc1",rownames(m_VDJ_BCR))]
# ids_TCR_filter = rownames(m_VDJ_TCR)[grep("Kidney_met_panc1",rownames(m_VDJ_TCR))]
# m_VDJ_BCR = m_VDJ_BCR[ids_BCR_filter,]
# m_VDJ_TCR = m_VDJ_TCR[ids_TCR_filter,]
# 


### ensure no counts are zeroed

m_VDJ_TCR[which(m_VDJ_TCR[,"n_umis1"]=='-'),"n_umis1"] = 0
m_VDJ_TCR[which(m_VDJ_TCR[,"n_umis2"]=='-'),"n_umis2"] = 0
m_VDJ_BCR[which(m_VDJ_BCR[,"n_umis1"]=='-'),"n_umis1"] = 0
m_VDJ_BCR[which(m_VDJ_BCR[,"n_umis2"]=='-'),"n_umis2"] = 0
m_VDJ_TCR[which(m_VDJ_TCR[,"mixed_contig_n_umis2"]=='-'),"mixed_contig_n_umis2"] = 0
m_VDJ_BCR[which(m_VDJ_BCR[,"mixed_contig_n_umis1"]=='-'),"mixed_contig_n_umis1"] = 0

VDJ_ids = rownames(m_VDJ_TCR)


#Final pbmc that matches the VDJ data and has doublets
pbmc <- subset(pbmc,cell = VDJ_ids)
#probable doublet+garbage cells
doublet_with_ND_cells <- rownames(pbmc@meta.data[pbmc@meta.data$broad_cell_type =="-",])


threshold_umis_all = c(1:3) ## this can be edited
#High Confidence Doublets
#High confidence T-B cell doublets: BCR+TCR+ (sum of >=3 umis for both)

names_cols = c("homotypic_doublets_B",
               "homotypic_doublets_T",
               "heterotypic_doublets",
               "triplets_2B_T",
               "triplets_B_2T","quadruplets","other")

mat_hc = t(matrix(data = 0, nrow = length(threshold_umis_all), ncol = length(names_cols), dimnames = c(list(threshold_umis_all), list(names_cols))))
mat_lc = t(matrix(data = 0, nrow = length(threshold_umis_all), ncol = length(names_cols), dimnames = c(list(threshold_umis_all), list(names_cols))))

n_hc_singlets = NULL
classification_list = NULL
i=3


hc_list = c(rep(list(c()), length(names_cols)))
lc_list = c(rep(list(c()), length(names_cols)))
names(hc_list) = names_cols
names(lc_list) = names_cols

threshold_umis = threshold_umis_all[i]
h_c_TCR_ids = VDJ_ids[which((as.numeric(m_VDJ_TCR[VDJ_ids,"n_umis1"]) + as.numeric(m_VDJ_TCR[VDJ_ids,"n_umis2"]))>=threshold_umis)]
h_c_BCR_ids = VDJ_ids[which((as.numeric(m_VDJ_BCR[VDJ_ids,"n_umis1"]) + as.numeric(m_VDJ_BCR[VDJ_ids,"n_umis2"]))>=threshold_umis)]
second_TCR = VDJ_ids[which(as.numeric(m_VDJ_TCR[VDJ_ids,"mixed_contig_n_umis2"])>=threshold_umis)]
second_BCR = VDJ_ids[which(as.numeric(m_VDJ_BCR[VDJ_ids,"mixed_contig_n_umis1"])>=threshold_umis)]

hc_list[["quadruplets"]] = intersect(second_TCR, second_BCR)
hc_list[["triplets_2B_T"]] = setdiff(intersect(h_c_TCR_ids, second_BCR), unlist(hc_list))
hc_list[["triplets_B_2T"]] = setdiff(intersect(h_c_BCR_ids, second_TCR), unlist(hc_list))
hc_list[["homotypic_doublets_B"]] = setdiff(second_BCR, unlist(hc_list))
hc_list[["homotypic_doublets_T"]] = setdiff(second_TCR, unlist(hc_list))
hc_list[["heterotypic_doublets"]] = setdiff(intersect(h_c_TCR_ids, h_c_BCR_ids), unlist(hc_list))

all_hc = unlist(hc_list)
names(all_hc) <- NULL

threshold_umis = 1 ## for the low confidence droplets
l_c_TCR_ids = VDJ_ids[which((as.numeric(m_VDJ_TCR[VDJ_ids,"n_umis1"]) + as.numeric(m_VDJ_TCR[VDJ_ids,"n_umis2"]))>=threshold_umis)]
l_c_BCR_ids = VDJ_ids[which((as.numeric(m_VDJ_BCR[VDJ_ids,"n_umis1"]) + as.numeric(m_VDJ_BCR[VDJ_ids,"n_umis2"]))>=threshold_umis)]
second_TCR = VDJ_ids[which(as.numeric(m_VDJ_TCR[VDJ_ids,"mixed_contig_n_umis2"])>=threshold_umis)]
second_BCR = VDJ_ids[which(as.numeric(m_VDJ_BCR[VDJ_ids,"mixed_contig_n_umis1"])>=threshold_umis)]

lc_list[["quadruplets"]] = setdiff(intersect(second_TCR, second_BCR), all_hc)
lc_list[["triplets_2B_T"]] = setdiff(intersect(l_c_TCR_ids, second_BCR), c(unlist(lc_list), all_hc))
lc_list[["triplets_B_2T"]] = setdiff(intersect(l_c_BCR_ids, second_TCR), c(unlist(lc_list), all_hc))
lc_list[["homotypic_doublets_B"]] = setdiff(second_BCR, c(unlist(lc_list), all_hc))
lc_list[["homotypic_doublets_T"]] = setdiff(second_TCR, c(unlist(lc_list), all_hc))
lc_list[["heterotypic_doublets"]] = setdiff(intersect(h_c_TCR_ids, h_c_BCR_ids), c(unlist(lc_list), all_hc))
ND_cells <- setdiff(doublet_with_ND_cells,c(unlist(hc_list),unlist(lc_list)))
lc_list[["other"]] = ND_cells
hc_singlets = setdiff(VDJ_ids, c(unlist(hc_list),unlist(lc_list), ND_cells ))

## (1) get nUMI outliers from singlets: Singlet_outlier function

#outlier replacement function 
outlierreplacement <- function(dataframe){
  dataframe %>%          
    map_if(is.numeric, ~ replace(.x, .x %in% boxplot.stats(.x)$out, NA)) %>%
    bind_cols
}

library(RColorBrewer)
Singlet_outlier <- function(pbmc, hc_singlets){
  
  pbmc_singlet <- subset(pbmc, cells = hc_singlets)
  x<- pbmc_singlet@meta.data
  pdf("Before Removal Singlet.pdf", height = 5, width = 5)
  p<- ggplot(x, aes(x = cell_refined_annotation, y = nCount_RNA , fill = cell_refined_annotation))+
    geom_boxplot()+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
          axis.text.y = element_text(vjust = 0.5,
                                     hjust = 0.5)) +
    labs(x = "Cell Type",y = "nUMI" )+
    scale_fill_manual(values = cols)
  print(p)
  dev.off()
  
  x$cell_id <- rownames(x)
  m<- x[,c("nCount_RNA","cell_refined_annotation","cell_id")]
  library(tidyr)
  library(dplyr)
  library(tidyr)
  library(purrr)
  y<- spread(m, key =  cell_refined_annotation, value = nCount_RNA)
  k<- outlierreplacement(y)
  z<- gather(k, key = "cell_refined_annotation", value = "nCount_RNA",2:ncol(k) )
  z<- z[!is.na(z$nCount_RNA),]
  outliers <- setdiff(k$cell_id, z$cell_id)
  hc_singlets = setdiff(hc_singlets, outliers)
  pbmc_singlet <- subset(pbmc, cells = hc_singlets)
  x<- pbmc_singlet@meta.data
  pdf("after Removal Singlet.pdf", height = 5, width = 5)
  p<- ggplot(x, aes(x = cell_refined_annotation, y = nCount_RNA , fill = cell_refined_annotation))+
    geom_boxplot()+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
          axis.text.y = element_text(vjust = 0.5,
                                     hjust = 0.5)) +
    labs(x = "Cell Type",y = "nUMI" )+
    scale_fill_manual(values = cols)
  print(p)
  dev.off()
  
  return(outliers)
  
}
singlet_outliers <- Singlet_outlier(pbmc,hc_singlets)
lc_list[["other"]] = append(lc_list[["other"]],singlet_outliers)
hc_singlets = setdiff(hc_singlets, singlet_outliers)

## (2) get MLplet prediction from doublets
hc_doublets = hc_list[["heterotypic_doublets"]]

### Identify high confidence singlets removing cells that express B cell markers in T cell or T cell markers in B cell
pbmc <- SetIdent(pbmc,value = "cell_refined_annotation")

bmarkers <- c("CD19", "CD22", "CD79A", "CD79B", "MS4A1", "BLK", "BANK1", 
              "TCL1A","CD52")

tmarkers <- c("CD8A","CD8B","CD28","CCR7","TCF7","CD3D","CD4","CD3E",
              "TRAC","TRBC1","TRBC2","TIGIT","GZMK","IL7R","LCK")

bcell_markers <- intersect(bmarkers, rownames(pbmc@assays$RNA@counts))
tcell_markers <- intersect(tmarkers , rownames(pbmc@assays$RNA@counts))

pbmc <- AddModuleScore(pbmc, features = list(bcell_markers), 
                       name = "B_cell_score", 
                       assay = "RNA", 
                       overwrite = TRUE)

pbmc <- AddModuleScore(pbmc, features = list(tcell_markers), 
                       name = "t_cell_score", 
                       assay = "RNA", 
                       overwrite = TRUE)
pdf("T_cell_module_score.pdf", height = 5, width = 5)
FeaturePlot(pbmc,
            features = "t_cell_score1", label = TRUE, repel = TRUE) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
dev.off()
pdf("B_cell_module_score.pdf", height = 5, width = 5)
FeaturePlot(pbmc,
            features = "B_cell_score1", label = TRUE, repel = TRUE) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
dev.off()


#Subset B cell and T cell
bcell <- subset(pbmc, subset = broad_cell_type== "B cell")
tcell <- subset(pbmc, subset = broad_cell_type== "T cell")


#Find High Confidence B cell and High Confidence T cell
iqr_B <- IQR(bcell@meta.data$t_cell_score1)
qB<- quantile(bcell@meta.data$t_cell_score1, probs = 0.75)

qn_b <- qB + 1.5*iqr_B
hc_bcells <- colnames(subset(bcell, subset = t_cell_score1<qn_b))
lc_bcell <- setdiff(colnames(bcell), hc_bcells)



#T cell
iqr_T <- IQR(tcell@meta.data$B_cell_score1)
qT<- quantile(tcell@meta.data$B_cell_score1, probs = 0.75)
qn_t <- qT + 1.5*iqr_T
hc_tcells <- colnames(subset(tcell, subset = B_cell_score1<qn_t))
lc_tcell <- setdiff(colnames(tcell), hc_tcells)

hc_list

#Get the low confidence cells from the hc_singlets

lc_singlets <- setdiff(hc_singlets, c(hc_tcells,hc_bcells))
lc_list[["other"]] = append(lc_list[["other"]],lc_singlets)
hc_singlets <- setdiff(hc_singlets, lc_singlets)






double_type_prediction<- function(pbmc, hc_doublets, hc_singlets){
  
  
  
  library("ggsci")
  library(Seurat)
  library(ggplot2)
  library(ggpubr)
  #set the Working Directory
  print("---------------Part1------------------Successful")
  #Calculate nUMI (this is the same as nCount_RNA which some datasets might include instead)
  nUMI = Matrix::colSums(pbmc@assays$RNA@counts)
  pbmc = AddMetaData(object = pbmc, metadata = as.data.frame(nUMI) ,col.name = "nUMI")
  
  #Separate the Singlet and Doublets
  pbmc_singlet <- subset(pbmc, cells = hc_singlets)
  pbmc_doublet <- subset(pbmc, cells = hc_doublets)
  pbmc_singlet
  ##Check to identify the proper names for each of the variables
  pbmc_singlet@meta.data$cell_refined_annotation <- make.names(pbmc_singlet@meta.data$cell_refined_annotation)
  colnames(pbmc_singlet@meta.data)
  colnames(pbmc_singlet@meta.data)[which(names(pbmc_singlet@meta.data)=="cell_refined_annotation")] <- "cell_type_overall"
  
  
  # get information from the main singlet dataset for simulating the doublets
  
  cell_ids = rownames(pbmc_singlet@meta.data)
  raw_data = pbmc_singlet@assays$RNA@counts 
  
  #Get gene names
  genes = rownames(raw_data)
  
  #Get original IDs i.e. just sample name, no sequence barcode
  orig.ident = pbmc_singlet@meta.data$orig.ident
  cell_type = pbmc_singlet@meta.data$cell_type_overall
  #Attach cell ids to cell type data and original IDs
  names(cell_type) = cell_ids
  names(orig.ident) = cell_ids
  
  cell_types = sort(unique(cell_type))
  cell_types
  
  clusters = cell_type #also clusters=cell_types
  
  orig.idents = sort(unique(orig.ident))
  #calculating the mito_ribo ratio of the model. I havent used it in the model but there is an interesting 
  #chance that I would start to do such soon in the model
  
  mitrib <- grep(pattern = "^MT-|^RP[SL]", x = rownames(x = raw_data), value = TRUE, ignore.case = TRUE)
  mito <- grep(pattern = "^MT-", x = rownames(x = raw_data), value = TRUE, ignore.case = TRUE)
  mitribcounts<- raw_data[which(rownames(raw_data) %in% mitrib), ]
  mitoribo_ratio <- Matrix::colSums(mitribcounts[mito, , drop = FALSE])/Matrix::colSums(mitribcounts)
  pbmc_singlet <- AddMetaData(object = pbmc_singlet, metadata = as.data.frame(mitoribo_ratio) ,col.name = "mito.ribo_ratio")
  
  #Collect the names of the Doublets and Singlets
  true_doublets = colnames(pbmc_doublet)
  true_singlets <- colnames(pbmc_singlet)
  #cell_ids
  
  #### Prepare for simulating the Doublets dataset
  #cell_types
  #Get the data for the cells we are using
  #This part I might delete later as I am now working with one sample each case 
  #This will lead to easy interpretation of the data and I will not be faced with excess coding
  true_doublets_sample = true_doublets
  true_singlets_sample = setdiff(cell_ids,true_doublets)
  
  use = true_singlets_sample[which(true_singlets_sample %in% names(cluster))]
  
  orig.ident1 = orig.ident[use]
  
  #Addition of the SHM data along with the model
  #using the v_mm1, v_mm2 for the as a feature for the machine learning model
  
  x <- m_VDJ_BCR  
  h_SHM = as.numeric(x[,"V_mm1"])
  l_SHM = as.numeric(x[,"V_mm2"])
  h_SHM[which(is.na(h_SHM))] = -1
  l_SHM[which(is.na(l_SHM))] = -1
  shm <- ifelse(h_SHM>3 | l_SHM>3,"high_SHM", "low_SHM")#check
  shm[intersect(which(h_SHM==-1), which(l_SHM==-1))] = NA
  names(shm) <- rownames(x)
  
  ### Actual Simulation of Doublets from the singlet data
  full_SHM <- NULL
  
  doublet_true_names = NULL
  doublet_expression = NULL
  doublet_sample_origin = NULL
  doublet_meta_data = NULL
  doublet_raw_data = NULL
  doublet_orig.idents = NULL
  
  triplet_true_names = NULL
  triplet_expression = NULL
  triplet_sample_origin = NULL
  triplet_meta_data = NULL
  triplet_raw_data = NULL
  triplet_orig.idents = NULL
  
  #clusters 
  #cell_types
  cluster1 <- clusters
  
  #Name of the cell types
  B_cell <- pbmc_singlet@meta.data[pbmc_singlet@meta.data$broad_cell_type == "B cell",]
  B_cell <- unique(B_cell$cell_type_overall)
  T_cell <- pbmc_singlet@meta.data[pbmc_singlet@meta.data$broad_cell_type == "T cell",]
  T_cell <- unique(T_cell$cell_type_overall)
  
  clusters1<- T_cell
  clusters2 <- B_cell
  
  w_choose_sample = cell_ids
  
  for(c1 in c(1:length(clusters1))){
    for(c2 in c(1:length(clusters2))){
      clust1 = clusters1[c1]
      clust2 = clusters2[c2]
      
      w1 = intersect(names(which(cluster1==clust1)), w_choose_sample)
      w2 = intersect(names(which(cluster1==clust2)), w_choose_sample)
      
      rand = cbind(sample(w1, 200, replace = T), sample(w2, 200, replace = T))
      rand = unique(rand) #Remove any duplicates
      
      #Get a list of which samples the doublets are from
      ident = apply(cbind(rand[,1], rand[,2]), 1, paste, collapse = '..')
      doublet_sample_origin = c(doublet_sample_origin, ident)
      print(paste(clust1, clust2, sep = ".."))
      
      #Code for SHM information adding
      a<- rand[,2]
      each_SHM <- shm[a]
      names(each_SHM) <- ident 
      full_SHM <- c(full_SHM, each_SHM)
      
      #Orig.idents  
      orig_ident = orig.ident[rand[,1]]
      #print(head(orig.ident))
      names(orig_ident) = ident
      #orig.ident
      
      doublet_orig.idents = c(doublet_orig.idents, orig_ident) #Add each loop on to the vector 
      #print(head(doublet_orig.idents))
      
      #Get cell expression values for the randomised cell list for each cluster
      #combine the cell expression proportions 
      #with alpha% coming from one cluster and 1-alpha% from the other cluster
      
      simul_raw_data = (raw_data[,rand[,1]])+ (raw_data[,rand[,2]])
      colnames(simul_raw_data) = ident
      
      variables = pbmc_singlet@meta.data[,c("percent.mito", "nUMI","mito.ribo_ratio")]
      
      #Calculate metadata for each doublet
      variables_doublets = data.frame(matrix(nrow=nrow(rand),ncol=0))
      variables_doublets$nUMI = (variables[rand[,1],]$nUMI)+ (variables[rand[,2],]$nUMI)
      rownames(variables_doublets) = ident
      #head(variables_doublets)
      
      #Get a list with the cell type combination of each doublet cell
      t_names = apply(cbind(cluster1[rand[,1]], cluster1[rand[,2]]), 1, paste, collapse = '..')
      names(t_names) = ident
      doublet_true_names = c(doublet_true_names, t_names) #Add each loop on to the vector 
      
      #doublet,triplet or quadruplet?
      type_doublet = rep('doublet',length(doublet_sample_origin))
      names(type_doublet) = doublet_sample_origin
      
      #Put together doublet meta data
      if(length(doublet_meta_data)==0){doublet_meta_data = variables_doublets #Add the first loop as the first entry
      }else{doublet_meta_data = rbind(doublet_meta_data ,variables_doublets)} #rbind subsequent loops to the first loop
      
      
      #Put together the doublet raw data
      if(length(doublet_raw_data)==0){doublet_raw_data = simul_raw_data #Add the first loop as the first entry
      }else{doublet_raw_data = cbind(doublet_raw_data ,simul_raw_data)} #rbind subsequent loops to the first loop
      print (dim(doublet_raw_data)) #Print the dimensions of the new expression matrix after every loop
    }
  }
  ##Find out how many cell_types we need to simulate for the triplets to equal the doublets
  library(gtools)
  cell_no<- round(ncol(doublet_raw_data)/nrow(combinations(length(cell_types),3, cell_types, repeats = TRUE)))
  
  
  
  for(c1 in c(1:length(cell_types))){
    for(c2 in c(c1:length(cell_types))){
      for(c3 in c(c2:length(cell_types))){
        if(c1<=c2 & c2<=c3){
          
          clust1 = cell_types[c1]
          clust2 = cell_types[c2]
          clust3 = cell_types[c3]
          
          w1 = intersect(names(which(cluster1==clust1)), w_choose_sample)
          w2 = intersect(names(which(cluster1==clust2)), w_choose_sample)
          w3 = intersect(names(which(cluster1==clust3)), w_choose_sample)
          
          ### randomly sample XX cells from each cluster
          rand = cbind(sample(w1, cell_no, replace = T), sample(w2,cell_no, replace = T), sample(w3,cell_no, replace = T))
          rand = unique(rand) #Remove any duplicates
          
          #Get a list of which samples the doublets are from
          ident = apply(cbind(rand[,1], rand[,2], rand[,3]), 1, paste, collapse = '..')
          triplet_sample_origin = c(triplet_sample_origin, ident)
          print(paste(clust1, clust2, clust3, sep =".."))
          
          #Orig.idents  
          orig_ident = orig.ident[rand[,1]]
          names(orig_ident) = ident
          triplet_orig.idents = c(triplet_orig.idents, orig_ident) #Add each loop on to the vector 
          
          
          #Get cell expression values for the randomised cell list for each cluster
          #combine the cell expression proportions 
          #simulated = data.frame(((1/3)*cell_exp_matrix[,rand[,1]])+ ((1/3)*cell_exp_matrix[,rand[,2]])+ ((1/3)*cell_exp_matrix[,rand[,3]]))
          
          simul_raw_data = (raw_data[,rand[,1]])+ (raw_data[,rand[,2]])+ (raw_data[,rand[,3]])
          colnames(simul_raw_data) = ident
          
          variables = pbmc_singlet@meta.data[,c("percent.mito", "nUMI","mito.ribo_ratio")]
          #Calculate metadata for each triplet      
          
          variables_triplets = data.frame(matrix(nrow=nrow(rand),ncol=0))
          variables_triplets$nUMI = (variables[rand[,1],]$nUMI)+ (variables[rand[,2],]$nUMI)+ (variables[rand[,3],]$nUMI)
          rownames(variables_triplets) = ident
          
          #Get a list with the cell type combination of each doublet cell
          t_names = apply(cbind(cluster1[rand[,1]], cluster1[rand[,2]], cluster1[rand[,3]]), 1, paste, collapse = '--')
          names(t_names) = ident
          triplet_true_names = c(triplet_true_names, t_names) #Add each loop on to the vector 
          
          type_triplet = rep('triplet',length(triplet_sample_origin))
          names(type_triplet) = triplet_sample_origin
          
          #Put together doublet meta data
          if(length(triplet_meta_data)==0){triplet_meta_data = variables_triplets #Add the first loop as the first entry
          }else{triplet_meta_data = rbind(triplet_meta_data ,variables_triplets)} #rbind subsequent loops to the first loop
          
          #Put together the doublet raw data
          if(length(triplet_raw_data)==0){triplet_raw_data = simul_raw_data #Add the first loop as the first entry
          }else{triplet_raw_data = cbind(triplet_raw_data ,simul_raw_data)} #rbind subsequent loops to the first loop
          print (dim(triplet_raw_data)) #Print the dimensions of the new expression matrix after every loop
        }}}}
  
  
  ncol(doublet_raw_data)
  ncol(triplet_raw_data)
  print("---------------Part2------------------Successful")
  #Put everything together
  multiplet_true_names = c(doublet_true_names,triplet_true_names)#,quadruplet_true_names)
  multiplet_raw_data = cbind(doublet_raw_data, triplet_raw_data)
  multiplet_meta_data = rbind(doublet_meta_data,triplet_meta_data)
  multiplet_sample_origin = c(doublet_sample_origin,triplet_sample_origin)#,quadruplet_sample_origin)
  multiplet_type = c(type_doublet,type_triplet)#,type_quadruplet)
  multiplet_orig.idents = c(doublet_orig.idents,triplet_orig.idents)
  
  
  
  
  ## Calculate the mitoribo ratio of the simulated doublets
  #Calculate mito.ribo_ratio
  #Calculate mito.ribo_ratio
  mitrib <- grep(pattern = "^MT-|^RP[SL]", x = rownames(x = multiplet_raw_data), value = TRUE, ignore.case = TRUE)
  mito <- grep(pattern = "^MT-", x = rownames(x = multiplet_raw_data), value = TRUE, ignore.case = TRUE)
  mitribcounts<- multiplet_raw_data[which(rownames(multiplet_raw_data) %in% mitrib), ]
  mitoribo_ratio <- Matrix::colSums(mitribcounts[mito, , drop = FALSE])/Matrix::colSums(mitribcounts)
  multiplet_meta_data$mito.ribo_ratio <- mitoribo_ratio
  
  simulated_multiplets = c(list(multiplet_true_names),  list(multiplet_raw_data), 
                           list(multiplet_meta_data), list(multiplet_sample_origin), list(multiplet_type),list(multiplet_orig.idents))
  names(simulated_multiplets) = c("multiplet_true_names", "multiplet_raw_data", 
                                  "multiplet_meta_data", "multiplet_sample_origin", "multiplet_type","multiplet_orig.idents")
  
  
  ##Store the expression data and other meta data in one object for future references
  
  
  ### Part 2
  
  
  #Load the real Doublet data and add the simulated and real doublet data together
  
  #loading
  raw_data_doublet = pbmc_doublet@assays$RNA@counts
  cell_ids_doublet = rownames(pbmc_doublet@meta.data)
  
  
  # Start with merging the simulated and raw data
  
  #Merge raw data
  simulated_doublet_raw_data = simulated_multiplets$multiplet_raw_data
  true_doublet_raw_data = raw_data_doublet
  doublet_raw_data_total = cbind(simulated_doublet_raw_data, true_doublet_raw_data)
  
  #Adding meta data of simulated and raw doublet data together
  simulated_doublet_meta_data = simulated_multiplets$multiplet_meta_data
  
  
  #adding the mito ribo ratio to the doublet dataset
  mitrib_d <- grep(pattern = "^MT-|^RP[SL]", x = rownames(x = raw_data_doublet), value = TRUE, ignore.case = TRUE)
  mito_d <- grep(pattern = "^MT-", x = rownames(x = raw_data_doublet), value = TRUE, ignore.case = TRUE)
  mitribcounts_d<- raw_data_doublet[which(rownames(raw_data_doublet) %in% mitrib_d), ]
  mitoribo_ratio_d <- Matrix::colSums(mitribcounts_d[mito_d, , drop = FALSE])/Matrix::colSums(mitribcounts_d)
  pbmc_doublet <- AddMetaData(object = pbmc_doublet, metadata = as.data.frame(mitoribo_ratio_d) ,col.name = "mito.ribo_ratio")
  
  # **** percent.mito in place of percent.mt
  variables_doublet <- pbmc_doublet@meta.data[,c("percent.mito", "nUMI","mito.ribo_ratio")]
  true_doublet_meta_data = variables_doublet[true_doublets_sample, which(colnames(variables_doublet) %in% colnames(simulated_doublet_meta_data))]
  doublet_meta_data_total = rbind(simulated_doublet_meta_data,true_doublet_meta_data)
  
  
  #Add label to say whether true doublet
  names_true_doublets = rep("true doublet", length(true_doublets_sample))
  names(names_true_doublets) = true_doublets_sample
  
  #Get all names and origins
  true_doublets_sample <- pbmc_doublet@meta.data$orig.ident
  names(true_doublets_sample) <- true_doublets
  
  doublet_true_names_all = append(simulated_multiplets$multiplet_true_names, names_true_doublets)
  names(doublet_true_names_all) = rownames(doublet_meta_data_total)
  
  doublet_orig.ident_all = c(simulated_multiplets$multiplet_orig.idents, true_doublets_sample)
  names(doublet_orig.ident_all) = names(doublet_true_names_all)
  
  multiplet_types_all = c(simulated_multiplets$multiplet_type,names_true_doublets)
  names(multiplet_types_all) <- names(doublet_true_names_all)
  
  
  #Put it all together and save
  doublet_sample = c(list(doublet_raw_data_total),
                     list(doublet_true_names_all),list(multiplet_types_all),list(doublet_meta_data_total),list(doublet_orig.ident_all))
  
  
  #Renaming the same doublet objects again. I can change it after I am done with the work first
  #Here it is doing the same thing, just the variables have been renamed.
  
  doublet_raw_data = doublet_sample[[1]]
  doublet_true_names_all = doublet_sample[[2]]
  multiplet_types_all = doublet_sample[[3]]
  doublet_meta_data = doublet_sample[[4]]
  doublet_orig.idents = doublet_sample[[5]]
  
  
  #Get groups- true doublets and singlets, and simulated doublets and multiplets.
  singlets = colnames(pbmc_singlet)[!colnames(pbmc_singlet) %in% true_doublets_sample]
  doublets = names(multiplet_types_all)[multiplet_types_all=='true doublet']
  doublets_sim = names(multiplet_types_all)[multiplet_types_all=='doublet']
  multiplets_sim = names(multiplet_types_all)[multiplet_types_all=='triplet']
  #combine singlet and simulated raw data
  raw_data_all = cbind(raw_data[,singlets],doublet_raw_data)
  print("---------------Part3------------------Successful")
  
  #### Part - 3
  # Creating Seurat downstream processing on the data
  
  # Move the data into "Seurat Object"
  data <-  CreateSeuratObject(counts = raw_data_all, project = "all")
  
  #There can be additional step regarding Quality control and filtering 
  #Although this might be troublesome due to we want the extra feature containing cells
  
  # Normalisation
  data <- NormalizeData(data)
  
  #Dimensionality Reduction
  #We start with finding highly variable genes. But the interesting part is that we are not using the
  # variable genes from this dataset. Rather we are using the variable genes from the main pbmc dataset
  #here in this case it is the singlet dataset
  
  data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
  
  #ALso we want to exclude genes of cell cycle, TCR, BCR , MT and RP specific genes from the analysis. So that they do not bias the result
  
  vf = data@ assays$ RNA@ var.features
  
  genes = rownames(data)
  
  IGHV = genes [grep("IGHV",genes)]
  
  IGKLV = c(genes [grep("IGKV",genes)], genes [grep("IGLV",genes)])
  
  variable = c(genes [grep("^MT",genes)], genes [grep("^RP",genes)])
  
  TCG = c(genes [grep("^TRA",genes)], genes [grep("^TRB",genes)],genes [grep("^TRD",genes)],genes [grep("^TRG",genes)])
  
  cell_cycle =c("STMN1","CDC20","ANP32E","CKS1B","NUF2","ASPM","UBE2T","CENPF","RRM2","SGOL2","SGOL1",
                
                "SMC4","CENPE","CCNA2","CCNB1","PTTG1","HMMR","MXD3","HIST1H4C","CENPW","H2AFV","CKS2","SMC2",
                
                "ZWINT","CDK1","MKI67","C12orf75","CDKN3","NUSAP1","CCNB2","KIAA0101","PRC1","ARL6IP1","PLK1",
                
                "CENPN","AURKB","TOP2A","KPNA2","HN1","TK1","BIRC5","TYMS","TPX2","UBE2C","CALM3","UBE2S","GTSE1",
                
                "CLSPN","CDCA7","CENPU","GMNN","MCM7","CCDC34","CDCA5","HELLS","TMEM106C","IDH2","SNRNP25","CDT1","PCNA",
                
                "UHRF1","ASF1B")
  
  noisy = c("MALAT1", "JCHAIN", "XIST")
  
  variable = c(variable, cell_cycle, noisy)
  
  #Then exclude these from your variable gene list:
  
  added_genes <- c("CD8A","CD8B","CD4","CD27","IGHD","CCR7")
  data@ assays$ RNA@ var.features=unique(c(data@ assays$ RNA@ var.features,added_genes))
  
  vf = setdiff(data@ assays$ RNA@ var.features, c(variable, IGHV, IGKLV, TCG))
  #sort for repeat
  vf = sort(unique(vf))
  
  ## Keep only the variable genes from the new created dataset.
  data = data[vf,]
  
  #Add labels for singlets, true doublets and simulated doublets
  class <- c(rep('singlet',length(singlets)),multiplet_types_all)
  
  names(class) = colnames(raw_data_all)
  
  #Adding class to meta data of data seurat object
  
  data$class = class
  
  Idents(data) <- "class"
  
  # Follow through with the normal seurat procedure
  
  #Finding Highly Variable genes
  data <- FindVariableFeatures(data,
                               selection.method = "vst",
                               nfeatures = 2000)
  
  #Principal Component Analysis
  
  data <- ScaleData(data)
  
  #Run PCA and set the number of PC you want to include in the machine learning model and see how many PC explains the most Variability
  print("---------------Part4------------------Successful")
  data <- RunPCA(data,approx=F, npc = 50)
  
  #Run the UMAP
  data <- RunUMAP(data, dims = 1:20)
  
  # We now want to check how the UMAP projection looks like compared to the simulated vs True Doublet dataset
  #Check the easy exploration of the primary source of heterogeneity in a dataset
  pdf("heatmap for PCS.pdf")
  p<- DimHeatmap(data, dims = 1:20, cells = 500, balanced = TRUE)
  print(p)
  dev.off()
  #First we need to subset the dataset to only include the simulated and true doublet and remove the singlets from the dataset
  
  data1 = subset(data,cells=names(doublet_true_names_all)) 
  #create a new column with the information of simulated and true doublets dataset
  data1$cell_type = doublet_true_names_all
  Idents(data1) = 'cell_type'
  
  
  ## Organize the Meta data for the singlet, simulated doublets and true doublets
  #check the nUMI count for the dataset
  
  ##There are two points here to make. One is to use log of nCountRNA or non log
  #Another is to use the raw data or to use the normalised gene removed data
  
  ## Part 4
  #This works in two ways. One would be to remove the multiplets from doublets
  
  ## Creating Machine Learning Model for Analysis
  
  ## Creating the dataframe required for the machine learning model
  ## Subset the singlets from the dataset
  
  
  data_subset <- subset(data, subset = class !="singlet")
  
  
  
  # Prepare the dataframe for the machine learning model
  
  data2 <- as.data.frame(Embeddings(data_subset, reduction = "pca")[,1:20])
  
  ## From here there are two things to do
  ## First, We wanted to see whether there are additional features that we can use in along with the PC for the model development
  ## So added few more columns with the PC containing data
  
  #Adding additional columns for the machine learning model
  data2$nCount_RNA <- data_subset@meta.data$nCount_RNA
  data2$nFeature_RNA <- data_subset@meta.data$nFeature_RNA
  data2$nFeature_RNA <- as.numeric(data2$nFeature_RNA)
  #Use specific genes for the model
  data3 <- as.data.frame(GetAssayData(object = data_subset, slot = "counts"))
  
  added_genes <- c("CD8A","CD8B","CD4","CD27","IGHD","CCR7")
  genes_to_add <- intersect(added_genes, rownames(data3))
  print(genes_to_add)
  data3 <- data3[genes_to_add,]
  data3 <- t(data3)
  
  #Adding the columns to the main data
  data2 <- cbind(data2, data3)
  #removed some of the columns as not needed or would bias the model
  data2$nFeature_RNA <- NULL
  data2$nCount_RNA <- NULL
  #data2$CD4 <- NULL
  #data2$CD8A <- NULL
  #data2$CD8B <- NULL
  #Add class column
  #Check if the columns are right
  table(rownames(data_subset@meta.data) == rownames(data2))
  data2$class <- data_subset@meta.data$class
  nrow(data2)
  nrow(data_subset@meta.data)
  
  table(data2$class)
  
  real_doublet_data <- data2[data2$class =="true doublet",]
  data2 <- data2[data2$class !="true doublet",]
  class <- data2$class
  data4 <-  data2
  data4$class <- NULL
  #install.packages("caret")
  
  library(caret)
  testindex <- createDataPartition(class, p = .30, list = F)
  
  #testindex <- sort(sample(c(1:length(w_use)), trunc(length(w_use)/3))) #Random sample for test set
  testset   <- data4[testindex,] #Get test set
  trainset  <- data4[-testindex,] #Get training set
  
  class_test = class[testindex] #Get cell type names of training cells
  class_train= class[-testindex] #Get cell type names of test cells
  
  #k-fold cross validation for random forest model
  trainset$class_train <- class_train
  testset$class_test <- class_test
  trcontrol = trainControl(method='cv', number=10, savePredictions = T,
                           classProbs = TRUE,returnResamp="all")
  library(randomForest)
  library(caret)
  library(pROC)
  print("---------------Part5------------------Successful")
  model = train(class_train ~ ., data=trainset, method = "rf", trControl = trcontrol) 
  
  model$resample
  
  #Predict the test set
  
  predict_prob_rf_cv <- predict(model, testset,type='prob')
  #Get lists of predicitons and models
  #predict_probs<-list(predict_prob_rf_cv)
  fits<-model
  #Check the AUC of the model
  pred1 = colnames(predict_prob_rf_cv) [max.col(predict_prob_rf_cv, 'first')]
  names(pred1) = names(class_test)
  roc.multi = multiclass.roc(class_test,as.numeric(as.factor(pred1)))
  print(paste("The AUC for Multiplet", roc.multi$auc))
  print(roc.multi$auc)
  conMat_Multiplet = confusionMatrix(as.factor(pred1),as.factor(class_test))
  print("Confusion Matrix for Multiplet Detection")
  print(conMat_Multiplet)
  pro_multiplet <- prop.table(conMat_Multiplet$table,2)*100 
  
  #Extract true doublets
  #xy_simul1 = xy_simul[which(doublet_true_names_all=="true doublet" ),] 
  
  name_true <- real_doublet_data[,"class"]
  xy_simul1 = real_doublet_data
  xy_simul1$class <- NULL
  #Identify true doublets from multiplets
  
  Identify_doublets<-function(xy_simul1, fits){
    library(e1071)
    data = as.data.frame(xy_simul1)
    names(data) = make.names(names(data))
    
    fit = fits #Extract the model
    
    predict_prob <- predict(fit, data,type='prob') #Predict the true doublet identities
    
    rf.pred = colnames(predict_prob) [max.col(predict_prob, 'first')]
    prob_max = apply(predict_prob, 1,max)
    
    names(rf.pred) <- rownames(xy_simul1)
    
    predicted_interactions_doublets = c(list(predict_prob), list(rf.pred), list(prob_max))
    names(predicted_interactions_doublets) = c("predict_prob", "rf.pred","prob_max")
    return(predicted_interactions_doublets)
    
  }
  
  #Identify the true doublets
  identify_true_doublets = Identify_doublets(xy_simul1, fits) 
  
  #Compare true and simulated for each patient using DGE
  
  #Get predictions
  preds_sim_ml <- pred1
  preds_true_ml <- identify_true_doublets$rf.pred
  preds_ml = c(preds_sim_ml,preds_true_ml)
  
  #Get simulated and true names
  sim = names(preds_sim_ml)
  true = names(preds_true_ml)
  
  #Remove totalC cells from doublet_raw_data
  doublet_raw_data1 = doublet_raw_data[!grepl("TotalC|TotalSeqC", rownames(doublet_raw_data)),]
  
  
  #preparing the predicted doublet type data into dataframe
  #make the preds_true vector a dataframe and add a cell with the name
  preds_true_ml <- as.data.frame(preds_true_ml)
  preds_true_ml$cells <-rownames(preds_true_ml)
  #Changing the factor data to character variable
  preds_true_ml$preds_true_ml <- as.character(preds_true_ml$preds_true_ml)
  
  table(preds_true_ml$preds_true_ml)
  
  final_doublet_cell <- rownames(preds_true_ml[preds_true_ml$preds_true_ml =="doublet",]) 
  final_triplet_cell <- rownames(preds_true_ml[preds_true_ml$preds_true_ml =="triplet",]) 
  hc_list[["heterotypic_doublets"]] <- final_doublet_cell
  hc_list[["other"]] <- final_triplet_cell
  hc_list <<- hc_list
  
  if(length(hc_list[["heterotypic_doublets"]]) == 0){
    print("No High Confidence Doublet in the dataset to work with")
  }
  #Outlier Detection in the Doublet Dataset
  
  ##Finding outliers in the Data
  #First get the doublet that we get after removal of triplet from the data
  
  pbmc_doublet <- subset(pbmc_doublet, cells = final_doublet_cell)
  
  
  x<- pbmc_doublet@meta.data
  library(tidyr)
  colnames(x)
  x$cell_id <- rownames(x)
  
  m<- x[,c("nCount_RNA","B_Tcell","cell_id")]
  head(m)
  rownames(m) <- m$cell_id
  y<- spread(m, key =  B_Tcell, value = nCount_RNA)
  
  outlierreplacement <- function(dataframe){
    dataframe %>%          
      map_if(is.numeric, ~ replace(.x, .x %in% boxplot.stats(.x)$out, NA)) %>%
      bind_cols 
    
    
    
  }
  k<- outlierreplacement(y)
  z<- gather(k, key = "B/T/D_cell", value = "nCount_RNA",2:ncol(k) )
  ##get the final doublets from the dataset
  
  final_doublet <- z[!is.na(z$nCount_RNA),]$cell_id
  final_triplet <- z[is.na(z$nCount_RNA),]$cell_id
  
  #Doublet outlier Figure
  
  x$outlier <- "Before"
  y<- x[final_doublet,]
  y$outlier <- "After"
  final_file <- rbind(x,y)
  final_file$outlier <- factor(final_file$outlier, levels = c("Before","After"))
  pdf("Doublet outlier.pdf", height = 3, width = 3)
  
  p <- ggplot(final_file, aes(x = outlier, y = nCount_RNA, fill = outlier)) +
    geom_boxplot()+
    theme_bw()+
    scale_fill_npg()+
    theme(axis.text.x = element_text(vjust = 0.5),
          axis.text.y = element_text(hjust = 0.5)) +
    labs(x = "Outlier", y = "nUMI")+
    scale_fill_manual(values = cols)
  print(p)
  dev.off()
  
  hc_list[["heterotypic_doublets"]] <- final_doublet
  hc_list[["other"]] <- append(hc_list[["other"]],final_triplet )
  
  hc_list <<- hc_list
  
  pbmc_doublet <- subset(pbmc_doublet, cells = final_doublet)
  #Now that we have found the Triplets in the Dataset. Now the final work would be to identify
  #the doublets and use them to predict the interaction between two cell types
  
  #Now run the same machine learning model for the doublet dataset 
  
  # Prepare the dataframe for the machine learning model
  
  data2 <- as.data.frame(Embeddings(data_subset, reduction = "pca")[,1:20])
  
  ## From here there are two things to do
  ## First, We wanted to see whether there are additional features that we can use in along with the PC for the model development
  ## So added few more columns with the PC containing data
  
  #Adding additional columns for the machine learning model
  data2$nCount_RNA <- data_subset@meta.data$nCount_RNA
  data2$nFeature_RNA <- data_subset@meta.data$nFeature_RNA
  data2$nFeature_RNA <- as.numeric(data2$nFeature_RNA)
  #Use specific genes for the model
  data3 <- as.data.frame(GetAssayData(object = data_subset, slot = "counts"))
  data3 <- data3[genes_to_add,]
  data3 <- t(data3)
  
  #Adding the columns to the main data
  data2 <- cbind(data2, data3)
  
  #Get doublet cell type names
  true_names_simul<-as.numeric(as.factor(doublet_true_names_all))
  
  n_type_doublet = length(unique(doublet_true_names_all[which(multiplet_types_all=='doublet')]))
  w_use = which(doublet_true_names_all %in% unique(doublet_true_names_all)[c(1:n_type_doublet)])
  
  #First we create the training and test dataset. This will be useful for later actual model building as well.
  data4 <- data2[w_use,]
  
  #Added SHM column
  data4$SHM <- as.factor(full_SHM)
  if(any(is.na(data4$SHM)) == TRUE) {
    position <- which(is.na(data4$SHM))
    class = factor(doublet_true_names_all[w_use])
    class = class[-position]
    
    data4 <- na.omit(data4)
    print("NA value was there and was removed")
  }else{
    #First we create the training and test dataset. This will be useful for later actual model building as well.
    data4 <- data2[w_use,]
    
    #Added SHM column
    data4$SHM <- as.factor(full_SHM)
    
    #data = as.data.frame(xy_simul[w_use,]) #Get simulated doublets
    class = factor(doublet_true_names_all[w_use]) #Get cell type combination names of simulated doublets
    #names(data2) = make.names(names(data2)) #Make cell type names
    print("No NA value was there")
  }
  
  #names(data2) = make.names(names(data2)) #Make cell type names
  
  testindex <- createDataPartition(class, p = .30, list = F)
  
  #testindex <- sort(sample(c(1:length(w_use)), trunc(length(w_use)/3))) #Random sample for test set
  testset   <- data4[testindex,] #Get test set
  trainset  <- data4[-testindex,] #Get training set
  
  class_test = class[testindex] #Get cell type names of training cells
  class_train= class[-testindex] #Get cell type names of test cells
  ## From here we are separating into two parts. One part is to do the recursive feature selection and
  ## Second part is doing the actual machine learning model.
  # Define the control using a random forest selection function
  print("It is doing Recursive feature Selection and takes a long time")
  control <- rfeControl(functions = rfFuncs, # random forest
                        method = "repeatedcv", # repeated cv
                        repeats = 5, # number of repeats
                        number = 10) # number of folds

  #
  # # Run RFE
  result_rfe2 <- rfe(x = trainset,
                     y = class_train,
                     sizes = c(1:5), # 20 features in total
                     rfeControl = control)

  # # Print the results
  result_rfe2

  # Variable importance
  varimp_data <- data.frame(feature = row.names(varImp(result_rfe2))[1:ncol(trainset)],
                            importance = varImp(result_rfe2)[1:ncol(trainset), 1])

  pdf("Top_features for model.pdf", height = 5, width =5)
  p<- ggplot(data = varimp_data,
             aes(x = reorder(feature, -importance), y = importance, fill = feature)) +
    geom_bar(stat="identity") + labs(x = "Features", y = "Variable Importance") +
    geom_text(aes(label = round(importance, 1)), vjust=1.6, color="black", size =1.5) +
    theme_bw() + theme(legend.position = "none")+
    pdf("Top_features for model.pdf", height = 7, width =7)
  p<- ggplot(data = varimp_data,
             aes(x = reorder(feature, -importance), y = importance, fill = feature)) +
    geom_bar(stat="identity") + labs(x = "Features", y = "Variable Importance") +
    #geom_text(aes(label = round(importance, 1)), vjust=1.6, color="black", size =2) +
    theme_bw() + theme(legend.position = "none")+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5,
                                     size = 16,
                                     hjust = 1),
          axis.text.y = element_text(vjust = 0.5,
                                     hjust = 1,
                                     size = 16),
          legend.text = element_text(size=8))+
    scale_fill_manual(values = my_colors)
  
  print(p)
  dev.off()

  my_colors <- c(
    "#FF5733",  # Red
             "#3366FF",  # Blue
             "#33FF33",  # Green
             "#000000",  # Black
             "#FFA500",  # Orange
             "#800080",  # Purple
             "#FFFF00",  # Yellow
             "#FF1493",  # Pink
             "#8B4513",  # Brown
             "#00CED1",  # Cyan
             "#008080",  # Teal
             "#808080",  # Gray
             "#E6E6FA",  # Lavender
             "#FF00FF",  # Magenta
             "#808000",  # Olive
             "#800000",  # Maroon
             "#FFD700",  # Gold
             "#C0C0C0",  # Silver
             "#4B0082",  # Indigo
             "#FF6F61",  # Coral
             "#40E0D0",  # Turquoise
             "#DDA0DD",  # Plum
             "#F5F5DC",  # Beige
             "#F0E68C",  # Khaki
             "#EE82EE",  # Violet
             "#708090",  # Slate
             "#36454F",  # Charcoal
             "#AA4069",  # Ruby
             "#F0FFFF",  # Azure
             "#98FF98"   # Mint
  )
  
  #k-fold cross validation for random forest model
  trainset$class_train <- class_train
  testset$class_test <- class_test
  trcontrol = trainControl(method='cv', number=10, savePredictions = T,
                           classProbs = TRUE,returnResamp="all")
  
  print("---------------Part6------------------Successful")
  model = train(class_train ~ ., data=trainset, method = "rf", trControl = trcontrol) 
  
  model$resample
  
  #Predict the test set
  
  predict_prob_rf_cv <- predict(model, testset,type='prob')
  #Get lists of predicitons and models
  predict_probs<-list(predict_prob_rf_cv)
  fits<-list(model)
  
  
  #Set up storage for loop
  predicted_interactions = vector('list',length=1)
  names(predicted_interactions)<-c('rf')
  library(caret)
  library(pROC)
  #Extract results of model testing
  for ( i in c(1:length(predict_probs))){
    
    predict_prob = predict_probs[[i]]
    fit = fits[[i]]
    predict_prob
    
    pred = colnames(predict_prob) [max.col(predict_prob, 'first')] #Extract the predictions
    names(pred) = names(class_test)
    roc.multi = multiclass.roc(class_test,as.numeric(as.factor(pred)))
    conMat = confusionMatrix(as.factor(pred),class_test)
    table_correct = table(pred, class_test) #Not essential, this is in confusion matrix
    percentage_correct = length(which(pred== class_test))*100/length(class_test) #What proportion the svm got right
    prob_max = apply(predict_prob, 1,max) #Get probability for each cell
    w = which(prob_max>0.5) #Get those with probability greater than 0.5
    percentage_correct_high_prob = length(which(pred[w]== class_test[w]))*100/length(w) #The percentage that are correct with high probability
    percentage_high_prob = length(w)*100/length(class_test) #Percentage of those with a high probability
    
    predicted_interactions[[i]] = c(list(fit), list(table_correct), list(percentage_correct), 
                                    list(prob_max), list(percentage_correct_high_prob), 
                                    list(percentage_high_prob),list(pred),list(roc.multi),list(conMat))
    
    names(predicted_interactions[[i]]) = c('model', "table_correct","percentage_correct","prob_max", "percentage_correct_high_prob", 
                                           "percentage_high_prob","pred",'roc.multi','confusion_matrix')
    
  }
  
  print(predicted_interactions)
  #Extract stats for comparison
  predicted_interactions_compare<-data.frame(model=c('rf'),
                                             AUC=c(predicted_interactions$rf$roc.multi$auc),
                                             percentage_correct=c(predicted_interactions$rf$percentage_correct),
                                             percentage_correct_high_prob=c(predicted_interactions$rf$percentage_correct_high_prob),
                                             sensitivity=c(mean(predicted_interactions$rf$confusion_matrix$byClass[,1])),
                                             specificity=c(mean(predicted_interactions$rf$confusion_matrix$byClass[,2])))
  
  print(predicted_interactions_compare)
  
  write.csv(file='predicted_interactions_compare_multiplets_removed.csv',predicted_interactions_compare)
  
  #Get the stats on a class level
  library(tidyverse)
  sens_spec_by_class<-data.frame(sensitivity=cbind(predicted_interactions$rf$confusion_matrix$byClass[,1]),
                                 specificity=cbind(predicted_interactions$rf$confusion_matrix$byClass[,2]))
  
  rownames(sens_spec_by_class)<-gsub("Class: ","",rownames(sens_spec_by_class))
  sens_spec_by_class<-sens_spec_by_class %>% rownames_to_column('cell_type')
  sens_spec_by_class$cell_type<-as.factor(sens_spec_by_class$cell_type)
  
  print(sens_spec_by_class)
  
  #Random forest stats
  rf_stats_by_class <- predicted_interactions$rf$confusion_matrix$byClass
  
  write.csv(rf_stats_by_class,'rf_stats_by_class.csv')
  
  #Get prediction accuracy by class
  accuracy_by_class_rf<-prop.table(predicted_interactions$rf$confusion_matrix$table,2)*100
  write.csv(accuracy_by_class_rf,'accuracy_by_class_rf_lung.csv')
  
  
  xy_simul1 = data2[final_doublet,]
  
  xy_simul1$SHM <- as.factor(shm[ which(names(shm)%in% rownames(xy_simul1))])
  
  #Classify true doublets
  
  Classify_doublets<-function(xy_simul1, predicted_interactions){
    library(e1071)
    data = as.data.frame(xy_simul1)
    names(data) = make.names(names(data))
    
    fit = predicted_interactions$rf$model #Extract the model
    
    predict_prob <- predict(fit, data,type='prob') #Predict the true doublet identities
    
    rf.pred = colnames(predict_prob) [max.col(predict_prob, 'first')]
    prob_max = apply(predict_prob, 1,max)
    
    names(rf.pred) <- rownames(xy_simul1)
    
    predicted_interactions_doublets = c(list(predict_prob), list(rf.pred), list(prob_max))
    names(predicted_interactions_doublets) = c("predict_prob", "rf.pred","prob_max")
    return(predicted_interactions_doublets)
    
  }
  #Classify the true doublets
  predicted_interactions_doublets = Classify_doublets(xy_simul1, predicted_interactions) 
  
  #Save doublet predictions
  prediction_doublets_overall_rf = c(list(predicted_interactions), list(predicted_interactions_doublets))
  names(prediction_doublets_overall_rf) = c("predicted_interactions","predicted_interactions_doublets")
  
  
  #Read in the doublet predictions
  predicted_doublet_types <- prediction_doublets_overall_rf
  
  #Get predictions
  preds_sim <- predicted_doublet_types$predicted_interactions$rf$pred
  preds_true <- predicted_doublet_types$predicted_interactions_doublets$rf.pred
  preds = c(preds_sim,preds_true)
  
  #Get simulated and true names
  sim = names(preds_sim)
  true = names(preds_true)
  
  
  #preparing the predicted doublet type data into dataframe
  #make the preds_true vector a dataframe and add a cell with the name
  preds_true <- as.data.frame(preds_true)
  preds_true$cells <-rownames(preds_true)
  #Changing the factor data to character variable
  preds_true$preds_true <- as.character(preds_true$preds_true)
  ## Get the hc list and lc list and name them
  write.csv(preds_true, "doublet_types_with_cell_name.csv")
  ## Part 5
  
  ## Check cellular interaction for the model and also differential gene expression plus the pathways that are active.
  
  #Create seurat object
  pbmc_doub <- CreateSeuratObject(counts = doublet_raw_data1[,names(preds)], project = "all")
  
  #Label whether simulated or true doublet
  type =c(rep('sim',length(sim)),rep('true',length(true)))
  names(type) = names(preds)
  
  #Get original IDs
  orig.idents_doub = doublet_orig.idents[names(preds)]
  
  #Add meta data
  pbmc_doub = AddMetaData(pbmc_doub,type,'type')
  pbmc_doub = AddMetaData(pbmc_doub,orig.idents_doub,'orig.ident')
  
  
  
  
  #Change Idents
  Idents(pbmc_doub)<-'type'
  
  #Find markers
  #For this part we are going to subset the seurat object into individual cell type and find the marker on them
  markers_all = list()
  
  #Check whether you have more than 3 cells for each subtype for the DEG to work
  if(any(as.vector(table(preds_true$preds_true))>=3)){
    for(i in unique(preds)){
      print(i)
      use <- colnames(pbmc_doub)[preds == i]
      pbmc_use <- subset(pbmc_doub, cells = use) #subset based on each cell type
      table_use <- c(length(which(pbmc_use@meta.data$type=='sim')),length(which(pbmc_use@meta.data$type=='true'))) 
      #Find the simulated and true doublet data number
      #If we want to work with this part we need to have more than 3 cells to say that this is a significant result we have found
      #Therefore we would omit when we get cell number which is lower than 3
      if(table_use[1]>=3 & table_use[2]>=3){
        lvels = levels(x = pbmc_use)
        markers <- list(FindMarkers(pbmc_use, ident.1 = lvels[2], ident.2 = lvels[1],min.pct = 0.1,logfc.threshold = 0.1,test.use = "poisson"))
        names(markers) = i
        markers_all = append(markers_all, markers)
        
      }else{
        next
      }
    }
    print(markers_all)
    marker_csv_file <- markers_all
    #Organise data
    for (i in seq_along(markers_all)){
      markers_all[[i]]=markers_all[[i]][markers_all[[i]]$p_val_adj<=0.05,] #Remove non-significant p-vals
    }
    markers_all = Filter(function(x) dim(x)[1] >= 1, markers_all)
    for (i in seq_along(markers_all)){
      if(nrow(markers_all[[i]]) >20){
        markers_all[[i]]$sig=c(rep('1',20),rep('0',nrow(markers_all[[i]])-20)) #Mark top 20 most significant 
        markers_all[[i]]$title=rep(names(markers_all)[[i]]) #Add titles for plots
        markers_all[[i]]$gene=rownames(markers_all[[i]]) #Add gene column
        markers_all[[i]]$log10_p_val = -log10(markers_all[[i]]$p_val_adj)
        print(paste("1",i))
      }else{
        
        markers_all[[i]]$sig=c(rep('1',nrow(markers_all[[i]]))) #Mark top 10 most significant 
        markers_all[[i]]$title=rep(names(markers_all)[[i]]) #Add titles for plots
        markers_all[[i]]$gene=rownames(markers_all[[i]]) #Add gene column
        markers_all[[i]]$log10_p_val = -log10(markers_all[[i]]$p_val_adj)
        print(paste("2",i))
      }
      
    }
    library(ggthemes)
    library(ggrepel)
    #Function for volcano plots
    plot_func <- function(df){
      g = ggplot(df, aes(avg_logFC,log10_p_val)) +
        geom_point(aes(col=sig))+
        scale_color_manual(values=c("black", "red"))+
        theme_base()+
        theme(legend.position = "none")+
        geom_vline(xintercept = 0, linetype="dashed",color = "grey")+
        geom_text_repel(data=df[df$sig==1,], 
                        aes(avg_logFC, log10_p_val,label=rownames(df[df$sig==1,])),
                        min.segment.length = Inf, seed = 42, box.padding = 0.5,cex=5)+
        ggtitle(df$title[1])
      return(g)
    }
    
    plot_list=lapply(markers_all,plot_func)
    
    pdf('Differential_Expression_true_doublet_types.pdf',width = 20, height = 20)
    print(ggarrange(plotlist = plot_list,ncol = 2, nrow = 2))
    dev.off()
    
  }else{
    print("Non enough cells to do DEG analysis")
  }
  
  
  print("-------------New entry Successful-------------")
  #Plot the Different Doublet types identified and the different proportions identified
  doublet_type <- preds_true %>%
    group_by(preds_true) %>%
    summarise(cnt = n()) %>%
    mutate(freq = round(cnt /sum(cnt),3)) %>%
    arrange(desc(freq))
  
  doublet_type$preds_true <- gsub('..'," - ",doublet_type$preds_true, fixed = T)
  doublet_type$preds_true <- gsub('.'," ",doublet_type$preds_true, fixed = T)
  
  write.csv(doublet_type, "Doublet Number and Proportions.csv")
  
  pdf("Number of Different Doublet cell type.pdf", height = 5, width =5)
  p<- ggplot(doublet_type, aes(x = preds_true, y =cnt, fill = preds_true))+
    geom_bar(stat = "identity") + 
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
          axis.text.y = element_text(vjust = 0.5,
                                     hjust = 0.5))+
    theme(legend.position="none")+
    labs(x = "Cell Type", y = "Number of Cells")+
    scale_fill_manual(values = cols)
  print(p)
  dev.off()
  
  pdf("Proportion of Different Doublet cell type.pdf", height = 5, width =5)
  p<- ggplot(doublet_type, aes(x = preds_true, y =freq, fill = preds_true))+
    geom_bar(stat = "identity") + 
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
          axis.text.y = element_text(vjust = 0.5,
                                     hjust = 0.5))+
    theme(legend.position="none")+
    labs(x = "Cell Type", y = "Frequency of Cells")+
    scale_fill_manual(values = cols)
  print(p)
  dev.off()
  
  
  ## Get the hc list and lc list and name them
  
  
  #Change the name to match for each subtype you want
  library(reshape2)
  
  HC<- melt(hc_list)
  HC$L1 <- paste("hc", HC$L1, sep="_")
  LC <- melt(lc_list)
  LC$L1 <- paste("lc", LC$L1, sep="_")
  
  #Adding the droplet information
  
  pbmc@meta.data$Droplet_type <- "NA"
  pbmc@meta.data[as.character(HC$value),]$Droplet_type <- HC$L1
  pbmc@meta.data[as.character(LC$value),]$Droplet_type <- LC$L1
  
  pbmc@meta.data[hc_singlets,]$Droplet_type <- "hc_singlet"
  
  #Adding the Doublet Cell type information to the pbmc dataset
  
  pbmc@meta.data$Doublet_cell_type<- "NA"
  pbmc@meta.data[preds_true$cells,]$Doublet_cell_type <- preds_true$preds_true
  
  #Change the name for proper format
  pbmc@meta.data$Droplet_type <- str_replace_all(pbmc@meta.data$Droplet_type, pattern = "_",replacement = " " )
  print("---------------Part8------------------Successful")
  
  library(stringr)         
  library(ggsci)
  
  dis <- as.data.frame(pbmc@meta.data)
  dis2<- dis[grep("hc", dis$Droplet_type),]

  dis2$Droplet_type<- ifelse(dis2$Droplet_type == "hc heterotypic doublets","HC Doublet B-T cell",
                             ifelse(dis2$Droplet_type =="hc homotypic doublets B", "HC Doublet B-B cell",
                                    ifelse(dis2$Droplet_type =="hc homotypic doublets T", "HC Doublet T-T cell",
                                           ifelse(dis2$Droplet_type =="hc other","HC Other",
                                                  ifelse(dis2$Droplet_type =="hc quadruplets", "HC Quadruplet 2B-2T cell",
                                                         ifelse(dis2$Droplet_type =="hc triplets 2B T","HC Triplet 2B-T cell",
                                                                ifelse(dis2$Droplet_type =="hc triplets B 2T","HC Triplet B-2T cell",
                                                                       ifelse(dis2$Droplet_type == "hc singlet", "HC Singlet","NA"))))))))
  dis2$Droplet_type <-  factor(dis2$Droplet_type, levels = c("HC Singlet","HC Doublet B-T cell","HC Doublet B-B cell","HC Doublet T-T cell","HC Triplet 2B-T cell","HC Triplet B-2T cell","HC Quadruplet 2B-2T cell","HC Other"))
  pdf("Distribution on nUMI per droplet type.pdf", height = 7, width =7)
  p<-ggplot(dis2, aes(y = nUMI, x = Droplet_type, fill = Droplet_type)) +
    geom_boxplot()+
    theme_bw()+
    theme(axis.title.x=element_text(
                                    hjust = 0.5,
                                    size=24),  # X axis title
          axis.title.y=element_text(size=24),  # Y axis title
          axis.text.x=element_text(size=24, 
                                   angle = 90,
                                   vjust=.5,
                                   hjust = 1),  # X axis text
          axis.text.y=element_text(size=24)) + # Y axis text  
    scale_fill_manual(values = cols)+
    theme(legend.position="None")
  print(p)
  dev.off()
  
  print("---------------Part9------------------Successful")
  #Draw a table of Droplets per type
  library(gridExtra)
  library(grid)
  doublet_freq<- as.data.frame(table(pbmc@meta.data$Droplet_type))
  colnames(doublet_freq) <- c("Droplet Type","Frequency")
  write.csv(doublet_freq,"Droplet_freq_table.csv")
  
  pdf("droplet frequency Table.pdf")
  g <- tableGrob(doublet_freq, theme = ttheme_minimal())
  separators <- replicate(ncol(g) - 2,
                          segmentsGrob(x1 = unit(0, "npc"), gp=gpar(lty=2)),
                          simplify=FALSE)
  ## add vertical lines on the left side of columns (after 2nd)
  g <- gtable::gtable_add_grob(g, grobs = separators,
                               t = 2, b = nrow(g), l = seq_len(ncol(g)-2)+2)
  grid.draw(g)
  dev.off()
  
  print("---------------Part10------------------Successful")
  ##Mixed Gene Expression
  
  pbmc_ge <- subset(pbmc, subset = Droplet_type =="hc heterotypic doublets" | Droplet_type =="hc singlet")
  pbmc_ge@meta.data$broad_cell_type <- ifelse(pbmc_ge@meta.data$broad_cell_type == "-","HC B-T Doublets", pbmc_ge@meta.data$broad_cell_type)
  

  
  
  pdf("Expression of B and T cell marker gene.pdf", width  = 7, height =7 )
#  p<-VlnPlot(pbmc_ge,features = c("CD19", "IGKC","MS4A1","CD3E", "CD3G", "CD8A"), group.by = "broad_cell_type", pt.size = 0.01,do.return = TRUE, return.plotlist = TRUE)
  p<-   VlnPlot(pbmc_ge,features = c("CD19", "IGKC","MS4A1","CD3E", "CD3G", "CD8A"), pt.size = 0.001)
  print(p)
  dev.off()
  pdf("NON_MARKER_Singlet vs Doublet.pdf",width  = 7, height = 7)
  p<- VlnPlot(pbmc_ge,features = c("CD14","CCL17","PPBP","FCER1A","CST3","MS4A7"), group.by = "broad_cell_type")
  print(p)
  dev.off()
  print("---------------Part11------------------Successful")
  #Plotting Accuracy by class
  accuracy_by_class_rf <- as.matrix(accuracy_by_class_rf)
  
  
  # Transform the matrix in long format
  df <- reshape::melt(accuracy_by_class_rf)
  class(df)
  colnames(df) <- c("x", "y", "value")
  
  library(dplyr)
  df<- df %>% 
    mutate_if(is.numeric, round, digits = 2)
  
  df$x <- gsub('..',"-",df$x, fixed = T)
  df$x <- gsub('.'," ",df$x, fixed = T)
  
  
  df$y <- gsub('..',"-",df$y, fixed = T)
  df$y <- gsub('.'," ",df$y, fixed = T)
  
  
  pdf("Classification Accuracy.pdf", width =16, height = 16)
  p<-ggplot(df, aes(x = x, y = y, fill = value)) +
    geom_tile(color = "black") +
    ## geom_text(aes(label = value), color = "black") +
    scale_fill_gradient(low = "white", high = "red") +
    coord_fixed()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5,
                                     size = 36,
                                     hjust = 1),
          axis.text.y = element_text(vjust = 0.5,
                                     hjust = 1,
                                     size = 36),
          legend.text = element_text(size=20))
  print(p)
  dev.off()
  print("---------------Part12------------------Successful")
  ## UMAP overlapping between simulated vs HC B-T doublet
  
  # Plot the UMAP of the doublets and the simulation
  
  table(data1@meta.data$new_class)
  data1@meta.data$new_class <- data1@meta.data$class
  data1@meta.data$new_class[rownames(data1@meta.data) %in% final_doublet] <- "HC Heterotypic Doublets"
  data1 <- subset(data1, subset = new_class != "true doublet")
  data1@meta.data$new_class <- ifelse(data1@meta.data$new_class =="doublet","Simulated Doublet",ifelse(data1@meta.data$new_class =="HC Heterotypic Doublets","HC Heterotypic Doublets","Simulated Multiplet"))
  table(data1@meta.data$new_class)
  pdf('doublet_multiplet_umap_cell_types.pdf',width=7,height=6)
  
  p<- DimPlot(data1,reduction = "umap",group.by = 'new_class', pt.size = 0.5)+ #scale_color_manual(values = cols)+
    theme(legend.position="bottom")+
    theme(axis.title.x=element_text(
      hjust = 0.5,
      size=24),  # X axis title
      axis.title.y=element_text(size=24),  # Y axis title
      axis.text.x=element_text(size=24, 
                               angle = 90,
                               vjust=.5,
                               hjust = 1),  # X axis text
      axis.text.y=element_text(size=24),
      legend.text = element_text(size=12))
    
  
  print(p)
  dev.off()
  
  pdf('doublet_Simulated Doublets.pdf',width=8,height=8)
  
  p<-DimPlot(subset(data1,cells=c(doublets_sim, final_doublet)),reduction = "umap",group.by = 'new_class', pt.size = 0.5)+ scale_color_manual(values = cols)+
    theme(legend.position="bottom")+
    theme(axis.title.x=element_text(
      hjust = 0.5,
      size=24),  # X axis title
      axis.title.y=element_text(size=24),  # Y axis title
      axis.text.x=element_text(size=24, 
                               angle = 90,
                               vjust=.5,
                               hjust = 1),  # X axis text
      axis.text.y=element_text(size=24),
      legend.text = element_text(size=12))
  print(p)
  dev.off()
  pdf('grey and red doublet_Simulated Doublets.pdf',width=5,height=5)
  
  p<-DimPlot(subset(data1,cells=c(doublets_sim, final_doublet)),reduction = "umap",group.by = 'new_class', pt.size = 0.5,cols = c("darkred","grey"))+
    theme(legend.position="bottom")+
    theme(axis.title.x=element_text(
      hjust = 0.5,
      size=24),  # X axis title
      axis.title.y=element_text(size=24),  # Y axis title
      axis.text.x=element_text(size=24, 
                               angle = 90,
                               vjust=.5,
                               hjust = 1),  # X axis text
      axis.text.y=element_text(size=24),
      legend.text = element_text(size=12))
  print(p)
  dev.off()
  pdf('Simulated Doublets all labels.pdf',width=8,height=5)
  
  p<-DimPlot(subset(data1,cells=c(doublets_sim)),reduction = "umap", pt.size = 0.3,label=F, repel = TRUE)+
    theme(legend.position="right")+
    theme(axis.title.x=element_text(
      hjust = 0.5,
      size=24),  # X axis title
      axis.title.y=element_text(size=24),  # Y axis title
      axis.text.x=element_text(size=24, 
                               angle = 90,
                               vjust=.5,
                               hjust = 1),  # X axis text
      axis.text.y=element_text(size=24),
      legend.text = element_text(size=16))
  print(p)
  dev.off()
  print("---------------Part13------------------Successful")
  
  
  fileout1=concat(c("Enrichment_of_B_T_cell_doublet_network_",".pdf"))
  pdf(file=fileout1)
  BCR_TCR_sub = unlist(hc_list["heterotypic_doublets"])
  clone_BCR = apply(cbind(m_VDJ_BCR[BCR_TCR_sub,"clone1"], m_VDJ_BCR[BCR_TCR_sub,"clone2"]), 1, paste, collapse = "_")
  clone_TCR = apply(cbind(m_VDJ_TCR[BCR_TCR_sub,"clone1"], m_VDJ_TCR[BCR_TCR_sub,"clone2"]), 1, paste, collapse = "_")
  
  string_sample <- strsplit(unique(orig.ident),"_")[[1]][2]
  clone_BCR = gsub(concat(paste("||","RCC-metP", sep ="")),"B",clone_BCR, fixed = T)
  clone_TCR = gsub(concat(paste("||","RCC-metP", sep ="")),"T",clone_TCR, fixed = T)
  
  t_clones = table(c(clone_BCR, clone_TCR))
  clone_ids = unique(c(clone_BCR, clone_TCR))
  library(igraph)
  g <- graph.empty( n=0, directed=FALSE)
  g <- igraph::add.vertices(g, length(clone_ids), name= clone_ids) 
  names <- V(g)$name
  ids <- 1:length(names)
  names(ids) <- names
  
  x = cbind(clone_BCR, clone_TCR)
  x_unique = unique(x)
  from = NULL
  to = NULL
  edge_strength = NULL
  for(i in c(1:length(x_unique[,1]))){
    from = c(from, x_unique[i,1])
    to = c(to, x_unique[i,2])
    w = intersect(which(x[,1]==x_unique[i,1]), which(x[,2]==x_unique[i,2]))
    edge_strength = c(edge_strength, length(w))
  }
  edges <- matrix(c(ids[from], ids[to]), nc=2)
  g <- add.edges(g, t(edges), weight= edge_strength)
  
  V(g)$label <- V(g)$name
  V(g)$size<-1
  V(g)$label.cex<-0.0001
  
  #layout1 =layout_with_mds(g1, dist = NULL, dim = 2, options = arpack_defaults)
  #layout = layout_with_graphopt(g1,niter = 100,start = layout1, spring.length = 1)
  layout = layout_with_graphopt(g,niter = 600)
  # layout1 =layout_in_circle(g)
  # layout1 = layout_with_graphopt(g,niter = 500)
  rownames(layout) = V(g)$name
  
  library(yarrr)
  library(RColorBrewer)
  library(prettyGraphs)
  
  sizes = rep(0, length(clone_ids))
  names(sizes) = clone_ids
  tx = table(clone_BCR)
  sizes[names(tx)]= tx
  tx = table(clone_TCR)
  sizes[names(tx)]= tx
  
  sizes1 = sizes^0.5
  thresh1 = 50^0.5
  sizes1 = (sizes1*20/thresh1)-1
  V(g)$size<-sizes1
  cols=add.alpha(colorRampPalette(c("pink",'red', "darkred"))(100), alpha = 0.7)
  m1 = edge_strength*100/max(edge_strength)
  m1 = round(m1, digits = 0)
  range(m1)
  E(g)$color = cols[m1]
  edge_strength_plot = edge_strength ^0.5
  edge_strength_plot = (edge_strength_plot*7/3.7)+0.05
  
  col1 =  add.alpha (brewer.pal(8, "Dark2")[c(1,6)], alpha = 0.95)
  cols_node = rep(col1[1], length(V(g)$name))
  cols_node[grep("B",V(g)$name)] = col1[2]
  V(g)$color = cols_node
  # V(g)[grep("BCR",V(g)$name)]$shape <- "square"
  # V(g)[grep("TCR",V(g)$name)]$shape <- "circle"
  V(g)$label.cex <- rep(0.0000001, length(V(g)$name))
  # w = which(sizes>=4)
  # V(g)$label.cex[w] <- 0.8
  plot(g, layout=layout, edge.color= E(g)$color, main="", edge.width= edge_strength_plot,vertex.label.family="sans", edge.lty = 1,xlim = c(-1, 1.5), ylim = c(-1, 1),cex.main =3) # vertex.label.cex = 0.4
  # title(main = samples)
  layout2 = layout
  layout2[,1] = layout2[,1]/max(abs(layout2[,1]))
  layout2[,2] = (layout2[,2]/max(abs(layout2[,2])))
  mult = 1.025
  layout2 = layout2*mult
  w = which(sizes>=4)
  if(length(w)>1){
    plot_ids = names(w)
    xy = layout2[plot_ids,]
    x_lab_pos= 1.3
    order = order(xy[,2])
    y_lab_pos = seq(from = -0.9, to = 0.7, length = length(plot_ids))
    plot_ids1 =gsub("_BCR,","-", plot_ids)
    plot_ids1 =gsub("_TCR,","-", plot_ids1)
    text(x_lab_pos, y = y_lab_pos, label = plot_ids1[order], cex = 1.2, col = "darkblue")
    # text(layout2[plot_ids,], label = plot_ids)		
    for(i in c(1:length(order))){
      segments(xy[order[i],1], xy[order[i],2], x_lab_pos-0.2, y_lab_pos[i], col = add.alpha("darkblue",alpha = 0.5), lty = 3, lwd = 1)
    }
  }
  
  legend("topleft",  legend= c("T cell","B cell"), pch= 21, bty="n",cex= 1.2, title="", pt.bg = col1,pt.cex = 1.1)
  scale = c(1,max(sizes))
  scale1 = scale ^0.5
  thresh1 = 90^0.5
  scale1 = (scale1*20/thresh1)-1
  legend("bottomleft",  legend= as.character(scale), pch= 21, bty="n",cex= 1, title="", pt.bg = "grey",pt.cex = scale1*0.4, col = "grey",y.intersp=2, xjust=0.5, x.intersp = 2)
  
  dev.off()
  print("---------------Part14------------------Successful")
  
  #chi square test for the expected and observed counts
  orig_cells<- pbmc@meta.data %>% 
    group_by(cell_refined_annotation, broad_cell_type) %>%
    summarize(num_cells =n())
  
  orig_cells <- as.data.frame(orig_cells[orig_cells$cell_refined_annotation != "-",])
  
  # Define the cells and their frequencies
  b_cells <- filter(orig_cells, broad_cell_type == "B cell") 
  
  t_cells <-filter(orig_cells, broad_cell_type == "T cell") 
  b_cells_name <- make.names(b_cells$cell_refined_annotation)
  t_cells_name <- make.names(t_cells$cell_refined_annotation)
  freqs_b <- b_cells$num_cells
  freqs_t <- t_cells$num_cells
  
  # Define the observed frequencies
  observed_cells <- as.data.frame(table(preds_true$preds_true))
  observed <- observed_cells$Freq
  names(observed) <- observed_cells$Var1
  
  # Calculate the total number of combinations
  n_combinations <- sum(observed)
  # Calculate the expected frequencies for each combination
  probs_b <- freqs_b / sum(freqs_b)
  
  probs_t <- freqs_t / sum(freqs_t)
  
  expected <- numeric(length(probs_t*probs_t))
  expected
  idx <- 1
  name <- c()
  for (i in 1:length(t_cells_name)) {
    for (j in 1:length(b_cells_name)) {
      print(paste(t_cells_name[i],b_cells_name[j]))
      expected[idx] <- probs_t[i] * probs_b[j] * n_combinations
      cell_name <- paste(t_cells_name[i],b_cells_name[j], sep = "..")
      name <- append(name,cell_name)
      print(expected[idx])
      idx <- idx + 1
    }
  }
  expected <- expected/sum(expected)
  sum(expected)
  names(expected) <- name
  a<- rep(0, length(expected))
  names(a) <- name
  a[names(observed)] <- observed
  a
  observed <- a
  observed
  # Create a table of observed and expected frequencies
  table <- data.frame(Observed = observed, Expected = expected)
  table
  rownames(table) <- names(observed)
  
  # Perform the chi-square test
  chisq <- chisq.test(observed, p = expected)
  
  # Print the table
  print(table)
  
  table$Expected <- table$Expected * n_combinations
  
  
  library(ggplot2)
  library(dplyr)
  
  library(reshape2)
  library(ggsci)
  table$cell_type <- rownames(table)
  table<-gather(table, key = "NAME", value = "Value",  -cell_type)
  table
  pdf("chisquare test.pdf", width = 5, height =7)
  p<- ggplot(table, aes(x = cell_type, y = Value, fill = NAME)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(x = "Combination of two cells", y = "Frequency",
         caption = paste("p value = ",chisq$p.value)) +
    scale_fill_npg()+
    theme_classic()+
    theme(axis.text.x = element_text(vjust = 0.5, angle = 90, size =14),
          axis.text.y = element_text(hjust = 0.5, size = 14))
  print(p)
  
  dev.off()
  
  
  #new entry-----------------
  #Mito ribo ratio and mito percentage
  mitrib <- grep(pattern = "^MT-|^RP[SL]", x = rownames(x = pbmc@assays$RNA@counts), value = TRUE, ignore.case = TRUE)
  mito <- grep(pattern = "^MT-", x = rownames(x = pbmc@assays$RNA@counts), value = TRUE, ignore.case = TRUE)
  mitribcounts<- pbmc@assays$RNA@counts[which(rownames(pbmc@assays$RNA@counts) %in% mitrib), ]
  mitoribo_ratio <- Matrix::colSums(mitribcounts[mito, , drop = FALSE])/Matrix::colSums(mitribcounts)
  pbmc <- AddMetaData(object = pbmc, metadata = as.data.frame(mitoribo_ratio) ,col.name = "mito.ribo_ratio")
  
  
  pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
  
  
  
  # Plot both the data
  #mito percentage
  
  x<- pbmc@meta.data
  
  x1<- x[x$Droplet_type == "hc heterotypic doublets" |x$Droplet_type == "hc singlet", ]
  
  pdf("mitopercentage.pdf")
  
  p<-ggplot(x1, aes(x = Droplet_type, y = percent.mt, fill = Droplet_type)) +
    geom_boxplot()+
    theme_bw()
  
  print(p)
  dev.off()
  
  pdf("mito_ribo_ratio.pdf")
  
  p<- ggplot(x1, aes(x = Droplet_type, y = mito.ribo_ratio, fill = Droplet_type)) +
    geom_boxplot()+
    theme_bw()
  
  print(p)
  dev.off()
  
  
  write.csv(pbmc@meta.data, "metadata.csv")
  saveRDS(pbmc, "Tissue2_with_mito_percentage.RDS")
  
  print("final success")
  return(pbmc)
  
  
  
}

potential_interaction = double_type_prediction(pbmc, hc_doublets, hc_singlets)

write.csv(potential_interaction, "tissue6653.csv")
write.csv(potential_interaction, "tissue6100.csv")
write.csv(potential_interaction, "blood.csv")




#Calculating mito ribo ratio and mito percentage
#calculating the mito_ribo ratio of the model. I havent used it in the model but there is an interesting 
#chance that I would start to do such soon in the model

mitrib <- grep(pattern = "^MT-|^RP[SL]", x = rownames(x = pbmc@assays$RNA@counts), value = TRUE, ignore.case = TRUE)
mito <- grep(pattern = "^MT-", x = rownames(x = pbmc@assays$RNA@counts), value = TRUE, ignore.case = TRUE)
mitribcounts<- pbmc@assays$RNA@counts[which(rownames(pbmc@assays$RNA@counts) %in% mitrib), ]
mitoribo_ratio <- Matrix::colSums(mitribcounts[mito, , drop = FALSE])/Matrix::colSums(mitribcounts)
pbmc <- AddMetaData(object = pbmc, metadata = as.data.frame(mitoribo_ratio) ,col.name = "mito.ribo_ratio")


pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")



# Plot both the data
#mito percentage

x<- pbmc@meta.data

x1<- x[x$Droplet_type == "hc heterotypic doublets" |x$Droplet_type == "hc singlet", ]

pdf("mitopercentage.pdf")

p<-ggplot(x1, aes(x = Droplet_type, y = percent.mt, fill = Droplet_type)) +
  geom_boxplot()+
  theme_bw()

print(p)
dev.off()

pdf("mito_ribo_ratio.pdf")

p<- ggplot(x1, aes(x = Droplet_type, y = mito.ribo_ratio, fill = Droplet_type)) +
  geom_boxplot()+
  theme_bw()

print(p)
dev.off()

saveRDS(pbmc, "RCC_with_mito_percentage_and_mito_ribo_ratio.rds")

hc_list[["heterotypic_doublets"]]=setdiff(hc_list[["heterotypic_doublets"]], potential_multiplets)
lc_list[["other"]] = potential_multiplets

## (3) get nUMI outliers from doublets: doublets outlier function
hc_doublets = hc_list[["heterotypic_doublets"]]
doublet_outliers = Doublet_outlier(pbmc, hc_doublets)
hc_list[["heterotypic_doublets"]]=setdiff(hc_list[["heterotypic_doublets"]], doublet_outliers)
lc_list[["other"]] = doublet_outliers

for(c in c(1:length(names_cols))){
  mat_hc[c,i] = length(hc_list[[names_cols[c]]])
  mat_lc[c,i] = length(lc_list[[names_cols[c]]])
}
n_hc_singlets = c(n_hc_singlets, length(hc_singlets))

droplet_type = rep("-", length(VDJ_ids))
names(droplet_type) = VDJ_ids
for(c in c(1:length(names_cols))){
  droplet_type[hc_list[[names_cols[c]]]] = concat(c("HC_", names_cols[c]))
  droplet_type[lc_list[[names_cols[c]]]] = concat(c("LC_", names_cols[c]))
}
droplet_type[hc_singlets] = "HC singlets"
#check everything is named (i.e. no "-"s)
table(droplet_type)
classification_list = c(classification_list, list(droplet_type))
}


print(n_hc_doublets)

fileout1=concat(c("/well/immune-rep/shared/10X_GENOMICS/RCC_WORKING_DATA/FINAL/Tmp_TB_doublet_count_nUMI_threshold.pdf"))
w=2.5
pdf(file=fileout1, height=w*3, width=w*2)
par(mfrow= c(3,2), mar = c(4,4,4,4))
library(RColorBrewer)
cols <- brewer.pal(8, "Dark2")

## high confidence
range = range(mat_hc)
plot(range(threshold_umis_all), range, pch = 21, col = "white",bg = "white", ylim =c(0,max(range)), ylab = "n droplets", xlab = "threshold")
for(i in c(1:length(names_cols))){
  points(threshold_umis_all,mat_hc[i,], type = "l", col = cols[i], lwd = 1.5)
  points(threshold_umis_all,mat_hc[i,], pch = 21, col = cols[i],bg = cols[i], cex =0.5)
}
plot(range, pch = 21, col = "white",bg = "white", ylab = "", xlab = "", axes = F)
labs = apply(cbind("HC ", names_cols), 1, paste, collapse = "")
legend("topleft",  legend= labs, pch= 21, bty="n",cex= 0.9, title="", pt.bg = cols,pt.cex = 0.9, col = cols)

## low confidence
range = range(mat_lc)
plot(range(threshold_umis_all), range, pch = 21, col = "white",bg = "white", ylim =c(0,max(range)), ylab = "n droplets", xlab = "threshold")
for(i in c(1:length(names_cols))){
  points(threshold_umis_all,mat_lc[i,], type = "l", col = cols[i], lwd = 1.5)
  points(threshold_umis_all,mat_lc[i,], pch = 21, col = cols[i],bg = cols[i], cex =0.5)
}
plot(range, pch = 21, col = "white",bg = "white", ylab = "", xlab = "", axes = F)
labs = apply(cbind("LC ", names_cols), 1, paste, collapse = "")
legend("topleft",  legend= labs, pch= 21, bty="n",cex= 0.9, title="", pt.bg = cols,pt.cex = 0.9, col = cols)

## singlets
range = range(n_hc_singlets)
plot(range(threshold_umis_all), range, pch = 21, col = "white",bg = "white", ylim =c((range)), ylab = "n droplets", xlab = "threshold")
i = 1
points(threshold_umis_all,n_hc_singlets, type = "l", col = cols[i], lwd = 1.5)
points(threshold_umis_all,n_hc_singlets, pch = 21, col = cols[i],bg = cols[i], cex =0.5)

plot(range, pch = 21, col = "white",bg = "white", ylab = "", xlab = "", axes = F)
labs = "n_hc_singlets"
legend("topleft",  legend= labs, pch= 21, bty="n",cex= 0.9, title="", pt.bg = cols,pt.cex = 0.9, col = cols)

dev.off()



#######################---------------------- <- <- <- <- ############################






tissue6653 <- read.csv("/Users/mohammadumersharifshohan/Desktop/Everything Oxford/Lab Work/Starter Sharif/B T cell Doublet/New Algorithm/New Algorithm from Beginning/Plan for Writing /Figure/Figure for all/Tissue 6653 sample run/doublet_types_with_cell_name.csv")
tissue6100 <- read.csv("/Users/mohammadumersharifshohan/Desktop/Everything Oxford/Lab Work/Starter Sharif/B T cell Doublet/New Algorithm/New Algorithm from Beginning/Plan for Writing /Figure/Figure for all/tissue 6100/doublet_types_with_cell_name.csv")
blood <- read.csv("/Users/mohammadumersharifshohan/Desktop/Everything Oxford/Lab Work/Starter Sharif/B T cell Doublet/New Algorithm/New Algorithm from Beginning/Plan for Writing /Figure/Figure for all/Blood/doublet_types_with_cell_name.csv")
table(tissue6653$preds_true)
tissue6653 <- read.csv("/Users/mohammads/Downloads/Everything Oxford/Lab Work/B T cell Doublet/New Algorithm/New Algorithm from Beginning/PDAC_everything/PDAC3/Tissue_1/doublet_types_with_cell_name.csv")
tissue6100 <- read.csv("/Users/mohammads/Downloads/Everything Oxford/Lab Work/B T cell Doublet/New Algorithm/New Algorithm from Beginning/PDAC_everything/PDAC3/Tissue_2/doublet_types_with_cell_name.csv")
blood <- read.csv("/Users/mohammads/Downloads/Everything Oxford/Lab Work/B T cell Doublet/New Algorithm/New Algorithm from Beginning/PDAC_everything/PDAC3/Tissue_combined/doublet_types_with_cell_name.csv")
x<- c("a","b","c","a","b","c","a","b","c","a","b","a","b","a","b","c","a","c","a","c","a")


library(dplyr)

x<- tissue6653 %>%
  group_by(preds_true) %>%
  summarise(cnt = n()) %>%
  mutate(freq = round(cnt /sum(cnt),3)) %>%
  arrange(desc(freq))
y<- tissue6100 %>%
  group_by(preds_true) %>%
  summarise(cnt = n()) %>%
  mutate(freq = round(cnt /sum(cnt),3)) %>%
  arrange(desc(freq))
z<- blood %>%
  group_by(preds_true) %>%
  summarise(cnt = n()) %>%
  mutate(freq = round(cnt /sum(cnt),3)) %>%
  arrange(desc(freq))
z




merged_tissue<- merge(x, y, by="preds_true", all.x = TRUE)
merged_tissue[is.na(merged_tissue)] <- 0

merged_tissue
install.packages("ggpubr")
library(ggpubr)
colnames(merged_tissue) <- c("Cell_Type","cntx","Tissue_1", "cnty", "Tissue_2")
pdf("Tissue1_vs_Tissue2.pdf")
ggplot(merged_tissue, aes(Tissue_1, Tissue_2 ))+
  stat_cor() +
  geom_point(size=3, aes(col = Cell_Type))+
  geom_smooth(method='lm', se=FALSE, color='black') +
  theme_bw() +
  labs(title='Linear Regression Plot') +
  theme(plot.title = element_text(hjust=0.5, size=20, face='bold'))
dev.off()     

colnames(z)[1] <- "Cell_Type"
merged_b_t <- merge(merged_tissue, z, by = "Cell_Type", all.x = TRUE)   
merged_b_t

merged_b_t[is.na(merged_b_t)] <- 0

colnames(merged_b_t) <- c("Cell_Type","cntx","Tissue_1", "cnty", "Tissue_2", "cntz","Combined")
library(tidyr)
a<- merged_b_t %>%
  gather(key = "Sample_type", value = "Proportion", -c(cntx,cnty,cntz,Cell_Type))
a
str(a)
a$Cell_Type<- gsub(".."," - ", a$Cell_Type, fixed =T)
a$Cell_Type<- gsub("."," ", a$Cell_Type, fixed = T)
library(ggsci)

pdf("proportions_of_Doublets.pdf")
ggplot(a, aes(x = Cell_Type, y  = Proportion, fill = Sample_type)) +
  geom_bar(colour="black", stat="identity", position = position_dodge(width = 0.8), width=0.7) +
  theme_bw()+
  scale_fill_npg(name = "Sample")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 15),
        axis.text.y = element_text(size = 15,hjust = 0.5))+
  labs(x = "Cell Type", y = "Proportions")

dev.off()      


#Droplet proportion
tissue6653 <- read.csv("/Users/mohammadumersharifshohan/Desktop/Everything Oxford/Lab Work/Starter Sharif/B T cell Doublet/New Algorithm/New Algorithm from Beginning/Plan for Writing /Figure/Figure for all/Tissue 6653 sample run/Droplet_freq_table.csv")
tissue6100 <- read.csv("/Users/mohammadumersharifshohan/Desktop/Everything Oxford/Lab Work/Starter Sharif/B T cell Doublet/New Algorithm/New Algorithm from Beginning/Plan for Writing /Figure/Figure for all/tissue 6100/Droplet_freq_table.csv")
blood <- read.csv("/Users/mohammadumersharifshohan/Desktop/Everything Oxford/Lab Work/Starter Sharif/B T cell Doublet/New Algorithm/New Algorithm from Beginning/Plan for Writing /Figure/Figure for all/Blood/Droplet_freq_table.csv")

tissue6653 <- read.csv("/Users/mohammads/Downloads/Everything Oxford/Lab Work/B T cell Doublet/New Algorithm/New Algorithm from Beginning/PDAC_everything/PDAC3/Tissue_1/Droplet_freq_table.csv")
tissue6100 <- read.csv("/Users/mohammads/Downloads/Everything Oxford/Lab Work/B T cell Doublet/New Algorithm/New Algorithm from Beginning/PDAC_everything/PDAC3/Tissue_2/Droplet_freq_table.csv")
blood <- read.csv("/Users/mohammads/Downloads/Everything Oxford/Lab Work/B T cell Doublet/New Algorithm/New Algorithm from Beginning/PDAC_everything/PDAC3/Tissue_combined/Droplet_freq_table.csv")
x<- c("a","b","c","a","b","c","a","b","c","a","b","a","b","a","b","c","a","c","a","c","a")





library(dplyr)
tissue6100
x<- tissue6653 %>%
  mutate(freq = round(Frequency /sum(Frequency),3)) %>%
  arrange(desc(freq))
y<- tissue6100  %>%
  mutate(freq = round(Frequency /sum(Frequency),3)) %>%
  arrange(desc(freq))
z<- blood %>%
  mutate(freq = round(Frequency /sum(Frequency),3)) %>%
  arrange(desc(freq))
x$Sample <- "Tissue_1"
y$Sample <- "Tissue_2"
z$Sample <- "Combined"

all_sample <- rbind(x,y)
all_sample <- rbind(all_sample,z)
head(all_sample)
all_sample <- all_sample[all_sample$Droplet.Type !="hc singlet",]
library(ggplot2)
library(ggsci)
scale_color_n
all_sample$droplet<-c("HC singlets","LC other","HC 2T doublets","HC B-T doublets",
                      "HC other","HC B 2T triplets","LC 2T doublets","HC 2B 2T quadruplets",
                      "HC 2B T triplets","HC singlets","LC other","HC 2T doublets",
                      "HC B-T doublets","LC 2T doublets","HC other","HC B 2T triplets",
                      "HC 2B doublets","LC 2B doublets","LC 2B T triplets","LC B 2T triplets",
                      "HC singlets","LC other","HC B-T doublets","HC 2T doublets",
                      "HC other","HC B 2T triplets","HC 2B doublets","LC 2T doublets",
                      "HC 2B 2T quadruplets","HC 2B T triplets") 
pdf("Proportions of Droplet Type.pdf", width = 4, height = 4)
ggplot(all_sample, aes(x = droplet, y = freq, shape = Sample, col = Sample)) +
  geom_point(position = position_dodge(width = 0.6), size = 1)+
  scale_color_npg()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        axis.text.y = element_text(hjust = 0.5))+
  labs(x = "Droplet Type", y = "Proportion")
dev.off()

pdf("Proportions of Droplet Type.pdf")
ggplot(all_sample, aes(x = Droplet.Type, y = freq, shape = Sample, col = Sample)) +
  geom_point(position = position_dodge(width = 0.6), size = 1.5)+
  scale_color_npg()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        axis.text.y = element_text(hjust = 0.5))+
  labs(x = "Droplet Type", y = "Proportion")
dev.off()





#### Distrub=ution across tissues or blood vs tissue

#Droplet proportion
tissue6653 <- read.csv("/Users/mohammadumersharifshohan/Desktop/Everything Oxford/Lab Work/Starter Sharif/B T cell Doublet/New Algorithm/New Algorithm from Beginning/Plan for Writing /Figure/Figure for all/Tissue 6653 sample run/Droplet_freq_table.csv")
tissue6100 <- read.csv("/Users/mohammadumersharifshohan/Desktop/Everything Oxford/Lab Work/Starter Sharif/B T cell Doublet/New Algorithm/New Algorithm from Beginning/Plan for Writing /Figure/Figure for all/tissue 6100/Droplet_freq_table.csv")
blood <- read.csv("/Users/mohammadumersharifshohan/Desktop/Everything Oxford/Lab Work/Starter Sharif/B T cell Doublet/New Algorithm/New Algorithm from Beginning/Plan for Writing /Figure/Figure for all/Blood/Droplet_freq_table.csv")

library(dplyr)
tissue6100
x<- tissue6653 %>%
  mutate(freq = round(Frequency /sum(Frequency),4)) %>%
  arrange(desc(freq))
y<- tissue6100  %>%
  mutate(freq = round(Frequency /sum(Frequency),4)) %>%
  arrange(desc(freq))
z<- blood %>%
  mutate(freq = round(Frequency /sum(Frequency),4)) %>%
  arrange(desc(freq))




merged_tissue<- merge(x, y, by="Droplet.Type", all.x = TRUE)
merged_tissue[is.na(merged_tissue)] <- 0


merged_tissue
library(ggpubr)
colnames(merged_tissue) <- c("Droplet Type","x","cntx","Tissue_1","y", "cnty", "Tissue_2")
merged_tissue$droplet<-c("HC B-T doublets","HC 2T doublets","HC other","HC 2B 2T quadruplets",
                         "HC singlets","HC 2B T triplets","HC B 2T triplets","LC 2T doublets",
                         "LC other")
merged_tissue2<- merged_tissue
merged_tissue2 <- merged_tissue2[merged_tissue2$droplet !="HC singlets",]
pdf("Tissue1_vs_Tissue2.pdf" ,height =5, width  = 5)
ggplot(merged_tissue2, aes(Tissue_1, Tissue_2 ))+
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),r.accuracy = 0.0001) +
  geom_point(size=2, aes(col = droplet))+
  geom_smooth(method='lm', se=FALSE, color='black') +
  theme_bw() +
  labs(title='Linear Regression Plot') +
  theme(plot.title = element_text(hjust=0.5, face='bold'))
dev.off()  

merged_tissue



x<- tissue6653 %>%
  mutate(freq = round(Frequency /sum(Frequency),4)) %>%
  arrange(desc(freq))
y<- tissue6100  %>%
  mutate(freq = round(Frequency /sum(Frequency),4)) %>%
  arrange(desc(freq))
z<- blood %>%
  mutate(freq = round(Frequency /sum(Frequency),4)) %>%
  arrange(desc(freq))




merged_tissue<- merge(y, z, by="Droplet.Type", all.x = TRUE)
merged_tissue[is.na(merged_tissue)] <- 0
colnames(merged_tissue) <- c("Droplet Type","x","cntx","Tissue_2","z", "cntz", "Blood")
merged_tissue
merged_tissue$droplet<-c("HC B-T doublets","HC 2B doublets","HC 2T doublets","HC other",
                         "HC singlets","HC 2B T triplets","LC 2B doublets","LC 2T doublets",
                         "LC other","LC 2B T triplets","LC B 2T triplets")

pdf("Tissue2_vs_Blood.pdf", height = 5, width =5)
merged_tissue2<- merged_tissue
merged_tissue2 <- merged_tissue2[merged_tissue2$droplet !="HC singlets",]
ggplot(merged_tissue2, aes(Tissue_2, Blood ))+
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),r.accuracy = 0.0001) +
  geom_point(size=2, aes(col = droplet))+
  geom_smooth(method='lm', se=FALSE, color='black') +
  theme_bw() +
  labs(title='Linear Regression Plot') +
  theme(plot.title = element_text(hjust=0.5, face='bold'))
dev.off()  
merged_tissue



ggplot(merged_tissue, aes(x = Tissue_2, y  = Blood, col =droplet))+
  geom_point()
