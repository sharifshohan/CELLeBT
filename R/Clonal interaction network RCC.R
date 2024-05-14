
Draw_box_plot<-function(box,x,width,c,lwd,line_col){
  segments(x, box[2], x, box[3], col = line_col,lwd =lwd)
  segments(x-(width/2), box[2], x+(width/2), box[2], col = line_col,lwd =lwd)
  segments(x-(width/2), box[3], x+(width/2), box[3], col = line_col,lwd =lwd)
  rect(x-width, box[4], x+width, box[5], col = c,lwd =lwd, border = line_col)
  segments(x-width, box[1], x+width, box[1], col = line_col,lwd=2*lwd)}

Means_factor = function(factor, x){
  m = NULL
  for(i1 in c(1:length(levels(factor)))){
    x1 = x[which(factor==levels(factor)[i1])]
    x1 = x1[which(x1!=-1)]
    m = c(m, mean(x1))}
  return(m)}

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


type = "ALL"
analysis = "DOUBLETS"
batch = "RCCmPan"
outdir = "~/Google_Drive/Projects/Single cell analysis/InfiltrImmune/PROJECT_CITE-seq CC/Doublet characterisation project/Doublet analysis/"

############ read doublet classification files
###########
files = c(concat(c(outdir,"Blood_with_mito_percentage.RDS")), 
          concat(c(outdir,"Tissue1_with_mito_percentage.RDS")),
          concat(c(outdir,"Tissue2_with_mito_percentage.RDS")))
names(files) = c("RCC blood", "RCC tissue 1","RCC tissue 2")

library(Seurat)

information = NULL
plotting_info = NULL
for(f in c(1:length(files))){
  data = readRDS(file = files[f])
  droplet_type = data@meta.data[,"Droplet_type"]
  data1 = data@meta.data[which(droplet_type=="hc heterotypic doublets"),]
  information = c(information, list(data1))
  print(f)
}
names(information) = names(files)

all_ids = NULL
for(f in c(1:length(files))){
  all_ids = c(all_ids, rownames(information[[f]]))
}

## combine tissue 1 and 2
information1 = c(list(information[[1]]), list(rbind(information[[2]], information[[3]])))
names(information1) = c("blood", "tumour")

######## create network

############ match VDJ
files_VDJ = c("~/Google_Drive/Projects/Single cell analysis/InfiltrImmune/PROJECT_RCC_metP/Plots/BCR_TCR/BCR_information_all_RCCmetP_RCC.txt", 
              "~/Google_Drive/Projects/Single cell analysis/InfiltrImmune/PROJECT_RCC_metP/Plots/BCR_TCR/TCR_information_all_RCCmetP_RCC.txt")
names(files_VDJ) = c("BCR","TCR")
mVDJ = NULL
mVDJ_all = NULL
for(f in c(1:length(files_VDJ))){
  p <- as.matrix(read.csv(files_VDJ[f], head=T, sep="\t"))
  id = strsplit(rownames(p), "||", fixed = T)
  for(i in c(1:length(id))){
    id[i] = concat(c(id[[i]][2], "_", id[[i]][1]))
  }
  id = unlist(id)
  length(intersect(id, all_ids))
  length(all_ids)
  rownames(p) = id
  mVDJ_all = c(mVDJ_all, list(p))
  p1 = p[all_ids,]
  mVDJ = c(mVDJ, list(p1))
  print (f)
}
names(mVDJ) = names(files_VDJ)
names(mVDJ_all) = names(files_VDJ)

### clone size versus doublet singlet

m_VDJ_BCR = mVDJ$BCR
m_VDJ_TCR = mVDJ$TCR

## here just shortening the names
sample_name = c("04_Kidney_met_panc1_blood","Kidney_met_panc1_biopsy_CD45")
clone_BCR = apply(cbind(mVDJ_all$BCR[,"clone1"], mVDJ_all$BCR[,"clone2"]), 1, paste, collapse = "_")
clone_TCR = apply(cbind(mVDJ_all$TCR[,"clone1"], mVDJ_all$TCR[,"clone2"]), 1, paste, collapse = "_")
clone_BCR = gsub(concat(c("||RCC-metP")),"B",clone_BCR, fixed = T)
clone_TCR = gsub(concat(c("||RCC-metP")),"T",clone_TCR, fixed = T)

fileout1=concat(c(outdir, "Clone_sizes_singlet_doublet_", batch,".pdf"))
w=2.7
pdf(file=fileout1, height=w*1.1*2, width=w*1.1*2)
par(mfrow= c(2,2), mar = c(5,5,5,5))

for(s in c(1:length(information1))){
  sample = names(information1)[s]
  ids_sample = names(clone_BCR)[grep(sample_name[s], names(clone_BCR))]
  tb = table(clone_BCR[ids_sample])
  tt = table(clone_TCR[ids_sample])
  tb = tb[grep("-", names(tb), invert = T)]
  tt = tt[grep("-", names(tt), invert = T)]
  BCR_TCR_sub = rownames(information1[[s]])
  doublet_bclone= unique(clone_BCR[BCR_TCR_sub])
  doublet_tclone= unique(clone_TCR[BCR_TCR_sub])
  doublet_bclone = doublet_bclone[grep("-", doublet_bclone, invert = T)]
  doublet_tclone = doublet_tclone[grep("-", doublet_tclone, invert = T)]
  singlet_bclone = setdiff(names(tb), doublet_bclone)
  singlet_tclone = setdiff(names(tt), doublet_tclone)
  singlet_bclone = singlet_bclone[grep("-", singlet_bclone, invert = T)]
  singlet_tclone = singlet_tclone[grep("-", singlet_tclone, invert = T)]
  
  t_B = cbind(tb*0, tb*0)
  t_B1 = cbind(tb*0, tb*0)
  colnames(t_B) = c("singlet","doublet")
  colnames(t_B1) = c("singlet","doublet")
  doublet_bclone= table(clone_BCR[BCR_TCR_sub])
  doublet_bclone = doublet_bclone[grep("-", names(doublet_bclone), invert = T)]
  t_B[names(doublet_bclone),"doublet"] = doublet_bclone*100/sum(doublet_bclone)
  t_B1[names(doublet_bclone),"doublet"] = doublet_bclone
  singlet_bclone= table(clone_BCR[setdiff(ids_sample,BCR_TCR_sub)])
  singlet_bclone = singlet_bclone[grep("-", names(singlet_bclone), invert = T)]
  t_B[names(singlet_bclone),"singlet"] = singlet_bclone*100/sum(singlet_bclone)
  t_B1[names(singlet_bclone),"singlet"] = singlet_bclone
  
  t_Bu = unique(t_B1)
  count = NULL
  p_value = NULL
  for(i in c(1:length(t_Bu[,1]))){
    w = intersect(which(t_B1[,1]==t_Bu[i,1]), which(t_B1[,2]==t_Bu[i,2]))
    count = c(count, length(w))
    x =cbind(t_Bu[i,], colSums(t_Bu[-i,]))
    p_value = c(p_value, fisher.test(x, alternative = "two")$p.value)
  }
  sort(p_value)
  
  t_Bu = unique(t_B)
  count = NULL
  for(i in c(1:length(t_Bu[,1]))){
    w = intersect(which(t_B[,1]==t_Bu[i,1]), which(t_B[,2]==t_Bu[i,2]))
    count = c(count, length(w))
  }
  
  size = count^0.2
  size = size*3/max(size)
  main = concat(c(names(information1)[s]," ", "B cells"))
  cols = rep("blue", length(p_value))
  cols[which(p_value<0.05)] = "red"
  plot(log10(t_Bu+1), cex = size, pch = 21, bg = add.alpha(cols, alpha = 0.5), col = "NA", xlim = log10(range(t_Bu)+1), ylim = log10(range(t_Bu)+1), xlab = "% clone size in singlets (log10+1)",ylab = "% clone size in singlets (log10+1)", main= main)
  segments(1e-1000,1e-1000,1000,1000,col = add.alpha("grey", alpha = 0.5), lwd = 2, lty = 2)
  
  
  t_T = cbind(tt*0, tt*0)
  t_T1 = cbind(tt*0, tt*0)
  colnames(t_T) = c("singlet","doublet")
  colnames(t_T1) = c("singlet","doublet")
  doublet_tclone= table(clone_TCR[BCR_TCR_sub])
  doublet_tclone = doublet_tclone[grep("-", names(doublet_tclone), invert = T)]
  t_T[names(doublet_tclone),"doublet"] = doublet_tclone*100/sum(doublet_tclone)
  t_T1[names(doublet_tclone),"doublet"] = doublet_tclone
  singlet_tclone= table(clone_TCR[setdiff(ids_sample,BCR_TCR_sub)])
  singlet_tclone = singlet_tclone[grep("-", names(singlet_tclone), invert = T)]
  t_T[names(singlet_tclone),"singlet"] = singlet_tclone*100/sum(singlet_tclone)
  t_T1[names(singlet_tclone),"singlet"] = singlet_tclone
  
  t_Tu = unique(t_T1)
  count = NULL
  p_value = NULL
  for(i in c(1:length(t_Tu[,1]))){
    w = intersect(which(t_T1[,1]==t_Tu[i,1]), which(t_T1[,2]==t_Tu[i,2]))
    count = c(count, length(w))
    x =cbind(t_Tu[i,], colSums(t_Tu[-i,]))
    p_value = c(p_value, fisher.test(x, alternative = "two")$p.value)
  }
  sort(p_value)
  
  t_Tu = unique(t_T)
  count = NULL
  for(i in c(1:length(t_Tu[,1]))){
    w = intersect(which(t_T[,1]==t_Tu[i,1]), which(t_T[,2]==t_Tu[i,2]))
    count = c(count, length(w))
  }
  size = count^0.2
  size = size*3/max(size)
  main = concat(c(names(information1)[s]," ", "T cells"))
  cols = rep("blue", length(p_value))
  cols[which(p_value<0.05)] = "red"
  plot(log10(t_Tu+1), cex = size, pch = 21, bg = add.alpha(cols, alpha = 0.5), col = "NA", xlim = log10(range(t_Tu)+1), ylim = log10(range(t_Tu)+1), xlab = "% clone size in singlets (log10+1)",ylab = "% clone size in singlets (log10+1)", main= main)
  segments(1e-1000,1e-1000,1000,1000,col = add.alpha("grey", alpha = 0.5), lwd = 2, lty = 2)

}
dev.off()


### for potting out the networks and 
Plot_networks<-function(mVDJ, information1, outdir){
  m_VDJ_BCR = mVDJ$BCR
  m_VDJ_TCR = mVDJ$TCR
  parameters_degree_centrality = NULL
  size_list = NULL
  parameters_mean_sizes_adjacent= NULL
  parameters_max_sizes_adjacent= NULL
  fileout1=concat(c(outdir, "Enrichment_of_B_T_cell_doublet_network_", batch,".pdf"))
  w=4.5
  pdf(file=fileout1, height=w*1.1*1, width=w*1.3*3)
  par(mfrow= c(1,3), mar = c(12,5,5,3))
  mult = 1.025
  for(s in c(1:length(information1))){
    sample = names(information1)[s]
    print (s)
    BCR_TCR_sub = rownames(information1[[s]])
    clone_BCR = apply(cbind(m_VDJ_BCR[BCR_TCR_sub,"clone1"], m_VDJ_BCR[BCR_TCR_sub,"clone2"]), 1, paste, collapse = "_")
    clone_TCR = apply(cbind(m_VDJ_TCR[BCR_TCR_sub,"clone1"], m_VDJ_TCR[BCR_TCR_sub,"clone2"]), 1, paste, collapse = "_")
    
    clone_BCR = gsub(concat(c("||RCC-metP")),"B",clone_BCR, fixed = T)
    clone_TCR = gsub(concat(c("||RCC-metP")),"T",clone_TCR, fixed = T)
    clone_BCR = gsub(concat(c("B_")),":",clone_BCR, fixed = T)
    clone_TCR = gsub(concat(c("T_")),":",clone_TCR, fixed = T)
    w = intersect(grep("-", clone_BCR, invert = T), grep("-", clone_TCR, invert = T))
    clone_BCR = clone_BCR[w]
    clone_TCR = clone_TCR[w]
    
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
    
   
    library(yarrr)
    library(RColorBrewer)
    
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
    
    layout =layout_with_mds(g, dist = NULL, dim = 2, options = arpack_defaults)
    layout = layout_with_graphopt(g,niter = 200,start = layout, spring.length = 1)
    #layout = layout_with_graphopt(g,niter = 600)
    # layout1 =layout_in_circle(g)
    # layout1 = layout_with_graphopt(g,niter = 500)
    rownames(layout) = V(g)$name
    
    plot(g, layout=layout, edge.color= E(g)$color, main="", edge.width= edge_strength_plot,vertex.label.family="sans", edge.lty = 1, main = sample,xlim = c(-1, 1.5), ylim = c(-1, 1),cex.main =3) # vertex.label.cex = 0.4
    title(main = sample)
    layout2 = layout
    layout2[,1] = layout2[,1]/max(abs(layout2[,1]))
    layout2[,2] = (layout2[,2]/max(abs(layout2[,2])))
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
      text(x_lab_pos+0.15, y = y_lab_pos, label = plot_ids1[order], cex = 1.2, col = "darkblue")
      # text(layout2[plot_ids,], label = plot_ids)		
      for(i in c(1:length(order))){
        shift = 0
        if(s==2){shift = 0.02}
        segments(xy[order[i],1], xy[order[i],2]+shift, x_lab_pos-0.2, y_lab_pos[i], col = add.alpha("darkblue",alpha = 0.5), lty = 3, lwd = 1)
      }
    }
    
    legend("topleft",  legend= c("T cell","B cell"), pch= 21, bty="n",cex= 1.2, title="", pt.bg = col1,pt.cex = 1.1)
    scale = c(1,max(sizes))
    scale1 = scale ^0.5
    thresh1 = 90^0.5
    scale1 = (scale1*20/thresh1)-1
    legend("bottomleft",  legend= as.character(scale), pch= 21, bty="n",cex= 1, title="", pt.bg = "grey",pt.cex = scale1*0.4, col = "grey",y.intersp=2, xjust=0.5, x.intersp = 2)
    
    degree = igraph::degree(g, mode = "all")
    size_list = c(size_list, list(sizes))
    parameters_degree_centrality = c(parameters_degree_centrality, list(degree))
    
    mean_sizes_adjacent = NULL
    max_sizes_adjacent = NULL
    for(i in c(1:length(V(g)$name))){
      a = adjacent_vertices(g, V(g)$name[i], mode = "all")[[1]]
      mean_sizes_adjacent = c(mean_sizes_adjacent, mean(sizes[a]))
      max_sizes_adjacent = c(max_sizes_adjacent, max(sizes[a]))
    }
    parameters_mean_sizes_adjacent = c(parameters_mean_sizes_adjacent, list(mean_sizes_adjacent))
    parameters_max_sizes_adjacent = c(parameters_max_sizes_adjacent, list(max_sizes_adjacent))
    
  }
  dev.off()
  
  samples = names(information1)
  names(parameters_mean_sizes_adjacent) = samples
  names(parameters_max_sizes_adjacent) = samples
  names(parameters_degree_centrality) = samples
  names(size_list) = samples
  
 
  ###### centrality between B and T cells?
  b_centrality = NULL
  t_centrality = NULL
  centrality_all = NULL
  size_all = NULL
  b_proportion_high_centrality = NULL
  t_proportion_high_centrality = NULL
  for(i in c(1:length(parameters_degree_centrality))){
    c = parameters_degree_centrality[[i]]
    centrality_all = c(centrality_all, c)
    size_all = c(size_all, size_list[[i]])
    w = grep("B", names(c))
    b_centrality =c(b_centrality, c[w])
    b_proportion_high_centrality = c(b_proportion_high_centrality, length(which(c[w]>1))/length(c))
    w = grep("T", names(c))
    t_centrality =c(t_centrality, c[w])
    t_proportion_high_centrality = c(t_proportion_high_centrality, length(which(c[w]>1))/length(c))
  }
  
  library(yarrr)
  library(RColorBrewer)
  cols1 =  add.alpha (piratepal(palette = "google"), alpha = 0.95)
  cols =  add.alpha (cols1, alpha = 0.5)
  
  fileout1=concat(c(outdir, "Enrichment_of_B_T_cell_doublet_network_", batch,"_2.pdf"))
  w=2.7
  pdf(file=fileout1, height=w*1.1*2, width=w*1.1*2)
  par(mfrow= c(2,2), mar = c(5,5,5,5))
  library(psych)
  
  for(i in c(1:length(parameters_degree_centrality))){
    w =  grep("B",names(size_list[[i]]))
    
    x =size_list[[i]][w]
    y = parameters_degree_centrality[[i]][w]
    fit1 <- lm( y~x)
    p_value = signif(summary.lm(fit1)$ coefficients["x",4], digits = 3)
    cortest = corr.test(cbind(x,y),adjust="holm")
    pval = cortest $ p[1,2]
    rval = signif(cortest $ r[1,2], digits = 3)
    
    plot(x,y, pch = 21, col = cols1[1], bg = cols1[1], main = concat(c(names(parameters_degree_centrality)[i], "\nB cell clones")),
         xlab = "Clone size", ylab = "Clone centrality", log = "xy")
    segments(-10, -10, 10000, 10000, col = add.alpha("grey", alpha = 0.5), lty = 3)
    
    new = seq(from= min(x), to = max(x)*3, length = 1000)
    predicted.intervals <- predict(fit1, data.frame(x= new),interval='prediction',level = 0.85)
    cols2 =  add.alpha ("grey", alpha = 0.85)
    points(new,predicted.intervals[,1],col= cols2,lwd=4, type = "l")
    polygon(x = c(new, rev(new)),y=c(predicted.intervals[,2],rev(predicted.intervals[,3]) ),  col= add.alpha(cols2,alpha = 0.25),border = NA)
    points(x,y, pch = 21, bg=cols1[1], cex = 1, col= add.alpha("white", alpha = 0))	
    legend("topleft", concat(c("p= ",p_value,"\nr2= ",rval,"\nn=",length(w))), pch = 21,cex= 0.8, bty="n", pt.bg = "white", col = "white", pt.lwd = 2, text.font = 2)
    
    
    w =  grep("T",names(size_list[[i]]))
    x =size_list[[i]][w]
    y = parameters_degree_centrality[[i]][w]
    fit1 <- lm( y~x)
    p_value = signif(summary.lm(fit1)$ coefficients["x",4], digits = 3)
    cortest = corr.test(cbind(x,y),adjust="holm")
    pval = cortest $ p[1,2]
    rval = signif(cortest $ r[1,2], digits = 3)
   
    plot(x,y, pch = 21, col = cols1[1], bg = cols1[1], main = concat(c(names(parameters_degree_centrality)[i], "\nT cell clones")),
         xlab = "Clone size", ylab = "Clone centrality", log = "xy")
    segments(-10, -10, 10000, 10000, col = add.alpha("grey", alpha = 0.5), lty = 3)
    
    new = seq(from= min(x)-1, to = max(x)*3, length = 1000)
    predicted.intervals <- predict(fit1, data.frame(x= new),interval='prediction',level = 0.85)
    cols2 =  add.alpha ("grey", alpha = 0.85)
    points(new,predicted.intervals[,1],col= cols2,lwd=4, type = "l")
    polygon(x = c(new, rev(new)),y=c(predicted.intervals[,2],rev(predicted.intervals[,3]) ),  col= add.alpha(cols2,alpha = 0.25),border = NA)
    points(x,y, pch = 21, bg=cols1[1], cex = 1, col= add.alpha("white", alpha = 0))	
    legend("topleft", concat(c("p= ",p_value,"\nr2= ",rval,"\nn=",length(w))), pch = 21,cex= 0.8, bty="n", pt.bg = "white", col = "white", pt.lwd = 2, text.font = 2)
    
  }
  
  dev.off()
  
}


Plot_networks(mVDJ, information1, outdir)

