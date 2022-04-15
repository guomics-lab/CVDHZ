req.pcg.install <- function(pcg) {
  new <- pcg[!(pcg %in% installed.packages()[, "Package"])]
  if (length(new))
    install.packages(new, dependencies = T)
  sapply(pcg, require, ch = T)
  
}
req.pcg.install(c(    "readr",
                      "plyr",
                      "readxl",
                      "stringr",
                      "magrittr",
                      "magrittr",
                      "ggplot2",
                      "BiocManager",
                      "pheatmap",
                      "RColorBrewer",
                      "umap",
                      "Rtsne",
                      "vioplot"))

ge.na.ratio <- function(x){
  sum(is.na(x))/dim(x)[1]/dim(x)[2]
}
#as.data.frame(lapply(datexpr,as.numeric))
ge.convert.geneprot <- function(vector,type="uniprot"){
  keytypes(org.Hs.eg.db) 
  protID = bitr(vector, fromType="UNIPROT", toType=c("SYMBOL", "UNIPROT"),drop = F, OrgDb="org.Hs.eg.db")
  return(paste0(vector,"_",protID$SYMBOL[match(vector,protID$UNIPROT)]))
}

ge.convert.number <- function(data){
  return(as.data.frame(lapply(data, as.numeric)))
}

ge.color <- function(n){
  if(n<10){
    set.seed(21)
    m <- sample(9,n)
    color <- brewer.pal(9,"Set1")[m]
  }else{
    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    set.seed(10)
    color=sample(col_vector, n)
  }
  return(color)
}


ge.split <- function(data,split,which=1,start=NULL){
  if(is.null(start)){
    sapply(data,function(v){strsplit(v,split)[[1]][which]})
  }else{
    tmp <- sapply(data,function(v){strsplit(v,split)[[1]][1]})
    sapply(tmp,function(v){strsplit(v,start)[[1]][2]})
  }
}


ge.read <- function(data,colname = T){
  suffix <- strsplit(data,"\\.")[[1]][length(strsplit(data,"\\.")[[1]])]
  data2 <- switch(suffix,
                  txt=read.table(data,sep = "\t",header = colname,stringsAsFactors = F),
                  csv=read.csv(data,header = colname),
                  tsv=read_tsv(data,col_names = colname),
                  xlsx=read_xlsx(data,col_names = colname)
  )
  return(data2)
}

ge.readtable <- function(data,sep = "\t",header = T){
  read.table(data,sep = sep,header = header,stringsAsFactors = F)
}

ge.writetable <- function(data,filename ,sep = "\t",col.names = T,row.names = T,quote = F){
  write.table(data,filename,sep=sep,col.names = col.names,row.names = row.names,quote = quote)
}

NA_threshold_table <- function(matrix) {
  na_ratio_in_each_prot = apply(matrix, 1, function(x) {
    sum(is.na(x))/ncol(matrix)
  })
  
  temp = data.frame(sample = names(na_ratio_in_each_prot),
                    na_ratio = na_ratio_in_each_prot,
                    stringsAsFactors = F)
  
  Table1_na = sapply(10:1, function(x) {
    threshold = x/10
    
    prot_choose = temp$na_ratio <= threshold
    prot_num = sum(prot_choose)
    na = ge.na.ratio(matrix[prot_choose,]) %>% round(4)
    
    return(c(threshold = paste0(10*x,"%"),
             protein_num = prot_num,
             NA_ratio = paste0(100*na,"%")))
  }) %>% t()
  
  return(Table1_na)
}

ge.plot.density <- function(data){
  plot(density(na.omit(unlist(data))),main="")
  #axis(1,cex.axis = 3)
}

ge.remove.techrep <- function(data,pattern="_repB",method="mean"){
  repB <- names(data)[grepl(pattern, names(data))]
  for (i in repB) {
    repa <- str_split(i,pattern)[[1]][1]
    df1 <- data[,which(names(data) %in% c(repa,i))]
    data <- data[,-which(names(data) %in% c(repa,i))]
    new_mean <- apply(df1, 1, function(x){ifelse(sum(is.na(x))==ncol(df1),NA, mean(as.numeric(x),na.rm=T))} )
    data <- cbind(data,new_mean)
    names(data)[ncol(data)] <- repa
  }
  return(data)
}

ge.plot.techrep.correlation <- function(cor1,cor2,name="pearson_correlation"){
  pdf(paste0(name,".pdf"))
  r <- cor(cor1, cor2, use = "pairwise.complete.obs")   
  smoothScatter(cor1, cor2, nrpoints = 100,cex = 2,
                colramp = colorRampPalette(c(blues9,"orange", "red")),
                main = name, xlab = "repA", ylab = "repB")
  abline(lm(cor1 ~ cor2), col="red", lwd=2, lty=2)
  text(min(cor1,na.rm = T)*1.3,max(cor2,na.rm = T)*0.8,labels =paste0( "r =", as.character(round(r,4))),cex = 1.2)
  dev.off()
}

ge.plot.pool.correlation <- function(data,name="bio_cor",method="circle",height = 7, width = 7,cor=1){
  library(corrplot)
  df_cor <- data.frame(data)
  pdf(paste0(name,".pdf"),height = height, width = width)
  mycor=cor(df_cor, use = "pairwise.complete.obs")
  if(cor==2){corrplot(mycor,
           method = "circle",
           type = "upper",
           is.corr = T,addCoef.col = "grey",
           tl.col = "black",
           tl.cex = 1.3,
           cl.cex = 1.2,
           number.cex = 1.5,
           mar=c(2,0,5,2))}else{corrplot(mycor, method=method,type = "upper",tl.col = "black",tl.srt = 45, tl.cex = 1.5)}
  dev.off()
}




ge.plot.boxplot <- function(data,x,y,type,filename,title="boxplot"){
  a <- ggplot(data=data, aes(x =x, y =y ,color=type,group=type)) +
    geom_jitter(alpha = 0.3,size=3) +
    geom_boxplot(alpha = .5,size=1)+
    labs(x="sample",y="value",fill= "type")+
    ggtitle(title)+
    theme_bw() + 
    theme(panel.border = element_blank())+
    theme(axis.line = element_line(size=1, colour = "black")) +
    theme(panel.grid =element_blank())+  
    theme(axis.text = element_text(size = 15,colour = "black"),text = element_text(size = 15,colour = "black"))+
    theme(axis.text.x = element_text( hjust = 1,angle = 45))
  ggsave(paste0(filename, ".pdf"),plot=a,width=8,height=8)
}

# p <- ggboxplot(df1, x="dose", y="len", color = "dose", 
#                palette = c("#00AFBB", "#E7B800", "#FC4E07"), 
#                add = "jitter", shape="dose")
# my_comparisons <- list(c("0.5", "1"), c("1", "2"), c("0.5", "2"))
# p+stat_compare_means(comparisons = my_comparisons)+
#   stat_compare_means(label.y = 50)



ge.plot.pca <- function(data,type,title="",ptColors=NULL,label2=NULL,width=12,height=8,quan=NULL){ 
  M <- t(data)
  M[is.na(M)] <- 0
  M <- apply(M,2,function(v){(v-mean(v,na.rm=T))/sd(v,na.rm=T)})
  clnames <- row.names(data)
  set.seed(10)
  m1 <- prcomp(M);
  Y  <- scale(M, m1$center, m1$scale) %*% m1$rotation 
  Y  <- Y[,c(1,2)]
  
  Y <- data.frame(Y,type);
  colnames(Y) <- c("PC1","PC2","label")
  eigs <- m1$sdev^2
  percentages <- eigs[1:2] / sum(eigs)
  p <- ggplot(Y, aes(x=PC1, y=PC2, colour=label)) + geom_point(size=4)
  p <- p + theme(  panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   axis.text = element_text(size = 30,color = "black"),
                   panel.border = element_blank(),
                   axis.line.x = element_line(color="black", size = 0.25),
                   axis.line.y = element_line(color="black", size = 0.25),
                   plot.title   = element_text(size=30),
                   axis.title   =element_text(size=30),
                   panel.background = element_blank())
  
  strLabx <- sprintf("PC1(%4.2f%%)",percentages[1]*100)
  p <- p +  labs(x =strLabx,y = sprintf("PC2(%4.2f%%)",percentages[2]*100),
                 title =sprintf("PCA:%d features",length(clnames)))
  if(!is.null(ptColors)){
    p <- p +   scale_colour_manual(values=ptColors)
  }
  if(!is.null(label2)){
    p <- p +   geom_text(aes(label=label2,vjust = -0.8, hjust = 0.5,size=0.5),show.legend = FALSE)
  }
  # 添加椭圆
  if(!is.null(quan)){
  p<-p+stat_ellipse(level = 0.95, show.legend = F)
  }
  ggsave(paste0(title,"_PCA.pdf"),plot =p ,width=width,height=height,device = NULL)
}


ge.plot.tsne <- function(data,type,title="",label=NULL){
  col2=ge.color(length(unique(type)))
  cl <- data.frame(col2,row.names =unique(type),stringsAsFactors = F)
  cl2 <- cl[match(type,row.names(cl)),1]
  
  df10 <- data
  df10[is.na(df10)] <- 0
  df10 <- t(apply(df10, 1, scale))
  colnames(df10) <- type
  set.seed(10)
  df11.tsne <- Rtsne(t(df10), dims = 2, perplexity = (ncol(data)-1)/3-1, verbose = T , check_duplicates = FALSE)
  pdf(paste0(title,"_TNSE.pdf"))
  if(is.null(label)){
  plot(df11.tsne$Y,col=cl2, main = "tsne", pch = 20,cex=2,cex.axis=2,cex.lab=2)
  legend("topright",legend=row.names(cl), fill=cl$col2, lty=1,lwd=1)
  }else{
    plot(df11.tsne$Y,col=cl2, main = "tsne", pch = 20,cex=2,cex.axis=2,cex.lab=2)
    text(df11.tsne$Y,pos = 1, labels = label, col= "DimGrey",cex = 0.8)
    legend("topright",legend=row.names(cl), fill=cl$col2, lty=1,lwd=1)
  }
  dev.off()
}


ge.plot.umap<- function(data,type,title="",label=NULL){
  col2=ge.color(length(unique(type)))
  cl <- data.frame(col2,row.names =unique(type),stringsAsFactors = F)
  cl2 <- cl[match(type,row.names(cl)),1]
  
  df10 <- data
  df10[is.na(df10)] <- 0
  df10 <- t(apply(df10, 1, scale))
  colnames(df10) <- type
  set.seed(10)
  df.umap <- umap(t(df10),n_neighbors=ncol(data)-1)
  pdf(paste0(title,"_UMAP.pdf"))
  if(is.null(label)){
  plot(df.umap$layout,col = cl2, main = "umap", pch = 20,cex=2,cex.axis=2,cex.lab=2)
    
  legend("topright",legend=row.names(cl), fill=cl$col2, lty=1,lwd=1)
  }else{
    plot(df.umap$layout,col = cl2, main = "umap", pch = 20,cex=2,cex.axis=2,cex.lab=2)
    text(df.umap$layout, pos = 1, labels = label, col= "DimGrey",cex = 0.8)
    legend("topright",legend=row.names(cl), fill=cl$col2, lty=1,lwd=1)
  }
  dev.off()
}



ge.plot.volcano <- function(data, group1, group2, fc= 1, pvalue = 0.05, str1= "grp1",str2= "grp2",pair=F,adjust.bool=T) {
  fc <- as.numeric(fc)
  pvalue <- as.numeric(pvalue)
  len.grp1 <- length(group1)
  len.grp2 <- length(group2)
  data <- data[,c(group1,group2)]
  data <- data[apply(data, 1, function(x){sum(!is.na(x))})>0,]
  df8 <- 2^data
  df8[is.na(df8)] <- min(df8,na.rm = T)*0.8
  df8$fd <-
    apply(df8, 1, function(x)
      log2((mean(x[1:len.grp1], na.rm = T) / mean(x[(len.grp1+1):ncol(data)], na.rm = T))))
  x <- c(0.0, 0.0)
  df9 <- data
  df9[is.na(df9)] <- min(df9,na.rm = T)*0.8
  
  df8$P_value <- apply(df9,1,function(y) {p_try = tryCatch(t.test(y[1:len.grp1],
                                             y[(len.grp1+1):ncol(data)],
                                             paired = pair,
                                             var.equal = F)$p.value,
                                      error = function(x) NA)})
  
  
  if(adjust.bool){
    df8$P_value_adjust <- p.adjust(df8$P_value, method = "BH")
    df8$P <- df8$P_value_adjust
    y.lab <- "-log10 (adjust P)"
  }else{
    df8$P <- df8$P_value
    y.lab <- "-log10 (P_value)"
  }
  pdf(paste0(str1, "_", str2, "_volcano.pdf"),
      width = 4,
      height = 4,
  )
  plot(
    df8$fd,
    -log10(df8$P),
    col = "#00000033",
    pch = 19,
    xlab = paste("log2 (fold change)"),
    ylab = y.lab,
    #xlim = c(-4, 4),
    main = paste0(str1, " vs ", str2)
  )
  
  up <- subset(df8, df8$P < pvalue & df8$fd > fc)
  down <- subset(df8, df8$P < pvalue & df8$fd < -1*fc)
  write.csv(up, file = paste0(str1, "_", str2, "_up_volcano.csv"))
  write.csv(down, file = paste0(str1, "_", str2, "_dw_volcano.csv"))
  write.csv(df8, file = paste0(str1, "_", str2, "_all_volcano.csv"))
  
  points(
    up$fd,
    -log10(up$P),
    col = 1,
    bg = brewer.pal(9, "YlOrRd")[6],
    pch = 21,
    cex = 1.5
  )
  points(
    down$fd,
    -log10(down$P),
    col = 1,
    bg = brewer.pal(11, "RdBu")[9],
    pch = 21,
    cex = 1.5
  )
  abline(
    h = -log10(pvalue),
    v = c(-1*fc, fc),
    lty = 2,
    lwd = 1
  )
  dev.off()
  print(paste0("up:",nrow(up)))
  print(paste0("down:",nrow(down)))
  return(c(row.names(up),row.names(down)))
}

ge.plot.bar <- function(data,sample,value,group=group,title="",xlab="sample",ylab="value"){
  a <- ggplot(data,aes(x=sample,y=value,group=group))+ 
    geom_bar(position = "dodge",stat = "identity",width =0.8,alpha=0.8,aes(fill=group))+
    ggtitle(paste0(title,"_barplot"))+
    xlab(xlab)+
    ylab(ylab)+
    theme(legend.text = element_text(size = 15,color = "black"),legend.position = 'top',
          legend.title = element_text(size=15,color="black") ,
          panel.grid.major =element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))+
    theme(panel.grid =element_blank())+
    theme(axis.text = element_text(size = 10,color = "black"))+
    theme(axis.text.x = element_text( hjust = 1,angle = 45))+
    theme(plot.subtitle=element_text(size=30, hjust=0, color="black"))+
    theme(axis.title.x=element_text(size=17, hjust=0.5, color="black"))+
    theme(axis.title.y=element_text(size=17, hjust=0.5, color="black")) + geom_text(aes(x=sample,label=value,vjust = -0.8, hjust = 0.5),position = "dodge",stat = "identity",show.legend = FALSE)
  ggsave(paste0(title,"_barplot.pdf"),plot=a,width=10,height=8)
}


ge.plot.line <- function(data,sample,value,group=group,title="",xlab="sample",ylab="value"){
  a <- ggplot(data,aes(x=sample,y=value,group=group,color=group))+ 
    geom_line()+
    geom_point()+
    ggtitle(paste0(title,"_lineplot"))+
    xlab(xlab)+
    ylab(ylab)+
    theme(legend.text = element_text(size = 15,color = "black"),legend.position = 'top',
          legend.title = element_text(size=15,color="black") ,
          panel.grid.major =element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))+
    theme(panel.grid =element_blank())+
    theme(axis.text = element_text(size = 10,color = "black"))+
    theme(axis.text.x = element_text( hjust = 1,angle = 45))+
    theme(plot.subtitle=element_text(size=30, hjust=0, color="black"))+
    theme(axis.title.x=element_text(size=17, hjust=0.5, color="black"))+
    theme(axis.title.y=element_text(size=17, hjust=0.5, color="black"))+ geom_text(aes(label=group,vjust = -0.8, hjust = 0.5),show.legend = FALSE)
  ggsave(paste0(title,"_lineplot.png"),plot=a,width=10,height=8)
}

ge.plot.alluvial <- function(data,axis1,axis2,title="",axis3 = NULL,axis4 = NULL,axis5 = NULL,fill=axis1){
  a <- ggplot(data = data,
              aes(axis1 = axis1, axis2 = axis2, axis3 = axis3, axis4 = axis4, axis5 = axis5) )+
    scale_x_discrete(limits = c("axis1", "axis2"), expand = c(.1, .05)) +
    geom_alluvium(aes(fill = factor(fill))) +
    geom_stratum( ) + 
    geom_text(stat = "stratum",aes(label = after_stat(stratum))) + # 
    theme_minimal() +
    labs(fill="group")+
    theme(panel.grid =element_blank()) +   
    theme(axis.text.y  = element_blank()) +   
    theme(axis.ticks = element_blank())+   
    scale_fill_brewer(palette = "Dark2") 
  #scale_color_brewer(type = "qual", palette = "Set3") 
  ggsave(paste0(title,"_alluvial.pdf"),plot=a,width=10,height=8)
}


ge.plot.vioplot <- function(data.list,title="",xlab="sample",ylab="value",color=NULL,width = 7,height = 7){
  df2 <-data.list  
  pdf(paste0(title, "_violin.pdf"),width = width,height = height)
  if(is.null(color)){
    vioplot(df2 ,
            areaEqual=FALSE,
            lineCol=rep("black",length(df2)),
            border=rep("black",length(df2)),
            names=paste0(names(df2),"\n",lapply(df2, function(x){round(median(x,na.rm = T),4)})),
            main=title, xlab=xlab, ylab=ylab,plotCentre = "point")
  }else{
    vioplot(df2 ,
            areaEqual=FALSE,
            rectCol= color, 
            col= color,
            lineCol=rep("black",length(df2)),
            border=rep("black",length(df2)),
            names=paste0(names(df2),"\n",lapply(df2, function(x){round(median(x,na.rm = T),4)})),
            main=title, xlab=xlab, ylab=ylab,plotCentre = "point")
  }
  dev.off()
}



ge.mfuzz.cselection <- function(data,range=seq(5,50,5),repeats = 5){
  df3a<-as.matrix(data)
  df3Ex<- ExpressionSet(assayData = df3a)
  if(interactive()){
    df3F <- filter.NA(df3Ex)
    df3F <- fill.NA(df3F)
    df3F <- standardise(df3F)
  }
  
  df3F <- filter.NA(df3F)
  m<-mestimate(df3F)
  cselection(df3F,m=m,crange = range,repeats = repeats,visu = T)
  return(df3F)
}

ge.mfuzz.getresult <- function(expressionSet, pic,time.label,filename,anova=F,alldata=NULL,type=NULL){
  set.seed(10)
  cl <- mfuzz(expressionSet,c=pic,m=1.25)
  dir.create(path=filename,recursive = TRUE)
  pdf(paste0(filename,".pdf"),width = 10, height  = 8)
  mfuzz.plot2(expressionSet, cl=cl,time.labels=time.label,mfrow=c(4,4),centre=TRUE,x11=F,centre.lwd=0.2)#min.mem=0.99
  dev.off()
  
  for(i in 1:pic){
    potname<-names(cl$cluster[unname(cl$cluster)==i])
    write.csv(cl[[4]][potname,i],paste0(filename,"/mfuzz_",i,".csv"))
  }
  if(anova){
    for(ii in 1:pic){
      potname<-names(cl$cluster[unname(cl$cluster)==ii])
      tmp <- data.frame(label=as.factor(type),t(alldata[potname,]))
      anova <- c()
      for (n in 1:length(potname)) {
        aov<- summary(aov(tmp[,n+1] ~ label,data=tmp))[[1]]$`Pr(>F)`[1] %>% format(digits = 3, scientific = FALSE)
        anova <- c(anova,aov)
      }
      anova.adjust <-p.adjust(anova, method="BH")
      newdf <- data.frame(prot=names(tmp)[-1],anova,anova.adjust)
      newdf2 <- newdf[newdf$anova.adjust<0.05,]
      write.csv(newdf,paste0(filename,"/mfuzz_anova_",ii,".csv"),row.names = F)
    }
  }
}


ge.upset <- function(){
  listinput <- list(AD_NC = prot.AD_NC,
                    ADmild_NC = prot.ADmild_NC,
                    ALS_NC = prot.ALS_NC,
                    HD_NC = prot.HD_NC,
                    MCI_NC = prot.MCI_NC)
  
  
  library(UpSetR)
  pdf(file='venn.pdf',height = 8,width = 8)
  upset(fromList(listinput),nsets = 5, order.by = "freq")
  dev.off()
  
  
  for (i in 2:length(listinput)) {
    com <- combn(length(listinput),i)
    for (n in 1:ncol(com)) {
      all <- row.names(df2)
      fl <- ""
      for (m in 1:nrow(com)) {
        all <-  intersect(all,listinput[[com[m,n]]])
        fl <- paste0(fl,"_",names(listinput)[[com[m,n]]])
      }
      write.csv(all,paste0("CSF/overlap/",fl,".csv"))
    }
  }
}



ge.plot.heatmap <- function(data,df.lable,scale = "row", cluster_rows = T, cluster_cols = F,show_rownames = T, show_colnames = T, filename = "heatmap.pdf",width=10,height=20,fontsize_row=NULL){
  
    ann_col <- df.lable
    
    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    colors <- list()
    for (n in 1:ncol(ann_col)) {
      nm <- length(levels(factor(ann_col[,n])))
      if(nm<11){
        color <- c(brewer.pal(10,"Paired")[1:nm])
      }else{
        set.seed(10)
        color <- sample(col_vector, nm)
      }
      colors[[n]] <- color
      names(colors[[n]]) <- levels(factor(ann_col[,n]))
    }
    names(colors) <- names(ann_col)

    a <- pheatmap(data, color = c(brewer.pal(11,"RdYlBu")[9:7],"azure1",brewer.pal(11,"RdYlBu")[4:1]), #fontsize_col = 8,
                  annotation_col = ann_col,fontsize_row=fontsize_row,
                  annotation_colors = colors,na_col = "grey60",scale = scale,
                  cluster_rows = cluster_rows, cluster_cols = cluster_cols,show_rownames = show_rownames, show_colnames = show_colnames, 
                  filename = filename,width=width,height=height)
    
}


ge.plot.enrich <- function(genelist,type="GO",keyType = 'ENTREZID',organism = 'hsa',pvalueCutoff = 0.05,qvalueCutoff = 0.2,filename=""){
  if(type=="GO"){
   go <- enrichGO(genelist, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = pvalueCutoff,  qvalueCutoff = qvalueCutoff,keyType = keyType)
   
  pdf(paste0(filename,"_GO_barplot.pdf"))
  barplot(go,showCategory=20,drop=T)
  dev.off()
  
  pdf(paste0(filename,"_GO_dotplot.pdf"))
  dotplot(go,showCategory=50)
  dev.off()
  }else if(type=="KEGG"){
   kegg <- enrichKEGG(genelist, organism = organism, keyType = 'kegg', pvalueCutoff = pvalueCutoff,pAdjustMethod = 'BH', minGSSize = 10,maxGSSize = 500,qvalueCutoff = qvalueCutoff,use_internal_data = FALSE)
   
   pdf(paste0(filename,"_KEGG_dotplot.pdf"))
   dotplot(kegg, showCategory=30)
   dev.off()
  }
}

auto_preprocess <-
  function(filename = "peptides.txt",
           tech_rep_f = NULL,
           batchf = NULL,
           psep = "\t",
           tsep = "\t",
           pheader = TRUE,
           theader = FALSE,
           bheader = TRUE,
           bsep = "\t",
           lr_top3 = "top3") {
    t1 <- proc.time()
    pep.data <-
      read.table(
        filename,
        header = pheader,
        sep = psep,
        check.names = F,
        quote = "",
        fill = TRUE
      )
    
    pep.data[is.na(pep.data)]<-NA
    pep.data<-pep.data[complete.cases(pep.data[,1]),]
    pep.data <- pep.data[!grepl("^1/CON", pep.data[, 2], fixed = F), ]
    
    pep.data[pep.data == 0] <- NA
    pep <- as.vector(as.matrix(pep.data[, 3:ncol(pep.data)]))
    #print(paste("missing rate is: ",sum(is.na(pep))/length(pep),sep=""))
    
    #log2 transform
    pep.data.log2 <- pep.data
    rownames(pep.data.log2) <- pep.data.log2[, 1]
    pep.data.log2 <- pep.data.log2[, -1]
    pep.data.log2[, 2:ncol(pep.data.log2)] <-
      log2(as.matrix(pep.data.log2[, 2:ncol(pep.data.log2)]))
    if(ncol(pep.data.log2)==2){
      pep.data.log2.group<-split(pep.data.log2,pep.data.log2[,1])
      prot.t<-lapply(pep.data.log2.group, function(x){
        x<-x[order(x[,2],decreasing = T),]
        y<-0
        if(nrow(x)>3){
          y<-mean(x[1:3,2],na.rm = T)
        }else{
          y<-mean(x[1:nrow(x),2],na.rm = T) 
        }
        return(y)
      })
      prot.t.d<-do.call(rbind,prot.t)
      prot.d<-data.frame(prot=rownames(prot.t.d),intensity=prot.t.d[,1])
      colnames(prot.d)[2]<-colnames(pep.data)[3]
      return(prot.d)
    }
    #R preprocessCore normalize.quantiles()
    
    pep.data.log2.qn = normalize.quantiles(as.matrix(pep.data.log2[, 2:ncol(pep.data.log2)]))
    colnames(pep.data.log2.qn) <- colnames(pep.data.log2)[-1]
    rownames(pep.data.log2.qn) <- rownames(pep.data.log2)
    #write.table(pep.data.log2.qn, "data/qn_log2_pep.txt",col.names=T,row.names=T,quote = F,sep = "\t",na = "NA")
    ####technical imputation
    #no technical imputaiton
    if(is.null(tech_rep_f)){
      data.tech.rep<-cbind(pep.data[,1:2],pep.data.log2.qn)
    }
    ##technical imputation
    else{    
      #read annotation
      anno <-
        read.table(tech_rep_f,
                   header = theader,
                   sep = tsep,
                   check.names = F)
      if(nrow(anno)!=ncol(pep.data.log2.qn))
        return("Number of samples do not match in peptide file and technical replicate file!")
      else if(length(intersect(unlist(anno[,1]),colnames(pep.data.log2.qn)))!=nrow(anno)){
        return("Samples do not match in peptide file and technical replicate file!")
      }
      else{
        
      }
      colnames(anno) <- c("V1", "V2")
      rownames(anno) <- anno$V1
      anno <- anno[colnames(pep.data.log2.qn), ]
      replicates <- split(as.vector(anno$V1), as.factor(anno$V2))
      replicates_impute <- function(x) {
        index <- which(is.na(x))
        x.median <- median(x, na.rm = T)
        x[index] <- x.median
        return(x)
      }
      data <- pep.data.log2.qn
      
      data.reps <- data.frame()
      for (reps in replicates) {
        if (length(reps) > 1) {
          d <- apply(data[, reps], 1, function(x) {
            return(replicates_impute(x))
          })
        }
        else {
          d <- data[, reps]
          d <- t(as.matrix(d))
          rownames(d) <- reps
        }
        data.reps <- rbind(data.reps, d)
      }
      data.tech.rep <- t(data.reps)
      data.tech.rep[is.nan(data.tech.rep)] <- NA
      data.tech.rep <- cbind(pep.data[, 1:2], data.tech.rep)
      #write.table(data.tech.rep, "data/data.tech.rep.txt",col.names=T,row.names=F,quote = F,sep = "\t",na = "NA")
    }
    
    if (!is.null(batchf)) {
      #print(batchf)
      
      batch <- read.table(batchf, header = bheader, sep = bsep)
      if(nrow(batch)!=ncol(data.tech.rep)-2)
        return("Number of samples do not match in peptide file and batch file!")
      else if(length(intersect(unlist(batch[,1]),colnames(data.tech.rep[,-c(1:2)])))!=nrow(batch)){
        return("Samples do not match in peptide file and batch file!")
      }else{
        data.tech.rep <- mycombat(data.tech.rep, batch)
        temp<-data.tech.rep
        temp<-data.tech.rep[,-c(1,2)]
        temp[temp<0]<-NA
        data.tech.rep[,-c(1,2)]<-temp
        rm(temp)
      }
      
      
    }
    
    ###order
    
    data <- data.tech.rep
    colnames(data)[1:2] <- c("tg", "prot")
    
    n = ncol(data)
    pep2 <- apply(data[, -c(1, 2)], 1, function(x) {
      NAs <- length(which(is.na(x)))
      meanexpr1 <- sum(as.numeric(x), na.rm = TRUE) / (n - NAs)
      meanexpr2 <- sum(as.numeric(x), na.rm = TRUE) / n
      d <- c(NAs, meanexpr1, meanexpr2)
      return(d)
    })
    pep2 <- t(pep2)
    colnames(pep2) <- c("NAs", "meanexpr1", "meanexpr2")
    pep_expr = cbind(data[, 1], pep2, data[, c(-1)])
    ########################################threee order methods ##################just choose only one
    ##order by pg ,#NA,intesity
    pep_order = pep_expr[order(pep_expr[, 5], pep_expr[, 2], -pep_expr[, 3]), ]
    colnames(pep_order)[1] <- "tg"
    ##order by pg intensity,#NA(excluding NAs)
    # pep_order=pep_expr[order(pep_expr[,5],-pep_expr[,3],pep_expr[,2]),]
    # colnames(pep_order)[1]<-"tg"
    # ###order by pg intensity,#NA(all samples)
    # pep_order=pep_expr[order(pep_expr[,5],-pep_expr[,4],pep_expr[,2]),]
    # colnames(pep_order)[1]<-"tg"
    
    #######################################################################################################
    pep_order2 <- pep_order[, c(-2, -3, -4)]
    
    ###select top 3 pep
    pre.prot = ""
    same.count = 0
    pep_order2.top3 <- data.frame()
    for (i in 1:nrow(pep_order2)) {
      if (pre.prot == as.vector(pep_order2[i, "prot"])) {
        same.count = same.count + 1
      }
      else {
        pre.prot = as.vector(pep_order2[i, "prot"])
        same.count = 1
      }
      if (same.count <= 3) {
        pep_order2.top3 <- rbind(pep_order2.top3, pep_order2[i, ])
      }
    }
    ##protetypic proteins are saved
    # pep_order2.top3 <-
    # pep_order2.top3[grep("^1/", pep_order2.top3$prot), ]
    
    pep_order2.top3 <-
      pep_order2.top3[c("prot", "tg", colnames(pep_order2.top3)[3:ncol(pep_order2.top3)])]
    pep_order2.top3[pep_order2.top3 == 0] <- NA
    if(lr_top3=="top3"){
      
      top3.prot.group<-split(pep_order2.top3,pep_order2.top3$prot)
      top3.mean<-lapply(top3.prot.group, function(l){
        #apply(l[,3:ncol(l)],2,mean,na.rm=T)
        apply(l[,3:ncol(l)],2,function(x){round(mean(x,na.rm=T),2)})
        
        # apply(l[,3:ncol(l)],2,function(v){
        #   if(is.nan(mean(v,na.rm=T)))
        #     return('')
        #   else return(mean(v,na.rm=T))
        # })
      })
      top3.mean<-do.call(rbind,top3.mean)
      top3.mean<-data.frame(top3.mean)
      top3.mean[is.na(top3.mean)]<-''
      top3.mean<-cbind(prot=rownames(top3.mean),top3.mean)
      
      return(top3.mean)
    }else{
      #write.table(pep_order2.top3, paste("data/",Sys.Date(),"pep.top3.txt",sep = ""),row.names = F,  quote = F,sep = "\t",na = "NA")
      
      #############lr for pep2prot
      
      prot.matrix <- pep2prot(pep_order2.top3)
      prot.matrix <- prot.matrix[, -2]
      
      return(prot.matrix)
    }
    
  }


library(speedglm)
mylm <- function(x, y) {
  mylr <- speedlm(y ~ x + 1, data = data.frame(y, x))
  mylr.summary <- summary(mylr)
  p=10
  if(nrow(mylr.summary$coefficients)>1)
    p = mylr.summary$coefficients[2, "p.value"]
  
  return(
    list(
      mylr = mylr,
      rSqured = mylr.summary$r.squared,
      fstatistic = mylr.summary$fstatistic[1],
      p=p
    )
  )
}
pep2prot <- function(filename) {
  top3 <- auto_preprocess(filename)
  prot.group <- split(top3, top3$prot)
  prot.matrix <- lapply(prot.group, function(d) {
    Y <- d[1, ]
    if (is.null(d)) {
    }
    else if (nrow(d) < 2) {
      
    }
    else{
      for (i in 2:nrow(d)) {
        x <- d[i, -c(1:2)]
        y <- d[1, -c(1:2)]
        
        xy.index <- which((!is.na(x)) & (!is.na(y)))
        xy0.index <- which((!is.na(x)) & (is.na(y)))
        if (length(xy.index) < 4 | all(unlist(x[xy.index]) == unlist(y[xy.index]))) {
          next
        }
        lr.results <- mylm(unlist(x[xy.index]), unlist(y[xy.index]))
        if(is.na(lr.results$p)|is.na(lr.results$rSqured)){
          next
        }else if(lr.results$rSqured > 0.36 &
                 lr.results$p < 0.05 & length(xy0.index) > 0 & length(xy.index) > 2) {
          tryCatch(
            y.predict <-
              predict(lr.results$mylr, new = data.frame(x = unlist(x[xy0.index]))),
            error = function (e) {
              print(d[i, ])
              stop()
            }
          )
          
          y.predict[y.predict < 4] <- NA
          Y[xy0.index + 2] <- y.predict
        }
      }
    }
    return(Y)
  })
  
  prot.matrix <- do.call(rbind, prot.matrix)
  return(prot.matrix)
}
