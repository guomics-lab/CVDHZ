# pacman::p_load(readxl,magrittr,
#                RColorBrewer,
#                Rtsne, umap,ggthemes,
#                pheatmap,vioplot,
#                ggpubr, ggplot2,
#                corrplot, stringr,
#                barplot)

#################################################################matrix process#################################################################
linda.summary.matrix <- function(matrix){
   cat("Sample size: ",
     ncol(matrix),
     "\nNumber of protein: ",
     nrow(matrix),
     "\nMissing rate: ",
     round(sum(is.na(matrix))/dim(matrix)[1]/dim(matrix)[2],4),
     "\nMaximum value: ",
     round(fivenum(as.matrix(matrix))[5],4),
     "\nMinimum value: ",
     round(fivenum(as.matrix(matrix))[1],4),
     "\nIQR: ",
     round(IQR(as.matrix(matrix),na.rm = T),4) 
   )
 
}
####################################################read data###############################
ge.read <- function(data,colname = T){
  data=input
  suffix <- strsplit(data,"\\.")[[1]][length(strsplit(data,"\\.")[[1]])]
  data2 <- switch(suffix,
                  txt=read.table(data,sep = "\t",header = colname,stringsAsFactors = F),
                  csv=read.csv(data,header = colname),
                  tsv=read_tsv(data,col_names = colname),
                  xlsx=read_xlsx(data,col_names = colname)
  )
  return(data2)
}
############################remove complete NA rows or column####################
linda.removeblanksamples <- function(matrix){
  newR <- matrix[apply(matrix, 1, function(x) {sum(!is.na(x))>0}),];
  newC <- newR[,apply(newR, 2, function(x) {sum(!is.na(x))>0})];
  return(newC)
}
#####################################clear DIA/DDA matrix##################################
linda.ClearRawmat <- function(matrix){   
  newR <- matrix[apply(matrix, 1, function(x) {sum(!is.na(x))>0}),];
  newC <- newR[,apply(newR, 2, function(x) {sum(!is.na(x))>0})];
  Dat<- newC[!grepl(pattern = ";",newC[,1]),]
  Dat<- Dat[!grepl(pattern = "^SWISS|^TRE|^H-|ENSEMBL",Dat[,1]),]
  Dat<- data.frame(Dat[,-1],row.names=Dat[,1])
}
##############################################################0-1标准化##############################################################
norm.mat <- function(matrix){
  apply(matrix, 1, function(x){
    (x-min(x,na.rm=T))/(max(x,na.rm = T)-min(x,na.rm=T))
  })
}
#################################split string#######################################################################
linda.split <- function(string,
                        pattern="_",
                        character_pos=1
){
  as.character(lapply(string, function(x) {str_split(x,pattern)[[1]][character_pos]}))
}
#############################################color###########################################################################
ge.color <- function(n){
  if(n<10){
    set.seed(21)
    m <- sample(9,n)
    color <- brewer.pal(9,"Set1")[m]
  }else{
    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))#74种颜色
    set.seed(10)
    color=sample(col_vector, n)
  }
  return(color)
}
##################################################################QC--densityplot######################################################
linda.plot.density <- function(data,outlier=F){
  if(!outlier){
    plot(density(na.omit(unlist(data))),main="density default",mgp = c(2, 0.7, 0),cex.lab=0.9,xlab = "Abundance")      #多个type ggplot?
  }else{
    df2.matrix <- as.matrix(data)
    outline1 <- (fivenum(df2.matrix)[4]-fivenum(df2.matrix)[2])*2+fivenum(df2.matrix)[4]
    outline2 <- fivenum(df2.matrix)[2]-(fivenum(df2.matrix)[4]-fivenum(df2.matrix)[2])*2
    data[data>outline1] <- outline1
    data[data<outline2] <- outline2
    plot(density(na.omit(unlist(data))),main="density default",mgp = c(2, 0.7, 0),cex.lab=0.9,xlab = "Abundance")
  }
}
##########################################################barplot count proteins of each sample#######################################################################################
linda.samplefile.statistic <- function(matrix,cex.axis = 1.5,label =F,col = brewer.pal(12,"Set3")[1]){
  df_sample3=matrix
  count <- apply(df_sample3, 2, function(x) sum(!is.na(x)))
  count_sample <- data.frame(file=names(df_sample3),proteins=count)
  write.csv(count_sample,file = "count_sample.csv",row.names = F)
  data <- count_sample[order(count_sample[,2]),]
  
  # pdf("file_proteins_statistic.pdf",width = 30,height = 10)
  # barplot(data$proteins,
  #         names.arg = data$file,cex.names = 0.8,
  #         xlab = "sample",ylab = "protein number",
  #         cex.axis = cex.axis, 
  #         col = col
  #         # main = "Number of proteins in each sample"
  # )
  # title(main = "Number of proteins in each sample", font.main = 15)
  # dev.off()
  p<- ggbarplot(data, 
                x = "file",
                y = "proteins",
                xlab ="sample",
                width = 0.8,
                ylab = "protein number",
                fill = col,               # 分组调整填充颜色 
                color = "black",            # 轮廓颜色为白色 
                palette =col,
                label = label,   #显示每个bar的数值
                title = "Number of proteins in each sample"           
                # sort.val = "asc",           # 按升序对值进行排序
                # sort.by.groups = T,         # 组内进行排序
                #x.text.angle = 90           # 旋转垂直x轴文本
  )+theme(axis.title.x = element_text(size = 15),axis.title.y = element_text(size = 15))
  # + theme(axis.text.x = element_blank()
  
  ggsave("file.stat.pdf",plot = p,width = 7,height = 5) 
}

linda.samplefile.statbardot <- function(matrix,cex.axis = 1.5,label =NULL,col = c("#377EB8","#984EA3","#A65628"),
                                        dot.size = "missing",
                                        shape = 1,
                                        font.label=list(size = 7, color = "black",vjust=0.5),
                                        legend.position = "none",
                                        plot.width = 10,plot.height = 7
                                        ){
  count <- apply(matrix, 2, function(x) sum(!is.na(x)))
  na_ratio <- apply(matrix, 2, function(x) sum(!is.na(x))/length(x)) 
  count_sample <- data.frame(file=names(matrix),proteins=count,missing=na_ratio)
  data <- count_sample[order(count_sample[,2]),]
  data$batch <- str_extract(data$file,"[Fb][0-9]*")
  write.csv(count_sample,file = "count_sample.csv",row.names = F)
  
  p <- ggdotchart(data,x="file",y="proteins",color = "batch",
                  palette = col,
                  sorting = "desc",#组别内降序排列
                  add = "segment",
                  add.params = list(color = "batch",size=1),#线段的颜色
                  group = "batch",
                  #dot.size = 8,
                  dot.size = dot.size,
                  shape =  shape,
                  label = label,#dot的标签data$proteins
                  #label.rectangle=T,
                  font.label =font.label,#dot的标签颜色大小距离dot的调整距离等
                  xlab = "",ylab = "Protein number",
                  ggtheme = theme_pubr()
                  
  )+
    theme(axis.title.y = element_text(face = 'bold',color = 'black',size = 15),
          #axis.title.x = element_text(face = 'bold',color = 'black',size = 10),
          axis.text.y = element_text(color = 'black',size = 12),
          axis.text.x = element_text(color = 'black',size = 12,angle = 80,vjust = 1), #x轴标签偏转45°，并下降0.1
          # axis.text.x = element_blank(),
          #axis.line = element_line(colour = "black"),
          panel.grid = element_blank(),
          #panel.background = element_blank(),
          #legend.text = element_text(face = 'bold',color = 'black',size = 10)),
          legend.position = legend.position
    )
  ggsave("file.stat_bardotplot.pdf",p,width = plot.width,height = plot.height)
  
}
###########################################################pool correlation corrplot############################################################################
linda.pool_corrplot <- function(df_pool=df_pool,title=paste0("pool","_","correlation.pdf"),
                                r.method = "pearson",#alternative:"kendall", "spearman"
                                plot.method = "circle",#alternative:square', 'ellipse', 'number', 'pie', 'shade' and 'color'
                                plot.type="lower",#c("full", "lower", "upper")
                                color=c(brewer.pal(11,'BrBG')[10:7],"azure1",brewer.pal(11,'BrBG')[5:2]),
                                corr.bool=TRUE,
                                corr.color="grey",
                                tl.col = "black",
                                tl.cex = 0.5,
                                cl.cex = 0.6,
                                cl.ratio = 0.1,
                                number.cex = 1.5
){
  pdf(title,width = 10,height = 10)
  corrplot(cor(df_pool,method = r.method,use = "pairwise.complete.obs"),
           method =  plot.method,
           type = plot.type,
           col = color,
           is.corr = corr.bool,
           addCoef.col = corr.color,
           tl.col = tl.col,
           tl.cex = tl.cex,
           cl.cex = cl.cex,
           cl.ratio = cl.ratio,
           number.cex = number.cex,
           mar=c(2,0,5,2))
  dev.off()
}
###########################################################bioreplivate CV violin############################################################################
linda.violin_biorep.cv <- function(matrix,type,col=col,title=title,width = 8, height = 8){
  data=matrix
  tmp <- data.frame(t(data))
  tmp1 <- tmp[,apply(tmp, 2, function(x) {sum(!is.na(x))>0})]
  
  tmp1 <- data.frame(group=type,tmp1)
  
  tmp.cv <- aggregate(tmp1[,-c(1)],by = list(tmp1$group),FUN = function(x) {sd(x,na.rm = T)/mean(x,na.rm =T)})
  
  box.cv <- melt(tmp.cv,measure.vars = names(tmp.cv)[-1],variable.name = "protein",value.name = "cv")
  
  median <- round(tapply(box.cv$cv,box.cv$Group.1,median,na.rm=T),3)
  
  box.cv$mid <- median[match(box.cv$Group.1,names(median))]
  
  p <- ggplot(box.cv,aes(x = Group.1, y=cv,color= Group.1)) + 
    geom_violin(trim=F)+
    ylab("Coefficient of variation")+
    scale_fill_manual(values=col)+
    geom_boxplot(width=0.2)+
    xlab("")+
    # stat_compare_means(method = "wilcox.test",comparisons = my_comparisons)+
    # stat_compare_means(label.y = 1.5,label.x = 2.2)+
    #coord_trans(x = "identity", y = "identity", xlim = NULL, ylim = c(0,0.8))+
    geom_hline(yintercept = 0.1,linetype=2,show.legend = T)+
    theme(legend.direction = 'horizontal',legend.position = 'top',panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))+
    theme(axis.text = element_text(size = 15,colour = "black"),text = element_text(size = 15,colour = "black"))+
    theme(axis.text.x = element_text(hjust = 1,angle = 45))+
    theme(legend.position = "none")+ #reomove legend
    theme(plot.margin = unit(rep(2,4),'lines'))+ #调整边缘距离
    geom_text(aes(y=0.6,label=paste0("Median=",box.cv$mid)),hjust=0.5,vjust=1,size=5) #调整文字位置
  ggsave(paste0(title,"_violin_CV.pdf"),plot =p ,device = NULL,width = width, height = height) 
}
###########################################################bioreplivate CV half_violin###########################################################
linda.halfviolin_biorep.cv <-  function(matrix,#待计算的全距阵
                                        info,#File.Name列（matrix的列名）和Biology_replicate列(sample name名）是必须的
                                        dup,#生物重复的sample名，非重复非实验名
                                        col,title,width = 8, height = 8){
  
  bio.cv <- c()
  Bio_CV_Median <- c()
  violin.sum <- c()
  for (i in dup) {
    tmp <-
      matrix[, which(names(matrix) %in% info$File.Name[which(info$Biology_replicate== i)])]
    tmp.cv <- apply(tmp, 1 , function(x) {
      # sd(x) / mean(x)
      sd(x,na.rm = T) / mean(x,na.rm=T)
    })
    bio.cv <- c(bio.cv,tmp.cv)
    Bio_CV_Median <- c(Bio_CV_Median,median(bio.cv,na.rm = T))
  }
  violin.sum <- rbind(violin.sum,data.frame(Group=rep(dup,each=nrow(matrix)),cv=as.numeric(bio.cv),median=rep(Bio_CV_Median,each=nrow(matrix))))
  
  stat_data <-  violin.sum[,-2] %>% unique()
  p <- ggplot(violin.sum,aes(x = Group, y=cv,fill= Group)) + 
    geom_half_violin(position = position_nudge(x = 0.2),side=2,alpha = 0.8)+
    scale_fill_manual(values=col)+
    geom_jitter(aes(y=cv, color = Group), 
                position = position_jitter(width = 0.12),#鎵�鏈夌偣鍗犵殑瀹藉害
                size = 0.2,#鐐圭殑澶у皬
                alpha = 0.5)+
    scale_color_manual(values=col)+
    geom_boxplot(width=0.05,outlier.shape = NA,fill="white",size=0.5,position=position_nudge(x=0.2))+ 
    geom_text(data=stat_data,aes(y=median-0.02,x=Group,label=paste0("Median=",round(median,3)),
                                 vjust = 0.5, #绔栫洿涓婁笅璋?
                                 hjust = 0.1#姘村钩宸﹀彸璋冩暣
    ),size=5)+
    ylab("Coefficient of variation")+
    xlab("")+
    coord_flip()+
    theme_bw()+ #鍔犲洓杈规
    theme(panel.grid =element_blank(), 
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))+
    theme(axis.text = element_text(size = 12,colour = "black"),
          text = element_text(size = 10,colour = "black"),#璋冩暣label text鐨勫瓧浣撳ぇ灏?
          axis.title.x = element_text(size = 15,colour = "black")
    )+
    theme(legend.position = "none")+ #reomove legend
    theme(plot.margin = unit(rep(2,4),'lines')) #璋冩暣杈圭紭璺濈
  ggsave(paste0(title,"_violin_CV.pdf"),plot =p ,device = NULL,width = width, height = height) 
}
######################################################correlation smooth scatter density plot#########################################
library(ggpointdensity)
library(viridis) #色系
linda.techcor.scatter <- function(corrdata,method.cor="pearson",point.viridiscolor="D",smooth.line="lm",smooth.col="#e13c50",text.size=7,title="",xlab="rep1",ylab="rep2"){
  cor1 = corrdata[,1]
  cor2 = corrdata[,2]
  r=cor.test(cor1,cor2,method = method.cor)
  p <- ggplot(corrdata,aes(get(names(corrdata)[1]),get(names(corrdata)[2])))+
    geom_pointdensity(adjust=0.2,
                      #size = 3
                      #shape = 17
    )+
    scale_color_viridis(option = point.viridiscolor)+
    geom_smooth(method = smooth.line,col=smooth.col,lty=2,bg="white")+
    geom_text(aes(x=min(cor1,na.rm = T)+0.5,y=max(cor2,na.rm = T)-0.5,label=paste0("r: ",round(r$estimate,2))),size = text.size)+
    geom_text(aes(x=min(cor1,na.rm = T)+0.5,y=max(cor2,na.rm = T)-0.8,label=paste0("pvalue: ",round(r$p.value,2))),size = text.size)+
    theme_base()+
    theme(legend.position = "none")+
    xlab(xlab)+
    ylab(ylab)
  ggsave(paste0(title,"Correlation.scatter_density.pdf"),plot = p)
}
###########################################################pca############################################################################
ge.plot.pca <- function(data,type,title="",ptColors=NULL,label2=NULL,width=12,height=8){ 
  M <- t(data)
  M[is.na(M)] <- 0
  M <- apply(M,2,function(v){(v-mean(v,na.rm=T))/sd(v,na.rm=T)})
  clnames <- row.names(data)
  set.seed(200)
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
  
  ggsave(paste0(title,"_PCA.pdf"),plot =p ,width=width,height=height,device = NULL)
}
################################################################umap######################################################################
library(umap)
ge.plot.umap<- function(data,type,label1,batch_color=batch_color,title="",label=NULL){
  # col2=ge.color(length(unique(type)))
  col2= batch_color[c(1:length(unique(type)))]
  cl <- data.frame(col2,row.names =unique(type),stringsAsFactors = F)
  cl2 <- cl[match(type,row.names(cl)),1]
  
  df10 <- data
  df10[is.na(df10)] <- 0
  df10 <- t(apply(df10, 1, scale))
  colnames(df10) <- label1
  set.seed(10)
  df.umap <- umap(t(df10),n_neighbors=ncol(data)-1)
  
  umap_data <- df.umap$layout
  colnames(umap_data) <- c("umap1","umap2")
  pdf(paste0(title,"_UMAP.pdf"))
  if(is.null(label)){
    # plot(df.umap$layout,col = cl2, main = "umap", pch = 20,cex=2,cex.axis=2,cex.lab=2)
    plot(umap_data,col = cl2, main = "umap", pch = 20,cex=2,cex.axis=2,cex.lab=2)
    
    legend("topright",legend=row.names(cl), fill=cl$col2, lty=1,lwd=1)
  }else{
    # plot(df.umap$layout,col = cl2, main = "umap", pch = 20,cex=2,cex.axis=2,cex.lab=2)
    plot(umap_data,col = cl2, main = "umap", pch = 20,cex=2,cex.axis=2,cex.lab=2)
    
    text(umap_data, pos =2, labels = label, col= "DimGrey",cex = 0.4)
    legend("topright",legend=row.names(cl), fill=cl$col2, lty=1,lwd=1)
  }
  dev.off()
  return(df.umap)
}
# #################################################################计算T-test & foldchange##################################################
linda.Ttest <- function(matrix,
                        group1,#postion number of group1 samples in clumns
                        group2,#postion number of group2 samples in clumns
                        paired_bool = F,
                        var.equal_bool=F,
                        fold_change=1.5,
                        pvalue=0.05,
                        title=title){
  df1 <- as.data.frame(matrix)
  df1$`log2(foldchange)` <- as.numeric(apply(df1,1, function(x) log2(abs((mean(x[group1],na.rm = T)/mean(x[group2],na.rm = T))))))
  ##可能会有负数
  P_value <- apply(df1,
                   1,
                   function(y) {
                     
                     p_try = tryCatch(t.test(y[group1],
                                             y[group2],
                                             paired = paired_bool,
                                             var.equal = var.equal_bool)$p.value,
                                      error = function(x) NA)
                   })  
  
  #给data增加p-value列
  df1$P_value<- P_value
  df1$P_value_adjust<-p.adjust(df1$P_value, method="BH")
  up <- subset(df1, df1$P_value < pvalue & df1$`log2(foldchange)` > log2(fold_change))
  down <- subset(df1, df1$P_value < pvalue & df1$`log2(foldchange)` < -log2(fold_change))
  print(paste0("down:",nrow(down)))
  print(paste0("up:",nrow(up)))
  write.csv(up,file = paste0(title,"_upPvalue0.05_volcano.csv"))
  write.csv(down,file = paste0(title,"_dwPvalue0.05_volcano.csv"))
  write.csv(df1,file = paste0(title, "_all_volcano.csv"))
}
# #################################################################不填充NA的t-TEST多################################################################################
linda.volcano.plot <- function(data,#data.frame
                               logtransformation=T,
                         group1,
                         group2,
                         label1,
                         label2,
                         paired = F,
                         pvalue= 0.05,
                         P.adjust=T,
                         var.equal = F,
                         fold_change = 2,
                         up_col="#EC9B8F", #dw_col = "#9370DB",up_col = "#FF7256"
                         dw_col="#A7DEE8"){
  df1 <- as.data.frame(data)
  #填充NA
  #df1[is.na(df1)]=min(df1,na.rm = T)*0.8
  #计算foldchange
  if(logtransformation){
    df1$`log2(foldchange)` <- apply(2^df1,1, function(x) log2((mean(x[group1],na.rm = T)/mean(x[group2],na.rm = T))))
  }else{
    df1$`log2(foldchange)` <- apply(df1,1, function(x) log2((mean(x[group1],na.rm = T)/mean(x[group2],na.rm = T))))
    }
  
  if(paired){
    P_value <- apply(df1,
                     1,
                     function(y) {
                       
                       p_try = tryCatch(t.test(y[group1],
                                               y[group2],
                                               paired = paired,
                                               var.equal = var.equal)$p.value,
                                        error = function(x) NA)
                     })
  }else{
    P_value <- apply(df1,
                     1,
                     function(y) {
                       
                       p_try = tryCatch(t.test(y[group1],
                                               y[group2],
                                               paired = F,
                                               var.equal = var.equal)$p.value,
                                        error = function(x) NA)
                     })  
  }
  #给data增加p-value列
  df1$P_value<- P_value
  df1$P_value_adjust<-p.adjust(df1$P_value, method="BH")
  
  if(P.adjust){
    #存储上调和下调蛋白信息
    up <- subset(df1, df1$P_value_adjust < pvalue & df1$`log2(foldchange)` > log2(fold_change))
    down <- subset(df1, df1$P_value_adjust < pvalue & df1$`log2(foldchange)` < -log2(fold_change))
    print(paste0("down:",nrow(down)))
    print(paste0("up:",nrow(up)))
    write.csv(up,file = paste0(label1,"_",label2, "_up_volcano.csv"))
    write.csv(down,file = paste0(label1,"_",label2, "_dw_volcano.csv"))
    write.csv(df1,file = paste0(label1,"_",label2, "_all_volcano.csv"))
    differprot <- c(row.names(up),row.names(down))
    
    #画图
    pdf(paste0(label1,"_",label2, "_volcano.pdf"),width=4, height=4)
    
    plot(df1$`log2(foldchange)`, -log10(df1$P_value_adjust), col="#00000033", pch=19,
         xlab=paste("log2 (fold change)"),
         ylab="-log10 (adjust P value)",
         xlim= c(-round(max(abs(df1$`log2(foldchange)`),na.rm=T),1),round(max(abs(df1$`log2(foldchange)`),na.rm=T),1)),#na.rm=T只能过滤NA值如果有-lnf等也不可以计算最值
         main=paste0(label1," : ",label2))
    points(up$`log2(foldchange)`, -log10(up$P_value_adjust), col=1, 
           #bg = brewer.pal(9, "YlOrRd")[6], 
           bg=up_col,
           pch=21, cex=1)
    points(down$`log2(foldchange)`, -log10(down$P_value_adjust), col = 1, 
           #bg = brewer.pal(11,"RdBu")[9], 
           bg=dw_col,
           pch = 21,cex=1)
    abline(h=-log10(pvalue), v=c(-log2(fold_change),log2(fold_change)), lty=2,lwd=1)
    dev.off()}else{
      #存储上调和下调蛋白信息
      up <- subset(df1, df1$P_value < pvalue & df1$`log2(foldchange)` > log2(fold_change))
      down <- subset(df1, df1$P_value < pvalue & df1$`log2(foldchange)` < -log2(fold_change))
      print(paste0("down:",nrow(down)))
      print(paste0("up:",nrow(up)))
      write.csv(up,file = paste0(label1,"_",label2, "_up_volcano.csv"))
      write.csv(down,file = paste0(label1,"_",label2, "_dw_volcano.csv"))
      write.csv(df1,file = paste0(label1,"_",label2, "_all_volcano.csv"))
      differprot <- c(row.names(up),row.names(down))
      #画图
      pdf(paste0(label1,"_",label2, "_volcano.pdf"),width=4, height=4)
      
      plot(df1$`log2(foldchange)`, -log10(df1$P_value), col="#00000033", pch=19,
           xlab=paste("log2 (fold change)"),
           ylab="-log10 (P value)",xlim=c(-round(max(abs(df1$`log2(foldchange)`),na.rm=T),1),round(max(abs(df1$`log2(foldchange)`),na.rm=T),1)),
           main=paste0(label1," : ",label2))

      points(up$`log2(foldchange)`, -log10(up$P_value), col=1, 
            # bg = brewer.pal(9, "YlOrRd")[6],
            bg=up_col,
             pch=21, cex=1.2)
      points(down$`log2(foldchange)`, -log10(down$P_value), col = 1, 
            # bg = brewer.pal(11,"RdBu")[9],
             bg=dw_col,
             pch = 21,cex=1.2)
      abline(h=-log10(pvalue), v=c(-log2(fold_change),log2(fold_change)), lty=2,lwd=1)
      dev.off()
    }
  return(differprot)
}
#############################################label volcano##########################################

linda.volcano.plot2 <- function(data,#data.frame
                                logtransformation=F,
                                group1,
                                group2,
                                label1,
                                label2,
                                paired = F,
                                pvalue= 0.05,
                                P.adjust=T,
                                var.equal = F,
                                fold_change = 2){
  df1 <- as.data.frame(data)
  
  if(logtransformation){
    df1$`log2(foldchange)` <- apply(2^df1,1, function(x) log2((mean(x[group1],na.rm = T)/mean(x[group2],na.rm = T))))
  }else{
    df1$`log2(foldchange)` <- apply(df1,1, function(x) log2((mean(x[group1],na.rm = T)/mean(x[group2],na.rm = T))))
  }
  
  if(paired){
    P_value <- apply(df1,
                     1,
                     function(y) {
                       
                       p_try = tryCatch(t.test(y[group1],
                                               y[group2],
                                               paired = paired,
                                               var.equal = var.equal)$p.value,
                                        error = function(x) NA)
                     })
  }else{
    P_value <- apply(df1,
                     1,
                     function(y) {
                       
                       p_try = tryCatch(t.test(y[group1],
                                               y[group2],
                                               paired = F,
                                               var.equal = var.equal)$p.value,
                                        error = function(x) NA)
                     })  
  }
  #给data增加p-value列
  df1$P_value<- P_value
  df1$P_value_adjust<-p.adjust(df1$P_value, method="BH")
  
  if(P.adjust){
    
    #存储上调和下调蛋白信息
    up <- subset(df1, df1$P_value_adjust < pvalue & df1$`log2(foldchange)` > log2(fold_change))
    down <- subset(df1, df1$P_value_adjust < pvalue & df1$`log2(foldchange)` < -log2(fold_change))
    print(paste0("down:",nrow(down)))
    print(paste0("up:",nrow(up)))
    write.csv(up,file = paste0(label1,"_",label2, "_up_volcano.csv"))
    write.csv(down,file = paste0(label1,"_",label2, "_dw_volcano.csv"))
    write.csv(df1,file = paste0(label1,"_",label2, "_all_volcano.csv"))
    differprot <- c(row.names(up),row.names(down))
    
    #画图
    pdf(paste0(label1,"_",label2, "_volcano.pdf"),width=4, height=4)
    
    plot(df1$`log2(foldchange)`, -log10(df1$P_value_adjust), col="#00000033", pch=19,cex=0.5,
         xlab=paste("log2 (fold change)"),
         ylab="-log10 (adjust P value)",
         xlim= c(-round(max(abs(df1$`log2(foldchange)`),na.rm=T),1),#na.rm=T只能过滤NA值如果log2(foldchange)有非法字符如-lnf等也不可以计算最值
                 round(max(abs(df1$`log2(foldchange)`),na.rm=T),1)),
         main=paste0(label1," : ",label2))
  
    points(up$`log2(foldchange)`, -log10(up$P_value_adjust), col="#FFA500", bg = "#FFA500", pch=21,cex=0.6)
    points(down$`log2(foldchange)`, -log10(down$P_value_adjust), col = "#32CD32", bg = "#32CD32", pch = 21,cex=0.6)
    abline(h=-log10(pvalue), v=c(-log2(fold_change),log2(fold_change)), lty=2,lwd=1)
    if(length(row.names(down))!=0 & length(row.names(up))!=0){
      text(down$`log2(foldchange)`, -log10(down$P_value_adjust),labels= row.names(down),adj = c(1,1),cex = 0.2)
      text(up$`log2(foldchange)`, -log10(up$P_value_adjust),labels= row.names(up),adj = c(1,1),cex = 0.2)
    }else if(length(row.names(down))==0 & length(row.names(up))!=0){
      text(up$`log2(foldchange)`, -log10(up$P_value_adjust),labels= row.names(up),adj = c(1,1),cex = 0.2)
    }else if(length(row.names(down))!=0 & length(row.names(up))==0){
      text(down$`log2(foldchange)`, -log10(down$P_value_adjust),labels= row.names(down),adj = c(1,1),cex = 0.2)
    }
    dev.off()}else{
      #存储上调和下调蛋白信息
      up <- subset(df1, df1$P_value < pvalue & df1$`log2(foldchange)` > log2(fold_change))
      down <- subset(df1, df1$P_value < pvalue & df1$`log2(foldchange)` < -log2(fold_change))
      print(paste0("down:",nrow(down)))
      print(paste0("up:",nrow(up)))
      write.csv(up,file = paste0(label1,"_",label2, "_up_volcano.csv"))
      write.csv(down,file = paste0(label1,"_",label2, "_dw_volcano.csv"))
      write.csv(df1,file = paste0(label1,"_",label2, "_all_volcano.csv"))
      differprot <- c(row.names(up),row.names(down))
      #画图
      pdf(paste0(label1,"_",label2, "_volcano.pdf"),width=4, height=4)
      
      plot(df1$`log2(foldchange)`, -log10(df1$P_value), col="#00000033", pch=19,cex=0.5,
           xlab=paste("log2 (fold change)"),
           ylab="-log10 (P value)",xlim=c(-round(max(abs(df1$`log2(foldchange)`),na.rm=T),1),round(max(abs(df1$`log2(foldchange)`),na.rm=T),1)),
           main=paste0(label1," : ",label2))
      points(up$`log2(foldchange)`, -log10(up$P_value), col="#FFA500", bg = "#FFA500", pch=21, cex=0.6)
      points(down$`log2(foldchange)`, -log10(down$P_value), col = "#32CD32", bg = "#32CD32" , pch = 21,cex=0.6)
      abline(h=-log10(pvalue), v=c(-log2(fold_change),log2(fold_change)), lty=2,lwd=1)
      if(length(row.names(down))!=0 & length(row.names(up))!=0){
        text(down$`log2(foldchange)`, -log10(down$P_value),labels= row.names(down),adj = c(1,1),cex = 0.2)
        text(up$`log2(foldchange)`, -log10(up$P_value),labels= row.names(up),adj = c(1,1),cex = 0.2)
      }else if(length(row.names(down))==0 & length(row.names(up))!=0){
        text(up$`log2(foldchange)`, -log10(up$P_value),labels= row.names(up),adj = c(1,1),cex = 0.2)
      }else if(length(row.names(down))!=0 & length(row.names(up))==0){
        text(down$`log2(foldchange)`, -log10(down$P_value),labels= row.names(down),adj = c(1,1),cex = 0.2)
      }
      dev.off()
    }
  return(differprot)
}
##############################renew: solve the problem of xlim infinite plot by ggpubr and ggplot2#################################################
linda.volcano.plot3 <- function(data,#data.frame
                                logtransformation=F,
                                group1,
                                group2,
                                label1,
                                label2,
                                paired = F,
                                pvalue= 0.05,
                                P.adjust=T,
                                var.equal = F,
                                fold_change = 2,
                                select_gene = NULL,
                                palette =  c("#ff7e00","#00beff","#00000033") #修改颜色
){
  data <- data[,c(group1,group2)]
  data <- data[apply(data, 1, function(x){sum(!is.na(x))})>0,]
  
  df1 <- as.data.frame(data)
  
  if(logtransformation){
    df1$`log2(foldchange)` <- apply(2^df1,1, function(x) log2((mean(x[group1],na.rm = T)/mean(x[group2],na.rm = T))))
  }else{
    df1$`log2(foldchange)` <- apply(df1,1, function(x) log2((mean(x[group1],na.rm = T)/mean(x[group2],na.rm = T))))
  }
  
  if(paired){
    P_value <- apply(df1,
                     1,
                     function(y) {
                       
                       p_try = tryCatch(t.test(y[group1],
                                               y[group2],
                                               paired = paired,
                                               var.equal = var.equal)$p.value,
                                        error = function(x) NA)
                     })
  }else{
    P_value <- apply(df1,
                     1,
                     function(y) {
                       
                       p_try = tryCatch(t.test(y[group1],
                                               y[group2],
                                               paired = F,
                                               var.equal = var.equal)$p.value,
                                        error = function(x) NA)
                     })  
  }
  #给data增加p-value列
  df1$P_value<- P_value
  df1$P_value_adjust<- p.adjust(df1$P_value, method="BH")
  
  if(P.adjust){
    df1$P <- df1$P_value_adjust
    y.lab <- "-log10 (adjust P_value)"
  }else{
    df1$P <- df1$P_value
    y.lab <- "-log10 (unadjust P_value)"
  }
  
  df1$Class = factor(ifelse(df1$P < pvalue & (df1$`log2(foldchange)` > log2(fold_change) | df1$`log2(foldchange)` < -log2(fold_change)),ifelse(df1$`log2(foldchange)` > log2(fold_change),'Up','Down'),'NoSignifi'),levels=c('Up','Down','NoSignifi'))
  df1$pro <-  rownames(df1)
  up <- subset(df1, df1$P < pvalue & df1$`log2(foldchange)` > log2(fold_change))
  down <- subset(df1, df1$P < pvalue & df1$`log2(foldchange)` < -log2(fold_change))
  print(paste0("down:",nrow(down)))
  print(paste0("up:",nrow(up)))
  write.csv(up,file = paste0(label1,"_",label2, "_up_volcano.csv"))
  write.csv(down,file = paste0(label1,"_",label2, "_dw_volcano.csv"))
  write.csv(df1,file = paste0(label1,"_",label2, "_all_volcano.csv"))
  differprot <- c(row.names(up),row.names(down))
  
  df1$`-log10(P)` <- -log10(df1$P)
  
  # if(!is.null(select_gene)) {
  #   df1$sign <- ifelse(rownames(df1)%in%select_gene,rownames(df1)[rownames(df1)%in%select_gene],NA)
  # }
  # else {df1$sign <- NA}
  df1$sign <- NA
  df1$sign[which(rownames(df1)%in%as.character(select_gene))] <- as.character(select_gene)
  
  p <- ggscatter(df1,x= "log2(foldchange)",y="-log10(P)",
                 
                 xlab = "log2 (foldchange)",
                 ylab = y.lab,
                 color="Class",
                 palette = palette, #修改颜色
                 label = "sign",repel = T, #显示每个点表示的样本名
                 #ellipse = F, #分组圈图
                 font.label = c(10, "bold","black"),#label的样式,如果不指定颜色则默认用点的颜色代替
                 shape = 19,
                 #size = 1.5,
                 size = 2.5,
                 show.legend.text = F,
                 # title="DLE vs. NC" #不知道怎么调大小和位置
                 title=paste(label1,"vs.",label2)
  )+
    geom_vline(xintercept = c(-log2(fold_change),log2(fold_change)),linetype=2)+
    geom_hline(yintercept = -log10(pvalue),linetype=2)+
    geom_text(aes(x=-max(df1$`log2(foldchange)`[is.finite(df1$`log2(foldchange)`)])/2,y=2,label=paste0("down: ",sum(na.omit(df1$`log2(foldchange)`< (-log2(fold_change)) & df1$`-log10(P)`> (-log10(pvalue)))))),size = 5)+
    geom_text(aes(x=max(df1$`log2(foldchange)`[is.finite(df1$`log2(foldchange)`)])/2,y=2,label=paste0("up: ",sum(na.omit(df1$`log2(foldchange)`> (log2(fold_change)) & df1$`-log10(P)`> (-log10(pvalue)))))),size = 5)+
    theme(axis.title.x =element_text(size=15), axis.title.y=element_text(size=15),
          axis.text = element_text(colour ="black", size = 10),
          axis.line.x = element_line(size = 0.5),
          axis.line.y = element_line(size = 0.5),
          legend.position = "none",
          panel.background = element_blank(),
          panel.grid =element_blank(),
          axis.ticks = element_blank())
  # theme_base()#调整legend位置在图的右侧
  ggsave(paste0(label1,"_",label2, "_volcano_id.pdf"),plot =p,device = NULL)
  
  return(differprot)
}
################label
# library(ggrepel)
tan.volcano.plot <- function(data,#data.frame
                             group1,
                             group2,
                             label1,
                             label2,
                             paired = F,
                             pvalue= 0.05,
                             P.adjust= F,
                             var.equal = F,
                             select_gene = NULL,
                             fold_change1 = 1.2,
                             fold_change2 = 1.2,
                             label_size = 3){
  df1 <- data
  #濉厖NA
  df1[is.na(df1)]=min(df1,na.rm = T)
  #璁＄畻foldchange
  df1$`log2(foldchange)` <- apply(df1,1, function(x) log2((mean(x[group1],na.rm = T)/mean(x[group2],na.rm = T))))
  
  if(paired){
    P_value <- apply(df1,
                     1,
                     function(y) {
                       
                       p_try = tryCatch(t.test(y[group1],
                                               y[group2],
                                               paired = paired,
                                               var.equal = var.equal)$p.value,
                                        error = function(x) NA)
                     })
  }else{
    P_value <- apply(df1,
                     1,
                     function(y) {
                       
                       p_try = tryCatch(t.test(y[group1],
                                               y[group2],
                                               paired = F,
                                               var.equal = var.equal)$p.value,
                                        error = function(x) NA)
                     })  
  }
  #缁檇ata澧炲姞p-value鍒?
  df1$P_value<- P_value
  df1$P_value_adjust<-p.adjust(df1$P_value, method="BH")
  
  
  if(P.adjust){
    df1$threshold = factor(ifelse(df1$P_value_adjust < pvalue & (df1$`log2(foldchange)` > log2(fold_change1) | df1$`log2(foldchange)` < -log2(fold_change2)),ifelse(df1$`log2(foldchange)` > log2(fold_change1),'Up','Down'),'NoSignifi'),levels=c('Up','Down','NoSignifi'))
    
    df1$pro <-  rownames(df1)
    up <- subset(df1, df1$P_value_adjust < pvalue & df1$`log2(foldchange)` > log2(fold_change1))
    down <- subset(df1, df1$P_value_adjust < pvalue & df1$`log2(foldchange)` < -log2(fold_change2))
    print(paste0("down:",nrow(down)))
    print(paste0("up:",nrow(up)))
    write.csv(up,file = paste0(label1,"_",label2, "_up_volcano.csv"))
    write.csv(down,file = paste0(label1,"_",label2, "_dw_volcano.csv"))
    write.csv(df1,file = paste0(label1,"_",label2, "_all_volcano.csv"))
    differprot <- list(up=row.names(up),down=row.names(down))
    
    if(!is.null(select_gene)) {data_text = df1[select_gene,]}
    else {data_text = subset(df1, P_value_adjust < pvalue & (df1$`log2(foldchange)` > log2(fold_change1) | df1$`log2(foldchange)` < -log2(fold_change2)))}
    
    pp <- ggplot(df1,aes(x=`log2(foldchange)`,y=-log10(P_value_adjust),color=threshold))+
      geom_point()+
      scale_color_manual(values=c('Up'="#DC143C",'Down'="#0A9731",'NoSignifi'="#808080"))+#纭畾鐐圭殑棰滆壊
      geom_text_repel(
        # data = subset(df1, P_value_adjust < pvalue & (df1$`log2(foldchange)` > log2(fold_change1) | df1$`log2(foldchange)` < -log2(fold_change2))),
        data = data_text,
        aes(label = pro),
        colour = "black",  ###鏍囨敞鍑烘潵鐨勫熀鍥犵殑棰滆壊
        size = label_size,
        max.iter=2500,
        segment.color = "black", show.legend = FALSE )+
      theme_bw()+#淇敼鍥剧墖鑳屾櫙
      theme(
        legend.title = element_blank()#涓嶆樉绀哄浘渚嬫爣棰?
      )+
      ylab("-log10 (adjust P value)")+#淇敼y杞村悕绉?
      xlab(paste("log2 (fold change)"))+#淇敼x杞村悕绉?
      geom_vline(xintercept=c(-log2(fold_change2),log2(fold_change1)),lty=3,col="black",lwd=0.5) +#娣诲姞妯嚎|FoldChange|>2
      geom_hline(yintercept = -log10(pvalue),lty=3,col="black",lwd=0.5)#娣诲姞绔栫嚎padj<0.05
    #dev.off()
    ggsave(paste0(label1,"_",label2, "_volcano_id.pdf"),plot =pp ,width=8,height=6,device = NULL)
  }else{
    df1$threshold = factor(ifelse(df1$P_value < pvalue & (df1$`log2(foldchange)` > log2(fold_change1) | df1$`log2(foldchange)` < -log2(fold_change2)),ifelse(df1$`log2(foldchange)` > log2(fold_change1),'Up','Down'),'NoSignifi'),levels=c('Up','Down','NoSignifi'))
    df1$pro <-  rownames(df1)
    #鐢诲浘
    if(!is.null(select_gene)) {data_text = df1[select_gene,]}
    else {data_text = subset(df1, P_value_adjust < pvalue & (df1$`log2(foldchange)` > log2(fold_change1) | df1$`log2(foldchange)` < -log2(fold_change2)))}
    
    pp <- ggplot(df1,aes(x=`log2(foldchange)`,y=-log10(P_value),color=threshold))+
      
      geom_point()+
      scale_color_manual(values=c('Up'="#DC143C",'Down'="#0A9731",'NoSignifi'="#808080"))+#纭畾鐐圭殑棰滆壊
      geom_text_repel(
        data = data_text,
        # data = subset(df1, P_value < pvalue & (df1$`log2(foldchange)` > log2(fold_change1) | df1$`log2(foldchange)` < -log2(fold_change2))),
        aes(label = pro),
        colour = "black", ###鏍囨敞鍑烘潵鐨勫熀鍥犵殑棰滆壊
        size = label_size,
        max.iter=3000,
        segment.color = "black", show.legend = FALSE )+
      
      theme_bw()+#淇敼鍥剧墖鑳屾櫙
      theme(
        legend.title = element_blank()#涓嶆樉绀哄浘渚嬫爣棰?
      )+
      ylab("-log10 (P value)")+#淇敼y杞村悕绉?
      xlab(paste("log2 (fold change)"))+#淇敼x杞村悕绉?
      geom_vline(xintercept=c(-log2(fold_change2),log2(fold_change1)),lty=3,col="black",lwd=0.5) +#娣诲姞妯嚎|FoldChange|>2
      geom_hline(yintercept = -log10(pvalue),lty=3,col="black",lwd=0.5)#娣诲姞绔栫嚎padj<0.05
    ggsave(paste0(label1,"_",label2, "_volcano_id.pdf"),plot =pp ,width=12,height=8,device = NULL)
    
    #瀛樺偍涓婅皟鍜屼笅璋冭泲鐧戒俊鎭?
    up <- subset(df1, df1$P_value < pvalue & df1$`log2(foldchange)` > log2(fold_change1))
    down <- subset(df1, df1$P_value < pvalue & df1$`log2(foldchange)` < -log2(fold_change2))
    print(paste0("down:",nrow(down)))
    print(paste0("up:",nrow(up)))
    write.csv(up,file = paste0(label1,"_",label2, "_up_volcano.csv"))
    write.csv(down,file = paste0(label1,"_",label2, "_dw_volcano.csv"))
    write.csv(df1,file = paste0(label1,"_",label2, "_all_volcano.csv"))
    differprot <- list(up=row.names(up),down=row.names(down))
    
  }
  return(differprot)
}
##################################################################填充NA的t-TEST多################################################################################

########################################annova###########################
linda.annova2 <- function(data,label,adjust.p=T,pvalue=0.05,title=""){
  df6 =data
  df.allanov <- data.frame(deal=label,t(df6))
  row.names(df.allanov) <-names(df6)
  p <- c()
  for(i in 2:ncol(df.allanov)){
    name <- df.allanov[,colnames(df.allanov)[i]]
    anova_p = tryCatch(summary(aov(name ~ df.allanov$deal))[[1]]$`Pr(>F)`[1] %>%format(digits = 3, scientific = FALSE),error = function(x) NA)#保证组间自由度大于等于1，至少有2组label有数据否则会error
    p <- append(p,anova_p)
  }
  
  adp <- p.adjust(p, method = "BH")%>%format(digits = 3, scientific = FALSE)
  
  allanov_result <- data.frame(df6,pvalue=p,adjust.pvalue = adp)
  write.csv(allanov_result,file = paste0(title,"_allprot_ajustp0.05.csv"),row.names = T)
  if(adjust.p){
    differ_result <- allanov_result[(allanov_result$adjust.pvalue<0.05),] 
    write.csv(differ_result,file = paste0(title,"_differprot_ajustp0.05.csv"),row.names = T)
  }else{
    differ_result <- allanov_result[(allanov_result$pvalue<0.05),] 
    write.csv(differ_result,file = paste0(title,"_differprot_pvalue0.05.csv"),row.names = T)
  }
}
###################################mfuzz############################
ge.mfuzz.cselection <- function(data,range=seq(5,50,5),repeats = 5){
  df3a<-as.matrix(data)
  df3Ex<- ExpressionSet(assayData = df3a)
  if(interactive()){
    df3F <- filter.NA(df3Ex)
    df3F <- fill.NA(df3F)
    df3F <- standardise(df3F) #标准化分析
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
    write.csv(cl[[4]][potname,i],paste0(filename,"/mfuzz_",i,".csv"))  #存储每个cluster的蛋白名
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
      write.csv(newdf,paste0(filename,"/mfuzz_anova_",ii,".csv"),row.names = F)  #给每个cluster在time label上加上一个annova分析，并存储结果
    }
  }
}
#############################################示例：linda.tsne(matrix,differ_grade,"Gradenew")################################################################
linda.tsne <- function(data,type,shape=24,label=F,ellipse=F,col=c("#33709F","#CC6C5E"),size=5,width=10,height=7,title=""){
  df10 <- data
  df10[is.na(df10)] <- 0
  df10 <- t(apply(df10, 1, scale))
  df10[is.na(df10)] <- 0
  df10<-data.frame(df10)
  colnames(df10) <- type
  #set.seed(10)
  df11.tsne <- Rtsne(t(df10), 
                     dims = 2,
                     perplexity = (ncol(data)-1)/3-1, 
                     verbose = T , 
                     check_duplicates = FALSE)
  #
  colnames(df11.tsne$Y) <- c("TSNE1","TSNE2")
  tsne.data <- data.frame(sample=colnames(data),Type=type,df11.tsne$Y,samplename=names(data))
  
  if(label){
  p <- ggscatter(tsne.data,x= "TSNE1",y="TSNE2",
                 color="Type", 
                 palette = col,
                 # palette = c("blue","red"), 修改颜色
                 #label = "samplename",repel = T, #显示每个点表示的样本名
                 label = "samplename",repel = T,
                 ellipse = ellipse, #分组圈图
                 shape = shape, #形状为三角形24,圆是19
                 fill = "Type", #三角形的填充颜色
                 #fill = c("#33709F","#CC6C5E"),
                 size = size,
                 main=paste0(nrow(df10),"features"))+
    theme(axis.title.x =element_text(size=15), axis.title.y=element_text(size=15),
          axis.text = element_text(colour ="black", size = 10),
          axis.line.x = element_line(size = 0.5),
          axis.line.y = element_line(size = 0.5),
          panel.background = element_blank(),
          panel.grid =element_blank(),
          axis.ticks = element_blank())+
    theme_base()#调整图里位置在图的右侧，没有theme_base() 则在上方且图形没有闭合的黑框（上和右）
  }else{
    p <- ggscatter(tsne.data,x= "TSNE1",y="TSNE2",
                   color="Type", 
                   palette = col,
                   # palette = c("blue","red"), 修改颜色
                   #label = "samplename",repel = T, #显示每个点表示的样本名
                   ellipse = F, #分组圈图
                   shape = shape, #形状为三角形24,圆是19
                   fill = "Type", #三角形的填充颜色
                   #fill = c("#33709F","#CC6C5E"),
                   size = size,
                   main="")+
      theme(axis.title.x =element_text(size=15), axis.title.y=element_text(size=15),
            axis.text = element_text(colour ="black", size = 10),
            axis.line.x = element_line(size = 0.5),
            axis.line.y = element_line(size = 0.5),
            panel.background = element_blank(),
            panel.grid =element_blank(),
            axis.ticks = element_blank())+
      theme_base()
  }
  ggsave(paste0(title,"_tsne.pdf"),plot = p,width=width,height=height,device = NULL)
} 
#Rtsne参数说明:
### dims降维后的维度，perplexity困惑度参数的取值必须小于(nrow(data) - 1 )/ 3; theta参数取值越大，结果的准确度越低，默认值为0.5，max_iter参数设置最大迭代次数。
### pca=F,参数表示是否对输入的原始数据进行PCA分析，然后使用PCA得到的topN主成分进行后续分析，t-SNE算法的计算量是特别大的，对于维度较高的数据数据，先采用PCA降维可以有效提高运行的效率，
### 默认采用top50的主成分进行后续分析，当然也可以通过initial_dims参数修改这个值

