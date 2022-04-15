##CIBERSORT##
##wangyingrui
##20211215  order by celltype (B cell, T cell...)
##CVDHZ impute with minimal of each protein

#install.packages("remotes")
#remotes::install_github("icbi-lab/immunedeconv")
setwd("D:/lab-wyr/lab/groups/CVDHZ/data_analysis/immunedeconv/20cells")
source("E:/PCSHA-DIANN/ProteomeCommonCode_20201202.R",encoding = 'UTF-8')
library(EPIC)
library(immunedeconv)
library(dplyr)
library(ggplot2)
library(reshape2)
library(tidyr)
library(immunedeconv)
library(tibble)
library(readxl)
library(ggpubr)
library(stringr)
library(viridisLite)
library(viridis)
library(pheatmap)
source("CIBERSORT.R")
library(e1071)
library(parallel)
library(preprocessCore)
library(RColorBrewer)
source("E:/proteomics_library_linda.R")
set_cibersort_binary("/path/to/CIBERSORT.R")
set_cibersort_mat("/path/to/LM20.txt")
df_PBMC3 <- read.csv("input/df_PBMC3_fill_with_min.csv",row.names = 1) 
df_PBMC3 <- df_PBMC3[apply(df_PBMC3, 1, function(x) sum(!is.na(x))>0),]
df_PBMC3 =df_PBMC3[apply(df_PBMC3, 1, function(x) {(sum(is.na(x))/length(x)) <0.9}),]  ##7342*492
linda.plot.density(df_PBMC3,outlier = T)

#outlier 值处理
df2.matrix <- as.matrix(df_PBMC3)
outline1 <- (fivenum(df2.matrix)[4]-fivenum(df2.matrix)[2])*2+fivenum(df2.matrix)[4]
outline2 <- fivenum(df2.matrix)[2]-(fivenum(df2.matrix)[4]-fivenum(df2.matrix)[2])*2
df2.matrix[df2.matrix>outline1] <- outline1
df2.matrix[df2.matrix<outline2] <- outline2
linda.plot.density(df2.matrix)

df <- as.data.frame(df2.matrix)
tdf <- as.data.frame(t(df))
tdf$MS_ID <- row.names(tdf)

##input sample info and arrange by group and TIME 
info_PBMCsample2 <- read.csv("input/info_PBMCsample2.csv")
dfall <- inner_join(info_PBMCsample2,tdf)  ##492*7356

dfall <- arrange(dfall,dfall$group,dfall$TIME)
row.names(dfall) <- dfall$MS_ID
repn <- dfall[dfall$replicate=="rep",]
repnn <- row.names(repn)

dfall2 <- as.data.frame(t(dfall))
infonew <- dfall[,1:14]
write.csv(infonew,file ="output/infonew1215_m.csv")
write.csv(dfall2,file="output/dfall_with_infonew1215_m.csv")
dfall2norep = dfall[dfall$replicate==c("no"),]
dfall22norep = as.data.frame(t(dfall2norep))  ##7356*478
write.csv(dfall22norep ,file="output/dfall_with_infonew1215_norep_m.csv")

df_m <- dfall22norep[15:nrow(dfall2),]
write.csv(df_m,file="output/pmatrix_all_arrange_by_group_time_1215_m.csv")

##change Uniprot ID to HGNC symbol
df_m = read.csv(file="output/pmatrix_all_arrange_by_group_time_1215_m.csv")
HGNC_ID = read.table(file="input/PBMC_HGNC_ID.txt",header =T)
HGNC_symbol = read.table(file="input/PBMC_HGNC_symbol.txt",sep ="\t")
names(HGNC_symbol) = c("To","Symbol")
res_mcp_counter = deconvolute(dataset_racle$expr_mat, "mcp_counter")
names(HGNC_ID) = c("X","To")
HGNC <- inner_join(HGNC_ID,HGNC_symbol) ###6331 proteins to 6359 proteins
df_mname<- inner_join(HGNC,df_m)
df_mname= df_mname[!duplicated(df_mname$Symbol),]
row.names(df_mname) <- df_mname$Symbol  #6352*481
write.csv(df_mname,file="output/allname_pmatrix_6352_1215_m.csv")
df <- df_mname[,4:ncol(df_mname)]  
write.table(df,file="pmatrix_6352DEPs_3ttest_norep_m.txt",sep ="\t")
dfname <- df_mname[,1:3]


#####ttest 并集deps###############################################################3
##input DEPs，最后再筛， 3 times group12 vs group0
tt1DEPs <- read.csv("input/Time1group12_group0_all_volcano.csv",row.names = 1)
tt1 <- tt1DEPs[,(ncol(tt1DEPs)-2):ncol(tt1DEPs)]
tt1$DEPs <- row.names(tt1)
tt2DEPs <- read.csv("input/Time2group12_group0_all_volcano.csv",row.names = 1)
tt2 <- tt2DEPs[,(ncol(tt2DEPs)-2):ncol(tt2DEPs)]
tt2$DEPs <- row.names(tt2)
tt3DEPs <- read.csv("input/Time3group12_group0_all_volcano.csv",row.names = 1)
tt3 <- tt3DEPs[,(ncol(tt3DEPs)-2):ncol(tt3DEPs)]
tt3$DEPs <- row.names(tt3)
tt123DEPs = rbind(tt1,tt2,tt3) ##12162
#write.csv(tt123DEPs,file="output/tt123DEPS_allttest_for_IPA.csv")
tt123DEPs <- tt123DEPs[tt123DEPs$P_value < 0.05,] ##818
tt123DEPs <- tt123DEPs[!is.na(tt123DEPs$P_value),]  ##818
tt123DEPs = tt123DEPs[!duplicated(tt123DEPs$DEPs),]   ##749

dfDEPs <- dfall2[15:nrow(dfall2),]
dfDEPs$DEPs = row.names(dfDEPs)
pDEPs <- inner_join(tt123DEPs,dfDEPs)
row.names(pDEPs) <- pDEPs$DEPs
pDEPs = pDEPs[,5:ncol(pDEPs)]  ##746*492

pDEPs2 = as.matrix(pDEPs)
min2 = min(pDEPs2,na.rm = TRUE) #0.01

mat2 <- pDEPs
mat2$X = row.names(mat2)
df3t = left_join(mat2,dfname)
which(is.na(df3t$Symbol))  #601
df3t2 = df3t[-601,]
which(is.na(df3t2$Symbol))
row.names(df3t2) = df3t2$Symbol
df3t2 = df3t2[,1:492]

df3tnorep <- df3t2[,-which(names(df3t2) %in% repnn)]  ##delete replicates
write.table(df3tnorep,file="pmatrix_746DEPs_3ttest_norep0311_20.txt",sep ="\t")
infonew <- infonew[-which(row.names(infonew) %in% repnn),]

library(RColorBrewer)
results_depso <- CIBERSORT('LM20.txt',"pmatrix_746DEPs_3ttest_norep0311_20.txt", perm=1000, QN=F)
re <- results_depso[,-(21:23)]
re.matrix = as.matrix(re)

mypalette <- colorRampPalette(brewer.pal(8,"Set1"))

dat <- re %>% as.data.frame() %>%
  rownames_to_column("Sample") %>%
  gather(key = Cell_type,value = Proportion,-Sample)
x <- factor(as.integer(rownames(dat)),labels= dat$Sample)


#给 cell_type加了一个顺序，按照cibersort输出细胞类型的顺序
plotOrder <- c(colnames(re))
dat$Cell_type <- factor(dat$Cell_type,levels = plotOrder) 

p1 <-  
  ggplot(data = dat,aes(x,Proportion,fill = Cell_type))+ 
  geom_bar(stat ="identity",position = 'stack')+
  labs(fill = "Cell_type",x = "",y = "Estiamted Proportion") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") +
  scale_y_continuous(expand = c(0.01,0)) +
  scale_fill_manual(values = mypalette(22)) +
  theme(axis.text.x = element_text(
    angle = 100,
    hjust =1,
    vjust =0.5
    
  ))
p1 

ggsave("output/3ttest746DEPs_all_samples_20.pdf",p1,width = 10,height=6,units = c("in"))

infonew$Sample = infonew$MS_ID
dat = left_join(dat,infonew)
dat$Group = dat$group
dat$newlabel=gsub("group1","group2",dat$Group)

library(ggsignif)
library(ggpubr)
my_comparisons <- list(c("group0","group1"), c("group1","group2"),c("group0", "group2"))
bartlett.test(dat$Proportion ~ dat$Group)
#二分类 3time

p22 <- ggbarplot(dat, x = "Cell_type", y = "Proportion", add = "mean_se",
                 fill = "newlabel", palette = mypalette(22)[c(4,13)], 
                 position = position_dodge(0.8))+
  theme(axis.text.x = element_text(angle=80,vjust = 0.5))
p22 <- p22 +
  stat_compare_means(aes(group = newlabel, label = ..p.signif..),hide.ns = T,method = "wilcox.test", label.y = 0.2) 
p22
ggsave("output/bar_mean_3ttest746DEPs_2groups_1215_wilcoxtest_20.pdf",p22,width = 10,height=6,units = c("in"))

#三分类 barplot

p33 <- ggbarplot(dat, x = "Cell_type", y = "Proportion", add = "mean_se",
                 fill = "Group", palette = mypalette(22)[c(4,13,10)], 
                 position = position_dodge(0.8))+
  theme(axis.text.x = element_text(angle=80,vjust = 0.5))
p33 <- p33 +
  stat_compare_means(aes(group = newlabel, label = ..p.signif..),hide.ns = T,method = "kruskal.test", label.y = 0.2) 
p33
ggsave("output/bar_mean_3ttest746DEPs_3groups_1215_kruskal.test_20.pdf",p33,width = 10,height=6,units = c("in"))


#6types 区分时间点and group
dat$type6 = paste0(dat$TIME,"_",dat$newlabel) 
dat <- arrange(dat,dat$TIME,dat$group)

#6分类

my_comparisons <- list(c("time1_group0","time1_group2"), c("time2_group0","time2_group2"),
                       c("time3_group0", "time3_group2"))
options(repr.plot.width=4, repr.plot.height=4)

##6 colors
pba2 <- ggbarplot(dat, x = "type6", y = "Proportion", add = "mean_se",
                fill = "type6", facet.by = "Cell_type",palette =  c("#84BADB","#FFB600FF","#408DC3","#FF9200FF","#105BA2","#FF6D00FF")
                , short.panel.labs = FALSE)+
  #geom_bar()+
  geom_smooth(aes(group= newlabel,color = newlabel),method="loess",se=T,span=0.5, size=1,
             alpha=0.4,level = 0.95)+
  labs(title = "22 cells")+
  stat_compare_means(comparisons = my_comparisons,label.y = c(0.12,0.14,0.16,0.18,0.2, 0.22, 0.24,0.25,0.3),
                     hide.ns = T,label="p.signif",
                     method = "wilcox.test"
  )+
  theme_bw()+
  #theme_classic()+
  theme(axis.text = element_text(size = 15,colour = "black"),
        text = element_text(size = 15,colour = "black"),#调整label text的字体大小
        axis.title.x = element_text(size = 15,colour = "black")
  )+
  theme(axis.text.x = element_text(hjust = 1,angle = 45))+
  theme(panel.grid =element_blank()) +
  theme(axis.ticks = element_blank())
pba2

ggsave("output/barplot_20_3ttest746DEPs_6types_1215_wilcox6_20.pdf",pba2,width = 20,height=25,units = c("in"))


#9分类
dat$type9 = paste0(dat$TIME,"_",dat$group)
write.csv(dat,file="output/746DEPs_results_1215.csv")
dat = read.csv(file="output/746DEPs_results_1215.csv",row.names = 1)
my_comparisons <- list(c("time1_group0","time1_group1","time1_group2"), c("time2_group0","time2_group1","time2_group2"),
                       c("time3_group0", "time3_group1","time3_group2"))

##堆积图 6types

infonew$type <- paste0(infonew$TIME,"_",infonew$group)
infonew$sampleID <- row.names(infonew)
infonew$group2 <- infonew$group  
infonew$newlabel = gsub("[12]","1+2",infonew$group2) ##替换group1，group2 to group12 []or
infonew$type6 <- paste0(infonew$TIME,"_",infonew$newlabel)
results_depso 
tresults_depso <- as.data.frame(t(results_depso))

names(tresults_depso) = infonew$type6[match(names(tresults_depso),infonew$sampleID)]
df.mfu <- aggregate(t(tresults_depso),by = list(names(tresults_depso)),mean,na.rm=T)
df.mfu2 <- data.frame(t(df.mfu[,-1]))
names(df.mfu2) <- df.mfu[,1]
df.mfu2 <- df.mfu2[1:20,]

write.table(df.mfu2,file="PBMC_746deps_6types_1215_group_TIME_aftercibersort.txt", sep = "\t")

re <- t(df.mfu2)
mypalette <- colorRampPalette(brewer.pal(8,"Set1"))

dat <- re %>% as.data.frame() %>%
  rownames_to_column("Sample") %>% 
  gather(key = Cell_type,value = Proportion,-Sample)

plotOrder <- c(colnames(re))
dat$Cell_type <- factor(dat$Cell_type,levels = plotOrder)  


pp6<-ggplot(dat,aes(Sample,Proportion,fill = Cell_type)) + 
  geom_bar(stat = "identity",position = 'stack') +
  labs(fill = "Cell Type",x = "",y = "Estiamted Proportion") + 
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") + 
  scale_y_continuous(expand = c(0.01,0)) +
  scale_fill_manual(values = mypalette(22)) +
  theme(axis.text.x = element_text(
    angle = 100,
    hjust =1,
    vjust =0.5
    
  ))
pp6
ggsave("output/3ttest746DEPs_6types_1215_aftercibersort_order.pdf",pp6,width = 4,height=5,units = c("in"))
write.csv(re,file="output/746_DEPs_6types_results.csv")

##堆积图 9types

infonew$type9 <- paste0(infonew$TIME,"_",infonew$group)
infonew$sampleID <- row.names(infonew)
tresults_depso <- as.data.frame(t(results_depso))

names(tresults_depso) = infonew$type9[match(names(tresults_depso),infonew$sampleID)]
df.mfu <- aggregate(t(tresults_depso),by = list(names(tresults_depso)),mean,na.rm=T)
df.mfu3 <- data.frame(t(df.mfu[,-1]))
names(df.mfu3) <- df.mfu[,1]
df.mfu3 <- df.mfu3[1:20,]
write.table(df.mfu2,file="PBMC_746deps_9types_1215_group_TIME_aftercibersort.txt", sep = "\t")

re <- t(df.mfu3)
mypalette <- colorRampPalette(brewer.pal(8,"Set1"))

dat <- re %>% as.data.frame() %>%
  rownames_to_column("Sample") %>% 
  gather(key = Cell_type,value = Proportion,-Sample)

plotOrder <- c(colnames(re))
dat$Cell_type <- factor(dat$Cell_type,levels = plotOrder)  


pp9<-ggplot(dat,aes(Sample,Proportion,fill = Cell_type)) + 
  geom_bar(stat = "identity") +
  labs(fill = "Cell Type",x = "",y = "Estiamted Proportion") + 
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") + 
  scale_y_continuous(expand = c(0.01,0)) +
  scale_fill_manual(values = mypalette(22)) +
  theme(axis.text.x = element_text(
    angle = 100,
    hjust =1,
    vjust =0.5
    
  ))
pp9
ggsave("output/3ttest746DEPs_9types_1215_aftercibersort_order.pdf",pp9,width = 5,height=5,units = c("in"))


