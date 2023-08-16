#function
firstup <- function(x){
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}
compound_class <- function(x){
  if(x %in% phenolic_compounds){
    return("Phenolic")
  }else if(x %in% flavonoids_compounds){
    return("Flavonoids")
  }else if(x %in% anthraquinones_compounds){
    return("Anthraquinones")
  }
}
deg_fun <- function(x){
  if(x >= 2){
    return("up")
  }else if(x <= -2){
    return("down")
  }
}
expirment_fun <- function(x){
  a <- c()
  for(i in expirment){
    a = append(a, grep(i, x))
  }
  if(length(a) == 0){
    return(FALSE)
  }else{
    return(TRUE)
  }
}
get_tpm <- function(x){
  for(j in c(1:6)){
    
  }
  for(i in c(c:row.names(miRNA$miRNA.family))){
    
  }
}

#1、绘制小RNA fold图
micro_fold <- read.csv("./data/micro_fold.csv")
micro_fold$miRNA <- sapply(micro_fold$miRNA, function(x) gsub(" ", "", x))
micro_fold$miRNA <- sapply(micro_fold$miRNA, function(x) gsub("MI", "mi", x))
micro_fold$miRNA <- sapply(micro_fold$miRNA, function(x) gsub("Mi", "mi", x))
micro_fold$miRNA <- sapply(micro_fold$miRNA, function(x) gsub("hsa-", "", x))
micro_fold$miRNA <- sapply(micro_fold$miRNA, function(x) gsub("mmu-", "", x))
micro_fold$miRNA <- sapply(micro_fold$miRNA, function(x) gsub("<a8>C5p", "", x))
micro_fold$miRNA <- sapply(micro_fold$miRNA, function(x) gsub("L", "l", x))
micro_fold$miRNA <- sapply(micro_fold$miRNA, function(x) gsub("pre-", "", x))
micro_fold$miRNA <- sapply(micro_fold$miRNA, function(x) gsub("r", "R", x))
micro_fold$compound <- sapply(micro_fold$compound, function(x) firstup(x))

#plot
library(ggplot2)
library(ggsci)
ggplot(data = micro_fold, aes(x=Fold.chance..up.down., y=miRNA, colour=compound)) +
  geom_point(alpha = 0.7) +
  geom_vline(xintercept = c(-1,1), colour="red", size=0.5) +
  labs(color="") +
  xlab("FC") +
  theme_bw() +
  scale_color_igv()
ggsave("./picture/micro_fold.pdf", height=20, width=10)

#2、venn图
phenolic_compounds <- c("Syringic acid", "Vanillic acid", "Ferulic acid",
                        "Gallic acid", "Protocatechuic acid","Proanthocyanidins",
                        "Resveratrol")
flavonoids_compounds <- c("Quercetin", "Flavonol glycosides", "Hyperoside",
                          "Morin", "Kaempferol","Rutin")
anthraquinones_compounds <- c("Emodin")
micro_fold$class <- sapply(micro_fold$compound, function(x) compound_class(x))
micro_fold$deg <- sapply(micro_fold$Fold.chance..up.down., function(x) deg_fun(x) )
venn_micro <- micro_fold[micro_fold$deg == "down" | micro_fold$deg == "up", ]
#plot
library(VennDiagram)
library(RColorBrewer)
myCol <- brewer.pal(3, "Pastel2")
venn.diagram(
  x = list(venn_micro[venn_micro$class == "Phenolic",2], 
           venn_micro[venn_micro$class == "Flavonoids",2], 
           venn_micro[venn_micro$class == "Anthraquinones",2]),
  category.names = c("Phenolic" , "Flavonoids " , "Anthraquinones"),
  filename = './picture/venn_plot2.png',
  output=TRUE,
  
  # 设置输出：
  imagetype="png" ,
  height = 1000 , 
  width = 1000 , 
  resolution = 300,
  compression = "lzw",
  
  # 圆的调整：
  lwd = 2, # 描边粗细
  lty = "blank",  # 去掉描边
  fill = myCol,  # 填充颜色
  
  # 文字大小：
  cex = .5,  # 大小；
  fontface = "bold",  # 粗体
  fontfamily = "sans",  # 字体
  
  # 每个集合的名称：
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",  # 集合名称位置：outer -- 外部；text -- 内部
  cat.pos = c(180, 180, 180),  # 集合名称分别在圈圈的什么角度
  cat.dist = c(0.055, 0.055, 0.075),  # 外部多少距离
  cat.fontfamily = "sans",
  rotation = 1  # 旋转
)
#3、上调和下调的前6个miRNA对应的内源靶标基因的GO和KEGG分析
library("readxl")
downandup_gene <- c("miR-122", "miR-33", "miR-34a",
                    "miR-21", "miR-15a", "miR-9")
downup_micro <- micro_fold[sapply(micro_fold$miRNA, function(x) x %in% downandup_gene), ]
up_gene <- unique(downup_micro[downup_micro$deg == "up",]$miRNA)
up_gene <- paste("hsa-",up_gene, sep="")
for(i in up_gene){
  up_gene <- append(up_gene,paste(i,"-5p",sep=""))
  #up_gene <- append(up_gene,paste(i,"-3p",sep=""))
}
down_gene <- unique(downup_micro[downup_micro$deg == "down",]$miRNA)
down_gene <- paste("hsa-",down_gene, sep="")
for(i in down_gene){
  down_gene <- append(down_gene,paste(i,"-5p",sep=""))
  #down_gene <- append(down_gene,paste(i,"-3p",sep=""))
}
hsa <- read_excel("E:/workspace/R/case4_micro/data/hsa_MTI.xlsx")
#获取所有up_gene所对应的基因
expirment <- c("qRT-PCR", "Western blot")
hsa_up_gene <- hsa[sapply(hsa$miRNA, function(x) x %in% up_gene),]
hsa_up_gene <- hsa_up_gene[sapply(hsa_up_gene$Experiments, function(x) x %in% expirment),]
hsa_up_gene <- hsa_up_gene[!duplicated(hsa_up_gene$`Target Gene`),]
write.csv(hsa_up_gene, "./data/hsa_up_gene.csv")
hsa_down_gene <- hsa[sapply(hsa$miRNA, function(x) x %in% down_gene),]
hsa_down_gene <- hsa_down_gene[sapply(hsa_down_gene$Experiments, function(x) x %in% expirment),]
hsa_down_gene <- hsa_down_gene[!duplicated(hsa_down_gene$`Target Gene`),]
write.csv(hsa_down_gene, "./data/hsa_down_gene.csv")

#plot
library(GOplot)
library(enrichplot)
#up
hsa_up_gene <- read.csv("./data/hsa_up_gene.csv", row.names = 1)
go_up_gene <- read.csv("E:/workspace/R/case4_micro/data/gProfiler_hsapiens_hsa_up_gene.csv")
go_up_gene <- go_up_gene[go_up_gene$source == "GO:MF" | go_up_gene$source == "KEGG", ]
up_plot <- hsa_up_gene[,c(4,2)]
for(i in c(1:length(row.names(go_up_gene)))){
  gene_list = strsplit(go_up_gene[i,11],",")[[1]]
  for(j in gene_list){
    micro_rna = hsa_up_gene[hsa_up_gene$Target.Gene == j,][1,2]
    up_plot = rbind(up_plot, c(go_up_gene[i,3], as.data.frame(micro_rna)[1,1]))
  }
}
up_plot <- up_plot[!duplicated(up_plot), ]
up_plot <- up_plot[!is.na(up_plot$miRNA), ]
up_plot$"logFC" = 2
colnames(up_plot) <- c("genes", "term", "logFC")
write.csv(up_plot, "data/up_plot.csv")
#down
hsa_down_gene <- read.csv("./data/hsa_down_gene.csv", row.names = 1)
go_down_gene <- read.csv("E:/workspace/R/case4_micro/data/gProfiler_hsapiens_hsa_down_gene.csv")
go_down_gene <- go_down_gene[go_down_gene$source == "GO:MF" | go_down_gene$source == "KEGG", ]
down_plot <- hsa_down_gene[,c(4,2)]
for(i in c(1:length(row.names(go_down_gene)))){
  gene_list = strsplit(go_down_gene[i,11],",")[[1]]
  for(j in gene_list){
    micro_rna = hsa_down_gene[hsa_down_gene$Target.Gene == j,][1,2]
    down_plot = rbind(down_plot, c(go_down_gene[i,3], as.data.frame(micro_rna)[1,1]))
  }
}
down_plot <- down_plot[!duplicated(down_plot), ]
down_plot <- down_plot[!is.na(down_plot$miRNA), ]
down_plot$"logFC" = 2
colnames(down_plot) <- c("genes", "term", "logFC")
write.csv(down_plot, "data/down_plot.csv")


#4、文献中报导的最多的荞麦衍生的miRNA
miRNA <- read.csv("E:/workspace/R/case4_micro/data/conserv_microrna.csv")
miRNA$sum <- apply(miRNA[,c(5:10)], 1, function(x) sum(x))
mi_sum <- aggregate(sum ~ miRNA.family, data=miRNA, sum)
mi_sum$percent <- (mi_sum[,2]/sum(mi_sum[,2]))*100
mi_sum_8 <- mi_sum[order(mi_sum$percent, decreasing = TRUE),][c(1:8),]
#miRNA$percent <- (miRNA[,11]/sum(miRNA[,11]))*100
# miRNA_copy <- miRNA
# miRNA <- miRNA[,c(1,4,5,6,7,8,9,10)]
#转换成tpm数据
# kb <- miRNA$Length..bp. / 1000
# miRNA[, c(3:8)] <- miRNA[, c(3:8)] / kb
# miRNA[, c(3:8)] <- t(t(miRNA[, c(3:8)])/colSums(miRNA[, c(3:8)])*1000000)
# miRNA$sum <- apply(miRNA[,c(5:10)], 1, function(x) sum(x))
#相同的miRNA 相加
#每个miRNA的占比
#plot柱状图
ggplot(mi_sum_8, aes(x=reorder(miRNA.family, percent, decreasing=TRUE), y=percent)) +
  geom_bar(stat="identity") +
  ylab("percent") +
  xlab("miRNA") +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5)) +
  labs(title = "miRNA abundance") +
  theme(plot.title = element_text(hjust = 0.5)) 
ggsave("./picture/mirna_abundance.pdf", height=5, width=5)
#5、miRNA 大于百分之二的miRNA
mi_sum <- mi_sum[mi_sum$percent > 0.3, ]
need_miRNA <- merge(mi_sum, miRNA, by="miRNA.family")
write.csv(need_miRNA, "./data/need_miRNA.csv")


#5、外源miRNA-mrna算法预测
library(stringr)
pita_results <- read.table("E:/workspace/R/case4_micro/out_miRNA/result_pita_results.tab",
                           sep = "\t", header = TRUE)
pita_results[,1] <- sapply(pita_results$UTR, function(x) str_extract(x, "Seq_.*? "))
pita_results[,1] <- sapply(pita_results$UTR, function(x) strsplit(x," ")[[1]][1])
pita_results[,1] <- sapply(pita_results$UTR, function(x) substr(x,5,20))
#得到蛋白名称
write.csv(unique(pita_results$UTR), "./out_miRNA/protein_name.csv")
gProfiler <- read.csv('./out_miRNA/gProfiler_hsapiens_2023-5-7 15-32-15.csv')
colnames(gProfiler)[1] <- "UTR"
a <- merge(pita_results[,c(1,2,13)],gProfiler[,c(1,2)], by="UTR")
b <- a
b <- b[b$converted_alias != "None",]
b <- b[!duplicated(b[,c(2,4)]),]
b <- b[grep("LOC", b$converted_alias, invert = TRUE), ]

pita_results_100 <- b[1,]
micro_rna <- unique(b$microRNA)
for(i in micro_rna){
  pita_results_100_test <- b[b$microRNA == i,][order(b[b$microRNA == i,]$ddG),]
  pita_results_100_test <- pita_results_100_test[c(1:100),]
  pita_results_100 <- rbind(pita_results_100, pita_results_100_test)
}
pita_results_100 <- pita_results_100[-1,]
write.csv(pita_results_100, "./out_miRNA/node.csv")
#看外源和内源的交集
inner_gene <- read.csv("./out_miRNA/inner_gene.csv")
intersect(unique(b$converted_alias),unique(inner_gene$gene))
write.csv(intersect(inner_gene$gene, b$converted_alias), "./out_miRNA/jiaojigene.csv")


#找出每个microRNA 得分最高的对应靶标
pita_results_max <- pita_results[1,]
for(i in micro_rna){
  pita_results_max_test <- pita_results[pita_results$microRNA == i,][which.min(pita_results[pita_results$microRNA == i,]$ddG),]
  pita_results_max <- rbind(pita_results_max, pita_results_max_test)
}
pita_results_max <- pita_results_max[-1,]
write.csv(pita_results_max, "./out_miRNA/pita_results_max.csv")

#plot富集弦图
library(GOplot)
gene_GO <- read.csv("E:/workspace/R/case4_micro/out_miRNA/nodefuji.csv")
gene_GO <- gene_GO[order(gene_GO$adjusted_p_value),][c(1:30),]
GOplotIn_BP<-gene_GO[,c(1,2,4,10)]
names(GOplotIn_BP)<-c('ID','Term','adj_pval','Genes')
GOplotIn_BP$Category = "BP"
#构建基因data数据
genedata <- data.frame(ID=c(), logFC=c())
for(i in c(1:length(row.names(GOplotIn_BP)))){
  IDs <- strsplit(GOplotIn_BP[i,4],",")
  log2FoldChange <- GOplotIn_BP[i,3]
  genedatas <- data.frame(ID=IDs[[1]], logFC=log2FoldChange)
  genedata <- rbind(genedata, genedatas)
}
genedata <- genedata[!duplicated(genedata$ID),]
circ_BP<-GOplot::circle_dat(GOplotIn_BP, genedata)
chord_BP<-chord_dat(data = circ_BP,genes = genedata)
chord_BP[,"logFC"] <- runif(length(row.names(chord_BP)), min=-2, max=2)
pdf(file="./picture/go_chordssss.pdf",width=30,height=30)
GOChord(data = chord_BP[c(1:100),],#弦图
        title = 'GO',space = 0.05,#GO Term间距
        limit = c(3,5),gene.order = 'logFC',gene.space = 0.25,gene.size = 5,
        border.size = 1,lfc.col = c('red','white','blue'), #上下调基因颜色
        process.label = 20) #GO Term字体大小
dev.off()


###########内源外源头目标基因venn图
library(VennDiagram)
library(RColorBrewer)
myCol <- c("#D7D7D7", "#A0A0A0")
venn.diagram(
  x = list(unique(b$converted_alias),unique(inner_gene$gene)),
  category.names = c("Exogenous" , "Endogenous"),
  filename = './picture/venn_in_out.png',
  output=TRUE,
  
  # 设置输出：
  imagetype="png" ,
  height = 1000 , 
  width = 1000 , 
  resolution = 300,
  compression = "lzw",
  
  # 圆的调整：
  lwd = 2, # 描边粗细
  lty = "blank",  # 去掉描边
  fill = myCol,  # 填充颜色
  
  # 文字大小：
  cex = .5,  # 大小；
  fontface = "bold",  # 粗体
  fontfamily = "sans",  # 字体
  
  # 每个集合的名称：
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",  # 集合名称位置：outer -- 外部；text -- 内部
  cat.pos = c(180, 180),  # 集合名称分别在圈圈的什么角度
  cat.dist = c(0.055, 0.055),  # 外部多少距离
  cat.fontfamily = "sans"
)
###########内源外源头目标基因GO富集绘图
GO_data <- read.csv("E:/workspace/R/case4_micro/out_miRNA/jiaojigofuji.csv")
BP <- GO_data[GO_data$source == "GO:BP",]
BP <- BP[order(BP$adjusted_p_value),][c(1:10),]
BP$adjusted_p_value <- -log10(BP$adjusted_p_value)
CC <- GO_data[GO_data$source == "GO:CC",]
CC <- CC[order(CC$adjusted_p_value),][c(1:10),]
CC$adjusted_p_value <- -log10(CC$adjusted_p_value)
MF <- GO_data[GO_data$source == "GO:MF",]
MF <- MF[order(MF$adjusted_p_value),][c(1:10),]
MF$adjusted_p_value <- -log10(MF$adjusted_p_value)
ggplot(BP, aes(x = adjusted_p_value, y = reorder(term_name, adjusted_p_value))) +
  geom_bar(stat="identity", position="identity", width = 0.5, fill = "#ABABAB") +
  labs(title = "BP", x = "-log10(Pval Adj)", y = "Numbers of genes") +
  theme(text = element_text(size = 20), plot.title = element_text(hjust = 0.5))
ggsave("bp.pdf", width=20,height=15)
ggplot(CC, aes(x = adjusted_p_value, y = reorder(term_name, adjusted_p_value))) +
  geom_bar(stat="identity", position="identity", width = 0.5, fill = "#ABABAB") +
  labs(title = "CC", x = "-log10(Pval Adj)", y = "Numbers of genes") +
  theme(text = element_text(size = 20), plot.title = element_text(hjust = 0.5))
ggsave("cc.pdf", width=20,height=15)
ggplot(MF, aes(x = adjusted_p_value, y = reorder(term_name, adjusted_p_value))) +
  geom_bar(stat="identity", position="identity", width = 0.5, fill = "#ABABAB") +
  labs(title = "MF", x = "-log10(Pval Adj)", y = "Numbers of genes") +
  theme(text = element_text(size = 20), plot.title = element_text(hjust = 0.5))
ggsave("mf.pdf", width=20,height=15)
