library(dplyr)
library(stringr)
library(UCSCXenaTools)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(randomForest)
library(tidyverse)
library(skimr)
library(DataExplorer)
library(caret)
library(pROC)
library(KEGGREST)
library(org.Hs.eg.db)
library(pheatmap)
library(tibble)
library(depmap)
library(ExperimentHub)
library(corrplot)
library(ggstatsplot)
library(ggside)
library(ggrepel)
library(survival)
library(survminer)
library(clusterProfiler)
library(fgsea)

#### 获取数据 ####
# 缺氧评分
hyScore <- read.table("~/FLAD1/HypoxiaScore/41588_2018_318_MOESM3_ESM.txt", header = T)  # 数据来源:https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-018-0318-2/MediaObjects/41588_2018_318_MOESM3_ESM.txt

# 肿瘤类型
cancerType <- unique(hyScore$tumour_type)
cancerType <- cancerType[-4]
cancerType[19:20] <- c("COAD", "READ")
cancerType

# 基因表达 fpkm
df_mRNA <- XenaData %>%
  filter(XenaHostNames == "gdcHub") %>%
  filter(Label == "HTSeq - FPKM")
df_mRNA$tumor <- gsub(".*\\((.*?)\\).*", "\\1", df_mRNA$XenaCohorts)
df_mRNA <- df_mRNA[df_mRNA$tumor %in% cancerType,]

df_mRNA <- df_mRNA %>%
  XenaGenerate() %>%  
  XenaQuery() %>% 
  XenaDownload() %>%
  XenaPrepare()

FLAD1_tumor_mRNA_score <- FLAD1_mRNA %>%
  dplyr::filter(sampletype == "tumor")
FLAD1_tumor_mRNA_score$id <- gsub(pattern = "-",replacement = ".", x = substr(row.names(FLAD1_tumor_mRNA_score),1,12) )   
FLAD1_tumor_mRNA_score <- FLAD1_tumor_mRNA_score[order(FLAD1_tumor_mRNA_score$id, -FLAD1_tumor_mRNA_score$FLAD1),]
FLAD1_tumor_mRNA_score <- FLAD1_tumor_mRNA_score[!duplicated(FLAD1_tumor_mRNA_score$id),]
FLAD1_tumor_mRNA_score <- merge(FLAD1_tumor_mRNA_score, hyScore, by.x = "id", by.y = "patient_id")

# 基因表达 counts
CancerType <- XenaData %>% filter(Label == "HTSeq - Counts")
CancerType <- CancerType[grep(pattern = "TCGA", x = CancerType$XenaDatasets),]
CancerType <- str_extract(string = str_extract(string = CancerType$XenaCohorts, 
                                               pattern = "\\([A-Z]+\\)"),
                          pattern = "[A-Z]+") 

list_counts_all <- list()
getcounts <- function(tumor){
  data("XenaData")
  counts <- XenaData %>% dplyr::filter(Label == "HTSeq - Counts")
  counts <- counts[grep(pattern = tumor, x = counts$XenaDatasets),]
  counts <- XenaGenerate(XenaData = counts) %>%
    XenaQuery() %>%
    XenaDownload() %>%
    XenaPrepare()
  return(counts)
}

for (i in 1:length(CancerType)) {
  list_counts_all[[i]] <- getcounts(tumor = CancerType[i])
  names(list_counts_all)[i] <- CancerType[i]
}

do_DEG <- function(counts, groupA, groupA_label, groupB, groupB_label){
  colnames(counts)[1] <- "Ensembl_ID"
  counts_rownames <- counts$Ensembl_ID 
  counts <- counts[,-1]
  counts <- 2^counts - 1 #log2(count+1)
  
  condition <- factor(x = c(rep(groupA_label, length(groupA)),
                            rep(groupB_label, length(groupB))),
                      levels = c(groupA_label, groupB_label))
  
  colData <- data.frame(row.names = c(groupA, groupB),
                        condition = condition)
  
  counts <- counts %>%
    dplyr::select(c(groupA, groupB))
  
  counts <- round(counts, digits = 0)
  counts <- as.data.frame(counts)
  row.names(counts) <- counts_rownames
  
  dds <- DESeqDataSetFromMatrix(countData = counts, 
                                colData = colData, 
                                design = ~ condition)
  dds <- DESeq(dds)
  res <- results(dds)
  
  return(res)
}

# FLAD1 fpkm
XenaData <- XenaData %>%
  filter(XenaHostNames == "gdcHub" & Label == "HTSeq - FPKM")

fetchfpkm <- function(tumor, gene){
  dataset <- paste0("TCGA-", tumor, ".htseq_fpkm.tsv")
  exp <- fetch_dense_values(host = "https://gdc.xenahubs.net", 
                            dataset =  dataset, 
                            identifiers = gene,
                            use_probeMap = TRUE)
  exp <- as.data.frame(t(exp))
  exp$type <- tumor
  
  return(exp)
} 

FLAD1_mRNA <- lapply(cancerType, function(tumor){
  fetchfpkm(tumor = tumor, gene = "FLAD1")
})

FLAD1_mRNA <- do.call("rbind", FLAD1_mRNA)
FLAD1_mRNA$sampletype <-ifelse(substr(row.names(FLAD1_mRNA),14,15) %in% c(paste0("0",1:9)), 
                               "tumor",
                               ifelse(substr(row.names(FLAD1_mRNA),14,15) %in% c(paste0("1",0:9)), 
                                      "normal", 
                                      "control")
)


# FLAD1 cnv
fetchcnv<- function(tumor, gene){
  dataset <- paste0("TCGA-", tumor, ".gistic.tsv")
  exp <- fetch_dense_values(host = "https://gdc.xenahubs.net", 
                            dataset =  dataset, 
                            identifiers = gene,
                            use_probeMap = TRUE)
  exp <- as.data.frame(t(exp))
  exp$type <- tumor
  return(exp)
} 

FLAD1_cnv <- lapply(cancerType, function(tumor){
  fetchcnv(tumor = tumor, gene = "FLAD1")
})
FLAD1_cnv <- do.call("rbind", FLAD1_cnv)
FLAD1_cnv$id <- row.names(FLAD1_cnv)
FLAD1_cnv <- FLAD1_cnv[which(substr(FLAD1_cnv$id, 14,15) %in% c(paste0("0",1:9))),]
FLAD1_cnv$id <- gsub("-",".",substr(FLAD1_cnv$id, 1,12)) 
length(unique(FLAD1_cnv$id)) 
FLAD1_cnv1 <- FLAD1_cnv[order(row.names(FLAD1_cnv)),]
FLAD1_cnv1$id[duplicated(FLAD1_cnv1$id)]
FLAD1_cnv1 <- FLAD1_cnv1[!duplicated(FLAD1_cnv1$id),]

FLAD1_cnv1 <- merge(FLAD1_cnv1, hyScore, by.x = "id", by.y = "patient_id")


#### figS3 ####
# * 1. 数据 ----
# 1. 代谢相关的基因集
hsa_path_gene  <- keggLink("pathway", "hsa") %>% 
  tibble(pathway = ., eg = sub("hsa:", "", names(.))) %>%
  mutate(
    symbol = mapIds(org.Hs.eg.db, eg, "SYMBOL", "ENTREZID")
  )

hsa_pathways <- keggList("pathway", "hsa") %>% 
  tibble(pathway = names(.), description = .)

hsa_path_gene$pathway <- gsub(pattern = "path:",
                              replacement = "",
                              x = hsa_path_gene$pathway)
path_gene <- merge(x = hsa_path_gene, y = hsa_pathways, 
                   by = "pathway",
                   all.x = T)

path_gene$description <- gsub(pattern = " - Homo sapiens \\(human\\)",
                              replacement = "",
                              x = path_gene$description)
length(path_gene$symbol[!duplicated(path_gene$symbol)])
path_gene$description[!duplicated(path_gene$description)]

# https://www.genome.jp/kegg/pathway.html
kegg_metabolism <- readxl::read_excel(path = "~/FLAD1/Metabolism/kegg_metabolism.xlsx",
                                      col_names = T)

kegg_metabolism <- merge(x = path_gene, y = kegg_metabolism, 
                         by.x = "description", by.y = "1.0 Global and overview maps")


# 2. 代谢基因fpkm
df_mRNA1 <- lapply(names(df_mRNA), function(name) {
  df <- df_mRNA[[name]]
  
  df <- merge(x = probe,y = df, by.x = "id", by.y = "Ensembl_ID" )
  
  df$mean <- apply(df[,-c(1,2)],1,mean)
  df <- df[order(df$gene, -df$mean),]
  df <- df[!duplicated(df$gene),]
  
  row.names(df) <- df$gene
  df <- df[,-which(colnames(df) %in% c("gene", "id", "mean"))]
  df <- as.data.frame(t(df))
  
  return(df)
})

df_mRNA2 <- do.call(rbind, df_mRNA1)
df_mRNA3 <- df_mRNA2[,colnames(df_mRNA2) %in% kegg_metabolism$symbol]

# 3. 肿瘤中上调基因
sample_type <- lapply(names(df_mRNA), function(name){
  df <- df_mRNA[[name]]
  df <- data.frame(
    sample = colnames(df)[-1],
    type = name
  )
  return(df)
})

sample_type <- do.call("rbind", sample_type)
sample_type$class <- ifelse(
  substr(sample_type$sample,14,15) %in% c(paste0("0",1:9)), "tumor",
  ifelse(
    substr(sample_type$sample,14,15) %in% c(paste0("1",0:9)), "normal",
    "other"
  )
)

unique(sample_type$type) %in% cancerType
cancerType1 <- table(sample_type$type, sample_type$class)
cancerType1 <- as.data.frame(cancerType1)
cancerType1 <- as.character(cancerType1[cancerType1$Var2 == "normal" & cancerType1$Freq < 4,]$Var1) # 去除normal 太少的样本
cancerType1 <- cancerType[!cancerType %in% cancerType1][-3]


DEG_NT <- lapply(cancerType1, function(type){
  tumorsample <- sample_type[sample_type$type == type & sample_type$class == "tumor",]$sample
  normalsample <- sample_type[sample_type$type == type & sample_type$class == "normal",]$sample
  
  tumorsample <- df_mRNA3[row.names(df_mRNA3) %in% tumorsample,]
  normalsample <- df_mRNA3[row.names(df_mRNA3) %in% normalsample,]
  
  df <- lapply(1:length(colnames(tumorsample)), function(x){
    data.frame(
      gene = colnames(tumorsample)[x],
      pval = wilcox.test(tumorsample[,x], normalsample[,x])$p.value,
      normalMean = mean(normalsample[,x]),
      tumorMean = mean(tumorsample[,x])
    )
  })
  df <- do.call("rbind", df)
  df$type <- type
  print(type)
  return(df)
})

DEG_NT <- do.call("rbind", DEG_NT)
DEG_NT$foldchange <- DEG_NT$tumorMean / DEG_NT$normalMean
DEG_NT$log2fc <- log(x = DEG_NT$foldchange,base = 2)

DEG_NT1 <- DEG_NT[DEG_NT$normalMean != 0,]

DEG_NT1 <- DEG_NT1[DEG_NT1$type != "PAAD",] 

res <- DEG_NT1 %>%
  group_by(gene) %>%  
  summarise(frequency = sum(pval < 0.05 & log2fc > 0), .groups = 'drop') 
res <- res[order(-res$frequency, res$gene),] # 肿瘤中显著上调的基因


# 4. 蛋白表达与mRNA水平一致基因
cor_RNA_protein  # 数据 https://doi.org/10.1016/j.cell.2020.08.036
res1 <- merge(x = res, y = cor_RNA_protein, by.x = "gene", by.y = "Gene", all.x = T)
res1 <- res1 %>%
  dplyr::filter(sig_indicator_for_Spearman_cor == "pos.sig") %>%
  dplyr::arrange(desc(frequency))

# 5. 合并缺氧得分与df_mRNA3
df_mRNA4 <- df_mRNA3
df_mRNA4$patient_id <- row.names(df_mRNA4)
tumorSample <- df_mRNA4$patient_id[which(substr(df_mRNA4$patient_id,14,15) %in% c(paste0("0",1:9)))]
df_mRNA4 <- df_mRNA4[df_mRNA4$patient_id %in% tumorSample,]
df_mRNA4$patient_id <- gsub("-",".",substr(df_mRNA4$patient_id, 1, 12))
df_mRNA4 <- merge(hyScore[,c(1,3)], df_mRNA4, by = "patient_id") # 选择 buffa/wintr/ragnum 缺氧评分

df_mRNA4 <- df_mRNA4[,-1]
colnames(df_mRNA4)[1]
colnames(df_mRNA4)[1] <- "hyScore"


# 6. 随机森林数据集df_mRNA5
df_mRNA5 <- df_mRNA4[,c(1,which(colnames(df_mRNA4) %in% res1$gene[c(1:49)]))]
df_mRNA5 <- scale(df_mRNA5)
df_mRNA5 <- as.data.frame(df_mRNA5)


# * 2. 数据拆分 ----
set.seed(42)
trains <- createDataPartition(y = df_mRNA5$hyScore, p = 0.75, list = F)
traindata <- df_mRNA5[trains,]
testdata <- df_mRNA5[-trains,]

hist(traindata$hyScore,breaks = 20)
hist(testdata$hyScore,breaks = 20)

colnames(traindata)
form_reg <- as.formula(
  paste0(
    "hyScore ~",
    paste(colnames(traindata)[2:length(colnames(traindata))], collapse = " + ")
  )
)
form_reg

# * 3. 训练模型 ----
set.seed(42)
fit_rf_reg <- randomForest(
  formula = form_reg,
  data = traindata,
  ntree = 500, 
  mtry = 6, 
  importance = T 
)

fit_rf_reg 
plot(fit_rf_reg, main = "ERROR & TREES")

a <- importance(fit_rf_reg)
a <- as.data.frame(a)
a$gene <- row.names(a)
a <- a[order(-a$IncNodePurity),]
a$gene <- factor(x = a$gene, levels = rev(a$gene))
# 分别使用winter, buffa, ragnum获得a_winter, a_buffa, a_ragnum

# * 4. 模型预测 ----
# 训练集
trainpred <- predict(fit_rf_reg, newdata = traindata)
defaultSummary(data.frame(obs = traindata$hyScore, pred = trainpred)) 

# 测试集预
testpred <- predict(fit_rf_reg, newdata = testdata)
defaultSummary(data.frame(obs = testdata$hyScore, pred = testpred))


# 训练集和测试集结果
preddresult <- data.frame(
  obs = c(traindata$hyScore, testdata$hyScore),
  pred = c(trainpred, testpred),
  group = c(rep("Train", length(trainpred)),
            rep("Test", length(testpred)))
)

# * 5. 结果展示 ----
# 1. 绘制相关性
cor_hypoxia <- hyScore[,2:4]
colnames(cor_hypoxia) <- c("Buffa", "Winter", "Ragnum")
corrplot.mixed(corr = cor(cor_hypoxia),
               lower = "number", upper = "pie",
               tl.pos = "d", tl.col = "black")

# 2. 绘制回归结果 
anno_train <- defaultSummary(data.frame(obs = traindata$hyScore, pred = trainpred))
anno_test <- defaultSummary(data.frame(obs = testdata$hyScore, pred = testpred))

anno_train <- paste0("Train: ","RMSE = ", round(anno_train[1],3), "  ",
                     "R2 = ", round(anno_train[2],3), "  ",
                     "MAE = ", round(anno_train[3],3))
anno_test <- paste0("Test:  ","RMSE = ", round(anno_test[1],3), "  ",
                    "R2 = ", round(anno_test[2],3), "  ",
                    "MAE = ", round(anno_test[3],3))

ggplot(preddresult,
       aes(x = obs, y = pred, color = group, fill = group)) +
  geom_point(size = 1,alpha = .6) +
  geom_smooth(method = "lm", se = T) +
  geom_abline(intercept = 0, slope = 1, size = 1.2 ) +
  geom_xsidedensity(aes(y=stat(density)), alpha = .1) +
  geom_ysidedensity(aes(x=stat(density)), alpha = .1) +
  scale_x_continuous(breaks = seq(-2,2,1)) +
  scale_color_manual(values = c(Train = "#FE817D", Test = "#81B8DF")) +
  scale_fill_manual(values = c(Train = "#FE817D", Test = "#81B8DF")) +
  scale_xsidey_continuous(breaks = seq(0,0.3,0.3)) +
  scale_ysidex_continuous(breaks = seq(0,0.3,0.3)) +
  xlab("Winter hypoxia score") + 
  ylab("Predicted hypoxia score") +
  annotate(geom = "text", x = -0.7, y = 2, 
           label = anno_train, 
           color = "#FE817D", size = 6)+
  annotate(geom = "text", x = -0.7, y = 1.7, 
           label = anno_test, 
           color = "#81B8DF", size = 6)+
  theme_bw() +
  theme(legend.position = "none",
        ggside.panel.scale = 0.2,
        axis.text.x = element_text(size = 30),
        axis.text.y = element_text(size = 30),
        axis.title = element_text(size = 30, face = "bold"),
  )  

# 3. 热图展示不同基因的在不同肿瘤中的表达情况
getdf_sample_gene_hypoxia <- function(genename){
  df <- df_mRNA3
  df$patient_id <- row.names(df)
  tumorSample <- df$patient_id[which(substr(df$patient_id,14,15) %in% c(paste0("0",1:9)))]
  df <- df[df$patient_id %in% tumorSample,]
 
  df$patient_id <- gsub("-",".",substr(df$patient_id, 1, 12))
  
  df <- merge(hyScore[,c(1,4)], df, by = "patient_id") # Ragnum缺氧评分
  colnames(df)[2] <- "hyScore"
  
  df <- df[order(-df$hyScore),]

  df <- df[,which(colnames(df) %in% genename)]
  
 
  return(df)
}


scale_to_range <- function(x, min_range, max_range){
  # 查找原始数据的最小值和最大值
  x_min <- min(x)
  x_max <- max(x)
  
  # 缩放到[min_range, max_range]
  scaled_x <-  (x - x_min) / (x_max - x_min) 
  
  return(scaled_x)
}

phe_sample_gene <- getdf_sample_gene_hypoxia(genename = unique(a_total$gene))
phe_sample_gene <- as.data.frame(lapply(phe_sample_gene, scale_to_range, min_range = -1, max_range = 1)) 

phe_sample_gene <- phe_sample_gene[,match(
  a_total[a_total$hyscore == "ragnum",]$gene,
  colnames(phe_sample_gene))]


pheatmap(phe_sample_gene,
         cluster_rows = F,
         cluster_cols = F,
         show_rownames = F,
         show_colnames = T,
         legend = T,
         fontsize = 20,
         angle_col = 90,
         main = "Ragnum")
skim(phe_sample_gene)

# 5. 基因列表(贡献度图)
a_total <- rbind(a_buffa[1:20,], a_winter[1:20,], a_ragnum[1:20,])

a_total <- a_total[a_total$gene %in% intersect(a_ragnum$gene[1:20], 
                                               intersect(a_buffa$gene[1:20], a_winter$gene[1:20])),]

a_total$hyscore <- factor(x = a_total$hyscore,
                          levels = c("ragnum", "winter", "buffa"))
a_total$gene <- factor(x = a_total$gene,
                       levels = rev(a_ragnum$gene))

ggplot(a_total, aes(x = gene, y = IncNodePurity)) +
  geom_segment(aes(x = gene, xend = gene, y = 0, yend = IncNodePurity), color="skyblue") +
  geom_point(aes(size = `%IncMSE`), color="blue", alpha=0.6) +
  theme_light() +
  coord_flip() +
  facet_grid(. ~ hyscore) +
  scale_y_continuous(breaks = seq(0,1000,500)) +
  theme_bw() +
  theme(
    legend.position="bottom",
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 20, face = "bold"),
    axis.title.y = element_blank(),
    strip.text = element_text(size = 20),
    strip.background = element_rect(fill = "white", color = "white")
  )



#### fig2  ####
# * 1. 计算基因与缺氧评分的相关性 ----
calcor <- function(tumor){
  df1 <- df_mRNA[[tumor]]
  colnames(df1)[1] <- "geneid"
  
  df1_row <- df1$geneid
  df1 <- df1[,-1]
  
  df1_col <- colnames(df1)[substr(colnames(df1),14,15) %in% c(paste0("0",1:9))] 
  df1_col <- df1_col[order(df1_col)] 
  df1 <- df1[,df1_col]
  
  colnames(df1) <- gsub(pattern = "-",replacement = ".", x = substr(colnames(df1),1,12)) 
  df1 <- df1[,!duplicated(colnames(df1))] 
  
  df1_col <- colnames(df1)
  df1_col <- intersect(df1_col, hyScore$patient_id) 
  df1 <- df1[,colnames(df1) %in% df1_col]
  df1 <- df1[,order(colnames(df1))]
  
  row.names(df1) <- df1_row
  
  df2 <- hyScore[hyScore$patient_id %in% df1_col,1:2] #选择buffa缺氧评分
  df2 <- df2[order(df2$patient_id),]
  colnames(df2)[2] <- "score"
  
  tmp_cor <- apply(df1, 1, function(x) cor(x, df2$score))
  res <- data.frame(gene = row.names(df1), cor = tmp_cor)
  res$type <- str_split_fixed(names(df_mRNA)[tumor],"[.]",5)[,2]
  res$num <- length(df2$patient_id)
  return(res)
}


res_tumor_HN <- lapply(1:20, function(x){
  calcor(tumor = x)
})

res_tumor_HN <- do.call("rbind", res_tumor_HN)

# * 2. 绘图 ----
# ** 1 tumor vs. normal ----
FLAD1_mRNA$type <- factor(x = FLAD1_mRNA$type, levels = sort(unique(FLAD1_mRNA$type)))


Num_tumor_and_normal <- as.data.frame(table(FLAD1_mRNA$type, FLAD1_mRNA$sampletype))
Num_tumor_and_normal <- dcast(Num_tumor_and_normal, Var1 ~ Var2, value.var = "Freq")
Num_tumor_and_normal <- Num_tumor_and_normal[order(Num_tumor_and_normal$Var1),]

ggboxplot(FLAD1_mRNA, x = "type", y = "FLAD1",
          fill = "sampletype",
          palette = c("#00AFBB", "#FC4E07"),
          outlier.shape = NA,
          ylim = c(1,6),
          ylab = "FLAD1 mRNA (log2(FPKM+1))")+
  stat_compare_means(aes(group=sampletype), label = "p.signif", 
                     label.y = 5.8, size = 6, angle = 90,
                     symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, Inf), 
                                        symbols = c("***", "**", "*", "ns"))
  ) +
  labs(fill = "Type") +
  scale_y_continuous(expand = c(0,0),
                     breaks = 1:6) +
  theme(axis.title.x = element_blank(),
        axis.title.y.left = element_text(size = 15, face = "bold"),
        axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 14),
        legend.position = "right",
        axis.text = element_text(size = 15, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        legend.text = element_text(face = "plain", size = 14),
        legend.title = element_text(face = "plain", size = 14),
        #axis.line = element_line(linetype=1,color="black",size=1.5)
  )+
  ggtitle("FLAD1 mRNA")+
  annotate(geom = "text", y=1.4, x=1:20, label=Num_tumor_and_normal$normal, size = 4.5)+
  annotate(geom = "text", y=1.2, x=1:20, label=Num_tumor_and_normal$tumor, size = 4.5)

# ** 2 hypoxia vs. normoxia ----
res_tumor_HN$type <- factor(x = res_tumor_HN$type,
                            levels = sort(unique(res_tumor_HN$type)))

FLAD1 <- res_tumor_HN[res_tumor_HN$gene == "ENSG00000160688.17",]
FLAD1 <- FLAD1[order(FLAD1$type),]


ggplot(res_tumor_HN, 
       aes(x = type, y = cor, fill = type)) +
  geom_boxplot(width = 0.3, outlier.shape = NA, color = "black", alpha = .5) +
  geom_point(data = FLAD1,
             aes(x = type, y = cor),
             shape = 23,   
             size = 3,     
             color = "red", 
             fill = "red" )+
  ylim(-0.5, 0.5)+
  ylab("Correlation")+
  annotate(geom = "text", y= -0.5, x=1:20, label=FLAD1$num, size = 4.5) +
  theme(panel.background = element_rect(color = "black", fill = "white"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold", size = 23),
        axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1),
        axis.text = element_text(face = "bold", size = 20),
        legend.position = "none")


# ** 3. 拷贝数 ----
getcnvbar <- function(a,score){
  
  num <- grep(score,colnames(a))
  
  a$group <- cut(
    a[, num],
    breaks = c(-Inf, quantile(a[, num], 0.3),
               quantile(a[, num], 0.7), Inf),
    labels = c("normoxia", "inter", "hypoxia")
  )
  
  
  a <- a %>%
    group_by(group, FLAD1) %>%
    count() %>%
    # as.data.frame() %>%
    # complete(group, FLAD1, fill = list(n = 0)) %>%
    pivot_wider(names_from = FLAD1, values_from = n, values_fill = 0) %>%
    as.data.frame() %>%
    filter(group != "inter")
  
  colnames(a)[2:4] <- c("del", "neutral", "amp")
  a$per_amp <- a$amp / (a$del + a$neutral + a$amp) * 100
  return(a)
}

plot_cnv_pan <- getcnvbar(a = FLAD1_cnv1, score = "Buffa_hypoxia_score_pan_cancer")
plot_cnv_pan

plot_cnv_LUAD <- getcnvbar(a = FLAD1_cnv1[FLAD1_cnv1$type == "LUAD",], 
                           score = "Buffa_hypoxia_score_pan_cancer")

plot_cnv_LUAD



getcnv_pval_fisher <- function(a){
  dt <- a[,1:4]
  
  dt$group <- as.character(dt$group)
  dt <- dt[order(dt$group),]
  
  dt$non_amp <- dt$del + dt$neutral
  dt <- dt[,c(5,4)]
  
  row.names(dt) <- c("hypoxia", "normoxia")
  dt <- as.matrix(dt)
  
  print(dt)
  return(fisher.test(dt)$p.value) 
}

getcnv_pval_fisher(a = plot_cnv_LUAD)
getcnv_pval_fisher(a = plot_cnv_pan)



plot_cnv_LUAD1 <- plot_cnv_LUAD[,1:4] %>% 
  summarise(across(where(is.numeric), sum)) %>%
  mutate(group = "total") %>%
  bind_rows(plot_cnv_LUAD[,1:4], .) %>%
  pivot_longer(cols = c(del:amp),
               names_to = "class",
               values_to = "number")

plot_cnv_LUAD1$group <- factor(plot_cnv_LUAD1$group, levels = c("total", "hypoxia", "normoxia"))
ggplot(plot_cnv_LUAD1, aes(x = group, y = number, fill = class)) +
  geom_bar(stat = "identity", position = "fill", width = .7) +
  ggtitle("LUAD FLAD1 cnv") +
  ylab("Percent (%)") +
  scale_y_continuous(breaks = c(seq(0,1,0.2))) +
  scale_fill_discrete(labels = c("amplification", "deletion", "neutral")) +
  theme_clean() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold", size = 12),
        axis.text = element_text(face = "bold", size = 10)) +
  annotate(geom = "text", x = 2.5, y = 1.1, label = "p value = 0.0013", size = 4)


plot_cnv_pan1 <- plot_cnv_pan[,1:4] %>% 
  summarise(across(where(is.numeric), sum)) %>%
  mutate(group = "total") %>%
  bind_rows(plot_cnv_pan[,1:4], .) %>%
  pivot_longer(cols = c(del:amp),
               names_to = "class",
               values_to = "number")

plot_cnv_pan1$group <- factor(plot_cnv_pan1$group, levels = c("total", "hypoxia", "normoxia"))
ggplot(plot_cnv_pan1, aes(x = group, y = number, fill = class)) +
  geom_bar(stat = "identity", position = "fill", width = 0.7) +
  ggtitle("pancancer FLAD1 CNV") +
  ylab("Percent (%)") +
  scale_y_continuous(breaks = c(seq(0,1,0.2))) +
  scale_fill_discrete(labels = c("amplification", "deletion", "neutral")) +
  theme_clean() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold", size = 12),
        axis.text = element_text(face = "bold", size = 10)) +
  annotate(geom = "text", x = 2.5, y = 1.1, label = "p value = 5.5e-21", size = 4)




#### fig S4 ####
Survival_plot <- function(Tumor, Gene){
  Sur <- fetch_dense_values(host = "https://gdc.xenahubs.net", 
                            dataset =  paste0("TCGA-", Tumor, ".htseq_fpkm.tsv"), 
                            identifiers = Gene,
                            use_probeMap = TRUE)
  Sur <- as.data.frame(t(Sur))
  Sur$sample <- row.names(Sur)
  
  data("XenaData")
  Sur_time <- XenaData %>% 
    dplyr::filter(XenaHostNames == "gdcHub") %>% 
    dplyr::filter(XenaDatasets == paste0("TCGA-", Tumor, ".survival.tsv")) %>%
    XenaGenerate() %>%
    XenaQuery() %>%
    XenaDownload() %>%
    XenaPrepare()
  
  Sur <- merge(Sur, Sur_time)
  Sur <- Sur %>% 
    dplyr::filter(substr(x = sample, start = 14, stop = 15) %in% c(paste0("0",1:9)) ) %>%
    dplyr::rename(time = OS.time,
                  status = OS)
  
  Sur$type <- Tumor
  
  return(Sur)
}


Sur <- Survival_plot(Tumor = "BLCA", Gene = c("FLAD1"))

Sur <- lapply(cancerType, function(x){
  Survival_plot(Tumor = x, Gene = c("FLAD1"))
})

Sur <- do.call("rbind", Sur)


Sur$class <- cut(x = Sur$FLAD1,
                 breaks = c(-Inf, quantile(Sur$FLAD1, 0.25), 
                            quantile(Sur$FLAD1, 0.75), Inf),
                 labels = c("Low", "mid", "High"))
Sur1 <- Sur[Sur$class != "mid",]

Sur1$class <- as.character(Sur1$class)


Sur1$months <- Sur1$time/30

Sur1 <- Sur1[Sur1$months <= 60,]
fit <- survfit(Surv(time = months, event = status) ~ class,
               data = Sur1)


ggsurvplot(fit = fit,
           conf.int = F,
           risk.table = T,
           pval = T,
           palette = c("red","blue"),
           surv.median.line = "hv",
           main = "Survival curves"
) + 
  theme_survminer(
    font.main = c(16, "bold"),
    font.x = c(14, "bold"),
    font.y = c(14, "bold"),
    font.tickslab = c(12, "bold"),
    legend = c(.85,.85),
    font.legend = c(12, "bold"),
  ) 



# HR森林图 
SurHR <- function(tumor){
  print(tumor)
  Sur <- Survival_plot(Tumor = tumor, Gene = c("FLAD1"))
  Sur$class <- cut(x = Sur$FLAD1,
                   breaks = c(-Inf, quantile(Sur$FLAD1, 0.25), 
                              quantile(Sur$FLAD1, 0.75), Inf),
                   labels = c("Low", "mid", "High"))
  Sur1 <- Sur[Sur$class != "mid",]
  Sur1$class <- as.character(Sur1$class)
  
  Sur1$months <- Sur1$time/30
  
  
  fit_cox <- coxph(Surv(time = months, event = status) ~ class,
                   data = Sur1)
  fit_cox
  
  forest_val <- summary(fit_cox)
  
  forest_hr <- forest_val$conf.int[,c(1,3,4)] %>%
    as.data.frame() %>%
    t()
  
  forest_p <- forest_val$coefficients[,c(2,5)] %>%
    as.data.frame() %>%
    t()
  
  forest_val <- merge(forest_hr, forest_p)
  
  forest_val$type <- tumor
  forest_val$number <- nrow(Sur1)
  forest_val$confint <- paste0("[",round(forest_val$`lower .95`,2)," ", 
                               round(forest_val$`upper .95`,2),"]")
  
  return(forest_val)
  
  
}


SurHR <- lapply(cancerType, function(x){
  SurHR(tumor = x)
})


SurHR1 <- do.call("rbind", SurHR)



SurHR1$pval <- ifelse(SurHR1$`Pr(>|z|)` < 0.001, "<0.001", 
                      round(SurHR1$`Pr(>|z|)`,3))
SurHR1$HR <- round(SurHR1$`exp(coef)`,2)
forestplot(SurHR1[,c(5,6,9,7,8)],
           mean = SurHR1$`exp(coef)`,
           lower = SurHR1$`lower .95`,
           upper = SurHR1$`upper .95`,
           zero = 1,
           boxsize = 0.3,
           graph.pos = 3,
           clip = c(0,2),
           hrzl_lines=list("1" = gpar(lwd=2, col="black"),
                           "21" = gpar(lwd=2, col="black")),
           xlab = "x"
)


#### fig S6 ####

# * 1. normal vs tumor ----
list_DEG_NT <- list()
for (i in 1:length(names(list_counts_all))) {
  
  counts <- list_counts_all[[i]]
  tumorSample <- colnames(counts)[which(substr(colnames(counts),14,15) %in% c(paste0("0",1:9)))]
  normalSample <- colnames(counts)[which(substr(colnames(counts),14,15) %in% c(paste0("1",0:9)))]
  
  if (length(normalSample) == 0) {
    next
  }
  
  DEG <- do_DEG(counts = counts, 
                groupA = normalSample, groupA_label = "normal",
                groupB = tumorSample, groupB_label = "tumor")
  
  DEG1 <- as.data.frame(DEG)
  
  DEG1$id <- row.names(DEG1)
  DEG1 <- merge(x = probe, y = DEG1, by = "id")
  DEG1$type <- names(list_counts_all)[i]
  
  list_DEG_NT[[(length(list_DEG_NT)+1)]] <- DEG1
  
}



# * 2. hypoxia vs normoxia ----

# 肿瘤类型
cancerType_HN <- unique(hyScore$tumour_type)

getsample <- function(hyScore, tumor, group, score, id){
  dt <- hyScore[hyScore$tumour_type == tumor,]
  
  score <- grep(score, colnames(hyScore))
  dt <- dt[,c(1,score,length(colnames(hyScore)))]
  
  dt$class <- cut(dt[,score], 
                  breaks = c(-Inf, quantile(dt[,score], 0.3), quantile(dt[,score], 0.7), Inf),
                  labels = c("normoxia", "inter", "hypoxia")
  )
  if(group == "hypoxia"){
    id_sim <- dt[dt$class == "hypoxia",]$patient_id
  }else{
    id_sim <- dt[dt$class == "normoxia",]$patient_id
  }
  id_sim  <- gsub("[.]", "-",id_sim)
  
  id_sim  <- lapply(id_sim, function(x){
    id[str_detect(id, x)]
  }) %>% 
    unlist()

  id_sim  <- id_sim[substr(id_sim,14,15)  %in% c(paste0("0",1:9)) & substr(id_sim,16,16) == "A"]
  
  return(id_sim)
  
  
}



list_DEG_HN <- list()

for (i in 1:length(cancerType_HN)) {
  
  if(cancerType_HN[i] == "COADREAD"){
    counts <-
      merge(list_counts_all[[which(names(list_counts_all) == "COAD")]],
            list_counts_all[[which(names(list_counts_all) == "READ")]],
            by = "Ensembl_ID")
  }else{
    counts <- list_counts_all[[grep(cancerType_HN[i], names(list_counts_all))]]
  }
  
  
  sample_hypoxia <- getsample(hyScore = hyScore, tumor = cancerType_HN[i], group = "hypoxia", 
                              score = "Buffa_hypoxia_score_pan_cancer",
                              id = colnames(counts) )
  
  sample_normoxia <- getsample(hyScore = hyScore, tumor = cancerType_HN[i], group = "normoxia", 
                               score = "Buffa_hypoxia_score_pan_cancer",
                               id = colnames(counts) )
  
  if (length(sample_hypoxia) == 0) {
    next
  }
  
  if (length(sample_normoxia) == 0) {
    next
  }
  
  DEG <- do_DEG(counts = counts, 
                groupA = sample_normoxia, groupA_label = "normoxia",
                groupB = sample_hypoxia, groupB_label = "hypoxia")
  
  print(head(DEG, 0)) # condition hypoxia vs normoxia
  
  DEG1 <- as.data.frame(DEG)
  
  DEG1$id <- row.names(DEG1)
  DEG1 <- merge(x = probe, y = DEG1, by = "id")
  DEG1$type <- cancerType_HN[i]
  
  list_DEG_HN[[(length(list_DEG_HN)+1)]] <- DEG1
  
}

# 绘图
pheatmap(plot_NT_ribo_reg,
         scale = "none",
         legend = F,
         cluster_rows = F, 
         cluster_cols = F,
         display_numbers = round(plot_NT_ribo_log2FC, 2),
         cellwidth = 30,
         cellheight = 30,
         color = colorRampPalette(c("#41b6e6", "white","#ff585d"))(3),
         border_color = "black",
         main = "Tumor VS. Normol",
         angle_col = 90,
         fontsize_row = 12, fontsize_col = 12, fontsize_number = 11)



pheatmap(plot_HN_ribo_reg,
         scale = "none",
         legend = F,
         cluster_rows = F, 
         cluster_cols = F,
         display_numbers = round(plot_HN_ribo_log2FC, 2),
         cellwidth = 30,
         cellheight = 30,
         color = colorRampPalette(c("#41b6e6", "white","#ff585d"))(3),
         border_color = "black",
         main = "Hypoxia VS. Normoxia",
         angle_col = 90,
         fontsize_row = 12, fontsize_col = 12, fontsize_number = 11)



ggboxplot(plot_LUAD_fpkm, x = "gene", y = "val",
          fill = "group_NT",
          palette = c("#00AFBB", "#FC4E07"),
          #add = "median_q1q3",
          outlier.shape = NA,
          ylim=c(0,10),
          ylab = "log2(FPKM + 1)")+
  stat_compare_means(aes(group=group_NT), label = "p.signif", 
                     label.y = 9, size=10,
                     symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, Inf), 
                                        symbols = c("***", "**", "*", "ns"))
  )+
  scale_y_continuous(breaks = c(seq(0,10,2)))+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 24, face = "bold"),
        axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 24, face = "bold"),
        axis.text.y = element_text(size = 24, face = "bold"),
        legend.position = "right",
        legend.text = element_text(face = "plain", size = 24),
        legend.title = element_text(face = "plain", size = 24))+
  ggtitle("Riboflavin Metabolism (tumor vs normal)")+
  labs(fill ="Type")



ggboxplot(plot_LUAD_fpkm[plot_LUAD_fpkm$group_HN != "inter",], x = "gene", y = "val",
          fill = "group_HN",
          palette = c("#00AFBB", "#FC4E07"),
          outlier.shape = NA,
          ylim=c(0,10),
          ylab = "log2(FPKM + 1)")+
  stat_compare_means(aes(group=group_HN), label = "p.signif", 
                     label.y = 9, size=10,
                     symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, Inf), 
                                        symbols = c("***", "**", "*", "ns"))
  )+
  scale_y_continuous(breaks = c(seq(0,10,2)))+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 24, face = "bold"),
        axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 24, face = "bold"),
        axis.text.y = element_text(size = 24, face = "bold"),
        legend.position = "right",
        legend.text = element_text(face = "plain", size = 24),
        legend.title = element_text(face = "plain", size = 24))+
  ggtitle("Riboflavin Metabolism (hypoxia vs. normoxia)")+
  labs(fill ="Type")

# * 3.缺氧程度与核黄素代谢水平 ----
calribo <- function(tumor, gene){
  exp <- fetchfpkm(tumor = tumor, gene = gene)#log2(fpkm+1)
  exp <- subset(exp, select= -type)
  exp <- 2^exp - 1
  exp <- scale(exp)
  
  exp <- data.frame(sample = row.names(exp), mean = apply(exp, 1, mean))
  exp$type <- tumor
  
  return(exp)
}


getsample_hyscore <- function(tumor, hyScore, score, dt_exp){
  dt <- hyScore[hyScore$tumour_type == tumor,]
  score <- grep(score, colnames(hyScore))
  dt <- dt[,c(1,score)]
  colnames(dt)[2] <- "hyScore"
  
  dt_exp <- dt_exp[substr(dt_exp$sample,14,15) %in% c(paste0("0",1:9)) & substr(dt_exp$sample,16,16) == "A",] 
  dt_exp$id_sim <- gsub("-",".",substr(x = dt_exp$sample,1,12))
  
  dt_exp <- merge(dt_exp, dt, by.x = "id_sim", by.y = "patient_id")
  
  return(dt_exp)
  
}

list_ribo_hypoxia <- lapply(cancerType_HN[-4], function(tumor){
  
  dt_exp <- calribo(tumor = tumor,
                    gene = ribogene$gene)
  
  dt_exp_score <- getsample_hyscore(tumor = tumor, 
                                    hyScore = hyScore,
                                    score = "Buffa_hypoxia_score_pan_cancer",
                                    dt_exp = dt_exp)
  
  dt_exp_score$class <- cut(dt_exp_score$hyScore, breaks = 5, labels = 1:5)
  dt_exp_score <- dt_exp_score %>%
    dplyr::select(mean, type, class) %>%
    dplyr::group_by(class, type) %>%
    dplyr::summarise(class_mean = mean(mean),.groups = "drop") %>%
    as.data.frame()
  
  print(tumor)
  return(dt_exp_score)
})

plot_ribo_hypoxia <- do.call("rbind", list_ribo_hypoxia)
plot_ribo_hypoxia <- plot_ribo_hypoxia %>%
  pivot_wider(names_from = "class",
              values_from = "class_mean") %>%
  as.data.frame()
row.names(plot_ribo_hypoxia) <- plot_ribo_hypoxia$type
plot_ribo_hypoxia <- plot_ribo_hypoxia[,-1]
colnames(plot_ribo_hypoxia) <- c("Very low", "Low", "Moderate", "High", "Very high")

pheatmap(plot_ribo_hypoxia,
         cluster_cols = F,
         cluster_rows = F,
         scale = "row",angle_col = 45, fontsize = 20)

#### fig S7 ####
# * 1. FLAD1 high vs FLAD1 low ----
fetchfpkm <- function(tumor, gene){
  dataset <- paste0("TCGA-", tumor, ".htseq_fpkm.tsv")
  exp <- fetch_dense_values(host = "https://gdc.xenahubs.net", 
                            dataset =  dataset, 
                            identifiers = gene,
                            use_probeMap = TRUE)
  exp <- as.data.frame(t(exp))
  exp$type <- tumor
  
  return(exp)
} 

getsample_FLAD1 <- function(tumor, class){
  dt_class <- fetchfpkm(tumor = tumor, gene = "FLAD1")
  
  dt_class$class <- cut(dt_class$FLAD1, 
                        breaks = c(-Inf, quantile(dt_class$FLAD1, 0.5), Inf),
                        labels = c("low", "high"))
  return(row.names(dt_class[dt_class$class == class,]))
}


list_DEG_FLAD1_HL <- list()
for (i in 1:length(names(list_counts_all))) {
  
  counts <- list_counts_all[[i]]
  
  sample_high <- getsample_FLAD1(tumor = names(list_counts_all)[i], class = "high")
  sample_low <- getsample_FLAD1(tumor = names(list_counts_all)[i], class = "low")
  
  
  if (length(sample_high) == 0) {
    next
  }
  
  if (length(sample_low) == 0) {
    next
  }
  
  DEG <- do_DEG(counts = counts, 
                groupA = sample_low, groupA_label = "low",
                groupB = sample_high, groupB_label = "high")
  
  print(head(DEG, 0)) 
  
  DEG1 <- as.data.frame(DEG)
  
  DEG1$id <- row.names(DEG1)
  DEG1 <- merge(x = probe, y = DEG1, by = "id")
  DEG1$type <- names(list_counts_all)[i]
  
  list_DEG_FLAD1_HL[[(length(list_DEG_FLAD1_HL)+1)]] <- DEG1
  
}

# * 2. GO分析FLAD1高表达组 ----
# LUAD
GO_LUAD_FLAD1_HL <- list_DEG_FLAD1_HL[[5]]
GO_LUAD_FLAD1_HL <- GO_LUAD_FLAD1_HL %>%
  dplyr::filter(log2FoldChange > 0.5) %>%
  dplyr::filter(padj < 0.05)

GO_LUAD_FLAD1_HL <- enrichGO(gene = GO_LUAD_FLAD1_HL$gene,
                             OrgDb = "org.Hs.eg.db", keyType = "SYMBOL",
                             ont = "CC", readable = T)

GO_LUAD_FLAD1_HL <- data.frame(GO_LUAD_FLAD1_HL)
GO_LUAD_FLAD1_HL <- GO_LUAD_FLAD1_HL[order(GO_LUAD_FLAD1_HL$p.adjust),]
GO_LUAD_FLAD1_HL <- GO_LUAD_FLAD1_HL[1:10,]
GO_LUAD_FLAD1_HL$Description <- str_to_title(GO_LUAD_FLAD1_HL$Description)

ggbarplot(GO_LUAD_FLAD1_HL, x = "Count", y = "Description", 
          fill = "p.adjust",
          ylab = "") +
  theme_clean() +
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 20))



#### fig S9 ####
# * 1. heatmap ----
pathways <- gmtPathways(gmt.file = "~/TCGA_data/GDCdata/c2.all.v2023.2.Hs.symbols.gmt") # 下载基因集：https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2023.2.Hs/c2.all.v2023.2.Hs.symbols.gmt

list_fgseaRes <- lapply(1:length(list_DEG_FLAD1_HL), function(x){
  DEG <- list_DEG_FLAD1_HL[[x]]
  
  ranks <- DEG %>%
    dplyr::select(gene, stat) %>%
    arrange(desc(stat))
  ranks <- deframe(ranks)
  ranks <- na.omit(ranks)
  
  fgseaRes <- fgsea(pathways = pathways, 
                    stats    = ranks,
                    minSize  = 5,
                    maxSize  = Inf,
                    eps = 0)
  fgseaRes$type <- DEG$type[1]
  return(fgseaRes)
})


kegg_met_path <- read.csv(file = "~/FLAD1/fig_data/KEGG_path1.csv",header = T) # KEGG gene set
kegg_met_path <- kegg_met_path$pathways[kegg_met_path$met == 1]

fgseaRes_met <- do.call("rbind", list_fgseaRes)
fgseaRes_met <- fgseaRes_met[fgseaRes_met$pathway %in% kegg_met_path,]

plot_fgseaRes_met_NES <- fgseaRes_met %>%
  dplyr::mutate(path = str_to_title(gsub("_"," ",gsub("KEGG_", "", fgseaRes_met$pathway)))) %>%
  dplyr::select(path, NES, type) %>%
  pivot_wider(names_from = type, values_from = NES) %>%
  as.data.frame()
row.names(plot_fgseaRes_met_NES) <- plot_fgseaRes_met_NES$path
plot_fgseaRes_met_NES <- plot_fgseaRes_met_NES[,2:ncol(plot_fgseaRes_met_NES)]

plot_fgseaRes_met_padj <- fgseaRes_met %>%
  dplyr::mutate(path = str_to_title(gsub("_"," ",gsub("KEGG_", "", fgseaRes_met$pathway)))) %>%
  dplyr::select(path, padj, type) %>%
  pivot_wider(names_from = type, values_from = padj) %>%
  as.data.frame()
row.names(plot_fgseaRes_met_padj) <- plot_fgseaRes_met_padj$path
plot_fgseaRes_met_padj <- plot_fgseaRes_met_padj[,2:ncol(plot_fgseaRes_met_padj)]

phe_num_met <- apply(plot_fgseaRes_met_NES, 1, mean)
phe_num_met <- data.frame(met = names(phe_num_met),
                          NES_mean = as.numeric(phe_num_met))
phe_num_met <- phe_num_met[order(-phe_num_met$NES_mean),]


plot_fgseaRes_met_NES <- plot_fgseaRes_met_NES[match(phe_num_met$met, row.names(plot_fgseaRes_met_NES) ),]
plot_fgseaRes_met_padj <- plot_fgseaRes_met_padj[match(phe_num_met$met, row.names(plot_fgseaRes_met_padj)),]
plot_fgseaRes_met_padj <- plot_fgseaRes_met_padj[1:20,] 

pheatmap(plot_fgseaRes_met_NES[1:20,],
         cluster_rows = F,
         cluster_cols = F,
         display_numbers = matrix(ifelse(plot_fgseaRes_met_padj < 0.05, "*", ""), nrow = nrow(plot_fgseaRes_met_padj)),
         fontsize_number = 13,
         angle_col = 45)

# * 2. 柱状图 氧化磷酸化 ----
gene_oxphos <- pathways[[which(names(pathways) == "KEGG_OXIDATIVE_PHOSPHORYLATION")]]
gene_oxphos <- gene_oxphos[gene_oxphos %in% probe$gene]

caloxphos <- function(tumor){
  exp <- fetchfpkm(tumor = tumor, gene = c(gene_oxphos, "FLAD1")) #log2(fpkm+1)
  exp <- subset(exp, select= -type)
  exp <- 2^exp - 1
  exp <- scale(exp)
  
  exp1 <- subset(exp, select = -FLAD1)
  exp1 <- data.frame(sample = row.names(exp1), mean = apply(exp1, 1, mean))
  
  exp <- subset(exp, select = FLAD1) 
  exp <- as.data.frame(exp)
  exp$sample <- row.names(exp)
  
  exp1 <- merge(exp1, exp, by = "sample")
  
  sample_high <- row.names(exp)[exp$FLAD1 > median(exp$FLAD1)]
  sample_low <- row.names(exp)[exp$FLAD1 <= median(exp$FLAD1)]
  
  exp1$class <- ifelse(exp1$sample %in% sample_high, "high",
                       ifelse(exp1$sample %in% sample_low, "low", "none"))
  exp1$type <- tumor
  
  return(exp1)
  
}


list_oxphos_HL <- lapply(CancerType, function(x){
  caloxphos(tumor = x)
})

list_oxphos_HL1 <- do.call("rbind", list_oxphos_HL)
list_oxphos_HL1 <- list_oxphos_HL1[list_oxphos_HL1$type != "UVM",]
list_oxphos_HL1$class <- factor(list_oxphos_HL1$class, levels = c("low", "high"))
list_oxphos_HL1 <- list_oxphos_HL1[order(list_oxphos_HL1$type),]

list_oxphos_HL_n <- list_oxphos_HL1 %>%
  dplyr::select(class, type) %>%
  group_by(type, class) %>%
  dplyr::count()


ggboxplot(list_oxphos_HL1, x = "type", y = "mean", 
          fill = "class", 
          palette = c( "#00AFBB","#FC4E07"),
          ylim = c(-1.5, 2),outlier.shape = NA,
          ylab = "Average gene expression level \n (normalized)") +
  stat_compare_means(aes(group=class), label = "p.signif", label.y = 1.8, 
                     size = 6, angle = 90,
                     symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, Inf), 
                                        symbols = c("***", "**", "*", "ns")))+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, face = "bold", size = 14),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 12, face = "plain"),
        legend.position = "right",
        legend.text = element_text(face = "plain", size = 10),
        legend.title = element_text(face = "plain", size = 14))+
  annotate(geom = "text", y=-1.15, x=1:32, label = list_oxphos_HL_n[list_oxphos_HL_n$class == "high",]$n)+
  annotate(geom = "text", y=-1.35, x=1:32, label = list_oxphos_HL_n[list_oxphos_HL_n$class == "low",]$n)

#* 3. LUAD OXPHOS GSEA ----
ranks <- list_DEG_FLAD1_HL[[grep("LUAD", CancerType)]] %>%
  dplyr::select(gene, stat) %>%
  arrange(desc(stat))
ranks <- deframe(ranks)
ranks <- na.omit(ranks)

plotEnrichment(pathways[["KEGG_OXIDATIVE_PHOSPHORYLATION"]],
               ranks) +  
  labs(title= "KEGG Oxidative Phosphorylation (LUAD)") +
  theme_bw() +
  ylab("Enrichment Score") +
  xlab("Rank in Ordered Dataset") +
  annotate(geom = "text", x = 50000, y = 0.6,
           label = "NES = 2.71") +
  annotate(geom = "text", x = 50000, y = 0.5,
           label = "p adj = 1.117e-17") +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 20),
        title = element_text(size = 20))

