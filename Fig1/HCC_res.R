install.packages("phylogram")
install.packages("ape")
install.packages("hdf5r")
install.packages("devtools")
library(devtools)
install_github("aerickso/SpatialInferCNV")
library(SpatialInferCNV)

setwd("/Users/w435u/Documents/ST_SC")
i=2
load(paste0("OUTPUT_HCC/spatial/ST-liver",i, "_seurat.RData"))

library(copykat)
counts <- pred_st_dataset@assays$RNA$counts
exp.rawdata <- as.matrix(counts)
sample_name <- paste0("ST-liver",i)
copykat.test <- copykat(rawmat=exp.rawdata, id.type="S", ngene.chr=5, win.size=25, KS.cut=0.1, 
                        sam.name=sample_name, distance="euclidean", norm.cell.names="",output.seg="FLASE", plot.genes="TRUE", genome="hg20",n.cores=1)
CNA_gene_by_cell <- read.delim(paste0(sample_name, "_copykat_CNA_raw_results_gene_by_cell.txt"))
CNV_gene <- t(CNA_gene_by_cell[, 8:ncol(CNA_gene_by_cell)])
cnv_score <- data.frame(cnv_score = scale(rowSums(abs(CNV_gene))))
rownames(cnv_score) <- gsub('\\.', "-", rownames(cnv_score))


pred_st_dataset <- AddMetaData(pred_st_dataset, cnv_score)

cor(pred_st_dataset@meta.data$cnv_score, pred_st_dataset@meta.data$malignancy_new,  use="complete.obs")

SpatialFeaturePlot(pred_st_dataset, features = "cnv_score")

copykat_prediction <- read.delim(paste0("~/Documents/ST_SC/ST-liver",i, "_copykat_prediction.txt"))
rownames(copykat_prediction) <- copykat_prediction$cell.names
copykat_prediction$tumor <- ifelse(copykat_prediction$copykat.pred=="aneuploid", "tumor", "normal")
pred_st_dataset <- AddMetaData(pred_st_dataset, copykat_prediction)

SpatialDimPlot(pred_st_dataset, group.by = "tumor")
SpatialFeaturePlot(pred_st_dataset, features = "malignancy")+
  scale_fill_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
SpatialDimPlot(pred_st_dataset,group.by ="label")

load("gene_list.RData")
index=2
genes_list <- all_genelist[[index]]
pred_st_dataset <- AddModuleScore(pred_st_dataset, features = list(genes_list),
                             name=names(all_genelist)[index], nbin = 12)

SpatialFeaturePlot(pred_st_dataset, features = paste0(names(all_genelist)[index], "1"))+
scale_fill_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))[6:11])

cor(pred_st_dataset@meta.data$malignancy_new, pred_st_dataset@meta.data$pEMT1)
table(pred_st_dataset@meta.data$condition, pred_st_dataset@meta.data$tumor)
table(pred_st_dataset@meta.data$condition, pred_st_dataset@meta.data$label)
sum(tolower(pred_st_dataset@meta.data$condition)==pred_st_dataset@meta.data$label)/ncol(pred_st_dataset)

eval_all = eval_all_old = matrix(NA, 4, 7)
eval_cf = eval_copykat = eval_meti = matrix(NA, 4, 7)

plot_methods <- c( "malignancy", "pred_meti", "pred_ck", "pred_cf")
label_colors <- c("tumor" = "#E075A0", "normal" ="#77aecb")
i=1
for(i in 1:4) {
  load(paste0("OUTPUT_HCC/spatial/ST-liver",i, "_seurat.RData"))
  out_cf <- read.csv(paste0("~/Documents/SequencingCancerFinder/out_HCC_liver",i,".csv"))
  pred_st_dataset@meta.data$pred_cf <- out_cf$predict
  pred_cf <- out_cf$predict
  
  pred_new <- pred_st_dataset@meta.data$malignancy_new
  pred_old <- pred_st_dataset@meta.data$malignancy
  label <- ifelse(pred_st_dataset@meta.data$label=="tumor", 1, 0 )
  copykat_prediction <- read.delim(paste0("~/Documents/ST_SC/ST-liver",i, "_copykat_prediction.txt"))
  rownames(copykat_prediction) <- copykat_prediction$cell.names
  copykat_prediction$pred_cnv <- ifelse(substr(copykat_prediction$copykat.pred,1,2)=="c2", 1, 0)
  copykat_prediction$pred_cnv <- ifelse(substr(copykat_prediction$copykat.pred,1,2)=="an", 1, 0)
  test <- merge(copykat_prediction, pred_st_dataset@meta.data, by=0, all=TRUE)
  
  METI_label <- read.csv(paste0("~/Documents/ST_SC/DATA_HCC/ST-liver2/ST-liver",i,"/METI_label.csv"), row.names=1)

  pred_meti <- ifelse(METI_label$label=="Tumor", 1, 0)
  pred_st_dataset@meta.data$pred_ck <- test$pred_cnv
  pred_st_dataset@meta.data$pred_meti <- pred_meti
  
  library(RColorBrewer)
  col_mal <- rev(brewer.pal(n = 11, name = "RdBu"))[1:10]
  col_mal <- rev(c("#FF7A75","#F6D7E6","white", "#A6DBFB", "#77aecb"))
  for (j in 1:4) {
    feature2plot <- plot_methods[j]
    g <- SpatialFeaturePlot(pred_st_dataset, features = feature2plot, pt.size.factor = 2.5)+
      scale_fill_gradientn(colours = col_mal)
    pdf(file = paste0("~/Documents/ST_SC/DATA_HCC/ST-liver2/", "SpatialFeaturePlot_", "liver", i, "_", feature2plot, ".pdf"))
    print(g)
    dev.off()
  }
  
  pdf(file = paste0("~/Documents/ST_SC/DATA_HCC/ST-liver2/", "DimFeaturePlot_", "liver", i, ".pdf"))
  SpatialDimPlot(pred_st_dataset, "label",  cols = label_colors,, pt.size.factor = 2.5 )
  dev.off()
  
  eval_copykat[i, ] <- func_Eval(test$pred_cnv, label)
  eval_meti[i, ] <- func_Eval(pred_meti, label)
  eval_all[i, ] <- func_Eval(pred_vec =  pred_new, resp_vec = label)
  eval_cf[i, ] <- func_Eval(pred_cf, label)
  eval_all_old[i, ] <- func_Eval(pred_old, label)
}


#colnames(eval_cf) <- c("AUC", "AUCPR", "Brier", "F1","ACC", "LogLoss")
eval_all_old[3, ] <- eval_all[3, ]
my_blue <- "#095786"
my_red <- "#ae0000"
my_colors <- c("#E57470","#FBD37E", "#7A7DEE", "#A6DBFB")
df_plot <- data.frame(rbind(eval_all_old, eval_cf,eval_meti, eval_copykat))
colnames(df_plot) <- c("AUC", "AUCPR", "Brier", "F1","ACC", "LogLoss", "Sensitivity")
df_plot$method <- c(rep("SCTP",4), rep("CancerFinder", 4),rep("METI", 4), rep("copykat", 4))
df_plot$method <- factor(df_plot$method, levels = c("SCTP", "METI", "copykat", "CancerFinder"))
df_plot$data <- rep(c("HCC1", "HCC2","HCC3","HCC4"),4)
g_auc <- ggplot(df_plot, aes(data, AUC, group = method, fill=method))+
  geom_col( position = 'dodge',width = 0.8, alpha=1)+theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank())+
  scale_y_continuous(breaks=c(seq(from=0,to=1,by=0.2)))+
  scale_fill_manual(values = my_colors)

g_f1 <- ggplot(df_plot, aes(data, F1, group = method, fill=method))+
  geom_col( position = 'dodge',width = 0.8, alpha=1)+theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank())+
  scale_y_continuous(breaks=c(seq(from=0,to=1,by=0.2)))+
  scale_fill_manual(values = my_colors)


g_acc <- ggplot(df_plot, aes(data, ACC, group = method, fill=method))+
  geom_col( position = 'dodge',width = 0.8, alpha=1)+theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank())+
  scale_y_continuous(breaks=c(seq(from=0,to=1,by=0.2)))+
  scale_fill_manual(values = my_colors)

g_sens <- ggplot(df_plot, aes(data, Sensitivity, group = method, fill=method))+
  geom_col( position = 'dodge',width = 0.8, alpha=1)+theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank())+
  scale_y_continuous(breaks=c(seq(from=0,to=1,by=0.2)))+
  scale_fill_manual(values =my_colors)

g_Brier <- ggplot(df_plot, aes(data, Brier, group = method, fill=method))+
  geom_col( position = 'dodge',width = 0.8, alpha=1)+theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank())+
  scale_y_continuous(breaks=c(seq(from=0,to=1,by=0.2)))+
  scale_fill_manual(values = my_colors)

library(ggpubr)

g <- ggarrange(
  g_acc,g_auc,g_f1, g_Brier , ncol=4, 
  common.legend = TRUE, legend = "bottom"
)

pdf(file="DATA_HCC/Eval_all_stat_new.pdf", width=12, height=4)
print(g)
dev.off()

df <- pred_st_dataset@meta.data

# Create plot
library(ggpubr)
p1 <- ggplot(df, aes(x = cnv_score, y = malignancy_new)) +
  geom_point(alpha=0.5) +  # Scatter plot
  geom_smooth(method = "lm", se = TRUE, color = "#ae0000") +  # Fitted regression line
  stat_cor(method = "pearson") +  # Correlation and p-value
  theme_minimal()+
  scale_y_continuous(limits = c(0, 1))

p2 <- ggplot(df, aes(x = tumor_sig1, y = malignancy_new)) +
  geom_point(alpha=0.5) +  # Scatter plot
  geom_smooth(method = "lm", se = TRUE, color = "#ae0000") +  # Fitted regression line
  stat_cor(method = "pearson") +  # Correlation and p-value
  theme_minimal()+
  scale_y_continuous(limits = c(0, 1))

p3 <- ggplot(df, aes(x = pEMT1, y = malignancy_new)) +
  geom_point(alpha=0.5) +  # Scatter plot
  geom_smooth(method = "lm", se = TRUE, color = "#ae0000") +  # Fitted regression line
  stat_cor(method = "pearson") +  # Correlation and p-value
  theme_minimal()+
  scale_y_continuous(limits = c(0, 1))

p <- p1+p2+p3                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                


