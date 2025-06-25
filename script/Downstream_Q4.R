setwd("/Users/w435u/Documents/ST_SC")
load("OUTPUT_CAO_new/Seurat_C4_subtumor.RData")
work_dir <- "OUTPUT_CAO_new"
name_st <- "C4"
library(presto)
library(ComplexHeatmap)
library(circlize)
library(Seurat)
library(dplyr)
library(tidyverse)
library(ggrepel)

markers<- presto::wilcoxauc(st_dataset, 'T12', assay = 'data')
markers<- top_markers(markers, n = 10, auc_min = 0.5, pct_in_min = 20, pct_out_max = 20)

markers <- markers[which(markers$padj<0.05), ]

tumor_enriched <- FindMarkers(st_dataset, ident.1 = "Positive", group.by = 'SCTP', verbose = FALSE)


out <- split( markers[,c("feature","logFC") ] , f = markers$group )


all_markers<- markers %>%
  dplyr::select(-rank) %>% 
  unclass() %>% 
  stack() %>%
  pull(values) %>%
  unique() %>%
  .[!is.na(.)]

out_name = paste0("OUTPUT_CAO_new/FIG/SpatialDimPlot_tj_",name_st, ".pdf")
pdf(out_name)
SpatialDimPlot(st_dataset, group.by = "T12")
dev.off()




tumor_enriched1 <- FindMarkers(st_dataset, ident.1 = "T2",ident.2 = "T1", group.by = 'T12', verbose = FALSE)
tumor_enriched2 <- FindMarkers(st_dataset, ident.1 = "T1",ident.2 = "Other", group.by = 'T12', verbose = FALSE)
tumor_enriched3 <- FindMarkers(st_dataset, ident.1 = "T2",ident.2 = "Other", group.by = 'T12', verbose = FALSE)

write.xlsx(tumor_enriched1, file="FIGURES/DEGs_FIG5E1.xlsx", sheetName="2 v.s. 3")
write.xlsx(tumor_enriched2, file="FIGURES/DEGs_FIG5E2.xlsx", sheetName="1 v.s. 3")
write.xlsx(tumor_enriched3, file="FIGURES/DEGs_FIG5E3.xlsx", sheetName="1 v.s. 2")


subgenes1 <- tumor_enriched1%>% filter(p_val<10e-14)%>% filter(abs(avg_log2FC)>1.5)
subgenes1$Group <- "T2 v.s T1"
subgenes2 <- tumor_enriched2%>% filter(p_val<10e-14)%>% filter(abs(avg_log2FC)>1.5)
subgenes2$Group <- "T1 v.s Other"
subgenes3 <- tumor_enriched3%>% filter(p_val<10e-14)%>% filter(abs(avg_log2FC)>1.5)
subgenes3$Group <- "T2 v.s Other"

suball <- rbind(subgenes1, subgenes2, subgenes3)
suball$gene <- rownames(suball)
suball$label <- ifelse(suball$avg_log2FC < (-0.5), '#0068B2',
                       ifelse(suball$avg_log2FC > 0.5, '#B2182B',
                              'black'))
dbar <- suball %>% group_by(Group) %>% summarise_all(list(min=min, max=max))%>%
  dplyr::select(Group, avg_log2FC_min, avg_log2FC_max)

my_col3 = c( "black","#e60012", "#00a0e9")


g <- ggplot()+
  geom_col(data = dbar,  # 绘制负向背景柱状图
           mapping = aes(x = Group,y = avg_log2FC_min),
           fill = "#dcdcdc",alpha = 0.6, width = 0.7) +
  geom_col(data = dbar, # 绘制正向背景柱状图
           mapping = aes(x = Group,y = avg_log2FC_max),
           fill = "#dcdcdc",alpha = 0.6, width = 0.7) +
  geom_jitter(data = suball, # 绘制所有数据点
              aes(x = Group, y = avg_log2FC, color = label),
              size = 0.85,
              width =0.3) +
  geom_tile(data = suball, # 绘制中心分组标记图
            aes(x = Group,
                y = 0,
                fill = Group),
            height=0.5,
            color = "black",
            alpha = 0.4,
            show.legend = F) +
  ggsci::scale_fill_npg(alpha=0.5) + # 自定义颜色
  scale_color_manual(values=rev(my_col3))+# 自定义颜色
  geom_text_repel(data = suball,  # 这里的filter很关键，筛选你想要标记的基因
                  aes(x = Group, y = avg_log2FC, label = gene),
                  size = 2, 
                  max.overlaps = getOption("ggrepel.max.overlaps", default = 15),
                  color = 'black',
                  force = 1.2,
                  arrow = arrow(length = unit(0.008, "npc"),
                                type = "open", ends = "last"))+
  labs(x="Cluster", y="Average logFC") +
  geom_text(data=suball, # 绘制中心分组标记图文本注释
            aes(x=Group, 
                y=0, 
                label=Group),
            size = 3,
            color ="white") +
  theme_minimal() + 
  theme(axis.title = element_text(size = 13,color = "black",face = "bold"),
        axis.line.y = element_line(color = "black",size = 1.2),
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        panel.grid = element_blank(),
        legend.position = "top",
        legend.direction = "vertical",
        legend.justification = c(1,0),
        legend.text = element_text(size = 13))


genes <- unique(c(rownames(tumor_enriched1)[1:25],rownames(tumor_enriched2)[1:25]))

for (i in 1:length(genes)) {
  gene = genes[i]
  gg <- SpatialFeaturePlot(st_dataset, features = gene)
  
  pdf(file=paste0("OUTPUT_CAO_new/FIG/FeaturePlot_", name_st, "_",gene, ".pdf"))
  print(gg)
  dev.off()
}

p <- DotPlot(object = st_dataset, features =genes, group.by = "T12")


df <- p$data
#df$features.plot <- rownames(df)
exp_mat<-df %>% 
  dplyr::select(-pct.exp, -avg.exp) %>%  
  pivot_wider(names_from = id, values_from = avg.exp.scaled) %>% 
  as.data.frame() 

exp_mat <- exp_mat[!is.na(exp_mat$features.plot), ]

row.names(exp_mat) <- exp_mat$features.plot  
exp_mat <- exp_mat[,-1] %>% data.matrix()

percent_mat<-df %>% 
  dplyr::select(-avg.exp, -avg.exp.scaled) %>%  
  pivot_wider(names_from = id, values_from = pct.exp) %>% 
  as.data.frame() 

percent_mat <- percent_mat[!is.na(percent_mat$features.plot), ]


row.names(percent_mat) <- percent_mat$features.plot  
percent_mat <- percent_mat[,-1] %>% data.matrix()

library(viridis)
col_fun = circlize::colorRamp2(c(-1, 0, 2), plasma(20)[c(1,10, 20)])
cell_fun = function(j, i, x, y, w, h, fill){
  grid.rect(x = x, y = y, width = w, height = h, 
            gp = gpar(col = NA, fill = NA))
  grid.circle(x=x,y=y,r= percent_mat[i, j]/100 * min(unit.c(w, h)),
              gp = gpar(fill = col_fun(exp_mat[i, j]), col = NA))}

write.csv(exp_mat, file="FIGURES/Expression_FIG5F.csv")


dot <- Heatmap(exp_mat,
               heatmap_legend_param=list(title="expression"),
               column_title = "clustered dotplot", 
               col=col_fun,
               rect_gp = gpar(type = "none"),
               cell_fun = cell_fun,
               row_names_gp = gpar(fontsize = 10),
               row_km = 3,
               border = "black",
               cluster_columns = FALSE,column_order = c("T2","T1","Other"))

dot

pdf("OUTPUT_CAO_new/FIG/Dotplot_C4_all_new.pdf")
dot
dev.off()

###### TCGA CRC Stages ---------
library(TCGAcrcmRNA)
data(TCGAmrna)

Pheno_CRC <- TCGAmrna@phenoData@data
Geno_CRC <- exprs(TCGAmrna)
annot <- TCGAmrna@featureData@data

rownames(COADREAD.clin.merged) <- COADREAD.clin.merged[,1]
pheno_stage <- data.frame(t(COADREAD.clin.merged[c("patient.bcr_patient_barcode","patient.stage_event.pathologic_stage"), ]))
pheno_stage <- pheno_stage[-1, ]
stages <- table(pheno_stage$patient.stage_event.pathologic_stage)
names(stages)
pheno_stage$stage <- rep(NA, nrow(pheno_stage))
pheno_stage$stage[which(pheno_stage$patient.stage_event.pathologic_stage%in%names(stages)[1:2])]="I"
pheno_stage$stage[which(pheno_stage$patient.stage_event.pathologic_stage%in%names(stages)[3:6])]="II"
pheno_stage$stage[which(pheno_stage$patient.stage_event.pathologic_stage%in%names(stages)[7:10])]="III"
pheno_stage$stage[which(pheno_stage$patient.stage_event.pathologic_stage%in%names(stages)[11:13])]="IV"
pheno_stage$stage0 <- rep(NA, nrow(pheno_stage))
pheno_stage$stage0[which(pheno_stage$patient.stage_event.pathologic_stage%in%names(stages)[1:2])]=1
pheno_stage$stage0[which(pheno_stage$patient.stage_event.pathologic_stage%in%names(stages)[3:6])]=2
pheno_stage$stage0[which(pheno_stage$patient.stage_event.pathologic_stage%in%names(stages)[7:10])]=3
pheno_stage$stage0[which(pheno_stage$patient.stage_event.pathologic_stage%in%names(stages)[11:13])]=4

pheno_stage <- pheno_stage[match(substr(Pheno_CRC$ID,1,12),toupper(pheno_stage$patient.bcr_patient_barcode)), ]
pheno_stage <- pheno_stage[!is.na(pheno_stage$patient.bcr_patient_barcode), ]
pheno_stage <- pheno_stage[!is.na(pheno_stage$patient.stage_event.pathologic_stage), ]

save(pheno_stage, Pheno_CRC, Geno_CRC, file = "OUTPUT_CAO_new/DATA_TCGA_CRC.RData")

load("OUTPUT_CAO_new/DATA_TCGA_CRC.RData")
my_cols <- c("#84BD00FF", "#FFCD00FF", "#ff7400", "#CC0C00FF")
index <- which(substr(Pheno_CRC$ID,1,12)%in%toupper(pheno_stage$patient.bcr_patient_barcode))
DEGS <- c("ADH1C","ATP1B1","C1R","CLCA1","DCN","FCGBP","TFF1", "TFF3","ITM2C","MT1E","MT1F","MT1G","MT1H",
          "MUC2","MMP2", "MYL9","PHGR1","SPINK4", "S100A10","SLC26A3")
DEGS <- c("MT1E","MMP2","MUC2","CEACAM5","PIGR")


#p_val <- c()
for (i in 1:length(DEGS)) {
  gene_name <- DEGS[i]
  #gene_name <- 'MUC2'
  GENE = which(annot$Symbol==gene_name)
  
  # #GENE = which(annot$Symbol=="MT1E")
  # pheno_stage$GENE <- Geno_CRC[GENE, index]
  # #fit = lm(stage0 ~ GENE, data = pheno_stage)
  # #summary(fit)
  # df_plot1 <- pheno_stage %>% group_by(stage) %>% summarise_at(vars(GENE),list(mean=mean,  Q1=~quantile(., probs = 0.25),Q3=~quantile(., probs = 0.75)))
  # 
  # GENE = which(annot$Symbol=="MUC2")
  # pheno_stage$GENE <- Geno_CRC[GENE, index]
  # df_plot2 <- pheno_stage %>% group_by(stage) %>% summarise_at(vars(GENE),list(mean=mean,  Q1=~quantile(., probs = 0.25),Q3=~quantile(., probs = 0.75)))
  # 
  # GENE = which(annot$Symbol=="MMP2")
  # pheno_stage$GENE <- Geno_CRC[GENE, index]
  # df_plot3 <- pheno_stage %>% group_by(stage) %>% summarise_at(vars(GENE),list(mean=mean,  Q1=~quantile(., probs = 0.25),Q3=~quantile(., probs = 0.75)))
  # 
  # df_plot <- rbind(df_plot1, df_plot2, df_plot3)
  # df_plot$gene <- rep(c("MT1E", "MUC2","PIGR"), each=4)
  # 
  # ggplot(df_plot, aes(x=stage, y=mean, group=gene, color=gene))+geom_errorbar(aes(ymin=Q1, ymax=Q3), width=.1) +
  #   geom_line() + theme_bw()+
  #   geom_point()+ theme(text = element_text(size=20), legend.position = "bottom", panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
  #                       legend.text = element_text(size=15),legend.title = element_text(size=20))
  # 
  #[which(pheno_stage$patient.bcr_patient_barcode!="tcga-aa-3695"),]
  # g <- ggplot(data=pheno_stage, aes(y=GENE,x=stage, fill=stage))+
  #   geom_boxplot()+theme_bw()+scale_fill_manual(values=alpha(my_cols,0.9))+
  #   theme(text = element_text(size=20),
  #         legend.text = element_text(size=15),
  #         legend.title=element_blank(), legend.position = 'bottom') +
  #   ylab("Gene expression")+xlab("Stage")+
  #   ylim(5, 12)+
  #   stat_summary(fun = median,
  #                geom = "line",
  #                aes(group = 1),
  #                col = "blue")+
  # 
  
  pheno_stage$GENE <- Geno_CRC[GENE, index]
  g <- ggplot(data=pheno_stage, aes(y=GENE,x=stage, color=stage))+
    geom_violin(linewidth = 1.5)+theme_bw()+scale_fill_manual(values=alpha(my_cols,0.9))+
    theme(text = element_text(size=20),
          legend.text = element_text(size=15),
          legend.title=element_blank(), legend.position = 'none') +
    ylab("Gene expression")+xlab("Stage")+
    scale_color_manual(values=alpha(my_cols,0.9))+
    stat_summary(fun = "mean",
                 geom = "line",
                 aes(group = 1),
                 col = "#00316E") +
    stat_summary(fun = "mean",
                 geom = "point",
                 color = "#b2182b")
  
  pdf(paste0("OUTPUT_CAO_new/Figures_C4/VIOLIN_GE_TCGA_", gene_name, ".pdf")  )
  print(g)
  dev.off()
 
  
}

CANCER_DEG <- read_excel("~/Documents/Supplemental_Table_S2.xlsx")

colnames(CANCER_DEG) <- CANCER_DEG[1, ]
CANCER_DEG <- CANCER_DEG[-1, ]
COAD_DEG <- CANCER_DEG[which(CANCER_DEG$Cancer=="COAD"), ]


DEGS[19]

GENE = which(annot$Symbol==DEGS[i])
pheno_stage$GENE <- Geno_CRC[GENE, index]

ggplot(data=pheno_stage, aes(y=GENE, fill=stage))+
  geom_boxplot()+theme_bw()



+ 
  ylab("GSVA score")+xlab("")+
  #geom_jitter(color="black", size=0.4, alpha=0.9) +
  stat_compare_means(method = "wilcox.test",label = "p.format",size = 8, label.y = 1.5)+
  theme(text = element_text(size=30),
        legend.text = element_text(size=20),
        legend.title=element_blank()) +
  scale_fill_manual(values=alpha(my_col2,1))
dev.off()


pdf(file=paste0("OUTPUT_CAO_new/FIG/FeaturePlot_", name_st, "_","PIGR", ".pdf"))
SpatialFeaturePlot(st_dataset, features = "PIGR")
dev.off()


##### Transfer to other ST data
my_cols <- c("Other"='#7C878EFF',
             "T1" ='#ff7400',
             "T2"='#CC0C00FF')
                              
c4_dataset <- st_dataset
name_st <- "C3"
load(file=paste0("DATA_STSC_CAO/Seurat_",name_st,".RData"))

anchors <- FindTransferAnchors(reference=c4_dataset, query = st_dataset, dims = 1:30, normalization.method = 'SCT')
predictions <- TransferData(anchorset = anchors, refdata = as.character(c4_dataset@meta.data$T12), dims = 1:30)

st_dataset@meta.data$deconv_c4 <- predictions$prediction.score.T2
st_dataset@meta.data$deconv_cat <- predictions$predicted.id
#st_dataset@meta.data$deconv_cat[which(st_dataset@meta.data$deconv_cat=="Other")] <- "PT"

png(paste0(work_dir, "/FIG/SpatialFeaturePlot_deconv_C4_", name_st, ".png"), width=600, height=450, unit='px')
SpatialFeaturePlot(st_dataset, features = "deconv_c4")+
  scale_fill_gradientn(colours = col[c(1:10)])+ theme(legend.position="none")
dev.off()


pdf(paste0(work_dir, "/FIG/SpatialFeaturePlot_deconv_C4_cat_", name_st, ".pdf"))
SpatialDimPlot(st_dataset, group.by = "deconv_cat",cols = my_cols)
dev.off()

