library(presto)
library(ComplexHeatmap)
library(circlize)
library(Seurat)
library(dplyr)
library(tidyverse)

setwd("/Users/w435u/Documents/ST_SC")
load("DATA_STSC_CAO/Seurat_C2.RData")

pred_st <- st_dataset@meta.data$pred_st
pred_cat <- as.character(cut(pred_st, 
                             breaks = c(-Inf, quantile(pred_st, c(0.25, 0.5, 0.75, 1))), 
                             labels = c("C1",  "C2", "C3", "C4")))
st_dataset <- AddMetaData(st_dataset, metadata = pred_cat, col.name = "predcat")

tumor_enriched <- FindMarkers(st_dataset, ident.1 = "C4",ident.2 = "C1", group.by = 'predcat', verbose = FALSE)
#tumor_enriched <- FindMarkers(st_dataset, ident.1 = "Positive", group.by = 'sign_C2', verbose = FALSE)
library(openxlsx)
write.xlsx(tumor_enriched, file="FIGURES/DEGs_P3.xlsx")

genes <- rownames(tumor_enriched)[1:40]

p <- DotPlot(object = st_dataset, features =genes, group.by = "predcat")
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


write.csv(exp_mat, file="FIGURES/Expression_FIG4H.csv")

dot <- Heatmap(exp_mat,
               heatmap_legend_param=list(title="expression"),
               column_title = "clustered dotplot", 
               col=col_fun,
               rect_gp = gpar(type = "none"),
               cell_fun = cell_fun,
               row_names_gp = gpar(fontsize = 10),
               row_km = 2,
               border = "black",
               cluster_columns = FALSE,column_order = c("C4","C3","C2", "C1"))



pdf("OUTPUT_CAO_new/FIG/Dotplot_C2_all_new.pdf")
print(dot)
dev.off()

#for (i in 1:length(genes)) {
#  gene= genes[i]
#  gg <- SpatialFeaturePlot(st_dataset, features = gene)
  
#  pdf(file=paste0("OUTPUT_CAO_new/FIG/FeaturePlot_", name_st, "_",gene, ".pdf"))
#  print(gg)
#  dev.off()
#}

