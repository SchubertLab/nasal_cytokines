#title: "Significance evaluation per cytokine expression"
#author: Juan Henao
#description: " Differential abundance analysis with KS and post-hoc Wilcoxon"


library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(openxlsx)
library(limma)
library(ggsignif)
library(reshape2)
library(ComplexHeatmap)
library(circlize)


# Set working directory and global variables:
root.dir <- dirname(dirname(rstudioapi::getActiveDocumentContext()$path))
print(root.dir)
setwd(root.dir)

results.dir <- paste0(root.dir, format(Sys.Date(), "/results/%Y%m%d/"))
prefix.tables<- paste0(results.dir,"/tables/")
prefix.plots <- paste0(results.dir,"/figures/")
dir.create(results.dir, recursive = TRUE, showWarnings = FALSE)
dir.create(prefix.tables, recursive = TRUE, showWarnings = FALSE)
dir.create(prefix.plots, recursive = TRUE, showWarnings = FALSE)


meta <- read.csv("./data/raw_data/meta.csv", sep=";", header=TRUE,stringsAsFactors = FALSE)
meta2 <- meta[meta["type"] == "measured",]
rownames(meta2)<-meta2$ID_pseudo

df_port <- read.table("./data/derivative_data/20220803_plate_normalized20220803_data_donwshift_imputed.csv", sep=" ", header=TRUE, stringsAsFactors = FALSE)

to.remove <- colnames(meta[,9:27])[!colnames(meta[,9:27]) %in% colnames(df_port)]
meta <- meta[,which(!colnames(meta) %in% to.remove)]

meta2$D_ASTHMA_SEVERITYGRADE_SCREEN[meta2$D_ASTHMA_SEVERITYGRADE_SCREEN == 9]<-0
meta2[rownames(meta2),colnames(df_port)] <- df_port[rownames(meta2),colnames(df_port)]

meta2$year <- as.factor(meta2$year)
meta2$plate_new <- as.factor(meta2$plate_new)
meta2$D_SEX<-as.factor(meta2$D_SEX)
meta2$D_ASTHMA_SEVERITYGRADE_SCREEN <-as.factor(meta2$D_ASTHMA_SEVERITYGRADE_SCREEN)
meta2$D_SmokingStatus<-as.factor(meta2$D_SmokingStatus)
meta2$OCS_mg <- ifelse(is.na(meta2$OCS_mg),0,meta2$OCS_mg)

meta2[, colnames(df_port)] <- removeBatchEffect(t(meta2[, colnames(df_port)]),batch = meta2$year) %>% t()


# Correlation

to.cor.plot <- meta2[,which(colnames(meta2) %in% colnames(df_port))] %>%
  cbind(meta2$MS_DIFF_EOS,
        meta2$MS_DIFF_NEUT,
        meta2$MB_DIFF_EOS_abs_corrected.FZB,
        meta2$MB_DIFF_NEUT_abs
  ) %>% na.omit()

colnames(to.cor.plot) <- c("CCL-26","G-CSF","IFN-g","IL-10",
                           "IL-13","IL-17","IL-1a","IL-37",
                           "IL-24","IL-33","IL-4","IL-5",
                           "IL-8","POSTN","SCGB1A1","TNF-a",
                           "Sputum Eosinophils","Sputum Neutrophils","Blood Eosinophils","Blood Neutrophils")


corrplot(cor(to.cor.plot),method="num",number.cex = 0.5,tl.col = 'black',type = 'upper',col=colorRampPalette(c("#313695", "white", "#A50026"))(200))

corr_matrix <- cor(to.cor.plot)
sparse_corr <- corr_matrix
sparse_corr[abs(sparse_corr) < 0.8] <- 0
sparse_corr <- sparse_corr[which(rowSums(sparse_corr) != 1.0),colSums(sparse_corr) != 1.0]

pdf(file = paste0(prefix.plots,"correlation.pdf"),width = 7/2.54,height = 7/2.54)
corrplot(sparse_corr,method="num",number.cex = 0.6,cl.cex=0.6,tl.cex=0.6,tl.col = 'black',type = 'upper',col=colorRampPalette(c("#313695", "white", "#A50026"))(200))
dev.off()

corrplot(sparse_corr,method="num",number.cex = 0.6,cl.cex=0.6,tl.cex=0.6,tl.col = 'black',type = 'upper',col=colorRampPalette(c("#313695", "white", "#A50026"))(200))

write.csv(x = cor(to.cor.plot),file = paste0(prefix.tables,"all_correlation_matrix.csv"))

corr <- cor(to.cor.plot) %>% 
  melt(id.vars="X")

corr <- corr[,c(2,1,3)]

red.corr <- corr[which(abs(corr[,3]) > 0.6),]

pairs <- paste0(red.corr$Var2,"--",red.corr$Var1)
to.keep <- c()
for(i in seq(nrow(red.corr))){
  if(paste0(red.corr[i,2],"--",red.corr[i,1]) %in% pairs & red.corr[i,2] != red.corr[i,1]){
    to.keep <- c(to.keep,i)
    pairs <- pairs[-c(which(paste0(red.corr[i,1],"--",red.corr[i,2]) == pairs))]
  }
}

red.corr <- red.corr[to.keep,]

red.corr$value <- round(red.corr$value,2)
colnames(red.corr) <- c("First variable","Second Variable","Correlation value")

red.corr <- red.corr[order(red.corr$`Correlation value`,decreasing = TRUE),]

write.csv(x = red.corr,file = paste0(prefix.plots,"ordered_correlation.csv"))


# Differential abundance Analysis

k <- c()
u <- c()
for(cyt in colnames(df_port)){
  k <- c(k,kruskal.test(meta2[,which(colnames(meta2) == cyt)]~meta2$D_ASTHMA_SEVERITYGRADE_SCREEN)$p.value)
  u <- c(u,kruskal.test(meta2[,which(colnames(meta2) == cyt)]~meta2$D_ASTHMA_SEVERITYGRADE_SCREEN)$statistic)
}

k.s <- p.adjust(k,method = "BH")

names(k.s) <- colnames(df_port)

stats <- data.frame(
  cyt=colnames(df_port),
  Kruskal=u,
  p.value=k,
  adj.p=k.s
)

rmarkdown::paged_table(stats)

write.csv(x = stats,file = paste0(prefix.tables,"kruskal_values.csv"),quote = FALSE)

#

p.h <- data.frame()

for(cyt in colnames(df_port)){
  pair <- pairwise.wilcox.test(meta2[,which(colnames(meta2) == cyt)],meta2$D_ASTHMA_SEVERITYGRADE_SCREEN,p.adjust.method = "BH")
  
  p.h <- rbind(p.h,
               data.frame("wilcoxon p.adj 0vs1" = pair$p.value[1,1],
                          "logFC 0vs1" = mean(unlist(subset(meta2,D_ASTHMA_SEVERITYGRADE_SCREEN==1)[cyt])) - mean(unlist(subset(meta2,D_ASTHMA_SEVERITYGRADE_SCREEN==0)[cyt])),
                          "wilcoxon p.adj 0vs2" = pair$p.value[2,1],
                          "logFC 0vs2" = mean(unlist(subset(meta2,D_ASTHMA_SEVERITYGRADE_SCREEN==2)[cyt])) - mean(unlist(subset(meta2,D_ASTHMA_SEVERITYGRADE_SCREEN==0)[cyt])),
                          "wilcoxon p.adj 0vs3" = pair$p.value[3,1],
                          "logFC 0vs3" = mean(unlist(subset(meta2,D_ASTHMA_SEVERITYGRADE_SCREEN==3)[cyt])) - mean(unlist(subset(meta2,D_ASTHMA_SEVERITYGRADE_SCREEN==0)[cyt])),
                          "wilcoxon p.adj 1vs2" = pair$p.value[2,2],
                          "logFC 1vs2" =mean(unlist(subset(meta2,D_ASTHMA_SEVERITYGRADE_SCREEN==2)[cyt])) - mean(unlist(subset(meta2,D_ASTHMA_SEVERITYGRADE_SCREEN==1)[cyt])),
                          "wilcoxon p.adj 1vs3" = pair$p.value[3,2],
                          "logFC 1vs3" = mean(unlist(subset(meta2,D_ASTHMA_SEVERITYGRADE_SCREEN==3)[cyt])) - mean(unlist(subset(meta2,D_ASTHMA_SEVERITYGRADE_SCREEN==1)[cyt])),
                          "wilcoxon p.adj 2vs3" = pair$p.value[3,3],
                          "logFC 2vs3" = mean(unlist(subset(meta2,D_ASTHMA_SEVERITYGRADE_SCREEN==3)[cyt])) - mean(unlist(subset(meta2,D_ASTHMA_SEVERITYGRADE_SCREEN==2)[cyt]))
)
  )
}

all.data <- cbind(k.s,p.h)


write.csv(x = all.data,file = paste0(prefix.tables,"wilcoxon_adj_pvalues_values.csv"),quote = FALSE)

to.plot <- as.matrix(all.data[, c(c(1), seq(from = 2, to = ncol(all.data), by = 2))])
colnames(to.plot) <- c("Kruskal-test","Healthy vs Mild","Healthy vs Moderate","Healthy vs Severe",
                       "Mild vs Moderate","Mild vs Severe", "Moderate vs Severe")
new.rows <- c("CCL-26","G-CSF","IFN-g","IL-10","IL-13","IL-17","IL-1a","IL-37","IL-24",
              "IL-33","IL-4","IL-5","IL-8","POSTN","SCGB1A1",expression(paste("TNF-", alpha)))
#

col_fun = colorRamp2(c(0, 0.5, 1), c("#A50026","white","#313695"))
#col_fun(seq(-3, 3))

pdf(file = paste0(prefix.plots,"significant_heatmap.pdf"),width = 21/2.54,height = 8/2.54)
Heatmap(to.plot, col = col_fun, show_column_dend = FALSE, column_order = colnames(to.plot), show_column_names = FALSE, column_title_gp = grid::gpar(fontsize = 8),
        bottom_annotation = HeatmapAnnotation(text = anno_text(colnames(to.plot),rot=360,offset = unit(0.5, "npc"), just = "center",gp=gpar(fontsize = 8))),
        row_names_gp = grid::gpar(fontsize = 8), show_heatmap_legend = FALSE,
        row_labels = new.rows,column_split = c("Kruskal-test",rep("Wilcoxon rank sum test (post-hoc)",6)),
        layer_fun = function(j, i, x, y, width, height, fill) {
          v = pindex(to.plot, i, j)
          l = v > 0
          grid.text(sprintf("%.3f", v[l]), x[l], y[l], gp = gpar(fontsize = 8))
        })
dev.off()

