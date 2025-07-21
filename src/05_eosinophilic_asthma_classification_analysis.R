library(corrplot)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(circlize)
library(ggpubr)
library(reshape2)
library(limma)
library(openxlsx)
library(caret)

# Set working directory and global variables:
root.dir <- dirname(dirname(rstudioapi::getActiveDocumentContext()$path))
print(root.dir)
setwd(root.dir)

results.dir <- paste0(root.dir, format(Sys.Date(), "/results/%Y%m%d/"))
results.tables<- paste0(results.dir,"/tables/")
results.figures <- paste0(results.dir,"/figures/")
dir.create(results.dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results.tables, recursive = TRUE, showWarnings = FALSE)
dir.create(results.figures, recursive = TRUE, showWarnings = FALSE)


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

meta2$ICS_YN <- ifelse(meta2$QMA_YN_ICS_Einzelinhalator == 1|meta2$QMA_YN_SystemicSteroids == 1,1,0) %>% as.factor()
meta2$OCS_YN <- meta2$QMA_YN_SystemicSteroids %>% as.factor()

meta2[, colnames(df_port)] <- removeBatchEffect(t(meta2[, colnames(df_port)]),batch = meta2$year) %>% t()






sputum <- ifelse(meta2$MS_DIFF_EOS >= 3,1,0)

bc150 <- ifelse(meta2$MB_DIFF_EOS_abs_corrected.FZB * 1000 >= 150,1,0)
bc300 <- ifelse(meta2$MB_DIFF_EOS_abs_corrected.FZB * 1000 >= 300,1,0)
bc400 <- ifelse(meta2$MB_DIFF_EOS_abs_corrected.FZB * 1000 >= 400,1,0)
perc2.70 <- ifelse(meta2$MB_DIFF_EOS_prec_corrected.FZB >= 2.70,1,0)
ENR0.05 <- ifelse((meta2$MB_DIFF_EOS_abs_corrected.FZB * 1000) / (meta2$MB_DIFF_NEUT_abs * 1000) >= 0.05,1,0)

bc150_FENO50 <- ifelse(meta2$MB_DIFF_EOS_abs_corrected.FZB * 1000 >= 150 & meta2$L_FeNO_Value >= 50,1,0)
bc300_FENO50 <- ifelse(meta2$MB_DIFF_EOS_abs_corrected.FZB * 1000 >= 300 & meta2$L_FeNO_Value >= 50,1,0)
bc400_FENO50 <- ifelse(meta2$MB_DIFF_EOS_abs_corrected.FZB * 1000 >= 400 & meta2$L_FeNO_Value >= 50,1,0)
perc2.70_FENO50 <- ifelse(meta2$MB_DIFF_EOS_prec_corrected.FZB >= 2.70 & meta2$L_FeNO_Value >= 50,1,0)
ENR0.05_FENO50 <- ifelse((meta2$MB_DIFF_EOS_abs_corrected.FZB * 1000) / (meta2$MB_DIFF_NEUT_abs * 1000) >= 0.05 & meta2$L_FeNO_Value >= 50,1,0)

predictors <- data.frame(
  sputum,
  bc150,
  bc150_FENO50,
  bc300,
  bc300_FENO50,
  bc400,
  bc400_FENO50,
  perc2.70,
  perc2.70_FENO50,
  ENR0.05,
  ENR0.05_FENO50
)

predictors <- predictors[which(!is.na(rowSums(predictors))),]

accuracy <- c()
accuracy.lower <- c()
accuracy.upper <- c()
F1 <- c()
false.positive <- c()
true.positive <- c()
false.negative <- c()
true.negative <- c()

for(i in seq(2,ncol(predictors))){
  
  dat <- confusionMatrix(as.factor(predictors[,i]),reference = as.factor(predictors[,1]),mode = "everything",)
  
  F1 <- c(F1,dat$byClass["F1"])
  
  accuracy <- c(accuracy,
                dat$overall[1])
  
  accuracy.lower <- c(accuracy.lower,
                      dat$overall[3])
  accuracy.upper <- c(accuracy.upper,
                      dat$overall[4])
  
  false.negative <- c(false.negative,
                      dat$table[1,2])
  true.negative <- c(true.negative,
                     dat$table[1,1])
  
  false.positive <- c(false.positive,
                      dat$table[2,1])
  true.positive <- c(true.positive,
                     dat$table[2,2])
  
}

test.results <- data.frame(
  Threshold = c("Blood ≥ 150 cells",
                "Blood ≥ 150 cells & FeNO ≥ 50 ppb",
                "Blood ≥ 300 cells",
                "Blood ≥ 300 cells & FeNO ≥ 50 ppb",
                "Blood ≥ 400 cells",
                "Blood ≥ 400 cells & FeNO ≥ 50 ppb",
                "Blood percentage ≥ 2.70",
                "Blood percentage ≥ 2.70 & FeNO ≥ 50 ppb",
                "Blood ENR ≥ 0.05",
                "Blood ENR ≥ 0.05 & FeNO ≥ 50 ppb"),
  F1=F1,
  Accuracy = accuracy,
  Accuracy.lower = accuracy.lower,
  accuracy.upper = accuracy.upper,
  FP = false.positive,
  FN = false.negative,
  TP = true.positive,
  TN = true.negative
)

write.csv(x = test.results,file = paste0(prefix.tables,"F1_scores.csv"),quote = FALSE)


to.plot <- melt(test.results[,c(1,2)])

p <- ggplot(to.plot, aes(x=value, y=Threshold, fill=variable)) + 
  geom_bar(position="dodge", stat="identity") +
  geom_text(aes(label=round(value,3)), hjust=1.6, vjust=0.5, color="white",
            position = position_dodge(0.9), size=3.5) +
  xlab("Score") +
  scale_fill_manual(values = c("#d8b365","#5ab4ac")) +
  scale_y_discrete(labels=c(quote(paste(Count >= 150," ",cells)),
                            quote(paste(Count >= 150," ",cells," ",and," ",FeNO >= 50," ",ppb)),
                            quote(paste(Count >= 300," ",cells)),
                            quote(paste(Count >= 300," ",cells," ",and," ",FeNO >= 50," ",ppb)),
                            quote(paste(Count >= 400," ",cells)),
                            quote(paste(Count >= 400," ",cells," ",and," ",FeNO >= 50," ",ppb)),
                            quote(paste(ENR >= 0.05)),
                            quote(paste(ENR >= 0.05," ",and," ",FeNO >= 50," ",ppb)),
                            quote(paste(Percentage >= 2.70)),
                            quote(paste(Percentage >= 2.70," ",and," ",FeNO >= 50," ",ppb))
  )) +
  theme_bw() + theme(axis.text.x = element_text(size=10),
                     axis.text.y = element_text(size=10),
                     axis.title.x = element_text(size=10,face="bold"),
                     axis.title.y = element_blank(), legend.position = "none", 
                     legend.text = element_text(size=10))

ggsave(filename = paste0(prefix.plots,"F1_scores.pdf"),plot = p,units = "cm",width = 12,height = 8)


###########################
#F1-score (bootstrapping)
###########################

meta3 <- meta2[which(!is.na(rowSums(meta2[,c("MB_DIFF_EOS_abs_corrected.FZB","MS_DIFF_EOS","MB_DIFF_EOS_prec_corrected.FZB","L_FeNO_Value")]))),]
boost.F1 <- data.frame()
result.list <- list()

set.seed(1234)
for(fk in seq(100)){
  samples <- sample(1:nrow(meta3), replace = TRUE)
  
  result.list[[fk]] <- meta3[samples,c(colnames(df_port),"MB_DIFF_EOS_abs_corrected.FZB","MS_DIFF_EOS","MB_DIFF_EOS_prec_corrected.FZB","L_FeNO_Value")]
  
  sputum <- ifelse(meta3$MS_DIFF_EOS[samples] >= 3,1,0)
  
  bc150 <- ifelse(meta3$MB_DIFF_EOS_abs_corrected.FZB[samples] * 1000 >= 150,1,0)
  bc300 <- ifelse(meta3$MB_DIFF_EOS_abs_corrected.FZB[samples] * 1000 >= 300,1,0)
  bc400 <- ifelse(meta3$MB_DIFF_EOS_abs_corrected.FZB[samples] * 1000 >= 400,1,0)
  perc2.70 <- ifelse(meta3$MB_DIFF_EOS_prec_corrected.FZB[samples] >= 2.70,1,0)
  ENR0.05 <- ifelse((meta3$MB_DIFF_EOS_abs_corrected.FZB[samples] * 1000) / (meta3$MB_DIFF_NEUT_abs[samples] * 1000) >= 0.05,1,0)
  
  bc150_FENO50 <- ifelse(meta3$MB_DIFF_EOS_abs_corrected.FZB[samples] * 1000 >= 150 & meta3$L_FeNO_Value[samples] >= 50,1,0)
  bc300_FENO50 <- ifelse(meta3$MB_DIFF_EOS_abs_corrected.FZB[samples] * 1000 >= 300 & meta3$L_FeNO_Value[samples] >= 50,1,0)
  bc400_FENO50 <- ifelse(meta3$MB_DIFF_EOS_abs_corrected.FZB[samples] * 1000 >= 400 & meta3$L_FeNO_Value[samples] >= 50,1,0)
  perc2.70_FENO50 <- ifelse(meta3$MB_DIFF_EOS_prec_corrected.FZB[samples] >= 2.70 & meta3$L_FeNO_Value[samples] >= 50,1,0)
  ENR0.05_FENO50 <- ifelse((meta3$MB_DIFF_EOS_abs_corrected.FZB[samples] * 1000) / (meta3$MB_DIFF_NEUT_abs[samples] * 1000) >= 0.05 & meta3$L_FeNO_Value[samples] >= 50,1,0)
  
  predictors <- data.frame(
    sputum,
    bc150,
    bc150_FENO50,
    bc300,
    bc300_FENO50,
    bc400,
    bc400_FENO50,
    perc2.70,
    perc2.70_FENO50,
    ENR0.05,
    ENR0.05_FENO50
  )
  
  F1 <- c()
  
  for(i in seq(2,ncol(predictors))){
    
    dat <- confusionMatrix(as.factor(predictors[,i]),reference = as.factor(predictors[,1]),mode = "everything",)
    
    F1 <- c(F1,dat$byClass["F1"])
  }
  
  boost.F1 <- rbind(boost.F1,F1)
}

colnames(boost.F1) <- colnames(predictors)[-c(1)]

to.plot <- melt(boost.F1)

p <- ggplot(to.plot, aes(x=variable, y=value)) + 
  geom_boxplot(aes(fill=variable),outlier.shape = NA) + 
  geom_jitter(size=0.4) +
  theme_bw() +
  labs(title = "F1-score distributions",
       subtitle = "Bootstraping 100 times",
       x = "Cutoff",
       y = "F1-score") +
  scale_fill_manual(values = RColorBrewer::brewer.pal(n = 10,name = "Paired")) +
  theme(axis.text.x = element_blank(), axis.text.y = element_text(size=8), axis.title=element_text(size=8), 
        legend.title = element_text(size=8), legend.position = "none", text = element_text(size=8))

ggsave(filename = paste0(prefix.plots,"F1_bootstraping.pdf"),plot = p,units = "cm",width = 12,height = 8)

#### SPUTUM

dummy.meta2 <- meta2[which(meta2$D_ASTHMA_SEVERITYGRADE_SCREEN != 0),]
dummy.meta2 <- dummy.meta2[which(!is.na(dummy.meta2$MS_DIFF_EOS)),]
dummy.meta2$EOS <- ifelse(dummy.meta2$MS_DIFF_EOS >= 3,"High","Low")

to.box.plot <- data.frame()
for (i in colnames(df_port)) {
  to.bind <-  data.frame(rep(i,nrow(dummy.meta2)),
                         dummy.meta2[,i],
                         dummy.meta2$EOS,row.names = NULL)
  colnames(to.bind) <- c('Protein','Expression','Condition')
  to.box.plot <- rbind(to.box.plot,to.bind)
}

to.box.plot$Condition <- factor(to.box.plot$Condition,levels = c("Low","High"))

to.sig <- c()
to.stats <- data.frame()
for(cyt in colnames(df_port)){
  to.stats <- rbind(
    to.stats,
    data.frame(
      cyt = cyt,
      Wilcoxon = wilcox.test(dummy.meta2[which(dummy.meta2$EOS == "High"), which(colnames(dummy.meta2) == cyt)],
                             dummy.meta2[which(dummy.meta2$EOS == "Low"), which(colnames(dummy.meta2) == cyt)])$statistic,
      log2FC = mean(dummy.meta2[which(dummy.meta2$EOS == "High"), which(colnames(dummy.meta2) == cyt)]) - mean(dummy.meta2[which(dummy.meta2$EOS == "Low"), which(colnames(dummy.meta2) == cyt)]),
      p = wilcox.test(dummy.meta2[which(dummy.meta2$EOS == "High"), which(colnames(dummy.meta2) == cyt)],
                      dummy.meta2[which(dummy.meta2$EOS == "Low"), which(colnames(dummy.meta2) == cyt)])$p.val
    )
  )
  to.sig <- c(to.sig,wilcox.test(dummy.meta2[which(dummy.meta2$EOS == "High"), which(colnames(dummy.meta2) == cyt)],
                                 dummy.meta2[which(dummy.meta2$EOS == "Low"), which(colnames(dummy.meta2) == cyt)])$p.val)
}

to.stats$adj.p <- to.sig <- p.adjust(to.sig,method = "BH")
names(to.sig) <- colnames(df_port)

sputum.fig <- to.box.plot
sputum.table <- to.sig

rmarkdown::paged_table(to.stats %>% arrange(adj.p)) %>% write.csv(file = paste0(prefix.tables,"sputum_wilcoxon.csv"),quote = FALSE)


### BLOOD ≥ 150 cells + FeNO >= 50 ppb

dummy.meta2 <- meta2[which(meta2$D_ASTHMA_SEVERITYGRADE_SCREEN != 0),]
dummy.meta2 <- dummy.meta2[which(!is.na(dummy.meta2$MB_DIFF_EOS_abs_corrected.FZB) & !is.na(dummy.meta2$L_FeNO_Value)),]
dummy.meta2$EOS <- ifelse(dummy.meta2$MB_DIFF_EOS_abs_corrected.FZB * 1000 >= 150 & dummy.meta2$L_FeNO_Value >= 50,"High","Low")

to.box.plot <- data.frame()
for (i in colnames(df_port)) {
  to.bind <-  data.frame(rep(i,nrow(dummy.meta2)),
                         dummy.meta2[,i],
                         dummy.meta2$EOS,row.names = NULL)
  colnames(to.bind) <- c('Protein','Expression','Condition')
  to.box.plot <- rbind(to.box.plot,to.bind)
}

to.box.plot$Condition <- factor(to.box.plot$Condition,levels = c("Low","High"))

to.sig <- c()
to.stats <- data.frame()
for(cyt in colnames(df_port)){
  to.stats <- rbind(
    to.stats,
    data.frame(
      cyt = cyt,
      Wilcoxon = wilcox.test(dummy.meta2[which(dummy.meta2$EOS == "High"), which(colnames(dummy.meta2) == cyt)],
                             dummy.meta2[which(dummy.meta2$EOS == "Low"), which(colnames(dummy.meta2) == cyt)])$statistic,
      log2FC = mean(dummy.meta2[which(dummy.meta2$EOS == "High"), which(colnames(dummy.meta2) == cyt)]) - mean(dummy.meta2[which(dummy.meta2$EOS == "Low"), which(colnames(dummy.meta2) == cyt)]),
      p = wilcox.test(dummy.meta2[which(dummy.meta2$EOS == "High"), which(colnames(dummy.meta2) == cyt)],
                      dummy.meta2[which(dummy.meta2$EOS == "Low"), which(colnames(dummy.meta2) == cyt)])$p.val
    )
  )
  to.sig <- c(to.sig,wilcox.test(dummy.meta2[which(dummy.meta2$EOS == "High"), which(colnames(dummy.meta2) == cyt)],
                                 dummy.meta2[which(dummy.meta2$EOS == "Low"), which(colnames(dummy.meta2) == cyt)])$p.val)
}

to.stats$adj.p <- to.sig <- p.adjust(to.sig,method = "BH")
names(to.sig) <- colnames(df_port)

blood.150.fig <- to.box.plot
blood.150.table <- to.sig

rmarkdown::paged_table(to.stats %>% arrange(adj.p)) %>% write.csv(file = paste0(prefix.tables,"blood150_feno50_wilcoxon.csv"),quote = FALSE)



### BLOOD ≥ 2.7% + FeNO >= 50 ppb

dummy.meta2 <- meta2[which(meta2$D_ASTHMA_SEVERITYGRADE_SCREEN != 0),]
dummy.meta2 <- dummy.meta2[which(!is.na(dummy.meta2$MB_DIFF_EOS_prec_corrected.FZB) & !is.na(dummy.meta2$L_FeNO_Value)),]
dummy.meta2$EOS <- ifelse(dummy.meta2$MB_DIFF_EOS_prec_corrected.FZB >= 2.70 & dummy.meta2$L_FeNO_Value >= 50,"High","Low")

to.box.plot <- data.frame()
for (i in colnames(df_port)) {
  to.bind <-  data.frame(rep(i,nrow(dummy.meta2)),
                         dummy.meta2[,i],
                         dummy.meta2$EOS,row.names = NULL)
  colnames(to.bind) <- c('Protein','Expression','Condition')
  to.box.plot <- rbind(to.box.plot,to.bind)
}

to.box.plot$Condition <- factor(to.box.plot$Condition,levels = c("Low","High"))

to.sig <- c()
to.stats <- data.frame()
for(cyt in colnames(df_port)){
  to.stats <- rbind(
    to.stats,
    data.frame(
      cyt = cyt,
      Wilcoxon = wilcox.test(dummy.meta2[which(dummy.meta2$EOS == "High"), which(colnames(dummy.meta2) == cyt)],
                             dummy.meta2[which(dummy.meta2$EOS == "Low"), which(colnames(dummy.meta2) == cyt)])$statistic,
      log2FC = mean(dummy.meta2[which(dummy.meta2$EOS == "High"), which(colnames(dummy.meta2) == cyt)]) - mean(dummy.meta2[which(dummy.meta2$EOS == "Low"), which(colnames(dummy.meta2) == cyt)]),
      p = wilcox.test(dummy.meta2[which(dummy.meta2$EOS == "High"), which(colnames(dummy.meta2) == cyt)],
                      dummy.meta2[which(dummy.meta2$EOS == "Low"), which(colnames(dummy.meta2) == cyt)])$p.val
    )
  )
  to.sig <- c(to.sig,wilcox.test(dummy.meta2[which(dummy.meta2$EOS == "High"), which(colnames(dummy.meta2) == cyt)],
                                 dummy.meta2[which(dummy.meta2$EOS == "Low"), which(colnames(dummy.meta2) == cyt)])$p.val)
}

to.stats$adj.p <- to.sig <- p.adjust(to.sig,method = "BH")
names(to.sig) <- colnames(df_port)

blood.27.fig <- to.box.plot
blood.27.table <- to.sig

rmarkdown::paged_table(to.stats %>% arrange(adj.p)) %>% write.csv(file = paste0(prefix.tables,"blood_percentage_feno50_wilcoxon.csv"),quote = FALSE)

##########################

cyt.table <- data.frame(
  Cytokine = names(sputum.table),
  Sputum = sputum.table,
  `Blood ≥ 150 + FeNO ≥ 50 ppb` = blood.150.table,
  `Blood ≥ 2.7% + FeNO ≥ 50 ppb` = blood.27.table
)


rmarkdown::paged_table(cyt.table)


sputum.fig$Type <- "Sputum ≥ 3%"
blood.150.fig$Type <- "Blood ≥ 150 cells/uL & FeNO ≥ 50 ppb"
blood.27.fig$Type <- "Blood ≥ 2.7% & FeNO ≥ 50 ppb"

to.box.plot <- rbind(sputum.fig,blood.150.fig)
to.box.plot$Type <- as.factor(to.box.plot$Type)
to.box.plot$Type <- relevel(x = to.box.plot$Type,ref = "Sputum ≥ 3%")

cyt <- 'Eotaxin.3'

p <- ggplot(data = subset(to.box.plot, (Protein == cyt)), aes(x = Type, y = Expression,fill = Condition)) + 
  geom_boxplot(width = 0.4,outlier.shape = NA) + theme_bw() +
  geom_jitter(size=0.7, position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.5)) +
  theme(axis.text.x = element_text(size=8), axis.text.y = element_text(size=8), axis.title=element_text(size=8), 
        legend.title = element_text(size=8), legend.position = "none", text = element_text(size=8)) +
  ylim(c(min(subset(to.box.plot, (Protein == cyt))$Expression) - 2,
         max(subset(to.box.plot, (Protein == cyt))$Expression)) + 2) +
  labs(
    title="Wilcoxon rank-sum - p.adj",
    subtitle="CCL-26",
    x="Cutoff",
    y="Normalized abundances"
  ) +
  scale_x_discrete(labels=c(
    quote(paste(Sputum >= 3,"%")),
    quote(paste(Blood >= 150," ",cells," & ",FeNO >= 50," ",ppb))
  )) +
  scale_fill_manual(name= "",labels =c("Low","High"), values= c("#fdbb84","#e34a33")) +
  geom_signif(y_position=max(subset(to.box.plot, (Protein == cyt & Type == "Sputum ≥ 3%"))$Expression) + 1, 
              xmin=c("Sputum ≥ 3%"), xmax=c("Sputum ≥ 3%"),annotation=c(formatC(cyt.table[which(cyt.table$Cytokine == cyt),c(2)],digits = 3)), tip_length=0,textsize = 3,linetype = "blank") +
  geom_signif(y_position=max(subset(to.box.plot, (Protein == cyt & Type == "Blood ≥ 150 cells/uL & FeNO ≥ 50 ppb"))$Expression) + 1, 
              xmin=c("Blood ≥ 150 cells/uL & FeNO ≥ 50 ppb"), xmax=c("Blood ≥ 150 cells/uL & FeNO ≥ 50 ppb"),annotation=c(formatC(cyt.table[which(cyt.table$Cytokine == cyt),c(3)],digits = 3)), tip_length=0,textsize = 3,linetype = "blank")

ggsave(filename = paste0(prefix.plots,"ccl26_wilcoxon.pdf"),plot = p,units = "cm",width = 10,height = 8)


# IL-5

cyt <- 'IL.5'

p <- ggplot(data = subset(to.box.plot, (Protein == cyt)), aes(x = Type, y = Expression,fill = Condition)) + 
  geom_boxplot(width = 0.4,outlier.shape = NA) + theme_bw() +
  geom_jitter(size=0.7, position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.5)) +
  theme(axis.text.x = element_text(size=8), axis.text.y = element_text(size=8), axis.title=element_text(size=8), 
        legend.title = element_text(size=8), legend.position = "none", text = element_text(size=8)) +
  ylim(c(min(subset(to.box.plot, (Protein == cyt))$Expression) - 2,
         max(subset(to.box.plot, (Protein == cyt))$Expression)) + 2) +
  labs(
    title="Wilcoxon rank-sum - p.adj",
    subtitle="IL-5",
    x="Cutoff",
    y="Normalized abundances"
  ) +
  scale_x_discrete(labels=c(
    quote(paste(Sputum >= 3,"%")),
    quote(paste(Blood >= 150," ",cells," & ",FeNO >= 50," ",ppb))
  )) +
  scale_fill_manual(name= "",labels =c("Low","High"), values= c("#fdbb84","#e34a33")) +
  geom_signif(y_position=max(subset(to.box.plot, (Protein == cyt & Type == "Sputum ≥ 3%"))$Expression) + 1, 
              xmin=c("Sputum ≥ 3%"), xmax=c("Sputum ≥ 3%"),annotation=c(formatC(cyt.table[which(cyt.table$Cytokine == cyt),c(2)],digits = 3)), tip_length=0,textsize = 3,linetype = "blank") +
  geom_signif(y_position=max(subset(to.box.plot, (Protein == cyt & Type == "Blood ≥ 150 cells/uL & FeNO ≥ 50 ppb"))$Expression) + 1, 
              xmin=c("Blood ≥ 150 cells/uL & FeNO ≥ 50 ppb"), xmax=c("Blood ≥ 150 cells/uL & FeNO ≥ 50 ppb"),annotation=c(formatC(cyt.table[which(cyt.table$Cytokine == cyt),c(3)],digits = 3)), tip_length=0,textsize = 3,linetype = "blank")

ggsave(filename = paste0(prefix.plots,"il5_wilcoxon.pdf"),plot = p,units = "cm",width = 10,height = 8)

p
