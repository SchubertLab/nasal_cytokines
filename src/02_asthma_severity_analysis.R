library(rstudioapi)
library(DirichletReg)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(openxlsx)
library(limma)
library(missRanger)
library(ComplexHeatmap)
library(M3C)
library(Rtsne)
library(ordinal)
library(effects)
library(AER)
library(MASS)
library(DescTools)
library(gridExtra)
library(sciplot)
library(ggsignif)
library(emmeans)
library(betareg)
library(openxlsx)
library(hdi)
library(glmnet)
library(compositions)
library(openxlsx)
library(brms)
library(bayestestR)
library(tidyr)
library(tidyverse)     



##### Function definitions:

clin_to_asthma_analysis <- function(dv, model, formula, data, wb, eps, dv.title, resultdir){
  
  fit <- model(formula, data = data)
  emm = emmeans(fit, specs = pairwise ~ D_ASTHMA_SEVERITYGRADE_SCREEN)
  pairwise_results <- summary(emm$contrasts, infer = TRUE)
  
  wb_sheet_title<-sub("\\]",")",sub("\\[", "(", sub("/", "_", dv.title)))
  addWorksheet(wb, wb_sheet_title)
  writeData(wb, wb_sheet_title, pairwise_results)
  
  
  base_y <- max(data[dv], na.rm = TRUE)
  buffer <- 0.05 * base_y
  
  if (!grepl("AQ5", dv.title)) {
    level<-c("healthy", "mild", "moderate", "severe")
    col <- c("#fceea9","#f3b55f", "#e95b3b", "#ad232c")
    k<-1
  } else {
    level<-c( "mild", "moderate", "severe")
    col <- c("#f3b55f", "#e95b3b", "#ad232c")
    k<-0
  }
  
  signif_data <- data.frame(
    xstart = as.numeric(sub("D_ASTHMA_SEVERITYGRADE_SCREEN", "", sub(" - .*", "", pairwise_results$contrast)))+k,
    xend = as.numeric(sub("D_ASTHMA_SEVERITYGRADE_SCREEN", "", sub(".* - ", "", pairwise_results$contrast)))+k,
    y = base_y + seq(from = buffer, by = buffer, length.out = length(pairwise_results$contrast)),
    p.value = pairwise_results$p.value
  )
  
  # Add asterisks based on significance levels
  signif_data$label <- ifelse(signif_data$p.value < 0.001, "",
                              ifelse(signif_data$p.value < 0.01, "",
                                     ifelse(signif_data$p.value < 0.05, "", "ns")))
  
  # Filter out non-significant comparisons
  signif_data <- signif_data[signif_data$p.value < 0.05, ]
  
  #create QQ-Sanity check
  qq<-ggplot() + geom_qq(aes(sample = rstandard(fit))) + geom_abline(color = "red") + coord_fixed() + ylab(dv.title)
  ggsave(paste0(resultdir,"/figures/","QQ_", dv, ".pdf"), qq, width = 6, height = 4, dpi = 300)  # Adjust width, height, and dpi
  
  # Create the boxplot
  p <- ggplot(data, aes_string(x ="D_ASTHMA_SEVERITYGRADE_SCREEN", y = dv)) +
    geom_boxplot(aes(fill = D_ASTHMA_SEVERITYGRADE_SCREEN), width = 0.4, outlier.shape = NA)+
    xlab("\nSeverity grade")+
    ylab(dv.title)+
    geom_jitter(color="black", size=0.4, alpha=0.9) +
    scale_x_discrete(labels=level, drop = TRUE)+
    scale_fill_manual(name= "Severity grade",labels=level, values=col, drop=T)+
    geom_signif(
      data = signif_data,
      aes(xmin = xstart, xmax = xend, annotations = label, y_position = y),
      manual = TRUE
    ) +
    theme_bw() +  # White background with black grid lines
    theme(
      axis.title.x = element_text(size = 12, face = "bold"),
      axis.title.y = element_text(size = 12, face = "bold"),
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10),
      axis.line = element_line(color = "black", size = 0.5),
      strip.text = element_text(size = 12),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      panel.grid.major = element_blank(),  # Remove major grid lines
      panel.grid.minor = element_blank(),  # Remove minor grid lines
      panel.border = element_blank(),  # Add black border
      plot.margin = margin(60, 10, 10, 10)  # Increase margins for publication
    )+coord_cartesian(clip = "off", ylim = c(0, max(data[dv], na.rm=T)))
  ggsave(paste0(resultdir,"/figures/", dv, ".pdf"), p, width = 6, height = 4, dpi = 300)  # Adjust width, height, and dpi 
}



################################################
################################################
################################################
################################################
########## ANALYSIS STARTS HERE ################

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


#Variables:
clin_2_model = list("L_FeNO_Value"=glm.nb,
                    "L_MBW_LCI"=lm,
                    "L_SPIRO_FEV1perc"=lm,
                    "L_SPIRO_FEV1_FVC"=lm, 
                    "L_SPIRO_FVCperc"=lm,
                    "L_BODY_RVperc"=lm,
                    "L_BODY_RV_TLC"=lm,
                    "L_BODY_TLCperc"=lm,
                    "S_ACQ_5_TOTAL"=lm,
                    "MS_DIFF_EOS"=lm,
                    "MB_DIFF_EOS_prec_corrected.FZB"=lm)



# read in data and do pre-processing:
meta <- read.csv( "./data/raw_data/meta.csv", sep=";", header=TRUE,stringsAsFactors = FALSE)
meta2 <- meta[meta["type"] == "measured",]
rownames(meta2)<-meta2$ID_pseudo

df_port <- read.table("./data/derivative_data/20220803_plate_normalized20220803_data_donwshift_imputed.csv", sep=" ", header=TRUE, stringsAsFactors = FALSE)

to.remove <- colnames(meta[,9:27])[!colnames(meta[,9:27]) %in% colnames(df_port)]
meta <- meta[,which(!colnames(meta) %in% to.remove)]

meta2$D_ASTHMA_SEVERITYGRADE_SCREEN[meta2$D_ASTHMA_SEVERITYGRADE_SCREEN == 9]<-0
meta2[rownames(meta2),colnames(df_port)] <- df_port[rownames(meta2),colnames(df_port)]

meta2$L_SPIRO_FEV1_FVC <- meta2$L_SPIRO_FEV1_FVC*100
meta2$L_BODY_RV_TLC<-meta2$L_BODY_RV_TLC*100
meta2$year <- as.factor(meta2$year)
meta2$plate_new <- as.factor(meta2$plate_new)
meta2$D_SEX<-as.factor(meta2$D_SEX)
meta2$D_ASTHMA_SEVERITYGRADE_SCREEN <-as.factor(meta2$D_ASTHMA_SEVERITYGRADE_SCREEN)
meta2$D_SmokingStatus<-as.factor(meta2$D_SmokingStatus)

meta2$L_FeNO_Value <- as.integer(meta2$L_FeNO_Value)
meta2$L_FeNO_LogValue <- log(meta2$L_FeNO_Value+0.01)

meta2$ICS_YN <-ifelse(as.numeric(meta2$D_ASTHMA_SEVERITYGRADE_SCREEN) > 2 | meta2$QMA_DOSE_Fluticasonequivalent > 0,1,0) %>% as.factor()
meta2$OCS_YN <- meta2$QMA_YN_SystemicSteroids %>% as.factor()


### Calculate Summary statistics:

meta2$QMA_Fluticasonequivalent_re2 <- as.factor(meta2$QMA_Fluticasonequivalent_re2)
meta2$ICShigh <- as.factor(meta2$QMA_DOSE_Fluticasonequivalent > 500)
meta2$D_SmokingStatus<-as.factor(meta2$D_SmokingStatus)

meta2$OCS_YN2 <- as.factor(meta2$OCS_mg > 0)
meta2$ICN_YN2 <- as.factor(meta2$QMA_DOSE_Fluticasonequivalent > 0)

clinical_summary <- meta2 %>%
  group_by(D_ASTHMA_SEVERITYGRADE_SCREEN) %>%
  summarise(
    across(
      where(is.numeric),
      list(
        mean = ~ mean(., na.rm = TRUE),
        sd = ~ sd(., na.rm = TRUE),
        median = ~ median(., na.rm = TRUE),
        iqrlow = ~ quantile(., na.rm = TRUE)[2],
        iqrhigh = ~ quantile(., na.rm = TRUE)[4]
      )
    ),
    across(
      where(is.factor),
      ~ list(table(.))
    )
  )


#### Clinical Variables Asthma Severity Analysis:
eps <- c(21, 0.5, 10, 5, 5, 10, 10, 10, 0.5, 5)
title<- c("FeNO [ppb]","LCI","FEV1 [% pred.]","FEV1/FVC [%]","FVC [% pred.]","RV [% pred.]",
          "RV/TLC [% pred.]", "TLC [% pred.]", "AQ5 score", "Sputum Eosinophiles [%]", "Blood Eosinophiles [%]")


formular_template <- "# ~ D_ASTHMA_SEVERITYGRADE_SCREEN + year+D_SEX+D_AGE + D_SmokingStatus + OCS_YN"

wb <- createWorkbook()

for (i in seq_along(clin_2_model)) {
  dv <- names(clin_2_model)[i]
  model <- clin_2_model[[i]]
  ep = eps[i]
  dvt= title[i]
  formula <- as.formula(sub("#", dv, formular_template))
  
  clin_to_asthma_analysis(dv, model, formula, droplevels(meta2[!is.na(meta2[dv]),]), wb, ep, dvt, results.dir)
  
}

saveWorkbook(wb, paste0(results.tables,"supplementary_table_1.xlsx"), overwrite = TRUE)

