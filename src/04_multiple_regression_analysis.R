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

cytokine_regression<-function(dv, data, wb, dv.title, seed=42){
  
  set.seed(seed)  # For reproducibility
  
  #create filtered data 
  fdata<-data[!is.na(data[dv]),]
  subsets<- list(all=c(0,1,2,3),
                 mild=c(1),
                 moderat=c(2),
                 severe=c(3))
  
  tmp<-names(subsets)
  tmp[1]<-dv.title
  names(subsets)<-tmp
  
  
  
  # subset analyses
  for(name in names(subsets)){
    
    s<- subsets[[name]]
    sfdata<- subset(fdata, D_ASTHMA_SEVERITYGRADE_SCREEN %in% s)
    
    if(name == "mild"){
      
      X <- model.matrix( ~ Eotaxin.3 + G.CSF + IFN.g + IL.10 + IL.13 + IL.17 
                         + IL.1alpha + IL.1F7 + IL.24 + IL.33 + IL.4 + IL.5 
                         + IL.8 + Periostin + SCGB1A.1 + TNF.alpha + year 
                         + D_SEX + D_AGE + D_SmokingStatus + ICS_YN - 1,  # -1 removes intercept (added automatically by ridge.proj)
                         data = sfdata)
      
    }else if(name %in% c("moderat", "severe") ){
      
      X <- model.matrix( ~ Eotaxin.3 + G.CSF + IFN.g + IL.10 + IL.13 + IL.17 
                         + IL.1alpha + IL.1F7 + IL.24 + IL.33 + IL.4 + IL.5 
                         + IL.8 + Periostin + SCGB1A.1 + TNF.alpha + year 
                         + D_SEX + D_AGE + D_SmokingStatus + OCS_YN  - 1,  # -1 removes intercept (added automatically by ridge.proj)
                         data = sfdata)
    }else{
      
      X <- model.matrix( ~ Eotaxin.3 + G.CSF + IFN.g + IL.10 + IL.13 + IL.17 
                         + IL.1alpha + IL.1F7 + IL.24 + IL.33 + IL.4 + IL.5 
                         + IL.8 + Periostin + SCGB1A.1 + TNF.alpha + year 
                         + D_SEX + D_AGE + D_SmokingStatus + OCS_YN + ICS_YN - 1,  # -1 removes intercept (added automatically by ridge.proj)
                         data = sfdata)
    }
    
    #scale data to cytokines match default setting of ridge.projection
    idx_scale <- append(1:16,19)
    X[,idx_scale]<-scale(X[,idx_scale])
    
    y <- unlist(sfdata[dv])
    
    ridge_fit <- ridge.proj(X,y,
                            standardize = F,    
                            family = "gaussian",  
                            betainit = "cv lasso", 
                            suppress.grouptesting = TRUE,
                            multiplecorr.method="BH"
    )
    
    conf_intervals <- confint(ridge_fit)
    
    #is.logtransformed <- rep(grepl( "Log", dv), length(ridge_fit$bhat))
    is.logtransformed<-rep(F, length(ridge_fit$bhat))
    
    # Format and display results table
    results <- data.frame(
      Predictor = colnames(X),
      Estimate = round(ifelse(is.logtransformed,exp(ridge_fit$bhat), ridge_fit$bhat), 4),
      Std.Error = round(ifelse(is.logtransformed, exp(ridge_fit$bhat)*exp(ridge_fit$se), ridge_fit$se), 4),
      CIlow =  round(ifelse(is.logtransformed, exp(conf_intervals[, 1]), conf_intervals[, 1]), 4), 
      CIhigh =  round(ifelse(is.logtransformed, exp(conf_intervals[, 2]), conf_intervals[, 2]), 4),
      p_value = format.pval(ridge_fit$pval, digits = 4)
    )
    
#    addWorksheet(wb, sub("/", "_", name))
#    writeData(wb, sub("/", "_", name), results)
    
    wb_sheet_title<-sub("\\]",")",sub("\\[", "(", sub("/", "_", name)))
    addWorksheet(wb, wb_sheet_title)
    writeData(wb, wb_sheet_title, results)
    
    
  }
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



###############################################################################
#### Post-Selection inference for multicolinear cytokines with ################
#### ridge-regression using hdi package by Bühlmann            ################
###############################################################################
clin_2_model = list("L_FeNO_LogValue"=lm,
                    "L_MBW_LCI"=lm,
                    "S_ACQ_5_TOTAL"=lm,
                    "L_SPIRO_FEV1_FVC"=lm,
                    "L_SPIRO_FEV1perc"=lm,
                    "L_SPIRO_FVCperc"=lm,
                    "L_BODY_RVperc"=lm,
                    "L_BODY_TLCperc"=lm,
                    "L_BODY_RV_TLC"=lm
)

title<- c("FeNO","LCI","AQ5 score","FEV1/FVC","FEV1","FVC","RV", "TLC", "RV_TLC") 


for (i in seq_along(clin_2_model)){
  wb <- createWorkbook()
  
  dv <- names(clin_2_model)[i]
  dvt= title[i]
  
  cytokine_regression(dv, meta2, wb, dvt, seed=57349)
  saveWorkbook(wb, paste0(results.tables,"supplementary_table_",i+1,".xlsx"), overwrite = TRUE)
}


###############################################################################
#### Compositional Data Analysis with CRL and ride projection  ################
###############################################################################

set.seed(123)

# 1. compositional data prep

continues_pred<-c("Eotaxin.3", "G.CSF", "IFN.g", "IL.10", "IL.13", "IL.17", "IL.1alpha", "IL.1F7", "IL.24",
                  "IL.33", "IL.4", "IL.5", "IL.8", "Periostin", "SCGB1A.1", "TNF.alpha", "D_AGE")
cyts <- c("Eotaxin.3", "G.CSF", "IFN.g", "IL.10", "IL.13", "IL.17", "IL.1alpha", "IL.1F7", "IL.24",
          "IL.33", "IL.4", "IL.5", "IL.8", "Periostin", "SCGB1A.1", "TNF.alpha")
cells     <- c("MS_DIFF_MONOS", "MS_DIFF_MAKROS", "MS_DIFF_NEUT", "MS_DIFF_EOS", "MS_DIFF_LYM", "MS_DIFF_FLIMMEREPITHEL")
D <- length(cells)

P_raw<-mutate_all(meta2[cells], function(x) as.numeric(as.character(x)))
comp <- meta2[complete.cases(P_raw),]
P_raw<-as.matrix(P_raw[complete.cases(P_raw),])

# 1.1 Compute CLR outcome for cell proportions
P_compos <- acomp(P_raw)
Y_clr    <- clr(P_compos)    # n x D matrix

# 1.2 Standardize cytokine predictors
comp[, continues_pred]<- scale(comp[,continues_pred])
X_scaled <- model.matrix(
  ~ Eotaxin.3 + G.CSF + IFN.g + IL.10 + IL.13 + IL.17 + IL.1alpha + IL.1F7 + IL.24 +
    IL.33 + IL.4 + IL.5 + IL.8 + Periostin + SCGB1A.1 + TNF.alpha +
    D_SEX + D_AGE + D_SmokingStatus + year + OCS_YN + ICS_YN,
  data = comp
)[,-1]

n <- nrow(X_scaled)
p<-ncol(X_scaled)


results_list <- vector("list", length = D)
names(results_list) <- colnames(Y_clr)
for(j in seq_len(D)) {
  results_list[[j]] <- ridge.proj(
    x                        = X_scaled,
    y                        = Y_clr[, j],
    family                   = "gaussian",
    standardize              = FALSE,
    suppress.grouptesting    = TRUE,
  )
}

# Extract and combine results 
summary_df <- map_dfr(names(results_list), function(cell) {
  out <- results_list[[cell]]
  
  conf_intervals <- confint(out)
  
  tibble(
    cell_type     = cell,
    predictor      = colnames(X_scaled),
    beta_debiased = out$bhat[,1],
    ci_low        = conf_intervals[, 1],
    ci_high       = conf_intervals[, 2],
    p_value       = out$pval
  )
})

summary_df <- summary_df %>%
  group_by(predictor) %>%
  mutate(
    significance = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      p_value < 0.1   ~ ".",
      TRUE            ~ ""
    )
  ) %>%
  ungroup()


write_csv(summary_df, paste0(results.tables, "compositional_debiased_estimates_summary.csv"))

