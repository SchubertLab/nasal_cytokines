library(openxlsx)
library(limma)
library(dplyr)

# Function for gaussian imputation

matequal <- function(x, y) {
  is.matrix(x) && is.matrix(y) && dim(x) == dim(y) && all(x == y)
}

impute_normal <- function(object, width=0.3, downshift=1.8, seed=100) {
  
  if (!is.matrix(object)) object <- as.matrix(object)
  mx <- max(object, na.rm=TRUE)
  mn <- min(object, na.rm=TRUE)
  if (mx - mn > 20) warning("Please make sure the values are log-transformed")
  
  set.seed(seed)
  object <- apply(object, 1, function(temp) {
    temp[!is.finite(temp)] <- NA
    temp_sd <- stats::sd(temp, na.rm=TRUE)
    temp_mean <- mean(temp, na.rm=TRUE)
    shrinked_sd <- width * temp_sd   # shrink sd width
    downshifted_mean <- temp_mean - downshift * temp_sd   # shift mean of imputed values
    n_missing <- sum(is.na(temp))
    temp[is.na(temp)] <- stats::rnorm(n_missing, mean=downshifted_mean, sd=shrinked_sd)
    temp
  })
  return(object)
}

# Important data extraction 

#set wokring directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


plate.correction <- read.xlsx('../data/raw_data/plate_correction.xlsx',rowNames = T)

raw.data <- read.table('../data/raw_data/nasal_cytokines_bl_data2020-07-01.csv',stringsAsFactors = F,
                       header = F,sep = ';')

filtered.data <- raw.data[which(raw.data[,7] == 'measured'),]
rownames(filtered.data) <- filtered.data[,1]
colnames(filtered.data) <- raw.data[1,]

normal.data <- filtered.data[,8:26]
colnames(normal.data) <- raw.data[1,8:26]

prot.conc <- as.numeric(filtered.data[,5])

normal.data <- as.data.frame(t(normal.data), stringsAsFactors = FALSE)
normal.data <- as.data.frame(sapply(normal.data, as.numeric),row.names = raw.data[1,8:26])
colnames(normal.data) <- row.names(filtered.data)

meta <- read.csv("../data/raw_data/meta.csv", sep=";", header=TRUE,stringsAsFactors = FALSE)
meta2 <- meta[meta["type"] == "measured",]
rownames(meta2)<-meta2$ID_pseudo

# Normalization

annot.sample <- filtered.data[,c(1,3,4,5,6,27,28,30,33,39,45)]
colnames(annot.sample) <- c('FullRunName','Year','Plate','Total_concentration','Year_Plate','Severity','Sex','Age',
                            'Condition','Asthma_severity','Smoke')
annot.sample$Year_Plate <- as.factor(annot.sample$Year_Plate)

for (i in rownames(plate.correction)) {
  for (j in colnames(plate.correction)) {
    normal.data[j,which(annot.sample$Year_Plate == i)] <- normal.data[j,which(annot.sample$Year_Plate == i)] / plate.correction[i,j]
  }
}

# Removal of cytokines with > 20% missing data

normal.data <- log2(normal.data)
tibble(normal.data)

tibble(data.frame(cytokine=names(rowMeans(ifelse(is.na(normal.data),1,0))),missingness=rowMeans(ifelse(is.na(normal.data),1,0))))
to.delete <- names(which(rowMeans(ifelse(is.na(normal.data),1,0)) < 0.2))
normal.data <- normal.data[to.delete,]

# Missing data imputation

rownames(normal.data) <- gsub("-",replacement = "\\.",x = rownames(normal.data))

width.set <- seq(0,5,by=0.1)
downshift.set <- seq(0,5,by=0.1)

normal.data <- as.data.frame(impute_normal(object = normal.data))
tibble(normal.data)

df_port <- t(normal.data)

write.csv(df_port, file=paste0("../data/derivative_data/", 
                               format(Sys.Date(), "%Y%m%d_"),
                               "plate_normalized_data_donwshift_imputed.csv")
  
)
            
            
