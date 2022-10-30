# cancer_origin_prediction
A deep learning approach to predict origin of cancer based on DNA methylation profile 

Data is coming from EWAS.

```R
load("C:/Users/qaedi/OneDrive - Queen's University/Documents/dl/sample_cancer.RData")

# EDA
print(paste0("There are " , dim(sample_cancer)[1], " samples in dataset with ", dim(sample_cancer)[2], " features." ))

print(paste0("The platform is ",names(table(sample_cancer$platform)), " used to assess ", table(sample_cancer$platform), " samples."))

# unique projects in the dataset:
print(paste0("Unique projects in the dataset:", length(unique(sample_cancer$project_id))))

# selecting tumor sample from projects with more than 100 samples
sample_cancer <- data.frame(sample_cancer)
sample_cancer$toSelect <- paste0(sample_cancer$project_id, "_", sample_cancer$sample_type)

# selecting
toSelect <- unique(sample_cancer$toSelect)
toSelect <- toSelect[table(sample_cancer$toSelect) >= 100]

#
metDat <- subSample[subSample$toSelect %in% toSelect,]

```
If you like to test on the TCGA data directly, her is how to download data:

```R
# DNA methylation aligned to hg38
library(TCGAbiolinks)
library(sesameData)
library(sesame)

projects <- TCGAbiolinks:::getGDCprojects()$project_id
projects <- projects[grepl('^TCGA',projects,perl=T)]

for(proj in 1:length(projects)){
query_met.hg38 <- GDCquery(project= projects[proj], 
                           data.category = "DNA Methylation",
                           platform = "Illumina Human Methylation 450",
                           data.type = "Methylation Beta Value")
GDCdownload(query_met.hg38)
data.hg38 <- GDCprepare(query_met.hg38)
saveRDS(data.hg38, paste0(projects[proj], "_450.RDS"))
}
```
To see type and count of samples in each cohort:

```R
rds_list <- list.files()[grepl('^TCGA',list.files(), perl = T )]

sample_count <- list()
for(f in 1:length(rds_list)){
  df <- readRDS(rds_list[f])
  sample_count[[f]] <- unclass(table(substr(colnames(df),14,15)))
}
names(sample_count) <- rds_list

```
The following cohorts wont be anymore considered because of < 100 samples, and need to be moved from the current directory:

`TCGA-ACC_450.RDS` , `TCGA-CHOL_450.RDS`, `TCGA-DLBC_450.RDS`, `TCGA-KICH_450.RDS`, `TCGA-MESO_450.RDS`, `TCGA-OV_450.RDS`, `TCGA-UCS_450.RDS`, `TCGA-UVM_450.RDS`, `TCGA-READ_450.RDS`. 

So the `rds_list` needs to be updated. The next step is to select tumor samples from each cohort. This can be done by looking at element 14-15 of the sample names, as "01" indicate primry tumor. 

```R
rds_list <- list.files()[grepl('^TCGA',list.files(), perl = T )]

for(f in 1:length(rds_list)){
  df <- readRDS(rds_list[f])
  df <- df[, substr(colnames(df),14,15) == "01"]
  saveRDS(df, paste0(substr(rds_list[f],1,length(rds_list[f])-8,, "_450_tumor_samples.RDS"))
}
```
Now we can merge the objects together to make a huge matrix which is needed for down stream analysis

```R
met_list <- list.files()[grepl('*_450_tumor_samples.RDS',list.files(), perl = T )]
met <- readRDS(met_list)[5]

```

## Working with methylation data from EWAS

```shell 
# running an interactive job on ComputeCanada
salloc --time=2:0:0 --ntasks=1 --account=def-gooding-ab --mem=150G
```

Then in ```R```:

```
library(BGData)
load.BGData("cancer_methylation_v1.RData", envir = parent.frame())

# subset of probes from Mc Intyre paper:
mcIn <- data.table::fread("main_FinalMcIntyreKeepProbesAndInfo.csv")

#subsetting the methylation matrix
metMat <- cancer_download[which(rownames(cancer_download) %in% mcIn$ID), ]
metMat <- metMat[, which(colnames(metMat) %in% metDat$sample_id)]

# drop probes with more than 80% NA;

r_idx <- which(rowMeans(!is.na(metMat)) > 0.80)
metMat <- metMat[r_idx, ]

rN <- rownames(metMat_2)[r_idx]

# impute values for NAs
metMat <- sapply(metMat, as.numeric)

r_idx <- which(rowMeans(!is.na(metMat)) > 0.80)
metMat <- metMat[r_idx, ]

rN <- rN[r_idx]

# writing huge table

fwrite(metMat, file = "metMat.csv", quote = "auto",sep = ",",  row.names = TRUE, col.names = TRUE)


# replacing Nas with row means
ind <- which(is.na(metMat), arr.ind=TRUE)
metMat[ind] <- rowMeans(metMat,  na.rm = TRUE)[ind[,1]]

# reverting rownames
rownames(metMat) <- rN



# splitting the dataframe into chuncks
chuncks_list <- str(split(metMat, (as.numeric(nrow(metMat))-1) %/% 10000))

df_list <- split(df, factor(sort(rank(row.names(df))%%50)))