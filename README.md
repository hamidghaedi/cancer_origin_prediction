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
subSample <- sample_cancer[sample_cancer$sample_type == "disease tissue", ]

# projects with more than 100 samples
#proj100 <- unique(subSample$project_id)[table(subSample$project_id) > 100]
project_count <- data.frame(unclass(table(subSample$project_id)))
project_count$project_id <- rownames(project_count)
names(project_count)[1] <- "count"
# finding projects with count > 100:
proj100 <- project_count$project_id[project_count$count >= 100]

# 
metDat <- subSample[subSample$project_id %in% proj100,]


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

```


