# cancer_origin_prediction
A deep learning approach to predict origin of cancer based on DNA methylation profile . The base line model reached ~97% of accuracy on predicting the origin of cancer. 


Data is coming from TCGA

```R
# DNA methylation aligned to hg38
library(TCGAbiolinks)
library(sesameData)
library(sesame)

projects <- TCGAbiolinks:::getGDCprojects()$project_id
projects <- projects[grepl('^TCGA',projects,perl=T)]

# check details of each project
pL <- list()

for(proj in 1:length(projects)){
  query_met.hg38 <- GDCquery(project= projects[proj], 
                             data.category = "DNA Methylation")
  pL[[proj]] <- getResults(query_met.hg38)
}

names(pL) <- projects

# filter retrived data 
pD<- do.call(rbind.data.frame, pL)

pD <- pD[pD$type == "methylation_beta_value" & 
           pD$sample_type == "Primary Tumor",]

# projects to drop because of < 100 tumor samples
tmp <- table(pD$project)
tmp <- names(tmp)[which(tmp > 100)]

pD <- pD[pD$project %in% tmp,]

# downloading data 
p450 <- pD[pD$platform == 'Illumina Human Methylation 450',]
p27 <-pD[pD$platform == 'Illumina Human Methylation 27',]

p27_projects <- unique(p27$project)
p450_projects <- unique(p450$project)

p27_metDownload <- function(proj){
  tmp <- p27[p27$project == proj,]
  barcode <- tmp$cases
  query_met.hg38 <- GDCquery(project= proj, 
                             data.category = "DNA Methylation",
                             platform = "Illumina Human Methylation 27",
                             data.type = "Methylation Beta Value",
                             barcode = barcode)
  GDCdownload(query_met.hg38)
  data.hg38 <- GDCprepare(query_met.hg38)
  saveRDS(data.hg38, paste0(proj, "_27.RDS"))
}

lapply(p27_projects, p27_metDownload)

# p450
p450_metDownload <- function(proj){
  tmp <- p450[p450$project == proj,]
  barcode <- tmp$cases
  query_met.hg38 <- GDCquery(project= proj, 
                             data.category = "DNA Methylation",
                             platform = "Illumina Human Methylation 450",
                             data.type = "Methylation Beta Value",
                             barcode = barcode)
  GDCdownload(query_met.hg38)
  data.hg38 <- GDCprepare(query_met.hg38)
  saveRDS(data.hg38, paste0(proj, "_450.RDS"))
}

lapply(p450_projects, p450_metDownload)

# combining the matrices

met27_files <- list.files()[grepl("27.RDS", list.files())]
metMat_27_list <- list()
for(f in 1:length(met27_files)){
  metMat_27_list[[f]] <- readRDS(met27_files[f])
  metMat_27_list[[f]] <- assay(metMat_27_list[[f]])
}

metMat27<- do.call(cbind.data.frame, metMat_27_list)


# combining 450
met450_files <- c(list.files()[grepl("450.RDS", list.files())],
                  list.files()[grepl("samples.RDS", list.files())])
metMat_450_list <- list()
for(f in 1:length(met450_files)){
  print(paste0(" file ", f, " of ", length(met450_files)))
  metMat_450_list[[f]] <- readRDS(met450_files[f])
  metMat_450_list[[f]] <- assay(metMat_450_list[[f]])
  metMat_450_list[[f]] <- metMat_450_list[[f]][rownames(metMat_450_list[[f]]) %in% rownames(metMat27),]
}

metMat450<- do.call(cbind.data.frame, metMat_450_list)


# merging files
metMat27 <- metMat27[rownames(metMat27) %in% rownames(metMat450),]

all(rownames(metMat27) %in% rownames(metMat450))

all(rownames(metMat27) == rownames(metMat450))

metMat27 <- metMat27[rownames(metMat450), ]
# reorder  metMat27 according to metMat450
metMat27 <- metMat27[match(rownames(metMat450), rownames(metMat27)),]
# to check
all(rownames(metMat27) == rownames(metMat450))

# combining dataframes
metMat <- rbind(metMat27, metMat450)

# filtering based on the mc Intyr list:
# link to article : https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-15-51

mcIntyr <- data.table::fread("/home/ghaedi/projects/def-gooding-ab/ghaedi/methyl/txt_files/main_FinalMcIntyreKeepProbesAndInfo.csv")

metMat <- metMat[rownames(metMat) %in% mcIntyr$ID,]

# dealing with NAs
r_idx <- which(rowMeans(!is.na(metMat)) > 0.80)
metMat <- metMat[r_idx, ]

# replacing Nas with row means
ind <- which(is.na(metMat), arr.ind=TRUE)
metMat[ind] <- rowMeans(metMat,  na.rm = TRUE)[ind[,1]]

# remove probes with snps and X and Y chromosome

library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)
#
# get the 450k annotation data
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
## remove probes that match to chromosome  X and Y
keep <- !(row.names(df) %in% ann450k$Name[ann450k$chr %in% c("chrX","chrY")])
table(keep)
#
#
# # probes without snp
no.snp.probe <- ann450k$Name[is.na(ann450k$Probe_rs)]
snp.probe <- ann450k[!is.na(ann450k$Probe_rs), ]
#snps with maf <= 0.05
snp5.probe <- snp.probe$Name[snp.probe$Probe_maf <= 0.05]
#
# filter dataset
# 
metMat <- metMat[row.names(metMat) %in% c(no.snp.probe, snp5.probe), ]

```
The rest of analysis is based on the what could be find in the "CNN_jupyter notebook" jupyter notebook in the repository.
