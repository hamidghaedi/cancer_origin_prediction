# Cancer Origin Prediction Using DNA Methylation Profiles

## Background

This project employs a convolutional neural network (CNN) to classify the origin of cancers based on DNA methylation profiles. The analysis achieves a high accuracy of ~97%, utilizing data from The Cancer Genome Atlas (TCGA). The project demonstrates the power of deep learning in cancer genomics, providing a robust framework for tumor origin classification and supporting precision oncology efforts.

### Data Processing and Preprocessing

1.	Data Source:
DNA methylation data aligned to the hg38 genome was retrieved using the TCGAbiolinks library. Two platforms were considered:
•	Illumina Human Methylation 450
•	Illumina Human Methylation 27

2.	Sample Selection and Filtering:
•	Primary tumor samples with at least 100 cases per project were selected.
•	Methylation probes were curated to exclude those on sex chromosomes (chrX and chrY) and probes with minor allele frequency (MAF) > 0.05.
•	Missing data was imputed using row means for robustness.

3.	Probe Annotation:
•	Annotations were incorporated using the IlluminaHumanMethylation450kanno.ilmn12.hg19 package to ensure the biological relevance of selected probes.


### Deep Learning Methodology

1.	Model Architecture:
A sequential neural network model was implemented using Keras. The structure includes:

•	Input Layer: Accepts 17,380 features corresponding to methylation probe values.
•	Hidden Layers: Three dense layers with ReLU activation to capture complex patterns.
•	Output Layer: Softmax activation outputs class probabilities for 25 cancer types.

3.	Training:
•	The model was trained using the categorical cross-entropy loss function and optimized with the Adam optimizer.
•	A 10-fold cross-validation strategy was applied to evaluate model performance.

4.	Performance Metrics:
•	Accuracy scores from cross-validation indicate the model's robust predictive capability with minimal variance.


### Next ups:

**Dealing with Batch Effects:** Apply methods like ComBat or quantile normalization to address systematic biases across datasets.

**Incorporating Additional Data:** Integrate clinical variables such as tumor stage and grade to improve model precision.

**Model Fine-Tuning:** Optimize model architecture and hyperparameters to enhance performance and reduce overfitting.

**Adding Another Data Modality:** Explore integrating RNA-seq or DNA-seq data to create a multi-omics model for cancer origin prediction.


## Analysis

### data acquisition and processing
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
### CNN development
Also available in the  "CNN_jupyter notebook" jupyter notebook in the repository.

```python
# importing libraries
import pandas as pd
import pyreadr
import pandas
from keras.models import Sequential
from keras.layers import Dense
from scikeras.wrappers import KerasClassifier
from keras.utils import np_utils
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import KFold
from sklearn.preprocessing import LabelEncoder
from sklearn.pipeline import Pipeline
```
```python
# reading data inot dataframes
metMat = pyreadr.read_r('../metMat.RDS') 
pD = pyreadr.read_r('../pD.RDS') 
metMat = metMat[None]
pD = pD[None]
# convert to array
met = metMat.values
p = pD[['project']].values.ravel()
# defining X and Y
X = met.astype(float)
Y = p
```
```python
#Encode the Output Variable

# encode class values as integers
encoder = LabelEncoder()
encoder.fit(Y)
encoded_Y = encoder.transform(Y)
# convert integers to dummy variables (i.e. one hot encoded)
dummy_y = np_utils.to_categorical(encoded_Y)
dummy_y
```
`
array([[0., 1., 0., ..., 0., 0., 0.],
       [0., 1., 0., ..., 0., 0., 0.],
       [0., 1., 0., ..., 0., 0., 0.],
       ...,
       [0., 0., 0., ..., 0., 0., 1.],
       [0., 0., 0., ..., 0., 0., 1.],
       [0., 0., 0., ..., 0., 0., 1.]], dtype=float32)
       `

```python
# define baseline model
def baseline_model():
    # create model
    model = Sequential()
    model.add(Dense(100, input_dim=17380, activation='relu'))
    model.add(Dense(64,  activation='relu'))
    model.add(Dense(64,  activation='relu'))
    model.add(Dense(25, activation='softmax'))
    # Compile model
    model.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['accuracy'])
    return model
```
```python
...
estimator = KerasClassifier(build_fn=baseline_model, epochs=200, batch_size=1000, verbose=1)
...
kfold = KFold(n_splits=10, shuffle=True)
...
results = cross_val_score(estimator, X, dummy_y, cv=kfold)
print("Baseline: %.2f%% (%.2f%%)" % (results.mean()*100, results.std()*100))
```

`
results
array([0.96586345, 0.94678715, 0.95582329, 0.9688755 , 0.96084337,
       0.95582329, 0.95883534, 0.96582915, 0.96582915, 0.96582915])
`
