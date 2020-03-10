if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("scDD")


library("scDD")

data(scDatEx)
class(scDatEx)
dim(scDatEx)

nDE <- 100
nDP <- 100
nDM <- 100
nDB <- 100
nEE <- 100
nEP <- 100
numSamples <- 100
seed <- 502

SD <- simulateSet(scDatEx, numSamples=numSamples,
                  nDE=nDE, nDP=nDP, nDM=nDM, nDB=nDB,
                  nEE=nEE, nEP=nEP, sd.range=c(2,2), modeFC=4, plots=FALSE,
                  random.seed=seed)

library(SingleCellExperiment)


head(rowData(SD))




