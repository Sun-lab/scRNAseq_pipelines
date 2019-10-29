library("splatter")
library("scater")
library("ggplot2")




#splat examples
#param default
params=newSplatParams()
params
#genes
params=setParam(params, "nGenes", 5000)
sim_splat=splatSimulate(params, batchCells=rep(200,20),verbose = FALSE)
counts(sim_splat)[1:10,1:10]


#BASiC examples
library(scater)
library("BASiCS")
data("sc_example_counts")
spike.info=data.frame(Name = rownames(sc_example_counts)[1:10],
                         Input = rnorm(10, 500, 200),
                         stringsAsFactors = FALSE)
params=BASiCSEstimate(sc_example_counts[1:100, 1:30],
                         spike.info)
params
sim_BASiCS=BASiCSSimulate(params)
counts(sim_BASiCS)[1:10,1:10]

#lun2 examples
library(scater)
data("sc_example_counts")
data("sc_example_cell_info")
plates=factor(sc_example_cell_info$Mutation_Status)
params=lun2Estimate(sc_example_counts, plates, min.size = 20)
params
sim_lun2=lun2Simulate()
counts(sim_lun2)[1:10,1:10]

#lun examples
library(scater)
data("sc_example_counts")
params=lunEstimate(sc_example_counts)
params
sim_lun=lunSimulate()
counts(sim_lun)[1:10,1:10]

#scDD examples
library(scater)
library("scDD")
data("sc_example_counts")
conditions=sample(1:2, ncol(sc_example_counts), replace = TRUE)
params=scDDEstimate(sc_example_counts, conditions = conditions)
params
sim_scDD=scDDSimulate()
counts(sim_scDD)[1:10,1:10]




