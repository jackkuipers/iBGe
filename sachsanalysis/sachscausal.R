
# preprocessing etc of the Sachs data
source("./sachssetup.R")

# load libraries
library(BiDAG)
library(Bestie)

# Use BiDAG with intervention scoring
insertSource("../usrscorefns.R", package = "BiDAG")
# Use Bestie with intervention scoring
insertSource("../usrparamfns.R", package = "Bestie")

# load causal pipeline taken and adapted from https://github.com/annlia/causalpipe
source("../Rfns/toyDAGfunctions.R")

library(data.table) # for last
library(DiagrammeR) # for making DAG plot
library(DiagrammeRsvg) ## for exporting svg for plotting to file
library(rsvg) ## for converting svg to png

inputData <- scale(data)
nDAGs <- 50
nSeeds <- 50
batch <- 100 + 1:nSeeds
labels4plot <- colnames(inputData) 
nNodes <- length(labels4plot)

for(seednumber in batch){
  timing <- proc.time()
  print(paste("Seed is", seednumber))
  sampleDAGs(inData=inputData, scoretype = "usr",
             usrpar = list(pctesttype = "bge", Imat = Imat, am = 0.1, bgremove = TRUE),
             nDigraphs=nDAGs, seed=seednumber)
  computeEffects(n=nNodes, seed=seednumber)
  print(proc.time() - timing)
}

data4plot <- loadsamples(seeds=batch, nn=nNodes)

graph2plot <- dagviz(data4plot$alldigraphs, rm_nodes = 1:6, style_mat = matrix(1, 11, 11), title_text = "")
rsvg_png(charToRaw(export_svg(graph2plot)), "SachsDAGs.png", width = 4000)

pdf("SachsEffects.pdf", width = 6, height = 6)
plotEffects(effects4plot = data4plot$alleffs, xmargs = c(0.1, 0.3), label_size = 1.5,
            sortlabs = 1:11, title_text = "")
dev.off()


