## ----global_options, include=FALSE---------------------------------------
library(knitr)
knitr::opts_chunk$set(dpi = 100, echo= TRUE, warning=FALSE, message=FALSE, fig.align = 'center', 
                      fig.show=TRUE, fig.keep = 'all', out.width = '50%') 

## ----message = FALSE-----------------------------------------------------
library(mixOmics)

## ------------------------------------------------------------------------
#data(stemcells)

input <- function(inputfile) {
  parameters <<- read.table(inputfile, as.is=T);
  rownames(parameters) <<- parameters[,1];
  dataset <<- parameters["data", 2]
  outcome <<- parameters["outcome", 2]
  studies <<- parameters["study", 2]
  X <<- read.csv(dataset, header=TRUE)
  rownames(X) <<- X[,1]
  X <<- X[,-1]
  Y <<- read.csv(outcome, header=TRUE)
  rownames(Y) <<- Y[,1]
  Y <<- Y[,-1]
  study <<- read.csv(studies, header=TRUE)
  rownames(study) <<- study[,1]
  study <<- study[,-1]
  study <<- as.factor(study)
  print(study)
}

run <- function() {}
#the combined data set X
#X = stemcells$gene
#dim(X) 

# the outcome vector Y:  
#Y = stemcells$celltype 
#length(Y) 
#summary(Y)

# the vector indicating each independent study
#study = stemcells$study
# number of samples per study:
#summary(study)

# experimental design
#table(Y,study)

output <- function(outputfile) {
## ------------------------------------------------------------------------
mint.plsda.res.perf = mint.plsda(X = X, Y = Y, study = study, ncomp = 5)

set.seed(2543)  # for reproducible result in this example
perf.mint.plsda.cell <- perf(mint.plsda.res.perf, validation = "Mfold", folds = 5, 
                  progressBar = FALSE, auc = TRUE) 

## ------------------------------------------------------------------------
plot(perf.mint.plsda.cell, col = color.mixo(5:7))

## ------------------------------------------------------------------------
write.csv(perf.mint.plsda.cell$global.error$BER, paste(outputfile, "error", "BER", "csv", sep="."))
write.csv(perf.mint.plsda.cell$global.error$overall, paste(outputfile, "error", "overall", "csv", sep="."))
write.csv(perf.mint.plsda.cell$global.error$error.rate.class$max.dist, paste(outputfile, "error", "max", "csv", sep="."))
write.csv(perf.mint.plsda.cell$global.error$error.rate.class$centroids.dist, paste(outputfile, "error", "centroids", "csv", sep="."))
write.csv(perf.mint.plsda.cell$global.error$error.rate.class$mahalanobis.dist, paste(outputfile, "error", "mahalanobis", "csv", sep="."))

## ------------------------------------------------------------------------
write.csv(perf.mint.plsda.cell$choice.ncomp, paste(outputfile, "optimalcomponents", "pretuned", "csv", sep="."))

## ------------------------------------------------------------------------
mint.plsda.res = mint.plsda(X = X, Y = Y, study = study, ncomp = 2)
#mint.plsda.res # lists the different functions
plotIndiv(mint.plsda.res, legend = TRUE, title = 'MINT PLS-DA', 
          subtitle = 'stem cell study', ellipse = T)

## ---- eval = TRUE, include = TRUE----------------------------------------
tune.mint = tune(X = X, Y = Y, study = study, ncomp = 2, test.keepX = seq(1, 100, 1), 
method = 'mint.splsda', dist = "max.dist", progressBar = FALSE)
# tune.mint   # lists the different types of outputs

# mean error rate per component and per tested keepX value
# tune.mint$error.rate

## ------------------------------------------------------------------------
# optimal number of components
#tune.mint$choice.ncomp #tune.mint$choice.ncomp # tell us again than ncomp=1 is sufficient
write.csv(tune.mint$choice.ncomp, paste(outputfile, "optimalcomponents", "posttuned", "csv", sep="."))

# optimal keepX
write.csv(tune.mint$choice.keepX, paste(outputfile, "optimalcomponents", "keep", "csv", sep="."))

plot(tune.mint, col = color.jet(2))

## ------------------------------------------------------------------------
mint.splsda.res = mint.splsda(X = X, Y = Y, study = study, ncomp = 2,  
                              keepX = tune.mint$choice.keepX)

#mint.splsda.res   # lists useful functions that can be used with a MINT object

## ------------------------------------------------------------------------
#selectVar(mint.splsda.res, comp = 1)
write.csv(selectVar(mint.splsda.res, comp = 1)$value, paste(outputfile, "importantvalues", "csv", sep="."))

## ------------------------------------------------------------------------
plotIndiv(mint.splsda.res, study = 'global', legend = TRUE, title = 'MINT sPLS-DA', 
          subtitle = 'Global', ellipse=T)

## ------------------------------------------------------------------------
plotIndiv(mint.splsda.res, study = 'all.partial',  title = 'MINT sPLS-DA', 
          subtitle = paste("Study",1:4))

## ------------------------------------------------------------------------
plotArrow(mint.splsda.res)

## ------------------------------------------------------------------------
plotVar(mint.splsda.res, cex = 4)

## ------------------------------------------------------------------------
cim(mint.splsda.res, comp = 1, margins=c(10,5), 
    row.sideColors = color.mixo(as.numeric(Y)), row.names = FALSE,
    title = "MINT sPLS-DA, component 1")

## ------------------------------------------------------------------------
network(mint.splsda.res, color.node = c(color.mixo(1), color.mixo(2)), comp = 1,
 shape.node = c("rectangle", "circle"),
 color.edge = color.jet(50),
 lty.edge = "solid", lwd.edge = 2,
 show.edge.labels = FALSE, interactive = FALSE,
 #,save = 'jpeg',    #uncomment the following if you experience margin issues with RStudio
#name.save = network
 )

## ------------------------------------------------------------------------
plotLoadings(mint.splsda.res, contrib="max", method = 'mean', comp=1, 
             study="all.partial", legend=FALSE, title="Contribution on comp 1", 
             subtitle = paste("Study",1:4))

## ------------------------------------------------------------------------
set.seed(123)  # for reproducibility of the results
perf.mint = perf(mint.splsda.res, progressBar = FALSE, dist = 'max.dist')


#perf.mint$global.error
write.csv(perf.mint$global.error$BER, paste(outputfile, "error", "finalmodel", "BER", "csv", sep="."))
write.csv(perf.mint$global.error$overall, paste(outputfile, "error", "finalmodel", "overall", "csv", sep="."))
write.csv(perf.mint$global.error$error.rate.class$max.dist, paste(outputfile, "error", "finalmodel", "max", "csv", sep="."))

## ------------------------------------------------------------------------
plot(perf.mint, col = color.mixo(5))

## ------------------------------------------------------------------------
# we predict on study 3
ind.test = which(study == "3")
test.predict <- predict(mint.splsda.res, newdata = X[ind.test, ], dist = "max.dist",
                        study.test = factor(study[ind.test]))
Prediction <- test.predict$class$max.dist[, 2]

# the confusion table compares the real subtypes with the predicted subtypes
write.csv(get.confusion_matrix(truth = Y[ind.test],
                     predicted = Prediction), paste(outputfile, "confusionmatrix", "csv", sep="."))

## ------------------------------------------------------------------------
auc.mint.splsda = auroc(mint.splsda.res, roc.comp = 2)

## ------------------------------------------------------------------------
auc.mint.splsda = auroc(mint.splsda.res, roc.comp = 2, roc.study = '2')
}
