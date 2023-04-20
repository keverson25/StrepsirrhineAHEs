library(hisse)
library(diversitree)

setwd("~/Documents/PostdocStuff/Strepsirrhine_AHE_Manuscript/Analyses/BAMM/")

StrepsirrhineTree<-read.tree("Strepsirrhini_MCMCTree/Strepsirrhine_TimeTree.tre")
StrepsirrhineTree<-reorder.phylo(StrepsirrhineTree)
is.ultrametric(StrepsirrhineTree)
is.binary(StrepsirrhineTree)
min(StrepsirrhineTree$edge.length)
StrepsirrhineTree<-ladderize(StrepsirrhineTree)
plot(StrepsirrhineTree)


# Set up MiSSE analysis ---------------------------------------------------

turnover <- c(1)
eps <- c(0)
one.rate <- MiSSE(StrepsirrhineTree, f=0.7, turnover=turnover, eps=eps)
turnover <- c(1,2)
eps <- c(0,0)
two.rate <- MiSSE(StrepsirrhineTree, f=0.7, turnover=turnover, eps=eps)
turnover <- c(1,2,3)
eps <- c(0,0,0)
three.rate <- MiSSE(StrepsirrhineTree, f=0.7, turnover=turnover, eps=eps)
#rate classes A:D
turnover <- c(1,2,3,4)
eps <- c(0,0,0,0)
four.rate <- MiSSE(StrepsirrhineTree, f=0.7, turnover=turnover, eps=eps)
#rate classes A:E
turnover <- c(1,2,3,4,5)
eps <- c(0,0,0,0,0)
five.rate <- MiSSE(StrepsirrhineTree, f=0.7, turnover=turnover, eps=eps)



# Plot the AIC values -----------------------------------------------------
one.rate
two.rate
three.rate
four.rate
plot(x=c(1,2,3,4),y=c(one.rate$AICc, two.rate$AICc, three.rate$AICc, four.rate$AICc), type="l")
GetAICWeights(list(one.rate, two.rate, three.rate, four.rate), criterion="AICc")


# Plot the top-ranked model -----------------------------------------------
two.rate.recon<-MarginReconMiSSE(phy=StrepsirrhineTree, f=0.7,  hidden.states=2, 
								pars=two.rate$solution, n.cores=2, AIC=two.rate$AIC)

plot.misse.states(two.rate.recon, rate.param="net.div", show.tip.label=TRUE, type="phylogram",
									fsize=.25, legend="none") #could also use rate.param="speciation"


#alternatively, we could plot the model-weighted average of multiple models
misse.results.list = list()
misse.results.list[[1]] = one.rate.recon
misse.results.list[[2]] = two.rate.recon
misse.results.list[[3]] = three.rate.recon
misse.results.list[[4]] = four.rate.recon
misse.results.list[[5]] = five.rate.recon
plot.misse.states(misse.results.list, rate.param="net.div", show.tip.label=TRUE, type="phylogram",
									fsize=.25, legend="none")





# HiSSE (trait-dependent diversification) liberal scheme; all missing taxa hybridize ---------------------------------
hybTrait<-read.delim(file="../MiSSE_HiSSE/HybTraits_v1.txt", header=F)
row.names(hybTrait)<-hybTrait[,1]

### "Dull null" model described in tutorial
turnover<-c(1,1)
extinction.fraction<-c(1,1)
f<-c(1,0.6)
trans.rate.null <- TransMatMaker.old(hidden.states=FALSE)
trans.rate.null
null <- hisse(phy=StrepsirrhineTree, data=hybTrait, f=f, turnover=turnover, eps=extinction.fraction, hidden.states=FALSE, trans.rate=trans.rate.null)

### CID-2 model like in tutorial (3 different transition rates and no double transitions allowed)
turnover<-c(1,1,2,2)
extinction.fraction <- rep(1, 4)
f<-c(1,0.6)
trans.rate.CID2<-TransMatMakerHiSSE(hidden.traits=1, make.null=TRUE)
trans.rate.CID2
CID2<-hisse(phy=StrepsirrhineTree, data=hybTrait, f=f, turnover=turnover, eps=extinction.fraction, hidden.states=T, trans.rate=trans.rate.CID2)

### CID-4 model from tutorial (3 different transition rates)
turnover <- c(1, 1, 2, 2, 3, 3, 4, 4)
extinction.fraction <- rep(1, 8) 
f<-c(1,0.6)
trans.rate.CID4 <- TransMatMakerHiSSE(hidden.traits=3, make.null=TRUE)
trans.rate.CID4
CID4<-hisse(phy=StrepsirrhineTree, data=hybTrait, f=f, turnover=turnover, eps=extinction.fraction, hidden.states=T, trans.rate=trans.rate.CID4)

### BiSSE model from tutorial
turnover <- c(1,2)
extinction.fraction<-c(1,1)
f<-c(1,0.6)
trans.rate.bisse <-  TransMatMakerHiSSE(hidden.traits=0)
trans.rate.bisse
bisse<-hisse(phy=StrepsirrhineTree, data=hybTrait, f=f, turnover=turnover, eps=extinction.fraction, hidden.states=F, trans.rate=trans.rate.bisse)

### HiSSE model (from tutorial)
turnover <- c(1,2,3,4)
extinction.fraction <- rep(1, 4)
f<-c(1,0.6)
trans.rate.hisse<-TransMatMakerHiSSE(hidden.traits=1)
trans.rate.hisse
hisse<-hisse(phy=StrepsirrhineTree, data=hybTrait, f=f, turnover=turnover, eps=extinction.fraction, hidden.states=T, trans.rate=trans.rate.hisse)


###plot the top ranking model
plot(x=seq(from=1, to=5, by=1),y=c(null$AICc, CID2$AICc, CID4$AICc, bisse$AICc, hisse$AICc), type="l")
c(null$AICc, CID2$AICc, CID4$AICc, bisse$AICc, hisse$AICc)
GetAICWeights(list(null, CID2, CID4, bisse, hisse), criterion="AICc")

bisse.recon<-MarginReconHiSSE(phy=bisse$phy, data=bisse$data, f=bisse$f,  hidden.states=1, 
															pars=bisse$solution, n.cores=2, AIC=bisse$AIC)
bisse.recon$rates.mat
plot.hisse.states(bisse.recon, rate.param="net.div", show.tip.label=TRUE, type="phylogram",
									fsize=.25, legend="none") #could also use rate.param="speciation"

save.image("~/Documents/PostdocStuff/Strepsirrhine_AHE_Manuscript/Analyses/MiSSE_HiSSE/LiberalAssignments_AllUnsampledHybridize.RData")


# HiSSE (trait-dependent diversification) liberal scheme; no missing taxa hybridize ---------------------------------
hybTrait<-read.delim(file="../MiSSE_HiSSE/HybTraits_v1.txt", header=F)
row.names(hybTrait)<-hybTrait[,1]

### "Dull null" model described in tutorial
turnover<-c(1,1)
extinction.fraction<-c(1,1)
f<-c(0.6,1)
trans.rate.null <- TransMatMaker.old(hidden.states=FALSE)
trans.rate.null
null <- hisse(phy=StrepsirrhineTree, data=hybTrait, f=f, turnover=turnover, eps=extinction.fraction, hidden.states=FALSE, trans.rate=trans.rate.null)

### CID-2 model like in tutorial (3 different transition rates and no double transitions allowed)
turnover<-c(1,1,2,2)
extinction.fraction <- rep(1, 4)
f<-c(0.6,1)
trans.rate.CID2<-TransMatMakerHiSSE(hidden.traits=1, make.null=TRUE)
trans.rate.CID2
CID2<-hisse(phy=StrepsirrhineTree, data=hybTrait, f=f, turnover=turnover, eps=extinction.fraction, hidden.states=T, trans.rate=trans.rate.CID2)

### CID-4 model from tutorial (3 different transition rates)
turnover <- c(1, 1, 2, 2, 3, 3, 4, 4)
extinction.fraction <- rep(1, 8) 
f<-c(0.6,1)
trans.rate.CID4 <- TransMatMakerHiSSE(hidden.traits=3, make.null=TRUE)
trans.rate.CID4
CID4<-hisse(phy=StrepsirrhineTree, data=hybTrait, f=f, turnover=turnover, eps=extinction.fraction, hidden.states=T, trans.rate=trans.rate.CID4)

### BiSSE model from tutorial
turnover <- c(1,2)
extinction.fraction<-c(1,1)
f<-c(0.6,1)
trans.rate.bisse <-  TransMatMakerHiSSE(hidden.traits=0)
trans.rate.bisse
bisse<-hisse(phy=StrepsirrhineTree, data=hybTrait, f=f, turnover=turnover, eps=extinction.fraction, hidden.states=F, trans.rate=trans.rate.bisse)

### HiSSE model (from tutorial)
turnover <- c(1,2,3,4)
extinction.fraction <- rep(1, 4)
f<-c(0.6,1)
trans.rate.hisse<-TransMatMakerHiSSE(hidden.traits=1)
trans.rate.hisse
hisse<-hisse(phy=StrepsirrhineTree, data=hybTrait, f=f, turnover=turnover, eps=extinction.fraction, hidden.states=T, trans.rate=trans.rate.hisse)


###plot the top ranking model
plot(x=seq(from=1, to=5, by=1),y=c(null$AICc, CID2$AICc, CID4$AICc, bisse$AICc, hisse$AICc), type="l")
c(null$AICc, CID2$AICc, CID4$AICc, bisse$AICc, hisse$AICc)
GetAICWeights(list(null, CID2, CID4, bisse, hisse), criterion="AICc")

bisse.recon<-MarginReconHiSSE(phy=bisse$phy, data=bisse$data, f=bisse$f,  hidden.states=1, 
															pars=bisse$solution, n.cores=2, AIC=bisse$AIC)
bisse.recon$rates.mat
plot.hisse.states(bisse.recon, rate.param="net.div", show.tip.label=TRUE, type="phylogram",
									fsize=.25, legend="none") #could also use rate.param="speciation"

save.image("~/Documents/PostdocStuff/Strepsirrhine_AHE_Manuscript/Analyses/MiSSE_HiSSE/LiberalAssignments_NoUnsampledHybridize.RData")





# HiSSE (trait-dependent diversification) conservative scheme; all missing taxa hybridize ---------------------------------
hybTrait<-read.delim(file="../MiSSE_HiSSE/HybTraits_conservative.txt", header=F)
row.names(hybTrait)<-hybTrait[,1]

### "Dull null" model described in tutorial
turnover<-c(1,1)
extinction.fraction<-c(1,1)
f<-c(1,0.6)
trans.rate.null <- TransMatMaker.old(hidden.states=FALSE)
trans.rate.null
null <- hisse(phy=StrepsirrhineTree, data=hybTrait, f=f, turnover=turnover, eps=extinction.fraction, hidden.states=FALSE, trans.rate=trans.rate.null)

### CID-2 model like in tutorial (3 different transition rates and no double transitions allowed)
turnover<-c(1,1,2,2)
extinction.fraction <- rep(1, 4)
f<-c(1,0.6)
trans.rate.CID2<-TransMatMakerHiSSE(hidden.traits=1, make.null=TRUE)
trans.rate.CID2
CID2<-hisse(phy=StrepsirrhineTree, data=hybTrait, f=f, turnover=turnover, eps=extinction.fraction, hidden.states=T, trans.rate=trans.rate.CID2)

### CID-4 model from tutorial (3 different transition rates)
turnover <- c(1, 1, 2, 2, 3, 3, 4, 4)
extinction.fraction <- rep(1, 8) 
f<-c(1,0.6)
trans.rate.CID4 <- TransMatMakerHiSSE(hidden.traits=3, make.null=TRUE)
trans.rate.CID4
CID4<-hisse(phy=StrepsirrhineTree, data=hybTrait, f=f, turnover=turnover, eps=extinction.fraction, hidden.states=T, trans.rate=trans.rate.CID4)

### BiSSE model from tutorial
turnover <- c(1,2)
extinction.fraction<-c(1,1)
f<-c(1,0.6)
trans.rate.bisse <-  TransMatMakerHiSSE(hidden.traits=0)
trans.rate.bisse
bisse<-hisse(phy=StrepsirrhineTree, data=hybTrait, f=f, turnover=turnover, eps=extinction.fraction, hidden.states=F, trans.rate=trans.rate.bisse)

### HiSSE model (from tutorial)
turnover <- c(1,2,3,4)
extinction.fraction <- rep(1, 4)
f<-c(1,0.6)
trans.rate.hisse<-TransMatMakerHiSSE(hidden.traits=1)
trans.rate.hisse
hisse<-hisse(phy=StrepsirrhineTree, data=hybTrait, f=f, turnover=turnover, eps=extinction.fraction, hidden.states=T, trans.rate=trans.rate.hisse)


###plot the top ranking model
plot(x=seq(from=1, to=5, by=1),y=c(null$AICc, CID2$AICc, CID4$AICc, bisse$AICc, hisse$AICc), type="l")
c(null$AICc, CID2$AICc, CID4$AICc, bisse$AICc, hisse$AICc)
GetAICWeights(list(null, CID2, CID4, bisse, hisse), criterion="AICc")

bisse.recon<-MarginReconHiSSE(phy=bisse$phy, data=bisse$data, f=bisse$f,  hidden.states=1, 
															pars=bisse$solution, n.cores=2, AIC=bisse$AIC)
bisse.recon$rates.mat
plot.hisse.states(bisse.recon, rate.param="net.div", show.tip.label=TRUE, type="phylogram",
									fsize=.25, legend="none") #could also use rate.param="speciation"

save.image("~/Documents/PostdocStuff/Strepsirrhine_AHE_Manuscript/Analyses/MiSSE_HiSSE/ConservativeAssignments_AllUnsampledHybridize.RData")


# HiSSE (trait-dependent diversification) conservative scheme; no missing taxa hybridize ---------------------------------
hybTrait<-read.delim(file="../MiSSE_HiSSE/HybTraits_conservative.txt", header=F)
row.names(hybTrait)<-hybTrait[,1]

### "Dull null" model described in tutorial
turnover<-c(1,1)
extinction.fraction<-c(1,1)
f<-c(0.6,1)
trans.rate.null <- TransMatMaker.old(hidden.states=FALSE)
trans.rate.null
null <- hisse(phy=StrepsirrhineTree, data=hybTrait, f=f, turnover=turnover, eps=extinction.fraction, hidden.states=FALSE, trans.rate=trans.rate.null)

### CID-2 model like in tutorial (3 different transition rates and no double transitions allowed)
turnover<-c(1,1,2,2)
extinction.fraction <- rep(1, 4)
f<-c(0.6,1)
trans.rate.CID2<-TransMatMakerHiSSE(hidden.traits=1, make.null=TRUE)
trans.rate.CID2
CID2<-hisse(phy=StrepsirrhineTree, data=hybTrait, f=f, turnover=turnover, eps=extinction.fraction, hidden.states=T, trans.rate=trans.rate.CID2)

### CID-4 model from tutorial (3 different transition rates)
turnover <- c(1, 1, 2, 2, 3, 3, 4, 4)
extinction.fraction <- rep(1, 8) 
f<-c(0.6,1)
trans.rate.CID4 <- TransMatMakerHiSSE(hidden.traits=3, make.null=TRUE)
trans.rate.CID4
CID4<-hisse(phy=StrepsirrhineTree, data=hybTrait, f=f, turnover=turnover, eps=extinction.fraction, hidden.states=T, trans.rate=trans.rate.CID4)

### BiSSE model from tutorial
turnover <- c(1,2)
extinction.fraction<-c(1,1)
f<-c(0.6,1)
trans.rate.bisse <-  TransMatMakerHiSSE(hidden.traits=0)
trans.rate.bisse
bisse<-hisse(phy=StrepsirrhineTree, data=hybTrait, f=f, turnover=turnover, eps=extinction.fraction, hidden.states=F, trans.rate=trans.rate.bisse)

### HiSSE model (from tutorial)
turnover <- c(1,2,3,4)
extinction.fraction <- rep(1, 4)
f<-c(0.6,1)
trans.rate.hisse<-TransMatMakerHiSSE(hidden.traits=1)
trans.rate.hisse
hisse<-hisse(phy=StrepsirrhineTree, data=hybTrait, f=f, turnover=turnover, eps=extinction.fraction, hidden.states=T, trans.rate=trans.rate.hisse)


###plot the top ranking model
plot(x=seq(from=1, to=5, by=1),y=c(null$AICc, CID2$AICc, CID4$AICc, bisse$AICc, hisse$AICc), type="l")
c(null$AICc, CID2$AICc, CID4$AICc, bisse$AICc, hisse$AICc)
GetAICWeights(list(null, CID2, CID4, bisse, hisse), criterion="AICc")

bisse.recon<-MarginReconHiSSE(phy=bisse$phy, data=bisse$data, f=bisse$f,  hidden.states=1, 
															pars=bisse$solution, n.cores=2, AIC=bisse$AIC)
bisse.recon$rates.mat
plot.hisse.states(bisse.recon, rate.param="net.div", show.tip.label=TRUE, type="phylogram",
									fsize=.25, legend="none") #could also use rate.param="speciation"

save.image("~/Documents/PostdocStuff/Strepsirrhine_AHE_Manuscript/Analyses/MiSSE_HiSSE/ConservativeAssignments_NoUnsampledHybridize.RData")





# Lemurs only: HiSSE (trait-dependent diversification) conservative scheme; no missing taxa hybridize ---------------------------------
hybTrait<-read.delim(file="../MiSSE_HiSSE/HybTraits_conservative_LemursOnly.txt", header=F)
row.names(hybTrait)<-hybTrait[,1]
StrepsirrhineTree<-drop.tip(StrepsirrhineTree, tip=c(74:90))

### "Dull null" model described in tutorial
turnover<-c(1,1)
extinction.fraction<-c(1,1)
f<-c(0.6,1)
trans.rate.null <- TransMatMaker.old(hidden.states=FALSE)
trans.rate.null
null <- hisse(phy=StrepsirrhineTree, data=hybTrait, f=f, turnover=turnover, eps=extinction.fraction, hidden.states=FALSE, trans.rate=trans.rate.null)

### CID-2 model like in tutorial (3 different transition rates and no double transitions allowed)
turnover<-c(1,1,2,2)
extinction.fraction <- rep(1, 4)
f<-c(0.6,1)
trans.rate.CID2<-TransMatMakerHiSSE(hidden.traits=1, make.null=TRUE)
trans.rate.CID2
CID2<-hisse(phy=StrepsirrhineTree, data=hybTrait, f=f, turnover=turnover, eps=extinction.fraction, hidden.states=T, trans.rate=trans.rate.CID2)

### CID-4 model from tutorial (3 different transition rates)
turnover <- c(1, 1, 2, 2, 3, 3, 4, 4)
extinction.fraction <- rep(1, 8) 
f<-c(0.6,1)
trans.rate.CID4 <- TransMatMakerHiSSE(hidden.traits=3, make.null=TRUE)
trans.rate.CID4
CID4<-hisse(phy=StrepsirrhineTree, data=hybTrait, f=f, turnover=turnover, eps=extinction.fraction, hidden.states=T, trans.rate=trans.rate.CID4)

### BiSSE model from tutorial
turnover <- c(1,2)
extinction.fraction<-c(1,1)
f<-c(0.6,1)
trans.rate.bisse <-  TransMatMakerHiSSE(hidden.traits=0)
trans.rate.bisse
bisse<-hisse(phy=StrepsirrhineTree, data=hybTrait, f=f, turnover=turnover, eps=extinction.fraction, hidden.states=F, trans.rate=trans.rate.bisse)

### HiSSE model (from tutorial)
turnover <- c(1,2,3,4)
extinction.fraction <- rep(1, 4)
f<-c(0.6,1)
trans.rate.hisse<-TransMatMakerHiSSE(hidden.traits=1)
trans.rate.hisse
hisse<-hisse(phy=StrepsirrhineTree, data=hybTrait, f=f, turnover=turnover, eps=extinction.fraction, hidden.states=T, trans.rate=trans.rate.hisse)


###plot the top ranking model
plot(x=seq(from=1, to=5, by=1),y=c(null$AICc, CID2$AICc, CID4$AICc, bisse$AICc, hisse$AICc), type="l")
c(null$AICc, CID2$AICc, CID4$AICc, bisse$AICc, hisse$AICc)
GetAICWeights(list(null, CID2, CID4, bisse, hisse), criterion="AICc")

bisse.recon<-MarginReconHiSSE(phy=bisse$phy, data=bisse$data, f=bisse$f,  hidden.states=1, 
															pars=bisse$solution, n.cores=2, AIC=bisse$AIC)
bisse.recon$rates.mat
plot.hisse.states(bisse.recon, rate.param="net.div", show.tip.label=TRUE, type="phylogram",
									fsize=.25, legend="none") #could also use rate.param="speciation"

save.image("~/Documents/PostdocStuff/Strepsirrhine_AHE_Manuscript/Analyses/MiSSE_HiSSE/ConservativeAssignments_LemursOnly_NoUnsampledHybridize.RData")


# Lemurs only: HiSSE (trait-dependent diversification) liberal scheme; no missing taxa hybridize---------------------------------
hybTrait<-read.delim(file="../MiSSE_HiSSE/HybTraits_Liberal_LemursOnly.txt", header=F)
row.names(hybTrait)<-hybTrait[,1]
StrepsirrhineTree<-drop.tip(StrepsirrhineTree, tip=c(74:90))

### "Dull null" model described in tutorial
turnover<-c(1,1)
extinction.fraction<-c(1,1)
f<-c(0.6,1)
trans.rate.null <- TransMatMaker.old(hidden.states=FALSE)
trans.rate.null
null <- hisse(phy=StrepsirrhineTree, data=hybTrait, f=f, turnover=turnover, eps=extinction.fraction, hidden.states=FALSE, trans.rate=trans.rate.null)

### CID-2 model like in tutorial (3 different transition rates and no double transitions allowed)
turnover<-c(1,1,2,2)
extinction.fraction <- rep(1, 4)
f<-c(0.6,1)
trans.rate.CID2<-TransMatMakerHiSSE(hidden.traits=1, make.null=TRUE)
trans.rate.CID2
CID2<-hisse(phy=StrepsirrhineTree, data=hybTrait, f=f, turnover=turnover, eps=extinction.fraction, hidden.states=T, trans.rate=trans.rate.CID2)

### CID-4 model from tutorial (3 different transition rates)
turnover <- c(1, 1, 2, 2, 3, 3, 4, 4)
extinction.fraction <- rep(1, 8) 
f<-c(0.6,1)
trans.rate.CID4 <- TransMatMakerHiSSE(hidden.traits=3, make.null=TRUE)
trans.rate.CID4
CID4<-hisse(phy=StrepsirrhineTree, data=hybTrait, f=f, turnover=turnover, eps=extinction.fraction, hidden.states=T, trans.rate=trans.rate.CID4)

### BiSSE model from tutorial
turnover <- c(1,2)
extinction.fraction<-c(1,1)
f<-c(0.6,1)
trans.rate.bisse <-  TransMatMakerHiSSE(hidden.traits=0)
trans.rate.bisse
bisse<-hisse(phy=StrepsirrhineTree, data=hybTrait, f=f, turnover=turnover, eps=extinction.fraction, hidden.states=F, trans.rate=trans.rate.bisse)

### HiSSE model (from tutorial)
turnover <- c(1,2,3,4)
extinction.fraction <- rep(1, 4)
f<-c(0.6,1)
trans.rate.hisse<-TransMatMakerHiSSE(hidden.traits=1)
trans.rate.hisse
hisse<-hisse(phy=StrepsirrhineTree, data=hybTrait, f=f, turnover=turnover, eps=extinction.fraction, hidden.states=T, trans.rate=trans.rate.hisse)


###plot the top ranking model
plot(x=seq(from=1, to=5, by=1),y=c(null$AICc, CID2$AICc, CID4$AICc, bisse$AICc, hisse$AICc), type="l")
c(null$AICc, CID2$AICc, CID4$AICc, bisse$AICc, hisse$AICc)
GetAICWeights(list(null, CID2, CID4, bisse, hisse), criterion="AICc")

bisse.recon<-MarginReconHiSSE(phy=bisse$phy, data=bisse$data, f=bisse$f,  hidden.states=1, 
															pars=bisse$solution, n.cores=2, AIC=bisse$AIC)
bisse.recon$rates.mat
plot.hisse.states(bisse.recon, rate.param="net.div", show.tip.label=TRUE, type="phylogram",
									fsize=.25, legend="none") #could also use rate.param="speciation"

save.image("~/Documents/PostdocStuff/Strepsirrhine_AHE_Manuscript/Analyses/MiSSE_HiSSE/LiberalAssignments_LemursOnly_NoUnsampledHybridize.RData")

# Lemurs only: HiSSE (trait-dependent diversification) liberal scheme; all missing taxa hybridize---------------------------------
hybTrait<-read.delim(file="../MiSSE_HiSSE/HybTraits_Liberal_LemursOnly.txt", header=F)
row.names(hybTrait)<-hybTrait[,1]
StrepsirrhineTree<-drop.tip(StrepsirrhineTree, tip=c(74:90))

### "Dull null" model described in tutorial
turnover<-c(1,1)
extinction.fraction<-c(1,1)
f<-c(1,0.6)
trans.rate.null <- TransMatMaker.old(hidden.states=FALSE)
trans.rate.null
null <- hisse(phy=StrepsirrhineTree, data=hybTrait, f=f, turnover=turnover, eps=extinction.fraction, hidden.states=FALSE, trans.rate=trans.rate.null)

### CID-2 model like in tutorial (3 different transition rates and no double transitions allowed)
turnover<-c(1,1,2,2)
extinction.fraction <- rep(1, 4)
f<-c(1,0.6)
trans.rate.CID2<-TransMatMakerHiSSE(hidden.traits=1, make.null=TRUE)
trans.rate.CID2
CID2<-hisse(phy=StrepsirrhineTree, data=hybTrait, f=f, turnover=turnover, eps=extinction.fraction, hidden.states=T, trans.rate=trans.rate.CID2)

### CID-4 model from tutorial (3 different transition rates)
turnover <- c(1, 1, 2, 2, 3, 3, 4, 4)
extinction.fraction <- rep(1, 8) 
f<-c(1,0.6)
trans.rate.CID4 <- TransMatMakerHiSSE(hidden.traits=3, make.null=TRUE)
trans.rate.CID4
CID4<-hisse(phy=StrepsirrhineTree, data=hybTrait, f=f, turnover=turnover, eps=extinction.fraction, hidden.states=T, trans.rate=trans.rate.CID4)

### BiSSE model from tutorial
turnover <- c(1,2)
extinction.fraction<-c(1,1)
f<-c(1,0.6)
trans.rate.bisse <-  TransMatMakerHiSSE(hidden.traits=0)
trans.rate.bisse
bisse<-hisse(phy=StrepsirrhineTree, data=hybTrait, f=f, turnover=turnover, eps=extinction.fraction, hidden.states=F, trans.rate=trans.rate.bisse)

### HiSSE model (from tutorial)
turnover <- c(1,2,3,4)
extinction.fraction <- rep(1, 4)
f<-c(1,0.6)
trans.rate.hisse<-TransMatMakerHiSSE(hidden.traits=1)
trans.rate.hisse
hisse<-hisse(phy=StrepsirrhineTree, data=hybTrait, f=f, turnover=turnover, eps=extinction.fraction, hidden.states=T, trans.rate=trans.rate.hisse)


###plot the top ranking model
plot(x=seq(from=1, to=5, by=1),y=c(null$AICc, CID2$AICc, CID4$AICc, bisse$AICc, hisse$AICc), type="l")
c(null$AICc, CID2$AICc, CID4$AICc, bisse$AICc, hisse$AICc)
GetAICWeights(list(null, CID2, CID4, bisse, hisse), criterion="AICc")

bisse.recon<-MarginReconHiSSE(phy=bisse$phy, data=bisse$data, f=bisse$f,  hidden.states=1, 
															pars=bisse$solution, n.cores=2, AIC=bisse$AIC)
bisse.recon$rates.mat
plot.hisse.states(bisse.recon, rate.param="net.div", show.tip.label=TRUE, type="phylogram",
									fsize=.25, legend="none") #could also use rate.param="speciation"

save.image("~/Documents/PostdocStuff/Strepsirrhine_AHE_Manuscript/Analyses/MiSSE_HiSSE/LiberalAssignments_LemursOnly_AllUnsampledHybridize.RData")