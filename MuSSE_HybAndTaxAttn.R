library(hisse)
library(diversitree)
library(phytools)

StrepsirrhineTree<-read.tree("Strepsirrhini_MCMCTree_result_combinedCalibrations_simplified.tre")

### Character dependent MuHiSSE model with no hidden traits
libTrait5<-read.delim(file="HybTraits_liberal_TaxAttn5.txt", header=F)
row.names(libTrait5)<-libTrait5[,1]

turnover <- c(1,2,3,4)
extinction.fraction<-c(0,0,0,0)
f<-c(1,1,1,1)
trans.rate.muhisse <-  TransMatMakerMuHiSSE(hidden.traits=0)
trans.rate.muhisse
muHisse<-MuHiSSE(phy=StrepsirrhineTree, data=libTrait5, f=f, turnover=turnover, eps=extinction.fraction, hidden.states=F, trans.rate=trans.rate.muhisse)
muHisse.recon<-MarginReconMuHiSSE(phy=muHisse$phy, data=muHisse$data, f=muHisse$f,  hidden.states=0, 
                              pars=muHisse$solution, n.cores=2, AIC=muHisse$AIC)
muHisse.recon$rates.mat


### Character dependent MuHiSSE model with no hidden traits, all unsampled hybridize
turnover <- c(1,2,3,4)
extinction.fraction<-c(0,0,0,0)
f<-c(1,1,0.6,0.6)
trans.rate.muhisse <-  TransMatMakerMuHiSSE(hidden.traits=0)
trans.rate.muhisse
muHisse_unsampHyb<-MuHiSSE(phy=StrepsirrhineTree, data=libTrait5, f=f, turnover=turnover, eps=extinction.fraction, hidden.states=F, trans.rate=trans.rate.muhisse)
muHisse_unsampHyb.recon<-MarginReconMuHiSSE(phy=muHisse_unsampHyb$phy, data=muHisse_unsampHyb$data, f=muHisse_unsampHyb$f,  hidden.states=0, 
                                  pars=muHisse_unsampHyb$solution, n.cores=2, AIC=muHisse_unsampHyb$AIC)
muHisse_unsampHyb.recon$rates.mat

### Character dependent MuHiSSE model with no hidden traits, no unsampled hybridize
turnover <- c(1,2,3,4)
extinction.fraction<-c(0,0,0,0)
f<-c(0.6,0.6,1,1)
trans.rate.muhisse <-  TransMatMakerMuHiSSE(hidden.traits=0)
trans.rate.muhisse
muHisse_unsampNo<-MuHiSSE(phy=StrepsirrhineTree, data=libTrait5, f=f, turnover=turnover, eps=extinction.fraction, hidden.states=F, trans.rate=trans.rate.muhisse)
muHisse_unsampNo.recon<-MarginReconMuHiSSE(phy=muHisse_unsampNo$phy, data=muHisse_unsampNo$data, f=muHisse_unsampNo$f,  hidden.states=0, 
                                            pars=muHisse_unsampNo$solution, n.cores=2, AIC=muHisse_unsampNo$AIC)
muHisse_unsampNo.recon$rates.mat


### Character dependent MuHiSSE model, lemurs only, with no hidden traits
LemurTree<-drop.tip(StrepsirrhineTree, tip=c(74:90))
plot(LemurTree)
lemLibTrait5<-read.delim(file="HybTraits_liberal_lemursOnly_TaxAttn5.txt", header=F)

turnover <- c(1,2,3,4)
extinction.fraction<-c(0,0,0,0)
f<-c(1,1,1,1)
trans.rate.muhisse <-  TransMatMakerMuHiSSE(hidden.traits=0)
trans.rate.muhisse
muHisse_lem<-MuHiSSE(phy=LemurTree, data=lemLibTrait5, f=f, turnover=turnover, eps=extinction.fraction, hidden.states=F, trans.rate=trans.rate.muhisse)
muHisse_lem.recon<-MarginReconMuHiSSE(phy=muHisse_lem$phy, data=muHisse_lem$data, f=muHisse_lem$f,  hidden.states=0, 
                                  pars=muHisse_lem$solution, n.cores=2, AIC=muHisse_lem$AIC)
muHisse_lem.recon$rates.mat


### Character dependent MuHiSSE model with no hidden traits, lemurs only, all unsampled hybridize
turnover <- c(1,2,3,4)
extinction.fraction<-c(0,0,0,0)
f<-c(1,1,0.6,0.6)
trans.rate.muhisse <-  TransMatMakerMuHiSSE(hidden.traits=0)
trans.rate.muhisse
muHisse_lem_unsampHyb<-MuHiSSE(phy=LemurTree, data=lemLibTrait5, f=f, turnover=turnover, eps=extinction.fraction, hidden.states=F, trans.rate=trans.rate.muhisse)
muHisse_lem_unsampHyb.recon<-MarginReconMuHiSSE(phy=muHisse_lem_unsampHyb$phy, data=muHisse_lem_unsampHyb$data, f=muHisse_lem_unsampHyb$f,  hidden.states=0, 
                                            pars=muHisse_lem_unsampHyb$solution, n.cores=2, AIC=muHisse_lem_unsampHyb$AIC)
muHisse_lem_unsampHyb.recon$rates.mat

### Character dependent MuHiSSE model with no hidden traits, lemurs only, no unsampled hybridize
turnover <- c(1,2,3,4)
extinction.fraction<-c(0,0,0,0)
f<-c(0.6,0.6,1,1)
trans.rate.muhisse <-  TransMatMakerMuHiSSE(hidden.traits=0)
trans.rate.muhisse
muHisse_lem_unsampNo<-MuHiSSE(phy=LemurTree, data=lemLibTrait5, f=f, turnover=turnover, eps=extinction.fraction, hidden.states=F, trans.rate=trans.rate.muhisse)
muHisse_lem_unsampNo.recon<-MarginReconMuHiSSE(phy=muHisse_lem_unsampNo$phy, data=muHisse_lem_unsampNo$data, f=muHisse_lem_unsampNo$f,  hidden.states=0, 
                                           pars=muHisse_lem_unsampNo$solution, n.cores=2, AIC=muHisse_lem_unsampNo$AIC)
muHisse_lem_unsampNo.recon$rates.mat