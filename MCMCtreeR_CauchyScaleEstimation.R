library(MCMCtreeR)

streps<-read.tree("SVDQ_fixedSpeciesTree_noLabels.tre")
species<-c("Tarsius_carlito_syrichta","Saimiri_sciureus","Macaca_mulatta","Gorilla_gorilla","Homo_sapiens","Pan_troglodytes","Loris_tardigradus","Galago_moholi")
pruned.tree<-drop.tip(streps,streps$tip.label[-match(species, streps$tip.label)])

plot(pruned.tree)
nodelabels()

monophyleticGroups.user <- tipDes(pruned.tree, seq(from=9, to=15, by=1))
monophyleticGroups.user

HomoPan <- estimateCauchy(minAge = .05, maxAge = .10, 
                                  monoGroups = monophyleticGroups.user[[6]], phy = pruned.tree, 
                                  offset = 0.01, estimateScale=T,
                                  maxProb=0.95)[[1]]
HomoPan

HomoGorilla <- estimateCauchy(minAge = .10, maxAge = .18, 
                          monoGroups = monophyleticGroups.user[[5]], phy = pruned.tree, 
                          offset = 0.01, estimateScale=T,
                          maxProb=0.95)[[1]]
HomoGorilla

Catarrhini <- estimateCauchy(minAge = .21, maxAge = .34, 
                              monoGroups = monophyleticGroups.user[[4]], phy = pruned.tree, 
                              offset = 0.01, estimateScale=T,
                              maxProb=0.95)[[1]]
Catarrhini

Anthropoidea <- estimateCauchy(minAge = .34, maxAge = .56, 
                             monoGroups = monophyleticGroups.user[[3]], phy = pruned.tree, 
                             offset = 0.01, estimateScale=T,
                             maxProb=0.95)[[1]]
Anthropoidea

Lorisiformes <- estimateCauchy(minAge = .369, maxAge = .47, 
                               monoGroups = monophyleticGroups.user[[7]], phy = pruned.tree, 
                               offset = 0.01, estimateScale=T,
                               maxProb=0.95)[[1]]
Lorisiformes

Haplorhini <- estimateCauchy(minAge = .558, maxAge = .70, 
                               monoGroups = monophyleticGroups.user[[2]], phy = pruned.tree, 
                               offset = 0.01, estimateScale=T,
                               maxProb=0.95)[[1]]
Haplorhini

Primates <- estimateCauchy(minAge = .55935, maxAge = .66095, 
                             monoGroups = monophyleticGroups.user[[1]], phy = pruned.tree, 
                             offset = 0.01, estimateScale=T,
                             maxProb=0.95)[[1]]
Primates
