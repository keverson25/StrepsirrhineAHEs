library(RPANDA)

setwd("~/Documents/PostdocStuff/Strepsirrhine_AHE_Manuscript/Analyses/BAMM/")


# Example Dataset ---------------------------------------------------------

data("Caprimulgidae_ClaDS2") 
# plot the mcmc chains 
par(mar=c(1,1,1,1))
plot_ClaDS_chains(Caprimulgidae_ClaDS2$sampler) 
# extract the Maxima A Posteriori for each parameter 
maps = getMAPS_ClaDS(Caprimulgidae_ClaDS2$sampler, thin = 1) 
print(paste0("sigma = ", maps[1], " ; alpha = ", maps[2], " ; epsilon = ", maps[3], " ; l_0 = ", maps[4] )) 
# plot the infered branch specific speciation rates 
plot_ClaDS_phylo(Caprimulgidae_ClaDS2$tree, maps[-(1:4)])


# Read in data ------------------------------------------------------------

StrepsirrhineTree<-read.tree("Strepsirrhini_MCMCTree/Strepsirrhine_TimeTree.tre")
plot(StrepsirrhineTree)
StrepsirrhineTree<-reorder.phylo(StrepsirrhineTree)
is.ultrametric(StrepsirrhineTree)
is.binary(StrepsirrhineTree)
min(StrepsirrhineTree$edge.length)



# Fit a clads2 model -----------------------------------------

#StrepClads0<-fit_ClaDS0(tree=StrepsirrhineTree, name=NULL)
StrepClads<-fit_ClaDS(tree=StrepsirrhineTree, sample_fraction=0.7, iterations = 100000)
plot_ClaDS_chains(StrepClads)
StrepMaps = getMAPS_ClaDS(StrepClads, thin = 1) 
print(paste0("sigma = ", StrepMaps[1], " ; alpha = ", StrepMaps[2], " ; epsilon = ", StrepMaps[3], " ; l_0 = ", StrepMaps[4] )) 
plot_ClaDS_phylo(StrepClads$tree, StrepMaps[-(1:4)])


