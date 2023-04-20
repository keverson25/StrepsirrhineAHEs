library(phytools)
library(ape)

setwd("~/Documents/PostdocStuff/Strepsirrhine_AHE_Manuscript/Analyses/TACT/")
StrepsirrhineTrees<-read.tree("Strepsirrhini.tacted.AllTrees.newick.tre")

# DR metric / inverse equal splits
DRstat <- function(tree) {
	
	spRate <- function(sp, tree) {
		#get branch lengths from root to tip
		edges <- vector()
		daughterNode <- match(sp, tree$tip.label)
		while (daughterNode != (length(tree$tip.label) + 1)) {
			parentNode <- tree$edge[which(tree$edge[,2] == daughterNode), 1]
			edges <- c(edges, tree$edge.length[which(tree$edge[,1] == parentNode & tree$edge[,2] == daughterNode)])
			daughterNode <- parentNode
		}
		
		res <- sum(sapply(1:length(edges), function(x) edges[x] * (1/(2 ^ (x-1)))))
		res <- res ^ (-1)
		
		return(res)
	}
	
	rates <- unlist(lapply(tree$tip.label, function(x) spRate(x, tree)))
	names(rates) <- tree$tip.label
	
	return(rates)
}


# Loop through TACT tree replicates calculating tip DR for each sp --------

allTreesDR<-DRstat(StrepsirrhineTrees[[1]])
for(i in 2:length(StrepsirrhineTrees)){
	allTreesDR<-Map(c, allTreesDR, DRstat(StrepsirrhineTrees[[i]]))
	}
allTreesDRdf<-as.data.frame(allTreesDR)

max(colMeans(allTreesDRdf))

# Merge raw species stats into clade stats --------------------------------

Lorisidae<-c("Arctocebus_aureus","Arctocebus_calabarensis","Perodicticus_potto","Perodicticus_ibeanus","Perodicticus_edwardsi","Loris_lydekkerianus","Loris_tardigradus","Nycticebus_menagensis","Nycticebus_javanicus","Nycticebus_pygmaeus","Nycticebus_bancanus","Nycticebus_bengaliensis","Nycticebus_kayan","Nycticebus_coucang","Nycticebus_hilleri","Nycticebus_borneanus")
Galagidae<-c("Euoticus_pallidus","Euoticus_elegantulus","Galagoides_thomasi","Galagoides_kumbirensis","Galagoides_demidoff","Otolemur_garnettii","Otolemur_crassicaudatus","Otolemur_monteiri","Sciurocheirus_makandensis","Sciurocheirus_alleni","Sciurocheirus_gabonensis","Sciurocheirus_cameronensis","Galago_matschiei","Galago_moholi","Galago_gallarum","Galago_senegalensis","Paragalago_orinus","Paragalago_granti","Paragalago_rondoensis","Paragalago_zanzibaricus","Paragalago_nyasae","Paragalago_cocos")
Daubentoniidae<-c("Daubentonia_madagascariensis")
Lemuridae<-c("Varecia_rubra","Varecia_variegata","Lemur_catta","Prolemur_simus","Hapalemur_gilberti","Hapalemur_meridionalis","Hapalemur_griseus","Hapalemur_occidentalis","Hapalemur_alaotrensis","Hapalemur_aureus","Eulemur_rubriventer","Eulemur_mongoz","Eulemur_coronatus","Eulemur_macaco","Eulemur_flavifrons","Eulemur_sanfordi","Eulemur_albifrons","Eulemur_fulvus","Eulemur_cinereiceps","Eulemur_collaris","Eulemur_rufifrons","Eulemur_rufus")
Lepilemuridae<-c("Lepilemur_seali","Lepilemur_hollandorum","Lepilemur_wrightae","Lepilemur_mustelinus","Lepilemur_sahamalazensis","Lepilemur_betsileo","Lepilemur_grewcockorum","Lepilemur_otto","Lepilemur_scottorum","Lepilemur_fleuretae","Lepilemur_petteri","Lepilemur_ruficaudatus","Lepilemur_hubbardorum","Lepilemur_aeeclis","Lepilemur_randrianasoloi","Lepilemur_mittermeieri","Lepilemur_septentrionalis","Lepilemur_microdon","Lepilemur_milanoii","Lepilemur_jamesorum","Lepilemur_ankaranensis","Lepilemur_edwardsi","Lepilemur_leucopus","Lepilemur_ahmansonorum","Lepilemur_tymerlachsoni","Lepilemur_dorsalis")
Indriidae<-c("Indri_indri","Avahi_unicolor","Avahi_occidentalis","Avahi_cleesei","Avahi_laniger","Avahi_peyrierasi","Avahi_betsileo","Avahi_mooreorum","Avahi_ramanantsoavanai","Avahi_meridionalis","Propithecus_edwardsi","Propithecus_diadema","Propithecus_perrieri","Propithecus_verreauxi","Propithecus_coronatus","Propithecus_tattersalli","Propithecus_coquereli","Propithecus_candidus","Propithecus_deckenii")
Cheirogaleidae<-c("Cheirogaleus_sibreei","Cheirogaleus_major","Cheirogaleus_minusculus","Cheirogaleus_crossleyi","Cheirogaleus_medius","Cheirogaleus_shethi","Cheirogaleus_andysabini","Cheirogaleus_lavasoensis","Cheirogaleus_grovesi","Mirza_zaza","Mirza_coquereli","Allocebus_trichotis","Phaner_electromontis","Phaner_pallescens","Phaner_parienti","Phaner_furcifer","Microcebus_macarthurii","Microcebus_jollyae","Microcebus_ravelobensis","Microcebus_manitatra","Microcebus_griseorufus","Microcebus_murinus","Microcebus_ganzhorni","Microcebus_bongolavensis","Microcebus_danfossi","Microcebus_marohita","Microcebus_gerpi","Microcebus_boraha","Microcebus_mittermeieri","Microcebus_lehilahytsara","Microcebus_mamiratra","Microcebus_myoxinus","Microcebus_simmonsi","Microcebus_rufus","Microcebus_tanosi","Microcebus_tavaratra","Microcebus_arnholdi","Microcebus_sambiranensis","Microcebus_berthae","Microcebus_jonahi","Microcebus_margotmarshae")

LorisidaeDRs<-as.data.frame(unlist(allTreesDRdf[,Lorisidae]))
LorisidaeDRs<-cbind(LorisidaeDRs, rep("Lorisidae",times=length(LorisidaeDRs)))
names(LorisidaeDRs)<-c("DR","names")
ggplot(LorisidaeDRs, aes(x=DR))+geom_density(fill="gray", adjust=4)

GalagidaeDRs<-as.data.frame(unlist(allTreesDRdf[,Galagidae]))
GalagidaeDRs<-cbind(GalagidaeDRs, rep("Galagidae",times=length(GalagidaeDRs)))
names(GalagidaeDRs)<-c("DR","names")
ggplot(GalagidaeDRs, aes(x=DR))+geom_density(fill="gray", adjust=4)

LorisiformDRs<-rbind(LorisidaeDRs, GalagidaeDRs)
ggplot(LorisiformDRs, aes(x=names))+geom_density(fill="gray", adjust=4)

ggplot(LorisiformDRs, aes(x=DR, y=names, fill=names))+
	geom_density_ridges(adjust=4)+
	theme_ridges()+
	theme(
		legend.position="none",
		panel.spacing=unit(0.1, "lines"),
		strip.text.x=element_text(size=8)
	)

LemuridaeDRs<-as.data.frame(unlist(allTreesDRdf[,Lemuridae]))
LemuridaeDRs<-cbind(LemuridaeDRs, rep("Lemuridae",times=length(LemuridaeDRs)))
names(LemuridaeDRs)<-c("DR","names")
ggplot(LemuridaeDRs, aes(x=DR))+geom_density(fill="gray", adjust=4)

IndriidaeDRs<-as.data.frame(unlist(allTreesDRdf[,Indriidae]))
IndriidaeDRs<-cbind(IndriidaeDRs, rep("Indriidae",times=length(IndriidaeDRs)))
names(IndriidaeDRs)<-c("DR","names")
ggplot(IndriidaeDRs, aes(x=DR))+geom_density(fill="gray", adjust=4)

LepilemuridaeDRs<-as.data.frame(unlist(allTreesDRdf[,Lepilemuridae]))
LepilemuridaeDRs<-cbind(LepilemuridaeDRs, rep("Lepilemuridae",times=length(LepilemuridaeDRs)))
names(LepilemuridaeDRs)<-c("DR","names")
ggplot(LepilemuridaeDRs, aes(x=DR))+geom_density(fill="gray", adjust=4)

CheirogaleidaeDRs<-as.data.frame(unlist(allTreesDRdf[,Cheirogaleidae]))
CheirogaleidaeDRs<-cbind(CheirogaleidaeDRs, rep("Cheirogaleidae",times=length(CheirogaleidaeDRs)))
names(CheirogaleidaeDRs)<-c("DR","names")
ggplot(CheirogaleidaeDRs, aes(x=DR))+geom_density(fill="gray", adjust=4)

DaubentoniidaeDRs<-as.data.frame(unlist(allTreesDRdf[,Daubentoniidae]))
DaubentoniidaeDRs<-cbind(DaubentoniidaeDRs, rep("Daubentoniidae",times=length(DaubentoniidaeDRs)))
names(DaubentoniidaeDRs)<-c("DR","names")

LemuriformDRs<-rbind(LemuridaeDRs, IndriidaeDRs, LepilemuridaeDRs, CheirogaleidaeDRs)
ggplot(LemuriformDRs, aes(x=DR, y=names, fill=names))+
	geom_density_ridges(adjust=4)+
	theme_ridges()+
	theme(
		legend.position="none",
		panel.spacing=unit(0.1, "lines"),
		strip.text.x=element_text(size=8)
	)

AllStrepDRs<-rbind(GalagidaeDRs, LorisidaeDRs, LemuridaeDRs, IndriidaeDRs, LepilemuridaeDRs, CheirogaleidaeDRs, DaubentoniidaeDRs)
ggplot(AllStrepDRs, aes(x=DR, y=names, fill=names))+
	geom_density_ridges(bandwidth=5, scale=2)+
	xlim(-20,100)+
	theme_ridges()+
	theme(
		legend.position="none",
		panel.spacing=unit(0.1, "lines"),
		strip.text.x=element_text(size=8)
	)


NoAyeAyeDRs<-rbind(GalagidaeDRs, LorisidaeDRs, LemuridaeDRs, IndriidaeDRs, LepilemuridaeDRs, CheirogaleidaeDRs)
ggplot(NoAyeAyeDRs, aes(x=DR, y=names, fill=names))+
	geom_density_ridges(bandwidth=3, scale=2)+
	xlim(-10,100)+
	theme_ridges()+
	theme(
		legend.position="none",
		panel.spacing=unit(0.1, "lines"),
		strip.text.x=element_text(size=8)
	)



# Simulated Trees for Expected DR Distributions ---------------------------

#Prune the 1000 trees (with tips stochastically added w/TACT) to only Lorisidae
LorisidaeTrees<-drop.tip.multiPhylo(StrepsirrhineTrees, tip=c(Galagidae, Lemuridae, Daubentoniidae, Indriidae, Cheirogaleidae, Lepilemuridae))
#Estimate Gamma for each Lorisidae subtree in the set of TACT trees
LorisidaeGamma[1]<-yule(LorisidaeTrees[[1]])$lambda
for(i in 2:1000){
	LorisidaeGamma[i]<-yule(LorisidaeTrees[[i]])$lambda
}
LorisidaeGamma<-unlist(LorisidaeGamma)
#Simulate new trees using the same gamma and n(spp) from each Lorisidae subtree
LorisidaeSimTrees<-NULL
for(i in 1:1000){
	LorisidaeSimTrees[[i]]<-pbtree(b=LorisidaeGamma[i], d=0, n=length(Lorisidae))
}
#Estimate tip DR on each simulated tree
LorisidaeSimTreesDR<-DRstat(LorisidaeSimTrees[[1]])
for(i in 2:length(LorisidaeSimTrees)){
	LorisidaeSimTreesDR<-Map(c, LorisidaeSimTreesDR, DRstat(LorisidaeSimTrees[[i]]))
}
LorisidaeSimTreesDRdf<-as.data.frame(LorisidaeSimTreesDR)
#Create a density plot of the tip DR values for the simulated Lorisidae Trees
LorisidaeSimTreesDR<-as.data.frame(unlist(LorisidaeSimTreesDRdf))
LorisidaeSimTreesDR<-cbind(LorisidaeSimTreesDR, rep("LorisidaeSims",times=length(LorisidaeSimTreesDR)))
names(LorisidaeSimTreesDR)<-c("DR","names")
ggplot(LorisidaeSimTreesDR, aes(x=DR))+geom_density(fill="gray", adjust=4)


#Prune the 1000 trees (with tips stochastically added w/TACT) to only Galagidae
GalagidaeTrees<-drop.tip.multiPhylo(StrepsirrhineTrees, tip=c(Lorisidae, Lemuridae, Daubentoniidae, Indriidae, Cheirogaleidae, Lepilemuridae))
#Estimate Gamma for each Galagidae subtree in the set of TACT trees
GalagidaeGamma<-yule(GalagidaeTrees[[1]])$lambda
for(i in 2:1000){
	GalagidaeGamma[i]<-yule(GalagidaeTrees[[i]])$lambda
}
GalagidaeGamma<-unlist(GalagidaeGamma)
#Simulate new trees using the same gamma and n(spp) from each Galagidae subtree
GalagidaeSimTrees<-NULL
for(i in 1:1000){
	GalagidaeSimTrees[[i]]<-pbtree(b=GalagidaeGamma[i], d=0, n=length(Galagidae))
}
#Estimate tip DR on each simulated tree
GalagidaeSimTreesDR<-DRstat(GalagidaeSimTrees[[1]])
for(i in 2:length(GalagidaeSimTrees)){
	GalagidaeSimTreesDR<-Map(c, GalagidaeSimTreesDR, DRstat(GalagidaeSimTrees[[i]]))
}
GalagidaeSimTreesDRdf<-as.data.frame(GalagidaeSimTreesDR)
#Create a density plot of the tip DR values for the simulated Galagidae Trees
GalagidaeSimTreesDR<-as.data.frame(unlist(GalagidaeSimTreesDRdf))
GalagidaeSimTreesDR<-cbind(GalagidaeSimTreesDR, rep("GalagidaeSims",times=length(GalagidaeSimTreesDR)))
names(GalagidaeSimTreesDR)<-c("DR","names")
ggplot(GalagidaeSimTreesDR, aes(x=DR))+geom_density(fill="gray", adjust=4)


#Prune the 1000 trees (with tips stochastically added w/TACT) to only Lemuridae
LemuridaeTrees<-drop.tip.multiPhylo(StrepsirrhineTrees, tip=c(Lorisidae, Galagidae, Daubentoniidae, Indriidae, Cheirogaleidae, Lepilemuridae))
#Estimate Gamma for each Lemuridae subtree in the set of TACT trees
LemuridaeGamma<-yule(LemuridaeTrees[[1]])$lambda
for(i in 2:1000){
	LemuridaeGamma[i]<-yule(LemuridaeTrees[[i]])$lambda
}
LemuridaeGamma<-unlist(LemuridaeGamma)
#Simulate new trees using the same gamma and n(spp) from each Lemuridae subtree
LemuridaeSimTrees<-NULL
for(i in 1:1000){
	LemuridaeSimTrees[[i]]<-pbtree(b=LemuridaeGamma[i], d=0, n=length(Lemuridae))
}
#Estimate tip DR on each simulated tree
LemuridaeSimTreesDR<-DRstat(LemuridaeSimTrees[[1]])
for(i in 2:length(LemuridaeSimTrees)){
	LemuridaeSimTreesDR<-Map(c, LemuridaeSimTreesDR, DRstat(LemuridaeSimTrees[[i]]))
}
LemuridaeSimTreesDRdf<-as.data.frame(LemuridaeSimTreesDR)
#Create a density plot of the tip DR values for the simulated Lemuridae Trees
LemuridaeSimTreesDR<-as.data.frame(unlist(LemuridaeSimTreesDRdf))
LemuridaeSimTreesDR<-cbind(LemuridaeSimTreesDR, rep("LemuridaeSims",times=length(LemuridaeSimTreesDR)))
names(LemuridaeSimTreesDR)<-c("DR","names")
ggplot(LemuridaeSimTreesDR, aes(x=DR))+geom_density(fill="gray", adjust=4)

#Prune the 1000 trees (with tips stochastically added w/TACT) to only Indriidae
IndriidaeTrees<-drop.tip.multiPhylo(StrepsirrhineTrees, tip=c(Lorisidae, Galagidae, Daubentoniidae, Lemuridae, Cheirogaleidae, Lepilemuridae))
#Estimate Gamma for each Indriidae subtree in the set of TACT trees
IndriidaeGamma<-yule(IndriidaeTrees[[1]])$lambda
for(i in 2:1000){
	IndriidaeGamma[i]<-yule(IndriidaeTrees[[i]])$lambda
}
IndriidaeGamma<-unlist(IndriidaeGamma)
#Simulate new trees using the same gamma and n(spp) from each Indriidae subtree
IndriidaeSimTrees<-NULL
for(i in 1:1000){
	IndriidaeSimTrees[[i]]<-pbtree(b=IndriidaeGamma[i], d=0, n=length(Indriidae))
}
#Estimate tip DR on each simulated tree
IndriidaeSimTreesDR<-DRstat(IndriidaeSimTrees[[1]])
for(i in 2:length(IndriidaeSimTrees)){
	IndriidaeSimTreesDR<-Map(c, IndriidaeSimTreesDR, DRstat(IndriidaeSimTrees[[i]]))
}
IndriidaeSimTreesDRdf<-as.data.frame(IndriidaeSimTreesDR)
#Create a density plot of the tip DR values for the simulated Indriidae Trees
IndriidaeSimTreesDR<-as.data.frame(unlist(IndriidaeSimTreesDRdf))
IndriidaeSimTreesDR<-cbind(IndriidaeSimTreesDR, rep("IndriidaeSims",times=length(IndriidaeSimTreesDR)))
names(IndriidaeSimTreesDR)<-c("DR","names")
ggplot(IndriidaeSimTreesDR, aes(x=DR))+geom_density(fill="gray", adjust=4)

#Prune the 1000 trees (with tips stochastically added w/TACT) to only Lepilemuridae
LepilemuridaeTrees<-drop.tip.multiPhylo(StrepsirrhineTrees, tip=c(Lorisidae, Galagidae, Daubentoniidae, Lemuridae, Cheirogaleidae, Indriidae))
#Estimate Gamma for each Lepilemuridae subtree in the set of TACT trees
LepilemuridaeGamma<-yule(LepilemuridaeTrees[[1]])$lambda
for(i in 2:1000){
	LepilemuridaeGamma[i]<-yule(LepilemuridaeTrees[[i]])$lambda
}
LepilemuridaeGamma<-unlist(LepilemuridaeGamma)
#Simulate new trees using the same gamma and n(spp) from each Lepilemuridae subtree
LepilemuridaeSimTrees<-NULL
for(i in 1:1000){
	LepilemuridaeSimTrees[[i]]<-pbtree(b=LepilemuridaeGamma[i], d=0, n=length(Lepilemuridae))
}
#Estimate tip DR on each simulated tree
LepilemuridaeSimTreesDR<-DRstat(LepilemuridaeSimTrees[[1]])
for(i in 2:length(LepilemuridaeSimTrees)){
	LepilemuridaeSimTreesDR<-Map(c, LepilemuridaeSimTreesDR, DRstat(LepilemuridaeSimTrees[[i]]))
}
LepilemuridaeSimTreesDRdf<-as.data.frame(LepilemuridaeSimTreesDR)
#Create a density plot of the tip DR values for the simulated Lepilemuridae Trees
LepilemuridaeSimTreesDR<-as.data.frame(unlist(LepilemuridaeSimTreesDRdf))
LepilemuridaeSimTreesDR<-cbind(LepilemuridaeSimTreesDR, rep("LepilemuridaeSims",times=length(LepilemuridaeSimTreesDR)))
names(LepilemuridaeSimTreesDR)<-c("DR","names")
ggplot(LepilemuridaeSimTreesDR, aes(x=DR))+geom_density(fill="gray", adjust=4)

#Prune the 1000 trees (with tips stochastically added w/TACT) to only Cheirogaleidae
CheirogaleidaeTrees<-drop.tip.multiPhylo(StrepsirrhineTrees, tip=c(Lorisidae, Galagidae, Daubentoniidae, Lemuridae, Cheirogaleidae, Lepilemuridae))
#Estimate Gamma for each Cheirogaleidae subtree in the set of TACT trees
CheirogaleidaeGamma<-yule(CheirogaleidaeTrees[[1]])$lambda
for(i in 2:1000){
	CheirogaleidaeGamma[i]<-yule(CheirogaleidaeTrees[[i]])$lambda
}
CheirogaleidaeGamma<-unlist(CheirogaleidaeGamma)
#Simulate new trees using the same gamma and n(spp) from each Cheirogaleidae subtree
CheirogaleidaeSimTrees<-NULL
for(i in 1:1000){
	CheirogaleidaeSimTrees[[i]]<-pbtree(b=CheirogaleidaeGamma[i], d=0, n=length(Cheirogaleidae))
}
#Estimate tip DR on each simulated tree
CheirogaleidaeSimTreesDR<-DRstat(CheirogaleidaeSimTrees[[1]])
for(i in 2:length(CheirogaleidaeSimTrees)){
	CheirogaleidaeSimTreesDR<-Map(c, CheirogaleidaeSimTreesDR, DRstat(CheirogaleidaeSimTrees[[i]]))
}
CheirogaleidaeSimTreesDRdf<-as.data.frame(CheirogaleidaeSimTreesDR)
#Create a density plot of the tip DR values for the simulated Cheirogaleidae Trees
CheirogaleidaeSimTreesDR<-as.data.frame(unlist(CheirogaleidaeSimTreesDRdf))
CheirogaleidaeSimTreesDR<-cbind(CheirogaleidaeSimTreesDR, rep("CheirogaleidaeSims",times=length(CheirogaleidaeSimTreesDR)))
names(CheirogaleidaeSimTreesDR)<-c("DR","names")
ggplot(CheirogaleidaeSimTreesDR, aes(x=DR))+geom_density(fill="gray", adjust=4)

EmpSimDRs<-rbind(GalagidaeDRs, LorisidaeDRs, LemuridaeDRs, IndriidaeDRs, LepilemuridaeDRs, CheirogaleidaeDRs, GalagidaeSimTreesDR, LorisidaeSimTreesDR, LemuridaeSimTreesDR, IndriidaeSimTreesDR, LepilemuridaeSimTreesDR, CheirogaleidaeSimTreesDR)
ggplot(EmpSimDRs, aes(x=DR, y=names, fill=names))+
	geom_density_ridges(bandwidth=2, scale=2)+
	xlim(-10,100)+
	theme_ridges()+
	theme(
		legend.position="none",
		panel.spacing=unit(0.1, "lines"),
		strip.text.x=element_text(size=8)
	)

