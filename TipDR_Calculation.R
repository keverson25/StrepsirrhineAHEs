library(phytools)
library(ape)
library(ggplot2)
library(ggridges)
library(epm)

StrepsirrhineTrees<-read.tree("6000TactTrees.trees")

for(i in 1:length(StrepsirrhineTrees)){
  print(max(nodeHeights(StrepsirrhineTrees[[i]])))
}

#fix scaling of edge lengths
for(i in 1:length(StrepsirrhineTrees)){
  StrepsirrhineTrees[[i]]$edge.length<-(StrepsirrhineTrees[[i]]$edge.length)*100
}
plot(StrepsirrhineTrees[[1]])

# Create lists of taxa for lemurs and lorisiforms for later analyses
Lorisiforms<-c("X.Arctocebus.aureus.","X.Arctocebus.calabarensis.","X.Perodicticus.potto.","X.Perodicticus.ibeanus.","X.Perodicticus.edwardsi.","X.Loris.lydekkerianus.","X.Loris.tardigradus.","X.Nycticebus.bengaliensis.","X.Nycticebus.bancanus.","X.Nycticebus.borneanus.","X.Nycticebus.javanicus.","X.Nycticebus.menagensis.","X.Nycticebus.pygmaeus.","X.Nycticebus.coucang.","X.Nycticebus.hilleri.","X.Nycticebus.kayan.","X.Euoticus.pallidus.","X.Euoticus.elegantulus.","X.Galagoides.thomasi.","X.Galagoides.kumbirensis.","Galagoides_demidoff","X.Otolemur.monteiri.","X.Otolemur.crassicaudatus.","X.Otolemur.garnettii.","X.Sciurocheirus.gabonensis.","X.Sciurocheirus.cameronensis.","X.Sciurocheirus.makandensis.","X.Sciurocheirus.alleni.","Galago_moholi","X.Galago.senegalensis.","Galago_gallarum","X.Galago.matschiei.","X.Paragalago.nyasae.","X.Paragalago.orinus.","X.Paragalago.granti.","Paragalago_cocos","X.Paragalago.rondoensis.","X.Paragalago.zanzibaricus.")
Lemuriforms<-c("X.Daubentonia.madagascariensis.","Varecia_rubra","X.Varecia.variegata.","X.Lemur.catta.","Prolemur_simus","X.Hapalemur.gilberti.","Hapalemur_griseus","X.Hapalemur.meridionalis.","Hapalemur_aureus","X.Hapalemur.alaotrensis.","X.Hapalemur.occidentalis.","X.Eulemur.mongoz.","X.Eulemur.rubriventer.","X.Eulemur.coronatus.","X.Eulemur.flavifrons.","Eulemur_macaco","X.Eulemur.albifrons.","X.Eulemur.sanfordi.","X.Eulemur.rufifrons.","Eulemur_fulvus","Eulemur_rufus","Eulemur_collaris","X.Eulemur.cinereiceps.","X.Indri.indri.","X.Propithecus.edwardsi.","X.Propithecus.diadema.","X.Propithecus.perrieri.","X.Propithecus.coquereli.","X.Propithecus.candidus.","X.Propithecus.tattersalli.","X.Propithecus.deckenii.","X.Propithecus.verreauxi.","X.Propithecus.coronatus.","Avahi_cleesei","X.Avahi.occidentalis.","X.Avahi.unicolor.","X.Avahi.betsileo.","Avahi_mooreorum","X.Avahi.laniger.","Avahi_peyrierasi","X.Avahi.ramanantsoavanai.","X.Avahi.meridionalis.","X.Lepilemur.scottorum.","Lepilemur_seali","X.Lepilemur.wrightae.","X.Lepilemur.fleuretae.","X.Lepilemur.mittermeieri.","X.Lepilemur.betsileo.","Lepilemur_leucopus","X.Lepilemur.mustelinus.","X.Lepilemur.septentrionalis.","X.Lepilemur.ankaranensis.","X.Lepilemur.milanoii.","X.Lepilemur.tymerlachsoni.","Lepilemur_dorsalis","Lepilemur_edwardsi","X.Lepilemur.ahmansonorum.","X.Lepilemur.otto.","X.Lepilemur.microdon.","Lepilemur_grewcockorum","X.Lepilemur.petteri.","X.Lepilemur.hollandorum.","X.Lepilemur.sahamalazensis.","X.Lepilemur.ruficaudatus.","Lepilemur_hubbardorum","X.Lepilemur.randrianasoloi.","Lepilemur_jamesorum","Lepilemur_aeeclis","X.Cheirogaleus.lavasoensis.","Cheirogaleus_sibreei","Cheirogaleus_major","Cheirogaleus_grovesi","Cheirogaleus_crossleyi","X.Cheirogaleus.minusculus.","X.Cheirogaleus.shethi.","X.Cheirogaleus.andysabini.","Cheirogaleus_medius","X.Allocebus.trichotis.","Mirza_coquereli","Mirza_zaza","X.Phaner.pallescens.","X.Phaner.electromontis.","X.Phaner.parienti.","X.Phaner.furcifer.","X.Microcebus.macarthurii.","X.Microcebus.ravelobensis.","X.Microcebus.mamiratra.","Microcebus_jollyae","Microcebus_griseorufus","X.Microcebus.murinus.","X.Microcebus.danfossi.","X.Microcebus.ganzhorni.","X.Microcebus.manitatra.","Microcebus_boraha","X.Microcebus.marohita.","Microcebus_gerpi","X.Microcebus.bongolavensis.","X.Microcebus.tavaratra.","X.Microcebus.arnholdi.","X.Microcebus.margotmarshae.","X.Microcebus.sambiranensis.","X.Microcebus.berthae.","X.Microcebus.mittermeieri.","Microcebus_rufus","X.Microcebus.myoxinus.","X.Microcebus.lehilahytsara.","X.Microcebus.tanosi.","X.Microcebus.jonahi.","X.Microcebus.simmonsi.")

# Loop through TACT tree replicates calculating tip DR for each sp --------

allTreesDR<-DRstat(StrepsirrhineTrees[[1]])
for(i in 2:length(StrepsirrhineTrees)){
	allTreesDR<-Map(c, allTreesDR, DRstat(StrepsirrhineTrees[[i]]))
	}
allTreesDRdf<-as.data.frame(allTreesDR)

max(colMeans(allTreesDRdf))
max((allTreesDRdf))

# Merge raw species stats into clade stats --------------------------------

LorisiformsDRs<-as.data.frame(unlist(allTreesDRdf[,Lorisiforms]))
LorisiformsDRs<-cbind(LorisiformsDRs, rep("Lorisiforms",times=length(LorisiformsDRs)))
names(LorisiformsDRs)<-c("DR","names")
ggplot(LorisiformsDRs, aes(x=DR))+geom_density(fill="gray", adjust=4)+xlim(0,3)+ylim(0,7)

LemuriformsDRs<-as.data.frame(unlist(allTreesDRdf[,Lemuriforms]))
LemuriformsDRs<-cbind(LemuriformsDRs, rep("Lemuriforms",times=length(LemuriformsDRs)))
names(LemuriformsDRs)<-c("DR","names")
ggplot(LemuriformsDRs, aes(x=DR))+geom_density(fill="gray", adjust=4)+xlim(0,3)+ylim(0,3)


# Simulated Trees for Expected DR Distributions ---------------------------

LorisiformTips<-c("'Arctocebus aureus'","'Arctocebus calabarensis'","'Perodicticus potto'","'Perodicticus ibeanus'","'Perodicticus edwardsi'","'Loris lydekkerianus'","'Loris tardigradus'","'Nycticebus bengaliensis'","'Nycticebus bancanus'","'Nycticebus borneanus'","'Nycticebus javanicus'","'Nycticebus menagensis'","'Nycticebus pygmaeus'","'Nycticebus coucang'","'Nycticebus hilleri'","'Nycticebus kayan'","'Euoticus pallidus'","'Euoticus elegantulus'","'Galagoides thomasi'","'Galagoides kumbirensis'","Galagoides_demidoff","'Otolemur monteiri'","'Otolemur crassicaudatus'","'Otolemur garnettii'","'Sciurocheirus gabonensis'","'Sciurocheirus cameronensis'","'Sciurocheirus makandensis'","'Sciurocheirus alleni'","Galago_moholi","'Galago senegalensis'","Galago_gallarum","'Galago matschiei'","'Paragalago nyasae'","'Paragalago orinus'","'Paragalago granti'","Paragalago_cocos","'Paragalago rondoensis'","'Paragalago zanzibaricus'")
LemuriformTips<-c("'Daubentonia madagascariensis'","Varecia_rubra","'Varecia variegata'","'Lemur catta'","Prolemur_simus","'Hapalemur gilberti'","Hapalemur_griseus","'Hapalemur meridionalis'","Hapalemur_aureus","'Hapalemur alaotrensis'","'Hapalemur occidentalis'","'Eulemur mongoz'","'Eulemur rubriventer'","'Eulemur coronatus'","'Eulemur flavifrons'","Eulemur_macaco","'Eulemur albifrons'","'Eulemur sanfordi'","'Eulemur rufifrons'","Eulemur_fulvus","Eulemur_rufus","Eulemur_collaris","'Eulemur cinereiceps'","'Indri indri'","'Propithecus edwardsi'","'Propithecus diadema'","'Propithecus perrieri'","'Propithecus coquereli'","'Propithecus candidus'","'Propithecus tattersalli'","'Propithecus deckenii'","'Propithecus verreauxi'","'Propithecus coronatus'","Avahi_cleesei","'Avahi occidentalis'","'Avahi unicolor'","'Avahi betsileo'","Avahi_mooreorum","'Avahi laniger'","Avahi_peyrierasi","'Avahi ramanantsoavanai'","'Avahi meridionalis'","'Lepilemur scottorum'","Lepilemur_seali","'Lepilemur wrightae'","'Lepilemur fleuretae'","'Lepilemur mittermeieri'","'Lepilemur betsileo'","Lepilemur_leucopus","'Lepilemur mustelinus'","'Lepilemur septentrionalis'","'Lepilemur ankaranensis'","'Lepilemur milanoii'","'Lepilemur tymerlachsoni'","Lepilemur_dorsalis","Lepilemur_edwardsi","'Lepilemur ahmansonorum'","'Lepilemur otto'","'Lepilemur microdon'","Lepilemur_grewcockorum","'Lepilemur petteri'","'Lepilemur hollandorum'","'Lepilemur sahamalazensis'","'Lepilemur ruficaudatus'","Lepilemur_hubbardorum","'Lepilemur randrianasoloi'","Lepilemur_jamesorum","Lepilemur_aeeclis","'Cheirogaleus lavasoensis'","Cheirogaleus_sibreei","Cheirogaleus_major","Cheirogaleus_grovesi","Cheirogaleus_crossleyi","'Cheirogaleus minusculus'","'Cheirogaleus shethi'","'Cheirogaleus andysabini'","Cheirogaleus_medius","'Allocebus trichotis'","Mirza_coquereli","Mirza_zaza","'Phaner pallescens'","'Phaner electromontis'","'Phaner parienti'","'Phaner furcifer'","'Microcebus macarthurii'","'Microcebus ravelobensis'","'Microcebus mamiratra'","Microcebus_jollyae","Microcebus_griseorufus","'Microcebus murinus'","'Microcebus danfossi'","'Microcebus ganzhorni'","'Microcebus manitatra'","Microcebus_boraha","'Microcebus marohita'","Microcebus_gerpi","'Microcebus bongolavensis'","'Microcebus tavaratra'","'Microcebus arnholdi'","'Microcebus margotmarshae'","'Microcebus sambiranensis'","'Microcebus berthae'","'Microcebus mittermeieri'","Microcebus_rufus","'Microcebus myoxinus'","'Microcebus lehilahytsara'","'Microcebus tanosi'","'Microcebus jonahi'","'Microcebus simmonsi'")

#create Lorisiform and Lemuriform subtrees
LorisiformTrees<-keep.tip(StrepsirrhineTrees, tip=LorisiformTips)
LemuriformsTrees<-keep.tip(StrepsirrhineTrees, tip=LemuriformTips)

#Estimate Gamma for each Lorisiform subtree in the set of TACT trees
LorisiformGamma<-NULL
LorisiformGamma[1]<-yule(LorisiformTrees[[1]])$lambda
for(i in 2:6000){
  LorisiformGamma[i]<-yule(LorisiformTrees[[i]])$lambda
}
LorisiformGamma<-unlist(LorisiformGamma)
#Simulate new trees using the same gamma and n(spp) from each Lorisiform subtree
LorisiformSimTrees<-NULL
for(i in 1:6000){
  LorisiformSimTrees[[i]]<-pbtree(b=LorisiformGamma[i], d=0, n=length(Lorisiforms))
}
#Estimate tip DR on each simulated tree
LorisiformSimTreesDR<-DRstat(LorisiformSimTrees[[1]])
for(i in 2:length(LorisiformSimTrees)){
  LorisiformSimTreesDR<-Map(c, LorisiformSimTreesDR, DRstat(LorisiformSimTrees[[i]]))
}
LorisiformSimTreesDRdf<-as.data.frame(LorisiformSimTreesDR)
#Create a density plot of the tip DR values for the simulated Lorisiform Trees
LorisiformSimTreesDR<-as.data.frame(unlist(LorisiformSimTreesDRdf))
LorisiformSimTreesDR<-cbind(LorisiformSimTreesDR, rep("LorisiformSims",times=length(LorisiformSimTreesDR)))
names(LorisiformSimTreesDR)<-c("DR","names")
ggplot(LorisiformSimTreesDR, aes(x=DR))+geom_density(fill="gray", adjust=4)+xlim(0,3)+ylim(0,7)


#Estimate Gamma for each Lemuriform subtree in the set of TACT trees
LemuriformGamma<-NULL
LemuriformGamma[1]<-yule(LemuriformsTrees[[1]])$lambda
for(i in 2:6000){
  LemuriformGamma[i]<-yule(LemuriformsTrees[[i]])$lambda
}
LemuriformGamma<-unlist(LemuriformGamma)
#Simulate new trees using the same gamma and n(spp) from each Lemuriform subtree
LemuriformSimTrees<-NULL
for(i in 1:6000){
  LemuriformSimTrees[[i]]<-pbtree(b=LemuriformGamma[i], d=0, n=length(Lemuriforms))
}
#Estimate tip DR on each simulated tree
LemuriformSimTreesDR<-DRstat(LemuriformSimTrees[[1]])
for(i in 2:length(LemuriformSimTrees)){
  LemuriformSimTreesDR<-Map(c, LemuriformSimTreesDR, DRstat(LemuriformSimTrees[[i]]))
}
LemuriformSimTreesDRdf<-as.data.frame(LemuriformSimTreesDR)
#Create a density plot of the tip DR values for the simulated Lemuriform Trees
LemuriformSimTreesDR<-as.data.frame(unlist(LemuriformSimTreesDRdf))
LemuriformSimTreesDR<-cbind(LemuriformSimTreesDR, rep("LemuriformSims",times=length(LemuriformSimTreesDR)))
names(LemuriformSimTreesDR)<-c("DR","names")
ggplot(LemuriformSimTreesDR, aes(x=DR))+geom_density(fill="gray", adjust=4)+xlim(0,3)+ylim(0,3)


# LineagesThroughTimePlots ------------------------------------------------

for(i in 1:length(LorisiformTrees)){
  LorisiformTrees[[i]]<-force.ultrametric(LorisiformTrees[[i]])
}

for(i in 1:length(LemuriformsTrees)){
  LemuriformsTrees[[i]]<-force.ultrametric(LemuriformsTrees[[i]])
}

#lineages through time
mccr<-mccr(ltt(LorisiformTrees[[i]], log.lineages=T), rho=.66, nsim=500)
LorisiformGammas<-mccr$gamma
LorisiformSimGammas<-mccr$null.gamma
for(i in 2:6000){
  mccr<-mccr(ltt(LorisiformTrees[[i]], log.lineages=T), rho=.66, nsim=500)
  LorisiformGammas[i]<-mccr$gamma
  LorisiformSimGammas<-c(LorisiformSimGammas,mccr$null.gamma)
}
hist(LorisiformSimGammas, breaks=20, xlim=c(-4,4))
hist(LorisiformGammas, breaks=20, xlim=c(-4,4))
hist(LorisiformSimGammas, ylim=c(0,1000), xlim=c(-5,8))
hist(LorisiformGammas, add=T, col="blue")

mccr<-mccr(ltt(LemuriformsTrees[[i]], log.lineages=T), rho=.66, nsim=500)
LemurGammas<-mccr$gamma
LemurSimGammas<-mccr$null.gamma
for(i in 2:6000){
  mccr<-mccr(ltt(LemuriformsTrees[[i]], log.lineages=T), rho=.66, nsim=500)
  LemurGammas[i]<-mccr$gamma
  LemurSimGammas<-c(LemurSimGammas,mccr$null.gamma)
}
hist(LemurSimGammas)
hist(LemurGammas)
hist(LemurSimGammas, ylim=c(0,500), xlim=c(-5,8))
hist(LemurGammas, add=T, col="blue")

# print DR medians per species
DRsPerSpecies<-as.data.frame(sapply(allTreesDRdf, median))
write.csv(DRsPerSpecies, file = "DRsPerSpecies_6000TactTrees.txt")

