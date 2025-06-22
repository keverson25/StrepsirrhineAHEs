library(hisse)
library(diversitree)
library(phytools)

StrepsirrhineTrees<-read.tree("6000TactTrees.trees")


# Pybus' Gamma on Lemur Tree with fewer lemurs (leaving at least 1 sp. per genus) --------
LemurTree<-drop.tip(StrepsirrhineTrees, tip=c(1:38))
plot(LemurTree)
plot(ltt(LemurTree[[1]]), rho=.71))

DroppableLemurs<-c("'Avahi betsileo'","'Avahi meridionalis'","'Avahi unicolor'","'Cheirogaleus shethi'","'Eulemur albifrons'","'Eulemur cinereiceps'","'Eulemur coronatus'","'Eulemur mongoz'","'Eulemur rubriventer'","'Eulemur rufifrons'","'Eulemur sanfordi'","'Hapalemur gilberti'","'Hapalemur meridionalis'","'Hapalemur occidentalis'","'Lepilemur ahmansonorum'","'Lepilemur ankaranensis'","'Lepilemur betsileo'","'Lepilemur milanoii'","'Lepilemur mustelinus'","'Lepilemur otto'","'Lepilemur petteri'","'Lepilemur randrianasoloi'","'Lepilemur ruficaudatus'","'Lepilemur septentrionalis'","'Lepilemur tymerlachsoni'","'Microcebus arnholdi'","'Microcebus berthae'","'Microcebus ganzhorni'","'Microcebus lehilahytsara'","'Microcebus margotmarshae'","'Microcebus marohita'","'Microcebus mittermeieri'","'Microcebus murinus'","'Microcebus myoxinus'","'Microcebus ravelobensis'","'Microcebus sambiranensis'","'Microcebus tavaratra'","'Propithecus coronatus'","'Propithecus diadema'","'Propithecus edwardsi'","'Propithecus perrieri'","'Propithecus tattersalli'","'Propithecus verreauxi'","Avahi_cleesei","Avahi_peyrierasi","Cheirogaleus_crossleyi","Cheirogaleus_major","Cheirogaleus_sibreei","Eulemur_collaris","Eulemur_fulvus","Eulemur_macaco","Eulemur_rufus","Lepilemur_grewcockorum","Lepilemur_hubbardorum","Lepilemur_seali","Microcebus_boraha","Microcebus_gerpi","Microcebus_griseorufus","Mirza_zaza","Lepilemur_aeeclis")


# drop 5 ------------------------------------------------------------------
LemurTree5<-drop.tip(LemurTree[[1]], tip=sample(DroppableLemurs, 5, replace=F))
LemurGamma5<-ltt(LemurTree5)$gamma
for(i in 2:6000){
  LemurTree5<-drop.tip(LemurTree[[i]], tip=sample(DroppableLemurs, 5, replace=F))
  LemurGamma5[i]<-ltt(LemurTree5)$gamma
}
median(LemurGamma5)
hist(LemurGamma5)


# drop 10 -----------------------------------------------------------------
LemurTree10<-drop.tip(LemurTree[[1]], tip=sample(DroppableLemurs, 10, replace=F))
LemurGamma10<-ltt(LemurTree10)$gamma
for(i in 2:6000){
  LemurTree10<-drop.tip(LemurTree[[i]], tip=sample(DroppableLemurs, 10, replace=F))
  LemurGamma10[i]<-ltt(LemurTree10)$gamma
}
median(LemurGamma10)
hist(LemurGamma10)

# drop 15 -----------------------------------------------------------------

LemurTree15<-drop.tip(LemurTree[[1]], tip=sample(DroppableLemurs, 15, replace=F))
LemurGamma15<-ltt(LemurTree15)$gamma
for(i in 2:6000){
  LemurTree15<-drop.tip(LemurTree[[i]], tip=sample(DroppableLemurs, 15, replace=F))
  LemurGamma15[i]<-ltt(LemurTree15)$gamma
}
median(LemurGamma15)
hist(LemurGamma15)

# drop 20 -----------------------------------------------------------------

LemurTree20<-drop.tip(LemurTree[[1]], tip=sample(DroppableLemurs, 20, replace=F))
LemurGamma20<-ltt(LemurTree20)$gamma
for(i in 2:6000){
  LemurTree20<-drop.tip(LemurTree[[i]], tip=sample(DroppableLemurs, 20, replace=F))
  LemurGamma20[i]<-ltt(LemurTree20)$gamma
}
median(LemurGamma20)
hist(LemurGamma20)


# Pybus' Gamma if you add 5 random lorisiforms with splits < 10 my ----------------------------
LorisiformTree<-drop.tip(StrepsirrhineTrees, tip=c(39:147))
for(i in 1:length(LorisiformTree)){
  LorisiformTree[[i]]$edge.length<-(LorisiformTree[[i]]$edge.length)*100
}

for(i in 1:6000){
  j=1
  success <- FALSE  # Variable to track whether the iteration was successful
  while (!success) {
    result <- try({
      LorisiformTree5<-bind.tip(LorisiformTree[[i]], tip.label=paste("n.s.",j), where=which(LorisiformTree[[i]]$tip.label==sample(x=LorisiformTree[[i]]$tip.label, size=1)), position=runif(n=1, min=0.005, max=1))
    }, silent = TRUE)
    if (inherits(result, "try-error")) {
      cat("Error encountered on iteration", j, "- retrying...\n")
    } else {
      success <- TRUE  # The iteration was successful, exit the while loop
    }
  }
      for (j in 2:5) {
    success <- FALSE  # Variable to track whether the iteration was successful
    while (!success) {
      result <- try({
        # Your code here (e.g., something that might fail)
        LorisiformTree5<-bind.tip(LorisiformTree5, tip.label=paste("n.s.",j), where=which(LorisiformTree5$tip.label==sample(x=LorisiformTree5$tip.label, size=1)), position=runif(n=1, min=0.005, max=1))
      }, silent = TRUE)
      
      if (inherits(result, "try-error")) {
        cat("Error encountered on iteration", j, "- retrying...\n")
      } else {
        success <- TRUE  # The iteration was successful, exit the while loop
      }
    }
  }
  LorisiformGamma5[i]<-ltt(LorisiformTree5)$gamma
}
median(LorisiformGamma5)
hist(LorisiformGamma5)
plot(LorisiformTree5)

# 10 lorisiforms < 10 MYA -------------------------------------------------
LorisiformGamma10<-NULL
for(i in 1:6000){
  j=1
  success <- FALSE  # Variable to track whether the iteration was successful
  while (!success) {
    result <- try({
      LorisiformTree10<-bind.tip(LorisiformTree[[i]], tip.label=paste("n.s.",j), where=which(LorisiformTree[[i]]$tip.label==sample(x=LorisiformTree[[i]]$tip.label, size=1)), position=runif(n=1, min=0.005, max=1))
    }, silent = TRUE)
    if (inherits(result, "try-error")) {
      cat("Error encountered on iteration", j, "- retrying...\n")
    } else {
      success <- TRUE  # The iteration was successful, exit the while loop
    }
  }
  for (j in 2:10) {
    success <- FALSE  # Variable to track whether the iteration was successful
    while (!success) {
      result <- try({
        # Your code here (e.g., something that might fail)
        LorisiformTree10<-bind.tip(LorisiformTree10, tip.label=paste("n.s.",j), where=which(LorisiformTree10$tip.label==sample(x=LorisiformTree10$tip.label, size=1)), position=runif(n=1, min=0.005, max=1))
      }, silent = TRUE)
      
      if (inherits(result, "try-error")) {
        cat("Error encountered on iteration", j, "- retrying...\n")
      } else {
        success <- TRUE  # The iteration was successful, exit the while loop
      }
    }
  }
  LorisiformGamma10[i]<-ltt(LorisiformTree10)$gamma
}
median(LorisiformGamma10)
hist(LorisiformGamma10)
plot(LorisiformTree10)

# 15 lorisiforms < 10 MYA -------------------------------------------------
LorisiformGamma15<-NULL
for(i in 1:6000){
  j=1
  success <- FALSE  # Variable to track whether the iteration was successful
  while (!success) {
    result <- try({
      LorisiformTree15<-bind.tip(LorisiformTree[[i]], tip.label=paste("n.s.",j), where=which(LorisiformTree[[i]]$tip.label==sample(x=LorisiformTree[[i]]$tip.label, size=1)), position=runif(n=1, min=0.005, max=1))
    }, silent = TRUE)
    if (inherits(result, "try-error")) {
      cat("Error encountered on iteration", j, "- retrying...\n")
    } else {
      success <- TRUE  # The iteration was successful, exit the while loop
    }
  }
  for (j in 2:15) {
    success <- FALSE  # Variable to track whether the iteration was successful
    while (!success) {
      result <- try({
        # Your code here (e.g., something that might fail)
        LorisiformTree15<-bind.tip(LorisiformTree15, tip.label=paste("n.s.",j), where=which(LorisiformTree15$tip.label==sample(x=LorisiformTree15$tip.label, size=1)), position=runif(n=1, min=0.005, max=1))
      }, silent = TRUE)
      
      if (inherits(result, "try-error")) {
        cat("Error encountered on iteration", j, "- retrying...\n")
      } else {
        success <- TRUE  # The iteration was successful, exit the while loop
      }
    }
  }
  LorisiformGamma15[i]<-ltt(LorisiformTree15)$gamma
}
median(LorisiformGamma15)
hist(LorisiformGamma15)
plot(LorisiformTree15)

# 20 lorisiforms < 10 MYA -------------------------------------------------
LorisiformGamma20<-NULL
for(i in 1:6000){
  j=1
  success <- FALSE  # Variable to track whether the iteration was successful
  while (!success) {
    result <- try({
      LorisiformTree20<-bind.tip(LorisiformTree[[i]], tip.label=paste("n.s.",j), where=which(LorisiformTree[[i]]$tip.label==sample(x=LorisiformTree[[i]]$tip.label, size=1)), position=runif(n=1, min=0.005, max=1))
    }, silent = TRUE)
    if (inherits(result, "try-error")) {
      cat("Error encountered on iteration", j, "- retrying...\n")
    } else {
      success <- TRUE  # The iteration was successful, exit the while loop
    }
  }
  for (j in 2:20) {
    success <- FALSE  # Variable to track whether the iteration was successful
    while (!success) {
      result <- try({
        # Your code here (e.g., something that might fail)
        LorisiformTree20<-bind.tip(LorisiformTree20, tip.label=paste("n.s.",j), where=which(LorisiformTree20$tip.label==sample(x=LorisiformTree20$tip.label, size=1)), position=runif(n=1, min=0.005, max=1))
      }, silent = TRUE)
      
      if (inherits(result, "try-error")) {
        cat("Error encountered on iteration", j, "- retrying...\n")
      } else {
        success <- TRUE  # The iteration was successful, exit the while loop
      }
    }
  }
  LorisiformGamma20[i]<-ltt(LorisiformTree20)$gamma
}
median(LorisiformGamma20)
hist(LorisiformGamma20)
plot(LorisiformTree20)


