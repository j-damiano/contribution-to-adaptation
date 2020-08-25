# ROBUSTESSE PREMIER ORDRE

library(rlist)

fitness <- function(s,R) 1-R/180*s

rotation <- function(pheno,ligne,newpheno){ #pour s'affranchir des effets de bord + passer d'une ligne à l'autre de la mat
  if (newpheno<pheno){
    delta <- pheno-newpheno
    ligne <- c(ligne[-c(1:delta)],ligne[1:delta])
  }
  if (newpheno>pheno){
    delta <- newpheno-pheno
    ligne <- c(ligne[(length(ligne)-delta+1):length(ligne)],ligne[-c((length(ligne)-delta+1):length(ligne))])
  }
  return(ligne)
}


robustness <- function(individu,Pheno_value){   
  line <- individu[[2]]   
  position <- which(line != 0)  #list of accessible phenotypes
  rob <- 0
  for (i in 1:length(position)){
    rob <- rob + abs(Pheno_value[individu[[1]]]-Pheno_value[position[i]]) 
  }
  return (rob/length(position)) #robustness = mean of the distances 
  #between the current phenotype and the accessible phenotypes
}

initialisation_line <- function(npheno){
  l <- c(rep(1,18),rep(0,162))
  l <- sample(l,npheno)
  return(l)
}

individu_initialisation <- function (optimum,numberP){  #sigma : standard deviation, optimum : value of the optimum phenotype, 
  #numberP : phenotypes number
  phenotype <- 90   # optimum phenotype
  l <- initialisation_line(numberP)
  return (list(phenotype,l))
}

pop_initialisation <- function(N,sigma,optimum,numberP,Pheno_value){ #N :population size
  pop <- list()
  for (i in 1:N){
    ind <- individu_initialisation(optimum,numberP)
    pop <- list.append(pop,ind)
  }
  statable<-rbind(linetab(0,pop,Pheno_value))  #table initialisation
  statable <- as.data.frame(statable)
  colnames(statable) <- c("Generation","Mean_phenotype", "Sd_phenotype", "Median_phenotype",
                          "Quantile_0.25_phenotype", "Quantile_0.75_phenotype",
                          "Quantile_0.05_phenotype", "Quantile_0.95_phenotype","Mean_robustness", 
                          "Sd_robustness", "Median_robustness", "Quantile_0.25_robustness", 
                          "Quantile_0.75_robustness","Quantile_0.05_robustness", "Quantile_0.95_robustness")
  return(list(pop,statable))
}

linetab <- function(generation,pop,Pheno_value){
  pheno <- c()
  robus <- c()
  for (i in 1:length(pop)){
    pheno <- c(pheno,Pheno_value[pop[[i]][[1]]])
    robus <- c(robus,robustness(pop[[i]],Pheno_value))
  }
  line <- c(generation, mean(pheno),sd(pheno),median(pheno),quantile(pheno,0.25),quantile(pheno,0.75),
            quantile(pheno,0.05),quantile(pheno,0.95),mean(robus),sd(robus),median(robus),
            quantile(robus,0.25),quantile(robus,0.75),quantile(robus,0.05),quantile(robus,0.95))
  return (line)  
}


mutation_line <- function(l){
  a <- sample(1:length(l),1)
  if (l[a]==1){ #if it's 1
    l[a] <- 0 #1 becomes 0
    a <- sample(1:length(l),1) #sample to a 0
    while (l[a] != 0){
      a <- sample(1:length(l),1)
    }
    l[a] <- 1 #0 becomes 1
  }
  if (l[a]==0){
    l[a] <- 1 
    a <- sample(1:length(l),1) 
    while (l[a] != 1){
      a <- sample(1:length(l),1)
    }
    l[a] <- 0 
  }
  return(l)
}

mutation_ind <- function (l) {
  nbmut <- rbinom(1,180,1e-05)
  while (nbmut>0){
    nbmut <- nbmut-1
    l<-mutation_line(l)
  }
  return(l)
}


somme_fitness <- function(pop,Pv,s){   # fitness sum of the entire population (W : fitness list)
  sigma=0
  for (i in 1:length(pop)){
    sigma <- sigma+ fitness(s,robustness(pop[[i]],Pv))
  }
  return(sigma)
}


reproduction <- function(pop,Pv,s){
  ftot <- somme_fitness(pop,Pv,s)
  proba <- list()
  for (i in 1:length(pop)){
    proba <- list.append(proba,fitness(s,robustness(pop[[i]],Pv))/ftot)   #list of the individuals fitness
  }
  newpop <- sample(pop, size=length(pop), replace = T, prob=proba)  #Wright Fisher
  return (newpop)
}

simulation <- function (n,T,s,nom1,nom2){    #n : number of indiviuals, T : number of generation
  t=0                   #mu : mutation rate sigma : standard deviation of fitness function
  T1<-Sys.time()
  a <- s*180
  Pv <- (-89:90) #list of traits/phenotypes value 
  pop <- pop_initialisation(n,10,90,180,Pv)
  statable <- as.data.frame(pop[[2]])
  pop <-pop[[1]]
  while(t<T){
    t <- t+1
    pop <- reproduction(pop,Pv,a)  #reproduction
    for (i in 1:n){   #mutations
      pop[[i]][[2]] <- mutation_ind(pop[[i]][[2]])
    }
    if (t%%200 == 0){                             #tables
      statable<-rbind(statable,linetab(t,pop,Pv))
    }
    if (t%%25000 == 0){                             
      write.table(statable,file=nom1)
    }
  }
  T2<-Sys.time()
  print(difftime(T2,T1))
  saveRDS(object = pop, file = nom2)
  return(statable)
}
