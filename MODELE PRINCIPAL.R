library(rlist)

fitness_1 <- function(P,sigma) exp(-(P/sigma)**2) #one peak


rotation <- function(pheno,ligne,newpheno){ #to be free from side effects
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
  line <- rotation(individu[[1]],line,90) #centered on the optimum
  position <- which(line != 0)  #list of accessible phenotypes
  rob <- 0
  for (i in 1:length(position)){
    rob <- rob + abs(Pheno_value[90]-Pheno_value[position[i]]) 
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
                          "Quantile_0.05_phenotype", "Quantile_0.95_phenotype","Mean_distance",
                          "Sd_distance", "Median_distance", "Quantile_0.25_distance",
                          "Quantile_0.75_distance","Quantile_0.05_distance", "Quantile_0.95_distance")
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



mutation_pheno <- function (mu,individu){   #mu :mutation rate
  if (rbinom(1,1,mu) == 1){
    ligne <- individu[[2]]  # line of the individual phenotype
    position <- which(ligne != 0)   # list of accessible phenotypes
    p <- position[sample(1:length(position),1)]  #phenotype randomly selected (uniform law)
    individu[[2]]<- rotation(individu[[1]],individu[[2]],p)
    individu[[1]]<- p
  }
  return (individu)
}

somme_fitness <- function(pop,W){   # fitness sum of the entire population (W : fitness list)
  s=0
  for (i in 1:length(pop)){
    s <- s+ W[pop[[i]][[1]]]
  }
  return(s)
}

reproduction <- function(pop,W){
  ftot <- somme_fitness(pop,W)
  proba <- list()
  for (i in 1:length(pop)){
    proba <- list.append(proba,W[pop[[i]][[1]]]/ftot)   #list of the individuals fitness
  }
  newpop <- sample(pop, size=length(pop), replace = T, prob=proba)  #Wright Fisher
  return (newpop)
}

simulation <- function (n,T,mu,sigma,nom1,nom2){    #n : number of indiviuals, T : number of generation
  t=0                   #mu : mutation rate sigma : standard deviation of fitness function
  T1<-Sys.time()
  Pv <- (-89:90) #list of traits/phenotypes value 
  W <- c(fitness_1(Pv,sigma)) # list of the fitness for each phenotype
  pop <- pop_initialisation(n,10,90,180,Pv)
  statable <- as.data.frame(pop[[2]])
  pop <-pop[[1]]
  while(t<T){
    t <- t+1
    pop <- reproduction(pop,W)  #reproduction
    for (i in 1:n){   #mutations
      pop[[i]]<-mutation_pheno(mu,pop[[i]])
      pop[[i]][[2]] <- mutation_ind(pop[[i]][[2]]) 
    }
    if (t%%200 == 0){                             #tables
      statable<-rbind(statable,linetab(t,pop,Pv))
    }
    if (t%%25000 == 0){           #backs up every 25,000 generations                  
      write.table(statable,file=nom1)
    }
  }
  T2<-Sys.time()
  print(difftime(T2,T1))
  saveRDS(object = pop, file = nom2)
}


## to continue a simulation because the equilibrium is not reached
simulation1 <- function (n,T,TT,mu,sigma,nom1,nom2,pop,statable){    #n : number of indiviuals, T : number of generation now, TT : new number of generation                  #mu : mutation rate sigma : standard deviation of fitness function
  T1<-Sys.time()
  Pv <- (-89:90) #list of traits/phenotypes value (≠ of phenotypes index) here 200 (mettre dans les paramètres ?)
  W <- c(fitness_1(Pv,sigma)) # list of the fitness for each phenotype
  while(T<TT){
    T <- T+1
    pop <- reproduction(pop,W)  #reproduction
    for (i in 1:n){   #mutations
      pop[[i]]<-mutation_pheno(mu,pop[[i]])
      pop[[i]][[2]] <- mutation_line(pop[[i]][[2]]) 
    }
    if (T%%200 == 0){                             #tables
      statable<-rbind(statable,linetab(T,pop,Pv))
    }
    if (T%%25000 == 0){           #backs up every 25,000 generations                  
      write.table(statable,file=nom1)
    }
  }
  T2<-Sys.time()
  print(difftime(T2,T1))
  saveRDS(object = pop, file = nom2)
  return(statable)
}

## two peaks

fitness_2 <- function(P,distance,alpha,sigma){
  peak1 <- fitness_1(P,sigma)  
  peak2 <- fitness_1(P,sigma)*alpha
  peak1<-rotation(90,peak1,90-0.5*distance)
  peak2<-rotation(90,peak2,90+0.5*distance)
  l <-c()
  for (i in 1:length(peak1)){
    l<-c(l,peak1[i]+peak2[i])
  }
  return(l)
} 

newpheno <- function (pheno,distance) {
  if ((pheno-0.5*distance)>=1) {
    p <- pheno-0.5*distance
  }
  else{
    p <- 180+(pheno-0.5*distance)
  }
  return(p)
}



simulation2 <- function (n,mu,sigma,distance, nom1,nom2,pop){    #n : number of indiviuals, T : number of generation now, TT : new number of generation                  #mu : mutation rate sigma : standard deviation of fitness function
  T1<-Sys.time()
  Pv <- (-89:90) #list of traits/phenotypes value (≠ of phenotypes index) 
  W <- fitness_2(Pv,distance,2,sigma) # list of the fitness for each phenotype
  for (i in 1:n){   #peak 1 "moves" with the arrival of the second peak
    p <- newpheno(pop[[i]][[1]],distance)
    pop[[i]][[2]]<-rotation(pop[[i]][[1]],pop[[i]][[2]],p)
    pop[[i]][[1]]<-p
  } 
  statable<-rbind(linetab(0,pop,Pv))  #table initialisation
  statable <- as.data.frame(statable)
  colnames(statable) <- c("Generation","Mean_phenotype", "Sd_phenotype", "Median_phenotype",
                          "Quantile_0.25_phenotype", "Quantile_0.75_phenotype",
                          "Quantile_0.05_phenotype", "Quantile_0.95_phenotype","Mean_distance",
                          "Sd_distance", "Median_distance", "Quantile_0.25_distance",
                          "Quantile_0.75_distance","Quantile_0.05_distance", "Quantile_0.95_distance")
  t=0
  pheno <- c()
  for (i in 1:length(pop)){
    pheno <- c(pheno,Pv[pop[[i]][[1]]])
  }
  med <- median(pheno)
  while(med < 0.95*distance/2 ){  #median >= 95% optimum
    t <- t+1
    pop <- reproduction(pop,W)  #reproduction
    for (i in 1:n){   #mutations
      pop[[i]]<-mutation_pheno(mu,pop[[i]])
      pop[[i]][[2]] <- mutation_line(pop[[i]][[2]]) 
    }
    pheno <- c()
    for (i in 1:length(pop)){
      pheno <- c(pheno,Pv[pop[[i]][[1]]])
    }
    med <- median(pheno)
    if (t%%10 == 0){                             #tables
      statable<-rbind(statable,linetab(t,pop,Pv))
    }
    if (t%%100 == 0){                           
      write.table(statable,file=nom1)
    }
  }
  statable<-rbind(statable,linetab(t,pop,Pv))
  TT <- t+1000
  while(t<TT){
    t <- t+1
    pop <- reproduction(pop,W)  #reproduction
    for (i in 1:n){   #mutations
      pop[[i]]<-mutation_pheno(mu,pop[[i]])
      pop[[i]][[2]] <- mutation_ind(pop[[i]][[2]]) 
    }
    if (t%%10 == 0){                             #tables
      statable<-rbind(statable,linetab(t,pop,Pv))
    }
    if (t%%100 == 0){                             
      write.table(statable,file=nom1)
    }
  }
  T2<-Sys.time()
  print(difftime(T2,T1))
  write.table(statable,file=nom1)
  saveRDS(object = pop, file = nom2)
}
