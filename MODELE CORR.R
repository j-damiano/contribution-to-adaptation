library(igraph)
library(rlist)
library(ggplot2)


# Fitness

fitness_1 <- function(P,sigma) exp(-(P/sigma)**2) #one peak
fitness_2 <- function(P,distance,alpha,sigma) exp(-((P+0.5*distance)/sigma)**2) + alpha * exp(-((P-0.5*distance)/sigma)**2) #two peaks
#P :phenotype, sigma : standard deviation, alpha : height of the second peak , distance : distance between the two peaks


# Robustness
robustness <- function(individu,Pheno_value){   
  transition_matrix <- individu[[3]]   
  line <- transition_matrix[individu[[1]],]  # line of the individual phenotype
  position <- which(line != 0)  #list of accessible phenotypes
  rob <- 0
  for (i in 1:length(position)){
    rob <- rob + abs(Pheno_value[individu[[1]]]-Pheno_value[position[i]]) 
  }
  return (rob/length(position)) #robustness = mean of the distances 
                                #between the current phenotype and the accessible phenotypes
}

# Initialisation
individu_initialisation <- function (sigma,optimum,numberP){  #sigma : standard deviation, optimum : value of the optimum phenotype, numberP : phenotypes number
  phenotype <- sample(c((optimum-sigma):(optimum-sigma/2),(optimum+sigma/2):(optimum+sigma)),1)   # random phenotype
  transition <- sample_degseq(rep(numberP/10,numberP))  # random network with 1/10 of 1
  transition_matrix<- as.matrix(as_adjacency_matrix(transition)) #network matrix
  return (list(phenotype, transition,transition_matrix))
}

pop_initialisation <- function(N,sigma,optimum,numberP,Pheno_value){ #N :population size
  pop <- list()
  for (i in 1:N){
    ind <- individu_initialisation(sigma,optimum,numberP)
    pop <- list.append(pop,ind)
  }
  statable<-rbind(linetab(0,pop,Pheno_value))  #table initialisation
  statable <- as.data.frame(statable)
  colnames(statable) <- c("Generation","Mean_phenotype", "Sd_phenotype", "Median_phenotype",
                          "Quantile_0.05_phenotype", "Quantile_0.95_phenotype","Mean_robustness", 
                          "Sd_robustness", "Median_robustness", "Quantile_0.05_robustness", "Quantile_0.95_robustness")
  phenotable<-rbind(pheno_counter(pop,0,numberP)) # phenotypic composition of the population table
  phenotable <- as.data.frame(phenotable)
  colnames(phenotable) <- c("Generation", 1:numberP)
  return(list(pop,statable,phenotable))
}

# Table
linetab <- function(generation,pop,Pheno_value){
  pheno <- c()
  robus <- c()
  for (i in 1:length(pop)){
    pheno <- c(pheno,Pheno_value[pop[[i]][[1]]])
    robus <- c(robus,robustness(pop[[i]],Pheno_value))
  }
  line <- c(generation, mean(pheno),sd(pheno),median(pheno),quantile(pheno,0.05),quantile(pheno,0.95),
            mean(robus),sd(robus),median(robus),quantile(robus,0.05),quantile(robus,0.95))
  return (line)  
}

pheno_counter <- function(pop,generation, numberP){
  compteur <- rep(0,numberP)
  for (i in 1:length(pop)){
    compteur[ pop[[i]][[1]][1]] <- compteur[ pop[[i]][[1]][1]] + 1
  }
  return(c(generation,compteur)  )
}

# Mutation of transition network

mutation_matrix <- function (individu, numberP){   
  nbmut <- rbinom(1,numberP^2,1e-06)  #number of mutations on the matrix (40 000 values with a mutation rate of 1e-06 )
  if (nbmut >= 1){
    individu[[2]] <-rewire(individu[[2]],keeping_degseq(niter=nbmut)) #nbmut changes in the network
    individu[[3]] <- as.matrix(as_adjacency_matrix(individu[[2]]))  
  }
  return (individu)
}



# Phenotypic mutations
mutation_pheno <- function (mu,individu){   #mu :mutation rate
  if (rbinom(1,1,mu) == 1){
    ligne <- individu[[3]][individu[[1]],]  # line of the individual phenotype
    position <- which(ligne != 0)   # list of accessible phenotypes
    p <- position[sample(1:length(position),1)]  #phenotype randomly selected (uniform law)
    individu[[1]]<- p
  }
  return (individu)
}

# Reproduction
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

# Simulation

simulation <- function (n,T,mu,sigma,nom1,nom2){    #n : number of indiviuals, T : number of generation
  t=0                   #mu : mutation rate sigma : standard deviation of fitness function
  T1<-Sys.time()
  Pv <- (-89:90) #list of traits/phenotypes value (≠ of phenotypes index) 
  W <- c(fitness_1(Pv,sigma)) # list of the fitness for each phenotype
  pop <- pop_initialisation(n,10,90,180,Pv)
  statable <- as.data.frame(pop[[2]])
  phenotable <- as.data.frame(pop[[3]])
  pop <-pop[[1]]
  while(t<T){
    t <- t+1
    pop <- reproduction(pop,W)  #reproduction
    for (i in 1:n){   #mutations
      pop[[i]]<-mutation_pheno(mu,pop[[i]])
      pop[[i]] <- mutation_matrix(pop[[i]],200) 
    }
    if (t%%200 == 0){                             #tables
      statable<-rbind(statable,linetab(t,pop,Pv))
    }
    if (t%%5000==0){
      phenotable<-rbind(phenotable,pheno_counter(pop,t,200))
    }
  }
  T2<-Sys.time()
  print(difftime(T2,T1))
  write.table(statable,file=nom1)
  write.table(phenotable,file=nom2)
}

 


