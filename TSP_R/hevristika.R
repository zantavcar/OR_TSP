### HEVRISTICNE METODE TSP
simetricna <- function(n){ #matrika, ce predpostavimo simetricni TSP, n je st. vozlisc
  A <- matrix(rep(Inf,n*n),ncol = n,nrow = n)
  for (i in (1:(n-1))){
    for (j in ((i+1):n)){
      A[i,j] <- sample(1:10,1)
      A[j,i] <- A[i,j]
    }
  }
  return(A)
}

metrika <- function(n){ #evklidska metrika, random tocke
  A <- matrix(rep(Inf,n*n),ncol = n,nrow = n)
  for (i in (1:(n-1))){
    for (j in ((i+1):n)){
      x <- sample(1:10,2)
      y <- sample(1:10,2)
      razdalja <- sqrt((x[1]-y[1])^2+(x[2]-y[2])^2)
      A[i,j] <- razdalja
      A[j,i] <- A[i,j]
    }
  }
  return(A)
}

teza <- function(A,pot){#funkcija vrne tezo poti podano kot vektor 
  #(npr. 1231 je pot 1->2,2->3,3->1), A je matrika tez povezav
  s <- 0
  n <- length(pot)
  for (i in 1:(n-1)){
    s <- s+A[pot[i],pot[i+1]]
  }
  return(s)
}
### METODA NAJBLIŽJIH SOSEDOV

najblizji_sosed <- function(A,z){ #TSP po metodi najblizjih sosedov, funkcija vrne pot in dolzino poti
  #A je matrika tez povezav, z pa vozlisce kjer zacnemo
  Q <- 1:dim(A)[1] #mnozica tock
  zacetek <- z
  pot <- c()
  while (!length(Q)==0){
    Q <- setdiff(Q,z)
    pot <- append(pot,z)
    k <- z
    z <- intersect(which(A[k,]==min(A[k,Q]),arr.ind = TRUE),Q)[1]
  }
  
  pot <- append(pot,zacetek)
  return(pot)
}

### Najblizji sosed za razlicne zacetke

najblizji <- function(A){ #funkcija najblizjih sosedov za vse mozne zacetke, vrne matriko v prvem stolpcu je zacetek,
  #v drugem prepotovana pot, v tretjem pa teza poti
  n <- dim(A)[1]
  Resitev <- matrix(rep(NA,3*n),ncol=3,nrow=n)
  for (i in 1:n){
    Resitev[i,1] <- i
    Resitev[i,2] <- paste0(as.character(najblizji_sosed(A,i)), collapse = "-")
    Resitev[i,3] <- teza(A,najblizji_sosed(A,i))
  }
  colnames(Resitev) <- c("Zacetek","Pot","Teza")
  return(Resitev)
}

#Optimalen zacetek - potrebovali bomo pri kombinirani metodi

opt_zacetek <- function(A){
  n <- dim(A)[1]
  razdalje <- c()
  for (i in 1:n){
    razdalje <- append(razdalje,teza(A,najblizji_sosed(A,i)))
    
  }
  return(which.min(razdalje))
}

### Obrati