require(dplyr)
require(lpSolve) #resevanje Madzarske metode
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

najblizji_sosed <- function(A,z){ #TSP po metodi najblizjih sosedov, funkcija vrne pot
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
  zacetek <- c()
  pot <- c()
  teze <- c()
  for (i in 1:n){
    zacetek <- append(zacetek,i)
    pot <- append(pot,paste0(as.character(najblizji_sosed(A,i)), collapse = "-"))
    teze <- append(teze,teza(A,najblizji_sosed(A,i)))
  }
  tabela <- data_frame(zacetek,pot,teze)
  return(tabela)
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

### SUBTOUR REVERSAL HEURISTIC

reversal <- function(A,zacetna=najblizji_sosed(A,opt_zacetek(A))){ #kombinirana metoda, zacnemo s potjo, ki 
  #je optimalna pri metodi najblizjih sosedov
  n <- length(zacetna)
  opt_pot <- zacetna
  opt_teza <- teza(A,zacetna)
  
  #vrednosti vnesene v tabelo
  rotacija <- c(" ")
  pot <- paste0(as.character(zacetna), collapse = "-")
  teza <- c(opt_teza)
  
  for (i in 2:(n-1)){ #koliko jih rotiramo, zacetek in konec je fiksen
    for (j in 2:(n-i)){# zacetek "drsece" rotacije, odvisne od i, tj. koliko jih v koraku obrnemo
      nova_pot <- c(opt_pot[1:(j-1)],rev(opt_pot[j:(j+i-1)]),opt_pot[(j+i):n])
      nova_teza <- teza(A,nova_pot)
      if (nova_teza < opt_teza && !is.na(nova_teza)){
      opt_pot <- nova_pot
      opt_teza <- nova_teza
      #Vnasamo v tabelo
      rotacija <- append(rotacija,paste0(as.character(rev(opt_pot[j:(j+i-1)])),collapse="-"))
      pot <- append(pot,paste0(as.character(opt_pot), collapse = "-"))
      teza <- append(teza,opt_teza)
      }
      else {
        break #da ne mece errorjev na koncu
      }
    } 
    
  }
  tabela <- data_frame(rotacija,pot,teza)
  return(tabela)
}

### EKSAKTNI ALGORITMI (BB algoritem)
permutacije <- function(A){ #vhod je matrika velikosti nx2 (rešitev nxn assignment problema), izhod pa vsi cikli
  n=nrow(A)
  perm <- list()
  Q=A[,1]
  i=1 #stevec da dodajamo v list
  while (!length(Q)==0){
    cikel <- c()
    zacetek <- Q[1]
    cikel <- append(cikel,zacetek)
    Q <- setdiff(Q,zacetek)
    indeks <- which(A[,1]==zacetek)
    while (!A[indeks,2]==zacetek){
      j <- A[indeks,2]
      cikel <- append(cikel,j)
      Q <- setdiff(Q,j)
      indeks <- which(A[,1]==j)
    }
    cikel <- as.integer(append(cikel,zacetek))
    perm[[i]] <- cikel
    i <- i+1
  }
  return(perm)
}

