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
permutacije <- function(A){ #vhod je matrika velikosti nx2 (rešitev nxn assignment problema), izhod
  # pa list ciklov
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

indeks_cikla <- function(l){ #indeks minimalnega podcikla v list
  n <- length(l)
  dolzine <- c()
  for (i in 1:n){
    dolzine <- append(dolzine,length(l[[i]]))
  }
  return(which.min(dolzine))
  
  
}

teza_permutacij <- function(A,resitev){ #vhod je lista permutacij
  n <- length(resitev)
  cena <- 0
  for (i in 1:n){
    cikel <- resitev[[i]]
    cena <- cena + teza(A,cikel)
  }
  return(cena)
}

veja <- function(A){ #resevanje enega nivoja za eno vejo BB algoritma, vrne bodisi resitev, bodisi podproblem
  A[!is.finite(A)]<-100 #nastavimo inf na velike vrednosti, lpSolve ne zna delati z inf
  #IN NITI Z "PRE"VELIKIMI VREDNOSTMI!!!
  C <-lp.assign(A,direction = "min")$solution
  resitev <- permutacije(which(!C == 0,arr.ind = TRUE))
  upper <- reversal(A)$teza
  if (length(resitev)==1){ #ce imamo en sam cikel, smo dobili obhod, vrnemo pot in ceno poti
    cena <- teza_permutacij(A,resitev)
    pot <- paste0(as.character(resitev[[1]]),collapse = "-")
    if (cena <= upper){
      upper <- cena #NIZANJE UPPER LIMIT
    return(c(pot,cena))
    }
    else{
      return(c(""))
    }
  }
  else{ #sicer razbivamo najmanjsi cikel, dobimo vec manjsih podproblemov
  cikel <- resitev[[indeks_cikla(resitev)]]
  lista_matrik <- list()
  for (i in (1:(length(cikel)-1))){
    B <- A
    B[cikel[i],cikel[i+1]] <- 100
    B1 <- lp.assign(B)$solution
    resitev1 <- cena <- teza_permutacij(A,permutacije(which(!B1 == 0,arr.ind = TRUE)))
    if (resitev1 < upper){
    lista_matrik[[i]]<-B  
    }
    else{
      lista_matrik[[i]] <- c("") #ZAPIRANJE VEJE, CE JE LE TA ZE NAD RESITVIJO IZ ALGORITMA REVERSAL
    }
    
  
  }
  return(lista_matrik)
  }
}

nivo <- function(A){#vhod je nek list, v katerem imamo class character in matrix.
  n <- length(A) #dolzina nivoja
  lista <- list()
  for (i in (1:n)){
    if (is.character(A[[i]])){
      lista[[i]] <- A[[i]] #kar smo že rešili samo prepišemo
    }
    else{
      resitev <- veja(A[[i]])
      if (is.character(resitev)){
        lista[[i]] <- resitev
      }
      else{ # dodajanje matrike v resitev
        a=length(resitev) #na koliko vej se cepi
        for (k in (0:(a-1))){
          lista[[i+k]]<-resitev[[k+1]]
        }
        n <- n + (a-1)
      }
    }
  }
  return(lista)
}

#BB ALGORITEM - dodamo upper_limit + ostale omejitve,
#ki pospesijo algoritem, ni nam treba racunati celotnega drevesa

bb <- function(A){
  k <- veja(A)
  logicne <- c(FALSE)
  while (!all(logicne)){
    k<-nivo(k)
    n <- length(k)
    logicne <- c()
  for (i in 1:n){
    logicne[i] <- is.character(k[[i]])
  }

  }
  pot <- c()
  cena <- c()
  min <- reversal(A)$teza
  for (j in 1:n){
    if ((!is.na(as.integer(k[[j]][2]))) && (as.integer(k[[j]][2])<=min)){
      cena <- as.integer(k[[j]][2])
      min <- cena
      pot <- k[[j]][1]
    }
    else{
      next
    }
    
  }
 return(c(pot,cena))
}
  
###PRIMERJAVE METOD


primerjava_najblizji <- function(m){#primerjajmo kako dobra je metoda najbliznjih
  #sosedov, naredimo m simulacij in gledamo relativno napako
  relativna <- c()
  for (i in 1:m){
    A <- simetricna(10)
    nakljucni_zacetek <- sample(1:10,1)
    prava_vrednost <- as.numeric(bb(A)[2])
    pot_sosed <- najblizji_sosed(A,nakljucni_zacetek)
    pot_teza <- teza(A,pot_sosed)
    napaka <- abs(pot_teza-prava_vrednost)/prava_vrednost
    relativna <- append(relativna,napaka)
  }
  povp <- mean(relativna)
  histogram <- hist(relativna,probability = TRUE,main = paste0("Relativna napaka metode najbližjih sosedov"),
                    xlab = "Relativna napaka",ylab = "Gostota",breaks=30)
  histogram <- abline(v=mean(relativna),col="red",lwd=2)
  histogram <- legend("topright",col=c("red"), lwd=2,
                      legend=c(paste0("Povprečje = ",round(povp,2))), cex=0.8)
  
  return(histogram)
}

primerjava_optnajblizji <- function(m){#primerjajmo kako dobra je metoda najbliznjih
  #sosedov, ce zacnemo v vozliscu, ki je optimalno. naredimo m simulacij in gledamo relativno napako
  relativna <- c()
  for (i in 1:m){
    A <- simetricna(10)
    zacetek <- opt_zacetek(A)
    prava_vrednost <- as.numeric(bb(A)[2])
    pot_sosed <- najblizji_sosed(A,zacetek)
    pot_teza <- teza(A,pot_sosed)
    napaka <- abs(pot_teza-prava_vrednost)/prava_vrednost
    relativna <- append(relativna,napaka)
  }
  povp <- mean(relativna)
  histogram <- hist(relativna,probability = TRUE,main = paste0("Relativna napaka opt. metode najbližjih sosedov"),
                    xlab = "Relativna napaka",ylab = "Gostota",breaks = 30)
  histogram <- abline(v=mean(relativna),col="red",lwd=2)
  histogram <- legend("topright",col=c("red"), lwd=2,
                      legend=c(paste0("Povprečje = ",round(povp,2))), cex=0.8)
  
  return(histogram)
}

primerjava_reversal <- function(m){#primerjajmo kako dobra je metoda reversal
  #naredimo m simulacij in gledamo relativno napako
  relativna <- c()
  for (i in 1:m){
    A <- simetricna(10)
    prava_vrednost <- as.numeric(bb(A)[2])
    pot_teza <- as.numeric(reversal(A)$teza)
    napaka <- abs(pot_teza-prava_vrednost)/prava_vrednost
    relativna <- append(relativna,napaka)
  }
  povp <- mean(relativna)
  histogram <- hist(relativna,probability = TRUE,main = paste0("Relativna napaka metode obratov"),
                    xlab = "Relativna napaka",ylab = "Gostota",breaks = 30)
  histogram <- abline(v=mean(relativna),col="red",lwd=2)
  histogram <- legend("topright",col=c("red"), lwd=2,
                      legend=c(paste0("Povprečje = ",round(povp,2))), cex=0.8)
  
  return(histogram)
}