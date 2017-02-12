require(dplyr)
require(lpSolve) #resevanje Madzarske metode
require(xtable)
require(ggplot2)
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
  teza <- teza(A,pot)
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
permutacije <- function(A){ #vhod je matrika, izhod pa permutacije (rešitev assingment problem)
  # pa list ciklov
  A[!is.finite(A)] <- 1000 #nastavimo dovolj visoko, saj lp.assign ne zna delati z Inf, ne pa preveč,
  #ker pride do težav iz neznanih razlogov
  C <-lp.assign(A,direction = "min")$solution
  resitev <- which(!C == 0,arr.ind = TRUE)
  n=nrow(resitev)
  perm <- list()
  Q=resitev[,1]
  i=1 #stevec da dodajamo v list
  while (!length(Q)==0){
    cikel <- c()
    zacetek <- Q[1]
    cikel <- append(cikel,zacetek) 
    Q <- setdiff(Q,zacetek) #vozlisce odstranimo iz Q
    indeks <- which(resitev[,1]==zacetek)
    while (!resitev[indeks,2]==zacetek){
      j <- resitev[indeks,2]
      cikel <- append(cikel,j)
      Q <- setdiff(Q,j)
      indeks <- which(resitev[,1]==j)
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

teza_permutacij <- function(A){ #vhod je matrika
  resitev <- permutacije(A)
  n <- length(resitev)
  cena <- 0
  for (i in 1:n){
    cikel <- resitev[[i]]
    cena <- cena + teza(A,cikel)
  }
  return(cena)
}

veja <- function(A){ #resevanje enega nivoja za eno vejo BB algoritma, vrne bodisi resitev, bodisi podproblem
  resitev <- permutacije(A)
  lista_matrik <- list()
  if (length(resitev)==1){ #ce imamo en sam cikel, smo dobili obhod, vrnemo pot in ceno poti
    cena <- teza_permutacij(A)
    pot <- paste0(as.character(resitev[[1]]),collapse = "-")
    lista_matrik[[1]]<-c(cena,pot)
  }
  else{ #sicer razbivamo najmanjsi cikel, dobimo vec manjsih podproblemov
  cikel <- resitev[[indeks_cikla(resitev)]] #določitev cikla z najmnajšo dolžino
  for (i in (1:(length(cikel)-1))){
    B <- A
    B[cikel[i],cikel[i+1]] <- Inf
    lista_matrik[[i]] <- B
  }
  }
  return(lista_matrik)
}

nivo <- function(A){#vhod je nek list, v katerem imamo class character (že rešeno) in matrix.
  n <- length(A) #dolzina nivoja
  lista <- list()
  zgornja <- 10000
  for (i in (1:n)){
    if (is.character(A[[i]])){
      lista[[i]] <- A[[i]]
    }
    else{
      resitev <- veja(A[[i]])
      a <- length(resitev)
      for (j in 1:a){
        if (is.character(resitev[[j]])){
          lista[[i+j-1]] <- resitev[[j]]
        }
        else{
          if (zgornja > teza_permutacij(resitev[[j]])){
          lista[[i+j-1]] <- resitev[[j]]
          zgornja <- teza_permutacij(resitev[[j]])
          }
          else {
            lista[[i+j-1]] <- c("")
          }
        }
        n <- n + a - 1
      }
    }
  }
  return(lista)
}

#BB ALGORITEM - dodamo upper_limit + ostale omejitve,
#ki pospesijo algoritem, ni nam treba racunati celotnega drevesa

bb <- function(A){
  k <- veja(A)
  if (length(k)==1){ #že v štartu dobimo cikel in ni treba delat vej naprej
    return(c(k[[1]][1],k[[1]][2]))
  }
  logicne <- c(FALSE)
  while (!all(logicne)){
    k<-nivo(k)
    n <- length(k)
    logicne <- c()

  for (i in 1:n){
    logicne[i] <- is.character(k[[i]][1])
  }

  }
  pot <- c()
  cena <- c()
  for (j in 1:n){
  pot <- append(pot,k[[j]][2])
  cena <- append(cena,as.numeric(k[[j]][1]))
  }
  indeks <- which.min(cena)
 return(c(pot[indeks],cena[indeks]))
}
  
###PRIMERJAVE METOD


primerjava <- function(m){#primerjajmo kako dobra je metoda najblizjih
  #sosedov, naredimo m simulacij in gledamo relativno napako, zacnemo v random tocki
  relativna1 <- c()
  relativna2 <- c()
  relativna3 <- c()
  for (i in 1:m){
    
    A <- simetricna(10)
    nakljucni_zacetek <- sample(1:10,1)
    prava_vrednost <- as.numeric(bb(A)[2])
    pot_sosed1 <- najblizji_sosed(A,nakljucni_zacetek)
    pot_teza1 <- teza(A,pot_sosed1)
    napaka1 <- abs(pot_teza1-prava_vrednost)/prava_vrednost
    relativna1 <- append(relativna1,napaka1)

    zacetek <- opt_zacetek(A)
    pot_sosed2 <- najblizji_sosed(A,zacetek)
    pot_teza2 <- teza(A,pot_sosed2)
    napaka2 <- abs(pot_teza2-prava_vrednost)/prava_vrednost
    relativna2 <- append(relativna2,napaka2)

    pot_teza3 <- as.numeric(reversal(A)$teza)
    napaka3 <- abs(pot_teza3-prava_vrednost)/prava_vrednost
    relativna3 <- append(relativna3,napaka3)

  }
  povp1 <- mean(relativna1)
  histogram1 <- hist(relativna1,probability = TRUE,main = paste0("Relativna napaka metode najbližjih sosedov"),
                     xlab = "Relativna napaka",ylab = "Gostota",breaks=30)
  histogram1 <- abline(v=povp1,col="red",lwd=2)
  histogram1 <- legend("topleft",col=c("red"), lwd=2,bty="n",xjust=1, seg.len=0.5,
                       legend=c(paste0("Povprečje = ",round(povp1,4))), cex=0.8)
  
  povp2 <- mean(relativna2)
  histogram2 <- hist(relativna2,probability = TRUE,main = paste0("Relativna napaka opt. metode najbližjih sosedov"),
                     xlab = "Relativna napaka",ylab = "Gostota",breaks = 30)
  histogram2 <- abline(v=povp2,col="red",lwd=2)
  histogram2 <- legend("topleft",col=c("red"), lwd=2,bty="n",xjust=1, seg.len=0.5,
                       legend=c(paste0("Povprečje = ",round(povp2,4))), cex=0.8)
  
  povp3 <- mean(relativna3)
  histogram3 <- hist(relativna3,probability = TRUE,main = paste0("Relativna napaka metode obratov"),
                     xlab = "Relativna napaka",ylab = "Gostota",breaks = 30)
  histogram3 <- abline(v=povp3,col="red",lwd=2)
  histogram3 <- legend("topleft",col=c("red"), lwd=2,bty="n",xjust=1, seg.len=0.5,
                       legend=c(paste0("Povprečje = ",round(povp3,4))), cex=0.8)
  
  par(mfrow=c(3,1))
  histogram1
  histogram2
  histogram3
  return()
}

### RAČUNANJE ČASA ZA IZVEDBO PROGRAMA PRI RAZLIČNIH VELIKOSTIH MATRIKE

casovna_zahtevnost <- function(n){ #gledali bomo, koliko casa porabimo za izvedbo
  #TSP z razlicnimi algoritmi in kako se cas spreminja z vecanjem velikosti matrike
  # pri vsaki velikosti bomo z vsako metodo (sosedi,reversal,BB) po 10-krat izracunali
  #cas izvajanja ter vzeli povprecnega. Velikost matrike bo tekla od 1x1 do 20x20
  #gradili bomo tabelo, za kasnjejse delo z dplyr in lazjo graficno predstavitev
  #Tidy data (velikost matrike, tip resevanja, povp cas)
  velikost <- c()
  tip <- c() #nacin resevanja ("sosed","rev","bb")
  cas_tabela <- c()
  
  for (s in 2:n){
    cas_sosed_sim <- c()
    cas_reversal_sim <- c()
    cas_bb_sim <- c()
    for (j in 1:100){
    A <- simetricna(s) #simulacija matrike
    z <- opt_zacetek(A)
      #NAJBLIZJI SOSED
      start1 <- proc.time()
      k <- najblizji_sosed(A,z)
      konec1 <- proc.time()
      cas1 <- (konec1-start1)[1]
      cas_sosed_sim <- append(cas_sosed_sim,cas1)
      #REVERSAL
      start2 <- proc.time()
      reversal(A,k)
      konec2 <- proc.time()
      cas2 <- (konec2-start2)[1]
      cas_reversal_sim <- append(cas_reversal_sim,cas2+cas1) #kaj je optimalna pot,
      #smo ze poracunali, vendar je to vhod pri reversal metodi, zato pristejemo tudi cas1
      #BB ALGORITEM
      start3 <- proc.time()
      bb(A)
      konec3 <- proc.time()
      cas3 <- (konec3-start3)[1]
      cas_bb_sim <- append(cas_bb_sim,cas3)
    }
    velikost <- append(velikost,rep(s,3))
    tip <- append(tip,c("sosed","rev","bb"))
    cas_tabela <- append(cas_tabela,c(mean(cas_sosed_sim),mean(cas_reversal_sim),mean(cas_bb_sim)))
  }
  return(data_frame(velikost,tip,cas_tabela))
}

#graf <- ggplot(data = tabela2,aes(x =velikost,y=cas_tabela,colour=tip)) + geom_point() + geom_line() + xlab("Velikost matrike")+
#  ylab("Časovna zahtevnost") + ggtitle("Časovna zahtevnost v odvisnosti od velikosti matrike")
