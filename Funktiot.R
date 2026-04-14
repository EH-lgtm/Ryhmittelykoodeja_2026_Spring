#gradukoodi, Essi Heiskanen, Funktiot

#####PAKETIT#####
library(DescTools)
library(MASS)
library(cluster)
library(pheatmap)
library(mclust)
library(arrow)
library(fossil)
library(aricode)
library(dendextend)
library(ggplot2)
library(reshape2)
library(tidyverse)

set.seed(11)

#options(max.print=2500) Tätä tullaan tarvitsemaan

###############################

#FUNKTIOT:

###############################

##indeksien laskenta

indeksilaskuri <- function(DNAclust, PiirreClust){
  #[0:1] #rand indeksi
  R <- rand.index(DNAclust, PiirreClust)
  #[-1:1] #adjusted rand ind
  R_A <- adjustedRandIndex(DNAclust, PiirreClust)
  #[0:1] #normalized mutual information
  #N_M_I <- NMI(DNAclust, PiirreClust) 
  #[0:1] #fowlkes-mallows index FMI (hierarkkisille, ei tietoa toimiiko muille)
  F_M <- FM_index_R(DNAclust, PiirreClust)
  A_M_I <-AMI(DNAclust, PiirreClust)
  #print("Rand");print(R)
  #print("AdjRand");print(R_A)
  #print("NMI");print(N_M_I)
  #print("FM");print(F_M)
  return(list(R =R,RA = R_A, a_m_i=A_M_I, fm=F_M))
  
}

indeksiyhd_0 <- function(L0,L1,L2,L3,L4){
  R <- list(L0$R,L1$R,L2$R,L3$R,L4$R)
  RA <- list(L0$RA,L1$RA,L2$RA,L3$RA,L4$RA)
  ami <- list(L0$a_m_i,L1$a_m_i,L2$a_m_i,L3$a_m_i,L4$a_m_i)
  fm <- list(L0$fm,L1$fm,L2$fm,L3$fm,L4$fm)
  
  R<- list(round(min(unlist(R)),2),round(mean(unlist(R)),2),round(max(unlist(R)),2))
  RA<-list(round(min(unlist(RA)),2),round(mean(unlist(RA)),2),round(max(unlist(RA)),2))
  a_m_i<-list(round(min(unlist(ami)),2),round(mean(unlist(ami)),2),round(max(unlist(ami)),2))
  fm<-list(round(min(unlist(fm)),2),round(mean(unlist(fm)),2),round(max(unlist(fm)),2))
  return(list(R =R, RA=RA,ami=a_m_i,fm=fm))
}


########### 
# HCLUST

#FAMILY FUNKTIO: HCLUST
hclustlaskentaF <- function(sampleDNAtaiPiirre, sampleMeta, k){
  dendro0a <- hclust(d= sampleDNAtaiPiirre, method = 'ward.D') 
  dendro0b <- hclust(d= sampleDNAtaiPiirre, method = 'ward.D2')
  dendro0c <- hclust(d= sampleDNAtaiPiirre, method = 'complete')
  
  plot(dendro0a, labels = F)
  rect.hclust(dendro0a, k = k)
  plot(dendro0b, labels = F)
  rect.hclust(dendro0b, k = k)
  plot(dendro0c, labels = F)
  rect.hclust(dendro0c, k = k)
  
  print(table(cutree(dendro0a,k = k)))
  print(table(cutree(dendro0a,k = k), sampleMeta$family ))
  apu1Fam <- table(cutree(dendro0a,k = k), sampleMeta$family )
  print(table(sampleMeta$family))
  print("valmis a, b alkaa")
  print(table(cutree(dendro0b,k = k)))
  print(table(cutree(dendro0b,k = k), sampleMeta$family ))
  apu2Fam <- table(cutree(dendro0b,k = k), sampleMeta$family )
  print(table(sampleMeta$family))
  print("valmis b, c alkaa")
  print(table(cutree(dendro0c,k = k)))
  print(table(cutree(dendro0c,k = k), sampleMeta$family ))
  apu3Fam <- table(cutree(dendro0c,k = k), sampleMeta$family )
  print(table(sampleMeta$family))
  print("valmis c")
  
  A <- cutree(dendro0a,k = k)
  B <- cutree(dendro0b,k = k)
  C <- cutree(dendro0c,k = k)
  
  return(list(a = A,b = B,c = C, D = apu1Fam, E = apu2Fam, f = apu3Fam))
  
}

#ORDER FUNKTIO: HCLUST
hclustlaskentaOR <- function(sampleDNAtaiPiirre, sampleMeta, k){
  dendro0a <- hclust(d= sampleDNAtaiPiirre, method = 'ward.D') 
  dendro0b <- hclust(d= sampleDNAtaiPiirre, method = 'ward.D2')
  dendro0c <- hclust(d= sampleDNAtaiPiirre, method = 'complete')
  
  plot(dendro0a, labels = F)
  rect.hclust(dendro0a, k = k)
  plot(dendro0b, labels = F)
  rect.hclust(dendro0b, k = k)
  plot(dendro0c, labels = F)
  rect.hclust(dendro0c, k = k)
  
  print(table(cutree(dendro0a,k = k)))
  print(table(cutree(dendro0a,k = k), sampleMeta$order ))
  apu1Ord <- table(cutree(dendro0a,k = k), sampleMeta$order )
  print(table(sampleMeta$order))
  print("valmis a, b alkaa")
  print(table(cutree(dendro0b,k = k)))
  print(table(cutree(dendro0b,k = k), sampleMeta$order ))
  apu2Ord <- table(cutree(dendro0b,k = k), sampleMeta$order )
  print(table(sampleMeta$order))
  print("valmis b, c alkaa")
  print(table(cutree(dendro0c,k = k)))
  print(table(cutree(dendro0c,k = k), sampleMeta$order ))
  apu3Ord <- table(cutree(dendro0c,k = k), sampleMeta$order )
  print(table(sampleMeta$order))
  print("valmis c")
  
  A <- cutree(dendro0a,k = k)
  B <- cutree(dendro0b,k = k)
  C <- cutree(dendro0c,k = k)
  
  return(list(a = A,b = B,c = C, D = apu1Ord, E = apu2Ord, f = apu3Ord))
}


##########

#KMEDOIDS: PAM

pamLaskentaF <- function(sampleDNAtaiPiirre, samplemeta, k){
  pam_resA <- pam(sampleDNAtaiPiirre, diss = T, k = k, metric = "euclidean") #ota diss pois jos ei syötetä dissimilaritya (tällöin myös euc ja man tekevät erilaiset ryhmittelyt)
  print("keskukset A (10kpl)"); print(pam_resA$medoids[1:10])
  pam_resB <- pam(sampleDNAtaiPiirre, diss =T, k = k, metric = "manhattan")
  print("keskukset B (10kpl)"); print(pam_resB$medoids[1:10])
  
  print(table(pam_resA$clustering))
  print(table(pam_resA$clustering, samplemeta$family))
  apu1<- table(pam_resA$clustering, samplemeta$family)
  
  print(table(samplemeta$family))
  print(table(pam_resB$clustering))
  print(table(pam_resB$clustering, samplemeta$family))
  apu2<- table(pam_resB$clustering, samplemeta$family)
  
  print(table(samplemeta$family))
  
  A <- pam_resA$clustering
  B <- pam_resB$clustering
  return(list(a = A,b = B, pam1F = apu1, pam2F = apu2))
}


pamLaskentaOR <- function(sampleDNAtaiPiirre, samplemeta, k){
  pam_resA <- pam(sampleDNAtaiPiirre,diss = T, k = k, metric = "euclidean")
  print("keskukset A (10kpl)"); print(pam_resA$medoids[1:10])
  pam_resB <- pam(sampleDNAtaiPiirre,diss = T, k = k, metric = "manhattan")
  print("keskukset B (10kpl)"); print(pam_resB$medoids[1:10])
  
  print(table(pam_resA$clustering))
  print(table(pam_resA$clustering, samplemeta$order))
  apu1<- table(pam_resA$clustering, samplemeta$order)
  
  print(table(samplemeta$order))
  print(table(pam_resB$clustering))
  print(table(pam_resB$clustering, samplemeta$order))
  apu2 <- table(pam_resB$clustering, samplemeta$order)
  print(table(samplemeta$order))
  
  A <- pam_resA$clustering
  B <- pam_resB$clusterin
  return(list(a = A,b = B, pam1Or = apu1, pam2Or = apu2))
}

##############

#hierarkkinen jakava: DIANA

#family
dianalaskentaF <- function(sampleDNAtaiPiirre, samplemeta, k){
  diana_resA <- diana(sampleDNAtaiPiirre, diss = T, metric ="euclidean") #metric tekee asioita vain jos diss=F, mutta nyt käsitellään vain diss otuksia
  print("A dc:");print(diana_resA$dc)
  diana_resB <- diana(sampleDNAtaiPiirre, diss = T, metric ="manhattan")
  print("B dc:");print(diana_resB$dc)
  
  A <- cutree(as.hclust(diana_resA),k = k)
  plot(as.hclust(diana_resA), labels= F, main = "Diana dendro")
  print("A"); print(table(A))
  print(table(A, samplemeta$family ))
  print(table(samplemeta$family))
  rect.hclust(as.hclust(diana_resA), k=k) 
  apu1 <- table(A, samplemeta$family)
  
  B <- cutree(as.hclust(diana_resB),k = k)
  plot(as.hclust(diana_resA), labels= F, main = "Diana dendro")
  print("B");print(table(B))
  print(table(B, samplemeta$family ))
  print(table(samplemeta$family))
  rect.hclust(as.hclust(diana_resB), k=k)
  apu2 <- table(B, samplemeta$family)
  
  return(list(a = A, b = B, diana1F =apu1,diana2F = apu2))
}


#order
dianalaskentaOR <- function(sampleDNAtaiPiirre, samplemeta, k){
  diana_resA <- diana(sampleDNAtaiPiirre, diss = T, metric ="euclidean")
  print("A dc:");print(diana_resA$dc)
  diana_resB <- diana(sampleDNAtaiPiirre, diss = T, metric ="manhattan")
  print("B dc:");print(diana_resB$dc)
  
  A <- cutree(as.hclust(diana_resA),k = k)
  plot(as.hclust(diana_resA), labels= F, main = "Diana dendro")
  print("A"); print(table(A))
  print(table(A, samplemeta$order ))
  print(table(samplemeta$order))
  rect.hclust(as.hclust(diana_resA), k=k) 
  apu1 <- table(A, samplemeta$order)
  
  B <- cutree(as.hclust(diana_resB),k = k)
  plot(as.hclust(diana_resA), labels= F, main = "Diana dendro")
  print("B");print(table(B))
  print(table(B, samplemeta$order ))
  print(table(samplemeta$order))
  rect.hclust(as.hclust(diana_resB), k=k)
  apu2 <- table(B, samplemeta$order)
  
  return(list(a = A, b = B, diana1Or =apu1,diana2Or = apu2))
}


#####################
#Tuloksien luvun automatisointi

#automata = painottamaton kaava
Automata <- function(df, threshold = 0) {
  
  df <- as.matrix(df)
  rs <- rowSums(df)
  
  values <- sapply(seq_len(ncol(df)), function(j) {
    a <- df[, j]
    idx <- a > 0 & rs > 0
    
    if (!any(idx)) return(0)
    
    incorrect <- sum(rs[idx] - a[idx])
    correct   <- sum(a[idx])
    
    (incorrect / correct)
  })
  
  names(values) <- colnames(df)
  
  # valitaan col, jotka alle rajan
  hits <- values[values < threshold]
  if(length(hits)==0){print("Ei sarakkeita alle rajan (painottamaton)")}
  else {print("Sarakkeet alle rajan (painottamaton):"); print(hits);
    invisible(hits)}
}

##automata3 = painotettu kaava

Automata3 <- function(df, threshold = 0) {
  
  df <- as.matrix(df)
  rs <- rowSums(df)
  
  values <- numeric(ncol(df))
  names(values) <- colnames(df)
  
  for (j in seq_len(ncol(df))) {
    a <- df[, j]
    idx <- a > 0 & rs > 0
    
    if (!any(idx)) {
      values[j] <- 0
    } else {
      w <- a[idx]/sum(a[idx])
      values[j] <- sum((w*(rs[idx]-a[idx]))/sum(a[idx]))
      
    }
  }
  
  # valitaan col, jotka alle rajan
  hits <- values[values < threshold]
  if(length(hits)==0){print("Ei sarakkeita alle rajan (painotettu)")}
  else {print("Sarakkeet alle rajan (painotettu):");print(hits);
    invisible(hits)}
}


Yhdistäjä <- function(DNAtulos,PiirreTulos){
  nimet1 <- names(DNAtulos)
  nimet2 <- names(PiirreTulos)
  #if(is.null(intersect(nimet1, nimet2))){print("Ei yhteisiä :<")}
  if(length(intersect(nimet1, nimet2))== 0){print("Ei yhteisiä")}
  else {print("Löytyy molemmista:");intersect(nimet1, nimet2)}
  
}


Laskija <- function(H1, H2, H3,P,D){
  dt <- c(H1, H2, H3,P,D)
  A <- sort(table(dt),decreasing = T)
  print("Tulos:");print(A)
}

LaskijaSuper <- function(lista) {
  
  # Get all unique names
  nimet <- unique(unlist(lapply(lista, names)))
  
  # kasvatetaan jokainen table sisältämään kaikki nimet
  aligned <- lapply(lista, function(tab) {
    out <- setNames(rep(0, length(nimet)), nimet)
    out[names(tab)] <- tab
    out
  })
  
  # yhdistys
  A <- Reduce("+", aligned)
  A <- sort(A,decreasing = T)
  print("Tulos:");print(A)
}
