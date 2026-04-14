#gradu koodi lepido
library(MASS)
library(cluster)
library(pheatmap)
library(mclust)
library(arrow)
library(fossil)
library(aricode)
library(dendextend)

#UUSITTU
#####ANALYYSI ALKAA#####

#k = 43, 40, 46, 37, 41 kun pienet mukana FAMILY
#k = 46, 51, 42, 47, 48, kun pienet mukana SUBFAMILY (analyysi eri tiedostossa)



#KMEDOIDS: PAM

#dna family
PAMajo0f80 <- pamLaskentaF(lepi0F80,lepi0Fmeta, k = 43)
PAMajo1f80 <- pamLaskentaF(lepi1F80,lepi1Fmeta, k = 40)
PAMajo2f80 <- pamLaskentaF(lepi2F80,lepi2Fmeta, k = 46)
PAMajo3f80 <- pamLaskentaF(lepi3F80,lepi3Fmeta, k = 37)
PAMajo4f80 <- pamLaskentaF(lepi4F80,lepi4Fmeta, k = 41)

#piirteet family
PAMajo0f <- pamLaskentaF(lepi0F,lepi0Fmeta, k = 43)
PAMajo1f <- pamLaskentaF(lepi1F,lepi1Fmeta, k = 40)
PAMajo2f <- pamLaskentaF(lepi2F,lepi2Fmeta, k = 46)
PAMajo3f <- pamLaskentaF(lepi3F,lepi3Fmeta, k = 37)
PAMajo4f <- pamLaskentaF(lepi4F,lepi4Fmeta, k = 41)

#automaatinen tulosten luku uusifu -tiedostossa

#indeksit family
indeksilaskuri(PAMajo0f80$a, PAMajo0f$a);print("######### B ALKAA #########"); indeksilaskuri(PAMajo0f80$b, PAMajo0f$b)
indeksilaskuri(PAMajo1f80$a, PAMajo1f$a);print("######### B ALKAA #########"); indeksilaskuri(PAMajo1f80$b, PAMajo1f$b)
indeksilaskuri(PAMajo2f80$a, PAMajo2f$a);print("######### B ALKAA #########"); indeksilaskuri(PAMajo2f80$b, PAMajo2f$b)
indeksilaskuri(PAMajo3f80$a, PAMajo3f$a);print("######### B ALKAA #########"); indeksilaskuri(PAMajo3f80$b, PAMajo3f$b)
indeksilaskuri(PAMajo4f80$a, PAMajo4f$a);print("######### B ALKAA #########"); indeksilaskuri(PAMajo4f80$b, PAMajo4f$b)

#####################################################################

#hierarkkinen jakava: DIANA

#dna family
Dia0F80 <- dianalaskentaF(lepi0F80,lepi0Fmeta, k = 43)
Dia1F80 <- dianalaskentaF(lepi1F80,lepi1Fmeta, k = 40)
Dia2F80 <- dianalaskentaF(lepi2F80,lepi2Fmeta, k = 46)
Dia3F80 <- dianalaskentaF(lepi3F80,lepi3Fmeta, k = 37)
Dia4F80 <- dianalaskentaF(lepi4F80,lepi4Fmeta, k = 41)

#piirteet family
Dia0F <- dianalaskentaF(lepi0F,lepi0Fmeta, k = 43)
Dia1F <- dianalaskentaF(lepi1F,lepi1Fmeta, k = 40)
Dia2F <- dianalaskentaF(lepi2F,lepi2Fmeta, k = 46)
Dia3F <- dianalaskentaF(lepi3F,lepi3Fmeta, k = 37)
Dia4F <- dianalaskentaF(lepi4F,lepi4Fmeta, k = 41)

#automaatinen tulosten luku uusifu -tiedostossa

#indeksit
indeksilaskuri(Dia0F80$a,Dia0F$a);print("######### B ALKAA #########");indeksilaskuri(Dia0F80$b,Dia0F$b)
indeksilaskuri(Dia1F80$a,Dia1F$a);print("######### B ALKAA #########");indeksilaskuri(Dia1F80$b,Dia1F$b)
indeksilaskuri(Dia2F80$a,Dia2F$a);print("######### B ALKAA #########");indeksilaskuri(Dia2F80$b,Dia2F$b)
indeksilaskuri(Dia3F80$a,Dia3F$a);print("######### B ALKAA #########");indeksilaskuri(Dia3F80$b,Dia3F$b)
indeksilaskuri(Dia4F80$a,Dia4F$a);print("######### B ALKAA #########");indeksilaskuri(Dia4F80$b,Dia4F$b)

#hierarkkinen yhdistävä: HCLUST
options(max.print=2500)
#DNA

D080 <- as.dist(lepi0F80)
D180 <- as.dist(lepi1F80)
D280 <- as.dist(lepi2F80)
D380 <- as.dist(lepi3F80)
D480 <- as.dist(lepi4F80)

ajo0_80 <- hclustlaskentaF(D080,lepi0Fmeta, k = 43 )
ajo0_80
ajo1_80 <- hclustlaskentaF(D180,lepi1Fmeta, k = 40 )
ajo2_80 <-hclustlaskentaF(D280,lepi2Fmeta, k = 46 )
ajo3_80 <-hclustlaskentaF(D380,lepi3Fmeta, k = 37 )
ajo4_80 <-hclustlaskentaF(D480,lepi4Fmeta, k = 41 )

#PIIRTEET
lepi0F <- as.dist(lepi0F)
lepi1F <- as.dist(lepi1F)
lepi2F <- as.dist(lepi2F)
lepi3F <- as.dist(lepi3F)
lepi4F <- as.dist(lepi4F)

ajo0f <- hclustlaskentaF(lepi0F,lepi0Fmeta, k = 43 )
ajo1f <- hclustlaskentaF(lepi1F,lepi1Fmeta, k = 40 )
ajo2f <- hclustlaskentaF(lepi2F,lepi2Fmeta, k = 46 )
ajo3f <- hclustlaskentaF(lepi3F,lepi3Fmeta, k = 37 )
ajo4f <- hclustlaskentaF(lepi4F,lepi4Fmeta, k = 41 )

#automaatinen tulosten luku uusifu -tiedostossa

#family
indeksilaskuri(ajo0_80$a, ajo0f$a);print("######### B ALKAA #########");indeksilaskuri(ajo0_80$b, ajo0f$b);print("######### C ALKAA #########");indeksilaskuri(ajo0_80$c, ajo0f$c)
indeksilaskuri(ajo1_80$a, ajo1f$a);print("######### B ALKAA #########");indeksilaskuri(ajo1_80$b, ajo1f$b);print("######### C ALKAA #########");indeksilaskuri(ajo1_80$c, ajo1f$c)
indeksilaskuri(ajo2_80$a, ajo2f$a);print("######### B ALKAA #########");indeksilaskuri(ajo2_80$b, ajo2f$b);print("######### C ALKAA #########");indeksilaskuri(ajo2_80$c, ajo2f$c)
indeksilaskuri(ajo3_80$a, ajo3f$a);print("######### B ALKAA #########");indeksilaskuri(ajo3_80$b, ajo3f$b);print("######### C ALKAA #########");indeksilaskuri(ajo3_80$c, ajo3f$c)
indeksilaskuri(ajo4_80$a, ajo4f$a);print("######### B ALKAA #########");indeksilaskuri(ajo4_80$b, ajo4f$b);print("######### C ALKAA #########");indeksilaskuri(ajo4_80$c, ajo4f$c)

#####