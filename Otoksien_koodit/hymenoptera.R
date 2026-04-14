#gradu koodi hymenoptera
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

#####ANALYYSI ALKAA#####

#k = 33, 37, 38, 38, 34 FAMILY


#KMEDOIDS: PAM

#dna family
PAMajo0f80 <- pamLaskentaF(hyme0F80,hyme0Fmeta, k = 33)
PAMajo1f80 <- pamLaskentaF(hyme1F80,hyme1Fmeta, k = 37)
PAMajo2f80 <- pamLaskentaF(hyme2F80,hyme2Fmeta, k = 38)
PAMajo3f80 <- pamLaskentaF(hyme3F80,hyme3Fmeta, k = 38)
PAMajo4f80 <- pamLaskentaF(hyme4F80,hyme4Fmeta, k = 34)

#piirteet family
PAMajo0f <- pamLaskentaF(hyme0F,hyme0Fmeta, k = 33)
PAMajo1f <- pamLaskentaF(hyme1F,hyme1Fmeta, k = 37)
PAMajo2f <- pamLaskentaF(hyme2F,hyme2Fmeta, k = 38)
PAMajo3f <- pamLaskentaF(hyme3F,hyme3Fmeta, k = 38)
PAMajo4f <- pamLaskentaF(hyme4F,hyme4Fmeta, k = 34)

#tulosten automaattinen luku uusifu:ssa

#indeksit family
indeksilaskuri(PAMajo0f80$a, PAMajo0f$a);print("######### B ALKAA #########"); indeksilaskuri(PAMajo0f80$b, PAMajo0f$b)
indeksilaskuri(PAMajo1f80$a, PAMajo1f$a);print("######### B ALKAA #########"); indeksilaskuri(PAMajo1f80$b, PAMajo1f$b)
indeksilaskuri(PAMajo2f80$a, PAMajo2f$a);print("######### B ALKAA #########"); indeksilaskuri(PAMajo2f80$b, PAMajo2f$b)
indeksilaskuri(PAMajo3f80$a, PAMajo3f$a);print("######### B ALKAA #########"); indeksilaskuri(PAMajo3f80$b, PAMajo3f$b)
indeksilaskuri(PAMajo4f80$a, PAMajo4f$a);print("######### B ALKAA #########"); indeksilaskuri(PAMajo4f80$b, PAMajo4f$b)

#hierarkkinen jakava: DIANA

#dna family
Dia0F80 <- dianalaskentaF(hyme0F80,hyme0Fmeta, k = 33)
Dia1F80 <- dianalaskentaF(hyme1F80,hyme1Fmeta, k = 37)
Dia2F80 <- dianalaskentaF(hyme2F80,hyme2Fmeta, k = 38)
Dia3F80 <- dianalaskentaF(hyme3F80,hyme3Fmeta, k = 38)
Dia4F80 <- dianalaskentaF(hyme4F80,hyme4Fmeta, k = 34)

#piirteet family
Dia0F <- dianalaskentaF(hyme0F,hyme0Fmeta, k = 33)
Dia1F <- dianalaskentaF(hyme1F,hyme1Fmeta, k = 37)
Dia2F <- dianalaskentaF(hyme2F,hyme2Fmeta, k = 38)
Dia3F <- dianalaskentaF(hyme3F,hyme3Fmeta, k = 38)
Dia4F <- dianalaskentaF(hyme4F,hyme4Fmeta, k = 34)

#tulosten automaattinen luku uusifu:ssa

#indeksit
indeksilaskuri(Dia0F80$a,Dia0F$a);print("######### B ALKAA #########");indeksilaskuri(Dia0F80$b,Dia0F$b)
indeksilaskuri(Dia1F80$a,Dia1F$a);print("######### B ALKAA #########");indeksilaskuri(Dia1F80$b,Dia1F$b)
indeksilaskuri(Dia2F80$a,Dia2F$a);print("######### B ALKAA #########");indeksilaskuri(Dia2F80$b,Dia2F$b)
indeksilaskuri(Dia3F80$a,Dia3F$a);print("######### B ALKAA #########");indeksilaskuri(Dia3F80$b,Dia3F$b)
indeksilaskuri(Dia4F80$a,Dia4F$a);print("######### B ALKAA #########");indeksilaskuri(Dia4F80$b,Dia4F$b)

#hierarkkinen yhdistävä: HCLUST
options(max.print=2500)
#DNA

D080 <- as.dist(hyme0F80)
D180 <- as.dist(hyme1F80)
D280 <- as.dist(hyme2F80)
D380 <- as.dist(hyme3F80)
D480 <- as.dist(hyme4F80)

ajo0_80 <- hclustlaskentaF(D080,hyme0Fmeta, k = 33 )
#ajo0_80$D
ajo1_80 <- hclustlaskentaF(D180,hyme1Fmeta, k = 37 )
ajo2_80 <-hclustlaskentaF(D280,hyme2Fmeta, k = 38 )
ajo3_80 <-hclustlaskentaF(D380,hyme3Fmeta, k = 38 )
ajo4_80 <-hclustlaskentaF(D480,hyme4Fmeta, k = 34 )

#PIIRTEET
hyme0F <- as.dist(hyme0F)
hyme1F <- as.dist(hyme1F)
hyme2F <- as.dist(hyme2F)
hyme3F <- as.dist(hyme3F)
hyme4F <- as.dist(hyme4F)

ajo0f <- hclustlaskentaF(hyme0F,hyme0Fmeta, k = 33 )
ajo1f <- hclustlaskentaF(hyme1F,hyme1Fmeta, k = 37 )
ajo2f <- hclustlaskentaF(hyme2F,hyme2Fmeta, k = 38 )
ajo3f <- hclustlaskentaF(hyme3F,hyme3Fmeta, k = 38 )
ajo4f <- hclustlaskentaF(hyme4F,hyme4Fmeta, k = 34 )

#tulosten automaattinen luku uusifu:ssa

#family
indeksilaskuri(ajo0_80$a, ajo0f$a);print("######### B ALKAA #########");indeksilaskuri(ajo0_80$b, ajo0f$b);print("######### C ALKAA #########");indeksilaskuri(ajo0_80$c, ajo0f$c)
indeksilaskuri(ajo1_80$a, ajo1f$a);print("######### B ALKAA #########");indeksilaskuri(ajo1_80$b, ajo1f$b);print("######### C ALKAA #########");indeksilaskuri(ajo1_80$c, ajo1f$c)
indeksilaskuri(ajo2_80$a, ajo2f$a);print("######### B ALKAA #########");indeksilaskuri(ajo2_80$b, ajo2f$b);print("######### C ALKAA #########");indeksilaskuri(ajo2_80$c, ajo2f$c)
indeksilaskuri(ajo3_80$a, ajo3f$a);print("######### B ALKAA #########");indeksilaskuri(ajo3_80$b, ajo3f$b);print("######### C ALKAA #########");indeksilaskuri(ajo3_80$c, ajo3f$c)
indeksilaskuri(ajo4_80$a, ajo4f$a);print("######### B ALKAA #########");indeksilaskuri(ajo4_80$b, ajo4f$b);print("######### C ALKAA #########");indeksilaskuri(ajo4_80$c, ajo4f$c)

