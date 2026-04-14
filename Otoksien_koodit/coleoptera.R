#gradu koodi coleoptera
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

#56, 54, 58, 61, 59
#KMEDOIDS: PAM

#dna family
PAMajo0f80 <- pamLaskentaF(cole0F80,cole0Fmeta, k = 56)
PAMajo1f80 <- pamLaskentaF(cole1F80,cole1Fmeta, k = 54)
PAMajo2f80 <- pamLaskentaF(cole2F80,cole2Fmeta, k = 58)
PAMajo3f80 <- pamLaskentaF(cole3F80,cole3Fmeta, k = 61)
PAMajo4f80 <- pamLaskentaF(cole4F80,cole4Fmeta, k = 59)

#piirteet family
PAMajo0f <- pamLaskentaF(cole0F,cole0Fmeta, k = 56)
PAMajo1f <- pamLaskentaF(cole1F,cole1Fmeta, k = 54)
PAMajo2f <- pamLaskentaF(cole2F,cole2Fmeta, k = 58)
PAMajo3f <- pamLaskentaF(cole3F,cole3Fmeta, k = 61)
PAMajo4f <- pamLaskentaF(cole4F,cole4Fmeta, k = 59)

#tulosten automaattiluku uusifu:ssa

#indeksit family
indeksilaskuri(PAMajo0f80$a, PAMajo0f$a);print("######### B ALKAA #########"); indeksilaskuri(PAMajo0f80$b, PAMajo0f$b)
indeksilaskuri(PAMajo1f80$a, PAMajo1f$a);print("######### B ALKAA #########"); indeksilaskuri(PAMajo1f80$b, PAMajo1f$b)
indeksilaskuri(PAMajo2f80$a, PAMajo2f$a);print("######### B ALKAA #########"); indeksilaskuri(PAMajo2f80$b, PAMajo2f$b)
indeksilaskuri(PAMajo3f80$a, PAMajo3f$a);print("######### B ALKAA #########"); indeksilaskuri(PAMajo3f80$b, PAMajo3f$b)
indeksilaskuri(PAMajo4f80$a, PAMajo4f$a);print("######### B ALKAA #########"); indeksilaskuri(PAMajo4f80$b, PAMajo4f$b)


#hierarkkinen jakava: DIANA

#dna family
Dia0F80 <- dianalaskentaF(cole0F80,cole0Fmeta, k = 56)
Dia1F80 <- dianalaskentaF(cole1F80,cole1Fmeta, k = 54)
Dia2F80 <- dianalaskentaF(cole2F80,cole2Fmeta, k = 58)
Dia3F80 <- dianalaskentaF(cole3F80,cole3Fmeta, k = 61)
Dia4F80 <- dianalaskentaF(cole4F80,cole4Fmeta, k = 59)

#piirteet family
Dia0F <- dianalaskentaF(cole0F,cole0Fmeta, k = 56)
Dia1F <- dianalaskentaF(cole1F,cole1Fmeta, k = 54)
Dia2F <- dianalaskentaF(cole2F,cole2Fmeta, k = 58)
Dia3F <- dianalaskentaF(cole3F,cole3Fmeta, k = 61)
Dia4F <- dianalaskentaF(cole4F,cole4Fmeta, k = 59)

#tulosten automaattiluku uusifu:ssa

#indeksit
indeksilaskuri(Dia0F80$a,Dia0F$a);print("######### B ALKAA #########");indeksilaskuri(Dia0F80$b,Dia0F$b)
indeksilaskuri(Dia1F80$a,Dia1F$a);print("######### B ALKAA #########");indeksilaskuri(Dia1F80$b,Dia1F$b)
indeksilaskuri(Dia2F80$a,Dia2F$a);print("######### B ALKAA #########");indeksilaskuri(Dia2F80$b,Dia2F$b)
indeksilaskuri(Dia3F80$a,Dia3F$a);print("######### B ALKAA #########");indeksilaskuri(Dia3F80$b,Dia3F$b)
indeksilaskuri(Dia4F80$a,Dia4F$a);print("######### B ALKAA #########");indeksilaskuri(Dia4F80$b,Dia4F$b)


#hierarkkinen yhdistävä: HCLUST
options(max.print=2500)
#DNA

D080 <- as.dist(cole0F80)
D180 <- as.dist(cole1F80)
D280 <- as.dist(cole2F80)
D380 <- as.dist(cole3F80)
D480 <- as.dist(cole4F80)

ajo0_80 <- hclustlaskentaF(D080,cole0Fmeta, k = 56 )
#ajo0_80
ajo1_80 <- hclustlaskentaF(D180,cole1Fmeta, k = 54 )
ajo2_80 <-hclustlaskentaF(D280,cole2Fmeta, k = 58 )
ajo3_80 <-hclustlaskentaF(D380,cole3Fmeta, k = 61 )
ajo4_80 <-hclustlaskentaF(D480,cole4Fmeta, k = 59 )

#PIIRTEET
cole0F <- as.dist(cole0F)
cole1F <- as.dist(cole1F)
cole2F <- as.dist(cole2F)
cole3F <- as.dist(cole3F)
cole4F <- as.dist(cole4F)

ajo0f <- hclustlaskentaF(cole0F,cole0Fmeta, k = 56 )
ajo1f <- hclustlaskentaF(cole1F,cole1Fmeta, k = 54 )
ajo2f <- hclustlaskentaF(cole2F,cole2Fmeta, k = 58 )
ajo3f <- hclustlaskentaF(cole3F,cole3Fmeta, k = 61 )
ajo4f <- hclustlaskentaF(cole4F,cole4Fmeta, k = 59 )

#tulosten automaattiluku uusifu:ssa

#family
indeksilaskuri(ajo0_80$a, ajo0f$a);print("######### B ALKAA #########");indeksilaskuri(ajo0_80$b, ajo0f$b);print("######### C ALKAA #########");indeksilaskuri(ajo0_80$c, ajo0f$c)
indeksilaskuri(ajo1_80$a, ajo1f$a);print("######### B ALKAA #########");indeksilaskuri(ajo1_80$b, ajo1f$b);print("######### C ALKAA #########");indeksilaskuri(ajo1_80$c, ajo1f$c)
indeksilaskuri(ajo2_80$a, ajo2f$a);print("######### B ALKAA #########");indeksilaskuri(ajo2_80$b, ajo2f$b);print("######### C ALKAA #########");indeksilaskuri(ajo2_80$c, ajo2f$c)
indeksilaskuri(ajo3_80$a, ajo3f$a);print("######### B ALKAA #########");indeksilaskuri(ajo3_80$b, ajo3f$b);print("######### C ALKAA #########");indeksilaskuri(ajo3_80$c, ajo3f$c)
indeksilaskuri(ajo4_80$a, ajo4f$a);print("######### B ALKAA #########");indeksilaskuri(ajo4_80$b, ajo4f$b);print("######### C ALKAA #########");indeksilaskuri(ajo4_80$c, ajo4f$c)

