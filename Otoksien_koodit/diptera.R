#gradu koodi diptera
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

#k = #32, 31, 30, 33, 32 eri tyyppiä samplessa 0-4 FAM 
#KMEDOIDS: PAM

#dna family
PAMajo0f80 <- pamLaskentaF(dip0F80,dip0Fmeta, k = 32)
PAMajo1f80 <- pamLaskentaF(dip1F80,dip1Fmeta, k = 31)
PAMajo2f80 <- pamLaskentaF(dip2F80,dip2Fmeta, k = 30)
PAMajo3f80 <- pamLaskentaF(dip3F80,dip3Fmeta, k = 33)
PAMajo4f80 <- pamLaskentaF(dip4F80,dip4Fmeta, k = 32)

#piirteet family
PAMajo0f <- pamLaskentaF(dip0F,dip0Fmeta, k = 32)
PAMajo1f <- pamLaskentaF(dip1F,dip1Fmeta, k = 31)
PAMajo2f <- pamLaskentaF(dip2F,dip2Fmeta, k = 30)
PAMajo3f <- pamLaskentaF(dip3F,dip3Fmeta, k = 33)
PAMajo4f <- pamLaskentaF(dip4F,dip4Fmeta, k = 32)

#automaattinen tulosten luku uusifu:ssa

#indeksit family
indeksilaskuri(PAMajo0f80$a, PAMajo0f$a);print("######### B ALKAA #########"); indeksilaskuri(PAMajo0f80$b, PAMajo0f$b)
indeksilaskuri(PAMajo1f80$a, PAMajo1f$a);print("######### B ALKAA #########"); indeksilaskuri(PAMajo1f80$b, PAMajo1f$b)
indeksilaskuri(PAMajo2f80$a, PAMajo2f$a);print("######### B ALKAA #########"); indeksilaskuri(PAMajo2f80$b, PAMajo2f$b)
indeksilaskuri(PAMajo3f80$a, PAMajo3f$a);print("######### B ALKAA #########"); indeksilaskuri(PAMajo3f80$b, PAMajo3f$b)
indeksilaskuri(PAMajo4f80$a, PAMajo4f$a);print("######### B ALKAA #########"); indeksilaskuri(PAMajo4f80$b, PAMajo4f$b)

#hierarkkinen jakava: DIANA

#dna family
Dia0F80 <- dianalaskentaF(dip0F80,dip0Fmeta, k = 32)
Dia1F80 <- dianalaskentaF(dip1F80,dip1Fmeta, k = 31)
Dia2F80 <- dianalaskentaF(dip2F80,dip2Fmeta, k = 30)
Dia3F80 <- dianalaskentaF(dip3F80,dip3Fmeta, k = 33)
Dia4F80 <- dianalaskentaF(dip4F80,dip4Fmeta, k = 32)

#piirteet family
Dia0F <- dianalaskentaF(dip0F,dip0Fmeta, k = 32)
Dia1F <- dianalaskentaF(dip1F,dip1Fmeta, k = 31)
Dia2F <- dianalaskentaF(dip2F,dip2Fmeta, k = 30)
Dia3F <- dianalaskentaF(dip3F,dip3Fmeta, k = 33)
Dia4F <- dianalaskentaF(dip4F,dip4Fmeta, k = 32)

#automaattinen tulosten luku uusifu:ssa

#indeksit
indeksilaskuri(Dia0F80$a,Dia0F$a);print("######### B ALKAA #########");indeksilaskuri(Dia0F80$b,Dia0F$b)
indeksilaskuri(Dia1F80$a,Dia1F$a);print("######### B ALKAA #########");indeksilaskuri(Dia1F80$b,Dia1F$b)
indeksilaskuri(Dia2F80$a,Dia2F$a);print("######### B ALKAA #########");indeksilaskuri(Dia2F80$b,Dia2F$b)
indeksilaskuri(Dia3F80$a,Dia3F$a);print("######### B ALKAA #########");indeksilaskuri(Dia3F80$b,Dia3F$b)
indeksilaskuri(Dia4F80$a,Dia4F$a);print("######### B ALKAA #########");indeksilaskuri(Dia4F80$b,Dia4F$b)



#hierarkkinen yhdistävä: HCLUST
options(max.print=2500)
#DNA

D080 <- as.dist(dip0F80)
D180 <- as.dist(dip1F80)
D280 <- as.dist(dip2F80)
D380 <- as.dist(dip3F80)
D480 <- as.dist(dip4F80)

ajo0_80 <- hclustlaskentaF(D080,dip0Fmeta, k = 32 )
#ajo0_80
ajo1_80 <- hclustlaskentaF(D180,dip1Fmeta, k = 31 )
ajo2_80 <-hclustlaskentaF(D280,dip2Fmeta, k = 30 )
ajo3_80 <-hclustlaskentaF(D380,dip3Fmeta, k = 33 )
ajo4_80 <-hclustlaskentaF(D480,dip4Fmeta, k = 32 )

#PIIRTEET
dip0F <- as.dist(dip0F)
dip1F <- as.dist(dip1F)
dip2F <- as.dist(dip2F)
dip3F <- as.dist(dip3F)
dip4F <- as.dist(dip4F)

ajo0f <- hclustlaskentaF(dip0F,dip0Fmeta, k = 32 )
ajo1f <- hclustlaskentaF(dip1F,dip1Fmeta, k = 31 )
ajo2f <- hclustlaskentaF(dip2F,dip2Fmeta, k = 30 )
ajo3f <- hclustlaskentaF(dip3F,dip3Fmeta, k = 33 )
ajo4f <- hclustlaskentaF(dip4F,dip4Fmeta, k = 32 )

#automaattinen tulosten luku uusifu:ssa

##indeksit
#family
indeksilaskuri(ajo0_80$a, ajo0f$a);print("######### B ALKAA #########");indeksilaskuri(ajo0_80$b, ajo0f$b);print("######### C ALKAA #########");indeksilaskuri(ajo0_80$c, ajo0f$c)
indeksilaskuri(ajo1_80$a, ajo1f$a);print("######### B ALKAA #########");indeksilaskuri(ajo1_80$b, ajo1f$b);print("######### C ALKAA #########");indeksilaskuri(ajo1_80$c, ajo1f$c)
indeksilaskuri(ajo2_80$a, ajo2f$a);print("######### B ALKAA #########");indeksilaskuri(ajo2_80$b, ajo2f$b);print("######### C ALKAA #########");indeksilaskuri(ajo2_80$c, ajo2f$c)
indeksilaskuri(ajo3_80$a, ajo3f$a);print("######### B ALKAA #########");indeksilaskuri(ajo3_80$b, ajo3f$b);print("######### C ALKAA #########");indeksilaskuri(ajo3_80$c, ajo3f$c)
indeksilaskuri(ajo4_80$a, ajo4f$a);print("######### B ALKAA #########");indeksilaskuri(ajo4_80$b, ajo4f$b);print("######### C ALKAA #########");indeksilaskuri(ajo4_80$c, ajo4f$c)
