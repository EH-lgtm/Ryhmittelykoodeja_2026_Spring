#gradu koodi hemiptera
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

#k = 36, 38, 36, 40, 33 kun pienet mukana
#KMEDOIDS: PAM

#dna family
PAMajo0f80 <- pamLaskentaF(hemi0F80,hemi0Fmeta, k = 36)
PAMajo1f80 <- pamLaskentaF(hemi1F80,hemi1Fmeta, k = 38)
PAMajo2f80 <- pamLaskentaF(hemi2F80,hemi2Fmeta, k = 36)
PAMajo3f80 <- pamLaskentaF(hemi3F80,hemi3Fmeta, k = 40)
PAMajo4f80 <- pamLaskentaF(hemi4F80,hemi4Fmeta, k = 33)

#piirteet family
PAMajo0f <- pamLaskentaF(hemi0F,hemi0Fmeta, k = 36)
PAMajo1f <- pamLaskentaF(hemi1F,hemi1Fmeta, k = 38)
PAMajo2f <- pamLaskentaF(hemi2F,hemi2Fmeta, k = 36)
PAMajo3f <- pamLaskentaF(hemi3F,hemi3Fmeta, k = 40)
PAMajo4f <- pamLaskentaF(hemi4F,hemi4Fmeta, k = 33)

#tulosten autoluku uusifu:ssa

#indeksit family
indeksilaskuri(PAMajo0f80$a, PAMajo0f$a);print("######### B ALKAA #########"); indeksilaskuri(PAMajo0f80$b, PAMajo0f$b)
indeksilaskuri(PAMajo1f80$a, PAMajo1f$a);print("######### B ALKAA #########"); indeksilaskuri(PAMajo1f80$b, PAMajo1f$b)
indeksilaskuri(PAMajo2f80$a, PAMajo2f$a);print("######### B ALKAA #########"); indeksilaskuri(PAMajo2f80$b, PAMajo2f$b)
indeksilaskuri(PAMajo3f80$a, PAMajo3f$a);print("######### B ALKAA #########"); indeksilaskuri(PAMajo3f80$b, PAMajo3f$b)
indeksilaskuri(PAMajo4f80$a, PAMajo4f$a);print("######### B ALKAA #########"); indeksilaskuri(PAMajo4f80$b, PAMajo4f$b)


#hierarkkinen jakava: DIANA

#dna family
Dia0F80 <- dianalaskentaF(hemi0F80,hemi0Fmeta, k = 36)
Dia1F80 <- dianalaskentaF(hemi1F80,hemi1Fmeta, k = 38)
Dia2F80 <- dianalaskentaF(hemi2F80,hemi2Fmeta, k = 36)
Dia3F80 <- dianalaskentaF(hemi3F80,hemi3Fmeta, k = 40)
Dia4F80 <- dianalaskentaF(hemi4F80,hemi4Fmeta, k = 33)

#piirteet family
Dia0F <- dianalaskentaF(hemi0F,hemi0Fmeta, k = 36)
Dia1F <- dianalaskentaF(hemi1F,hemi1Fmeta, k = 38)
Dia2F <- dianalaskentaF(hemi2F,hemi2Fmeta, k = 36)
Dia3F <- dianalaskentaF(hemi3F,hemi3Fmeta, k = 40)
Dia4F <- dianalaskentaF(hemi4F,hemi4Fmeta, k = 33)

#tulosten autoluku uusifu:ssa

#indeksit
indeksilaskuri(Dia0F80$a,Dia0F$a);print("######### B ALKAA #########");indeksilaskuri(Dia0F80$b,Dia0F$b)
indeksilaskuri(Dia1F80$a,Dia1F$a);print("######### B ALKAA #########");indeksilaskuri(Dia1F80$b,Dia1F$b)
indeksilaskuri(Dia2F80$a,Dia2F$a);print("######### B ALKAA #########");indeksilaskuri(Dia2F80$b,Dia2F$b)
indeksilaskuri(Dia3F80$a,Dia3F$a);print("######### B ALKAA #########");indeksilaskuri(Dia3F80$b,Dia3F$b)
indeksilaskuri(Dia4F80$a,Dia4F$a);print("######### B ALKAA #########");indeksilaskuri(Dia4F80$b,Dia4F$b)



#hierarkkinen yhdistävä: HCLUST
options(max.print=2500)
#DNA

D080 <- as.dist(hemi0F80)
D180 <- as.dist(hemi1F80)
D280 <- as.dist(hemi2F80)
D380 <- as.dist(hemi3F80)
D480 <- as.dist(hemi4F80)

ajo0_80 <- hclustlaskentaF(D080,hemi0Fmeta, k = 36 )
#ajo0_80$D
ajo1_80 <- hclustlaskentaF(D180,hemi1Fmeta, k = 38 )
ajo2_80 <-hclustlaskentaF(D280,hemi2Fmeta, k = 36 )
ajo3_80 <-hclustlaskentaF(D380,hemi3Fmeta, k = 40 )
ajo4_80 <-hclustlaskentaF(D480,hemi4Fmeta, k = 33 )

#PIIRTEET
hemi0F <- as.dist(hemi0F)
hemi1F <- as.dist(hemi1F)
hemi2F <- as.dist(hemi2F)
hemi3F <- as.dist(hemi3F)
hemi4F <- as.dist(hemi4F)

ajo0f <- hclustlaskentaF(hemi0F,hemi0Fmeta, k = 36 )
ajo1f <- hclustlaskentaF(hemi1F,hemi1Fmeta, k = 38 )
ajo2f <- hclustlaskentaF(hemi2F,hemi2Fmeta, k = 36 )
ajo3f <- hclustlaskentaF(hemi3F,hemi3Fmeta, k = 40 )
ajo4f <- hclustlaskentaF(hemi4F,hemi4Fmeta, k = 33 )

#tulosten autoluku uusifu:ssa

#family
indeksilaskuri(ajo0_80$a, ajo0f$a);print("######### B ALKAA #########");indeksilaskuri(ajo0_80$b, ajo0f$b);print("######### C ALKAA #########");indeksilaskuri(ajo0_80$c, ajo0f$c)
indeksilaskuri(ajo1_80$a, ajo1f$a);print("######### B ALKAA #########");indeksilaskuri(ajo1_80$b, ajo1f$b);print("######### C ALKAA #########");indeksilaskuri(ajo1_80$c, ajo1f$c)
indeksilaskuri(ajo2_80$a, ajo2f$a);print("######### B ALKAA #########");indeksilaskuri(ajo2_80$b, ajo2f$b);print("######### C ALKAA #########");indeksilaskuri(ajo2_80$c, ajo2f$c)
indeksilaskuri(ajo3_80$a, ajo3f$a);print("######### B ALKAA #########");indeksilaskuri(ajo3_80$b, ajo3f$b);print("######### C ALKAA #########");indeksilaskuri(ajo3_80$c, ajo3f$c)
indeksilaskuri(ajo4_80$a, ajo4f$a);print("######### B ALKAA #########");indeksilaskuri(ajo4_80$b, ajo4f$b);print("######### C ALKAA #########");indeksilaskuri(ajo4_80$c, ajo4f$c)

