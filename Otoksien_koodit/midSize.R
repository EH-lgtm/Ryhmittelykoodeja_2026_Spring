#gradu koodi MIDSIZE
library(MASS)
library(cluster)
library(pheatmap)
library(mclust)
library(arrow)
library(fossil)
library(aricode)
library(dendextend)

#####ANALYYSI ALKAA#####

#k = #51, 49, 49, 51, 55 kun pienet mukana, fam 


#KMEDOIDS: PAM
#51, 49, 49, 51, 55 
#dna family
PAMajo0f80 <- pamLaskentaF(mid0F80,mid0Fmeta, k = 51)
PAMajo1f80 <- pamLaskentaF(mid1F80,mid1Fmeta, k = 49)
PAMajo2f80 <- pamLaskentaF(mid2F80,mid2Fmeta, k = 49)
PAMajo3f80 <- pamLaskentaF(mid3F80,mid3Fmeta, k = 51)
PAMajo4f80 <- pamLaskentaF(mid4F80,mid4Fmeta, k = 55)

#piirteet family
PAMajo0f <- pamLaskentaF(mid0F,mid0Fmeta, k = 51)
PAMajo1f <- pamLaskentaF(mid1F,mid1Fmeta, k = 49)
PAMajo2f <- pamLaskentaF(mid2F,mid2Fmeta, k = 49)
PAMajo3f <- pamLaskentaF(mid3F,mid3Fmeta, k = 51)
PAMajo4f <- pamLaskentaF(mid4F,mid4Fmeta, k = 55)

#indeksit family
indeksilaskuri(PAMajo0f80$a, PAMajo0f$a);print("######### B ALKAA #########"); indeksilaskuri(PAMajo0f80$b, PAMajo0f$b)
indeksilaskuri(PAMajo1f80$a, PAMajo1f$a);print("######### B ALKAA #########"); indeksilaskuri(PAMajo1f80$b, PAMajo1f$b)
indeksilaskuri(PAMajo2f80$a, PAMajo2f$a);print("######### B ALKAA #########"); indeksilaskuri(PAMajo2f80$b, PAMajo2f$b)
indeksilaskuri(PAMajo3f80$a, PAMajo3f$a);print("######### B ALKAA #########"); indeksilaskuri(PAMajo3f80$b, PAMajo3f$b)
indeksilaskuri(PAMajo4f80$a, PAMajo4f$a);print("######### B ALKAA #########"); indeksilaskuri(PAMajo4f80$b, PAMajo4f$b)

#tulosten auto luku tiedostossa uusifu

#dna order
PAMajo0OR80 <- pamLaskentaOR(mid0OR80,mid0ORmeta, k = 6)
PAMajo1OR80 <- pamLaskentaOR(mid1OR80,mid1ORmeta, k = 6)
PAMajo2OR80 <- pamLaskentaOR(mid2OR80,mid2ORmeta, k = 6)
PAMajo3OR80 <- pamLaskentaOR(mid3OR80,mid3ORmeta, k = 6)
PAMajo4OR80 <- pamLaskentaOR(mid4OR80,mid4ORmeta, k = 6)

#piirteet
PAMajo0OR <- pamLaskentaOR(mid0OR,mid0ORmeta, k = 6)
PAMajo1OR <- pamLaskentaOR(mid1OR,mid1ORmeta, k = 6)
PAMajo2OR <- pamLaskentaOR(mid2OR,mid2ORmeta, k = 6)
PAMajo3OR <- pamLaskentaOR(mid3OR,mid3ORmeta, k = 6)
PAMajo4OR <- pamLaskentaOR(mid4OR,mid4ORmeta, k = 6)

#indeksit ###########################################################
#order
indeksilaskuri(PAMajo0OR80$a, PAMajo0OR$a);print("######### B ALKAA #########");indeksilaskuri(PAMajo0OR80$b, PAMajo0OR$b)
indeksilaskuri(PAMajo1OR80$a, PAMajo1OR$a);print("######### B ALKAA #########");indeksilaskuri(PAMajo1OR80$b, PAMajo1OR$b)
indeksilaskuri(PAMajo2OR80$a, PAMajo2OR$a);print("######### B ALKAA #########");indeksilaskuri(PAMajo2OR80$b, PAMajo2OR$b)
indeksilaskuri(PAMajo3OR80$a, PAMajo3OR$a);print("######### B ALKAA #########");indeksilaskuri(PAMajo3OR80$b, PAMajo3OR$b)
indeksilaskuri(PAMajo4OR80$a, PAMajo4OR$a);print("######### B ALKAA #########");indeksilaskuri(PAMajo4OR80$b, PAMajo4OR$b)

#tulosten auto luku tiedostossa uusifu


#hierarkkinen jakava: DIANA
#51, 49, 49, 51, 55 
#dna family
Dia0F80 <- dianalaskentaF(mid0F80,mid0Fmeta, k = 51)
Dia1F80 <- dianalaskentaF(mid1F80,mid1Fmeta, k = 49)
Dia2F80 <- dianalaskentaF(mid2F80,mid2Fmeta, k = 49)
Dia3F80 <- dianalaskentaF(mid3F80,mid3Fmeta, k = 51)
Dia4F80 <- dianalaskentaF(mid4F80,mid4Fmeta, k = 55)

#piirteet family
Dia0F <- dianalaskentaF(mid0F,mid0Fmeta, k = 51)
Dia1F <- dianalaskentaF(mid1F,mid1Fmeta, k = 49)
Dia2F <- dianalaskentaF(mid2F,mid2Fmeta, k = 49)
Dia3F <- dianalaskentaF(mid3F,mid3Fmeta, k = 51)
Dia4F <- dianalaskentaF(mid4F,mid4Fmeta, k = 55)

#indeksit
indeksilaskuri(Dia0F80$a,Dia0F$a);print("######### B ALKAA #########");indeksilaskuri(Dia0F80$b,Dia0F$b)
indeksilaskuri(Dia1F80$a,Dia1F$a);print("######### B ALKAA #########");indeksilaskuri(Dia1F80$b,Dia1F$b)
indeksilaskuri(Dia2F80$a,Dia2F$a);print("######### B ALKAA #########");indeksilaskuri(Dia2F80$b,Dia2F$b)
indeksilaskuri(Dia3F80$a,Dia3F$a);print("######### B ALKAA #########");indeksilaskuri(Dia3F80$b,Dia3F$b)
indeksilaskuri(Dia4F80$a,Dia4F$a);print("######### B ALKAA #########");indeksilaskuri(Dia4F80$b,Dia4F$b)

#tulosten auto luku tiedostossa uusifu

#order
#6,12 

#dna order
Dia0OR80 <- dianalaskentaOR(mid0OR80,mid0ORmeta, k = 6)
Dia1OR80 <- dianalaskentaOR(mid1OR80,mid1ORmeta, k = 6)
Dia2OR80 <- dianalaskentaOR(mid2OR80,mid2ORmeta, k = 6)
Dia3OR80 <- dianalaskentaOR(mid3OR80,mid3ORmeta, k = 6)
Dia4OR80 <- dianalaskentaOR(mid4OR80,mid4ORmeta, k = 6)

#piirteet order
Dia0OR <- dianalaskentaOR(mid0OR,mid0ORmeta, k = 6)
Dia1OR <- dianalaskentaOR(mid1OR,mid1ORmeta, k = 6)
Dia2OR <- dianalaskentaOR(mid2OR,mid2ORmeta, k = 6)
Dia3OR <- dianalaskentaOR(mid3OR,mid3ORmeta, k = 6)
Dia4OR <- dianalaskentaOR(mid4OR,mid4ORmeta, k = 6)

#indeksit order
indeksilaskuri(Dia0OR80$a,Dia0OR$a);print("######### B ALKAA #########");indeksilaskuri(Dia0OR80$b,Dia0OR$b)
indeksilaskuri(Dia1OR80$a,Dia1OR$a);print("######### B ALKAA #########");indeksilaskuri(Dia1OR80$b,Dia1OR$b)
indeksilaskuri(Dia2OR80$a,Dia2OR$a);print("######### B ALKAA #########");indeksilaskuri(Dia2OR80$b,Dia2OR$b)
indeksilaskuri(Dia3OR80$a,Dia3OR$a);print("######### B ALKAA #########");indeksilaskuri(Dia3OR80$b,Dia3OR$b)
indeksilaskuri(Dia4OR80$a,Dia4OR$a);print("######### B ALKAA #########");indeksilaskuri(Dia4OR80$b,Dia4OR$b)

#hierarkkinen yhdistävä: HCLUST
options(max.print=3000)
#DNA

D080 <- as.dist(mid0F80)
D180 <- as.dist(mid1F80)
D280 <- as.dist(mid2F80)
D380 <- as.dist(mid3F80)
D480 <- as.dist(mid4F80)

ajo0_80 <- hclustlaskentaF(D080,mid0Fmeta, k = 51 )
#ajo0_80
ajo1_80 <- hclustlaskentaF(D180,mid1Fmeta, k = 49 )
ajo2_80 <-hclustlaskentaF(D280,mid2Fmeta, k = 49 )
ajo3_80 <-hclustlaskentaF(D380,mid3Fmeta, k = 51 )
ajo4_80 <-hclustlaskentaF(D480,mid4Fmeta, k = 55 )

#PIIRTEET
mid0F <- as.dist(mid0F)
mid1F <- as.dist(mid1F)
mid2F <- as.dist(mid2F)
mid3F <- as.dist(mid3F)
mid4F <- as.dist(mid4F)

ajo0f <- hclustlaskentaF(mid0F,mid0Fmeta, k = 51 )
ajo1f <- hclustlaskentaF(mid1F,mid1Fmeta, k = 49 )
ajo2f <- hclustlaskentaF(mid2F,mid2Fmeta, k = 49 )
ajo3f <- hclustlaskentaF(mid3F,mid3Fmeta, k = 51 )
ajo4f <- hclustlaskentaF(mid4F,mid4Fmeta, k = 52 )

#tulosten auto luku tiedostossa uusifu

#order
D080OR <- as.dist(mid0OR80)
D180OR <- as.dist(mid1OR80)
D280OR <- as.dist(mid2OR80)
D380OR <- as.dist(mid3OR80)
D480OR <- as.dist(mid4OR80)

ajo080or <- hclustlaskentaOR(D080OR,mid0ORmeta, k = 6 )
ajo180or <- hclustlaskentaOR(D180OR,mid1ORmeta, k = 6 )
ajo280or <- hclustlaskentaOR(D280OR,mid2ORmeta, k = 6 )
ajo380or <- hclustlaskentaOR(D380OR,mid3ORmeta, k = 6 )
ajo480or <- hclustlaskentaOR(D480OR,mid4ORmeta, k = 6 )

#PIIRTEET
mid0OR <- as.dist(mid0OR)
mid1OR <- as.dist(mid1OR)
mid2OR <- as.dist(mid2OR)
mid3OR <- as.dist(mid3OR)
mid4OR <- as.dist(mid4OR)

ajo0OR <- hclustlaskentaOR(mid0OR,mid0ORmeta, k = 6 )
ajo1OR <- hclustlaskentaOR(mid1OR,mid1ORmeta, k = 6 )
ajo2OR <- hclustlaskentaOR(mid2OR,mid2ORmeta, k = 6 )
ajo3OR <- hclustlaskentaOR(mid3OR,mid3ORmeta, k = 6 )
ajo4OR <- hclustlaskentaOR(mid4OR,mid4ORmeta, k = 6 )

#tulosten auto luku tiedostossa uusifu


#family
indeksilaskuri(ajo0_80$a, ajo0f$a);print("######### B ALKAA #########");indeksilaskuri(ajo0_80$b, ajo0f$b);print("######### C ALKAA #########");indeksilaskuri(ajo0_80$c, ajo0f$c)
indeksilaskuri(ajo1_80$a, ajo1f$a);print("######### B ALKAA #########");indeksilaskuri(ajo1_80$b, ajo1f$b);print("######### C ALKAA #########");indeksilaskuri(ajo1_80$c, ajo1f$c)
indeksilaskuri(ajo2_80$a, ajo2f$a);print("######### B ALKAA #########");indeksilaskuri(ajo2_80$b, ajo2f$b);print("######### C ALKAA #########");indeksilaskuri(ajo2_80$c, ajo2f$c)
indeksilaskuri(ajo3_80$a, ajo3f$a);print("######### B ALKAA #########");indeksilaskuri(ajo3_80$b, ajo3f$b);print("######### C ALKAA #########");indeksilaskuri(ajo3_80$c, ajo3f$c)
indeksilaskuri(ajo4_80$a, ajo4f$a);print("######### B ALKAA #########");indeksilaskuri(ajo4_80$b, ajo4f$b);print("######### C ALKAA #########");indeksilaskuri(ajo4_80$c, ajo4f$c)
#order
indeksilaskuri(ajo080or$a, ajo0OR$a);print("######### B ALKAA #########");indeksilaskuri(ajo080or$b, ajo0OR$b);print("######### C ALKAA #########");indeksilaskuri(ajo080or$c, ajo0OR$c)
indeksilaskuri(ajo180or$a, ajo1OR$a);print("######### B ALKAA #########");indeksilaskuri(ajo180or$b, ajo1OR$b);print("######### C ALKAA #########");indeksilaskuri(ajo180or$c, ajo1OR$c)
indeksilaskuri(ajo280or$a, ajo2OR$a);print("######### B ALKAA #########");indeksilaskuri(ajo280or$b, ajo2OR$b);print("######### C ALKAA #########");indeksilaskuri(ajo280or$c, ajo2OR$c)
indeksilaskuri(ajo380or$a, ajo3OR$a);print("######### B ALKAA #########");indeksilaskuri(ajo380or$b, ajo3OR$b);print("######### C ALKAA #########");indeksilaskuri(ajo380or$c, ajo3OR$c)
indeksilaskuri(ajo480or$a, ajo4OR$a);print("######### B ALKAA #########");indeksilaskuri(ajo480or$b, ajo4OR$b);print("######### C ALKAA #########");indeksilaskuri(ajo480or$c, ajo4OR$c)


