#gradu koodi LARGESIZE
library(MASS)
library(cluster)
library(pheatmap)
library(mclust)
library(arrow)
library(fossil)
library(aricode)
library(dendextend)


#####ANALYYSI ALKAA#####

#KMEDOIDS: PAM
#119, 119, 117, 126, 117 Family
#dna family
PAMajo0f80 <- pamLaskentaF(large0F80,large0Fmeta, k = 119)
#PAMajo0f80$pam1F
#PAMajo0f80$pam2F
#PAMajo0f80$pam1F==PAMajo0f80$pam2F ihan samat
PAMajo1f80 <- pamLaskentaF(large1F80,large1Fmeta, k = 119)
PAMajo2f80 <- pamLaskentaF(large2F80,large2Fmeta, k = 117)
PAMajo3f80 <- pamLaskentaF(large3F80,large3Fmeta, k = 126)
PAMajo4f80 <- pamLaskentaF(large4F80,large4Fmeta, k = 117)

#piirteet family
PAMajo0f <- pamLaskentaF(large0F,large0Fmeta, k = 119)
PAMajo1f <- pamLaskentaF(large1F,large1Fmeta, k = 119)
PAMajo2f <- pamLaskentaF(large2F,large2Fmeta, k = 117)
PAMajo3f <- pamLaskentaF(large3F,large3Fmeta, k = 126)
PAMajo4f <- pamLaskentaF(large4F,large4Fmeta, k = 117)

#tuloksien auto luku tiedostossa uusifu 

#indeksit family
indeksilaskuri(PAMajo0f80$a, PAMajo0f$a);print("######### B ALKAA #########"); indeksilaskuri(PAMajo0f80$b, PAMajo0f$b)
indeksilaskuri(PAMajo1f80$a, PAMajo1f$a);print("######### B ALKAA #########"); indeksilaskuri(PAMajo1f80$b, PAMajo1f$b)
indeksilaskuri(PAMajo2f80$a, PAMajo2f$a);print("######### B ALKAA #########"); indeksilaskuri(PAMajo2f80$b, PAMajo2f$b)
indeksilaskuri(PAMajo3f80$a, PAMajo3f$a);print("######### B ALKAA #########"); indeksilaskuri(PAMajo3f80$b, PAMajo3f$b)
indeksilaskuri(PAMajo4f80$a, PAMajo4f$a);print("######### B ALKAA #########"); indeksilaskuri(PAMajo4f80$b, PAMajo4f$b)


###

#dna order
PAMajo0OR80 <- pamLaskentaOR(large0OR80,large0ORmeta, k = 5)
PAMajo1OR80 <- pamLaskentaOR(large1OR80,large1ORmeta, k = 5)
PAMajo2OR80 <- pamLaskentaOR(large2OR80,large2ORmeta, k = 5)
PAMajo3OR80 <- pamLaskentaOR(large3OR80,large3ORmeta, k = 5)
PAMajo4OR80 <- pamLaskentaOR(large4OR80,large4ORmeta, k = 5)
#PAMajo0OR80$pam1Or
#piirteet
PAMajo0OR <- pamLaskentaOR(large0OR,large0ORmeta, k = 5)
PAMajo1OR <- pamLaskentaOR(large1OR,large1ORmeta, k = 5)
PAMajo2OR <- pamLaskentaOR(large2OR,large2ORmeta, k = 5)
PAMajo3OR <- pamLaskentaOR(large3OR,large3ORmeta, k = 5)
PAMajo4OR <- pamLaskentaOR(large4OR,large4ORmeta, k = 5)

#tuloksien auto luku tiedostossa uusifu 

#indeksit ###########################################################
#order
indeksilaskuri(PAMajo0OR80$a, PAMajo0OR$a);print("######### B ALKAA #########");indeksilaskuri(PAMajo0OR80$b, PAMajo0OR$b)
indeksilaskuri(PAMajo1OR80$a, PAMajo1OR$a);print("######### B ALKAA #########");indeksilaskuri(PAMajo1OR80$b, PAMajo1OR$b)
indeksilaskuri(PAMajo2OR80$a, PAMajo2OR$a);print("######### B ALKAA #########");indeksilaskuri(PAMajo2OR80$b, PAMajo2OR$b)
indeksilaskuri(PAMajo3OR80$a, PAMajo3OR$a);print("######### B ALKAA #########");indeksilaskuri(PAMajo3OR80$b, PAMajo3OR$b)
indeksilaskuri(PAMajo4OR80$a, PAMajo4OR$a);print("######### B ALKAA #########");indeksilaskuri(PAMajo4OR80$b, PAMajo4OR$b)

####
#hierarkkinen jakava: DIANA
#119, 119, 117, 126, 117 Family
#dna family
Dia0F80 <- dianalaskentaF(large0F80,large0Fmeta, k = 119)
Dia1F80 <- dianalaskentaF(large1F80,large1Fmeta, k = 119)
Dia2F80 <- dianalaskentaF(large2F80,large2Fmeta, k = 117)
Dia3F80 <- dianalaskentaF(large3F80,large3Fmeta, k = 126)
Dia4F80 <- dianalaskentaF(large4F80,large4Fmeta, k = 117)
#Dia0F80$diana1F==Dia0F80$diana2F samat
#piirteet family
Dia0F <- dianalaskentaF(large0F,large0Fmeta, k = 119)
Dia1F <- dianalaskentaF(large1F,large1Fmeta, k = 119)
Dia2F <- dianalaskentaF(large2F,large2Fmeta, k = 117)
Dia3F <- dianalaskentaF(large3F,large3Fmeta, k = 126)
Dia4F <- dianalaskentaF(large4F,large4Fmeta, k = 117)

#tuloksien auto luku tiedostossa uusifu 

#indeksit
indeksilaskuri(Dia0F80$a,Dia0F$a);print("######### B ALKAA #########");indeksilaskuri(Dia0F80$b,Dia0F$b)
indeksilaskuri(Dia1F80$a,Dia1F$a);print("######### B ALKAA #########");indeksilaskuri(Dia1F80$b,Dia1F$b)
indeksilaskuri(Dia2F80$a,Dia2F$a);print("######### B ALKAA #########");indeksilaskuri(Dia2F80$b,Dia2F$b)
indeksilaskuri(Dia3F80$a,Dia3F$a);print("######### B ALKAA #########");indeksilaskuri(Dia3F80$b,Dia3F$b)
indeksilaskuri(Dia4F80$a,Dia4F$a);print("######### B ALKAA #########");indeksilaskuri(Dia4F80$b,Dia4F$b)

#####
#order

#dna order
Dia0OR80 <- dianalaskentaOR(large0OR80,large0ORmeta, k = 5)
Dia1OR80 <- dianalaskentaOR(large1OR80,large1ORmeta, k = 5)
Dia2OR80 <- dianalaskentaOR(large2OR80,large2ORmeta, k = 5)
Dia3OR80 <- dianalaskentaOR(large3OR80,large3ORmeta, k = 5)
Dia4OR80 <- dianalaskentaOR(large4OR80,large4ORmeta, k = 5)
#Dia0OR80$diana1Or==Dia0OR80$diana2Or
#piirteet order
Dia0OR <- dianalaskentaOR(large0OR,large0ORmeta, k = 5)
Dia1OR <- dianalaskentaOR(large1OR,large1ORmeta, k = 5)
Dia2OR <- dianalaskentaOR(large2OR,large2ORmeta, k = 5)
Dia3OR <- dianalaskentaOR(large3OR,large3ORmeta, k = 5)
Dia4OR <- dianalaskentaOR(large4OR,large4ORmeta, k = 5)

#tuloksien auto luku tiedostossa uusifu 

#indeksit order
indeksilaskuri(Dia0OR80$a,Dia0OR$a);print("######### B ALKAA #########");indeksilaskuri(Dia0OR80$b,Dia0OR$b)
indeksilaskuri(Dia1OR80$a,Dia1OR$a);print("######### B ALKAA #########");indeksilaskuri(Dia1OR80$b,Dia1OR$b)
indeksilaskuri(Dia2OR80$a,Dia2OR$a);print("######### B ALKAA #########");indeksilaskuri(Dia2OR80$b,Dia2OR$b)
indeksilaskuri(Dia3OR80$a,Dia3OR$a);print("######### B ALKAA #########");indeksilaskuri(Dia3OR80$b,Dia3OR$b)
indeksilaskuri(Dia4OR80$a,Dia4OR$a);print("######### B ALKAA #########");indeksilaskuri(Dia4OR80$b,Dia4OR$b)




#hierarkkinen yhdistävä: HCLUST
options(max.print=2500)
#DNA #119, 119, 117, 126, 117 Family

D080 <- as.dist(large0F80)
D180 <- as.dist(large1F80)
D280 <- as.dist(large2F80)
D380 <- as.dist(large3F80)
D480 <- as.dist(large4F80)

ajo0_80 <- hclustlaskentaF(D080,large0Fmeta, k = 119 )
#ajo0_80
ajo1_80 <- hclustlaskentaF(D180,large1Fmeta, k = 119 )
ajo2_80 <-hclustlaskentaF(D280,large2Fmeta, k = 117 )
ajo3_80 <-hclustlaskentaF(D380,large3Fmeta, k = 126 )
ajo4_80 <-hclustlaskentaF(D480,large4Fmeta, k = 117 )

#PIIRTEET
large0F <- as.dist(large0F)
large1F <- as.dist(large1F)
large2F <- as.dist(large2F)
large3F <- as.dist(large3F)
large4F <- as.dist(large4F)

ajo0f <- hclustlaskentaF(large0F,large0Fmeta, k = 119 )
ajo1f <- hclustlaskentaF(large1F,large1Fmeta, k = 119 )
ajo2f <- hclustlaskentaF(large2F,large2Fmeta, k = 117 )
ajo3f <- hclustlaskentaF(large3F,large3Fmeta, k = 126 )
ajo4f <- hclustlaskentaF(large4F,large4Fmeta, k = 117 )
######
#HUOMIO:
# tulosten autoluku tiedostossa: uusifu
######

#indeksit
#family
indeksilaskuri(ajo0_80$a, ajo0f$a);print("######### B ALKAA #########");indeksilaskuri(ajo0_80$b, ajo0f$b);print("######### C ALKAA #########");indeksilaskuri(ajo0_80$c, ajo0f$c)
indeksilaskuri(ajo1_80$a, ajo1f$a);print("######### B ALKAA #########");indeksilaskuri(ajo1_80$b, ajo1f$b);print("######### C ALKAA #########");indeksilaskuri(ajo1_80$c, ajo1f$c)
indeksilaskuri(ajo2_80$a, ajo2f$a);print("######### B ALKAA #########");indeksilaskuri(ajo2_80$b, ajo2f$b);print("######### C ALKAA #########");indeksilaskuri(ajo2_80$c, ajo2f$c)
indeksilaskuri(ajo3_80$a, ajo3f$a);print("######### B ALKAA #########");indeksilaskuri(ajo3_80$b, ajo3f$b);print("######### C ALKAA #########");indeksilaskuri(ajo3_80$c, ajo3f$c)
indeksilaskuri(ajo4_80$a, ajo4f$a);print("######### B ALKAA #########");indeksilaskuri(ajo4_80$b, ajo4f$b);print("######### C ALKAA #########");indeksilaskuri(ajo4_80$c, ajo4f$c)

#####

#DNA order #5 tai joku muu
D080OR <- as.dist(large0OR80)
D180OR <- as.dist(large1OR80)
D280OR <- as.dist(large2OR80)
D380OR <- as.dist(large3OR80)
D480OR <- as.dist(large4OR80)

ajo080or <- hclustlaskentaOR(D080OR,large0ORmeta, k = 5 )
ajo180or <- hclustlaskentaOR(D180OR,large1ORmeta, k = 5 )
ajo280or <- hclustlaskentaOR(D280OR,large2ORmeta, k = 5 )
ajo380or <- hclustlaskentaOR(D380OR,large3ORmeta, k = 5 )
ajo480or <- hclustlaskentaOR(D480OR,large4ORmeta, k = 5 )

#PIIRTEET
large0OR <- as.dist(large0OR)
large1OR <- as.dist(large1OR)
large2OR <- as.dist(large2OR)
large3OR <- as.dist(large3OR)
large4OR <- as.dist(large4OR)

ajo0OR <- hclustlaskentaOR(large0OR,large0ORmeta, k = 5 )
ajo1OR <- hclustlaskentaOR(large1OR,large1ORmeta, k = 5 )
ajo2OR <- hclustlaskentaOR(large2OR,large2ORmeta, k = 5 )
ajo3OR <- hclustlaskentaOR(large3OR,large3ORmeta, k = 5 )
ajo4OR <- hclustlaskentaOR(large4OR,large4ORmeta, k = 5 )

######
#HUOMIO:
# tulosten autoluku tiedostossa: uusifu
######

#order
indeksilaskuri(ajo080or$a, ajo0OR$a);print("######### B ALKAA #########");indeksilaskuri(ajo080or$b, ajo0OR$b);print("######### C ALKAA #########");indeksilaskuri(ajo080or$c, ajo0OR$c)
indeksilaskuri(ajo180or$a, ajo1OR$a);print("######### B ALKAA #########");indeksilaskuri(ajo180or$b, ajo1OR$b);print("######### C ALKAA #########");indeksilaskuri(ajo180or$c, ajo1OR$c)
indeksilaskuri(ajo280or$a, ajo2OR$a);print("######### B ALKAA #########");indeksilaskuri(ajo280or$b, ajo2OR$b);print("######### C ALKAA #########");indeksilaskuri(ajo280or$c, ajo2OR$c)
indeksilaskuri(ajo380or$a, ajo3OR$a);print("######### B ALKAA #########");indeksilaskuri(ajo380or$b, ajo3OR$b);print("######### C ALKAA #########");indeksilaskuri(ajo380or$c, ajo3OR$c)
indeksilaskuri(ajo480or$a, ajo4OR$a);print("######### B ALKAA #########");indeksilaskuri(ajo480or$b, ajo4OR$b);print("######### C ALKAA #########");indeksilaskuri(ajo480or$c, ajo4OR$c)

