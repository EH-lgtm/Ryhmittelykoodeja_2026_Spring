#gradu koodi alleven
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
##########
options(max.print=2500)
#KMEDOIDS: PAM

#dna family
PAMajo0f80 <- pamLaskentaF(allE0F80,allE0Fmeta, k = 142) 
#PAMajo0f80$pam1F == PAMajo0f80$pam2F #antaa eri tulokset jos diss asetus ei päällä
PAMajo1f80 <- pamLaskentaF(allE1F80,allE1Fmeta, k = 143)
PAMajo2f80 <- pamLaskentaF(allE2F80,allE2Fmeta, k = 144)
PAMajo3f80 <- pamLaskentaF(allE3F80,allE3Fmeta, k = 147)
PAMajo4f80 <- pamLaskentaF(allE4F80,allE4Fmeta, k = 142)

#piirteet family
PAMajo0f <- pamLaskentaF(allE0F,allE0Fmeta, k = 142)
#PAMajo0f$pam1F == PAMajo0f$pam2F
PAMajo1f <- pamLaskentaF(allE1F,allE1Fmeta, k = 143)
PAMajo2f <- pamLaskentaF(allE2F,allE2Fmeta, k = 144)
PAMajo3f <- pamLaskentaF(allE3F,allE3Fmeta, k = 147)
PAMajo4f <- pamLaskentaF(allE4F,allE4Fmeta, k = 142)

#tulosten autoluku uusifu:ssa

#indeksit family
indeksilaskuri(PAMajo0f80$a, PAMajo0f$a);print("######### B ALKAA #########"); indeksilaskuri(PAMajo0f80$b, PAMajo0f$b)
indeksilaskuri(PAMajo1f80$a, PAMajo1f$a);print("######### B ALKAA #########"); indeksilaskuri(PAMajo1f80$b, PAMajo1f$b)
indeksilaskuri(PAMajo2f80$a, PAMajo2f$a);print("######### B ALKAA #########"); indeksilaskuri(PAMajo2f80$b, PAMajo2f$b)
indeksilaskuri(PAMajo3f80$a, PAMajo3f$a);print("######### B ALKAA #########"); indeksilaskuri(PAMajo3f80$b, PAMajo3f$b)
indeksilaskuri(PAMajo4f80$a, PAMajo4f$a);print("######### B ALKAA #########"); indeksilaskuri(PAMajo4f80$b, PAMajo4f$b)

#dna order
PAMajo0OR80 <- pamLaskentaOR(allE0OR80,allE0ORmeta, k = 19)
PAMajo1OR80 <- pamLaskentaOR(allE1OR80,allE1ORmeta, k = 19)
PAMajo2OR80 <- pamLaskentaOR(allE2OR80,allE2ORmeta, k = 19)
PAMajo3OR80 <- pamLaskentaOR(allE3OR80,allE3ORmeta, k = 19)
PAMajo4OR80 <- pamLaskentaOR(allE4OR80,allE4ORmeta, k = 19)

#piirteet
PAMajo0OR <- pamLaskentaOR(allE0OR,allE0ORmeta, k = 19)
PAMajo1OR <- pamLaskentaOR(allE1OR,allE1ORmeta, k = 19)
PAMajo2OR <- pamLaskentaOR(allE2OR,allE2ORmeta, k = 19)
PAMajo3OR <- pamLaskentaOR(allE3OR,allE3ORmeta, k = 19)
PAMajo4OR <- pamLaskentaOR(allE4OR,allE4ORmeta, k = 19)

#tulosten autoluku uusifu:ssa

#indeksit 
#order
indeksilaskuri(PAMajo0OR80$a, PAMajo0OR$a);print("######### B ALKAA #########");indeksilaskuri(PAMajo0OR80$b, PAMajo0OR$b)
indeksilaskuri(PAMajo1OR80$a, PAMajo1OR$a);print("######### B ALKAA #########");indeksilaskuri(PAMajo1OR80$b, PAMajo1OR$b)
indeksilaskuri(PAMajo2OR80$a, PAMajo2OR$a);print("######### B ALKAA #########");indeksilaskuri(PAMajo2OR80$b, PAMajo2OR$b)
indeksilaskuri(PAMajo3OR80$a, PAMajo3OR$a);print("######### B ALKAA #########");indeksilaskuri(PAMajo3OR80$b, PAMajo3OR$b)
indeksilaskuri(PAMajo4OR80$a, PAMajo4OR$a);print("######### B ALKAA #########");indeksilaskuri(PAMajo4OR80$b, PAMajo4OR$b)

#################################################################

#hierarkkinen jakava: DIANA

#dna family
Dia0F80 <- dianalaskentaF(allE0F80,allE0Fmeta, k = 142) 
#Dia0F80$diana1F==Dia0F80$diana2F
Dia1F80 <- dianalaskentaF(allE1F80,allE1Fmeta, k = 143)
Dia2F80 <- dianalaskentaF(allE2F80,allE2Fmeta, k = 144)
Dia3F80 <- dianalaskentaF(allE3F80,allE3Fmeta, k = 147)
Dia4F80 <- dianalaskentaF(allE4F80,allE4Fmeta, k = 142)

#piirteet family
Dia0F <- dianalaskentaF(allE0F,allE0Fmeta, k = 142)
Dia1F <- dianalaskentaF(allE1F,allE1Fmeta, k = 143)
Dia2F <- dianalaskentaF(allE2F,allE2Fmeta, k = 144)
Dia3F <- dianalaskentaF(allE3F,allE3Fmeta, k = 147)
Dia4F <- dianalaskentaF(allE4F,allE4Fmeta, k = 142)

#tulosten autoluku uusifu:ssa

#indeksit
indeksilaskuri(Dia0F80$a,Dia0F$a);print("######### B ALKAA #########");indeksilaskuri(Dia0F80$b,Dia0F$b)
indeksilaskuri(Dia1F80$a,Dia1F$a);print("######### B ALKAA #########");indeksilaskuri(Dia1F80$b,Dia1F$b)
indeksilaskuri(Dia2F80$a,Dia2F$a);print("######### B ALKAA #########");indeksilaskuri(Dia2F80$b,Dia2F$b)
indeksilaskuri(Dia3F80$a,Dia3F$a);print("######### B ALKAA #########");indeksilaskuri(Dia3F80$b,Dia3F$b)
indeksilaskuri(Dia4F80$a,Dia4F$a);print("######### B ALKAA #########");indeksilaskuri(Dia4F80$b,Dia4F$b)

#order

#dna order
Dia0OR80 <- dianalaskentaOR(allE0OR80,allE0ORmeta, k = 19)
Dia1OR80 <- dianalaskentaOR(allE1OR80,allE1ORmeta, k = 19)
Dia2OR80 <- dianalaskentaOR(allE2OR80,allE2ORmeta, k = 19)
Dia3OR80 <- dianalaskentaOR(allE3OR80,allE3ORmeta, k = 19)
Dia4OR80 <- dianalaskentaOR(allE4OR80,allE4ORmeta, k = 19)

#piirteet order
Dia0OR <- dianalaskentaOR(allE0OR,allE0ORmeta, k = 19)
Dia1OR <- dianalaskentaOR(allE1OR,allE1ORmeta, k = 19)
Dia2OR <- dianalaskentaOR(allE2OR,allE2ORmeta, k = 19)
Dia3OR <- dianalaskentaOR(allE3OR,allE3ORmeta, k = 19)
Dia4OR <- dianalaskentaOR(allE4OR,allE4ORmeta, k = 19)

#tulosten autoluku uusifu:ssa

#indeksit order
indeksilaskuri(Dia0OR80$a,Dia0OR$a);print("######### B ALKAA #########");indeksilaskuri(Dia0OR80$b,Dia0OR$b)
indeksilaskuri(Dia1OR80$a,Dia1OR$a);print("######### B ALKAA #########");indeksilaskuri(Dia1OR80$b,Dia1OR$b)
indeksilaskuri(Dia2OR80$a,Dia2OR$a);print("######### B ALKAA #########");indeksilaskuri(Dia2OR80$b,Dia2OR$b)
indeksilaskuri(Dia3OR80$a,Dia3OR$a);print("######### B ALKAA #########");indeksilaskuri(Dia3OR80$b,Dia3OR$b)
indeksilaskuri(Dia4OR80$a,Dia4OR$a);print("######### B ALKAA #########");indeksilaskuri(Dia4OR80$b,Dia4OR$b)

##########################

#hierarkkinen yhdistävä: HCLUST
options(max.print=2500)
#DNA
D080 <- as.dist(allE0F80, diag = T) #käytetään as.dist, koska nämä ovat jo etäisyyksiä (olisi dist jos etäisyyksiä ei olisi laskettu)
D180 <- as.dist(allE1F80, diag = T)
D280 <- as.dist(allE2F80, diag = T)
D380 <- as.dist(allE3F80, diag = T)
D480 <- as.dist(allE4F80, diag = T)

ajo0_80 <- hclustlaskentaF(D080,allE0Fmeta, k = 142 )
#ajo0_80
ajo1_80 <- hclustlaskentaF(D180,allE1Fmeta, k = 143 )
ajo2_80 <-hclustlaskentaF(D280,allE2Fmeta, k = 144 )
ajo3_80 <-hclustlaskentaF(D380,allE3Fmeta, k = 147 )
ajo4_80 <-hclustlaskentaF(D480,allE4Fmeta, k = 142 )

#PIIRTEET
allE0F <- as.dist(allE0F)
allE1F <- as.dist(allE1F)
allE2F <- as.dist(allE2F)
allE3F <- as.dist(allE3F)
allE4F <- as.dist(allE4F)

ajo0f <- hclustlaskentaF(allE0F,allE0Fmeta, k = 142 )
ajo1f <- hclustlaskentaF(allE1F,allE1Fmeta, k = 143 )
ajo2f <- hclustlaskentaF(allE2F,allE2Fmeta, k = 144 )
ajo3f <- hclustlaskentaF(allE3F,allE3Fmeta, k = 147 )
ajo4f <- hclustlaskentaF(allE4F,allE4Fmeta, k = 142 )

#tulosten autoluku uusifu:ssa

##########
D080OR <- as.dist(allE0OR80)
D180OR <- as.dist(allE1OR80)
D280OR <- as.dist(allE2OR80)
D380OR <- as.dist(allE3OR80)
D480OR <- as.dist(allE4OR80)

ajo080or <- hclustlaskentaOR(D080OR,allE0ORmeta, k = 19 ) #19 oikea, voi kokeilla tuplaa
ajo180or <- hclustlaskentaOR(D180OR,allE1ORmeta, k = 19 )
ajo280or <- hclustlaskentaOR(D280OR,allE2ORmeta, k = 19 )
ajo380or <- hclustlaskentaOR(D380OR,allE3ORmeta, k = 19 )
ajo480or <- hclustlaskentaOR(D480OR,allE4ORmeta, k = 19 )

#PIIRTEET
allE0OR <- as.dist(allE0OR)
allE1OR <- as.dist(allE1OR)
allE2OR <- as.dist(allE2OR)
allE3OR <- as.dist(allE3OR)
allE4OR <- as.dist(allE4OR)

ajo0OR <- hclustlaskentaOR(allE0OR,allE0ORmeta, k = 19 )
ajo1OR <- hclustlaskentaOR(allE1OR,allE1ORmeta, k = 19 )
ajo2OR <- hclustlaskentaOR(allE2OR,allE2ORmeta, k = 19 )
ajo3OR <- hclustlaskentaOR(allE3OR,allE3ORmeta, k = 19 )
ajo4OR <- hclustlaskentaOR(allE4OR,allE4ORmeta, k = 19 )

#tulosten autoluku uusifu:ssa


##indeksit

#family
indeksilaskuri(ajo0_80$a, ajo0f$a);print("######### B ALKAA #########");indeksilaskuri(ajo0_80$b, ajo0f$b);print("######### C ALKAA #########");indeksilaskuri(ajo0_80$c, ajo0f$c)
indeksilaskuri(ajo1_80$a, ajo1f$a);print("######### B ALKAA #########");indeksilaskuri(ajo1_80$b, ajo1f$b);print("######### C ALKAA #########");indeksilaskuri(ajo1_80$c, ajo1f$c)
indeksilaskuri(ajo2_80$a, ajo2f$a);print("######### B ALKAA #########");indeksilaskuri(ajo2_80$b, ajo2f$b);print("######### C ALKAA #########");indeksilaskuri(ajo2_80$c, ajo2f$c)
indeksilaskuri(ajo3_80$a, ajo3f$a);print("######### B ALKAA #########");indeksilaskuri(ajo3_80$b, ajo3f$b);print("######### C ALKAA #########");indeksilaskuri(ajo3_80$c, ajo3f$c)
indeksilaskuri(ajo4_80$a, ajo4f$a);print("######### B ALKAA #########");indeksilaskuri(ajo4_80$b, ajo4f$b);print("######### C ALKAA #########");indeksilaskuri(ajo4_80$c, ajo4f$c)


#indeksien kootut tulokset indeksien_valit_laskennass

#order
indeksilaskuri(ajo080or$a, ajo0OR$a);print("######### B ALKAA #########");indeksilaskuri(ajo080or$b, ajo0OR$b);print("######### C ALKAA #########");indeksilaskuri(ajo080or$c, ajo0OR$c)
indeksilaskuri(ajo180or$a, ajo1OR$a);print("######### B ALKAA #########");indeksilaskuri(ajo180or$b, ajo1OR$b);print("######### C ALKAA #########");indeksilaskuri(ajo180or$c, ajo1OR$c)
indeksilaskuri(ajo280or$a, ajo2OR$a);print("######### B ALKAA #########");indeksilaskuri(ajo280or$b, ajo2OR$b);print("######### C ALKAA #########");indeksilaskuri(ajo280or$c, ajo2OR$c)
indeksilaskuri(ajo380or$a, ajo3OR$a);print("######### B ALKAA #########");indeksilaskuri(ajo380or$b, ajo3OR$b);print("######### C ALKAA #########");indeksilaskuri(ajo380or$c, ajo3OR$c)
indeksilaskuri(ajo480or$a, ajo4OR$a);print("######### B ALKAA #########");indeksilaskuri(ajo480or$b, ajo4OR$b);print("######### C ALKAA #########");indeksilaskuri(ajo480or$c, ajo4OR$c)

########## kuvia
plot_unweighted_noise(ajo0_80$D,"Coccinellidae")
plot_weighted_noise(ajo0f$D,"Coccinellidae")

