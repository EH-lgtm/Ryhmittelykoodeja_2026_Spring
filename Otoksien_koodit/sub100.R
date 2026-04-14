#gradu koodi SUB100
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


#####ANALYYSI ALKAA#####
#######################
#KMEDOIDS: PAM

#dna family
PAMajo0f80 <- pamLaskentaF(sub100F80,sub100Fmeta, k = 19)
PAMajo0f80$a
#piirteet family
PAMajo0f <- pamLaskentaF(sub100F,sub100Fmeta, k = 19)

#indeksit family
indeksilaskuri(PAMajo0f80$a, PAMajo0f$a);print("######### B ALKAA #########"); indeksilaskuri(PAMajo0f80$b, PAMajo0f$b)

## automaattinen tulosten lukija uusifu:ssa, sample0

##########
#dna order
PAMajo0OR80 <- pamLaskentaOR(sub100_80,metasub100, k = 8) #8,16
#piirteet
PAMajo0OR <- pamLaskentaOR(sub100,metasub100, k = 8) #8,16

#indeksit ###########################################################
#order
indeksilaskuri(PAMajo0OR80$a, PAMajo0OR$a);print("######### B ALKAA #########");indeksilaskuri(PAMajo0OR80$b, PAMajo0OR$b)
## automaattinen tulosten lukija uusifu:ssa, sample0



#hierarkkinen jakava: DIANA

#dna FAMILY
Dia0F80 <- dianalaskentaF(sub100F80,sub100Fmeta, k = 19)
#piirteet family
Dia0F <- dianalaskentaF(sub100F,sub100Fmeta, k = 19)
#indeksit
indeksilaskuri(Dia0F80$a,Dia0F$a);print("######### B ALKAA #########");indeksilaskuri(Dia0F80$b,Dia0F$b)
## automaattinen tulosten lukija uusifu:ssa, sample0


#ORDER
#dna order
Dia0OR80 <- dianalaskentaOR(sub100_80,metasub100, k = 8) #8, 16
table(metasub100$order)
#piirteet order
Dia0OR <- dianalaskentaOR(sub100,metasub100, k = 8) # 8, 16
#indeksit order
indeksilaskuri(Dia0OR80$a,Dia0OR$a);print("######### B ALKAA #########");indeksilaskuri(Dia0OR80$b,Dia0OR$b)
## automaattinen tulosten lukija uusifu:ssa, sample0


#hierarkkinen yhdistävä: HCLUST
options(max.print=2500)
#DNA
D080 <- as.dist(sub100F80)
#piirteet
sub100F <- as.dist(sub100F)

ajo0_80 <- hclustlaskentaF(D080,sub100Fmeta, k = 19 ) #FAM
ajo0f <- hclustlaskentaF(sub100F,sub100Fmeta, k = 19 )

## automaattinen tulosten lukija uusifu:ssa, sample0



#family
indeksilaskuri(ajo0_80$a, ajo0f$a);print("######### B ALKAA #########");indeksilaskuri(ajo0_80$b, ajo0f$b);print("######### C ALKAA #########");indeksilaskuri(ajo0_80$c, ajo0f$c)

########
#ORDER
D080OR <- as.dist(sub100_80)
ajo080or <- hclustlaskentaOR(D080OR,metasub100, k = 8 ) #ORDER k= 8, 16
sub100OR <- as.dist(sub100)
ajo0OR <- hclustlaskentaOR(sub100OR,metasub100, k = 8 ) #k= 8, 16

## automaattinen tulosten lukija uusifu:ssa, sample0
indeksilaskuri(ajo080or$a, ajo0OR$a);print("######### B ALKAA #########");indeksilaskuri(ajo080or$b, ajo0OR$b);print("######### C ALKAA #########");indeksilaskuri(ajo080or$c, ajo0OR$c)
############
