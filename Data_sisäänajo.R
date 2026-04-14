################################################################################
###############################

#DATAN SISÄÄNAJO

###############################
################################################################################

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



################################################################################

##All Even

################################################################################

#piirteet
allE0 <- read.csv("bioscan_bioclip_all_even_sample0_distances.csv")
allE1 <- read.csv("bioscan_bioclip_all_even_sample1_distances.csv")
allE2 <- read.csv("bioscan_bioclip_all_even_sample2_distances.csv")
allE3 <- read.csv("bioscan_bioclip_all_even_sample3_distances.csv")
allE4 <- read.csv("bioscan_bioclip_all_even_sample4_distances.csv")


#dna
allE080 <- read.csv("bioscan_all_even_sample0_k80dist.csv")
allE180 <- read.csv("bioscan_all_even_sample1_k80dist.csv")
allE280 <- read.csv("bioscan_all_even_sample2_k80dist.csv")
allE380 <- read.csv("bioscan_all_even_sample3_k80dist.csv")
allE480 <- read.csv("bioscan_all_even_sample4_k80dist.csv")


#meta
metaallE0<-read_parquet("bioscan_metadata_all_even_sample0.parquet")
metaallE1<-read_parquet("bioscan_metadata_all_even_sample1.parquet")
metaallE2<-read_parquet("bioscan_metadata_all_even_sample2.parquet")
metaallE3<-read_parquet("bioscan_metadata_all_even_sample3.parquet")
metaallE4<-read_parquet("bioscan_metadata_all_even_sample4.parquet")

#tarkastetaan taksonit
sum((metaallE0$subspecies)=="not_classified") #ei labeleitä, 1000 NC
sum((metaallE0$species)=="not_classified") #NC 859
sum((metaallE0$genus )=="not_classified") #NC 702
sum((metaallE0$tribe)=="not_classified") #NC 962
sum((metaallE0$subfamily)=="not_classified") #NC 694
sum((metaallE0$family)=="not_classified") #NC 44
sum((metaallE0$order)=="not_classified") #NC 0 = kaikki labelit paikassaan
#tarkastetaan family, order

length(unique(metaallE0$family)) 
length(unique(metaallE0$order))
table(metaallE4$family) #142, 143, 144, 147, 142 eri tyyppiä samplessa 0-4 UH OH Näitä on liikaa
table(metaallE0$order) #19,19,19,19,19 (kuten kuuluu)

#poistetaan NC, family

allE0F <- allE0[,2:1001] #otetaan nimi kolumni pois sotkemasta
allE1F <- allE1[,2:1001]
allE2F <- allE2[,2:1001]
allE3F <- allE3[,2:1001]
allE4F <- allE4[,2:1001]

#piirteet
allE0F <- allE0F[metaallE0$family != "not_classified",metaallE0$family != "not_classified"]
allE1F <- allE1F[metaallE1$family != "not_classified",metaallE1$family != "not_classified"]
allE2F <- allE2F[metaallE2$family != "not_classified",metaallE2$family != "not_classified"]
allE3F <- allE3F[metaallE3$family != "not_classified",metaallE3$family != "not_classified"]
allE4F <- allE4F[metaallE4$family != "not_classified",metaallE4$family != "not_classified"]
#dna
allE0F80 <- allE080[metaallE0$family != "not_classified",metaallE0$family != "not_classified"]
allE1F80 <- allE180[metaallE1$family != "not_classified",metaallE1$family != "not_classified"]
allE2F80 <- allE280[metaallE2$family != "not_classified",metaallE2$family != "not_classified"]
allE3F80 <- allE380[metaallE3$family != "not_classified",metaallE3$family != "not_classified"]
allE4F80 <- allE480[metaallE4$family != "not_classified",metaallE4$family != "not_classified"]

dim(allE0F)
dim(allE1F)
dim(allE2F)
dim(allE3F)
dim(allE4F)

dim(allE0F80);dim(allE1F80);dim(allE2F80);dim(allE3F80);dim(allE4F80)

#meta
allE0Fmeta <- metaallE0[metaallE0$family != "not_classified",]
allE1Fmeta <- metaallE1[metaallE1$family != "not_classified",]
allE2Fmeta <- metaallE2[metaallE2$family != "not_classified",]
allE3Fmeta <- metaallE3[metaallE3$family != "not_classified",]
allE4Fmeta <- metaallE4[metaallE4$family != "not_classified",]
sum(allE0Fmeta$family =="not_classified")
dim(allE0Fmeta)

# ORDER ei tarvitse poistoa (nimetään vain uudelleen)

#piirteet
allE0OR <- allE0[,2:1001] #otetaan nimi kolumni pois sotkemasta
allE1OR <- allE1[,2:1001]
allE2OR <- allE2[,2:1001]
allE3OR <- allE3[,2:1001]
allE4OR <- allE4[,2:1001]
#dna
allE0OR80 <- allE080
allE1OR80 <- allE180
allE2OR80 <- allE280
allE3OR80 <- allE380
allE4OR80 <- allE480

allE0ORmeta <- metaallE0
allE1ORmeta <- metaallE1
allE2ORmeta <- metaallE2
allE3ORmeta <- metaallE3
allE4ORmeta <- metaallE4

sum(allE0ORmeta$order =="not_classified")
dim(allE0ORmeta)

dim(allE0OR)
dim(allE0OR80)

################################################################################

##Coleoptera (family, olemassa sf, mutta ei oteta mukaan nyt)

################################################################################

#piirteet
cole0 <- read.csv("bioscan_bioclip_coleoptera_sample0_distances.csv")
cole1 <- read.csv("bioscan_bioclip_coleoptera_sample1_distances.csv")
cole2 <- read.csv("bioscan_bioclip_coleoptera_sample2_distances.csv")
cole3 <- read.csv("bioscan_bioclip_coleoptera_sample3_distances.csv")
cole4 <- read.csv("bioscan_bioclip_coleoptera_sample4_distances.csv")


#dna
cole080 <- read.csv("bioscan_coleoptera_sample0_k80dist.csv")
cole180 <- read.csv("bioscan_coleoptera_sample1_k80dist.csv")
cole280 <- read.csv("bioscan_coleoptera_sample2_k80dist.csv")
cole380 <- read.csv("bioscan_coleoptera_sample3_k80dist.csv")
cole480 <- read.csv("bioscan_coleoptera_sample4_k80dist.csv")


#meta
metacole0<-read_parquet("bioscan_metadata_coleoptera_sample0.parquet")
metacole1<-read_parquet("bioscan_metadata_coleoptera_sample1.parquet")
metacole2<-read_parquet("bioscan_metadata_coleoptera_sample2.parquet")
metacole3<-read_parquet("bioscan_metadata_coleoptera_sample3.parquet")
metacole4<-read_parquet("bioscan_metadata_coleoptera_sample4.parquet")

#tarkastetaan taksonit
sum((metacole0$subspecies)=="not_classified") #ei labeleitä, 1000 NC
sum((metacole0$species)=="not_classified") #NC 935
sum((metacole0$genus )=="not_classified") #NC 889
sum((metacole0$tribe)=="not_classified") #NC 903
sum((metacole0$subfamily)=="not_classified") #NC 593
sum((metacole0$family)=="not_classified") #NC 1
sum((metacole0$order)=="not_classified") #NC 0 = kaikki labelit paikassaan
#tarkastetaan family, subfamily

length(unique(metacole4$family)) 
length(unique(metacole4$subfamily))
table(metacole4$family) #56, 54, 58, 61, 59 eri tyyppiä samplessa 0-4
table(metacole0$subfamily) #49, 51, 51, 47, 52


#poistetaan NC, family

cole0F <- cole0[,2:1001] #otetaan nimi kolumni pois sotkemasta
cole1F <- cole1[,2:1001]
cole2F <- cole2[,2:1001]
cole3F <- cole3[,2:1001]
cole4F <- cole4[,2:1001]

#piirteet
cole0F <- cole0F[metacole0$family != "not_classified",metacole0$family != "not_classified"]
cole1F <- cole1F[metacole1$family != "not_classified",metacole1$family != "not_classified"]
cole2F <- cole2F[metacole2$family != "not_classified",metacole2$family != "not_classified"]
cole3F <- cole3F[metacole3$family != "not_classified",metacole3$family != "not_classified"]
cole4F <- cole4F[metacole4$family != "not_classified",metacole4$family != "not_classified"]
#dna
cole0F80 <- cole080[metacole0$family != "not_classified",metacole0$family != "not_classified"]
cole1F80 <- cole180[metacole1$family != "not_classified",metacole1$family != "not_classified"]
cole2F80 <- cole280[metacole2$family != "not_classified",metacole2$family != "not_classified"]
cole3F80 <- cole380[metacole3$family != "not_classified",metacole3$family != "not_classified"]
cole4F80 <- cole480[metacole4$family != "not_classified",metacole4$family != "not_classified"]

dim(cole0F)
dim(cole1F)
dim(cole2F)
dim(cole3F)
dim(cole4F)

dim(cole0F80);dim(cole1F80);dim(cole2F80);dim(cole3F80);dim(cole4F80)

#meta
cole0Fmeta <- metacole0[metacole0$family != "not_classified",]
cole1Fmeta <- metacole1[metacole1$family != "not_classified",]
cole2Fmeta <- metacole2[metacole2$family != "not_classified",]
cole3Fmeta <- metacole3[metacole3$family != "not_classified",]
cole4Fmeta <- metacole4[metacole4$family != "not_classified",]
sum(cole0Fmeta$family =="not_classified")
dim(cole0Fmeta)

################################################################################

##Diptera (family)

################################################################################

#piirteet
dip0 <- read.csv("bioscan_bioclip_diptera_sample0_distances.csv")
dip1 <- read.csv("bioscan_bioclip_diptera_sample1_distances.csv")
dip2 <- read.csv("bioscan_bioclip_diptera_sample2_distances.csv")
dip3 <- read.csv("bioscan_bioclip_diptera_sample3_distances.csv")
dip4 <- read.csv("bioscan_bioclip_diptera_sample4_distances.csv")


#dna
dip080 <- read.csv("bioscan_diptera_sample0_k80dist.csv")
dip180 <- read.csv("bioscan_diptera_sample1_k80dist.csv")
dip280 <- read.csv("bioscan_diptera_sample2_k80dist.csv")
dip380 <- read.csv("bioscan_diptera_sample3_k80dist.csv")
dip480 <- read.csv("bioscan_diptera_sample4_k80dist.csv")


#meta
metadip0<-read_parquet("bioscan_metadata_diptera_sample0.parquet")
metadip1<-read_parquet("bioscan_metadata_diptera_sample1.parquet")
metadip2<-read_parquet("bioscan_metadata_diptera_sample2.parquet")
metadip3<-read_parquet("bioscan_metadata_diptera_sample3.parquet")
metadip4<-read_parquet("bioscan_metadata_diptera_sample4.parquet")

#tarkastetaan taksonit
sum((metadip0$subspecies)=="not_classified") #ei labeleitä, 1000 NC
sum((metadip0$species)=="not_classified") #NC 955
sum((metadip0$genus )=="not_classified") #NC 800
sum((metadip0$tribe)=="not_classified") #NC 962
sum((metadip0$subfamily)=="not_classified") #NC 803
sum((metadip0$family)=="not_classified") #NC 1
sum((metadip0$order)=="not_classified") #NC 0 = kaikki labelit paikassaan
#tarkastetaan vain family

length(unique(metadip0$family)) 
length(unique(metadip0$subfamily))
table(metadip4$family) #32, 31, 30, 33, 32 eri tyyppiä samplessa 0-4
#table(metadip0$subfamily) #31, 34, 35, 31, 34


#poistetaan NC, family

dip0F <- dip0[,2:1001] #otetaan nimi kolumni pois sotkemasta
dip1F <- dip1[,2:1001]
dip2F <- dip2[,2:1001]
dip3F <- dip3[,2:1001]
dip4F <- dip4[,2:1001]

#piirteet
dip0F <- dip0F[metadip0$family != "not_classified",metadip0$family != "not_classified"]
dip1F <- dip1F[metadip1$family != "not_classified",metadip1$family != "not_classified"]
dip2F <- dip2F[metadip2$family != "not_classified",metadip2$family != "not_classified"]
dip3F <- dip3F[metadip3$family != "not_classified",metadip3$family != "not_classified"]
dip4F <- dip4F[metadip4$family != "not_classified",metadip4$family != "not_classified"]
#dna
dip0F80 <- dip080[metadip0$family != "not_classified",metadip0$family != "not_classified"]
dip1F80 <- dip180[metadip1$family != "not_classified",metadip1$family != "not_classified"]
dip2F80 <- dip280[metadip2$family != "not_classified",metadip2$family != "not_classified"]
dip3F80 <- dip380[metadip3$family != "not_classified",metadip3$family != "not_classified"]
dip4F80 <- dip480[metadip4$family != "not_classified",metadip4$family != "not_classified"]

dim(dip0F)
dim(dip1F)
dim(dip2F)
dim(dip3F)
dim(dip4F)
dim(dip0F80)

#meta
dip0Fmeta <- metadip0[metadip0$family != "not_classified",]
dip1Fmeta <- metadip1[metadip1$family != "not_classified",]
dip2Fmeta <- metadip2[metadip2$family != "not_classified",]
dip3Fmeta <- metadip3[metadip3$family != "not_classified",]
dip4Fmeta <- metadip4[metadip4$family != "not_classified",]
sum(dip0Fmeta$family =="not_classified")
dim(dip0Fmeta)

################################################################################

##Hemiptera (family)

################################################################################

#piirteet
hemi0 <- read.csv("bioscan_bioclip_hemiptera_sample0_distances.csv")
hemi1 <- read.csv("bioscan_bioclip_hemiptera_sample1_distances.csv")
hemi2 <- read.csv("bioscan_bioclip_hemiptera_sample2_distances.csv")
hemi3 <- read.csv("bioscan_bioclip_hemiptera_sample3_distances.csv")
hemi4 <- read.csv("bioscan_bioclip_hemiptera_sample4_distances.csv")


#dna
hemi080 <- read.csv("bioscan_hemiptera_sample0_k80dist.csv")
hemi180 <- read.csv("bioscan_hemiptera_sample1_k80dist.csv")
hemi280 <- read.csv("bioscan_hemiptera_sample2_k80dist.csv")
hemi380 <- read.csv("bioscan_hemiptera_sample3_k80dist.csv")
hemi480 <- read.csv("bioscan_hemiptera_sample4_k80dist.csv")


#meta
metahemi0<-read_parquet("bioscan_metadata_hemiptera_sample0.parquet")
metahemi1<-read_parquet("bioscan_metadata_hemiptera_sample1.parquet")
metahemi2<-read_parquet("bioscan_metadata_hemiptera_sample2.parquet")
metahemi3<-read_parquet("bioscan_metadata_hemiptera_sample3.parquet")
metahemi4<-read_parquet("bioscan_metadata_hemiptera_sample4.parquet")

#tarkastetaan taksonit
sum((metahemi0$subspecies)=="not_classified") #ei labeleitä, 1000 NC
sum((metahemi0$species)=="not_classified") #NC 935
sum((metahemi0$genus )=="not_classified") #NC 828
sum((metahemi0$tribe)=="not_classified") #NC 891
sum((metahemi0$subfamily)=="not_classified") #NC 762 HMMMM >:I
sum((metahemi0$family)=="not_classified") #NC 72
sum((metahemi0$order)=="not_classified") #NC 0 = kaikki labelit paikassaan
#tarkastellaan VAIN family:ä

#36, 38, 36, 40, 33 eri tyyppiä sampleissa 0-4
table(metahemi0$family) #ikävästi taas yksittäisiä


#poistetaan NC, family

hemi0F <- hemi0[,2:1001] #otetaan nimi kolumni pois sotkemasta
hemi1F <- hemi1[,2:1001]
hemi2F <- hemi2[,2:1001]
hemi3F <- hemi3[,2:1001]
hemi4F <- hemi4[,2:1001]

#piirteet
hemi0F <- hemi0F[metahemi0$family != "not_classified",metahemi0$family != "not_classified"]
hemi1F <- hemi1F[metahemi1$family != "not_classified",metahemi1$family != "not_classified"]
hemi2F <- hemi2F[metahemi2$family != "not_classified",metahemi2$family != "not_classified"]
hemi3F <- hemi3F[metahemi3$family != "not_classified",metahemi3$family != "not_classified"]
hemi4F <- hemi4F[metahemi4$family != "not_classified",metahemi4$family != "not_classified"]
#dna
hemi0F80 <- hemi080[metahemi0$family != "not_classified",metahemi0$family != "not_classified"]
hemi1F80 <- hemi180[metahemi1$family != "not_classified",metahemi1$family != "not_classified"]
hemi2F80 <- hemi280[metahemi2$family != "not_classified",metahemi2$family != "not_classified"]
hemi3F80 <- hemi380[metahemi3$family != "not_classified",metahemi3$family != "not_classified"]
hemi4F80 <- hemi480[metahemi4$family != "not_classified",metahemi4$family != "not_classified"]

dim(hemi0F)
dim(hemi1F)
dim(hemi2F)
dim(hemi3F)
dim(hemi4F)
dim(hemi0F80)

#meta
hemi0Fmeta <- metahemi0[metahemi0$family != "not_classified",]
hemi1Fmeta <- metahemi1[metahemi1$family != "not_classified",]
hemi2Fmeta <- metahemi2[metahemi2$family != "not_classified",]
hemi3Fmeta <- metahemi3[metahemi3$family != "not_classified",]
hemi4Fmeta <- metahemi4[metahemi4$family != "not_classified",]
sum(hemi0Fmeta$family =="not_classified")
dim(hemi0Fmeta)


################################################################################

##Hymenoptera (family - olemassa subf, mutta ei oteta tänne mukaan)

################################################################################

#piirteet
hyme0 <- read.csv("bioscan_bioclip_hymenoptera_sample0_distances.csv")
hyme1 <- read.csv("bioscan_bioclip_hymenoptera_sample1_distances.csv")
hyme2 <- read.csv("bioscan_bioclip_hymenoptera_sample2_distances.csv")
hyme3 <- read.csv("bioscan_bioclip_hymenoptera_sample3_distances.csv")
hyme4 <- read.csv("bioscan_bioclip_hymenoptera_sample4_distances.csv")


#dna
hyme080 <- read.csv("bioscan_hymenoptera_sample0_k80dist.csv")
hyme180 <- read.csv("bioscan_hymenoptera_sample1_k80dist.csv")
hyme280 <- read.csv("bioscan_hymenoptera_sample2_k80dist.csv")
hyme380 <- read.csv("bioscan_hymenoptera_sample3_k80dist.csv")
hyme480 <- read.csv("bioscan_hymenoptera_sample4_k80dist.csv")


#meta
metahyme0<-read_parquet("bioscan_metadata_hymenoptera_sample0.parquet")
metahyme1<-read_parquet("bioscan_metadata_hymenoptera_sample1.parquet")
metahyme2<-read_parquet("bioscan_metadata_hymenoptera_sample2.parquet")
metahyme3<-read_parquet("bioscan_metadata_hymenoptera_sample3.parquet")
metahyme4<-read_parquet("bioscan_metadata_hymenoptera_sample4.parquet")

#tarkastetaan taksonit
sum((metahyme0$subspecies)=="not_classified") #ei labeleitä, 1000 NC
sum((metahyme0$species)=="not_classified") #NC 795
sum((metahyme0$genus )=="not_classified") #NC 651
sum((metahyme0$tribe)=="not_classified") #NC 914
sum((metahyme0$subfamily)=="not_classified") #NC 621
sum((metahyme0$family)=="not_classified") #NC 5
sum((metahyme0$order)=="not_classified") #NC 0 = kaikki labelit paikassaan
#tarkastetaan family, subfamily

length(unique(metahyme0$family)) #33, 37, 38, 38, 34 eri tyyppiä samplessa 0-4
length(unique(metahyme0$subfamily))  #58, 56, 61, 61, 58

table(metahyme4$family) 
table(metahyme0$subfamily) #paljon paljon ykkösiä


#poistetaan NC, family

hyme0F <- hyme0[,2:1001] #otetaan nimi kolumni pois sotkemasta
hyme1F <- hyme1[,2:1001]
hyme2F <- hyme2[,2:1001]
hyme3F <- hyme3[,2:1001]
hyme4F <- hyme4[,2:1001]

#piirteet
hyme0F <- hyme0F[metahyme0$family != "not_classified",metahyme0$family != "not_classified"]
hyme1F <- hyme1F[metahyme1$family != "not_classified",metahyme1$family != "not_classified"]
hyme2F <- hyme2F[metahyme2$family != "not_classified",metahyme2$family != "not_classified"]
hyme3F <- hyme3F[metahyme3$family != "not_classified",metahyme3$family != "not_classified"]
hyme4F <- hyme4F[metahyme4$family != "not_classified",metahyme4$family != "not_classified"]
#dna
hyme0F80 <- hyme080[metahyme0$family != "not_classified",metahyme0$family != "not_classified"]
hyme1F80 <- hyme180[metahyme1$family != "not_classified",metahyme1$family != "not_classified"]
hyme2F80 <- hyme280[metahyme2$family != "not_classified",metahyme2$family != "not_classified"]
hyme3F80 <- hyme380[metahyme3$family != "not_classified",metahyme3$family != "not_classified"]
hyme4F80 <- hyme480[metahyme4$family != "not_classified",metahyme4$family != "not_classified"]

dim(hyme0F)
dim(hyme1F)
dim(hyme2F)
dim(hyme3F)
dim(hyme4F)
dim(hyme0F80)

#meta
hyme0Fmeta <- metahyme0[metahyme0$family != "not_classified",]
hyme1Fmeta <- metahyme1[metahyme1$family != "not_classified",]
hyme2Fmeta <- metahyme2[metahyme2$family != "not_classified",]
hyme3Fmeta <- metahyme3[metahyme3$family != "not_classified",]
hyme4Fmeta <- metahyme4[metahyme4$family != "not_classified",]
sum(hyme0Fmeta$family =="not_classified")
dim(hyme0Fmeta)


################################################################################

##Lepidoptera (family, sf olemassa)

################################################################################

#piirteet

lepi0 <- read.csv("bioscan_bioclip_lepidoptera_sample0_distances.csv")
lepi1 <- read.csv("bioscan_bioclip_lepidoptera_sample1_distances.csv")
lepi2 <- read.csv("bioscan_bioclip_lepidoptera_sample2_distances.csv")
lepi3 <- read.csv("bioscan_bioclip_lepidoptera_sample3_distances.csv")
lepi4 <- read.csv("bioscan_bioclip_lepidoptera_sample4_distances.csv")


#dna
lepi080 <- read.csv("bioscan_lepidoptera_sample0_k80dist.csv")
lepi180 <- read.csv("bioscan_lepidoptera_sample1_k80dist.csv")
lepi280 <- read.csv("bioscan_lepidoptera_sample2_k80dist.csv")
lepi380 <- read.csv("bioscan_lepidoptera_sample3_k80dist.csv")
lepi480 <- read.csv("bioscan_lepidoptera_sample4_k80dist.csv")


#meta
metalepi0<-read_parquet("bioscan_metadata_lepidoptera_sample0.parquet")
metalepi1<-read_parquet("bioscan_metadata_lepidoptera_sample1.parquet")
metalepi2<-read_parquet("bioscan_metadata_lepidoptera_sample2.parquet")
metalepi3<-read_parquet("bioscan_metadata_lepidoptera_sample3.parquet")
metalepi4<-read_parquet("bioscan_metadata_lepidoptera_sample4.parquet")

#tarkastetaan taksonit
sum((metalepi0$subspecies)=="not_classified") #ei labeleitä, 1000 NC
sum((metalepi0$species)=="not_classified") #NC 597
sum((metalepi0$genus )=="not_classified") #NC 530
sum((metalepi0$tribe)=="not_classified") #NC 969
sum((metalepi0$subfamily)=="not_classified") #NC 691 semi vähän jäljellä, mutta katotaan silti
sum((metalepi0$family)=="not_classified") #NC 233
sum((metalepi0$order)=="not_classified") #NC 0 = kaikki labelit paikassaan
#tarkastetaan family, subfamily


length(unique(metalepi4$species)) #209, 204, 193, 197, 198  NOPE
length(unique(metalepi4$genus)) #134, 131, 140, 136, 130 NOPE
length(unique(metalepi4$subfamily)) #46, 51, 42, 47, 48
length(unique(metalepi4$family)) #43, 40, 46, 37, 41

############

#poistetaan NC, family

lepi0F <- lepi0[,2:1001] #otetaan nimi kolumni pois sotkemasta
lepi1F <- lepi1[,2:1001]
lepi2F <- lepi2[,2:1001]
lepi3F <- lepi3[,2:1001]
lepi4F <- lepi4[,2:1001]

#piirteet
lepi0F <- lepi0F[metalepi0$family != "not_classified",metalepi0$family != "not_classified"]
lepi1F <- lepi1F[metalepi1$family != "not_classified",metalepi1$family != "not_classified"]
lepi2F <- lepi2F[metalepi2$family != "not_classified",metalepi2$family != "not_classified"]
lepi3F <- lepi3F[metalepi3$family != "not_classified",metalepi3$family != "not_classified"]
lepi4F <- lepi4F[metalepi4$family != "not_classified",metalepi4$family != "not_classified"]
#dna
lepi0F80 <- lepi080[metalepi0$family != "not_classified",metalepi0$family != "not_classified"]
lepi1F80 <- lepi180[metalepi1$family != "not_classified",metalepi1$family != "not_classified"]
lepi2F80 <- lepi280[metalepi2$family != "not_classified",metalepi2$family != "not_classified"]
lepi3F80 <- lepi380[metalepi3$family != "not_classified",metalepi3$family != "not_classified"]
lepi4F80 <- lepi480[metalepi4$family != "not_classified",metalepi4$family != "not_classified"]

dim(lepi0F)
dim(lepi1F)
dim(lepi2F)
dim(lepi3F)
dim(lepi4F)
dim(lepi0F80)

#meta
lepi0Fmeta <- metalepi0[metalepi0$family != "not_classified",]
lepi1Fmeta <- metalepi1[metalepi1$family != "not_classified",]
lepi2Fmeta <- metalepi2[metalepi2$family != "not_classified",]
lepi3Fmeta <- metalepi3[metalepi3$family != "not_classified",]
lepi4Fmeta <- metalepi4[metalepi4$family != "not_classified",]
sum(lepi0Fmeta$family =="not_classified")
dim(lepi0Fmeta)

#k = 43, 40, 46, 37, 41 kun pienet mukana FAMILY


################################################################################

##LargeSize, fam, order (kerro mitkä fam kun kirjotat)

################################################################################


#piirteet
large0 <- read.csv("bioscan_bioclip_largesize_sample0_distances.csv")
large1 <- read.csv("bioscan_bioclip_largesize_sample1_distances.csv")
large2 <- read.csv("bioscan_bioclip_largesize_sample2_distances.csv")
large3 <- read.csv("bioscan_bioclip_largesize_sample3_distances.csv")
large4 <- read.csv("bioscan_bioclip_largesize_sample4_distances.csv")


#dna
large080 <- read.csv("bioscan_largesize_sample0_k80dist.csv")
large180 <- read.csv("bioscan_largesize_sample1_k80dist.csv")
large280 <- read.csv("bioscan_largesize_sample2_k80dist.csv")
large380 <- read.csv("bioscan_largesize_sample3_k80dist.csv")
large480 <- read.csv("bioscan_largesize_sample4_k80dist.csv")


#meta
metalarge0<-read_parquet("bioscan_metadata_largesize_sample0.parquet")
metalarge1<-read_parquet("bioscan_metadata_largesize_sample1.parquet")
metalarge2<-read_parquet("bioscan_metadata_largesize_sample2.parquet")
metalarge3<-read_parquet("bioscan_metadata_largesize_sample3.parquet")
metalarge4<-read_parquet("bioscan_metadata_largesize_sample4.parquet")

#tarkastetaan taksonit
sum((metalarge0$subspecies)=="not_classified") #ei labeleitä, 1000
sum((metalarge0$species)=="not_classified") #NC 855
sum((metalarge0$genus )=="not_classified") #NC 741
sum((metalarge0$tribe)=="not_classified") #NC 919
sum((metalarge0$subfamily)=="not_classified") #NC 722
sum((metalarge0$family)=="not_classified") #NC 80
sum((metalarge0$order)=="not_classified") #NC 0 = kaikki labelit paikassaan
#tarkastetaan order, family

length(unique(metalarge4$family)) 
length(unique(metalarge0$subfamily))
table(metalarge4$family) #119, 119, 117, 126, 117 eri tyyppiä samplessa 0-4
table(metalarge0$order) #5,5,5,5,5 (valitse joku hyvä k)



#poistetaan NC, family

large0F <- large0[,2:1001] #otetaan nimi kolumni pois sotkemasta
large1F <- large1[,2:1001]
large2F <- large2[,2:1001]
large3F <- large3[,2:1001]
large4F <- large4[,2:1001]

#piirteet
large0F <- large0F[metalarge0$family != "not_classified",metalarge0$family != "not_classified"]
large1F <- large1F[metalarge1$family != "not_classified",metalarge1$family != "not_classified"]
large2F <- large2F[metalarge2$family != "not_classified",metalarge2$family != "not_classified"]
large3F <- large3F[metalarge3$family != "not_classified",metalarge3$family != "not_classified"]
large4F <- large4F[metalarge4$family != "not_classified",metalarge4$family != "not_classified"]
#dna
large0F80 <- large080[metalarge0$family != "not_classified",metalarge0$family != "not_classified"]
large1F80 <- large180[metalarge1$family != "not_classified",metalarge1$family != "not_classified"]
large2F80 <- large280[metalarge2$family != "not_classified",metalarge2$family != "not_classified"]
large3F80 <- large380[metalarge3$family != "not_classified",metalarge3$family != "not_classified"]
large4F80 <- large480[metalarge4$family != "not_classified",metalarge4$family != "not_classified"]

dim(large0F)
dim(large1F)
dim(large2F)
dim(large3F)
dim(large4F)
dim(large0F80)


#meta
large0Fmeta <- metalarge0[metalarge0$family != "not_classified",]
large1Fmeta <- metalarge1[metalarge1$family != "not_classified",]
large2Fmeta <- metalarge2[metalarge2$family != "not_classified",]
large3Fmeta <- metalarge3[metalarge3$family != "not_classified",]
large4Fmeta <- metalarge4[metalarge4$family != "not_classified",]
sum(large0Fmeta$family =="not_classified")
dim(large0Fmeta)

# ORDER ei tarvitse poistoa

#piirteet
large0OR <- large0[,2:1001] #otetaan nimi kolumni pois sotkemasta
large1OR <- large1[,2:1001]
large2OR <- large2[,2:1001]
large3OR <- large3[,2:1001]
large4OR <- large4[,2:1001]

large0OR80 <- large080
large1OR80 <- large180
large2OR80 <- large280
large3OR80 <- large380
large4OR80 <- large480

large0ORmeta <- metalarge0
large1ORmeta <- metalarge1
large2ORmeta <- metalarge2
large3ORmeta <- metalarge3
large4ORmeta <- metalarge4

sum(large0ORmeta$order =="not_classified")
dim(large0ORmeta)

dim(large0OR)
dim(large0OR80)

################################################################################

##MidSize

################################################################################

#piirteet
mid0 <- read.csv("bioscan_bioclip_midsize_sample0_distances.csv")
mid1 <- read.csv("bioscan_bioclip_midsize_sample1_distances.csv")
mid2 <- read.csv("bioscan_bioclip_midsize_sample2_distances.csv")
mid3 <- read.csv("bioscan_bioclip_midsize_sample3_distances.csv")
mid4 <- read.csv("bioscan_bioclip_midsize_sample4_distances.csv")


#dna
mid080 <- read.csv("bioscan_midsize_sample0_k80dist.csv")
mid180 <- read.csv("bioscan_midsize_sample1_k80dist.csv")
mid280 <- read.csv("bioscan_midsize_sample2_k80dist.csv")
mid380 <- read.csv("bioscan_midsize_sample3_k80dist.csv")
mid480 <- read.csv("bioscan_midsize_sample4_k80dist.csv")


#meta
metamid0<-read_parquet("bioscan_metadata_midsize_sample0.parquet")
metamid1<-read_parquet("bioscan_metadata_midsize_sample1.parquet")
metamid2<-read_parquet("bioscan_metadata_midsize_sample2.parquet")
metamid3<-read_parquet("bioscan_metadata_midsize_sample3.parquet")
metamid4<-read_parquet("bioscan_metadata_midsize_sample4.parquet")

#tarkastetaan taksonit
sum((metamid0$subspecies)=="not_classified") #ei labeleitä, 1200 NC #HUOM ERI KOKOINEN KUIN MUUT SAMPLET n = 1200
sum((metamid0$species)=="not_classified") #NC 1100
sum((metamid0$genus )=="not_classified") #NC 907
sum((metamid0$tribe)=="not_classified") #NC 1181
sum((metamid0$subfamily)=="not_classified") #NC 899
sum((metamid0$family)=="not_classified") #NC 25
sum((metamid0$order)=="not_classified") #NC 0 = kaikki labelit paikassaan
#tarkastetaan order, family

length(unique(metamid2$family)) 
length(unique(metamid1$order))
table(metamid4$family) #51, 49, 49, 51, 55 eri tyyppiä samplessa 0-4
table(metamid0$order) #6, 6, 6, 6, 6



#poistetaan NC, family

mid0F <- mid0[,2:1201] #otetaan nimi kolumni pois sotkemasta
mid1F <- mid1[,2:1201]
mid2F <- mid2[,2:1201]
mid3F <- mid3[,2:1201]
mid4F <- mid4[,2:1201]

#piirteet
mid0F <- mid0F[metamid0$family != "not_classified",metamid0$family != "not_classified"]
mid1F <- mid1F[metamid1$family != "not_classified",metamid1$family != "not_classified"]
mid2F <- mid2F[metamid2$family != "not_classified",metamid2$family != "not_classified"]
mid3F <- mid3F[metamid3$family != "not_classified",metamid3$family != "not_classified"]
mid4F <- mid4F[metamid4$family != "not_classified",metamid4$family != "not_classified"]
#dna
mid0F80 <- mid080[metamid0$family != "not_classified",metamid0$family != "not_classified"]
mid1F80 <- mid180[metamid1$family != "not_classified",metamid1$family != "not_classified"]
mid2F80 <- mid280[metamid2$family != "not_classified",metamid2$family != "not_classified"]
mid3F80 <- mid380[metamid3$family != "not_classified",metamid3$family != "not_classified"]
mid4F80 <- mid480[metamid4$family != "not_classified",metamid4$family != "not_classified"]

dim(mid0F)
dim(mid1F)
dim(mid2F)
dim(mid3F)
dim(mid4F)
dim(mid0F80)


#meta
mid0Fmeta <- metamid0[metamid0$family != "not_classified",]
mid1Fmeta <- metamid1[metamid1$family != "not_classified",]
mid2Fmeta <- metamid2[metamid2$family != "not_classified",]
mid3Fmeta <- metamid3[metamid3$family != "not_classified",]
mid4Fmeta <- metamid4[metamid4$family != "not_classified",]
sum(mid0Fmeta$family =="not_classified")
dim(mid0Fmeta)

# ORDER ei tarvitse poistoa

#piirteet
mid0OR <- mid0[,2:1201] #otetaan nimi kolumni pois sotkemasta
mid1OR <- mid1[,2:1201]
mid2OR <- mid2[,2:1201]
mid3OR <- mid3[,2:1201]
mid4OR <- mid4[,2:1201]

mid0OR80 <- mid080
mid1OR80 <- mid180
mid2OR80 <- mid280
mid3OR80 <- mid380
mid4OR80 <- mid480

mid0ORmeta <- metamid0
mid1ORmeta <- metamid1
mid2ORmeta <- metamid2
mid3ORmeta <- metamid3
mid4ORmeta <- metamid4

sum(mid0ORmeta$order =="not_classified")
dim(mid0ORmeta)

dim(mid0OR)
dim(mid0OR80)

#k = #51, 49, 49, 51, 55 kun pienet mukana, fam 

################################################################################

##Sub100 (fam, order)

################################################################################

#piirteet
sub100 <- read.csv("bioscan_bioclip_sub100_sample0_distances.csv")

#dna
sub100_80 <- read.csv("bioscan_sub100_sample0_k80dist.csv")

#meta
metasub100<-read_parquet("bioscan_metadata_sub100_sample0.parquet")


#tarkastetaan taksonit
sum((metasub100$subspecies)=="not_classified") #243 kaikki puuttuu
sum((metasub100$species)=="not_classified") #NC 174
sum((metasub100$genus )=="not_classified") #NC 116
sum((metasub100$tribe)=="not_classified") #NC 243
sum((metasub100$subfamily)=="not_classified") #NC 120
sum((metasub100$family)=="not_classified") #NC 17
sum((metasub100$order)=="not_classified") #NC 0 = kaikki labelit paikassaan
#tarkastetaan order, family

length(unique(metasub100$family))  #19
length(unique(metasub100$order))  #8


table(metasub100$family)
table(metasub100$order)


#poistetaan NC family

#piir
sub100F <- sub100[,2:244]
sub100F <- sub100F[metasub100$family != "not_classified",metasub100$family != "not_classified"]
#dna
sub100F80 <- sub100_80[metasub100$family != "not_classified",metasub100$family != "not_classified"]
#meta
sub100Fmeta <- metasub100[metasub100$family != "not_classified",]

dim(sub100F)
dim(sub100F80)
dim(sub100Fmeta)


#orderista ei tarvitse poistaa NC, mutta otetaan sille oma piirre (kaikki on paikallaan)

sub100 <- sub100[,2:244]

