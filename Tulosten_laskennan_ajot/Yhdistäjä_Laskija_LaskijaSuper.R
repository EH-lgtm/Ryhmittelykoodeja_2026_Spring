#tulosten luenta osa-aineistoittain

#LASKIJA ORDER
#AUTOMATA
raja <- 0.05
raja <- 0.1
#sample0
AE0H1<-Yhdistäjä(Automata(ajo080or$D,raja),Automata(ajo0OR$D,raja))
AE0H2<-Yhdistäjä(Automata(ajo080or$E,raja),Automata(ajo0OR$E,raja))
AE0H3<-Yhdistäjä(Automata(ajo080or$f,raja),Automata(ajo0OR$f,raja))
AE0P<-Yhdistäjä(Automata(PAMajo0OR80$pam1Or,raja),Automata(PAMajo0OR$pam1Or ,raja))
AE0D<-Yhdistäjä(Automata(Dia0OR80$diana1Or,raja),Automata(Dia0OR$diana1Or,raja))

koko0 <-Laskija(AE0H1,AE0H2,AE0H3,AE0P,AE0D)
koko0
#sample1
AE1H1<-Yhdistäjä(Automata(ajo180or$D,raja),Automata(ajo1OR$D,raja))
AE1H2<-Yhdistäjä(Automata(ajo180or$E,raja),Automata(ajo1OR$E,raja))
AE1H3<-Yhdistäjä(Automata(ajo180or$f,raja),Automata(ajo1OR$f,raja))
AE1P<-Yhdistäjä(Automata(PAMajo1OR80$pam1Or,raja),Automata(PAMajo1OR$pam1Or,raja))
AE1D<-Yhdistäjä(Automata(Dia1OR80$diana1Or,raja),Automata(Dia1OR$diana1Or,raja))

koko1 <-Laskija(AE1H1,AE1H2,AE1H3,AE1P,AE1D)
#sample2
AE2H1<-Yhdistäjä(Automata(ajo280or$D,raja),Automata(ajo2OR$D,raja))
AE2H2<-Yhdistäjä(Automata(ajo280or$E,raja),Automata(ajo2OR$E,raja))
AE2H3<-Yhdistäjä(Automata(ajo280or$f,raja),Automata(ajo2OR$f,raja))
AE2P<-Yhdistäjä(Automata(PAMajo2OR80$pam1Or,raja),Automata(PAMajo2OR$pam1Or,raja))
AE2D<-Yhdistäjä(Automata(Dia2OR80$diana1Or,raja),Automata(Dia2OR$diana1Or,raja))

koko2 <-Laskija(AE2H1,AE2H2,AE2H3,AE2P,AE2D)
#sample3
AE3H1<-Yhdistäjä(Automata(ajo380or$D,raja),Automata(ajo3OR$D,raja))
AE3H2<-Yhdistäjä(Automata(ajo380or$E,raja),Automata(ajo3OR$E,raja))
AE3H3<-Yhdistäjä(Automata(ajo380or$f,raja),Automata(ajo3OR$f,raja))
AE3P<-Yhdistäjä(Automata(PAMajo3OR80$pam1Or,raja),Automata(PAMajo3OR$pam1Or,raja))
AE3D<-Yhdistäjä(Automata(Dia3OR80$diana1Or,raja),Automata(Dia3OR$diana1Or,raja))

koko3 <-Laskija(AE3H1,AE3H2,AE3H3,AE3P,AE3D)
#sample4
AE4H1<-Yhdistäjä(Automata(ajo480or$D,raja),Automata(ajo4OR$D,raja))
AE4H2<-Yhdistäjä(Automata(ajo480or$E,raja),Automata(ajo4OR$E,raja))
AE4H3<-Yhdistäjä(Automata(ajo480or$f,raja),Automata(ajo4OR$f,raja))
AE4P<-Yhdistäjä(Automata(PAMajo4OR80$pam1Or,raja),Automata(PAMajo4OR$pam1Or,raja))
AE4D<-Yhdistäjä(Automata(Dia4OR80$diana1Or,raja),Automata(Dia4OR$diana1Or,raja))

koko4 <-Laskija(AE4H1,AE4H2,AE4H3,AE4P,AE4D)

#yhteistulos
LaskijaSuper(list(koko1,koko2,koko3,koko4,koko0))


#LASKIJA
#AUTOMATA3

#sample0
A30H1<-Yhdistäjä(Automata3(ajo080or$D,raja),Automata3(ajo0OR$D,raja))
A30H2<-Yhdistäjä(Automata3(ajo080or$E,raja),Automata3(ajo0OR$E,raja))
A30H3<-Yhdistäjä(Automata3(ajo080or$f,raja),Automata3(ajo0OR$f,raja))
A30P<-Yhdistäjä(Automata3(PAMajo0OR80$pam1Or,raja),Automata3(PAMajo0OR$pam1Or,raja))
A30D<-Yhdistäjä(Automata3(Dia0OR80$diana1Or,raja),Automata3(Dia0OR$diana1Or,raja))

koko03 <-Laskija(A30H1,A30H2,A30H3,A30P,A30D)
koko03
#sample1
A31H1<-Yhdistäjä(Automata3(ajo180or$D,raja),Automata3(ajo1OR$D,raja))
A31H2<-Yhdistäjä(Automata3(ajo180or$E,raja),Automata3(ajo1OR$E,raja))
A31H3<-Yhdistäjä(Automata3(ajo180or$f,raja),Automata3(ajo1OR$f,raja))
A31P<-Yhdistäjä(Automata3(PAMajo1OR80$pam1Or,raja),Automata3(PAMajo1OR$pam1Or,raja))
A31D<-Yhdistäjä(Automata3(Dia1OR80$diana1Or,raja),Automata3(Dia1OR$diana1Or,raja))

koko13 <-Laskija(A31H1,A31H2,A31H3,A31P,A31D)
#sample2
A32H1<-Yhdistäjä(Automata3(ajo280or$D,raja),Automata3(ajo2OR$D,raja))
A32H2<-Yhdistäjä(Automata3(ajo280or$E,raja),Automata3(ajo2OR$E,raja))
A32H3<-Yhdistäjä(Automata3(ajo280or$f,raja),Automata3(ajo2OR$f,raja))
A32P<-Yhdistäjä(Automata3(PAMajo2OR80$pam1Or,raja),Automata3(PAMajo2OR$pam1Or,raja))
A32D<-Yhdistäjä(Automata3(Dia2OR80$diana1Or,raja),Automata3(Dia2OR$diana1Or,raja))

koko23 <-Laskija(A32H1,A32H2,A32H3,A32P,A32D)
#sample3
A33H1<-Yhdistäjä(Automata3(ajo380or$D,raja),Automata3(ajo3OR$D,raja))
A33H2<-Yhdistäjä(Automata3(ajo380or$E,raja),Automata3(ajo3OR$E,raja))
A33H3<-Yhdistäjä(Automata3(ajo380or$f,raja),Automata3(ajo3OR$f,raja))
A33P<-Yhdistäjä(Automata3(PAMajo3OR80$pam1Or,raja),Automata3(PAMajo3OR$pam1Or,raja))
A33D<-Yhdistäjä(Automata3(Dia3OR80$diana1Or,raja),Automata3(Dia3OR$diana1Or,raja))

koko33 <-Laskija(A33H1,A33H2,A33H3,A33P,A33D)
#sample4
A34H1<-Yhdistäjä(Automata3(ajo480or$D,raja),Automata3(ajo4OR$D,raja))
A34H2<-Yhdistäjä(Automata3(ajo480or$E,raja),Automata3(ajo4OR$E,raja))
A34H3<-Yhdistäjä(Automata3(ajo480or$f,raja),Automata3(ajo4OR$f,raja))
A34P<-Yhdistäjä(Automata3(PAMajo4OR80$pam1Or,raja),Automata3(PAMajo4OR$pam1Or,raja))
A34D<-Yhdistäjä(Automata3(Dia4OR80$diana1Or,raja),Automata3(Dia4OR$diana1Or,raja))

koko43 <-Laskija(A34H1,A34H2,A34H3,A34P,A34D)

LaskijaSuper(list(koko13,koko23,koko33,koko43,koko03))



#########################################
#LASKIJA FAMILY
#AUTOMATA
raja <- 0.05
raja <- 0.1
#sample0
AE0H1<-Yhdistäjä(Automata(ajo0_80$D,raja),Automata(ajo0f$D,raja))
AE0H2<-Yhdistäjä(Automata(ajo0_80$E,raja),Automata(ajo0f$E,raja))
AE0H3<-Yhdistäjä(Automata(ajo0_80$f,raja),Automata(ajo0f$f,raja))
AE0P<-Yhdistäjä(Automata(PAMajo0f80$pam1F,raja),Automata(PAMajo0f$pam1F,raja))
AE0D<-Yhdistäjä(Automata(Dia0F80$diana1F,raja),Automata(Dia0F$diana1F,raja))

koko0 <-Laskija(AE0H1,AE0H2,AE0H3,AE0P,AE0D)
koko0
#sample1
AE1H1<-Yhdistäjä(Automata(ajo1_80$D,raja),Automata(ajo1f$D,raja))
AE1H2<-Yhdistäjä(Automata(ajo1_80$E,raja),Automata(ajo1f$E,raja))
AE1H3<-Yhdistäjä(Automata(ajo1_80$f,raja),Automata(ajo1f$f,raja))
AE1P<-Yhdistäjä(Automata(PAMajo1f80$pam1F,raja),Automata(PAMajo1f$pam1F,raja))
AE1D<-Yhdistäjä(Automata(Dia1F80$diana1F,raja),Automata(Dia1F$diana1F,raja))

koko1 <-Laskija(AE1H1,AE1H2,AE1H3,AE1P,AE1D)
#sample2
AE2H1<-Yhdistäjä(Automata(ajo2_80$D,raja),Automata(ajo2f$D,raja))
AE2H2<-Yhdistäjä(Automata(ajo2_80$E,raja),Automata(ajo2f$E,raja))
AE2H3<-Yhdistäjä(Automata(ajo2_80$f,raja),Automata(ajo2f$f,raja))
AE2P<-Yhdistäjä(Automata(PAMajo2f80$pam1F,raja),Automata(PAMajo2f$pam1F,raja))
AE2D<-Yhdistäjä(Automata(Dia2F80$diana1F,raja),Automata(Dia2F$diana1F,raja))

koko2 <-Laskija(AE2H1,AE2H2,AE2H3,AE2P,AE2D)
#sample3
AE3H1<-Yhdistäjä(Automata(ajo3_80$D,raja),Automata(ajo3f$D,raja))
AE3H2<-Yhdistäjä(Automata(ajo3_80$E,raja),Automata(ajo3f$E,raja))
AE3H3<-Yhdistäjä(Automata(ajo3_80$f,raja),Automata(ajo3f$f,raja))
AE3P<-Yhdistäjä(Automata(PAMajo3f80$pam1F,raja),Automata(PAMajo3f$pam1F,raja))
AE3D<-Yhdistäjä(Automata(Dia3F80$diana1F,raja),Automata(Dia3F$diana1F,raja))

koko3 <-Laskija(AE3H1,AE3H2,AE3H3,AE3P,AE3D)
#sample4
AE4H1<-Yhdistäjä(Automata(ajo4_80$D,raja),Automata(ajo4f$D,raja))
AE4H2<-Yhdistäjä(Automata(ajo4_80$E,raja),Automata(ajo4f$E,raja))
AE4H3<-Yhdistäjä(Automata(ajo4_80$f,raja),Automata(ajo4f$f,raja))
AE4P<-Yhdistäjä(Automata(PAMajo4f80$pam1F,raja),Automata(PAMajo4f$pam1F,raja))
AE4D<-Yhdistäjä(Automata(Dia4F80$diana1F,raja),Automata(Dia4F$diana1F,raja))

koko4 <-Laskija(AE4H1,AE4H2,AE4H3,AE4P,AE4D)

#yhteistulos
LaskijaSuper(list(koko1,koko2,koko3,koko4,koko0))


#LASKIJA
#AUTOMATA3

#sample0
A30H1<-Yhdistäjä(Automata3(ajo0_80$D,raja),Automata3(ajo0f$D,raja))
A30H2<-Yhdistäjä(Automata3(ajo0_80$E,raja),Automata3(ajo0f$E,raja))
A30H3<-Yhdistäjä(Automata3(ajo0_80$f,raja),Automata3(ajo0f$f,raja))
A30P<-Yhdistäjä(Automata3(PAMajo0f80$pam1F,raja),Automata3(PAMajo0f$pam1F,raja))
A30D<-Yhdistäjä(Automata3(Dia0F80$diana1F,raja),Automata3(Dia0F$diana1F,raja))

koko03 <-Laskija(A30H1,A30H2,A30H3,A30P,A30D)
koko03
#sample1
A31H1<-Yhdistäjä(Automata3(ajo1_80$D,raja),Automata3(ajo1f$D,raja))
A31H2<-Yhdistäjä(Automata3(ajo1_80$E,raja),Automata3(ajo1f$E,raja))
A31H3<-Yhdistäjä(Automata3(ajo1_80$f,raja),Automata3(ajo1f$f,raja))
A31P<-Yhdistäjä(Automata3(PAMajo1f80$pam1F,raja),Automata3(PAMajo1f$pam1F,raja))
A31D<-Yhdistäjä(Automata3(Dia1F80$diana1F,raja),Automata3(Dia1F$diana1F,raja))

koko13 <-Laskija(A31H1,A31H2,A31H3,A31P,A31D)
#sample2
A32H1<-Yhdistäjä(Automata3(ajo2_80$D,raja),Automata3(ajo2f$D,raja))
A32H2<-Yhdistäjä(Automata3(ajo2_80$E,raja),Automata3(ajo2f$E,raja))
A32H3<-Yhdistäjä(Automata3(ajo2_80$f,raja),Automata3(ajo2f$f,raja))
A32P<-Yhdistäjä(Automata3(PAMajo2f80$pam1F,raja),Automata3(PAMajo2f$pam1F,raja))
A32D<-Yhdistäjä(Automata3(Dia2F80$diana1F,raja),Automata3(Dia2F$diana1F,raja))

koko23 <-Laskija(A32H1,A32H2,A32H3,A32P,A32D)
#sample3
A33H1<-Yhdistäjä(Automata3(ajo3_80$D,raja),Automata3(ajo3f$D,raja))
A33H2<-Yhdistäjä(Automata3(ajo3_80$E,raja),Automata3(ajo3f$E,raja))
A33H3<-Yhdistäjä(Automata3(ajo3_80$f,raja),Automata3(ajo3f$f,raja))
A33P<-Yhdistäjä(Automata3(PAMajo3f80$pam1F,raja),Automata3(PAMajo3f$pam1F,raja))
A33D<-Yhdistäjä(Automata3(Dia3F80$diana1F,raja),Automata3(Dia3F$diana1F,raja))

koko33 <-Laskija(A33H1,A33H2,A33H3,A33P,A33D)
#sample4
A34H1<-Yhdistäjä(Automata3(ajo4_80$D,raja),Automata3(ajo4f$D,raja))
A34H2<-Yhdistäjä(Automata3(ajo4_80$E,raja),Automata3(ajo4f$E,raja))
A34H3<-Yhdistäjä(Automata3(ajo4_80$f,raja),Automata3(ajo4f$f,raja))
A34P<-Yhdistäjä(Automata3(PAMajo4f80$pam1F,raja),Automata3(PAMajo4f$pam1F,raja))
A34D<-Yhdistäjä(Automata3(Dia4F80$diana1F,raja),Automata3(Dia4F$diana1F,raja))

koko43 <-Laskija(A34H1,A34H2,A34H3,A34P,A34D)

LaskijaSuper(list(koko13,koko23,koko33,koko43,koko03))

