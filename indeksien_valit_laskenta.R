#indeksien välit (min,mean,max)

#FAMILY
H1f<-indeksiyhd_0(indeksilaskuri(ajo0_80$a,ajo0f$a),indeksilaskuri(ajo1_80$a,ajo1f$a),indeksilaskuri(ajo2_80$a,ajo2f$a),indeksilaskuri(ajo3_80$a,ajo3f$a),indeksilaskuri(ajo4_80$a,ajo4f$a))
H1f
H2f<-indeksiyhd_0(indeksilaskuri(ajo0_80$b,ajo0f$b),indeksilaskuri(ajo1_80$b,ajo1f$b),indeksilaskuri(ajo2_80$b,ajo2f$b),indeksilaskuri(ajo3_80$b,ajo3f$b),indeksilaskuri(ajo4_80$b,ajo4f$b))
H2f
H3f<-indeksiyhd_0(indeksilaskuri(ajo0_80$c,ajo0f$c),indeksilaskuri(ajo1_80$c,ajo1f$c),indeksilaskuri(ajo2_80$c,ajo2f$c),indeksilaskuri(ajo3_80$c,ajo3f$c),indeksilaskuri(ajo4_80$c,ajo4f$c))
H3f
Pf<-indeksiyhd_0(indeksilaskuri(PAMajo0f80$a, PAMajo0f$a),indeksilaskuri(PAMajo1f80$a, PAMajo1f$a),indeksilaskuri(PAMajo2f80$a, PAMajo2f$a),indeksilaskuri(PAMajo3f80$a, PAMajo3f$a),indeksilaskuri(PAMajo4f80$a, PAMajo4f$a))
Pf
Df<-indeksiyhd_0(indeksilaskuri(Dia0F80$a,Dia0F$a),indeksilaskuri(Dia1F80$a,Dia1F$a),indeksilaskuri(Dia2F80$a,Dia2F$a),indeksilaskuri(Dia3F80$a,Dia3F$a),indeksilaskuri(Dia4F80$a,Dia4F$a))
Df

H1f
H2f
H3f
Pf
Df

#ORDER MUISTA VAIHTAA k ja ajaa ajot uudestaan, kun tarpeen!

H1o<-indeksiyhd_0(indeksilaskuri(ajo080or$a, ajo0OR$a),indeksilaskuri(ajo180or$a, ajo1OR$a),indeksilaskuri(ajo280or$a, ajo2OR$a),indeksilaskuri(ajo380or$a, ajo3OR$a),indeksilaskuri(ajo480or$a, ajo4OR$a))
H1o
H2o<-indeksiyhd_0(indeksilaskuri(ajo080or$b, ajo0OR$b),indeksilaskuri(ajo180or$b, ajo1OR$b),indeksilaskuri(ajo280or$b, ajo2OR$b),indeksilaskuri(ajo380or$b, ajo3OR$b),indeksilaskuri(ajo480or$b, ajo4OR$b))
H2o
H3o<-indeksiyhd_0(indeksilaskuri(ajo080or$c, ajo0OR$c),indeksilaskuri(ajo180or$c, ajo1OR$c),indeksilaskuri(ajo280or$c, ajo2OR$c),indeksilaskuri(ajo380or$c, ajo3OR$c),indeksilaskuri(ajo480or$c, ajo4OR$c))
H3o
Po<-indeksiyhd_0(indeksilaskuri(PAMajo0OR80$a, PAMajo0OR$a),indeksilaskuri(PAMajo1OR80$a, PAMajo1OR$a),indeksilaskuri(PAMajo2OR80$a, PAMajo2OR$a),indeksilaskuri(PAMajo3OR80$a, PAMajo3OR$a),indeksilaskuri(PAMajo4OR80$a, PAMajo4OR$a))
Po
Do<-indeksiyhd_0(indeksilaskuri(Dia0OR80$a,Dia0OR$a),indeksilaskuri(Dia1OR80$a,Dia1OR$a),indeksilaskuri(Dia2OR80$a,Dia2OR$a),indeksilaskuri(Dia3OR80$a,Dia3OR$a),indeksilaskuri(Dia4OR80$a,Dia4OR$a))
Do

H1o
H2o
H3o
Po
Do