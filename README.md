# Ryhmittelykoodeja_2026_Spring

Sisällöstä ja ohjeita lukijalle:

Tiedostot sisältävät koodeja, joilla olen ryhmitellyt hyönteisiä (BioSCAN-1M aineistosta lahkon/lahkojen perusteella palauttamatta arvottuja otoksia).

Osa-aineistot: 

- AllEven, LargeSize, MidSize, Sub100 ovat yhdistelmä-tyyppisiä eli sisältävät useita lahkoja (toteutettu sekä lahko-, että heimotason ryhmittelyt)

- Coleoptera, Diptera, Hemiptera, Hymenoptera, Lepidoptera sis. yksittäisiä lahkoja (toteutettu vain heimotason ryhmittelyt)

Ryhmittelymenetelmät:
Diana (jakava hierarkkinen), Pam (k-metoids), Hclust (agglomeratiivinen hierarkkinen) 3 linkillä

Koodien ajosta uusilla otoksilla kiinnostuneen pyydän huomioimaan koodin otosten oletetun koon (1000x1000 matriiseja paitsi midSize:lle 1200x1200 ja sub100:lle 243x243). Tämä repo ei sisällä (ainakaan vielä) itse käytettämiäni otoksia, vaan ne kiinnostunut saa laskea itse.

Sisältö:
- Data_sisäänajo.R näyttää itse sisäänajon lisäksi, miten puuttuvaa tietoa on käsitelty ja datan labeleiden vaihtelua.
- Funktiot.R sisältää (toivottavasti) kaikki käytetyt funktiot ja sen on suositeltavaa olla ajettuna aina.
- Otoksien_koodit -kansiossa ovat koodit analyysien ajamiseen osa-aineistoittain.
- Tulosten_laskennan_ajot -kansio sisältää koodit automaattilukijoiden tuloksien saamiseksi (otoksittain ja osa-aineistoittain) -- huom, osa-aineisto (sekä sisään_ajo että otoksien_koodit) pitää tosiaan ajaa tätä varten. 
- indeksien_valit_laskenta.R toimii indeksien laskentaan sille osa-aineistolle, jonka analyysit on viimeksi ajettu kokonaisuudessaan.

Pärjäilkää. E. H. 2026 Kevät
