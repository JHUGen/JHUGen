      integer cc0,cc1,cc2,cc00,cc11,cc12,cc22,
     . cc001,cc002,cc111,cc112,cc122,cc222,
     . cc0000,cc0011,cc0012,cc0022,cc1111,cc1112,cc1122,cc1222,cc2222,
     . cc00001,cc00002,cc00111,cc00112,cc00122,cc00222,cc11111,cc11112,
     . cc11122,cc11222,cc12222,cc22222,
     . cc000000,cc000011,cc000012,cc000022,
     . cc001111,cc001112,cc001122,cc001222,cc002222,
     . cc111111,cc111112,cc111122,cc111222,cc112222,cc122222,cc222222,
     . cc0000001,cc0000002,cc0000111,cc0000112,cc0000122,cc0000222,
     . cc0011111,cc0011112,cc0011122,cc0011222,cc0012222,cc0022222,
     . cc1111111,cc1111112,cc1111122,cc1111222,cc1112222,cc1122222,
     . cc1222222,cc2222222,
     . Ncc,Pcc,Ncmax

      parameter(cc0=1)
      parameter(cc1=2)
      parameter(cc2=3)
      parameter(cc00=4)
      parameter(cc11=5)
      parameter(cc12=6)
      parameter(cc22=7)
      parameter(cc001=8)
      parameter(cc002=9)
      parameter(cc111=10)
      parameter(cc112=11)
      parameter(cc122=12)
      parameter(cc222=13)
      parameter(cc0000=14)
      parameter(cc0011=15)
      parameter(cc0012=16)
      parameter(cc0022=17)
      parameter(cc1111=18)
      parameter(cc1112=19)
      parameter(cc1122=20)
      parameter(cc1222=21)
      parameter(cc2222=22)
      parameter(cc00001=23)
      parameter(cc00002=24)
      parameter(cc00111=25)
      parameter(cc00112=26)
      parameter(cc00122=27)
      parameter(cc00222=28)
      parameter(cc11111=29)
      parameter(cc11112=30)
      parameter(cc11122=31)
      parameter(cc11222=32)
      parameter(cc12222=33)
      parameter(cc22222=34)
      parameter(cc000000=35)
      parameter(cc000011=36)
      parameter(cc000012=37)
      parameter(cc000022=38)
      parameter(cc001111=39)
      parameter(cc001112=40)
      parameter(cc001122=41)
      parameter(cc001222=42)
      parameter(cc002222=43)
      parameter(cc111111=44)
      parameter(cc111112=45)
      parameter(cc111122=46)
      parameter(cc111222=47)
      parameter(cc112222=48)
      parameter(cc122222=49)
      parameter(cc222222=50)

      parameter(cc0000001=51)
      parameter(cc0000002=52)
      parameter(cc0000111=53)
      parameter(cc0000112=54)
      parameter(cc0000122=55)
      parameter(cc0000222=56)
      parameter(cc0011111=57)
      parameter(cc0011112=58)
      parameter(cc0011122=59)
      parameter(cc0011222=60)
      parameter(cc0012222=61)
      parameter(cc0022222=62)
      parameter(cc1111111=63)
      parameter(cc1111112=64)
      parameter(cc1111122=65)
      parameter(cc1111222=66)
      parameter(cc1112222=67)
      parameter(cc1122222=68)
      parameter(cc1222222=69)
      parameter(cc2222222=70)



C---Number of different coeffs
      parameter(Ncc=70)
C---Number of parameters for massless triangle
      parameter(Pcc=6)
C---Number of different points in cache
      parameter(Ncmax=100)
