
C   Dv(N+dd1,Dv(N+dd2,Dv(N+dd3,z)
c     do ep=-2,0
c     in(1,ep)=
c    & +f1*Dv(N+dd0,ep)
c    & +Cv(cc0+C134,ep)
c    & -Cv(cc0+C234,ep)
c     in(2,ep)=
c    & +f2*Dv(N+dd0,ep)
c    & +Cv(cc0+C124,ep)
c    & -Cv(cc0+C234,ep)
c     in(3,ep)=
c    & +f3*Dv(N+dd0,ep)
c    & +Cv(cc0+C123,ep)
c    & -Cv(cc0+C234,ep)



C---TWO

C   Dv(N+dd11,Dv(N+dd12,Dv(N+dd13,p)
c     do ep=-2,0
c     in(1,ep)=
c    & +f1*Dv(N+dd1,ep)
c    & +Csum0(ep)
c    & -2*Dv(N+dd00,ep)
c     in(2,ep)=
c    & +f2*Dv(N+dd1,ep)
c    & +Csum0(ep)
c    & +Cv(cc1+C124,ep)
c     in(3,ep)=
c    & +f3*Dv(N+dd1,ep)
c    & +Csum0(ep)
c    & +Cv(cc1+C123,ep)

c'B'Id,Csum0(P?,K?,L?,m1?,m2?,m3?)=
c'B'  +Cv(cc0+C234,ep)
c'B'  +Cv(cc1+C234,ep)
c'B'  +Cv(cc2+C234,ep)


C   Dv(N+dd12,Dv(N+dd22,Dv(N+dd23,k)
c     in(1,ep)=
c    & +f1*Dv(N+dd2,ep)
c    & +Cv(cc1+C134,ep)
c    & -Cv(cc1+C234,ep)

c     in(2,ep)=
c    & +f2*Dv(N+dd2,ep)
c    & -Cv(cc1+C234,ep)
c    & -2*Dv(N+dd00,ep)

c     in(3,ep)=
c    & +f3*Dv(N+dd2,ep)
c    & -Cv(cc1+C234,ep)
c    & +Cv(cc2+C123,ep)

c     enddo

c     call pvBackSubst(G, 3, perm, in)

c     do ep=-2,0
c     Dv(N+dd133333,ep)=in(1,ep)
c     Dv(N+dd233333,ep)=in(2,ep)
c     Dv(N+dd333333,ep)=in(3,ep)

C   Dv(N+dd13,Dv(N+dd23,Dv(N+dd33,l)
c     in(1,ep)=
c    & +f1*Dv(N+dd3,ep)
c    & +Cv(cc2+C134,ep)
c    & -Cv(cc2+C234,ep)

c     in(2,ep)=
c    & +f2*Dv(N+dd3,ep)
c    & +Cv(cc2+C124,ep)
c    & -Cv(cc2+C234,ep)
c     in(3,ep)=
c    & +f3*Dv(N+dd3,ep)
c    & -Cv(cc2+C234,ep)
c    & -2*Dv(N+dd00,ep)
c     enddo

c     call pvBackSubst(G, 3, perm, in)

c     do ep=-2,0
c     Dv(N+dd133333,ep)=in(1,ep)
c     Dv(N+dd233333,ep)=in(2,ep)
c     Dv(N+dd333333,ep)=in(3,ep)


C  THREE

C   Dv(N+dd111,Dv(N+dd112,Dv(N+dd113,pp)
c     in(1,ep)=
c    & +f1*Dv(N+dd11,ep)
c    & -Csum0(ep)
c    & -Csum1(ep)
c    & -Csum2(ep)
c    & -4*Dv(N+dd001,ep)

c     in(2,ep)=
c    & +f2*Dv(N+dd11,ep)
c    & -Csum0(ep)
c    & -Csum1(ep)
c    & -Csum2(ep)
c    & +Cv(cc11+C124,ep)

c     in(3,ep)=
c    & +f3*Dv(N+dd11,ep)
c    & -Csum0(ep)
c    & -Csum1(ep)
c    & -Csum2(ep)
c    & +Cv(cc11+C123,ep)

c     enddo

c     call pvBackSubst(G, 3, perm, in)

c     do ep=-2,0
c     Dv(N+dd133333,ep)=in(1,ep)
c     Dv(N+dd233333,ep)=in(2,ep)
c     Dv(N+dd333333,ep)=in(3,ep)
c'B'Id,Csum1(P?,K?,L?,m1?,m2?,m3?)=
c'B'    Cv(cc1+C234,ep)
c'B'   +Cv(cc11+C234,ep)
c'B'   +Cv(cc12+C234,ep)
c'B'Id,Csum2(P?,K?,L?,m1?,m2?,m3?)=
c'B'    Cv(cc2+C234,ep)
c'B'   +Cv(cc12+C234,ep)
c'B'   +Cv(cc22+C234,ep)
  

C   Dv(N+dd112,Dv(N+dd122,Dv(N+dd123,pk)
c     in(1,ep)=
c    & +f1*Dv(N+dd12,ep)
c    & +Csum1(ep)
c    & -2*Dv(N+dd002,ep)

c     in(2,ep)=
c    & +f2*Dv(N+dd12,ep)
c    & +Csum1(ep)
c    & -2*Dv(N+dd001,ep)
c     in(3,ep)=
c    & +f3*Dv(N+dd12,ep)
c    & +Csum1(ep)
c    & +Cv(cc12+C123,ep)
c     enddo

c     call pvBackSubst(G, 3, perm, in)

c     do ep=-2,0
c     Dv(N+dd133333,ep)=in(1,ep)
c     Dv(N+dd233333,ep)=in(2,ep)
c     Dv(N+dd333333,ep)=in(3,ep)


C   Dv(N+dd113,Dv(N+dd123,Dv(N+dd133,pl)
c     in(1,ep)=
c    & +f1*Dv(N+dd13,ep)
c    & +Csum2(ep)
c    & -2*Dv(N+dd003,ep)
c     in(2,ep)=
c    & +f2*Dv(N+dd13,ep)
c    & +Csum2(ep)
c    & +Cv(cc12+C124,ep)
c     in(3,ep)=
c    & +f3*Dv(N+dd13,ep)
c    & +Csum2(ep)
c    & -2*Dv(N+dd001,ep)
 
C   Dv(N+dd123,Dv(N+dd223,Dv(N+dd233,kl)
c     in(1,ep)=
c    & +f1*Dv(N+dd23,ep)
c    & +Cv(cc12+C134,ep)
c    & -Cv(cc12+C234,ep)
c     in(2,ep)=
c    & +f2*Dv(N+dd23,ep)
c    & -Cv(cc12+C234,ep)
c    & -2*Dv(N+dd003,ep)
c     in(3,ep)=
c    & +f3*Dv(N+dd23,ep)
c    & -Cv(cc12+C234,ep)
c    & -2*Dv(N+dd002,ep)
c     enddo

c     call pvBackSubst(G, 3, perm, in)

c     do ep=-2,0
c     Dv(N+dd133333,ep)=in(1,ep)
c     Dv(N+dd233333,ep)=in(2,ep)
c     Dv(N+dd333333,ep)=in(3,ep)


C   Dv(N+dd122,Dv(N+dd222,Dv(N+dd223,kk)
c     in(1,ep)=
c    & +f1*Dv(N+dd22,ep)
c    & +Cv(cc11+C134,ep)
c    & -Cv(cc11+C234,ep)

c     in(2,ep)=
c    & +f2*Dv(N+dd22,ep)
c    & -Cv(cc11+C234,ep)
c    & -4*Dv(N+dd002,ep)
c     in(3,ep)=
c    & +f3*Dv(N+dd22,ep)
c    & -Cv(cc11+C234,ep)
c    & +Cv(cc22+C123,ep)
c     enddo

c     call pvBackSubst(G, 3, perm, in)

c     do ep=-2,0
c     Dv(N+dd133333,ep)=in(1,ep)
c     Dv(N+dd233333,ep)=in(2,ep)
c     Dv(N+dd333333,ep)=in(3,ep)


C   Dv(N+dd133,Dv(N+dd233,Dv(N+dd333,ll)
c     in(1,ep)=
c    & +f1*Dv(N+dd33,ep)
c    & +Cv(cc22+C134,ep)
c    & -Cv(cc22+C234,ep)
c     in(2,ep)=
c    & +f2*Dv(N+dd33,ep)
c    & +Cv(cc22+C124,ep)
c    & -Cv(cc22+C234,ep)
c     in(3,ep)=
c    & +f3*Dv(N+dd33,ep)
c    & -Cv(cc22+C234,ep)
c    & -4*Dv(N+dd003,ep)
c     enddo

c     call pvBackSubst(G, 3, perm, in)

c     do ep=-2,0
c     Dv(N+dd133333,ep)=in(1,ep)
c     Dv(N+dd233333,ep)=in(2,ep)
c     Dv(N+dd333333,ep)=in(3,ep)

C   Dv(N+dd001,Dv(N+dd002,Dv(N+dd003,zz)
c     in(1,ep)=
c    & +f1*Dv(N+dd00,ep)
c    & +Cv(cc00+C134,ep)
c    & -Cv(cc00+C234,ep)
c     in(2,ep)=
c    & +f2*Dv(N+dd00,ep)
c    & +Cv(cc00+C124,ep)
c    & -Cv(cc00+C234,ep)
c     in(3,ep)=
c    & +f3*Dv(N+dd00,ep)
c    & +Cv(cc00+C123,ep)
c    & -Cv(cc00+C234,ep)
c     enddo

c     call pvBackSubst(G, 3, perm, in)

c     do ep=-2,0
c     Dv(N+dd133333,ep)=in(1,ep)
c     Dv(N+dd233333,ep)=in(2,ep)
c     Dv(N+dd333333,ep)=in(3,ep)







C FOUR
C   Dv(N+dd1111,Dv(N+dd1112,Dv(N+dd1113,ppp)
c     in(1,ep)=
c    & +f1*Dv(N+dd111,ep)
c    & +Csum0(ep)
c    & +2*Csum1(ep)
c    & +2*Csum2(ep)
c    & +2*Csum12(ep)
c    & +Csum11(ep)
c    & +Csum22(ep)
c    & -6*Dv(N+dd0011,ep)
c     in(2,ep)=
c    & +f2*Dv(N+dd111,ep)
c    & +Csum0(ep)
c    & +2*Csum1(ep)
c    & +2*Csum2(ep)
c    & +2*Csum12(ep)
c    & +Csum11(ep)
c    & +Csum22(ep)
c    & +Cv(cc111+C124,ep)
c     in(3,ep)=
c    & +f3*Dv(N+dd111,ep)
c    & +Csum0(ep)
c    & +2*Csum1(ep)
c    & +2*Csum2(ep)
c    & +2*Csum12(ep)
c    & +Csum11(ep)
c    & +Csum22(ep)
c    & +Cv(cc111+C123,ep)
c     enddo

c     call pvBackSubst(G, 3, perm, in)

c     do ep=-2,0
c     Dv(N+dd133333,ep)=in(1,ep)
c     Dv(N+dd233333,ep)=in(2,ep)
c     Dv(N+dd333333,ep)=in(3,ep)

C   Dv(N+dd1222,Dv(N+dd2222,Dv(N+dd2223,kkk)
c     in(1,ep)=
c    & +f1*Dv(N+dd222,ep)
c    & +Cv(cc111+C134,ep)
c    & -Cv(cc111+C234,ep)
c     in(2,ep)=
c    & +f2*Dv(N+dd222,ep)
c    & -Cv(cc111+C234,ep)
c    & -6*Dv(N+dd0022,ep)
c     in(3,ep)=
c    & +f3*Dv(N+dd222,ep)
c    & -Cv(cc111+C234,ep)
c    & +Cv(cc222+C123,ep)
c     enddo

c     call pvBackSubst(G, 3, perm, in)

c     do ep=-2,0
c     Dv(N+dd133333,ep)=in(1,ep)
c     Dv(N+dd233333,ep)=in(2,ep)
c     Dv(N+dd333333,ep)=in(3,ep)


C   Dv(N+dd1333,Dv(N+dd2333,Dv(N+dd3333,lll)
c     in(1,ep)=
c    & +f1*Dv(N+dd333,ep)
c    & +Cv(cc222+C134,ep)
c    & -Cv(cc222+C234,ep)
c     in(2,ep)=
c    & +f2*Dv(N+dd333,ep)
c    & +Cv(cc222+C124,ep)
c    & -Cv(cc222+C234,ep)
c     in(3,ep)=
c    & +f3*Dv(N+dd333,ep)
c    & -Cv(cc222+C234,ep)
c    & -6*Dv(N+dd0033,ep)
c     enddo

c     call pvBackSubst(G, 3, perm, in)

c     do ep=-2,0
c     Dv(N+dd133333,ep)=in(1,ep)
c     Dv(N+dd233333,ep)=in(2,ep)
c     Dv(N+dd333333,ep)=in(3,ep)


C   Dv(N+dd1112,Dv(N+dd1122,Dv(N+dd1123,ppk)
c     in(1,ep)=
c    & +f1*Dv(N+dd112,ep)
c    & -Csum1(ep)
c    & -Csum11(ep)
c    & -Csum12(ep)
c    & -4*Dv(N+dd0012,ep)
c     in(2,ep)=
c    & +f2*Dv(N+dd112,ep)
c    & -Csum1(ep)
c    & -Csum11(ep)
c    & -Csum12(ep)
c    & -2*Dv(N+dd0011,ep)
c     in(3,ep)=
c    & +f3*Dv(N+dd112,ep)
c    & -Csum1(ep)
c    & -Csum11(ep)
c    & -Csum12(ep)
c    & +Cv(cc112+C123,ep)
c     enddo

c     call pvBackSubst(G, 3, perm, in)

c     do ep=-2,0
c     Dv(N+dd133333,ep)=in(1,ep)
c     Dv(N+dd233333,ep)=in(2,ep)
c     Dv(N+dd333333,ep)=in(3,ep)


C   Dv(N+dd1113,Dv(N+dd1123,Dv(N+dd1133,ppl)
c     in(1,ep)=
c    & +f1*Dv(N+dd113,ep)
c    & -Csum2(ep)
c    & -Csum22(ep)
c    & -Csum12(ep)
c    & -4*Dv(N+dd0013,ep)
c     in(2,ep)=
c    & +f2*Dv(N+dd113,ep)
c    & -Csum2(ep)
c    & -Csum22(ep)
c    & -Csum12(ep)
c    & +Cv(cc112+C124,ep)
c     in(3,ep)=
c    & +f3*Dv(N+dd113,ep)
c    & -Csum2(ep)
c    & -Csum22(ep)
c    & -Csum12(ep)
c    & -2*Dv(N+dd0011,ep)
c     enddo

c     call pvBackSubst(G, 3, perm, in)

c     do ep=-2,0
c     Dv(N+dd133333,ep)=in(1,ep)
c     Dv(N+dd233333,ep)=in(2,ep)
c     Dv(N+dd333333,ep)=in(3,ep)


C   Dv(N+dd1122,Dv(N+dd1222,Dv(N+dd1223,pkk)
c     in(1,ep)=
c    & +f1*Dv(N+dd122,ep)
c    & +Csum11(ep)
c    & -2*Dv(N+dd0022,ep)

c     in(2,ep)=
c    & +f2*Dv(N+dd122,ep)
c    & +Csum11(ep)
c    & -4*Dv(N+dd0012,ep)
c     in(3,ep)=
c    & +f3*Dv(N+dd122,ep)
c    & +Csum11(ep)
c    & +Cv(cc122+C123,ep)
c     enddo

c     call pvBackSubst(G, 3, perm, in)

c     do ep=-2,0
c     Dv(N+dd133333,ep)=in(1,ep)
c     Dv(N+dd233333,ep)=in(2,ep)
c     Dv(N+dd333333,ep)=in(3,ep)

C   Dv(N+dd1133,Dv(N+dd1233,Dv(N+dd1333,pll)
c     in(1,ep)=
c    & +f1*Dv(N+dd133,ep)
c    & +Csum22(ep)
c    & -2*Dv(N+dd0033,ep)
c     in(2,ep)=
c    & +f2*Dv(N+dd133,ep)
c    & +Csum22(ep)
c    & +Cv(cc122+C124,ep)
c     in(3,ep)=
c    & +f3*Dv(N+dd133,ep)
c    & +Csum22(ep)
c    & -4*Dv(N+dd0013,ep)
c     enddo

c     call pvBackSubst(G, 3, perm, in)

c     do ep=-2,0
c     Dv(N+dd133333,ep)=in(1,ep)
c     Dv(N+dd233333,ep)=in(2,ep)
c     Dv(N+dd333333,ep)=in(3,ep)


C   Dv(N+dd1223,Dv(N+dd2223,Dv(N+dd2233,kkl)
c     in(1,ep)=
c    & +f1*Dv(N+dd223,ep)
c    & +Cv(cc112+C134,ep)
c    & -Cv(cc112+C234,ep)
c     in(2,ep)=
c    & +f2*Dv(N+dd223,ep)
c    & -Cv(cc112+C234,ep)
c    & -4*Dv(N+dd0023,ep)
c     in(3,ep)=
c    & +f3*Dv(N+dd223,ep)
c    & -Cv(cc112+C234,ep)
c    & -2*Dv(N+dd0022,ep)
c     enddo

c     call pvBackSubst(G, 3, perm, in)

c     do ep=-2,0
c     Dv(N+dd133333,ep)=in(1,ep)
c     Dv(N+dd233333,ep)=in(2,ep)
c     Dv(N+dd333333,ep)=in(3,ep)


C   Dv(N+dd1233,Dv(N+dd2233,Dv(N+dd2333,kll)
c     in(1,ep)=
c    & +f1*Dv(N+dd233,ep)
c    & +Cv(cc122+C134,ep)
c    & -Cv(cc122+C234,ep)

c    &
c     in(2,ep)=
c    & +f2*Dv(N+dd233,ep)
c    & -Cv(cc122+C234,ep)
c    & -2*Dv(N+dd0033,ep)
c     in(3,ep)=
c    & +f3*Dv(N+dd233,ep)
c    & -Cv(cc122+C234,ep)
c    & -4*Dv(N+dd0023,ep)
c     enddo

c     call pvBackSubst(G, 3, perm, in)

c     do ep=-2,0
c     Dv(N+dd133333,ep)=in(1,ep)
c     Dv(N+dd233333,ep)=in(2,ep)
c     Dv(N+dd333333,ep)=in(3,ep)


C   Dv(N+dd1123,Dv(N+dd1223,Dv(N+dd1233,pkl)
c     in(1,ep)=
c    & +f1*Dv(N+dd123,ep)
c    & +Csum12(ep)
c    & -2*Dv(N+dd0023,ep)
c     in(2,ep)=
c    & +f2*Dv(N+dd123,ep)
c    & +Csum12(ep)
c    & -2*Dv(N+dd0013,ep)
c     in(3,ep)=
c    & +f3*Dv(N+dd123,ep)
c    & +Csum12(ep)
c    & -2*Dv(N+dd0012,ep)
c     enddo

c     call pvBackSubst(G, 3, perm, in)

c     do ep=-2,0
c     Dv(N+dd133333,ep)=in(1,ep)
c     Dv(N+dd233333,ep)=in(2,ep)
c     Dv(N+dd333333,ep)=in(3,ep)



C   Dv(N+dd0011,Dv(N+dd0012,Dv(N+dd0013,zzp)
c     in(1,ep)=
c    & +f1*Dv(N+dd001,ep)
c    & +Csum00(ep)
c    & -2*Dv(N+dd0000,ep)
c     in(2,ep)=
c    & +f2*Dv(N+dd001,ep)
c    & +Csum00(ep)
c    & +Cv(cc001+C124,ep)
c     in(3,ep)=
c    & +f3*Dv(N+dd001,ep)
c    & +Csum00(ep)
c    & +Cv(cc001+C123,ep)
c     enddo

c     call pvBackSubst(G, 3, perm, in)

c     do ep=-2,0
c     Dv(N+dd133333,ep)=in(1,ep)
c     Dv(N+dd233333,ep)=in(2,ep)
c     Dv(N+dd333333,ep)=in(3,ep)

C  Id,Csum11(P?,K?,L?,m1?,m2?,m3?)=
Cc      & +Cv(cc11+C234,ep)
Cc      & +Cv(cc111+C234,ep)
Cc      & +Cv(cc112+C234,ep)
C  Id,Csum12(P?,K?,L?,m1?,m2?,m3?)=Cv(cc12+C234,ep)
Cc      &c    &   +Cv(cc112+C234,ep)
Cc      &c    &   +Cv(cc122+C234,ep)
C  Id,Csum22(P?,K?,L?,m1?,m2?,m3?)=Cv(cc22+C234,ep)
Cc      &c    &   +Cv(cc122+C234,ep)
Cc      &c    &   +Cv(cc222+C234,ep)
C  Id,Csum00(P?,K?,L?,m1?,m2?,m3?)=Cv(cc00+C234,ep)
Cc      &c    &   +Cv(cc001+C234,ep)
Cc      &c    &   +Cv(cc002+C234,ep)


C   Dv(N+dd0012,Dv(N+dd0022,Dv(N+dd0023,zzk)
c     in(1,ep)=
c    & +f1*Dv(N+dd002,ep)
c    & +Cv(cc001+C134,ep)
c    & -Cv(cc001+C234,ep)
c     in(2,ep)=
c    & +f2*Dv(N+dd002,ep)
c    & -Cv(cc001+C234,ep)
c    & -2*Dv(N+dd0000,ep)
c     in(3,ep)=
c    & +f3*Dv(N+dd002,ep)
c    & -Cv(cc001+C234,ep)
c    & +Cv(cc002+C123,ep)
c     enddo

c     call pvBackSubst(G, 3, perm, in)

c     do ep=-2,0
c     Dv(N+dd133333,ep)=in(1,ep)
c     Dv(N+dd233333,ep)=in(2,ep)
c     Dv(N+dd333333,ep)=in(3,ep)


C   Dv(N+dd0013,Dv(N+dd0023,Dv(N+dd0033,zzl)
c     in(1,ep)=
c    & +f1*Dv(N+dd003,ep)
c    & +Cv(cc002+C134,ep)
c    & -Cv(cc002+C234,ep)
c     in(2,ep)=
c    & +f2*Dv(N+dd003,ep)
c    & +Cv(cc002+C124,ep)
c    & -Cv(cc002+C234,ep)
c     in(3,ep)=
c    & +f3*Dv(N+dd003,ep)
c    & -Cv(cc002+C234,ep)
c    & -2*Dv(N+dd0000,ep)
c     enddo

c     call pvBackSubst(G, 3, perm, in)

c     do ep=-2,0
c     Dv(N+dd133333,ep)=in(1,ep)
c     Dv(N+dd233333,ep)=in(2,ep)
c     Dv(N+dd333333,ep)=in(3,ep)


C****  FIVE

C   Dv(N+dd11111,Dv(N+dd11112,Dv(N+dd11113,pppp)

c     in(1,ep)=
c    & +f1*Dv(N+dd1111,ep)
c    & -8*Dv(N+dd00111,ep)
c    & -3*Csum22(ep)
c    & -3*Csum11(ep)
c    & -3*Csum122(ep)
c    & -3*Csum112(ep)
c    & -Csum111(ep)
c    & -Csum222(ep)
c    & -6*Csum12(ep)
c    & -3*Csum2(ep)
c    & -3*Csum1(ep)
c    & -Csum0(ep)

c     in(2,ep)=
c    & +f2*Dv(N+dd1111,ep)
c    & +Cv(cc1111+C124,ep)
c    & -3*Csum22(ep)
c    & -3*Csum11(ep)
c    & -3*Csum122(ep)
c    & -3*Csum112(ep)
c    & -Csum111(ep)
c    & -Csum222(ep)
c    & -6*Csum12(ep)
c    & -3*Csum2(ep)
c    & -3*Csum1(ep)
c    & -Csum0(ep)

c     in(3,ep)=
c    & +f3*Dv(N+dd1111,ep)
c    & +Cv(cc1111+C123,ep)
c    & -3*Csum22(ep)
c    & -3*Csum11(ep)
c    & -3*Csum122(ep)
c    & -3*Csum112(ep)
c    & -Csum111(ep)
c    & -Csum222(ep)
c    & -6*Csum12(ep)
c    & -3*Csum2(ep)
c    & -3*Csum1(ep)
c    & -Csum0(ep)

c     enddo

c     call pvBackSubst(G, 3, perm, in)

c     do ep=-2,0
c     Dv(N+dd133333,ep)=in(1,ep)
c     Dv(N+dd233333,ep)=in(2,ep)
c     Dv(N+dd333333,ep)=in(3,ep)




C   Dv(N+dd11112,Dv(N+dd11122,Dv(N+dd11123,pppk)
c     in(1,ep)=
c    & +f1*Dv(N+dd1112,ep)
c    & -6*Dv(N+dd00112,ep)
c    & +Csum1(ep)
c    & +Csum111(ep)
c    & +Csum122(ep)
c    & +2*Csum112(ep)
c    & +2*Csum11(ep)
c    & +2*Csum12(ep)

c     in(2,ep)=
c    & +f2*Dv(N+dd1112,ep)
c    & -2*Dv(N+dd00111,ep)
c    & +Csum1(ep)
c    & +Csum111(ep)
c    & +Csum122(ep)
c    & +2*Csum112(ep)
c    & +2*Csum11(ep)
c    & +2*Csum12(ep)

c     in(3,ep)=
c    & +f3*Dv(N+dd1112,ep)
c    & +Cv(cc1112+C123,ep)
c    & +Csum1(ep)
c    & +Csum111(ep)
c    & +Csum122(ep)
c    & +2*Csum112(ep)
c    & +2*Csum11(ep)
c    & +2*Csum12(ep)
c     enddo

c     call pvBackSubst(G, 3, perm, in)

c     do ep=-2,0
c     Dv(N+dd133333,ep)=in(1,ep)
c     Dv(N+dd233333,ep)=in(2,ep)
c     Dv(N+dd333333,ep)=in(3,ep)


C   Dv(N+dd11113,Dv(N+dd11123,Dv(N+dd11133,pppl)
c     in(1,ep)=
c    & +f1*Dv(N+dd1113,ep)
c    & -6*Dv(N+dd00113,ep)
c    & +Csum2(ep)
c    & +Csum112(ep)
c    & +Csum222(ep)
c    & +2*Csum122(ep)
c    & +2*Csum12(ep)
c    & +2*Csum22(ep)


c     in(2,ep)=
c    & +f2*Dv(N+dd1113,ep)
c    & +Cv(cc1112+C124,ep)
c    & +Csum2(ep)
c    & +Csum112(ep)
c    & +Csum222(ep)
c    & +2*Csum122(ep)
c    & +2*Csum12(ep)
c    & +2*Csum22(ep)

c     in(3,ep)=
c    & +f3*Dv(N+dd1113,ep)
c    & -2*Dv(N+dd00111,ep)
c    & +Csum2(ep)
c    & +Csum112(ep)
c    & +Csum222(ep)
c    & +2*Csum122(ep)
c    & +2*Csum12(ep)
c    & +2*Csum22(ep)
c     enddo

c     call pvBackSubst(G, 3, perm, in)

c     do ep=-2,0
c     Dv(N+dd133333,ep)=in(1,ep)
c     Dv(N+dd233333,ep)=in(2,ep)
c     Dv(N+dd333333,ep)=in(3,ep)


C   Dv(N+dd11133,Dv(N+dd11233,Dv(N+dd11333,ppll)

c     in(1,ep)=
c    & +f1*Dv(N+dd1133,ep)
c    & -4*Dv(N+dd00133,ep)
c    & -Csum22(ep)
c    & -Csum122(ep)
c    & -Csum222(ep)
c     in(2,ep)=
c    & +f2*Dv(N+dd1133,ep)
c    & +Cv(cc1122+C124,ep)
c    & -Csum22(ep)
c    & -Csum122(ep)
c    & -Csum222(ep)
c     in(3,ep)=
c    & +f3*Dv(N+dd1133,ep)
c    & -4*Dv(N+dd00113,ep)
c    & -Csum22(ep)
c    & -Csum122(ep)
c    & -Csum222(ep)
c     enddo

c     call pvBackSubst(G, 3, perm, in)

c     do ep=-2,0
c     Dv(N+dd133333,ep)=in(1,ep)
c     Dv(N+dd233333,ep)=in(2,ep)
c     Dv(N+dd333333,ep)=in(3,ep)


C   Dv(N+dd11122,Dv(N+dd11222,Dv(N+dd11223,ppkk)
c     in(1,ep)=
c    & +f1*Dv(N+dd1122,ep)
c    & -Csum11(ep)
c    & -Csum111(ep)
c    & -Csum112(ep)
c    & -4*Dv(N+dd00122,ep)
c     in(2,ep)=
c    & +f2*Dv(N+dd1122,ep)
c    & -Csum11(ep)
c    & -Csum111(ep)
c    & -Csum112(ep)
c    & -4*Dv(N+dd00112,ep)
c     in(3,ep)=
c    & +f3*Dv(N+dd1122,ep)
c    & -Csum11(ep)
c    & -Csum111(ep)
c    & -Csum112(ep)
c    & +Cv(cc1122+C123,ep)
c     enddo

c     call pvBackSubst(G, 3, perm, in)

c     do ep=-2,0
c     Dv(N+dd133333,ep)=in(1,ep)
c     Dv(N+dd233333,ep)=in(2,ep)
c     Dv(N+dd333333,ep)=in(3,ep)

C   Dv(N+dd11233,Dv(N+dd12233,Dv(N+dd12333,pkll)
c     in(1,ep)=
c    & +f1*Dv(N+dd1233,ep)
c    & +Csum122(ep)
c    & -2*Dv(N+dd00233,ep)

c     in(2,ep)=
c    & +f2*Dv(N+dd1233,ep)
c    & +Csum122(ep)
c    & -2*Dv(N+dd00133,ep)

c     in(3,ep)=
c    & +f3*Dv(N+dd1233,ep)
c    & +Csum122(ep)
c    & -4*Dv(N+dd00123,ep)
c     enddo

c     call pvBackSubst(G, 3, perm, in)

c     do ep=-2,0
c     Dv(N+dd133333,ep)=in(1,ep)
c     Dv(N+dd233333,ep)=in(2,ep)
c     Dv(N+dd333333,ep)=in(3,ep)

C   Dv(N+dd11333,Dv(N+dd12333,Dv(N+dd13333,plll)
c     in(1,ep)=
c    & +f1*Dv(N+dd1333,ep)
c    & +Csum222(ep)
c    & -2*Dv(N+dd00333,ep)

c     in(2,ep)=
c    & +f2*Dv(N+dd1333,ep)
c    & +Csum222(ep)
c    & +Cv(cc1222+C124,ep)

c     in(3,ep)=
c    & +f3*Dv(N+dd1333,ep)
c    & +Csum222(ep)
c    & -6*Dv(N+dd00133,ep)
c     enddo

c     call pvBackSubst(G, 3, perm, in)

c     do ep=-2,0
c     Dv(N+dd133333,ep)=in(1,ep)
c     Dv(N+dd233333,ep)=in(2,ep)
c     Dv(N+dd333333,ep)=in(3,ep)

C   Dv(N+dd11222,Dv(N+dd12222,Dv(N+dd12223,pkkk)
c     in(1,ep)=
c    & +f1*Dv(N+dd1222,ep)
c    & +Csum111(ep)
c    & -2*Dv(N+dd00222,ep)


c     in(2,ep)=
c    & +f2*Dv(N+dd1222,ep)
c    & +Csum111(ep)
c    & -6*Dv(N+dd00122,ep)

c     in(3,ep)=
c    & +f3*Dv(N+dd1222,ep)
c    & +Csum111(ep)
c    & +Cv(cc1222+C123,ep)
c     enddo

c     call pvBackSubst(G, 3, perm, in)

c     do ep=-2,0
c     Dv(N+dd133333,ep)=in(1,ep)
c     Dv(N+dd233333,ep)=in(2,ep)
c     Dv(N+dd333333,ep)=in(3,ep)

C   Dv(N+dd12222,Dv(N+dd22222,Dv(N+dd22223,kkkk)
c     in(1,ep)=
c    & +f1*Dv(N+dd2222,ep)
c    & +Cv(cc1111+C134,ep)
c    & -Cv(cc1111+C234,ep)
c     in(2,ep)=
c    & +f2*Dv(N+dd2222,ep)
c    & -Cv(cc1111+C234,ep)
c    & -8*Dv(N+dd00222,ep)
c     in(3,ep)=
c    & +f3*Dv(N+dd2222,ep)
c    & -Cv(cc1111+C234,ep)
c    & +Cv(cc2222+C123,ep)
c     enddo

c     call pvBackSubst(G, 3, perm, in)

c     do ep=-2,0
c     Dv(N+dd133333,ep)=in(1,ep)
c     Dv(N+dd233333,ep)=in(2,ep)
c     Dv(N+dd333333,ep)=in(3,ep)



C   Dv(N+dd11223,Dv(N+dd12223,Dv(N+dd12233,pkkl)
c     in(1,ep)=
c    & +f1*Dv(N+dd1223,ep)
c    & +Csum112(ep)
c    & -2*Dv(N+dd00223,ep)
c     in(2,ep)=
c    & +f2*Dv(N+dd1223,ep)
c    & +Csum112(ep)
c    & -4*Dv(N+dd00123,ep)
c     in(3,ep)=
c    & +f3*Dv(N+dd1223,ep)
c    & +Csum112(ep)
c    & -2*Dv(N+dd00122,ep)
c     enddo

c     call pvBackSubst(G, 3, perm, in)

c     do ep=-2,0
c     Dv(N+dd133333,ep)=in(1,ep)
c     Dv(N+dd233333,ep)=in(2,ep)
c     Dv(N+dd333333,ep)=in(3,ep)


C   Dv(N+dd12223,Dv(N+dd22223,Dv(N+dd22233,kkkl)
c     in(1,ep)=
c    & +f1*Dv(N+dd2223,ep)
c    & +Cv(cc1112+C134,ep)
c    & -Cv(cc1112+C234,ep)
c     in(2,ep)=
c    & +f2*Dv(N+dd2223,ep)
c    & -Cv(cc1112+C234,ep)
c    & -6*Dv(N+dd00223,ep)

c     in(3,ep)=
c    & +f3*Dv(N+dd2223,ep)
c    & -Cv(cc1112+C234,ep)
c    & -2*Dv(N+dd00222,ep)
c     enddo

c     call pvBackSubst(G, 3, perm, in)

c     do ep=-2,0
c     Dv(N+dd133333,ep)=in(1,ep)
c     Dv(N+dd233333,ep)=in(2,ep)
c     Dv(N+dd333333,ep)=in(3,ep)


C   Dv(N+dd13333,Dv(N+dd23333,Dv(N+dd33333,llll)
c     in(1,ep)=
c    & +f1*Dv(N+dd3333,ep)
c    & +Cv(cc2222+C134,ep)
c    & -Cv(cc2222+C234,ep)
c     in(2,ep)=
c    & +f2*Dv(N+dd3333,ep)
c    & +Cv(cc2222+C124,ep)
c    & -Cv(cc2222+C234,ep)
c     in(3,ep)=
c    & +f3*Dv(N+dd3333,ep)
c    & -Cv(cc2222+C234,ep)
c    & -8*Dv(N+dd00333,ep)
c     enddo

c     call pvBackSubst(G, 3, perm, in)

c     do ep=-2,0
c     Dv(N+dd133333,ep)=in(1,ep)
c     Dv(N+dd233333,ep)=in(2,ep)
c     Dv(N+dd333333,ep)=in(3,ep)

C   Dv(N+dd11122,Dv(N+dd11222,Dv(N+dd11223,ppkk)




C   Dv(N+dd12333,Dv(N+dd22333,Dv(N+dd23333,klll)
c     in(1,ep)=
c    & +f1*Dv(N+dd2333,ep)
c    & +Cv(cc1222+C134,ep)
c    & -Cv(cc1222+C234,ep)
c     in(2,ep)=
c    & +f2*Dv(N+dd2333,ep)
c    & -Cv(cc1222+C234,ep)
c    & -2*Dv(N+dd00333,ep)

c     in(3,ep)=
c    & +f3*Dv(N+dd2333,ep)
c    & -Cv(cc1222+C234,ep)
c    & -6*Dv(N+dd00233,ep)

c     enddo

c     call pvBackSubst(G, 3, perm, in)

c     do ep=-2,0
c     Dv(N+dd133333,ep)=in(1,ep)
c     Dv(N+dd233333,ep)=in(2,ep)
c     Dv(N+dd333333,ep)=in(3,ep)
C   Dv(N+dd12233,Dv(N+dd22233,Dv(N+dd22333,kkll)
c     in(1,ep)=
c    & +f1*Dv(N+dd2233,ep)
c    & +Cv(cc1122+C134,ep)
c    & -Cv(cc1122+C234,ep)
c     in(2,ep)=
c    & +f2*Dv(N+dd2233,ep)
c    & -Cv(cc1122+C234,ep)
c    & -4*Dv(N+dd00233,ep)

c     in(3,ep)=
c    & +f3*Dv(N+dd2233,ep)
c    & -Cv(cc1122+C234,ep)
c    & -4*Dv(N+dd00223,ep)
c     enddo

c     call pvBackSubst(G, 3, perm, in)

c     do ep=-2,0
c     Dv(N+dd133333,ep)=in(1,ep)
c     Dv(N+dd233333,ep)=in(2,ep)
c     Dv(N+dd333333,ep)=in(3,ep)


C   Dv(N+dd00111,Dv(N+dd00112,Dv(N+dd00113,zzpp)

c     in(1,ep)=
c    & +f1*Dv(N+dd0011,ep)
c    & -Csum00(ep)
c    & -Csum001(ep)
c    & -Csum002(ep)
c    & -4*Dv(N+dd00001,ep)


c     in(2,ep)=
c    & +f2*Dv(N+dd0011,ep)
c    & -Csum00(ep)
c    & -Csum001(ep)
c    & -Csum002(ep)
c    & +Cv(cc0011+C124,ep)
c     in(3,ep)=
c    & +f3*Dv(N+dd0011,ep)
c    & -Csum00(ep)
c    & -Csum001(ep)
c    & -Csum002(ep)
c    & +Cv(cc0011+C123,ep)
c     enddo

c     call pvBackSubst(G, 3, perm, in)

c     do ep=-2,0
c     Dv(N+dd133333,ep)=in(1,ep)
c     Dv(N+dd233333,ep)=in(2,ep)
c     Dv(N+dd333333,ep)=in(3,ep)

C   Dv(N+dd00112,Dv(N+dd00122,Dv(N+dd00123,zzpk)

c     in(1,ep)=
c    & +f1*Dv(N+dd0012,ep)
c    & +Csum001(ep)
c    & -2*Dv(N+dd00002,ep)

c     in(2,ep)=
c    & +f2*Dv(N+dd0012,ep)
c    & +Csum001(ep)
c    & -2*Dv(N+dd00001,ep)
c     in(3,ep)=
c    & +f3*Dv(N+dd0012,ep)
c    & +Csum001(ep)
c    & +Cv(cc0012+C123,ep)
c     enddo

c     call pvBackSubst(G, 3, perm, in)

c     do ep=-2,0
c     Dv(N+dd133333,ep)=in(1,ep)
c     Dv(N+dd233333,ep)=in(2,ep)
c     Dv(N+dd333333,ep)=in(3,ep)

C   Dv(N+dd00122,Dv(N+dd00222,Dv(N+dd00223,zzkk)


c     in(1,ep)=
c    & +f1*Dv(N+dd0022,ep)
c    & +Cv(cc0011+C134,ep)
c    & -Cv(cc0011+C234,ep)

c     in(2,ep)=
c    & +f2*Dv(N+dd0022,ep)
c    & -Cv(cc0011+C234,ep)
c    & -4*Dv(N+dd00002,ep)

c     in(3,ep)=
c    & +f3*Dv(N+dd0022,ep)
c    & -Cv(cc0011+C234,ep)
c    & +Cv(cc0022+C123,ep)
c     enddo

c     call pvBackSubst(G, 3, perm, in)

c     do ep=-2,0
c     Dv(N+dd133333,ep)=in(1,ep)
c     Dv(N+dd233333,ep)=in(2,ep)
c     Dv(N+dd333333,ep)=in(3,ep)



C   Dv(N+dd00123,Dv(N+dd002223,Dv(N+dd00233,zzkl)

c     in(1,ep)=
c    & +f1*Dv(N+dd0023,ep)
c    & +Cv(cc0012+C134,ep)
c    & -Cv(cc0012+C234,ep)

c     in(2,ep)=
c    & +f2*Dv(N+dd0023,ep)
c    & -Cv(cc0012+C234,ep)
c    & -2*Dv(N+dd00003,ep)
c     in(3,ep)=
c    & +f3*Dv(N+dd0023,ep)
c    & -Cv(cc0012+C234,ep)
c    & -2*Dv(N+dd00002,ep)
c     enddo

c     call pvBackSubst(G, 3, perm, in)

c     do ep=-2,0
c     Dv(N+dd133333,ep)=in(1,ep)
c     Dv(N+dd233333,ep)=in(2,ep)
c     Dv(N+dd333333,ep)=in(3,ep)

C   Dv(N+dd00113,Dv(N+dd00123,Dv(N+dd00133,zzpl)


c     in(1,ep)=
c    & +f1*Dv(N+dd0013,ep)
c    & +Csum002(ep)
c    & -2*Dv(N+dd00003,ep)

c     in(2,ep)=
c    & +f2*Dv(N+dd0013,ep)
c    & +Csum002(ep)
c    & +Cv(cc0012+C124,ep)


c     in(3,ep)=
c    & +f3*Dv(N+dd0013,ep)
c    & +Csum002(ep)
c    & -2*Dv(N+dd00001,ep)
c     enddo

c     call pvBackSubst(G, 3, perm, in)

c     do ep=-2,0
c     Dv(N+dd133333,ep)=in(1,ep)
c     Dv(N+dd233333,ep)=in(2,ep)
c     Dv(N+dd333333,ep)=in(3,ep)

C   Dv(N+dd00133,Dv(N+dd00233,Dv(N+dd00333,zzll)

c     in(1,ep)=
c    & +f1*Dv(N+dd0033,ep)
c    & +Cv(cc0022+C134,ep)
c    & -Cv(cc0022+C234,ep)



c     in(2,ep)=
c    & +f2*Dv(N+dd0033,ep)
c    & +Cv(cc0022+C124,ep)
c    & -Cv(cc0022+C234,ep)


c     in(3,ep)=
c    & +f3*Dv(N+dd0033,ep)
c    & -Cv(cc0022+C234,ep)
c    & -4*Dv(N+dd00003,ep)
c     enddo

c     call pvBackSubst(G, 3, perm, in)

c     do ep=-2,0
c     Dv(N+dd133333,ep)=in(1,ep)
c     Dv(N+dd233333,ep)=in(2,ep)
c     Dv(N+dd333333,ep)=in(3,ep)


C   Dv(N+dd00123,Dv(N+dd00223,Dv(N+dd00233,ppkl)
c     in(1,ep)=
c    & +f1*Dv(N+dd1123,ep)
c    & -Csum12(ep)
c    & -Csum112(ep)
c    & -Csum122(ep)
c    & -4*Dv(N+dd00123,ep)

c     in(2,ep)=
c    & +f2*Dv(N+dd1123,ep)
c    & -Csum12(ep)
c    & -Csum112(ep)
c    & -Csum122(ep)
c    & -2*Dv(N+dd00113,ep)

c     in(3,ep)=
c    & +f3*Dv(N+dd1123,ep)
c    & -Csum12(ep)
c    & -Csum112(ep)
c    & -Csum122(ep)
c    & -2*Dv(N+dd00112,ep)
c     enddo

c     call pvBackSubst(G, 3, perm, in)

c     do ep=-2,0
c     Dv(N+dd133333,ep)=in(1,ep)
c     Dv(N+dd233333,ep)=in(2,ep)
c     Dv(N+dd333333,ep)=in(3,ep)

C  Id,Csum111(P?,K?,L?,m1?,m2?,m3?)=
Cc      & +Cv(cc111+C234,ep)
Cc      & +Cv(cc1111+C234,ep)
Cc      & +Cv(cc1112+C234,ep)
C  Id,Csum112(P?,K?,L?,m1?,m2?,m3?)=
Cc      & +Cv(cc112+C234,ep)
Cc      & +Cv(cc1112+C234,ep)
Cc      & +Cv(cc1122+C234,ep)
C  Id,Csum122(P?,K?,L?,m1?,m2?,m3?)=
Cc      & +Cv(cc122+C234,ep)
Cc      & +Cv(cc1122+C234,ep)
Cc      & +Cv(cc1222+C234,ep)
C  Id,Csum222(P?,K?,L?,m1?,m2?,m3?)=
Cc      & +Cv(cc222+C234,ep)
Cc      & +Cv(cc1222+C234,ep)
Cc      & +Cv(cc2222+C234,ep)
C  Id,Csum001(P?,K?,L?,m1?,m2?,m3?)=
Cc      & +Cv(cc001+C234,ep)
Cc      & +Cv(cc0011+C234,ep)
Cc      & +Cv(cc0012+C234,ep)
C  Id,Csum002(P?,K?,L?,m1?,m2?,m3?)=
Cc      & +Cv(cc002+C234,ep)
Cc      & +Cv(cc0012+C234,ep)
Cc      & +Cv(cc0022+C234,ep)


C   Dv(N+dd00001,Dv(N+dd00002,Dv(N+dd00003,zzzz)
c     in(1,ep)=
c    & +f1*Dv(N+dd0000,ep)
c    &  +Cv(cc0000+C134,ep)
c    & -Cv(cc0000+C234,ep)
c     in(2,ep)=
c    & +f2*Dv(N+dd0000,ep)
c    & +Cv(cc0000+C124,ep)
c    & -Cv(cc0000+C234,ep)
c     in(3,ep)=
c    & +f3*Dv(N+dd0000,ep)
c    & +Cv(cc0000+C123,ep)
c    & -Cv(cc0000+C234,ep)
c     enddo

c     call pvBackSubst(G, 3, perm, in)

c     do ep=-2,0
c     Dv(N+dd133333,ep)=in(1,ep)
c     Dv(N+dd233333,ep)=in(2,ep)
c     Dv(N+dd333333,ep)=in(3,ep)




*********SIX
C   Dv(N+dd111111,Dv(N+dd111112,Dv(N+dd111113,ppppp)
c     in(1,ep)=
c    & +f1*Dv(N+dd11111,ep)
c    & -10*Dv(N+dd001111,ep)
c    & +Csum0(ep)
c    & +4*Csum2(ep)
c    & +4*Csum1(ep)
c    & +4*Csum111(ep)
c    & +4*Csum222(ep)
c    & +6*Csum11(ep)
c    & +6*Csum22(ep)
c    & +12*Csum12(ep)
c    & +12*Csum112(ep)
c    & +12*Csum122(ep)
c    & +4*Csum1222(ep)
c    & +6*Csum1122(ep)
c    & +4*Csum1112(ep)
c    & +Csum1111(ep)
c    & +Csum2222(ep)

c     in(2,ep)=
c    & +f2*Dv(N+dd11111,ep)
c    & +Cv(cc11111+C124,ep)
c    & +Csum0(ep)
c    & +4*Csum2(ep)
c    & +4*Csum1(ep)
c    & +4*Csum111(ep)
c    & +4*Csum222(ep)
c    & +6*Csum11(ep)
c    & +6*Csum22(ep)
c    & +12*Csum12(ep)
c    & +12*Csum112(ep)
c    & +12*Csum122(ep)
c    & +4*Csum1222(ep)
c    & +6*Csum1122(ep)
c    & +4*Csum1112(ep)
c    & +Csum1111(ep)
c    & +Csum2222(ep)

c     in(3,ep)=
c    & +f3*Dv(N+dd11111,ep)
c    & +Cv(cc11111+C123,ep)
c    & +Csum0(ep)
c    & +4*Csum2(ep)
c    & +4*Csum1(ep)
c    & +4*Csum111(ep)
c    & +4*Csum222(ep)
c    & +6*Csum11(ep)
c    & +6*Csum22(ep)
c    & +12*Csum12(ep)
c    & +12*Csum112(ep)
c    & +12*Csum122(ep)
c    & +4*Csum1222(ep)
c    & +6*Csum1122(ep)
c    & +4*Csum1112(ep)
c    & +Csum1111(ep)
c    & +Csum2222(ep)
c     enddo

c     call pvBackSubst(G, 3, perm, in)

c     do ep=-2,0
c     Dv(N+dd133333,ep)=in(1,ep)
c     Dv(N+dd233333,ep)=in(2,ep)
c     Dv(N+dd333333,ep)=in(3,ep)


C   Dv(N+dd111112,Dv(N+dd111122,Dv(N+dd111123,ppppk)
c     in(1,ep)=
c    & +f1*Dv(N+dd11112,ep)
c    & -8*Dv(N+dd001112,ep)
c    & -Csum1(ep)
c    & -3*Csum11(ep)
c    & -3*Csum12(ep)
c    & -3*Csum122(ep)
c    & -6*Csum112(ep)
c    & -3*Csum111(ep)
c    & -3*Csum1112(ep)
c    & -3*Csum1122(ep)
c    & -Csum1111(ep)
c    & -Csum1222(ep)

c     in(2,ep)=
c    & +f2*Dv(N+dd11112,ep)
c    & -2*Dv(N+dd001111,ep)
c    & -Csum1(ep)
c    & -3*Csum11(ep)
c    & -3*Csum12(ep)
c    & -3*Csum122(ep)
c    & -6*Csum112(ep)
c    & -3*Csum111(ep)
c    & -3*Csum1112(ep)
c    & -3*Csum1122(ep)
c    & -Csum1111(ep)
c    & -Csum1222(ep)

c     in(3,ep)=
c    & +f3*Dv(N+dd11112,ep)
c    & +Cv(cc11112+C123,ep)
c    & -Csum1(ep)
c    & -3*Csum11(ep)
c    & -3*Csum12(ep)
c    & -3*Csum122(ep)
c    & -6*Csum112(ep)
c    & -3*Csum111(ep)
c    & -3*Csum1112(ep)
c    & -3*Csum1122(ep)
c    & -Csum1111(ep)
c    & -Csum1222(ep)
c     enddo

c     call pvBackSubst(G, 3, perm, in)

c     do ep=-2,0
c     Dv(N+dd133333,ep)=in(1,ep)
c     Dv(N+dd233333,ep)=in(2,ep)
c     Dv(N+dd333333,ep)=in(3,ep)


C   Dv(N+dd111113,Dv(N+dd111123,Dv(N+dd111133,ppppl)
c     in(1,ep)=
c    & +f1*Dv(N+dd11113,ep)
c    & -8*Dv(N+dd001113,ep)
c    & -Csum2(ep)
c    & -3*Csum12(ep)
c    & -3*Csum22(ep)
c    & -3*Csum222(ep)
c    & -6*Csum122(ep)
c    & -3*Csum112(ep)
c    & -3*Csum1122(ep)
c    & -3*Csum1222(ep)
c    & -Csum1112(ep)
c    & -Csum2222(ep)

c     in(2,ep)=
c    & +f2*Dv(N+dd11113,ep)
c    & +Cv(cc11112+C124,ep)
c    & -Csum2(ep)
c    & -3*Csum12(ep)
c    & -3*Csum22(ep)
c    & -3*Csum222(ep)
c    & -6*Csum122(ep)
c    & -3*Csum112(ep)
c    & -3*Csum1122(ep)
c    & -3*Csum1222(ep)
c    & -Csum1112(ep)
c    & -Csum2222(ep)
c     in(3,ep)=
c    & +f3*Dv(N+dd11113,ep)
c    & -2*Dv(N+dd001111,ep)
c    & -Csum2(ep)
c    & -3*Csum12(ep)
c    & -3*Csum22(ep)
c    & -3*Csum222(ep)
c    & -6*Csum122(ep)
c    & -3*Csum112(ep)
c    & -3*Csum1122(ep)
c    & -3*Csum1222(ep)
c    & -Csum1112(ep)
c    & -Csum2222(ep)
c     enddo

c     call pvBackSubst(G, 3, perm, in)

c     do ep=-2,0
c     Dv(N+dd133333,ep)=in(1,ep)
c     Dv(N+dd233333,ep)=in(2,ep)
c     Dv(N+dd333333,ep)=in(3,ep)

C   Dv(N+dd111122,Dv(N+dd111222,Dv(N+dd111223,pppkk)
c     in(1,ep)=
c    & +f1*Dv(N+dd11122,ep)
c    & -6*Dv(N+dd001122,ep)
c    & +Csum11(ep)
c    & +2*Csum111(ep)
c    & +2*Csum112(ep)
c    & +2*Csum1112(ep)
c    & +Csum1111(ep)
c    & +Csum1122(ep)
c     in(2,ep)=
c    & +f2*Dv(N+dd11122,ep)
c    & -4*Dv(N+dd001112,ep)
c    & +Csum11(ep)
c    & +2*Csum111(ep)
c    & +2*Csum112(ep)
c    & +2*Csum1112(ep)
c    & +Csum1111(ep)
c    & +Csum1122(ep)

c     in(3,ep)=
c    & +f3*Dv(N+dd11122,ep)
c    & +Cv(cc11122+C123,ep)
c    & +Csum11(ep)
c    & +2*Csum111(ep)
c    & +2*Csum112(ep)
c    & +2*Csum1112(ep)
c    & +Csum1111(ep)
c    & +Csum1122(ep)
c     enddo

c     call pvBackSubst(G, 3, perm, in)

c     do ep=-2,0
c     Dv(N+dd133333,ep)=in(1,ep)
c     Dv(N+dd233333,ep)=in(2,ep)
c     Dv(N+dd333333,ep)=in(3,ep)


C   Dv(N+dd111123,Dv(N+dd111223,Dv(N+dd111233,pppkl)
c     in(1,ep)=
c    & +f1*Dv(N+dd11123,ep)
c    & -6*Dv(N+dd001123,ep)
c    & +Csum12(ep)
c    & +2*Csum112(ep)
c    & +2*Csum122(ep)
c    & +2*Csum1122(ep)
c    & +Csum1112(ep)
c    & +Csum1222(ep)
c     in(2,ep)=
c    & +f2*Dv(N+dd11123,ep)
c    & -2*Dv(N+dd001113,ep)
c    & +Csum12(ep)
c    & +2*Csum112(ep)
c    & +2*Csum122(ep)
c    & +2*Csum1122(ep)
c    & +Csum1112(ep)
c    & +Csum1222(ep)

c     in(3,ep)=
c    & +f3*Dv(N+dd11123,ep)
c    & -2*Dv(N+dd001112,ep)
c    & +Csum12(ep)
c    & +2*Csum112(ep)
c    & +2*Csum122(ep)
c    & +2*Csum1122(ep)
c    & +Csum1112(ep)
c    & +Csum1222(ep)
c     enddo

c     call pvBackSubst(G, 3, perm, in)

c     do ep=-2,0
c     Dv(N+dd133333,ep)=in(1,ep)
c     Dv(N+dd233333,ep)=in(2,ep)
c     Dv(N+dd333333,ep)=in(3,ep)

C   Dv(N+dd111133,Dv(N+dd111233,Dv(N+dd111333,pppll)
c     in(1,ep)=
c    & +f1*Dv(N+dd11133,ep)
c    & -6*Dv(N+dd001133,ep)
c    & +Csum22(ep)
c    & +2*Csum122(ep)
c    & +2*Csum222(ep)
c    & +2*Csum1222(ep)
c    & +Csum1122(ep)
c    & +Csum2222(ep)
c     in(2,ep)=
c    & +f2*Dv(N+dd11133,ep)
c    & +Cv(cc11122+C124,ep)
c    & +Csum22(ep)
c    & +2*Csum122(ep)
c    & +2*Csum222(ep)
c    & +2*Csum1222(ep)
c    & +Csum1122(ep)
c    & +Csum2222(ep)
c     in(3,ep)=
c    & +f3*Dv(N+dd11133,ep)
c    & -4*Dv(N+dd001113,ep)
c    & +Csum22(ep)
c    & +2*Csum122(ep)
c    & +2*Csum222(ep)
c    & +2*Csum1222(ep)
c    & +Csum1122(ep)
c    & +Csum2222(ep)
c     enddo

c     call pvBackSubst(G, 3, perm, in)

c     do ep=-2,0
c     Dv(N+dd133333,ep)=in(1,ep)
c     Dv(N+dd233333,ep)=in(2,ep)
c     Dv(N+dd333333,ep)=in(3,ep)



C   Dv(N+dd111222,Dv(N+dd112222,Dv(N+dd112223,ppkkk)
c     in(1,ep)=
c    & +f1*Dv(N+dd11222,ep)
c    & -4*Dv(N+dd001222,ep)
c    & -Csum111(ep)
c    & -Csum1112(ep)
c    & -Csum1111(ep)
c     in(2,ep)=
c    & +f2*Dv(N+dd11222,ep)
c    & -6*Dv(N+dd001122,ep)
c    & -Csum111(ep)
c    & -Csum1112(ep)
c    & -Csum1111(ep)
c     in(3,ep)=
c    & +f3*Dv(N+dd11222,ep)
c    & +Cv(cc11222+C123,ep)
c    & -Csum111(ep)
c    & -Csum1112(ep)
c    & -Csum1111(ep)
c     enddo

c     call pvBackSubst(G, 3, perm, in)

c     do ep=-2,0
c     Dv(N+dd133333,ep)=in(1,ep)
c     Dv(N+dd233333,ep)=in(2,ep)
c     Dv(N+dd333333,ep)=in(3,ep)


C   Dv(N+dd112233,Dv(N+dd122233,Dv(N+dd122333,pkkll)
c     in(1,ep)=
c    & +f1*Dv(N+dd12233,ep)
c    & -2*Dv(N+dd002233,ep)
c    & +Csum1122(ep)
c     in(2,ep)=
c    & +f2*Dv(N+dd12233,ep)
c    & -4*Dv(N+dd001233,ep)
c    & +Csum1122(ep)
c     in(3,ep)=
c    & +f3*Dv(N+dd12233,ep)
c    & -4*Dv(N+dd001223,ep)
c    & +Csum1122(ep)
c     enddo

c     call pvBackSubst(G, 3, perm, in)

c     do ep=-2,0
c     Dv(N+dd133333,ep)=in(1,ep)
c     Dv(N+dd233333,ep)=in(2,ep)
c     Dv(N+dd333333,ep)=in(3,ep)

C   Dv(N+dd111223,Dv(N+dd112223,Dv(N+dd112233,ppkkl)
c     in(1,ep)=
c    & +f1*Dv(N+dd11223,ep)
c    & -4*Dv(N+dd001223,ep)
c    & -Csum112(ep)
c    & -Csum1122(ep)
c    & -Csum1112(ep)

c     in(2,ep)=
c    & +f2*Dv(N+dd11223,ep)
c    & -4*Dv(N+dd001123,ep)
c    & -Csum112(ep)
c    & -Csum1122(ep)
c    & -Csum1112(ep)
c     in(3,ep)=
c    & +f3*Dv(N+dd11223,ep)
c    & -2*Dv(N+dd001122,ep)
c    & -Csum112(ep)
c    & -Csum1122(ep)
c    & -Csum1112(ep)
c     enddo

c     call pvBackSubst(G, 3, perm, in)

c     do ep=-2,0
c     Dv(N+dd133333,ep)=in(1,ep)
c     Dv(N+dd233333,ep)=in(2,ep)
c     Dv(N+dd333333,ep)=in(3,ep)


C   Dv(N+dd111233,Dv(N+dd112233,Dv(N+dd112333,ppkll)
c     in(1,ep)=
c    & +f1*Dv(N+dd11233,ep)
c    & -4*Dv(N+dd001233,ep)
c    & -Csum122(ep)
c    & -Csum1222(ep)
c    & -Csum1122(ep)
c     in(2,ep)=
c    & +f2*Dv(N+dd11233,ep)
c    & -2*Dv(N+dd001133,ep)
c    & -Csum122(ep)
c    & -Csum1222(ep)
c    & -Csum1122(ep)
c     in(3,ep)=
c    & +f3*Dv(N+dd11233,ep)
c    & -4*Dv(N+dd001123,ep)
c    & -Csum122(ep)
c    & -Csum1222(ep)
c    & -Csum1122(ep)
c     enddo

c     call pvBackSubst(G, 3, perm, in)

c     do ep=-2,0
c     Dv(N+dd133333,ep)=in(1,ep)
c     Dv(N+dd233333,ep)=in(2,ep)
c     Dv(N+dd333333,ep)=in(3,ep)


C   Dv(N+dd111333,Dv(N+dd112333,Dv(N+dd113333,pplll)
c     in(1,ep)=
c    & +f1*Dv(N+dd11333,ep)
c    & -4*Dv(N+dd001333,ep)
c    & -Csum222(ep)
c    & -Csum2222(ep)
c    & -Csum1222(ep)
c     in(2,ep)=
c    & +f2*Dv(N+dd11333,ep)
c    & +Cv(cc11222+C124,ep)
c    & -Csum222(ep)
c    & -Csum2222(ep)
c    & -Csum1222(ep)
c     in(3,ep)=
c    & +f3*Dv(N+dd11333,ep)
c    & -6*Dv(N+dd001133,ep)
c    & -Csum222(ep)
c    & -Csum2222(ep)
c    & -Csum1222(ep)
c     enddo

c     call pvBackSubst(G, 3, perm, in)

c     do ep=-2,0
c     Dv(N+dd133333,ep)=in(1,ep)
c     Dv(N+dd233333,ep)=in(2,ep)
c     Dv(N+dd333333,ep)=in(3,ep)

C   Dv(N+dd112222,Dv(N+dd122222,Dv(N+dd122223,pkkkk)
c     in(1,ep)=
c    & +f1*Dv(N+dd12222,ep)
c    & -2*Dv(N+dd002222,ep)
c    & +Csum1111(ep)

c     in(2,ep)=
c    & +f2*Dv(N+dd12222,ep)
c    & -8*Dv(N+dd001222,ep)
c    & +Csum1111(ep)
c     in(3,ep)=
c    & +f3*Dv(N+dd12222,ep)
c    & +Csum1111(ep)
c    & +Cv(cc12222+C123,ep)
c     enddo

c     call pvBackSubst(G, 3, perm, in)

c     do ep=-2,0
c     Dv(N+dd133333,ep)=in(1,ep)
c     Dv(N+dd233333,ep)=in(2,ep)
c     Dv(N+dd333333,ep)=in(3,ep)




C   Dv(N+dd112223,Dv(N+dd122223,Dv(N+dd122233,pkkkl)
c     in(1,ep)=
c    & +f1*Dv(N+dd12223,ep)
c    & -2*Dv(N+dd002223,ep)
c    & +Csum1112(ep)
c     in(2,ep)=
c    & +f2*Dv(N+dd12223,ep)
c    & -6*Dv(N+dd001223,ep)
c    & +Csum1112(ep)
c     in(3,ep)=
c    & +f3*Dv(N+dd12223,ep)
c    & -2*Dv(N+dd001222,ep)
c    & +Csum1112(ep)
c     enddo

c     call pvBackSubst(G, 3, perm, in)

c     do ep=-2,0
c     Dv(N+dd133333,ep)=in(1,ep)
c     Dv(N+dd233333,ep)=in(2,ep)
c     Dv(N+dd333333,ep)=in(3,ep)


C   Dv(N+dd112333,Dv(N+dd122333,Dv(N+dd123333,pklll)
c     in(1,ep)=
c    & +f1*Dv(N+dd12333,ep)
c    & -2*Dv(N+dd002333,ep)
c    & +Csum1222(ep)
c     in(2,ep)=
c    & +f2*Dv(N+dd12333,ep)
c    & -2*Dv(N+dd001333,ep)
c    & +Csum1222(ep)
c     in(3,ep)=
c    & +f3*Dv(N+dd12333,ep)
c    & -6*Dv(N+dd001233,ep)
c    & +Csum1222(ep)
c     enddo

c     call pvBackSubst(G, 3, perm, in)

c     do ep=-2,0
c     Dv(N+dd133333,ep)=in(1,ep)
c     Dv(N+dd233333,ep)=in(2,ep)
c     Dv(N+dd333333,ep)=in(3,ep)

C   Dv(N+dd113333,Dv(N+dd123333,Dv(N+dd133333,pllll)
c     in(1,ep)=
c    & +f1*Dv(N+dd13333,ep)
c    & -2*Dv(N+dd003333,ep)
c    & +Csum2222(ep)

c     in(2,ep)=
c    & +f2*Dv(N+dd13333,ep)
c    & +Cv(cc12222+C124,ep)
c    & +Csum2222(ep)
c     in(3,ep)=
c    & +f3*Dv(N+dd13333,ep)
c    & -8*Dv(N+dd001333,ep)
c    & +Csum2222(ep)
c     enddo

c     call pvBackSubst(G, 3, perm, in)

c     do ep=-2,0
c     Dv(N+dd133333,ep)=in(1,ep)
c     Dv(N+dd233333,ep)=in(2,ep)
c     Dv(N+dd333333,ep)=in(3,ep)


C   Dv(N+dd122222,Dv(N+dd222222,Dv(N+dd222223,kkkkk)
c     in(1,ep)=
c    & +f1*Dv(N+dd22222,ep)
c    & +Cv(cc11111+C134,ep)
c    & -Cv(cc11111+C234,ep)
 
c     in(2,ep)=
c    & +f2*Dv(N+dd22222,ep)
c    & -10* Dv(N+dd002222,ep)
c    & -Cv(cc11111+C234,ep)
c     in(3,ep)=
c    & +f3*Dv(N+dd22222,ep)
c    & -Cv(cc11111+C234,ep)
c    & +Cv(cc22222+C123,ep)
c     enddo

c     call pvBackSubst(G, 3, perm, in)

c     do ep=-2,0
c     Dv(N+dd133333,ep)=in(1,ep)
c     Dv(N+dd233333,ep)=in(2,ep)
c     Dv(N+dd333333,ep)=in(3,ep)


C   Dv(N+dd122223,Dv(N+dd222223,Dv(N+dd222233,kkkkl)
c     in(1,ep)=
c    & +f1*Dv(N+dd22223,ep)
c    & +Cv(cc11112+C134,ep)
c    & -Cv(cc11112+C234,ep)
c     in(2,ep)=
c    & +f2*Dv(N+dd22223,ep)
c    & -8*Dv(N+dd002223,ep)
c    & -Cv(cc11112+C234,ep)
c     in(3,ep)=
c    & +f3*Dv(N+dd22223,ep)
c    & -2*Dv(N+dd002222,ep)
c    & -Cv(cc11112+C234,ep)
c      enddo

c     call pvBackSubst(G, 3, perm, in)

c     do ep=-2,0
c     Dv(N+dd133333,ep)=in(1,ep)
c     Dv(N+dd233333,ep)=in(2,ep)
c     Dv(N+dd333333,ep)=in(3,ep)


C   Dv(N+dd122233,Dv(N+dd222233,Dv(N+dd222333,kkkll)
c     in(1,ep)=
c    & +f1*Dv(N+dd22233,ep)
c    & +Cv(cc11122+C134,ep)
c    & -Cv(cc11122+C234,ep)
c     in(2,ep)=
c    & +f2*Dv(N+dd22233,ep)
c    & -6*Dv(N+dd002233,ep)
c    & -Cv(cc11122+C234,ep)
c     in(3,ep)=
c    &+f3*Dv(N+dd22233,ep)
c    & -4*Dv(N+dd002223,ep)
c    & -Cv(cc11122+C234,ep)
c     enddo

c     call pvBackSubst(G, 3, perm, in)

c     do ep=-2,0
c     Dv(N+dd133333,ep)=in(1,ep)
c     Dv(N+dd233333,ep)=in(2,ep)
c     Dv(N+dd333333,ep)=in(3,ep)

C   Dv(N+dd122333,Dv(N+dd222333,Dv(N+dd223333,kklll)
c     in(1,ep)=
c    & +f1*Dv(N+dd22333,ep)
c    & +Cv(cc11222+C134,ep)
c    & -Cv(cc11222+C234,ep)
c     in(2,ep)=
c    & +f2*Dv(N+dd22333,ep)
c    & -4*Dv(N+dd002333,ep)
c    & -Cv(cc11222+C234,ep)
c     in(3,ep)=
c    & +f3*Dv(N+dd22333,ep)
c    & -6*Dv(N+dd002233,ep)
c    & -Cv(cc11222+C234,ep)
c     enddo

c     call pvBackSubst(G, 3, perm, in)

c     do ep=-2,0
c     Dv(N+dd133333,ep)=in(1,ep)
c     Dv(N+dd233333,ep)=in(2,ep)
c     Dv(N+dd333333,ep)=in(3,ep)


C   Dv(N+dd123333,Dv(N+dd223333,Dv(N+dd233333,kllll)
c     in(1,ep)=
c    & +f1*Dv(N+dd23333,ep)
c    & +Cv(cc12222+C134,ep)
c    & -Cv(cc12222+C234,ep)
c     in(2,ep)=
c    & +f2*Dv(N+dd23333,ep)
c    & -2*Dv(N+dd003333,ep)
c    & -Cv(cc12222+C234,ep)
c     in(3,ep)=
c    & +f3*Dv(N+dd23333,ep)
c    & -8*Dv(N+dd002333,ep)
c    & -Cv(cc12222+C234,ep)
c     enddo

c     call pvBackSubst(G, 3, perm, in)

c     do ep=-2,0
c     Dv(N+dd133333,ep)=in(1,ep)
c     Dv(N+dd233333,ep)=in(2,ep)
c     Dv(N+dd333333,ep)=in(3,ep)


C   Dv(N+dd133333,Dv(N+dd233333,Dv(N+dd333333,lllll)
c     in(1,ep)=
c    & +f1*Dv(N+dd33333,ep)
c    & +Cv(cc22222+C134,ep)
c    & -Cv(cc22222+C234,ep)
c     in(2,ep)=
c    & +f2*Dv(N+dd33333,ep)
c    & +Cv(cc22222+C124,ep)
c    & -Cv(cc22222+C234,ep)
c     in(3,ep)=
c    & +f3*Dv(N+dd33333,ep)
c    & -10* Dv(N+dd003333,ep)
c    & -Cv(cc22222+C234,ep)
c     enddo

c     call pvBackSubst(G, 3, perm, in)

c     do ep=-2,0
c     Dv(N+dd133333,ep)=in(1,ep)
c     Dv(N+dd233333,ep)=in(2,ep)
c     Dv(N+dd333333,ep)=in(3,ep)




C   Dv(N+dd000011,Dv(N+dd000012,Dv(N+dd000013,zzzzp)
c     in(1,ep)=
c    & +f1*Dv(N+dd00001,ep)
c    & -2*Dv(N+dd000000,ep)
c    & +Csum0000(ep)
c     in(2,ep)=
c    & +f2*Dv(N+dd00001,ep)
c    & +Cv(cc00001+C124,ep)
c    & +Csum0000(ep)
c     in(3,ep)=
c    & +f3*Dv(N+dd00001,ep)
c    & +Cv(cc00001+C123,ep)
c    & +Csum0000(ep)
c     enddo

c     call pvBackSubst(G, 3, perm, in)

c     do ep=-2,0
c     Dv(N+dd133333,ep)=in(1,ep)
c     Dv(N+dd233333,ep)=in(2,ep)
c     Dv(N+dd333333,ep)=in(3,ep)

C  Id,Csum0000(P?,K?,L?,m1?,m2?,m3?)=
Cc      & +Cv(cc00001+C234,ep)
Cc      & +Cv(cc00002+C234,ep)
Cc      & +Cv(cc0000+C234,ep)

C   Dv(N+dd000012,Dv(N+dd000022,Dv(N+dd000023,zzzzk)
c     in(1,ep)=
c    & +f1*Dv(N+dd00002,ep)
c    & +Cv(cc00001+C134,ep)
c    & -Cv(cc00001+C234,ep)
c     in(2,ep)=
c    & +f2*Dv(N+dd00002,ep)
c    & -2*Dv(N+dd000000,ep)
c    & -Cv(cc00001+C234,ep)
c     in(3,ep)=
c    & +f3*Dv(N+dd00002,ep)
c    & -Cv(cc00001+C234,ep)
c    & +Cv(cc00002+C123,ep)
c     enddo

c     call pvBackSubst(G, 3, perm, in)

c     do ep=-2,0
c     Dv(N+dd133333,ep)=in(1,ep)
c     Dv(N+dd233333,ep)=in(2,ep)
c     Dv(N+dd333333,ep)=in(3,ep)


C   Dv(N+dd000013,Dv(N+dd000023,Dv(N+dd000033,zzzzl)
c     in(1,ep)=
c    & +f1*Dv(N+dd00003,ep)
c    & +Cv(cc00002+C134,ep)
c    & -Cv(cc00002+C234,ep)
c     in(2,ep)=
c    & +f2*Dv(N+dd00003,ep)
c    & +Cv(cc00002+C124,ep)
c    & -Cv(cc00002+C234,ep)
c     in(3,ep)=
c    & +f3*Dv(N+dd00003,ep)
c    & -2*Dv(N+dd000000,ep)
c    & -Cv(cc00002+C234,ep)
c     enddo

c     call pvBackSubst(G, 3, perm, in)

c     do ep=-2,0
c     Dv(N+dd133333,ep)=in(1,ep)
c     Dv(N+dd233333,ep)=in(2,ep)
c     Dv(N+dd333333,ep)=in(3,ep)



C   Dv(N+dd001111,Dv(N+dd001112,Dv(N+dd001113,zzppp)
c     in(1,ep)=
c    & +f1*Dv(N+dd00111,ep)
c    & -6*Dv(N+dd000011,ep)
c    & +Csum0011(ep)
c    & +Csum0022(ep)
c    & +2*Csum0012(ep)
c    & +2*Csum002(ep)
c    & +2*Csum001(ep)
c    & +Csum00(ep)
c     in(2,ep)=
c    &+f2*Dv(N+dd00111,ep)
c    & +Cv(cc00111+C124,ep)
c    & +Csum0011(ep)
c    & +Csum0022(ep)
c    & +2*Csum0012(ep)
c    & +2*Csum002(ep)
c    & +2*Csum001(ep)
c    & +Csum00(ep)
c     in(3,ep)=
c    & +f3*Dv(N+dd00111,ep)
c    & +Cv(cc00111+C123,ep)
c    & +Csum0011(ep)
c    & +Csum0022(ep)
c    & +2*Csum0012(ep)
c    & +2*Csum002(ep)
c    & +2*Csum001(ep)
c    & +Csum00(ep)
c     enddo

c     call pvBackSubst(G, 3, perm, in)

c     do ep=-2,0
c     Dv(N+dd133333,ep)=in(1,ep)
c     Dv(N+dd233333,ep)=in(2,ep)
c     Dv(N+dd333333,ep)=in(3,ep)


C   Dv(N+dd001112,Dv(N+dd001122,Dv(N+dd001123,zzppk)
c     in(1,ep)=
c    & +f1*Dv(N+dd00112,ep)
c    & -4*Dv(N+dd000012,ep)
c    &  -Csum001(ep)
c    &  -Csum0011(ep)
c    &  -Csum0012(ep)
c     in(2,ep)=
c    & +f2*Dv(N+dd00112,ep)
c    & -2*Dv(N+dd000011,ep)
c    &  -Csum001(ep)
c    &  -Csum0011(ep)
c    &  -Csum0012(ep)
c     in(3,ep)=
c    & +f3*Dv(N+dd00112,ep)
c    & +Cv(cc00112+C123,ep)
c    &  -Csum001(ep)
c    &  -Csum0011(ep)
c    &  -Csum0012(ep)
c     enddo

c     call pvBackSubst(G, 3, perm, in)

c     do ep=-2,0
c     Dv(N+dd133333,ep)=in(1,ep)
c     Dv(N+dd233333,ep)=in(2,ep)
c     Dv(N+dd333333,ep)=in(3,ep)

C   Dv(N+dd001113,Dv(N+dd001123,Dv(N+dd001133,zzppl)
c     in(1,ep)=
c    & +f1*Dv(N+dd00113,ep)
c    & -4*Dv(N+dd000013,ep)
c    &  -Csum0022(ep)
c    &  -Csum0012(ep)
c    &  -Csum002(ep)
c     in(2,ep)=
c    & +f2*Dv(N+dd00113,ep)
c    & +Cv(cc00112+C124,ep)
c    &  -Csum0022(ep)
c    &  -Csum0012(ep)
c    &  -Csum002(ep)
c     in(3,ep)=
c    & +f3*Dv(N+dd00113,ep)
c    & -2*Dv(N+dd000011,ep)
c    &  -Csum0022(ep)
c    &  -Csum0012(ep)
c    &  -Csum002(ep)
c     enddo

c     call pvBackSubst(G, 3, perm, in)

c     do ep=-2,0
c     Dv(N+dd133333,ep)=in(1,ep)
c     Dv(N+dd233333,ep)=in(2,ep)
c     Dv(N+dd333333,ep)=in(3,ep)

C   Dv(N+dd001122,Dv(N+dd001222,Dv(N+dd001223,zzpkk)
c     in(1,ep)=
c    & +f1*Dv(N+dd00122,ep)
c    & -2*Dv(N+dd000022,ep)
c    & +Csum0011(ep)

c     in(2,ep)=
c    & +f2*Dv(N+dd00122,ep)
c    & -4*Dv(N+dd000012,ep)
c    & +Csum0011(ep)
c     in(3,ep)=
c    & +f3*Dv(N+dd00122,ep)
c    & +Cv(cc00122+C123,ep)
c    & +Csum0011(ep)
c     enddo

c     call pvBackSubst(G, 3, perm, in)

c     do ep=-2,0
c     Dv(N+dd133333,ep)=in(1,ep)
c     Dv(N+dd233333,ep)=in(2,ep)
c     Dv(N+dd333333,ep)=in(3,ep)


C   Dv(N+dd001123,Dv(N+dd001223,Dv(N+dd001233,zzpkl)
c     in(1,ep)=
c    & +f1*Dv(N+dd00123,ep)
c    & -2*Dv(N+dd000023,ep)
c    & +Csum0012(ep)
c     in(2,ep)=
c    & +f2*Dv(N+dd00123,ep)
c    & -2*Dv(N+dd000013,ep)
c    & +Csum0012(ep)
c     in(3,ep)=
c    & +f3*Dv(N+dd00123,ep)
c    & -2*Dv(N+dd000012,ep)
c    & +Csum0012(ep)
c     enddo

c     call pvBackSubst(G, 3, perm, in)

c     do ep=-2,0
c     Dv(N+dd133333,ep)=in(1,ep)
c     Dv(N+dd233333,ep)=in(2,ep)
c     Dv(N+dd333333,ep)=in(3,ep)


C   Dv(N+dd001133,Dv(N+dd001233,Dv(N+dd001333,zzpll)
c     in(1,ep)=
c    & +f1*Dv(N+dd00133,ep)
c    & -2*Dv(N+dd000033,ep)
c    & +Csum0022(ep)
c     in(2,ep)=
c    & +f2*Dv(N+dd00133,ep)
c    & +Cv(cc00122+C124,ep)
c    & +Csum0022(ep)
c     in(3,ep)=
c    & +f3*Dv(N+dd00133,ep)
c    & -4*Dv(N+dd000013,ep)
c    & +Csum0022(ep)
c     enddo

c     call pvBackSubst(G, 3, perm, in)

c     do ep=-2,0
c     Dv(N+dd133333,ep)=in(1,ep)
c     Dv(N+dd233333,ep)=in(2,ep)
c     Dv(N+dd333333,ep)=in(3,ep)


C   Dv(N+dd001222,Dv(N+dd002222,Dv(N+dd002223,zzkkk)
c     in(1,ep)=
c    & +f1*Dv(N+dd00222,ep)
c    & +Cv(cc00111+C134,ep)
c    & -Cv(cc00111+C234,ep)
c     in(2,ep)=
c    & +f2*Dv(N+dd00222,ep)
c    & -6*Dv(N+dd000022,ep)
c    & -Cv(cc00111+C234,ep)
c     in(3,ep)=
c    & +f3*Dv(N+dd00222,ep)
c    & -Cv(cc00111+C234,ep)
c    & +Cv(cc00222+C123,ep)
c     enddo

c     call pvBackSubst(G, 3, perm, in)

c     do ep=-2,0
c     Dv(N+dd133333,ep)=in(1,ep)
c     Dv(N+dd233333,ep)=in(2,ep)
c     Dv(N+dd333333,ep)=in(3,ep)


C   Dv(N+dd001223,Dv(N+dd002223,Dv(N+dd002233,zzkkl)
c     in(1,ep)=
c    & +f1*Dv(N+dd00223,ep)
c    & +Cv(cc00112+C134,ep)
c    & -Cv(cc00112+C234,ep)
c     in(2,ep)=
c    & +f2*Dv(N+dd00223,ep)
c    & -4*Dv(N+dd000023,ep)
c    & -Cv(cc00112+C234,ep)
c     in(3,ep)=
c    & +f3*Dv(N+dd00223,ep)
c    & -2*Dv(N+dd000022,ep)
c    & -Cv(cc00112+C234,ep)
c     enddo

c     call pvBackSubst(G, 3, perm, in)

c     do ep=-2,0
c     Dv(N+dd133333,ep)=in(1,ep)
c     Dv(N+dd233333,ep)=in(2,ep)
c     Dv(N+dd333333,ep)=in(3,ep)


C   Dv(N+dd001233,Dv(N+dd002233,Dv(N+dd002333,zzkll)
c     in(1,ep)=
c    & +f1*Dv(N+dd00233,ep)
c    & +Cv(cc00122+C134,ep)
c    & -Cv(cc00122+C234,ep)
c     in(2,ep)=
c    & +f2*Dv(N+dd00233,ep)
c    & -2*Dv(N+dd000033,ep)
c    & -Cv(cc00122+C234,ep)
c     in(3,ep)=
c    & +f3*Dv(N+dd00233,ep)
c    & -4*Dv(N+dd000023,ep)
c    & -Cv(cc00122+C234,ep)
c     enddo

c     call pvBackSubst(G, 3, perm, in)

c     do ep=-2,0
c     Dv(N+dd133333,ep)=in(1,ep)
c     Dv(N+dd233333,ep)=in(2,ep)
c     Dv(N+dd333333,ep)=in(3,ep)

C   Dv(N+dd001333,Dv(N+dd002333,Dv(N+dd003333,zzlll)
c     in(1,ep)=
c    & +f1*Dv(N+dd00333,ep)
c    & +Cv(cc00222+C134,ep)
c    & -Cv(cc00222+C234,ep)
c     in(2,ep)=
c    & +f2*Dv(N+dd00333,ep)
c    & +Cv(cc00222+C124,ep)
c    & -Cv(cc00222+C234,ep)
c     in(3,ep)=
c    & +f3*Dv(N+dd00333,ep)
c    & -6*Dv(N+dd000033,ep)
c    & -Cv(cc00222+C234,ep)
c     enddo

c     call pvBackSubst(G, 3, perm, in)

c     do ep=-2,0
c     Dv(N+dd133333,ep)=in(1,ep)
c     Dv(N+dd233333,ep)=in(2,ep)
c     Dv(N+dd333333,ep)=in(3,ep)

C  Id,Csum0011(P?,K?,L?,m1?,m2?,m3?)=
Cc      & +Cv(cc0011+C234,ep)
Cc      & +Cv(cc00111+C234,ep)
Cc      & +Cv(cc00112+C234,ep)

C  Id,Csum0022(P?,K?,L?,m1?,m2?,m3?)=
Cc      & +Cv(cc0022+C234,ep)
Cc      & +Cv(cc00122+C234,ep)
Cc      & +Cv(cc00222+C234,ep)

C  Id,Csum0012(P?,K?,L?,m1?,m2?,m3?)=
Cc      & +Cv(cc0012+C234,ep)
Cc      & +Cv(cc00112+C234,ep)
Cc      & +Cv(cc00122+C234,ep)




*******SEVEN
C   Dv(N+dd0000001,Dv(N+dd0000002,Dv(N+dd0000003,zzzzzz)
      do ep=-2,0
      in(1,ep)=
     &    + f1*Dv(N+dd000000,ep)
     &    + Cv(cc000000+C134,ep)
     &    - Cv(cc000000+C234,ep)
      in(2,ep)=
     &    + f2*Dv(N+dd000000,ep)
     &    + Cv(cc000000+C124,ep)
     &    - Cv(cc000000+C234,ep)
      in(3,ep)=
     &    + f3*Dv(N+dd000000,ep)
     &    + Cv(cc000000+C123,ep)
     &    - Cv(cc000000+C234,ep)
      enddo

      call pvBackSubst(G, 3, perm, in)

      do ep=-2,0
      Dv(N+dd0000001,ep)=in(1,ep)
      Dv(N+dd0000002,ep)=in(2,ep)
      Dv(N+dd0000003,ep)=in(3,ep)
      enddo

      do ep=-2,0
C   Dv(N+dd0000111,Dv(N+dd0000112,Dv(N+dd0000113,zzzzpp)
      in(1,ep)=
     &    + f1*Dv(N+dd000011,ep)
     &    - Cv(cc000011+C234,ep)
     &    - 2*Cv(cc000012+C234,ep)
     &    - Cv(cc000022+C234,ep)
     &    - 2*Cv(cc00002+C234,ep)
     &    - 2*Cv(cc00001+C234,ep)
     &    - 4*Dv(N+dd0000001,ep)
     &    - Cv(cc0000+C234,ep)
      in(2,ep)=
     &    + f2*Dv(N+dd000011,ep)
     &    + Cv(cc000011+C124,ep)
     &    - Cv(cc000011+C234,ep)
     &    - 2*Cv(cc000012+C234,ep)
     &    - 2*Cv(cc00001+C234,ep)
     &    - Cv(cc000022+C234,ep)
     &    - 2*Cv(cc00002+C234,ep)
     &    - Cv(cc0000+C234,ep)
      in(3,ep)=
     &    + f3*Dv(N+dd000011,ep)
     &    + Cv(cc000011+C123,ep)
     &    - Cv(cc000011+C234,ep)
     &    - 2*Cv(cc000012+C234,ep)
     &    - 2*Cv(cc00001+C234,ep)
     &    - Cv(cc000022+C234,ep)
     &    - 2*Cv(cc00002+C234,ep)
     &    - Cv(cc0000+C234,ep)

      enddo

      call pvBackSubst(G, 3, perm, in)

      do ep=-2,0
      Dv(N+dd0000111,ep)=in(1,ep)
      Dv(N+dd0000112,ep)=in(2,ep)
      Dv(N+dd0000113,ep)=in(3,ep)


C   Dv(N+dd0000122,Dv(N+dd0000222,Dv(N+dd0000223,zzzzkk)
      in(1,ep)=
     &    + f1*Dv(N+dd000022,ep)
     &    + Cv(cc000011+C134,ep)
     &    - Cv(cc000011+C234,ep)
      in(2,ep)=
     &    + f2*Dv(N+dd000022,ep)
     &    - 4*Dv(N+dd0000002,ep)
     &    - Cv(cc000011+C234,ep)
      in(3,ep)=
     &    + f3*Dv(N+dd000022,ep)
     &    - Cv(cc000011+C234,ep)
     &    + Cv(cc000022+C123,ep)

      enddo

      call pvBackSubst(G, 3, perm, in)

      do ep=-2,0
      Dv(N+dd0000122,ep)=in(1,ep)
      Dv(N+dd0000222,ep)=in(2,ep)
      Dv(N+dd0000223,ep)=in(3,ep)



C   Dv(N+dd0000112,Dv(N+dd0000122,Dv(N+dd0000123,zzzzpk)

      in(1,ep)=
     &    + f1*Dv(N+dd000012,ep)
     &    - 2*Dv(N+dd0000002,ep)
     &    + Cv(cc000011+C234,ep)
     &    + Cv(cc000012+C234,ep)
     &    + Cv(cc00001+C234,ep)
      in(2,ep)=
     &    + f2*Dv(N+dd000012,ep)
     &    + Cv(cc000011+C234,ep)
     &    + Cv(cc000012+C234,ep)
     &    + Cv(cc00001+C234,ep)
     &    - 2*Dv(N+dd0000001,ep)
      in(3,ep)=
     &    + f3*Dv(N+dd000012,ep)
     &    + Cv(cc000011+C234,ep)
     &    + Cv(cc000012+C123,ep)
     &    + Cv(cc000012+C234,ep)
     &    + Cv(cc00001+C234,ep)

      enddo

      call pvBackSubst(G, 3, perm, in)

      do ep=-2,0
      Dv(N+dd0000112,ep)=in(1,ep)
      Dv(N+dd0000122,ep)=in(2,ep)
      Dv(N+dd0000123,ep)=in(3,ep)




C   Dv(N+dd0000113,Dv(N+dd0000123,Dv(N+dd0000133,zzzzpl)
      in(1,ep)=
     &    + f1*Dv(N+dd000013,ep)
     &    - 2*Dv(N+dd0000003,ep)
     &    + Cv(cc000012+C234,ep)
     &    + Cv(cc000022+C234,ep)
     &    + Cv(cc00002+C234,ep)
      in(2,ep)=
     &    + f2*Dv(N+dd000013,ep)
     &    + Cv(cc000012+C124,ep)
     &    + Cv(cc000012+C234,ep)
     &    + Cv(cc000022+C234,ep)
     &    + Cv(cc00002+C234,ep)
      in(3,ep)=
     &    + f3*Dv(N+dd000013,ep)
     &    + Cv(cc000012+C234,ep)
     &    + Cv(cc000022+C234,ep)
     &    + Cv(cc00002+C234,ep)
     &    - 2*Dv(N+dd0000001,ep)

      enddo

      call pvBackSubst(G, 3, perm, in)

      do ep=-2,0
      Dv(N+dd0000113,ep)=in(1,ep)
      Dv(N+dd0000123,ep)=in(2,ep)
      Dv(N+dd0000133,ep)=in(3,ep)


C   Dv(N+dd0000123,Dv(N+dd0000223,Dv(N+dd0000233,zzzzkl)
      in(1,ep)=
     &    + f1*Dv(N+dd000023,ep)
     &    + Cv(cc000012+C134,ep)
     &    - Cv(cc000012+C234,ep)
      in(2,ep)=
     &    + f2*Dv(N+dd000023,ep)
     &    - Cv(cc000012+C234,ep)
     &    - 2*Dv(N+dd0000003,ep)
      in(3,ep)=
     &    + f3*Dv(N+dd000023,ep)
     &    - 2*Dv(N+dd0000002,ep)
     &    - Cv(cc000012+C234,ep)

      enddo

      call pvBackSubst(G, 3, perm, in)

      do ep=-2,0
      Dv(N+dd0000123,ep)=in(1,ep)
      Dv(N+dd0000223,ep)=in(2,ep)
      Dv(N+dd0000233,ep)=in(3,ep)




C   Dv(N+dd0000133,Dv(N+dd0000233,Dv(N+dd0000333,zzzzll)

      in(1,ep)=
     &    + f1*Dv(N+dd000033,ep)
     &    + Cv(cc000022+C134,ep)
     &    - Cv(cc000022+C234,ep)
      in(2,ep)=
     &    + f2*Dv(N+dd000033,ep)
     &    + Cv(cc000022+C124,ep)
     &    - Cv(cc000022+C234,ep)
      in(3,ep)=
     &    + f3*Dv(N+dd000033,ep)
     &    - Cv(cc000022+C234,ep)
     &    - 4*Dv(N+dd0000003,ep)
      enddo

      call pvBackSubst(G, 3, perm, in)

      do ep=-2,0
      Dv(N+dd0000133,ep)=in(1,ep)
      Dv(N+dd0000233,ep)=in(2,ep)
      Dv(N+dd0000333,ep)=in(3,ep)
      enddo

C   Dv(N+dd0011111,Dv(N+dd0011112,Dv(N+dd0011113,zzpppp)
      do ep=-2,0
      in(1,ep)=
     &   - Cv(cc00+C234,ep)
     &   - 4*Cv(cc001+C234,ep)
     &   - 4*Cv(cc002+C234,ep)
     &   + f1*Dv(N+dd001111,ep)
     &   - 4*Cv(cc001112+C234,ep)
     &   - 6*Cv(cc001122+C234,ep)
     &   - 12*Cv(cc00112+C234,ep)
     &   - Cv(cc001111+C234,ep)
     &   - 4*Cv(cc00111+C234,ep)
     &   - 6*Cv(cc0011+C234,ep)
     &   - 4*Cv(cc001222+C234,ep)
     &   - 12*Cv(cc00122+C234,ep)
     &   - 12*Cv(cc0012+C234,ep)
     &   - Cv(cc002222+C234,ep)
     &   - 4*Cv(cc00222+C234,ep)
     &   - 6*Cv(cc0022+C234,ep)
     &   - 8*Dv(N+dd0000111,ep)

      in(2,ep)=
     &   - Cv(cc00+C234,ep)
     &   - 4*Cv(cc001+C234,ep)
     &   - 4*Cv(cc002+C234,ep)
     &   + f2*Dv(N+dd001111,ep)
     &   - 4*Cv(cc001112+C234,ep)
     &   - 6*Cv(cc001122+C234,ep)
     &   - 12*Cv(cc00112+C234,ep)
     &   + Cv(cc001111+C124,ep)
     &   - Cv(cc001111+C234,ep)
     &   - 4*Cv(cc00111+C234,ep)
     &   - 6*Cv(cc0011+C234,ep)
     &   - 4*Cv(cc001222+C234,ep)
     &   - 12*Cv(cc00122+C234,ep)
     &   - 12*Cv(cc0012+C234,ep)
     &   - Cv(cc002222+C234,ep)
     &   - 4*Cv(cc00222+C234,ep)
     &   - 6*Cv(cc0022+C234,ep)
      in(3,ep)=
     &   - Cv(cc00+C234,ep)
     &   - 4*Cv(cc001+C234,ep)
     &   - 4*Cv(cc002+C234,ep)
     &   + f3*Dv(N+dd001111,ep)
     &   - 4*Cv(cc001112+C234,ep)
     &   - 6*Cv(cc001122+C234,ep)
     &   - 12*Cv(cc00112+C234,ep)
     &   + Cv(cc001111+C123,ep)
     &   - Cv(cc001111+C234,ep)
     &   - 4*Cv(cc00111+C234,ep)
     &   - 6*Cv(cc0011+C234,ep)
     &   - 4*Cv(cc001222+C234,ep)
     &   - 12*Cv(cc00122+C234,ep)
     &   - 12*Cv(cc0012+C234,ep)
     &   - Cv(cc002222+C234,ep)
     &   - 4*Cv(cc00222+C234,ep)
     &   - 6*Cv(cc0022+C234,ep)

       enddo

       call pvBackSubst(G, 3, perm, in)

       do ep=-2,0
       Dv(N+dd0011111,ep)=in(1,ep)
       Dv(N+dd0011112,ep)=in(2,ep)
       Dv(N+dd0011113,ep)=in(3,ep)


C   Dv(N+dd0011112,Dv(N+dd0011122,Dv(N+dd0011123,zzpppk)
      in(1,ep)=
     &    + Cv(cc001+C234,ep)
     &    + f1*Dv(N+dd001112,ep)
     &    + 3*Cv(cc001112+C234,ep)
     &    + 3*Cv(cc001122+C234,ep)
     &    + 6*Cv(cc00112+C234,ep)
     &    + Cv(cc001111+C234,ep)
     &    + 3*Cv(cc00111+C234,ep)
     &    + 3*Cv(cc0011+C234,ep)
     &    + Cv(cc001222+C234,ep)
     &    + 3*Cv(cc00122+C234,ep)
     &    + 3*Cv(cc0012+C234,ep)
     &    - 6*Dv(N+dd0000112,ep)
      in(2,ep)=
     &    + Cv(cc001+C234,ep)
     &    + f2*Dv(N+dd001112,ep)
     &    + 3*Cv(cc001112+C234,ep)
     &    + 3*Cv(cc001122+C234,ep)
     &    + 6*Cv(cc00112+C234,ep)
     &    + Cv(cc001111+C234,ep)
     &    + 3*Cv(cc00111+C234,ep)
     &    + 3*Cv(cc0011+C234,ep)
     &    + Cv(cc001222+C234,ep)
     &    + 3*Cv(cc00122+C234,ep)
     &    + 3*Cv(cc0012+C234,ep)
     &    - 2*Dv(N+dd0000111,ep)
      in(3,ep)=
     &    + Cv(cc001+C234,ep)
     &    + f3*Dv(N+dd001112,ep)
     &    + Cv(cc001112+C123,ep)
     &    + 3*Cv(cc001112+C234,ep)
     &    + 3*Cv(cc001122+C234,ep)
     &    + 6*Cv(cc00112+C234,ep)
     &    + Cv(cc001111+C234,ep)
     &    + 3*Cv(cc00111+C234,ep)
     &    + 3*Cv(cc0011+C234,ep)
     &    + Cv(cc001222+C234,ep)
     &    + 3*Cv(cc00122+C234,ep)
     &    + 3*Cv(cc0012+C234,ep)
      enddo

      call pvBackSubst(G, 3, perm, in)

      do ep=-2,0
      Dv(N+dd0011112,ep)=in(1,ep)
      Dv(N+dd0011122,ep)=in(2,ep)
      Dv(N+dd0011123,ep)=in(3,ep)


C   Dv(N+dd0011113,Dv(N+dd0011123,Dv(N+dd0011133,zzpppl)
      in(1,ep)=
     &    + Cv(cc002+C234,ep)
     &    + f1*Dv(N+dd001113,ep)
     &    + Cv(cc001112+C234,ep)
     &    + 3*Cv(cc001122+C234,ep)
     &    + 3*Cv(cc00112+C234,ep)
     &    + 3*Cv(cc001222+C234,ep)
     &    + 6*Cv(cc00122+C234,ep)
     &    + 3*Cv(cc0012+C234,ep)
     &    + Cv(cc002222+C234,ep)
     &    + 3*Cv(cc00222+C234,ep)
     &    + 3*Cv(cc0022+C234,ep)
     &    - 6*Dv(N+dd0000113,ep)
      in(2,ep)=
     &   + Cv(cc002+C234,ep)
     &    + f2*Dv(N+dd001113,ep)
     &     + Cv(cc001112+C124,ep)
     &    + Cv(cc001112+C234,ep)
     &    + 3*Cv(cc001122+C234,ep)
     &    + 3*Cv(cc00112+C234,ep)
     &    + 3*Cv(cc001222+C234,ep)
     &    + 6*Cv(cc00122+C234,ep)
     &    + 3*Cv(cc0012+C234,ep)
     &    + Cv(cc002222+C234,ep)
     &    + 3*Cv(cc00222+C234,ep)
     &    + 3*Cv(cc0022+C234,ep)
      in(3,ep)=
     &    + Cv(cc002+C234,ep)
     &    + f3*Dv(N+dd001113,ep)
     &    + Cv(cc001112+C234,ep)
     &    + 3*Cv(cc001122+C234,ep)
     &    + 3*Cv(cc00112+C234,ep)
     &    + 3*Cv(cc001222+C234,ep)
     &    + 6*Cv(cc00122+C234,ep)
     &    + 3*Cv(cc0012+C234,ep)
     &    + Cv(cc002222+C234,ep)
     &    + 3*Cv(cc00222+C234,ep)
     &    + 3*Cv(cc0022+C234,ep)
     &    - 2*Dv(N+dd0000111,ep)

      enddo

      call pvBackSubst(G, 3, perm, in)

      do ep=-2,0
      Dv(N+dd0011113,ep)=in(1,ep)
      Dv(N+dd0011123,ep)=in(2,ep)
      Dv(N+dd0011133,ep)=in(3,ep)



C   Dv(N+dd0011122,Dv(N+dd0011222,Dv(N+dd0011223,zzppkk)
      in(1,ep)=
     &   + f1*Dv(N+dd001122,ep)
     &    - 2*Cv(cc001112+C234,ep)
     &    - Cv(cc001122+C234,ep)
     &    - 2*Cv(cc00112+C234,ep)
     &    - Cv(cc001111+C234,ep)
     &    - 2*Cv(cc00111+C234,ep)
     &    - 4*Dv(N+dd0000122,ep)
     &    - Cv(cc0011+C234,ep)
      in(2,ep)=
     &    + f2*Dv(N+dd001122,ep)
     &    - 2*Cv(cc001112+C234,ep)
     &    - Cv(cc001122+C234,ep)
     &    - 2*Cv(cc00112+C234,ep)
     &    - Cv(cc001111+C234,ep)
     &    - 2*Cv(cc00111+C234,ep)
     &    - Cv(cc0011+C234,ep)
     &    - 4*Dv(N+dd0000112,ep)
      in(3,ep)=
     &    + f3*Dv(N+dd001122,ep)
     &    - 2*Cv(cc001112+C234,ep)
     &    + Cv(cc001122+C123,ep)
     &    - Cv(cc001122+C234,ep)
     &    - 2*Cv(cc00112+C234,ep)
     &    - Cv(cc001111+C234,ep)
     &    - 2*Cv(cc00111+C234,ep)
     &    - Cv(cc0011+C234,ep)

      enddo

      call pvBackSubst(G, 3, perm, in)

      do ep=-2,0
      Dv(N+dd0011122,ep)=in(1,ep)
      Dv(N+dd0011222,ep)=in(2,ep)
      Dv(N+dd0011223,ep)=in(3,ep)



C   Dv(N+dd0011123,Dv(N+dd0011223,Dv(N+dd0011233,zzppkl)
      in(1,ep)=
     &    + f1*Dv(N+dd001123,ep)
     &    - Cv(cc001112+C234,ep)
     &    - 2*Cv(cc001122+C234,ep)
     &    - 2*Cv(cc00112+C234,ep)
     &    - Cv(cc001222+C234,ep)
     &    - 2*Cv(cc00122+C234,ep)
     &    - 4*Dv(N+dd0000123,ep)
     &    - Cv(cc0012+C234,ep)
      in(2,ep)=
     &    + f2*Dv(N+dd001123,ep)
     &    - Cv(cc001112+C234,ep)
     &    - 2*Cv(cc001122+C234,ep)
     &    - 2*Cv(cc00112+C234,ep)
     &    - Cv(cc001222+C234,ep)
     &    - 2*Cv(cc00122+C234,ep)
     &    - Cv(cc0012+C234,ep)
     &    - 2*Dv(N+dd0000113,ep)
      in(3,ep)=
     &    + f3*Dv(N+dd001123,ep)
     &    - Cv(cc001112+C234,ep)
     &    - 2*Cv(cc001122+C234,ep)
     &    - 2*Cv(cc00112+C234,ep)
     &    - Cv(cc001222+C234,ep)
     &    - 2*Cv(cc00122+C234,ep)
     &    - Cv(cc0012+C234,ep)
     &    - 2*Dv(N+dd0000112,ep)
      enddo

      call pvBackSubst(G, 3, perm, in)

      do ep=-2,0
      Dv(N+dd0011123,ep)=in(1,ep)
      Dv(N+dd0011223,ep)=in(2,ep)
      Dv(N+dd0011233,ep)=in(3,ep)


 
C   Dv(N+dd0011133,Dv(N+dd0011233,Dv(N+dd0011333,zzppll)
      in(1,ep)=
     &    + f1*Dv(N+dd001133,ep)
     &    - Cv(cc001122+C234,ep)
     &    - 2*Cv(cc001222+C234,ep)
     &    - 2*Cv(cc00122+C234,ep)
     &    - Cv(cc002222+C234,ep)
     &    - 2*Cv(cc00222+C234,ep)
     &    - 4*Dv(N+dd0000133,ep)
     &    - Cv(cc0022+C234,ep)
      in(2,ep)=
     &    + f2*Dv(N+dd001133,ep)
     &    + Cv(cc001122+C124,ep)
     &    - Cv(cc001122+C234,ep)
     &    - 2*Cv(cc001222+C234,ep)
     &    - 2*Cv(cc00122+C234,ep)
     &    - Cv(cc002222+C234,ep)
     &    - 2*Cv(cc00222+C234,ep)
     &    - Cv(cc0022+C234,ep)
      in(3,ep)=
     &    + f3*Dv(N+dd001133,ep)
     &    - Cv(cc001122+C234,ep)
     &    - 2*Cv(cc001222+C234,ep)
     &    - 2*Cv(cc00122+C234,ep)
     &    - Cv(cc002222+C234,ep)
     &    - 2*Cv(cc00222+C234,ep)
     &    - Cv(cc0022+C234,ep)
     &    - 4*Dv(N+dd0000113,ep)

      enddo

      call pvBackSubst(G, 3, perm, in)

      do ep=-2,0
      Dv(N+dd0011133,ep)=in(1,ep)
      Dv(N+dd0011233,ep)=in(2,ep)
      Dv(N+dd0011333,ep)=in(3,ep)


C   Dv(N+dd0011222,Dv(N+dd0012222,Dv(N+dd0012223,zzpkkk)
      in(1,ep)=
     &    + f1*Dv(N+dd001222,ep)
     &    + Cv(cc001112+C234,ep)
     &    + Cv(cc001111+C234,ep)
     &    - 2*Dv(N+dd0000222,ep)
     &    + Cv(cc00111+C234,ep)
      in(2,ep)=
     &    + f2*Dv(N+dd001222,ep)
     &    + Cv(cc001112+C234,ep)
     &    + Cv(cc001111+C234,ep)
     &    + Cv(cc00111+C234,ep)
     &    - 6*Dv(N+dd0000122,ep)
      in(3,ep)=
     &    + f3*Dv(N+dd001222,ep)
     &    + Cv(cc001112+C234,ep)
     &    + Cv(cc001111+C234,ep)
     &    + Cv(cc00111+C234,ep)
     &    + Cv(cc001222+C123,ep)
      enddo

      call pvBackSubst(G, 3, perm, in)

      do ep=-2,0
      Dv(N+dd0011222,ep)=in(1,ep)
      Dv(N+dd0012222,ep)=in(2,ep)
      Dv(N+dd0012223,ep)=in(3,ep)



C   Dv(N+dd0011223,Dv(N+dd0012223,Dv(N+dd0012233,zzpkkl)
      in(1,ep)=
     &    + f1*Dv(N+dd001223,ep)
     &    + Cv(cc001112+C234,ep)
     &    + Cv(cc001122+C234,ep)
     &    - 2*Dv(N+dd0000223,ep)
     &    + Cv(cc00112+C234,ep)
      in(2,ep)=
     &    + f2*Dv(N+dd001223,ep)
     &    + Cv(cc001112+C234,ep)
     &    + Cv(cc001122+C234,ep)
     &    + Cv(cc00112+C234,ep)
     &    - 4*Dv(N+dd0000123,ep)
      in(3,ep)=
     &    + f3*Dv(N+dd001223,ep)
     &    + Cv(cc001112+C234,ep)
     &    + Cv(cc001122+C234,ep)
     &    + Cv(cc00112+C234,ep)
     &    - 2*Dv(N+dd0000122,ep)
      enddo

      call pvBackSubst(G, 3, perm, in)

      do ep=-2,0
      Dv(N+dd0011223,ep)=in(1,ep)
      Dv(N+dd0012223,ep)=in(2,ep)
      Dv(N+dd0012233,ep)=in(3,ep)


C   Dv(N+dd0011233,Dv(N+dd0012233,Dv(N+dd0012333,zzpkll)
      in(1,ep)=
     &    + f1*Dv(N+dd001233,ep)
     &    + Cv(cc001222+C234,ep)
     &    + Cv(cc001122+C234,ep)
     &    - 2*Dv(N+dd0000233,ep)
     &    + Cv(cc00122+C234,ep)
      in(2,ep)=
     &    + f2*Dv(N+dd001233,ep)
     &    + Cv(cc001122+C234,ep)
     &    + Cv(cc001222+C234,ep)
     &    + Cv(cc00122+C234,ep)
     &    - 2*Dv(N+dd0000133,ep)
      in(3,ep)=
     &    + f3*Dv(N+dd001233,ep)
     &    + Cv(cc001122+C234,ep)
     &    + Cv(cc001222+C234,ep)
     &    + Cv(cc00122+C234,ep)
     &    - 4*Dv(N+dd0000123,ep)

      enddo

      call pvBackSubst(G, 3, perm, in)

      do ep=-2,0
      Dv(N+dd0011233,ep)=in(1,ep)
      Dv(N+dd0012233,ep)=in(2,ep)
      Dv(N+dd0012333,ep)=in(3,ep)



C   Dv(N+dd0011333,Dv(N+dd0012333,Dv(N+dd0013333,zzplll)
      in(1,ep)=
     &    + f1*Dv(N+dd001333,ep)
     &    + Cv(cc001222+C234,ep)
     &    + Cv(cc002222+C234,ep)
     &    - 2*Dv(N+dd0000333,ep)
     &    + Cv(cc00222+C234,ep)
      in(2,ep)=
     &    + f2*Dv(N+dd001333,ep)
     &    + Cv(cc001222+C124,ep)
     &    + Cv(cc001222+C234,ep)
     &    + Cv(cc002222+C234,ep)
     &    + Cv(cc00222+C234,ep)
      in(3,ep)=
     &    + f3*Dv(N+dd001333,ep)
     &    + Cv(cc001222+C234,ep)
     &    + Cv(cc002222+C234,ep)
     &    + Cv(cc00222+C234,ep)
     &    - 6*Dv(N+dd0000133,ep)
      enddo

      call pvBackSubst(G, 3, perm, in)

      do ep=-2,0
      Dv(N+dd0011333,ep)=in(1,ep)
      Dv(N+dd0012333,ep)=in(2,ep)
      Dv(N+dd0013333,ep)=in(3,ep)



C   Dv(N+dd0012222,Dv(N+dd0022222,Dv(N+dd0022223,zzkkkk)
      in(1,ep)=
     &    + f1*Dv(N+dd002222,ep)
     &    + Cv(cc001111+C134,ep)
     &    - Cv(cc001111+C234,ep)
      in(2,ep)=
     &    + f2*Dv(N+dd002222,ep)
     &    - Cv(cc001111+C234,ep)
     &    - 8*Dv(N+dd0000222,ep)
      in(3,ep)=
     &    + f3*Dv(N+dd002222,ep)
     &    - Cv(cc001111+C234,ep)
     &    + Cv(cc002222+C123,ep)
      enddo

      call pvBackSubst(G, 3, perm, in)

      do ep=-2,0
      Dv(N+dd0012222,ep)=in(1,ep)
      Dv(N+dd0022222,ep)=in(2,ep)
      Dv(N+dd0022223,ep)=in(3,ep)




C   Dv(N+dd0012223,Dv(N+dd0022223,Dv(N+dd0022233,zzkkkl)
      in(1,ep)=
     &    + f1*Dv(N+dd002223,ep)
     &    + Cv(cc001112+C134,ep)
     &    - Cv(cc001112+C234,ep)
      in(2,ep)=
     &    + f2*Dv(N+dd002223,ep)
     &    - Cv(cc001112+C234,ep)
     &    - 6*Dv(N+dd0000223,ep)
      in(3,ep)=
     &    + f3*Dv(N+dd002223,ep)
     &    - Cv(cc001112+C234,ep)
     &    - 2*Dv(N+dd0000222,ep)

      enddo

      call pvBackSubst(G, 3, perm, in)

      do ep=-2,0
      Dv(N+dd0012223,ep)=in(1,ep)
      Dv(N+dd0022223,ep)=in(2,ep)
      Dv(N+dd0022233,ep)=in(3,ep)





C   Dv(N+dd0012233,Dv(N+dd0022233,Dv(N+dd0022333,zzkkll)
      in(1,ep)=
     &    + f1*Dv(N+dd002233,ep)
     &    + Cv(cc001122+C134,ep)
     &    - Cv(cc001122+C234,ep)
      in(2,ep)=
     &    + f2*Dv(N+dd002233,ep)
     &    - Cv(cc001122+C234,ep)
     &    - 4*Dv(N+dd0000233,ep)
      in(3,ep)=
     &    + f3*Dv(N+dd002233,ep)
     &    - Cv(cc001122+C234,ep)
     &    - 4*Dv(N+dd0000223,ep)
      enddo

      call pvBackSubst(G, 3, perm, in)

      do ep=-2,0
      Dv(N+dd0012233,ep)=in(1,ep)
      Dv(N+dd0022233,ep)=in(2,ep)
      Dv(N+dd0022333,ep)=in(3,ep)



C   Dv(N+dd0012333,Dv(N+dd0022333,Dv(N+dd0023333,zzklll)
      in(1,ep)=
     &    + f1*Dv(N+dd002333,ep)
     &    + Cv(cc001222+C134,ep)
     &    - Cv(cc001222+C234,ep)
      in(2,ep)=
     &    + f2*Dv(N+dd002333,ep)
     &    - Cv(cc001222+C234,ep)
     &    - 2*Dv(N+dd0000333,ep)
      in(3,ep)=
     &    + f3*Dv(N+dd002333,ep)
     &    - Cv(cc001222+C234,ep)
     &    - 6*Dv(N+dd0000233,ep)
      enddo

      call pvBackSubst(G, 3, perm, in)

      do ep=-2,0
      Dv(N+dd0012333,ep)=in(1,ep)
      Dv(N+dd0022333,ep)=in(2,ep)
      Dv(N+dd0023333,ep)=in(3,ep)


C   Dv(N+dd0013333,Dv(N+dd0023333,Dv(N+dd0033333,zzllll)
      in(1,ep)=
     &    + f1*Dv(N+dd003333,ep)
     &    + Cv(cc002222+C134,ep)
     &    - Cv(cc002222+C234,ep)
      in(2,ep)=
     &    + f2*Dv(N+dd003333,ep)
     &    + Cv(cc002222+C124,ep)
     &    - Cv(cc002222+C234,ep)
      in(3,ep)=
     &    + f3*Dv(N+dd003333,ep)
     &    - Cv(cc002222+C234,ep)
     &    - 8*Dv(N+dd0000333,ep)

      enddo

      call pvBackSubst(G, 3, perm, in)

      do ep=-2,0
      Dv(N+dd0013333,ep)=in(1,ep)
      Dv(N+dd0023333,ep)=in(2,ep)
      Dv(N+dd0033333,ep)=in(3,ep)

      enddo

c      do ep=-2,0
c      do j=dd0000001,dd0033333
c      write(66,*) 'D:mx',j,ep,Dv(N+j,ep)
c      enddo
c      enddo
c      pause


C   Dv(N+dd1111111,Dv(N+dd1111112,Dv(N+dd1111113,pppppp)
c      in(1,ep)=
c     &    - Cv(cc0+C234,ep)
c     &    - 6*Cv(cc1+C234,ep)
c     &    - 6*Cv(cc2+C234,ep)
c     &    - 15*Cv(cc11+C234,ep)
c     &    - 30*Cv(cc12+C234,ep)
c     &    - 15*Cv(cc22+C234,ep)
c     &    - 20*Cv(cc111+C234,ep)
c     &    - 60*Cv(cc112+C234,ep)
c     &    - 60*Cv(cc122+C234,ep)
c     &    - 20*Cv(cc222+C234,ep)
c     &    - 15*Cv(cc1111+C234,ep)
c     &    - 60*Cv(cc1112+C234,ep)
c     &    - 90*Cv(cc1122+C234,ep)
c     &    - 60*Cv(cc1222+C234,ep)
c     &    - 15*Cv(cc2222+C234,ep)
c     &    + f1*Dv(N+dd111111,ep)
c     &    - 12*Dv(N+dd0011111,ep)
c     &    - Cv(cc111111+C234,ep)
c     &    - 6*Cv(cc111112+C234,ep)
c     &    - 15*Cv(cc111122+C234,ep)
c     &    - 20*Cv(cc111222+C234,ep)
c     &    - 15*Cv(cc112222+C234,ep)
c     &    - 6*Cv(cc122222+C234,ep)
c     &    - Cv(cc222222+C234,ep)
c     &    - 6*Cv(cc11111+C234,ep)
c     &    - 30*Cv(cc11112+C234,ep)
c     &    - 60*Cv(cc11122+C234,ep)
c     &    - 60*Cv(cc11222+C234,ep)
c     &    - 30*Cv(cc12222+C234,ep)
c     &    - 6*Cv(cc22222+C234,ep)
c      in(2,ep)=
c     &    - Cv(cc0+C234,ep)
c     &    - 6*Cv(cc1+C234,ep)
c     &    - 6*Cv(cc2+C234,ep)
c     &    - 15*Cv(cc11+C234,ep)
c     &    - 30*Cv(cc12+C234,ep)
c     &    - 15*Cv(cc22+C234,ep)
c     &    - 20*Cv(cc111+C234,ep)
c     &    - 60*Cv(cc112+C234,ep)
c     &    - 60*Cv(cc122+C234,ep)
c     &    - 20*Cv(cc222+C234,ep)
c     &    - 15*Cv(cc1111+C234,ep)
c     &    - 60*Cv(cc1112+C234,ep)
c     &    - 90*Cv(cc1122+C234,ep)
c     &    - 60*Cv(cc1222+C234,ep)
c     &    - 15*Cv(cc2222+C234,ep)
c     &    + f2*Dv(N+dd111111,ep)
c     &    - 15*Cv(cc111122+C234,ep)
c     &    - 15*Cv(cc112222+C234,ep)
c     &    - 20*Cv(cc111222+C234,ep)
c     &    - 60*Cv(cc11222+C234,ep)
c     &    - 60*Cv(cc11122+C234,ep)
c     &    - 6*Cv(cc111112+C234,ep)
c     &    + Cv(cc111111+C124,ep)
c     &    - Cv(cc111111+C234,ep)
c     &    - 6*Cv(cc11111+C234,ep)
c     &    - 30*Cv(cc11112+C234,ep)
c     &    - 6*Cv(cc122222+C234,ep)
c     &    - 30*Cv(cc12222+C234,ep)
c     &    - Cv(cc222222+C234,ep)
c     &    - 6*Cv(cc22222+C234,ep)

c      in(3,ep)=
c     &    - Cv(cc0+C234,ep)
c     &    - 6*Cv(cc1+C234,ep)
c     &    - 6*Cv(cc2+C234,ep)
c     &    - 15*Cv(cc11+C234,ep)
c     &    - 30*Cv(cc12+C234,ep)
c     &    - 15*Cv(cc22+C234,ep)
c     &    - 20*Cv(cc111+C234,ep)
c     &    - 60*Cv(cc112+C234,ep)
c     &    - 60*Cv(cc122+C234,ep)
c     &    - 20*Cv(cc222+C234,ep)
c     &    - 15*Cv(cc1111+C234,ep)
c     &    - 60*Cv(cc1112+C234,ep)
c     &    - 90*Cv(cc1122+C234,ep)
c     &    - 60*Cv(cc1222+C234,ep)
c     &    - 15*Cv(cc2222+C234,ep)
c     &    + f3*Dv(N+dd111111,ep)
c     &    - 15*Cv(cc111122+C234,ep)
c     &    - 15*Cv(cc112222+C234,ep)
c     &    - 20*Cv(cc111222+C234,ep)
c     &    - 60*Cv(cc11222+C234,ep)
c     &    - 60*Cv(cc11122+C234,ep)
c     &    - 6*Cv(cc111112+C234,ep)
c     &    + Cv(cc111111+C123,ep)
c     &    - Cv(cc111111+C234,ep)
c     &    - 6*Cv(cc11111+C234,ep)
c     &    - 30*Cv(cc11112+C234,ep)
c     &    - 6*Cv(cc122222+C234,ep)
c     &    - 30*Cv(cc12222+C234,ep)
c     &    - Cv(cc222222+C234,ep)
c     &    - 6*Cv(cc22222+C234,ep)
c      enddo

c      call pvBackSubst(G, 3, perm, in)

c      do ep=-2,0
c      Dv(N+dd1111111,ep)=in(1,ep)
c      Dv(N+dd1111112,ep)=in(2,ep)
c      Dv(N+dd1111113,ep)=in(3,ep)

C   Dv(N+dd1111112,Dv(N+dd1111122,Dv(N+dd1111123,pppppk)
c      in(1,ep)=
c     &    + Cv(cc1+C234,ep)
c     &    + 5*Cv(cc11+C234,ep)
c     &    + 5*Cv(cc12+C234,ep)
c     &    + 10*Cv(cc111+C234,ep)
c     &    + 20*Cv(cc112+C234,ep)
c     &    + 10*Cv(cc122+C234,ep)
c     &    + 10*Cv(cc1111+C234,ep)
c     &    + 30*Cv(cc1112+C234,ep)
c     &    + 30*Cv(cc1122+C234,ep)
c     &    + 10*Cv(cc1222+C234,ep)
c     &    + f1*Dv(N+dd111112,ep)
c     &    - 10*Dv(N+dd0011112,ep)
c     &    + Cv(cc111111+C234,ep)
c     &    + 5*Cv(cc111112+C234,ep)
c     &    + 10*Cv(cc111122+C234,ep)
c     &    + 10*Cv(cc111222+C234,ep)
c     &    + 5*Cv(cc112222+C234,ep)
c     &    + Cv(cc122222+C234,ep)
c     &    + 5*Cv(cc11111+C234,ep)
c     &    + 20*Cv(cc11112+C234,ep)
c     &    + 30*Cv(cc11122+C234,ep)
c     &    + 20*Cv(cc11222+C234,ep)
c     &    + 5*Cv(cc12222+C234,ep)
c      in(2,ep)=
c     &    + Cv(cc1+C234,ep)
c     &    + 5*Cv(cc11+C234,ep)
c     &    + 5*Cv(cc12+C234,ep)
c     &    + 10*Cv(cc111+C234,ep)
c     &    + 20*Cv(cc112+C234,ep)
c     &    + 10*Cv(cc122+C234,ep)
c     &    + 10*Cv(cc1111+C234,ep)
c     &    + 30*Cv(cc1112+C234,ep)
c     &    + 30*Cv(cc1122+C234,ep)
c     &    + 10*Cv(cc1222+C234,ep)
c     &    + f2*Dv(N+dd111112,ep)
c     &    + 10*Cv(cc111122+C234,ep)
c     &    + 5*Cv(cc112222+C234,ep)
c     &    + 10*Cv(cc111222+C234,ep)
c     &    + 20*Cv(cc11222+C234,ep)
c     &    + 30*Cv(cc11122+C234,ep)
c     &    + 5*Cv(cc111112+C234,ep)
c     &    + Cv(cc111111+C234,ep)
c     &    + 5*Cv(cc11111+C234,ep)
c     &    + 20*Cv(cc11112+C234,ep)
c     &    + Cv(cc122222+C234,ep)
c     &    + 5*Cv(cc12222+C234,ep)
c     &    - 2*Dv(N+dd0011111,ep)

c      in(3,ep)=
c     &   + Cv(cc1+C234,ep)
c     &    + 5*Cv(cc11+C234,ep)
c     &    + 5*Cv(cc12+C234,ep)
c     &    + 10*Cv(cc111+C234,ep)
c     &    + 20*Cv(cc112+C234,ep)
c     &    + 10*Cv(cc122+C234,ep)
c     &    + 10*Cv(cc1111+C234,ep)
c     &    + 30*Cv(cc1112+C234,ep)
c     &    + 30*Cv(cc1122+C234,ep)
c     &    + 10*Cv(cc1222+C234,ep)
c     &    + f3*Dv(N+dd111112,ep)
c     &c     + 10*Cv(cc111122+C234,ep)
c     &    + 5*Cv(cc112222+C234,ep)
c     &    + 10*Cv(cc111222+C234,ep)
c     &    + 20*Cv(cc11222+C234,ep)
c     &    + 30*Cv(cc11122+C234,ep)
c     &    + Cv(cc111112+C123,ep)
c     &    + 5*Cv(cc111112+C234,ep)
c     &    + Cv(cc111111+C234,ep)
c     &    + 5*Cv(cc11111+C234,ep)
c     &    + 20*Cv(cc11112+C234,ep)
c     &    + Cv(cc122222+C234,ep)
c     &    + 5*Cv(cc12222+C234,ep)
c      enddo

c      call pvBackSubst(G, 3, perm, in)

c      do ep=-2,0
c      Dv(N+dd1111112,ep)=in(1,ep)
c      Dv(N+dd1111122,ep)=in(2,ep)
c      Dv(N+dd1111123,ep)=in(3,ep)

C   Dv(N+dd1111113,Dv(N+dd1111123,Dv(N+dd1111133,pppppl)
c      in(1,ep)=
c     &    + Cv(cc2+C234,ep)
c     &    + 5*Cv(cc12+C234,ep)
c     &    + 5*Cv(cc22+C234,ep)
c     &    + 10*Cv(cc112+C234,ep)
c     &    + 20*Cv(cc122+C234,ep)
c     &    + 10*Cv(cc222+C234,ep)
c     &    + 10*Cv(cc1112+C234,ep)
c     &    + 30*Cv(cc1122+C234,ep)
c     &    + 30*Cv(cc1222+C234,ep)
c     &    + 10*Cv(cc2222+C234,ep)
c     &    + f1*Dv(N+dd111113,ep)
c     &    - 10*Dv(N+dd0011113,ep)
c     &    + Cv(cc111112+C234,ep)
c     &    + 5*Cv(cc111122+C234,ep)
c     &    + 10*Cv(cc111222+C234,ep)
c     &    + 10*Cv(cc112222+C234,ep)
c     &    + 5*Cv(cc122222+C234,ep)
c     &    + Cv(cc222222+C234,ep)
c     &    + 5*Cv(cc11112+C234,ep)
c     &    + 20*Cv(cc11122+C234,ep)
c     &    + 30*Cv(cc11222+C234,ep)
c     &    + 20*Cv(cc12222+C234,ep)
c     &    + 5*Cv(cc22222+C234,ep)
c      in(2,ep)=
c     &    + Cv(cc2+C234,ep)
c     &    + 5*Cv(cc12+C234,ep)
c     &    + 5*Cv(cc22+C234,ep)
c     &    + 10*Cv(cc112+C234,ep)
c     &    + 20*Cv(cc122+C234,ep)
c     &    + 10*Cv(cc222+C234,ep)
c     &    + 10*Cv(cc1112+C234,ep)
c     &    + 30*Cv(cc1122+C234,ep)
c     &    + 30*Cv(cc1222+C234,ep)
c     &    + 10*Cv(cc2222+C234,ep)
c     &    + f2*Dv(N+dd111113,ep)
c     &    + 5*Cv(cc111122+C234,ep)
c     &    + 10*Cv(cc112222+C234,ep)
c     &    + 10*Cv(cc111222+C234,ep)
c     &    + 30*Cv(cc11222+C234,ep)
c     &    + 20*Cv(cc11122+C234,ep)
c     &    + Cv(cc111112+C124,ep)
c     &    + Cv(cc111112+C234,ep)
c     &    + 5*Cv(cc11112+C234,ep)
c     &    + 5*Cv(cc122222+C234,ep)
c     &    + 20*Cv(cc12222+C234,ep)
c     &    + Cv(cc222222+C234,ep)
c     &    + 5*Cv(cc22222+C234,ep)
c      in(3,ep)=
c     &    + Cv(cc2+C234,ep)
c     &    + 5*Cv(cc12+C234,ep)
c     &    + 5*Cv(cc22+C234,ep)
c     &    + 10*Cv(cc112+C234,ep)
c     &    + 20*Cv(cc122+C234,ep)
c     &    + 10*Cv(cc222+C234,ep)
c     &    + 10*Cv(cc1112+C234,ep)
c     &    + 30*Cv(cc1122+C234,ep)
c     &    + 30*Cv(cc1222+C234,ep)
c     &    + 10*Cv(cc2222+C234,ep)
c     &    + f3*Dv(N+dd111113,ep)
c     &    + 5*Cv(cc111122+C234,ep)
c     &    + 10*Cv(cc112222+C234,ep)
c     &    + 10*Cv(cc111222+C234,ep)
c     &    + 30*Cv(cc11222+C234,ep)
c     &    + 20*Cv(cc11122+C234,ep)
c     &    + Cv(cc111112+C234,ep)
c     &    + 5*Cv(cc11112+C234,ep)
c     &    + 5*Cv(cc122222+C234,ep)
c     &    + 20*Cv(cc12222+C234,ep)
c     &    + Cv(cc222222+C234,ep)
c     &    + 5*Cv(cc22222+C234,ep)
c     &    - 2*Dv(N+dd0011111,ep)
c      enddo

c      call pvBackSubst(G, 3, perm, in)

c      do ep=-2,0
c      Dv(N+dd1111113,ep)=in(1,ep)
c      Dv(N+dd1111123,ep)=in(2,ep)
c      Dv(N+dd11111333,ep)=in(3,ep)


C   Dv(N+dd1111122,Dv(N+dd1111222,Dv(N+dd1111223,ppppkk)
c      in(1,ep)=
c     &    - Cv(cc11+C234,ep)
c     &    - 4*Cv(cc111+C234,ep)
c     &    - 4*Cv(cc112+C234,ep)
c     &    - 6*Cv(cc1111+C234,ep)
c     &    - 12*Cv(cc1112+C234,ep)
c     &    - 6*Cv(cc1122+C234,ep)
c     &    + f1*Dv(N+dd111122,ep)
c     &    - 8*Dv(N+dd0011122,ep)
c     &    - Cv(cc111111+C234,ep)
c     &    - 4*Cv(cc111112+C234,ep)
c     &    - 6*Cv(cc111122+C234,ep)
c     &    - 4*Cv(cc111222+C234,ep)
c     &    - Cv(cc112222+C234,ep)
c     &    - 4*Cv(cc11111+C234,ep)
c     &    - 12*Cv(cc11112+C234,ep)
c     &    - 12*Cv(cc11122+C234,ep)
c     &    - 4*Cv(cc11222+C234,ep)
c      in(2,ep)=
c     &    - Cv(cc11+C234,ep)
c     &    - 4*Cv(cc111+C234,ep)
c     &    - 4*Cv(cc112+C234,ep)
c     &    - 6*Cv(cc1111+C234,ep)
c     &    - 12*Cv(cc1112+C234,ep)
c     &    - 6*Cv(cc1122+C234,ep)
c     &    + f2*Dv(N+dd111122,ep)
c     &    - 6*Cv(cc111122+C234,ep)
c     &    - Cv(cc112222+C234,ep)
c     &    - 4*Cv(cc111222+C234,ep)
c     &    - 4*Cv(cc11222+C234,ep)
c     &    - 12*Cv(cc11122+C234,ep)
c     &    - 4*Cv(cc111112+C234,ep)
c     &    - Cv(cc111111+C234,ep)
c     &    - 4*Cv(cc11111+C234,ep)
c     &    - 12*Cv(cc11112+C234,ep)
c     &    - 4*Dv(N+dd0011112,ep)
c      in(3,ep)=
c     &    - Cv(cc11+C234,ep)
c     &    - 4*Cv(cc111+C234,ep)
c     &    - 4*Cv(cc112+C234,ep)
c     &    - 6*Cv(cc1111+C234,ep)
c     &    - 12*Cv(cc1112+C234,ep)
c     &    - 6*Cv(cc1122+C234,ep)
c     &    + f3*Dv(N+dd111122,ep)
c     &    + Cv(cc111122+C123,ep)
c     &    - 6*Cv(cc111122+C234,ep)
c     &    - Cv(cc112222+C234,ep)
c     &    - 4*Cv(cc111222+C234,ep)
c     &    - 4*Cv(cc11222+C234,ep)
c     &    - 12*Cv(cc11122+C234,ep)
c     &    - 4*Cv(cc111112+C234,ep)
c     &    - Cv(cc111111+C234,ep)
c     &    - 4*Cv(cc11111+C234,ep)
c     &    - 12*Cv(cc11112+C234,ep)
c      enddo

c      call pvBackSubst(G, 3, perm, in)

c      do ep=-2,0
c      Dv(N+dd1111122,ep)=in(1,ep)
c      Dv(N+dd1111222,ep)=in(2,ep)
c      Dv(N+dd1111223,ep)=in(3,ep)

C   Dv(N+dd1111123,Dv(N+dd1111223,Dv(N+dd1111233,ppppkl)
c      in(1,ep)=
c     &    - Cv(cc12+C234,ep)
c     &    - 4*Cv(cc112+C234,ep)
c     &    - 4*Cv(cc122+C234,ep)
c     &    - 6*Cv(cc1112+C234,ep)
c     &    - 12*Cv(cc1122+C234,ep)
c     &    - 6*Cv(cc1222+C234,ep)
c     &    + f1*Dv(N+dd111123,ep)
c     &    - 8*Dv(N+dd0011123,ep)
c     &    - Cv(cc111112+C234,ep)
c     &    - 4*Cv(cc111122+C234,ep)
c     &    - 6*Cv(cc111222+C234,ep)
c     &    - 4*Cv(cc112222+C234,ep)
c     &    - Cv(cc122222+C234,ep)
c     &    - 4*Cv(cc11112+C234,ep)
c     &    - 12*Cv(cc11122+C234,ep)
c     &    - 12*Cv(cc11222+C234,ep)
c     &    - 4*Cv(cc12222+C234,ep)
c      in(2,ep)=
c     &    - Cv(cc12+C234,ep)
c     &    - 4*Cv(cc112+C234,ep)
c     &    - 4*Cv(cc122+C234,ep)
c     &    - 6*Cv(cc1112+C234,ep)
c     &    - 12*Cv(cc1122+C234,ep)
c     &    - 6*Cv(cc1222+C234,ep)
c     &    + f2*Dv(N+dd111123,ep)
c     &    - 4*Cv(cc111122+C234,ep)
c     &    - 4*Cv(cc112222+C234,ep)
c     &    - 6*Cv(cc111222+C234,ep)
c     &    - 12*Cv(cc11222+C234,ep)
c     &    - 12*Cv(cc11122+C234,ep)
c     &    - Cv(cc111112+C234,ep)
c     &    - 4*Cv(cc11112+C234,ep)
c     &    - Cv(cc122222+C234,ep)
c     &    - 4*Cv(cc12222+C234,ep)
c     &    - 2*Dv(N+dd0011113,ep)
c      in(3,ep)=
c     &    - Cv(cc12+C234,ep)
c     &    - 4*Cv(cc112+C234,ep)
c     &    - 4*Cv(cc122+C234,ep)
c     &    - 6*Cv(cc1112+C234,ep)
c     &    - 12*Cv(cc1122+C234,ep)
c     &    - 6*Cv(cc1222+C234,ep)
c     &    + f3*Dv(N+dd111123,ep)
c     &    - 4*Cv(cc111122+C234,ep)
c     &    - 4*Cv(cc112222+C234,ep)
c     &    - 6*Cv(cc111222+C234,ep)
c     &    - 12*Cv(cc11222+C234,ep)
c     &    - 12*Cv(cc11122+C234,ep)
c     &    - Cv(cc111112+C234,ep)
c     &    - 4*Cv(cc11112+C234,ep)
c     &    - Cv(cc122222+C234,ep)
c     &    - 4*Cv(cc12222+C234,ep)
c     &    - 2*Dv(N+dd0011112,ep)
c      enddo

c      call pvBackSubst(G, 3, perm, in)

c      do ep=-2,0
c      Dv(N+dd1111123,ep)=in(1,ep)
c      Dv(N+dd1111223,ep)=in(2,ep)
c      Dv(N+dd1111233,ep)=in(3,ep)

C   Dv(N+dd1111133,Dv(N+dd1111233,Dv(N+dd1111333,ppppll)
c      in(1,ep)=
c     &    - Cv(cc22+C234,ep)
c     &    - 4*Cv(cc122+C234,ep)
c     &    - 4*Cv(cc222+C234,ep)
c     &    - 6*Cv(cc1122+C234,ep)
c     &    - 12*Cv(cc1222+C234,ep)
c     &    - 6*Cv(cc2222+C234,ep)
c     &    + f1*Dv(N+dd111133,ep)
c     &    - 8*Dv(N+dd0011133,ep)
c     &    - Cv(cc111122+C234,ep)
c     &    - 4*Cv(cc111222+C234,ep)
c     &    - 6*Cv(cc112222+C234,ep)
c     &    - 4*Cv(cc122222+C234,ep)
c     &    - Cv(cc222222+C234,ep)
c     &    - 4*Cv(cc11122+C234,ep)
c     &    - 12*Cv(cc11222+C234,ep)
c     &    - 12*Cv(cc12222+C234,ep)
c     &    - 4*Cv(cc22222+C234,ep)
c      in(2,ep)=
c     &    - Cv(cc22+C234,ep)
c     &    - 4*Cv(cc122+C234,ep)
c     &    - 4*Cv(cc222+C234,ep)
c     &    - 6*Cv(cc1122+C234,ep)
c     &    - 12*Cv(cc1222+C234,ep)
c     &    - 6*Cv(cc2222+C234,ep)
c     &    + f2*Dv(N+dd111133,ep)
c     &    + Cv(cc111122+C124,ep)
c     &    - Cv(cc111122+C234,ep)
c     &    - 6*Cv(cc112222+C234,ep)
c     &    - 4*Cv(cc111222+C234,ep)
c     &    - 12*Cv(cc11222+C234,ep)
c     &    - 4*Cv(cc11122+C234,ep)
c     &    - 4*Cv(cc122222+C234,ep)
c     &    - 12*Cv(cc12222+C234,ep)
c     &    - Cv(cc222222+C234,ep)
c     &    - 4*Cv(cc22222+C234,ep)
c      in(3,ep)=
c     &    - Cv(cc22+C234,ep)
c     &    - 4*Cv(cc122+C234,ep)
c     &    - 4*Cv(cc222+C234,ep)
c     &    - 6*Cv(cc1122+C234,ep)
c     &    - 12*Cv(cc1222+C234,ep)
c     &    - 6*Cv(cc2222+C234,ep)
c     &    + f3*Dv(N+dd111133,ep)
c     &    - Cv(cc111122+C234,ep)
c     &    - 6*Cv(cc112222+C234,ep)
c     &    - 4*Cv(cc111222+C234,ep)
c     &    - 12*Cv(cc11222+C234,ep)
c     &    - 4*Cv(cc11122+C234,ep)
c     &    - 4*Cv(cc122222+C234,ep)
c     &    - 12*Cv(cc12222+C234,ep)
c     &    - Cv(cc222222+C234,ep)
c     &    - 4*Cv(cc22222+C234,ep)
c     &    - 4*Dv(N+dd0011113,ep)
c      enddo

c      call pvBackSubst(G, 3, perm, in)

c      do ep=-2,0
c      Dv(N+dd1111133,ep)=in(1,ep)
c      Dv(N+dd1111233,ep)=in(2,ep)
c      Dv(N+dd1111333,ep)=in(3,ep)

C   Dv(N+dd1111222,Dv(N+dd1112222,Dv(N+dd1112223,pppkkk)
c      in(1,ep)=
c     &    + Cv(cc111+C234,ep)
c     &    + 3*Cv(cc1111+C234,ep)
c     &    + 3*Cv(cc1112+C234,ep)
c     &    + f1*Dv(N+dd111222,ep)
c     &    - 6*Dv(N+dd0011222,ep)
c     &    + Cv(cc111111+C234,ep)
c     &    + 3*Cv(cc111112+C234,ep)
c     &    + 3*Cv(cc111122+C234,ep)
c     &    + Cv(cc111222+C234,ep)
c     &    + 3*Cv(cc11111+C234,ep)
c     &    + 6*Cv(cc11112+C234,ep)
c     &    + 3*Cv(cc11122+C234,ep)
c      in(2,ep)=
c     &    + Cv(cc111+C234,ep)
c     &    + 3*Cv(cc1111+C234,ep)
c     &    + 3*Cv(cc1112+C234,ep)
c     &    + f2*Dv(N+dd111222,ep)
c     &    + 3*Cv(cc111122+C234,ep)
c     &    + Cv(cc111222+C234,ep)
c     &    + 3*Cv(cc11122+C234,ep)
c     &    + 3*Cv(cc111112+C234,ep)
c     &    + Cv(cc111111+C234,ep)
c     &    + 3*Cv(cc11111+C234,ep)
c     &    + 6*Cv(cc11112+C234,ep)
c     &    - 6*Dv(N+dd0011122,ep)
c      in(3,ep)=
c     &    + Cv(cc111+C234,ep)
c     &    + 3*Cv(cc1111+C234,ep)
c     &    + 3*Cv(cc1112+C234,ep)
c     &    + f3*Dv(N+dd111222,ep)
c     &    + 3*Cv(cc111122+C234,ep)
c     &    + Cv(cc111222+C123,ep)
c     &    + Cv(cc111222+C234,ep)
c     &    + 3*Cv(cc11122+C234,ep)
c     &    + 3*Cv(cc111112+C234,ep)
c     &    + Cv(cc111111+C234,ep)
c     &    + 3*Cv(cc11111+C234,ep)
c     &    + 6*Cv(cc11112+C234,ep)
c      enddo

c      call pvBackSubst(G, 3, perm, in)

c      do ep=-2,0
c      Dv(N+dd1111222,ep)=in(1,ep)
c      Dv(N+dd1112222,ep)=in(2,ep)
c      Dv(N+dd1112223,ep)=in(3,ep)

C   Dv(N+dd1111223,Dv(N+dd1112223,Dv(N+dd1112233,pppkkl)
c      in(1,ep)=
c     &    + Cv(cc112+C234,ep)
c     &    + 3*Cv(cc1112+C234,ep)
c     &    + 3*Cv(cc1122+C234,ep)
c     &    + f1*Dv(N+dd111223,ep)
c     &    - 6*Dv(N+dd0011223,ep)
c     &    + Cv(cc111112+C234,ep)
c     &    + 3*Cv(cc111122+C234,ep)
c     &    + 3*Cv(cc111222+C234,ep)
c     &    + Cv(cc112222+C234,ep)
c     &    + 3*Cv(cc11112+C234,ep)
c     &    + 6*Cv(cc11122+C234,ep)
c     &    + 3*Cv(cc11222+C234,ep)
c      in(2,ep)=
c     &    + Cv(cc112+C234,ep)
c     &    + 3*Cv(cc1112+C234,ep)
c     &    + 3*Cv(cc1122+C234,ep)
c     &    + f2*Dv(N+dd111223,ep)
c     &    + 3*Cv(cc111122+C234,ep)
c     &    + Cv(cc112222+C234,ep)
c     &    + 3*Cv(cc111222+C234,ep)
c     &    + 3*Cv(cc11222+C234,ep)
c     &    + 6*Cv(cc11122+C234,ep)
c     &    + Cv(cc111112+C234,ep)
c     &    + 3*Cv(cc11112+C234,ep)
c     &    - 4*Dv(N+dd0011123,ep)
c      in(3,ep)=
c     &    + Cv(cc112+C234,ep)
c     &    + 3*Cv(cc1112+C234,ep)
c     &    + 3*Cv(cc1122+C234,ep)
c     &    + f3*Dv(N+dd111223,ep)
c     &    + 3*Cv(cc111122+C234,ep)
c     &    + Cv(cc112222+C234,ep)
c     &    + 3*Cv(cc111222+C234,ep)
c     &    + 3*Cv(cc11222+C234,ep)
c     &    + 6*Cv(cc11122+C234,ep)
c     &    + Cv(cc111112+C234,ep)
c     &    + 3*Cv(cc11112+C234,ep)
c     &    - 2*Dv(N+dd0011122,ep)
c      enddo

c      call pvBackSubst(G, 3, perm, in)

c      do ep=-2,0
c      Dv(N+dd1111223,ep)=in(1,ep)
c      Dv(N+dd1112223,ep)=in(2,ep)
c      Dv(N+dd1112233,ep)=in(3,ep)

C   Dv(N+dd1111233,Dv(N+dd1112233,Dv(N+dd1112333,pppkll)
c      in(1,ep)=
c     &    + Cv(cc122+C234,ep)
c     &    + 3*Cv(cc1122+C234,ep)
c     &    + 3*Cv(cc1222+C234,ep)
c     &    + f1*Dv(N+dd111233,ep)
c     &    - 6*Dv(N+dd0011233,ep)
c     &    + Cv(cc111122+C234,ep)
c     &    + 3*Cv(cc111222+C234,ep)
c     &    + 3*Cv(cc112222+C234,ep)
c     &    + Cv(cc122222+C234,ep)
c     &    + 3*Cv(cc11122+C234,ep)
c     &    + 6*Cv(cc11222+C234,ep)
c     &    + 3*Cv(cc12222+C234,ep)
c      in(2,ep)=
c     &    + Cv(cc122+C234,ep)
c     &    + 3*Cv(cc1122+C234,ep)
c     &    + 3*Cv(cc1222+C234,ep)
c     &    + f2*Dv(N+dd111233,ep)
c     &    + Cv(cc111122+C234,ep)
c     &    + 3*Cv(cc112222+C234,ep)
c     &    + 3*Cv(cc111222+C234,ep)
c     &    + 6*Cv(cc11222+C234,ep)
c     &    + 3*Cv(cc11122+C234,ep)
c     &    + Cv(cc122222+C234,ep)
c     &    + 3*Cv(cc12222+C234,ep)
c     &    - 2*Dv(N+dd0011133,ep)
c      in(3,ep)=
c     &    + Cv(cc122+C234,ep)
c     &    + 3*Cv(cc1122+C234,ep)
c     &    + 3*Cv(cc1222+C234,ep)
c     &    + f3*Dv(N+dd111233,ep)
c     &    + Cv(cc111122+C234,ep)
c     &    + 3*Cv(cc112222+C234,ep)
c     &    + 3*Cv(cc111222+C234,ep)
c     &    + 6*Cv(cc11222+C234,ep)
c     &    + 3*Cv(cc11122+C234,ep)
c     &    + Cv(cc122222+C234,ep)
c     &    + 3*Cv(cc12222+C234,ep)
c     &    - 4*Dv(N+dd0011123,ep)
c      enddo

c      call pvBackSubst(G, 3, perm, in)

c      do ep=-2,0
c      Dv(N+dd1111233,ep)=in(1,ep)
c      Dv(N+dd1112233,ep)=in(2,ep)
c      Dv(N+dd1112333,ep)=in(3,ep)

C   Dv(N+dd1111333,Dv(N+dd1112333,Dv(N+dd1113333,ppplll)
c      in(1,ep)=
c     &    + Cv(cc222+C234,ep)
c     &    + 3*Cv(cc1222+C234,ep)
c     &    + 3*Cv(cc2222+C234,ep)
c     &    + f1*Dv(N+dd111333,ep)
c     &    - 6*Dv(N+dd0011333,ep)
c     &    + Cv(cc111222+C234,ep)
c     &    + 3*Cv(cc112222+C234,ep)
c     &    + 3*Cv(cc122222+C234,ep)
c     &    + Cv(cc222222+C234,ep)
c     &    + 3*Cv(cc11222+C234,ep)
c     &    + 6*Cv(cc12222+C234,ep)
c     &    + 3*Cv(cc22222+C234,ep)
c      in(2,ep)=
c     &    + Cv(cc222+C234,ep)
c     &    + 3*Cv(cc1222+C234,ep)
c     &    + 3*Cv(cc2222+C234,ep)
c     &    + f2*Dv(N+dd111333,ep)
c     &    + 3*Cv(cc112222+C234,ep)
c     &    + Cv(cc111222+C124,ep)
c     &    + Cv(cc111222+C234,ep)
c     &    + 3*Cv(cc11222+C234,ep)
c     &    + 3*Cv(cc122222+C234,ep)
c     &    + 6*Cv(cc12222+C234,ep)
c     &    + Cv(cc222222+C234,ep)
c     &    + 3*Cv(cc22222+C234,ep)
c      in(3,ep)=
c     &    + Cv(cc222+C234,ep)
c     &    + 3*Cv(cc1222+C234,ep)
c     &    + 3*Cv(cc2222+C234,ep)
c     &    + f3*Dv(N+dd111333,ep)
c     &    + 3*Cv(cc112222+C234,ep)
c     &    + Cv(cc111222+C234,ep)
c     &    + 3*Cv(cc11222+C234,ep)
c     &    + 3*Cv(cc122222+C234,ep)
c     &    + 6*Cv(cc12222+C234,ep)
c     &    + Cv(cc222222+C234,ep)
c     &    + 3*Cv(cc22222+C234,ep)
c     &    - 6*Dv(N+dd0011133,ep)
c      enddo

c      call pvBackSubst(G, 3, perm, in)

c      do ep=-2,0
c      Dv(N+dd1111333,ep)=in(1,ep)
c      Dv(N+dd1112333,ep)=in(2,ep)
c      Dv(N+dd1113333,ep)=in(3,ep)

C   Dv(N+dd1112222,Dv(N+dd1122222,Dv(N+dd1122223,ppkkkk)
c      in(1,ep)=
c     &    - Cv(cc1111+C234,ep)
c     &    + f1*Dv(N+dd112222,ep)
c     &    - 4*Dv(N+dd0012222,ep)
c     &    - Cv(cc111111+C234,ep)
c     &    - 2*Cv(cc111112+C234,ep)
c     &    - Cv(cc111122+C234,ep)
c     &    - 2*Cv(cc11111+C234,ep)
c     &    - 2*Cv(cc11112+C234,ep)
c      in(2,ep)=
c     &    - Cv(cc1111+C234,ep)
c     &    + f2*Dv(N+dd112222,ep)
c     &    - Cv(cc111122+C234,ep)
c     &    - 2*Cv(cc111112+C234,ep)
c     &    - Cv(cc111111+C234,ep)
c     &    - 2*Cv(cc11111+C234,ep)
c     &    - 2*Cv(cc11112+C234,ep)
c     &    - 8*Dv(N+dd0011222,ep)

c      in(3,ep)=
c     &    - Cv(cc1111+C234,ep)
c     &    + f3*Dv(N+dd112222,ep)
c     &    - Cv(cc111122+C234,ep)
c     &    + Cv(cc112222+C123,ep)
c     &    - 2*Cv(cc111112+C234,ep)
c     &    - Cv(cc111111+C234,ep)
c     &    - 2*Cv(cc11111+C234,ep)
c     &    - 2*Cv(cc11112+C234,ep)

c      enddo

c      call pvBackSubst(G, 3, perm, in)

c      do ep=-2,0
c      Dv(N+dd1112222,ep)=in(1,ep)
c      Dv(N+dd1122222,ep)=in(2,ep)
c      Dv(N+dd1122223,ep)=in(3,ep)

C   Dv(N+dd1112223,Dv(N+dd1122223,Dv(N+dd1122233,ppkkkl)
c      in(1,ep)=
c     &    - Cv(cc1112+C234,ep)
c     &    + f1*Dv(N+dd112223,ep)
c     &    - 4*Dv(N+dd0012223,ep)
c     &    - Cv(cc111112+C234,ep)
c     &    - 2*Cv(cc111122+C234,ep)
c     &    - Cv(cc111222+C234,ep)
c     &    - 2*Cv(cc11112+C234,ep)
c     &    - 2*Cv(cc11122+C234,ep)
c      in(2,ep)=
c     &    - Cv(cc1112+C234,ep)
c     &    + f2*Dv(N+dd112223,ep)
c     &    - 2*Cv(cc111122+C234,ep)
c     &    - Cv(cc111222+C234,ep)
c     &    - 2*Cv(cc11122+C234,ep)
c     &    - Cv(cc111112+C234,ep)
c     &    - 2*Cv(cc11112+C234,ep)
c     &    - 6*Dv(N+dd0011223,ep)
c      in(3,ep)=
c     &    - Cv(cc1112+C234,ep)
c     &    + f3*Dv(N+dd112223,ep)
c     &    - 2*Cv(cc111122+C234,ep)
c     &    - Cv(cc111222+C234,ep)
c     &    - 2*Cv(cc11122+C234,ep)
c     &    - Cv(cc111112+C234,ep)
c     &    - 2*Cv(cc11112+C234,ep)
c     &    - 2*Dv(N+dd0011222,ep)

c      enddo

c      call pvBackSubst(G, 3, perm, in)

c      do ep=-2,0
c      Dv(N+dd1112223,ep)=in(1,ep)
c      Dv(N+dd1122223,ep)=in(2,ep)
c      Dv(N+dd1122233,ep)=in(3,ep)



C   Dv(N+dd1112233,Dv(N+dd1122233,Dv(N+dd1122333,ppkkll)
c      in(1,ep)=
c     &    - Cv(cc1122+C234,ep)
c     &    + f1*Dv(N+dd112233,ep)
c     &    - Cv(cc111122+C234,ep)
c     &    - Cv(cc112222+C234,ep)
c     &    - 2*Cv(cc111222+C234,ep)
c     &    - 2*Cv(cc11222+C234,ep)
c     &    - 2*Cv(cc11122+C234,ep)
c     &    - 4*Dv(N+dd0012233,ep)

c      in(2,ep)=
c     &    - Cv(cc1122+C234,ep)
c     &    + f2*Dv(N+dd112233,ep)
c     &    - Cv(cc111122+C234,ep)
c     &    - Cv(cc112222+C234,ep)
c     &    - 2*Cv(cc111222+C234,ep)
c     &    - 2*Cv(cc11222+C234,ep)
c     &    - 2*Cv(cc11122+C234,ep)
c     &    - 4*Dv(N+dd0011233,ep)

c      in(3,ep)=
c     &    - Cv(cc1122+C234,ep)
c     &    + f3*Dv(N+dd112233,ep)
c     &    - Cv(cc111122+C234,ep)
c     &    - Cv(cc112222+C234,ep)
c     &    - 2*Cv(cc111222+C234,ep)
c     &    - 2*Cv(cc11222+C234,ep)
c     &    - 2*Cv(cc11122+C234,ep)
c     &    - 4*Dv(N+dd0011223,ep)

c      enddo

c      call pvBackSubst(G, 3, perm, in)

c      do ep=-2,0
c      Dv(N+dd1112233,ep)=in(1,ep)
c      Dv(N+dd1122233,ep)=in(2,ep)
c      Dv(N+dd1122333,ep)=in(3,ep)


C   Dv(N+dd1112333,Dv(N+dd1122333,Dv(N+dd1123333,ppklll)
c      in(1,ep)=
c     &    - Cv(cc1222+C234,ep)
c     &    + f1*Dv(N+dd112333,ep)
c     &    - 4*Dv(N+dd0012333,ep)
c     &    - Cv(cc111222+C234,ep)
c     &    - 2*Cv(cc112222+C234,ep)
c     &    - Cv(cc122222+C234,ep)
c     &    - 2*Cv(cc11222+C234,ep)
c     &    - 2*Cv(cc12222+C234,ep)

c      in(2,ep)=
c     &    - Cv(cc1222+C234,ep)
c     &    + f2*Dv(N+dd112333,ep)
c     &    - 2*Cv(cc112222+C234,ep)
c     &    - Cv(cc111222+C234,ep)
c     &    - 2*Cv(cc11222+C234,ep)
c     &    - Cv(cc122222+C234,ep)
c     &    - 2*Cv(cc12222+C234,ep)
c     &    - 2*Dv(N+dd0011333,ep)

c      in(3,ep)=
c     &    - Cv(cc1222+C234,ep)
c     &    + f3*Dv(N+dd112333,ep)
c     &    - 2*Cv(cc112222+C234,ep)
c     &    - Cv(cc111222+C234,ep)
c     &    - 2*Cv(cc11222+C234,ep)
c     &    - Cv(cc122222+C234,ep)
c     &    - 2*Cv(cc12222+C234,ep)
c     &    - 6*Dv(N+dd0011233,ep)

c      enddo

c      call pvBackSubst(G, 3, perm, in)

c      do ep=-2,0
c      Dv(N+dd1112333,ep)=in(1,ep)
c      Dv(N+dd1122333,ep)=in(2,ep)
c      Dv(N+dd1123333,ep)=in(3,ep)





C   Dv(N+dd1113333,Dv(N+dd1123333,Dv(N+dd1133333,ppllll)
c      in(1,ep)=
c     &    - Cv(cc2222+C234,ep)
c     &    + f1*Dv(N+dd113333,ep)
c     &    - 4*Dv(N+dd0013333,ep)
c     &    - Cv(cc112222+C234,ep)
c     &    - 2*Cv(cc122222+C234,ep)
c     &    - Cv(cc222222+C234,ep)
c     &    - 2*Cv(cc12222+C234,ep)
c     &    - 2*Cv(cc22222+C234,ep)
c      in(2,ep)=
c     &    - Cv(cc2222+C234,ep)
c     &    + f2*Dv(N+dd113333,ep)
c     &    + Cv(cc112222+C124,ep)
c     &    - Cv(cc112222+C234,ep)
c     &    - 2*Cv(cc122222+C234,ep)
c     &    - 2*Cv(cc12222+C234,ep)
c     &    - Cv(cc222222+C234,ep)
c     &    - 2*Cv(cc22222+C234,ep)

c      in(3,ep)=
c     &    - Cv(cc2222+C234,ep)
c     &    + f3*Dv(N+dd113333,ep)
c     &    - Cv(cc112222+C234,ep)
c     &    - 2*Cv(cc122222+C234,ep)
c     &    - 2*Cv(cc12222+C234,ep)
c     &    - Cv(cc222222+C234,ep)
c     &    - 2*Cv(cc22222+C234,ep)
c     &    - 8*Dv(N+dd0011333,ep)

c      enddo

c      call pvBackSubst(G, 3, perm, in)

c      do ep=-2,0
c      Dv(N+dd1113333,ep)=in(1,ep)
c      Dv(N+dd1123333,ep)=in(2,ep)
c      Dv(N+dd1133333,ep)=in(3,ep)


C   Dv(N+dd1122222,Dv(N+dd1222222,Dv(N+dd1222223,pkkkkk)
c      in(1,ep)=
c     &    + f1*Dv(N+dd122222,ep)
c     &    + Cv(cc111112+C234,ep)
c     &    + Cv(cc111111+C234,ep)
c     &    + Cv(cc11111+C234,ep)
c     &    - 2*Dv(N+dd0022222,ep)
c      in(2,ep)=
c     &    + f2*Dv(N+dd122222,ep)
c     &    + Cv(cc111112+C234,ep)
c     &    + Cv(cc111111+C234,ep)
c     &    + Cv(cc11111+C234,ep)
c     &    - 10*Dv(N+dd0012222,ep)
c      in(3,ep)=
c     &    + f3*Dv(N+dd122222,ep)
c     &    + Cv(cc111112+C234,ep)
c     &    + Cv(cc111111+C234,ep)
c     &    + Cv(cc11111+C234,ep)
c     &    + Cv(cc122222+C123,ep)

c      enddo

c      call pvBackSubst(G, 3, perm, in)

c      do ep=-2,0
c      Dv(N+dd1122222,ep)=in(1,ep)
c      Dv(N+dd1222222,ep)=in(2,ep)
c      Dv(N+dd1222223,ep)=in(3,ep)




C   Dv(N+dd1122223,Dv(N+dd1222223,Dv(N+dd1222233,pkkkkl)
c      in(1,ep)=
c     & + f1*Dv(N+dd122223,ep)
c     &    + Cv(cc111122+C234,ep)
c     &    + Cv(cc111112+C234,ep)
c     &    + Cv(cc11112+C234,ep)
c     &    - 2*Dv(N+dd0022223,ep)

c      in(2,ep)=
c     &    + f2*Dv(N+dd122223,ep)
c     &    + Cv(cc111122+C234,ep)
c     &    + Cv(cc111112+C234,ep)
c     &    + Cv(cc11112+C234,ep)
c     &    - 8*Dv(N+dd0012223,ep)

c      in(3,ep)=
c     &    + f3*Dv(N+dd122223,ep)
c     &    + Cv(cc111122+C234,ep)
c     &    + Cv(cc111112+C234,ep)
c     &    + Cv(cc11112+C234,ep)
c     &    - 2*Dv(N+dd0012222,ep)

c      enddo

c      call pvBackSubst(G, 3, perm, in)

c      do ep=-2,0
c      Dv(N+dd1122223,ep)=in(1,ep)
c      Dv(N+dd1222223,ep)=in(2,ep)
c      Dv(N+dd1222233,ep)=in(3,ep)





C   Dv(N+dd1122233,Dv(N+dd1222233,Dv(N+dd1222333,pkkkll)
c      in(1,ep)=
c     &    + f1*Dv(N+dd122233,ep)
c     &    + Cv(cc111122+C234,ep)
c     &    + Cv(cc111222+C234,ep)
c     &    + Cv(cc11122+C234,ep)
c     &    - 2*Dv(N+dd0022233,ep)

c      in(2,ep)=
c     &    + f2*Dv(N+dd122233,ep)
c     &    + Cv(cc111122+C234,ep)
c     &    + Cv(cc111222+C234,ep)
c     &    + Cv(cc11122+C234,ep)
c     &    - 6*Dv(N+dd0012233,ep)
c      in(3,ep)=
c     &    + f3*Dv(N+dd122233,ep)
c     &    + Cv(cc111122+C234,ep)
c     &    + Cv(cc111222+C234,ep)
c     &    + Cv(cc11122+C234,ep)
c     &    - 4*Dv(N+dd0012223,ep)

c      enddo

c      call pvBackSubst(G, 3, perm, in)

c      do ep=-2,0
c      Dv(N+dd1122233,ep)=in(1,ep)
c      Dv(N+dd1222233,ep)=in(2,ep)
c      Dv(N+dd1222333,ep)=in(3,ep)



C   Dv(N+dd1122333,Dv(N+dd1222333,Dv(N+dd1223333,pkklll)
c      in(1,ep)=
c     &    + f1*Dv(N+dd122333,ep)
c     &    + Cv(cc112222+C234,ep)
c     &    + Cv(cc111222+C234,ep)
c     &    + Cv(cc11222+C234,ep)
c     &    - 2*Dv(N+dd0022333,ep)
c      in(2,ep)=
c     &    + f2*Dv(N+dd122333,ep)
c     &    + Cv(cc112222+C234,ep)
c     &    + Cv(cc111222+C234,ep)
c     &    + Cv(cc11222+C234,ep)
c     &    - 4*Dv(N+dd0012333,ep)

c      in(3,ep)=
c     &    + f3*Dv(N+dd122333,ep)
c     &    + Cv(cc112222+C234,ep)
c     &    + Cv(cc111222+C234,ep)
c     &    + Cv(cc11222+C234,ep)
c     &    - 6*Dv(N+dd0012233,ep)

c      enddo

c      call pvBackSubst(G, 3, perm, in)

c      do ep=-2,0
c      Dv(N+dd1122333,ep)=in(1,ep)
c      Dv(N+dd1222333,ep)=in(2,ep)
c      Dv(N+dd1223333,ep)=in(3,ep)


C   Dv(N+dd1123333,Dv(N+dd1223333,Dv(N+dd1233333,pkllll)
c      in(1,ep)=
c     &    + f1*Dv(N+dd123333,ep)
c     &    + Cv(cc112222+C234,ep)
c     &    + Cv(cc122222+C234,ep)
c     &    + Cv(cc12222+C234,ep)
c     &    - 2*Dv(N+dd0023333,ep)
c      in(2,ep)=
c     &    + f2*Dv(N+dd123333,ep)
c     &    + Cv(cc112222+C234,ep)
c     &    + Cv(cc122222+C234,ep)
c     &    + Cv(cc12222+C234,ep)
c     &    - 2*Dv(N+dd0013333,ep)
c      in(3,ep)=
c     &    + f3*Dv(N+dd123333,ep)
c     &    + Cv(cc112222+C234,ep)
c     &    + Cv(cc122222+C234,ep)
c     &    + Cv(cc12222+C234,ep)
c     &    - 8*Dv(N+dd0012333,ep)

c      enddo

c      call pvBackSubst(G, 3, perm, in)

c      do ep=-2,0
c      Dv(N+dd1123333,ep)=in(1,ep)
c      Dv(N+dd1223333,ep)=in(2,ep)
c      Dv(N+dd1233333,ep)=in(3,ep)


C   Dv(N+dd1133333,Dv(N+dd1233333,Dv(N+dd1333333,plllll)
c      in(1,ep)=
c     &    + f1*Dv(N+dd133333,ep)
c     &    + Cv(cc122222+C234,ep)
c     &    + Cv(cc222222+C234,ep)
c     &    + Cv(cc22222+C234,ep)
c     &    - 2*Dv(N+dd0033333,ep)

c      in(2,ep)=
c     &    + f2*Dv(N+dd133333,ep)
c     &    + Cv(cc122222+C124,ep)
c     &    + Cv(cc122222+C234,ep)
c     &    + Cv(cc222222+C234,ep)
c     &    + Cv(cc22222+C234,ep)

c      in(3,ep)=
c     &    + f3*Dv(N+dd133333,ep)
c     &    + Cv(cc122222+C234,ep)
c     &    + Cv(cc222222+C234,ep)
c     &    + Cv(cc22222+C234,ep)
c     &    - 10*Dv(N+dd0013333,ep)

c      enddo

c      call pvBackSubst(G, 3, perm, in)

c      do ep=-2,0
c      Dv(N+dd1133333,ep)=in(1,ep)
c      Dv(N+dd1233333,ep)=in(2,ep)
c      Dv(N+dd1333333,ep)=in(3,ep)



C   Dv(N+dd1222222,Dv(N+dd2222222,Dv(N+dd2222223,kkkkkk)
c      in(1,ep)=
c     &    + f1*Dv(N+dd222222,ep)
c     &    + Cv(cc111111+C134,ep)
c     &    - Cv(cc111111+C234,ep)
c      in(2,ep)=
c     &    + f2*Dv(N+dd222222,ep)
c     &    - Cv(cc111111+C234,ep)
c     &    - 12*Dv(N+dd0022222,ep)
c      in(3,ep)=
c     &    + f3*Dv(N+dd222222,ep)
c     &    - Cv(cc111111+C234,ep)
c     &    + Cv(cc222222+C123,ep)

c      enddo

c      call pvBackSubst(G, 3, perm, in)

c      do ep=-2,0
c      Dv(N+dd1222222,ep)=in(1,ep)
c      Dv(N+dd2222222,ep)=in(2,ep)
c      Dv(N+dd2222223,ep)=in(3,ep)


C   Dv(N+dd1222223,Dv(N+dd2222223,Dv(N+dd2222233,kkkkkl)
c      in(1,ep)=
c     &    + f1*Dv(N+dd222223,ep)
c     &    + Cv(cc111112+C134,ep)
c     &    - Cv(cc111112+C234,ep)
c      in(2,ep)=
c     &    + f2*Dv(N+dd222223,ep)
c     &    - Cv(cc111112+C234,ep)
c     &    - 10*Dv(N+dd0022223,ep)
c      in(3,ep)=
c     &    + f3*Dv(N+dd222223,ep)
c     &    - Cv(cc111112+C234,ep)
c     &    - 2*Dv(N+dd0022222,ep)

c      enddo

c      call pvBackSubst(G, 3, perm, in)

c      do ep=-2,0
c      Dv(N+dd1222223,ep)=in(1,ep)
c      Dv(N+dd2222223,ep)=in(2,ep)
c      Dv(N+dd2222233,ep)=in(3,ep)


C   Dv(N+dd1222233,Dv(N+dd2222233,Dv(N+dd2222333,kkkkll)
c      in(1,ep)=
c     &    + f1*Dv(N+dd222233,ep)
c     &    + Cv(cc111122+C134,ep)
c     &    - Cv(cc111122+C234,ep)
c      in(2,ep)=
c     &    + f2*Dv(N+dd222233,ep)
c     &    - Cv(cc111122+C234,ep)
c     &    - 8*Dv(N+dd0022233,ep)
c      in(3,ep)=
c     &    + f3*Dv(N+dd222233,ep)
c     &    - Cv(cc111122+C234,ep)
c     &    - 4*Dv(N+dd0022223,ep)

c      enddo

c      call pvBackSubst(G, 3, perm, in)

c      do ep=-2,0
c      Dv(N+dd1222233,ep)=in(1,ep)
c      Dv(N+dd2222233,ep)=in(2,ep)
c      Dv(N+dd2222333,ep)=in(3,ep)


C   Dv(N+dd1222333,Dv(N+dd2222333,Dv(N+dd2223333,kkklll)
c      in(1,ep)=
c     &    + f1*Dv(N+dd222333,ep)
c     &    + Cv(cc111222+C134,ep)
c     &    - Cv(cc111222+C234,ep)
c      in(2,ep)=
c     &    + f2*Dv(N+dd222333,ep)
c     &    - Cv(cc111222+C234,ep)
c     &    - 6*Dv(N+dd0022333,ep)
c      in(3,ep)=
c     &    + f3*Dv(N+dd222333,ep)
c     &    - Cv(cc111222+C234,ep)
c     &    - 6*Dv(N+dd0022233,ep)

c      enddo

c      call pvBackSubst(G, 3, perm, in)

c      do ep=-2,0
c      Dv(N+dd1222333,ep)=in(1,ep)
c      Dv(N+dd2222333,ep)=in(2,ep)
c      Dv(N+dd2223333,ep)=in(3,ep)


C   Dv(N+dd1223333,Dv(N+dd2223333,Dv(N+dd2233333,kkllll)
c      in(1,ep)=
c     &    + f1*Dv(N+dd223333,ep)
c     &    + Cv(cc112222+C134,ep)
c     &    - Cv(cc112222+C234,ep)
c      in(2,ep)=
c     &    + f2*Dv(N+dd223333,ep)
c     &    - Cv(cc112222+C234,ep)
c     &    - 4*Dv(N+dd0023333,ep)
c      in(3,ep)=
c     &    + f3*Dv(N+dd223333,ep)
c     &    - Cv(cc112222+C234,ep)
c     &    - 8*Dv(N+dd0022333,ep)

c      enddo

c      call pvBackSubst(G, 3, perm, in)

c      do ep=-2,0
c      Dv(N+dd1223333,ep)=in(1,ep)
c      Dv(N+dd2223333,ep)=in(2,ep)
c      Dv(N+dd2233333,ep)=in(3,ep)


C   Dv(N+dd1233333,Dv(N+dd2233333,Dv(N+dd2333333,klllll)
c      in(1,ep)=
c     &    + f1*Dv(N+dd233333,ep)
c     &    + Cv(cc122222+C134,ep)
c     &    - Cv(cc122222+C234,ep)

c      in(2,ep)=
c     &    + f2*Dv(N+dd233333,ep)
c     &    - Cv(cc122222+C234,ep)
c     &    - 2*Dv(N+dd0033333,ep)
c      in(3,ep)=
c     &    + f3*Dv(N+dd233333,ep)
c     &    - Cv(cc122222+C234,ep)
c     &    - 10*Dv(N+dd0023333,ep)

c      enddo

c      call pvBackSubst(G, 3, perm, in)

c      do ep=-2,0
c      Dv(N+dd1233333,ep)=in(1,ep)
c      Dv(N+dd2233333,ep)=in(2,ep)
c      Dv(N+dd2333333,ep)=in(3,ep)


C   Dv(N+dd1333333,Dv(N+dd2333333,Dv(N+dd3333333,llllll)
c      in(1,ep)=
c     &   + f1*Dv(N+dd333333,ep)
c     &    + Cv(cc222222+C134,ep)
c     &    - Cv(cc222222+C234,ep)
c      in(2,ep)=
c     &    + f2*Dv(N+dd333333,ep)
c     &    + Cv(cc222222+C124,ep)
c     &    - Cv(cc222222+C234,ep)
c      in(3,ep)=
c     &    + f3*Dv(N+dd333333,ep)
c     &    - Cv(cc222222+C234,ep)
c     &    - 12*Dv(N+dd0033333,ep)

c      enddo

c      call pvBackSubst(G, 3, perm, in)

c      do ep=-2,0
c      Dv(N+dd1333333,ep)=in(1,ep)
c      Dv(N+dd2333333,ep)=in(2,ep)
c      Dv(N+dd3333333,ep)=in(3,ep)


c      enddo






