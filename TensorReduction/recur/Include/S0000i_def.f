      do ep=-2,0
      S0000i(1,ep)=-2d0*Csum00(ep)+2d0*m1*Dv(dd001+N0,ep)
      S0000i(2,ep)=+2d0*Cv(cc001+C234,ep)+2d0*m1*Dv(dd002+N0,ep)
      S0000i(3,ep)=+2d0*Cv(cc002+C234,ep)+2d0*m1*Dv(dd003+N0,ep)
      enddo
