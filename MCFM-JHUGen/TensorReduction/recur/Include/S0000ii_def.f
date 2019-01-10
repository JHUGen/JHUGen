      do ep=-2,0
      S0000ii(z2(1,1),ep)=+2d0*(Csum00(ep)+Csum001(ep)+Csum002(ep))
     .                    +2d0*m1*Dv(dd0011+N0,ep)
      S0000ii(z2(1,2),ep)=-2d0*Csum001(ep)
     .                    +2d0*m1*Dv(dd0012+N0,ep)
      S0000ii(z2(1,3),ep)=-2d0*Csum002(ep)
     .                    +2d0*m1*Dv(dd0013+N0,ep)
      S0000ii(z2(2,2),ep)=+2d0*Cv(cc0011+C234,ep)
     .                    +2d0*m1*Dv(dd0022+N0,ep)
      S0000ii(z2(2,3),ep)=+2d0*Cv(cc0012+C234,ep)
     .                    +2d0*m1*Dv(dd0023+N0,ep)
      S0000ii(z2(3,3),ep)=+2d0*Cv(cc0022+C234,ep)
     .                    +2d0*m1*Dv(dd0033+N0,ep)
      enddo
