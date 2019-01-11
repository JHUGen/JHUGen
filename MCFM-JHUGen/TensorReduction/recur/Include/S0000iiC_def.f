      do ep=-2,0
      S0000ii(z2(1,1),ep)=+2d0*(Bsum00(ep)+Bsum001(ep))
     .                    +2d0*m1*Cv(cc0011+N0,ep)
      S0000ii(z2(1,2),ep)=-2d0*Bsum001(ep)
     .                    +2d0*m1*Cv(cc0012+N0,ep)
      S0000ii(z2(2,2),ep)=+2d0*Bv(bb0011 + B23,ep)
     .                    +2d0*m1*Cv(cc0022+N0,ep)
      enddo
