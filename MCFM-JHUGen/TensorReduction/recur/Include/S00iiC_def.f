      do ep=-2,0
      S00ii(z2(1,1),ep)=+2d0*(Bsum0(ep)+Bsum1(ep))+2d0*m1*Cv(cc11+N0,ep)
      S00ii(z2(1,2),ep)=-2d0*Bsum1(ep)+2d0*m1*Cv(cc12+N0,ep)
      S00ii(z2(2,2),ep)=+2d0*Bv(bb11 + B23,ep)+2d0*m1*Cv(cc22+N0,ep)
      enddo
