      do ep=-2,0
      S00ii(z2(1,1),ep)=+2d0*(Csum0(ep)+Csum1(ep)+Csum2(ep))
     .                  +2d0*m1*Dv(dd11+N0,ep)
      S00ii(z2(1,2),ep)=-2d0*Csum1(ep)+2d0*m1*Dv(dd12+N0,ep)
      S00ii(z2(1,3),ep)=-2d0*Csum2(ep)+2d0*m1*Dv(dd13+N0,ep)
      S00ii(z2(2,2),ep)=+2d0*Cv(cc11+C234,ep)+2d0*m1*Dv(dd22+N0,ep)
      S00ii(z2(2,3),ep)=+2d0*Cv(cc12+C234,ep)+2d0*m1*Dv(dd23+N0,ep)
      S00ii(z2(3,3),ep)=+2d0*Cv(cc22+C234,ep)+2d0*m1*Dv(dd33+N0,ep)
      enddo
