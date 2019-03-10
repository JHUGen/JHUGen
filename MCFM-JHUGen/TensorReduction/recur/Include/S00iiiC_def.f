      do ep=-2,0
      S00iii(z3(1,1,1),ep)=-2d0*(Bsum0(ep)+2*Bsum1(ep)+Bsum11(ep))
     .                     +2d0*m1*Cv(cc111+N0,ep)
      S00iii(z3(1,1,2),ep)=+2d0*(Bsum1(ep)+Bsum11(ep))
     .                     +2d0*m1*Cv(cc112+N0,ep)
      S00iii(z3(1,2,2),ep)=-2d0*Bsum11(ep)
     .                     +2d0*m1*Cv(cc122+N0,ep)
      S00iii(z3(2,2,2),ep)=+2d0*Bv(bb111 + B23,ep)
     .                     +2d0*m1*Cv(cc222+N0,ep)
      enddo
