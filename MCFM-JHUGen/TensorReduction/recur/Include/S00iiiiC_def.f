      do ep=-2,0
      S00iiii(z4(1,1,1,1),ep)=+2d0*(Bsum0(ep)+3*Bsum1(ep)+3*Bsum11(ep)
     .           +Bsum111(ep))+2d0*m1*Cv(cc1111+N0,ep)
      S00iiii(z4(2,2,2,2),ep)=+2d0*Bv(bb1111 + B23,ep)
     .                        +2d0*m1*Cv(cc2222+N0,ep)
      S00iiii(z4(1,1,1,2),ep)=-2d0*(Bsum1(ep)+2*Bsum11(ep)+Bsum111(ep))
     .                        +2d0*m1*Cv(cc1112+N0,ep)
      S00iiii(z4(1,1,2,2),ep)=2d0*(Bsum11(ep)+Bsum111(ep))
     .                        +2d0*m1*Cv(cc1122+N0,ep)
      S00iiii(z4(1,2,2,2),ep)=-2d0*Bsum111(ep)
     .                        +2d0*m1*Cv(cc1222+N0,ep)
      enddo
