      do ep=-2,0
      S00i(1,ep)=-2d0*Csum0(ep)+2d0*m1*Dv(dd1+N0,ep) 
      S00i(2,ep)=+2d0*Cv(cc1+C234,ep)+2d0*m1*Dv(dd2+N0,ep)
      S00i(3,ep)=+2d0*Cv(cc2+C234,ep)+2d0*m1*Dv(dd3+N0,ep)
      enddo
