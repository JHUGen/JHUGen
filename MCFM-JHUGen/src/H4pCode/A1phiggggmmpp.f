c--- Results taken from:
c---   S.~D.~Badger, E.~W.~N.~Glover and K.~Risager,
c---   %``One-loop phigggg-MHV amplitudes using the unitarity bootstrap,''
c---   JHEP {\bf 0707}, 066 (2007)
c---   [arXiv:0704.3914 [hep-ph]].

c--- with the sign of all terms proportional to Np or b0 reversed

      double complex function A1phiggggmmpp(j1,j2,j3,j4,za,zb)
      implicit none
C     arXiv:0704.3914v3, Eq.(5.13)
      integer j1,j2,j3,j4
      include 'constants.f'
      include 'zprods_decl.f'
      double complex C4mmpphat,Rhat4mmpp

c---- version using "hatting" procedure
      A1phiggggmmpp=C4mmpphat(j1,j2,j3,j4,za,zb)
     &             +Rhat4mmpp(j1,j2,j3,j4,za,zb)

c      write(6,*) 'CHECK: ',A1phiggggmmpp

c--- note, we could equivalently have written
c      A1phiggggmmpp=   C4mmpp(j1,j2,j3,j4,za,zb)
c     &               +CR4mmpp(j1,j2,j3,j4,za,zb)
c     &             +Rhat4mmpp(j1,j2,j3,j4,za,zb)
c      write(6,*) 'CHECK: ',A1phiggggmmpp

      return
      end

      double complex function Rhat4mmpp(j1,j2,j3,j4,za,zb)
      implicit none
C     arXiv:0704.3914v3, Eq.(5.12) (factor of 16 pi^2 removed)
      include 'constants.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'nflav.f'
      integer j1,j2,j3,j4
      double complex A0phiggggmmpp,zab2
      double precision Np,s3
      zab2(j1,j2,j3,j4)=+za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
      s3(j1,j2,j3)=s(j1,j2)+s(j2,j3)+s(j3,j1)
      Np=2d0*(1d0-dfloat(nflav)/xn)

c--- Note: implicit use of A0(A,...)=A0phigggg(...)-A0phiggggdagger(...)
c---       and the fact that A0phiggggdagger=(c.c. of A0phigggg)
      Rhat4mmpp=2d0
     & *(A0phiggggmmpp(j1,j2,j3,j4,za,zb)
     &  -A0phiggggmmpp(j3,j4,j1,j2,zb,za))
     & +Np/6d0*zb(j4,j3)/za(j3,j4)
     & *(-za(j2,j3)*zab2(j1,j2,j4,j3)**2
     & /(za(j3,j4)*zb(j4,j3)*zb(j3,j2)*s3(j2,j3,j4))
     &  +za(j4,j1)*zab2(j3,j1,j2,j3)
     & /(za(j3,j4)*zb(j1,j2)*zb(j3,j2))
     &  -za(j1,j4)*zab2(j2,j1,j3,j4)**2
     & /(za(j3,j4)*zb(j4,j3)*zb(j4,j1)*s3(j3,j4,j1))
     &  +za(j3,j2)*zab2(j4,j1,j2,j4)
     & /(za(j3,j4)*zb(j1,j2)*zb(j4,j1))
     &  +za(j1,j2)**2/(za(j3,j4)*zb(j4,j3))-za(j1,j2)/zb(j1,j2)
     &  -za(j1,j2)*zab2(j2,j1,j3,j4)
     & /(2d0*zb(j4,j1)*s3(j3,j4,j1))
     &  +za(j1,j2)*zab2(j1,j2,j4,j3)
     & /(2d0*zb(j3,j2)*s3(j2,j3,j4))
     &  +za(j1,j2)**2/(2d0*s(j2,j3))+za(j1,j2)**2/(2d0*s(j4,j1)))

c--- MODIFIED: overall sign of Rhat
      Rhat4mmpp=-Rhat4mmpp
       
      return
      end


c--- This is the hatted second version of the
c--- function presented in the paper
      double complex function C4mmpphat(j1,j2,j3,j4,za,zb)
      implicit none
C     arXiv:0704.3914v3, Eq.(3.24) (factor of C_\Gamma removed)
      include 'constants.f'
      include 'sprods_com.f'
      include 'scprods_com.f'
      include 'zprods_decl.f'
      include 'nflav.f'
      integer j1,j2,j3,j4,i
      double complex A0phiggggmmpp,Born,F31m,F42me,F41m,sum,
     . BGRL3hat,BGRL1,zab2,zba2
      double precision Np,bb0
      integer,parameter::ii(7)=(/1,2,3,4,1,2,3/)

      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
      zba2(j1,j2,j3,j4)=zb(j1,j2)*za(j2,j4)+zb(j1,j3)*za(j3,j4)

c--- set up 's-comma' products
      sc(1,1)=zip
      sc(1,2)=s(j1,j2)
      sc(1,3)=s(j1,j2)+s(j1,j3)+s(j2,j3)
      sc(1,4)=s(j1,j2)+s(j1,j3)+s(j1,j4)+s(j2,j3)+s(j2,j4)+s(j3,j4)
      sc(2,1)=sc(1,4) !s(j1,j2)+s(j1,j3)+s(j1,j4)+s(j2,j3)+s(j2,j4)+s(j3,j4)
      sc(2,2)=zip
      sc(2,3)=s(j2,j3)
      sc(2,4)=s(j2,j3)+s(j2,j4)+s(j3,j4)
      sc(3,1)=s(j3,j4)+s(j3,j1)+s(j4,j1)
      sc(3,2)=sc(1,4) !s(j1,j2)+s(j1,j3)+s(j1,j4)+s(j2,j3)+s(j2,j4)+s(j3,j4)
      sc(3,3)=zip
      sc(3,4)=s(j3,j4)
      sc(4,1)=s(j4,j1)
      sc(4,2)=s(j4,j1)+s(j4,j2)+s(j1,j2)
      sc(4,3)=sc(1,4) !s(j1,j2)+s(j1,j3)+s(j1,j4)+s(j2,j3)+s(j2,j4)+s(j3,j4)
      sc(4,4)=zip
      
      Np=2d0*(1d0-dfloat(nflav)/xn)
      bb0=(11d0*xn-2d0*dfloat(nflav))/3d0
      Born=A0phiggggmmpp(j1,j2,j3,j4,za,zb)

      sum=czip
      do i=1,4      
c--- note: corrected third argument of F42me from
c---       s_(i,j+1) to s_(i,j-1) [with j=i+3 here]
      sum=sum
     . +F31m(sc(ii(i),ii(i+2)))-F31m(sc(ii(i),ii(i+3)))
     . -0.5d0*F42me(sc(ii(i),ii(i+3)),sc(ii(i+1),ii(i+2))
     .             ,sc(ii(i),ii(i+2)),sc(ii(i+1),ii(i+3)))
     . -0.5d0*F41m(sc(ii(i),ii(i+2)),sc(ii(i),ii(i+1)),
     .                               sc(ii(i+1),ii(i+2)))
      enddo

c--- implementation of Eq. (3.24)
      C4mmpphat=Born*sum-(
     . +1d0/(za(j2,j3)*za(j3,j4)*za(j4,j1))*(
     . Np/6d0*za(j1,j4)*zb(j4,j3)*za(j3,j2)
     . *za(j1,j3)*zba2(j3,j4,j1,j2)
     . *(za(j1,j3)*zba2(j3,j4,j1,j2)-za(j1,j4)*zb(j4,j3)*za(j3,j2))
     . *BGRL3hat(sc(3,1),sc(4,1))
     . +bb0/xn*za(j1,j2)**2*za(j1,j4)*zb(j4,j3)*za(j3,j2)
     . *BGRL1(sc(3,1),sc(4,1))
     . +Np/6d0
     . *za(j1,j4)*zb(j4,j3)*za(j3,j2)*zab2(j1,j2,j3,j4)*za(j4,j2)
     . *(zab2(j1,j2,j3,j4)*za(j4,j2)-za(j1,j4)*zb(j4,j3)*za(j3,j2))
     . *BGRL3hat(sc(2,4),sc(2,3))
     . +bb0/xn*za(j1,j2)**2*za(j1,j4)*zb(j4,j3)*za(j3,j2)
     . *BGRL1(sc(2,4),sc(2,3))))

      return
      end


c--- This is the original (i.e. unhatted) second version of the
c--- function presented in the paper, which is suitable for hatting
      double complex function C4mmpp(j1,j2,j3,j4,za,zb)
      implicit none
C     arXiv:0704.3914v3, Eq.(3.24) (factor of C_\Gamma removed)
      include 'constants.f'
      include 'sprods_com.f'
      include 'scprods_com.f'
      include 'zprods_decl.f'
      include 'nflav.f'
      integer j1,j2,j3,j4,i
      double complex A0phiggggmmpp,Born,F31m,F42me,F41m,sum,
     . BGRL3,BGRL1,zab2,zba2
      double precision Np,bb0
      integer,parameter::ii(7)=(/1,2,3,4,1,2,3/)

      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
      zba2(j1,j2,j3,j4)=zb(j1,j2)*za(j2,j4)+zb(j1,j3)*za(j3,j4)

c--- set up 's-comma' products
      sc(1,1)=zip
      sc(1,2)=s(j1,j2)
      sc(1,3)=s(j1,j2)+s(j1,j3)+s(j2,j3)
      sc(1,4)=s(j1,j2)+s(j1,j3)+s(j1,j4)+s(j2,j3)+s(j2,j4)+s(j3,j4)
      sc(2,1)=sc(1,4) !s(j1,j2)+s(j1,j3)+s(j1,j4)+s(j2,j3)+s(j2,j4)+s(j3,j4)
      sc(2,2)=zip
      sc(2,3)=s(j2,j3)
      sc(2,4)=s(j2,j3)+s(j2,j4)+s(j3,j4)
      sc(3,1)=s(j3,j4)+s(j3,j1)+s(j4,j1)
      sc(3,2)=sc(1,4) !s(j1,j2)+s(j1,j3)+s(j1,j4)+s(j2,j3)+s(j2,j4)+s(j3,j4)
      sc(3,3)=zip
      sc(3,4)=s(j3,j4)
      sc(4,1)=s(j4,j1)
      sc(4,2)=s(j4,j1)+s(j4,j2)+s(j1,j2)
      sc(4,3)=sc(1,4) !s(j1,j2)+s(j1,j3)+s(j1,j4)+s(j2,j3)+s(j2,j4)+s(j3,j4)
      sc(4,4)=zip
      
      Np=2d0*(1d0-dfloat(nflav)/xn)
      bb0=(11d0*xn-2d0*dfloat(nflav))/3d0
      Born=A0phiggggmmpp(j1,j2,j3,j4,za,zb)

      sum=czip
      do i=1,4
c--- note: I have changed third argument of F42me from
c---       s_(i,j+1) to s_(i,j-1) [with j=i+3 here]
      sum=sum
     . +F31m(sc(ii(i),ii(i+2)))-F31m(sc(ii(i),ii(i+3)))
     . -0.5d0*F42me(sc(ii(i),ii(i+3)),sc(ii(i+1),ii(i+2))
     .             ,sc(ii(i),ii(i+2)),sc(ii(i+1),ii(i+3)))
     . -0.5d0*F41m(sc(ii(i),ii(i+2)),sc(ii(i),ii(i+1)),
     .                               sc(ii(i+1),ii(i+2)))
      enddo

c--- implementation of Eq. (3.24)
      C4mmpp=Born*sum-(
     . +1d0/(za(j2,j3)*za(j3,j4)*za(j4,j1))*(
     . Np/6d0*za(j1,j4)*zb(j4,j3)*za(j3,j2)
     . *za(j1,j3)*zba2(j3,j4,j1,j2)
     . *(za(j1,j3)*zba2(j3,j4,j1,j2)-za(j1,j4)*zb(j4,j3)*za(j3,j2))
     . *BGRL3(sc(3,1),sc(4,1))
     . +bb0/xn*za(j1,j2)**2*za(j1,j4)*zb(j4,j3)*za(j3,j2)
     . *BGRL1(sc(3,1),sc(4,1))
     . +Np/6d0
     . *za(j1,j4)*zb(j4,j3)*za(j3,j2)*zab2(j1,j2,j3,j4)*za(j4,j2)
     . *(zab2(j1,j2,j3,j4)*za(j4,j2)-za(j1,j4)*zb(j4,j3)*za(j3,j2))
     . *BGRL3(sc(2,4),sc(2,3))
     . +bb0/xn*za(j1,j2)**2*za(j1,j4)*zb(j4,j3)*za(j3,j2)
     . *BGRL1(sc(2,4),sc(2,3))))

      return
      end


      double complex function CR4mmpp(j1,j2,j3,j4,za,zb)
      implicit none
c     arXiv:0704.3914v3, Eq.(4.21) (factor of C_\Gamma removed)
C---  not yet implemented; instead, included in C by hatting procedure
      include 'constants.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'scprods_com.f'
      include 'nflav.f'
      integer j1,j2,j3,j4,j5
      double complex sum,zab2,zba2,zab3,zba3
      double precision Np

      zab2(j1,j2,j3,j4)=+za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
      zba2(j1,j2,j3,j4)=+zb(j1,j2)*za(j2,j4)+zb(j1,j3)*za(j3,j4)
      zab3(j1,j2,j3,j4,j5)=
     . +za(j1,j2)*zb(j2,j5)+za(j1,j3)*zb(j3,j5)+za(j1,j4)*zb(j4,j5)
      zba3(j1,j2,j3,j4,j5)=
     . +zb(j1,j2)*za(j2,j5)+zb(j1,j3)*za(j3,j5)+zb(j1,j4)*za(j4,j5)
      
c--- set up 's-comma' products
      sc(1,1)=zip
      sc(1,2)=s(j1,j2)
      sc(1,3)=s(j1,j2)+s(j1,j3)+s(j2,j3)
      sc(1,4)=s(j1,j2)+s(j1,j3)+s(j1,j4)+s(j2,j3)+s(j2,j4)+s(j3,j4)
      sc(2,1)=sc(1,4) !s(j1,j2)+s(j1,j3)+s(j1,j4)+s(j2,j3)+s(j2,j4)+s(j3,j4)
      sc(2,2)=zip
      sc(2,3)=s(j2,j3)
      sc(2,4)=s(j2,j3)+s(j2,j4)+s(j3,j4)
      sc(3,1)=s(j3,j4)+s(j3,j1)+s(j4,j1)
      sc(3,2)=sc(1,4) !s(j1,j2)+s(j1,j3)+s(j1,j4)+s(j2,j3)+s(j2,j4)+s(j3,j4)
      sc(3,3)=zip
      sc(3,4)=s(j3,j4)
      sc(4,1)=s(j4,j1)
      sc(4,2)=s(j4,j1)+s(j4,j2)+s(j1,j2)
      sc(4,3)=sc(1,4) !s(j1,j2)+s(j1,j3)+s(j1,j4)+s(j2,j3)+s(j2,j4)+s(j3,j4)
      sc(4,4)=zip
      
      Np=2d0*(1d0-dfloat(nflav)/xn)
       
      sum=
c--- 1st line
     . za(j1,j4)*zb(j4,j3)*za(j3,j2)*zab2(j1,j2,j3,j4)*za(j4,j2)
     . *(zab2(j1,j2,j3,j4)*za(j4,j2)-za(j1,j4)*zb(j4,j3)*za(j3,j2))
     . /(s(j2,j3)*zab2(j4,j2,j3,j4)**2)
c--- 2nd line
     .+za(j1,j4)*zba3(j4,j2,j3,j4,j2)*zab3(j1,j2,j3,j4,j4)*za(j4,j2)
     . *(zab3(j1,j2,j3,j4,j4)*za(j4,j2)-za(j1,j4)*zba3(j4,j2,j3,j4,j2))
     . /(sc(2,4)*zab3(j4,j2,j3,j4,j4)**2) 
c--- 3rd line
     .+zab3(j1,j3,j4,j1,j3)*za(j3,j2)*za(j1,j3)*zba3(j3,j3,j4,j1,j2)
     . *(za(j1,j3)*zba3(j3,j3,j4,j1,j2)-zab3(j1,j3,j4,j1,j3)*za(j3,j2))
     . /(sc(3,1)*zab3(j3,j3,j4,j1,j3)**2)
c--- 4th line
     .+za(j1,j4)*zb(j4,j3)*za(j3,j2)*za(j1,j3)*zba2(j3,j4,j1,j2)
     . *(za(j1,j3)*zba2(j3,j4,j1,j2)-za(j1,j4)*zb(j4,j3)*za(j3,j2))
     . /(s(j4,j1)*zab2(j3,j4,j1,j3)**2)

      CR4mmpp=sum*Np/(12d0*za(j4,j1)*za(j2,j3)*za(j3,j4))

      return
      end







c--- NOTE: these functions are not used, but they may be useful
c---       for performing double-checks

c--- This is the first version of the function presented in the
c--- paper; it should agree with the implementation of Eq. (3.24)
      double complex function C4mmpp322(j1,j2,j3,j4,za,zb)
      implicit none
C     arXiv:0704.3914v3, Eq.(3.22) (factor of C_\Gamma removed)
      include 'constants.f'
      include 'sprods_com.f'
      include 'scprods_com.f'
      include 'zprods_decl.f'
      include 'nflav.f'
      integer j1,j2,j3,j4,i
      double complex A0phiggggmmpp,Born,F31m,F42me,F41m,sum,
     . BGRL3,BGRL2,BGRL1,trm,trm1432,trm2341
      double precision Np,bb0
      integer,parameter::ii(7)=(/1,2,3,4,1,2,3/)
      trm(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j3)*za(j3,j4)*zb(j4,j1)


c--- set up 's-comma' products
      sc(1,1)=zip
      sc(1,2)=s(j1,j2)
      sc(1,3)=s(j1,j2)+s(j1,j3)+s(j2,j3)
      sc(1,4)=s(j1,j2)+s(j1,j3)+s(j1,j4)+s(j2,j3)+s(j2,j4)+s(j3,j4)
      sc(2,1)=sc(1,4) !s(j1,j2)+s(j1,j3)+s(j1,j4)+s(j2,j3)+s(j2,j4)+s(j3,j4)
      sc(2,2)=zip
      sc(2,3)=s(j2,j3)
      sc(2,4)=s(j2,j3)+s(j2,j4)+s(j3,j4)
      sc(3,1)=s(j3,j4)+s(j3,j1)+s(j4,j1)
      sc(3,2)=sc(1,4) !s(j1,j2)+s(j1,j3)+s(j1,j4)+s(j2,j3)+s(j2,j4)+s(j3,j4)
      sc(3,3)=zip
      sc(3,4)=s(j3,j4)
      sc(4,1)=s(j4,j1)
      sc(4,2)=s(j4,j1)+s(j4,j2)+s(j1,j2)
      sc(4,3)=sc(1,4) !s(j1,j2)+s(j1,j3)+s(j1,j4)+s(j2,j3)+s(j2,j4)+s(j3,j4)
      sc(4,4)=zip
      
      trm1432=trm(j1,j4,j3,j2)
      trm2341=trm(j2,j3,j4,j1)
      Np=2d0*(1d0-dfloat(nflav)/xn)
      bb0=(11d0*xn-2d0*dfloat(nflav))/3d0
      Born=A0phiggggmmpp(j1,j2,j3,j4,za,zb)
      sum=czip
      do i=1,4
c--- note: I have changed third argument of F42me from
c---       s_(i,j+1) to s_(i,j-1) [with j=i+3 here]
      sum=sum
     . +F31m(sc(ii(i),ii(i+2)))-F31m(sc(ii(i),ii(i+3)))
     . -0.5d0*F42me(sc(ii(i),ii(i+3)),sc(ii(i+1),ii(i+2))
     .             ,sc(ii(i),ii(i+2)),sc(ii(i+1),ii(i+3)))
     . -0.5d0*F41m(sc(ii(i),ii(i+2)),sc(ii(i),ii(i+1)),
     .                               sc(ii(i+1),ii(i+2)))
      enddo

      sum=sum-(
     . +Np/3d0*(trm1432/s(j1,j2))**3*BGRL3(sc(3,1),sc(4,1))    
     . +Np/3d0*(trm2341/s(j1,j2))**3*BGRL3(sc(2,4),sc(2,3))    
     . -Np/2d0*(trm1432/s(j1,j2))**2*BGRL2(sc(3,1),sc(4,1))    
     . -Np/2d0*(trm2341/s(j1,j2))**2*BGRL2(sc(2,4),sc(2,3))    
     . +(Np/6d0+bb0/xn)*(trm1432/s(j1,j2))*BGRL1(sc(3,1),sc(4,1)) 
     . +(Np/6d0+bb0/xn)*(trm2341/s(j1,j2))*BGRL1(sc(2,4),sc(2,3))) 

      C4mmpp322=Born*sum

      return
      end



