c--- Results taken from an early draft of:
c---   S.~Badger, E.~W.~N.~Glover, P.~Mastrolia and C.~Williams,
c---   %``One-loop Higgs plus four gluon amplitudes: Full analytic results,''
c---   arXiv:0909.4475 [hep-ph].
c---   %%CITATION = ARXIV:0909.4475;%%

      double complex function A1Hggggpmmm(j1,j2,j3,j4,za,zb)
      implicit none
      integer j1,j2,j3,j4
      include 'constants.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      double complex A0phiggggpmmm,V4,FR4pmmm,F31m
      
      V4 = -A0phiggggpmmm(j1,j2,j3,j4,za,zb)*(
     &  +F31m(s(j1,j2))+F31m(s(j2,j3))+F31m(s(j3,j4))+F31m(s(j4,j1)))
     
      A1Hggggpmmm=V4-FR4pmmm(j1,j2,j3,j4,za,zb)

      return
      end

      double complex function W(mhsq,s234,s12,s23,s34,s14)
      implicit none
      double complex F41mF_BGMW,F42mhF_BGMW
      double precision mhsq,s234,s12,s23,s34,s14
c--- note reference to BGMW version of the finite piece of the
c--- box integrals, which is twice the definition in previous papers
c--- Note: second term has 3rd and 4th arguments switched
c---       compared with BGMW paper
      W=F41mF_BGMW(s234,s23,s34)
     . +F42mhF_BGMW(mhsq,s23,s14,s234)
     . +F42mhF_BGMW(mhsq,s34,s12,s234)

      return
      end

c--- these are just aliases to previously-defined finite pieces
c--- of the box functions, but with an additional factor of two
      double complex function F41mF_BGMW(psq,s,t)
      implicit none
      double complex F41mF
      double precision psq,s,t
      
      F41mF_BGMW=2d0*F41mF(psq,s,t)
      
      return
      end
      

      double complex function F42mhF_BGMW(psq,qsq,s,t)
      implicit none
      double complex F42mhF
      double precision psq,qsq,s,t
      
      F42mhF_BGMW=2d0*F42mhF(psq,qsq,s,t)

      return
      end
      

      double complex function FR4unsym(j1,j2,j3,j4,za,zb)
c--- This function is the sum of the finite cut-constructible
c--- and rational terms from the unpublished BGMW paper,
c--- i.e. it is the sum of F4 in Eq. (8.18) and R4 in Eq. (8.25)
      implicit none
      include 'constants.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'nflav.f'
      integer j1,j2,j3,j4,j
      double complex zab2,W,W1,W2,W3,
     . BGRL1,BGRL2hat,BGRL3hat,F41mF_BGMW,F33m,lnrat
      double precision mhsq,s123,s234,s134,s124,s12,s13,s14,s23,s24,s34
      double complex k1sq,k2sq,k1Dk2,gamma,factor,coef3mass,a1,a2,a3,a4
c----statement function
      zab2(j1,j2,j3,j4)=+za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
c----statement function

      mhsq=s(j1,j2)+s(j1,j3)+s(j1,j4)+s(j2,j3)+s(j2,j4)+s(j3,j4)
      s234=s(j2,j3)+s(j3,j4)+s(j2,j4)
      s134=s(j1,j3)+s(j3,j4)+s(j1,j4)
      s124=s(j1,j2)+s(j2,j4)+s(j1,j4)
      s123=s(j1,j2)+s(j2,j3)+s(j1,j3)
      s12=s(j1,j2)
      s13=s(j1,j3)
      s14=s(j1,j4)
      s23=s(j2,j3)
      s24=s(j2,j4)
      s34=s(j3,j4)
      W1=W(mhsq,s234,s12,s23,s34,s14)
      W2=W(mhsq,s134,s23,s34,s14,s12)
      W3=W(mhsq,s124,s34,s14,s12,s23)

C     p12flat=gamma/(gamma**2-k1sq*k2sq)*(gamma*(p1+p2)-k1sq*(p3+p4))
C----solve for gamma_+ and gamma_-

      k1sq=dcmplx(mhsq)
      k2sq=dcmplx(s12)
      k1Dk2=-dcmplx(s12+0.5d0*(s13+s14+s23+s24))
C-gamma+ = k1Dk2+sqrt(k1Dk2**2-k1sq*k2sq)
C-gamma- = k1Dk2-sqrt(k1Dk2**2-k1sq*k2sq)
      
      coef3mass=czip
      gamma=k1Dk2+sqrt(k1Dk2**2-k1sq*k2sq)      
      
      do j=1,2
      factor=gamma/(gamma**2-k1sq*k2sq)
      a1=factor*(-gamma-k1sq)
      a2=a1
      a3=factor*(-gamma)
      a4=a3
      
      coef3mass=coef3mass+mhsq**2*za(j3,j4)**3
     &  *(+a3*za(j2,j3)*zb(j3,j1)+a4*za(j2,j4)*zb(j4,j1))
     &  *(+a1*za(j2,j1)*zb(j1,j3)+a4*za(j2,j4)*zb(j4,j3))
     &  *(+a1*za(j2,j1)*zb(j1,j4)+a3*za(j2,j3)*zb(j3,j4))
     & /(gamma*(gamma+mhsq)*za(j1,j2)
     &  *(a2*s(j1,j2)+a3*s(j1,j3)+a4*s(j1,j4))
     &  *(a1*s(j1,j3)+a2*s(j2,j3)+a4*s(j3,j4))
     &  *(a1*s(j1,j4)+a2*s(j2,j4)+a3*s(j3,j4)))
      gamma=2d0*k1Dk2-gamma
      enddo

c---  Eq.(8.18)
      FR4unsym=-s234**3
     . /(4d0*zab2(j1,j3,j4,j2)*zab2(j1,j2,j3,j4)*zb(j2,j3)*zb(j3,j4))*W1
     .   *(-1d0)

     . +(zab2(j2,j3,j4,j1)**3
     . /(2d0*s134*zab2(j2,j1,j4,j3)*zb(j3,j4)*zb(j4,j1))
     .  +za(j3,j4)**3*mhsq**2
     .  /(2d0*s134*zab2(j1,j3,j4,j2)*zab2(j3,j1,j4,j2)*za(j4,j1)))*W2

     . +0.25d0/s124
     . *(zab2(j3,j2,j4,j1)**4
     . /(zab2(j3,j1,j4,j2)*zab2(j3,j1,j2,j4)*zb(j2,j1)*zb(j4,j1))
     . + za(j2,j4)**4*mhsq**2
     . /(za(j1,j2)*za(j1,j4)*zab2(j2,j1,j4,j3)*zab2(j4,j1,j2,j3)))*W3
     .   *(-1d0)*(1d0)
      
     . +coef3mass*F33m(mhsq,s12,s34)

     . +(1d0-dfloat(nflav)/4d0/xn)
     . *(2d0*zab2(j3,j2,j4,j1)**2/(s124*zb(j2,j4)**2)
     . *F41mF_BGMW(s124,s12,s14)*(-0.5d0)
     . +4d0*za(j2,j4)*zab2(j3,j2,j4,j1)**2/(s124*zb(j4,j2))
     . *BGRL1(s124,s12)
     . -4d0*za(j2,j3)*zab2(j4,j3,j2,j1)**2/(s123*zb(j3,j2))
     . *BGRL1(s123,s12))

     . +(1d0-dfloat(nflav)/xn)*(
     . zb(j1,j2)*zb(j4,j1)*zab2(j3,j1,j4,j2)*zab2(j3,j1,j2,j4)
     . /(2d0*s124*zb(j2,j4)**4)*F41mF_BGMW(s124,s12,s14)
     . +2d0*s124*za(j2,j4)*(za(j3,j4)*zb(j4,j1))**2/(3d0*zb(j4,j2))
     . *BGRL3hat(s124,s12)
     . +za(j3,j4)*zb(j4,j1)
     . *(3d0*s124*za(j3,j4)*zb(j4,j1)
     . -za(j2,j4)*zab2(j3,j2,j4,j1)*zb(j4,j2))
     . /(3d0*zb(j4,j2)**2)*BGRL2hat(s124,s12)

     . +(2d0*s124*(za(j3,j4)*zb(j4,j1))**2/(za(j2,j4)*zb(j4,j2)**3)
     . -za(j2,j4)*zab2(j3,j2,j4,j1)**2/(3d0*s124*zb(j4,j2)))
     . *BGRL1(s124,s12)

     . -zab2(j3,j2,j4,j1)*(4d0*s124*za(j3,j4)*zb(j4,j1)
     .  -zab2(j3,j2,j4,j1)*(2d0*s14+s24))/(s124*za(j2,j4)*zb(j4,j2)**3)
     . *lnrat(-s124,-s12)

     . -2d0*s123*za(j2,j3)*(za(j3,j4)*zb(j3,j1))**2/(3d0*zb(j3,j2))
     . *BGRL3hat(s123,s12)
     . -za(j2,j3)*za(j3,j4)*zb(j3,j1)*zab2(j4,j2,j3,j1)/(3d0*zb(j3,j2))
     . *BGRL2hat(s123,s12)

     . +za(j2,j3)*zab2(j4,j2,j3,j1)**2/(3d0*s123*zb(j3,j2))
     . *BGRL1(s123,s12))

C-----Add rational piece Eq.(8.25)
      FR4unsym=FR4unsym
     . -0.5d0*(1d0-dfloat(nflav)/xn)
     . *(-za(j2,j3)*za(j3,j4)*zab2(j4,j2,j3,j1)*zb(j3,j1) 
     . /(3d0*s123*za(j1,j2)*zb(j2,j1)*zb(j3,j2))
     . -zab2(j3,j2,j4,j1)**2/(s124*zb(j4,j2)**2)
     . -za(j2,j4)*za(j3,j4)*zab2(j3,j2,j4,j1)*zb(j4,j1)
     . /(3d0*s124*s12*zb(j4,j2))
     . -(zb(j1,j2)*za(j2,j3))**2/(s14*zb(j4,j2)**2)
     . -za(j2,j4)*(s23*s24+s23*s34+s24*s34)
     . /(3d0*za(j1,j2)*za(j1,j4)*zb(j2,j3)*zb(j3,j4)*zb(j4,j2))
     . +zab2(j2,j3,j4,j1)*zab2(j4,j2,j3,j1)
     . /(3d0*s234*zb(j2,j3)*zb(j3,j4))
     . -2d0*zb(j1,j2)*za(j2,j3)*zb(j3,j1)**2
     . /(3d0*zb(j2,j3)**2*zb(j4,j1)*zb(j3,j4)))

      return
      end
      


      double complex function FR4pmmm(j1,j2,j3,j4,za,zb)
      implicit none
      include 'constants.f'
      include 'zprods_decl.f'
      integer j1,j2,j3,j4
      double complex FR4unsym

      FR4pmmm=FR4unsym(j1,j2,j3,j4,za,zb)+FR4unsym(j1,j4,j3,j2,za,zb)

      return
      end
      

