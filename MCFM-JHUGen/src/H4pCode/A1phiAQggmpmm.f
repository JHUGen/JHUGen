      double complex function A1phiAQggmpmmL(k1,k2,k3,k4,za,zb)
c--- This is an implementation of Eq. (5.2) in
c---  S.~Badger, John.~M.~Campbell, R.~Keith Ellis and Ciaran Williams
c---  "Analytic results for the one-loop NMHV H-qbar-q-g-g amplitude."
c---   preprint DESY 09-180, FERMILAB-PUB-09-505-T, IPPP/09/86
c---   arXiv: 0910.4481 [hep-ph]
      implicit none
      include 'constants.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'scale.f'
      include 'epinv.f'
      include 'deltar.f'
      integer j,k1,k2,k3,k4
      double complex V1L,A0phiAQggmpmm,lnrat,zab2,Lsm1,Lsm1_2mht
      double complex sum,l23,l34,l41,l12,coef3m1234,coef3m1423,
     & S1,S2,K1DK2,a1,a2,a3,a4,gamma,factor,I3m,
     & BGRL1,BGRL2hat,BGRL3hat
      double precision s12,s13,s14,s23,s24,s34,s123,s234,s124,s134,mhsq
      zab2(k1,k2,k3,k4)=za(k1,k2)*zb(k2,k4)+za(k1,k3)*zb(k3,k4)
      s12=s(k1,k2)
      s13=s(k1,k3)
      s14=s(k1,k4)
      s23=s(k2,k3)
      s24=s(k2,k4)
      s34=s(k3,k4)
      mhsq=s12+s13+s14+s23+s24+s34
      s123=s12+s13+s23
      s124=s12+s14+s24
      s134=s13+s14+s34
      s234=s23+s24+s34

      l12=lnrat(musq,-s12)
      l23=lnrat(musq,-s23)
      l34=lnrat(musq,-s34)
      l41=lnrat(musq,-s14)


      V1L=
     . -epinv**2-epinv*l23-0.5d0*l23**2
     . -epinv**2-epinv*l34-0.5d0*l34**2
     . -epinv**2-epinv*l41-0.5d0*l41**2
     . +13d0/6d0*(epinv+l12)+119d0/18d0-deltar/6d0

      A1phiAQggmpmmL=A0phiAQggmpmm(k1,k2,k3,k4,za,zb)*V1L

      sum=
     & -s134**2/(zb(k4,k1)*zb(k3,k4)*zab2(k2,k1,k4,k3))
     & *Lsm1(-s14,-s134,-s34,-s134)

     & -zab2(k1,k3,k4,k2)**2/(zab2(k1,k2,k3,k4)*zb(k2,k3)*zb(k3,k4))
     & *Lsm1(-s34,-s234,-s23,-s234)

     &  +(mhsq**2*za(k1,k4)**2*za(k2,k4)
     & /(za(k1,k2)*zab2(k2,k1,k4,k3)*zab2(k4,k1,k2,k3)*s124)
     &   -zab2(k3,k1,k4,k2)**3
     & /(zb(k1,k2)*zb(k2,k4)*zab2(k3,k1,k2,k4)*s124))
     & *Lsm1(-s14,-s124,-s12,-s124)

     & +(zb(k2,k3)**2*zab2(k4,k2,k3,k1)**3
     & /(zb(k1,k2)*zb(k1,k3)**3*zab2(k4,k1,k2,k3)*s123)
     & -mhsq**2*za(k1,k3)**3
     & /(za(k1,k2)*zab2(k1,k2,k3,k4)*zab2(k3,k1,k2,k4)*s123))
     & *Lsm1(-s12,-s123,-s23,-s123)

     &   +za(k3,k4)*s134**2
     & /(zb(k3,k4)*za(k3,k4)*zb(k1,k4)*zab2(k2,k1,k4,k3))
     & *Lsm1_2mht(s12,s134,s34,mhsq)

     & -zab2(k1,k3,k4,k2)**2
     & /(zb(k3,k4)*zb(k2,k3)*zab2(k1,k2,k3,k4))
     & *Lsm1_2mht(s12,s234,s34,mhsq)

     & +(zab2(k4,k1,k3,k2)**3
     & /(zb(k1,k2)*zb(k2,k3)*zab2(k4,k1,k2,k3)*s123)
     & -mhsq**2*za(k1,k3)**3
     & /(za(k1,k2)*zab2(k1,k2,k3,k4)*zab2(k3,k1,k2,k4)*s123))
     & *Lsm1_2mht(s34,s123,s12,mhsq)

     & +(mhsq**2*za(k1,k4)**2*za(k2,k4)
     & /(za(k1,k2)*zab2(k2,k1,k4,k3)*zab2(k4,k1,k2,k3)*s124)
     & -zab2(k3,k2,k4,k1)*zab2(k3,k1,k4,k2)**2
     & /(zab2(k3,k1,k2,k4)*zb(k1,k4)*zb(k1,k2)*s124))
     & *Lsm1_2mht(s34,s124,s12,mhsq)

     &  +(zab2(k4,k1,k3,k2)**3
     & /(zb(k1,k2)*zb(k2,k3)*zab2(k4,k1,k2,k3)*s123)
     & -mhsq**2*za(k1,k3)**3
     & /(za(k1,k2)*zab2(k1,k2,k3,k4)*zab2(k3,k1,k2,k4)*s123))
     & *Lsm1_2mht(s14,s123,s23,mhsq)

     & +(mhsq**2*za(k2,k4)*za(k1,k4)**2
     & /(za(k1,k2)*zab2(k2,k1,k4,k3)*zab2(k4,k1,k2,k3)*s124)
     & -zab2(k3,k1,k4,k2)**2*zab2(k3,k2,k4,k1)
     & /(zb(k1,k2)*zb(k1,k4)*zab2(k3,k1,k2,k4)*s124))
     & *Lsm1_2mht(s23,s124,s14,mhsq)

C-----Now for three mass triangles
C-----Deal with the 12-34 case first
C-----K1=-(p1+p2+p3+p4)
C-----K2=-(p1+p2)
C     K1flat=gamma/(gamma**2-S1*S2)*(gamma*K1-S1*K2)
C     K1flat=gamma/(gamma**2-S1*S2)*((S1-gamma)*(p1+p2)-gamma*((p3+p4))
C----solve for gamma_+ and gamma_-
      S1=dcmplx(mhsq)
      S2=dcmplx(s12)
      K1DK2=dcmplx(s12+0.5d0*(s13+s14+s23+s24))

C-gamma+ = K1DK2+sqrt(K1DK2**2-S1*S2)
C-gamma- = K1DK2-sqrt(K1DK2**2-S1*S2)

      coef3m1234=czip
      gamma=K1DK2+sqrt(K1DK2**2-S1*S2)

      do j=1,2
C -- calculate the projections of K1 flat on k1,k2,k3,k4 called a1,a2,a3,a4
      factor=gamma/(gamma**2-S1*S2)
      a1=factor*(S1-gamma)
      a2=a1
      a3=-factor*gamma
      a4=a3

c---  The result for the 12-34 coefficient is
c       + iza(k1,k2)*iza(k3,k1f)*iza(k4,k1f) * (
c          + za(k1,k1f)^2*za(k3,k4)^3*ga^-1*S1^2*[ga-S1]^-1);
c     which we can rewrite as
c       + iza(k1,k2) * (
c          + 1/4*za(k3,k4)^3*zab2(k1,k1f,k3)*zab2(k1,k1f,k4)*k3.k1f^-1*
c         k4.k1f^-1*ga^-1*S1^2*[ga-S1]^-1);

      coef3m1234=coef3m1234
     & +mhsq**2*za(k3,k4)**3/(za(k1,k2)*gamma*(gamma-mhsq))
     & *(a2*za(k1,k2)*zb(k2,k3)+a4*za(k1,k4)*zb(k4,k3))   ! zab2(k1,k1f,k3)
     & *(a2*za(k1,k2)*zb(k2,k4)+a3*za(k1,k3)*zb(k3,k4))   ! zab2(k1,k1f,k4)
     * /(a1*s13+a2*s23+a4*s34)                            ! (2*k3.k1f)^-1*
     * /(a1*s14+a2*s24+a3*s34)                            ! (2*k4.k1f)^-1*

C----switch to other solution
      gamma=2*K1DK2-gamma

      enddo


C-----Now deal with the 14-23 case
C-----K1=-(p1+p2+p3+p4)
C-----K2=-(p1+p4)
C     K1flat=gamma/(gamma**2-S1*S2)*(gamma*K1-S1*K2)
C     K1flat=gamma/(gamma**2-S1*S2)*((S1-gamma)*(p1+p4)-gamma*((p2+p3))

C----solve for gamma_+ and gamma_-
      S1=dcmplx(mhsq)
      S2=dcmplx(s14)
      K1DK2=dcmplx(s14+0.5d0*(s12+s13+s24+s34))

C-gamma+ = K1DK2+sqrt(K1DK2**2-S1*S2)
C-gamma- = K1DK2-sqrt(K1DK2**2-S1*S2)

      coef3m1423=czip
      gamma=K1DK2+sqrt(K1DK2**2-S1*S2)

      do j=1,2
C -- calculate the projections of K1 flat on k1,k2,k3,k4 called a1,a2,a3,a4
      factor=gamma/(gamma**2-S1*S2)
      a1=factor*(S1-gamma)
      a4=a1
      a3=-factor*gamma
      a2=a3

C=--  The result for the 14-23 coefficient is
c       -za(k1,k4)^2*za(k3,k1f)^2*ga^-1*S1^2*[ga-S1]^-1
c           + iza(k1,k1f)*iza(k2,k1f) *
c     which we can rewrite as
c       -za(k1,k4)^2*zab2(k3,k1f,k1)*zab2(k3,k1f,k2)*S1^2
c         /(2*k1.k1f^-1)/*2*k2.k1f^-1*ga^-1*[ga-S1]^-1

      coef3m1423=coef3m1423
     &  -mhsq**2*za(k1,k4)**2/(2d0*gamma*(gamma-mhsq))
     & *(a2*za(k3,k2)*zb(k2,k1)+a4*za(k3,k4)*zb(k4,k1))   ! *zab2(k3,k1f,k1)
     & *(a1*za(k3,k1)*zb(k1,k2)+a4*za(k3,k4)*zb(k4,k2))   ! *zab2(k3,k1f,k2)
     * /(a2*s12+a3*s13+a4*s14)                            ! (2*k1.k1f)^-1*
     * /(a1*s12+a3*s23+a4*s24)                            ! (2*k2.k1f)^-1*


C----switch to other solution
      gamma=2*K1DK2-gamma
      enddo

      sum=sum-coef3m1234*I3m(mhsq,s12,s34)
     .       -coef3m1423*I3m(mhsq,s14,s23)

      sum=sum
     & -2d0/3d0*za(k1,k3)**2*za(k3,k4)*zab2(k4,k1,k2,k3)*zb(k1,k2)
     & *BGRL3hat(s123,s12)

     & +1d0/6d0*za(k3,k4)*za(k1,k3)
     & *(zab2(k4,k1,k3,k2)*zb(k1,k3)-3d0*zab2(k4,k2,k3,k1)*zb(k2,k3))
     & /zb(k1,k3)
     * *BGRL2hat(s123,s12)

     & +za(k1,k3)
     & *(0.5d0*zab2(k4,k1,k3,k2)*zab2(k4,k1,k2,k3)*zb(k1,k2)*zb(k1,k3)
     & -zab2(k4,k2,k3,k1)**2*zb(k2,k3)**2
     & -8d0/3d0*zab2(k4,k1,k3,k2)**2*zb(k1,k3)**2)
     & /(s123*zb(k1,k3)**2*zb(k2,k3))
     & *BGRL1(s123,s12)

     & -2d0/3d0*s124*za(k3,k4)**2*za(k1,k4)*zb(k4,k2)
     & *BGRL3hat(s124,s12)

     & +za(k3,k4)*za(k1,k4)
     & *(1d0/3d0*zab2(k3,k1,k4,k2)*zb(k1,k4)
     * -0.5d0*zab2(k3,k1,k2,k4)*zb(k1,k2))/zb(k1,k4)
     & *BGRL2hat(s124,s12)

     & +zab2(k3,k1,k4,k2)*(3d0/2d0*s124*za(k3,k4)
     & +11d0/3d0*zab2(k3,k1,k4,k2)*za(k4,k2))/(s124*zb(k1,k4))
     & *BGRL1(s124,s12)

     & +0.5d0*za(k1,k4)*za(k1,k3)*zab2(k4,k2,k3,k1)*zb(k1,k2)/zb(k3,k1)
     & *BGRL2hat(s123,s23)

     &-za(k1,k3)*zab2(k4,k2,k3,k1)*(3d0/2d0*zab2(k4,k1,k3,k2)*zb(k1,k3)
     &+zab2(k4,k2,k3,k1)*zb(k2,k3))/(s123*zb(k1,k3)**2)
     & *BGRL1(s123,s23)

     & +0.5d0*s234*za(k1,k4)*za(k3,k4)*zb(k4,k2)/zb(k4,k3)
     & *BGRL2hat(s234,s23)

     & +3d0/2d0*za(k3,k4)*zab2(k1,k3,k4,k2)/zb(k4,k3)
     & *BGRL1(s234,s23)

      A1phiAQggmpmmL=A1phiAQggmpmmL+sum

c--- now add the rational pieces
      sum=
     . za(k3,k4)*zab2(k3,k1,k4,k2)
     .  *(2d0*za(k2,k4)*zb(k4,k2)-za(k1,k2)*zb(k2,k1))
     .  /(12d0*s124*za(k1,k2)*zb(k2,k1)*zb(k4,k1))
     . +(za(k2,k3)*zab2(k4,k1,k3,k2)**2*(
     .    3d0*za(k1,k2)*zb(k2,k1)-2d0*za(k2,k3)*zb(k3,k2))
     .  -2d0*za(k1,k3)**2*za(k2,k4)*zab2(k4,k2,k3,k1)
     .      *zb(k2,k1)*zb(k3,k2))
     .  /(12d0*s123*za(k1,k2)*za(k2,k3)*zb(k2,k1)*zb(k3,k1)*zb(k3,k2))
     . +5d0*za(k3,k4)**2/(12d0*za(k2,k3)*zb(k3,k1))
     . +5d0*za(k3,k4)*zab2(k4,k1,k3,k2)
     .  /(6d0*za(k2,k3)*zb(k3,k1)*zb(k3,k2))
     . +zab2(k4,k1,k3,k2)**2
     .  /(6d0*za(k1,k2)*zb(k2,k1)*zb(k3,k1)*zb(k3,k2))
     . -za(k1,k3)*za(k1,k4)*za(k2,k4)*zb(k2,k1)
     .  /(3d0*za(k1,k2)*za(k2,k3)*zb(k3,k1)*zb(k3,k2))
     . -za(k1,k3)*za(k3,k4)/(12d0*za(k1,k2)*zb(k4,k1))
     . -za(k3,k4)**2*zb(k4,k2)/(6d0*za(k1,k2)*zb(k2,k1)*zb(k4,k1))
     . +za(k1,k3)*za(k2,k4)*zab2(k4,k1,k3,k4)
     .  /(4d0*za(k1,k2)*za(k2,k3)*zb(k3,k1)*zb(k4,k3))
     . -za(k1,k3)*zab2(k4,k1,k3,k4)/(3d0*za(k1,k2)*zb(k4,k1)*zb(k4,k3))
     . -5d0*za(k1,k4)**2*zb(k4,k1)/(12d0*za(k1,k2)*zb(k3,k1)*zb(k4,k3))
     . +za(k1,k4)**2*zb(k4,k2)/(6d0*za(k1,k2)*zb(k3,k2)*zb(k4,k3))

      A1phiAQggmpmmL=A1phiAQggmpmmL+sum

      return
      end


      double complex function A1phiAQggmpmmR(k1,k2,k3,k4,za,zb)
c--- This is an implementation of Eq. (5.10) in
c---  S.~Badger, John.~M.~Campbell, R.~Keith Ellis and Ciaran Williams
c---  "Analytic results for the one-loop NMHV H-qbar-q-g-g amplitude."
c---   preprint DESY 09-180, FERMILAB-PUB-09-505-T, IPPP/09/86
c---   arXiv: 0910.4481 [hep-ph]
      implicit none
      include 'constants.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'scale.f'
      include 'epinv.f'
      include 'deltar.f'
      integer j,k1,k2,k3,k4
      double complex VR,A0phiAQggmpmm,lnrat,zab2,Lsm1,Lsm1_2mht
      double complex sum,l12,coef3m1423,
     & S1,S2,K1DK2,a1,a2,a3,a4,gamma,factor,I3m,
     & BGRL1,BGRL2hat
      double precision s12,s13,s14,s23,s24,s34,s123,s234,s124,s134,mhsq
      zab2(k1,k2,k3,k4)=za(k1,k2)*zb(k2,k4)+za(k1,k3)*zb(k3,k4)
      s12=s(k1,k2)
      s13=s(k1,k3)
      s14=s(k1,k4)
      s23=s(k2,k3)
      s24=s(k2,k4)
      s34=s(k3,k4)
      mhsq=s12+s13+s14+s23+s24+s34
      s123=s12+s13+s23
      s124=s12+s14+s24
      s134=s13+s14+s34
      s234=s23+s24+s34

      l12=lnrat(musq,-s12)

      VR=-epinv**2-epinv*l12-0.5d0*l12**2
     .   -3d0/2d0*(epinv+l12)-7d0/2d0-deltar/2d0
      A1phiAQggmpmmR=A0phiAQggmpmm(k1,k2,k3,k4,za,zb)*VR

      sum=
     &  +zb(k1,k2)**2*zab2(k4,k1,k2,k3)**2/(zb(k1,k3)**3*zb(k2,k3)*s123)
     & *Lsm1(-s12,-s123,-s23,-s123)

     & +zab2(k3,k1,k4,k2)**2/(zb(k1,k4)*zb(k2,k4)*s124)
     & *Lsm1(-s14,-s124,-s12,-s124)

     & -zab2(k1,k3,k4,k2)**2/(zb(k2,k3)*zb(k3,k4)*zab2(k1,k2,k3,k4))
     & *Lsm1_2mht(s14,s234,s23,mhsq)

     & +s134**2/(zb(k1,k4)*zb(k3,k4)*zab2(k2,k1,k4,k3))
     & *Lsm1_2mht(s23,s134,s14,mhsq)

C-----Now deal with the 14-23 case
C-----K1=-(p1+p2+p3+p4)
C-----K2=-(p1+p4)
C     K1flat=gamma/(gamma**2-S1*S2)*(gamma*K1-S1*K2)
C     K1flat=gamma/(gamma**2-S1*S2)*((S1-gamma)*(p1+p4)-gamma*((p2+p3))

C----solve for gamma_+ and gamma_-
      S1=dcmplx(mhsq)
      S2=dcmplx(s14)
      K1DK2=dcmplx(s14+0.5d0*(s12+s13+s24+s34))

C-gamma+ = K1DK2+sqrt(K1DK2**2-S1*S2)
C-gamma- = K1DK2-sqrt(K1DK2**2-S1*S2)

      coef3m1423=czip
      gamma=K1DK2+sqrt(K1DK2**2-S1*S2)

      do j=1,2
C -- calculate the projections of K1 flat on k1,k2,k3,k4 called a1,a2,a3,a4
      factor=gamma/(gamma**2-S1*S2)
      a1=factor*(S1-gamma)
      a4=a1
      a3=-factor*gamma
      a2=a3

C=--  The result for the 14-23 coefficient is
c       -za(k1,k4)^2*za(k3,k1f)^2*ga^-1*S1^2*[ga-S1]^-1
c           + iza(k1,k1f)*iza(k2,k1f) *
c     which we can rewrite as
c       -za(k1,k4)^2*zab2(k3,k1f,k1)*zab2(k3,k1f,k2)*S1^2
c         /(2*k1.k1f^-1)/*2*k2.k1f^-1*ga^-1*[ga-S1]^-1
C---- NB   Factor of 1/2 added over and above the form, to get numerical agreement

      coef3m1423=coef3m1423
     &  -mhsq**2*za(k1,k4)**2/(2d0*gamma*(gamma-mhsq))
     & *(a2*za(k3,k2)*zb(k2,k1)+a4*za(k3,k4)*zb(k4,k1))   ! *zab2(k3,k1f,k1)
     & *(a1*za(k3,k1)*zb(k1,k2)+a4*za(k3,k4)*zb(k4,k2))   ! *zab2(k3,k1f,k2)
     * /(a2*s12+a3*s13+a4*s14)                            ! (2*k1.k1f)^-1*
     * /(a1*s12+a3*s23+a4*s24)                            ! (2*k2.k1f)^-1*


C----switch to other solution
      gamma=2*K1DK2-gamma
      enddo

      sum=sum-coef3m1423*I3m(mhsq,s14,s23)

      sum=sum
     & -0.5d0*(za(k1,k4)*zb(k1,k2)*zab2(k3,k1,k2,k4))**2
     & /(zb(k1,k4)*zb(k2,k4)*s124)
     & *BGRL2hat(s124,s12)

     & +2d0*za(k3,k4)*zab2(k3,k1,k4,k2)/zb(k1,k4)
     & *BGRL1(s124,s12)

     & +0.5d0*zab2(k3,k1,k4,k2)**2/(zb(k1,k4)*zb(k2,k4)*s124)
     & *lnrat(-s124,-s12)

     & -0.5d0*(za(k1,k4)*zb(k2,k4)*s234)**2
     & /(zb(k2,k3)*zb(k3,k4)*zab2(k1,k2,k3,k4))
     & *BGRL2hat(s234,s23)

     & -2d0*za(k3,k4)*zab2(k1,k3,k4,k2)/zb(k3,k4)
     & *BGRL1(s234,s23)

     & +0.5d0*zab2(k1,k3,k4,k2)**2
     & /(zb(k2,k3)*zb(k3,k4)*zab2(k1,k2,k3,k4))
     & *lnrat(-s234,-s23)

     &-0.5d0*(za(k1,k2)*zb(k1,k2)*zab2(k4,k2,k3,k1))**2*zb(k2,k3)
     & /(zb(k1,k3)**3*s123)
     & *BGRL2hat(s123,s23)

     & +2d0*za(k1,k3)*zb(k1,k2)*zab2(k4,k1,k2,k3)*zab2(k4,k2,k3,k1)
     & /(za(k2,k3)*zb(k1,k3)**2*zb(k2,k3))
     & *BGRL1(s123,s23)

     & +(-2d0*za(k1,k3)*zb(k1,k2)*zab2(k4,k1,k2,k3)*zab2(k4,k2,k3,k1)
     & /(zb(k1,k3)**2*za(k2,k3)*zb(k2,k3)*s123)
     & +0.5d0*zab2(k4,k2,k3,k1)**2*zb(k2,k3)
     & /(zb(k1,k3)**3*s123))
     & *lnrat(-s123,-s23)

     & -0.5d0*(za(k1,k3)*zb(k1,k2)*zab2(k4,k1,k2,k3))**2
     & /(zb(k1,k3)*zb(k2,k3)*s123)
     & *BGRL2hat(s123,s12)

     & +za(k3,k4)*zb(k1,k2)*zab2(k4,k1,k2,k3)
     & *(-2d0*za(k1,k3)*zb(k1,k3)-za(k2,k3)*zb(k2,k3))
     & /(za(k2,k3)*zb(k1,k3)**2*zb(k2,k3))
     & *BGRL1(s123,s12)

     & +zb(k1,k2)*zab2(k4,k1,k2,k3)
     & *(za(k2,k3)*zab2(k4,k1,k3,k2)+2d0*za(k1,k3)*zab2(k4,k2,k3,k1))
     & /(zb(k1,k3)**2*za(k2,k3)*zb(k2,k3)*s123)
     & *lnrat(-s123,-s12)

      A1phiAQggmpmmR=A1phiAQggmpmmR+sum

c--- now add the rational pieces
      sum=
     .-(za(k2,k4)**2*zb(k2,k1)**2)/(2d0*za(k2,k3)*zb(k3,k1)**3)
     .+(zab2(k4,k1,k2,k3)**2*zb(k2,k1)**2)
     .  /(2d0*s123*zb(k3,k1)**3*zb(k3,k2))
     .-(za(k1,k4)**2*zb(k2,k1))/(2*za(k1,k2)*zb(k3,k1)*zb(k3,k2))
     .+(zb(k2,k1)*(za(k1,k3)**2*za(k2,k3)
     .  *zab2(k4,k1,k2,k3)**2*zb(k3,k1)**2
     .+za(k1,k2)**3*zab2(k4,k2,k3,k1)**2*zb(k2,k1)*zb(k3,k2)))
     .  /(4d0*s123**2*za(k1,k2)*za(k2,k3)*zb(k3,k1)**3*zb(k3,k2))
     .+zab2(k3,k1,k4,k2)**2/(2d0*s124*zb(k4,k1)*zb(k4,k2))
     .-(za(k1,k3)**2*zb(k2,k1))/(2d0*za(k1,k2)*zb(k4,k1)*zb(k4,k2))
     .+(za(k1,k4)**2*zab2(k3,k1,k2,k4)**2*zb(k2,k1))
     .  /(4d0*s124**2*za(k1,k2)*zb(k4,k1)*zb(k4,k2))
     .-(za(k1,k3)*za(k1,k4)*zb(k4,k2))/(2d0*zab2(k1,k2,k3,k4)*zb(k4,k3))
     .-(s234*za(k1,k4)**2*zb(k4,k2)**2)
     .  /(4d0*za(k2,k3)*zab2(k1,k2,k3,k4)*zb(k3,k2)**2*zb(k4,k3))
     .-(za(k1,k4)**2*zb(k4,k2)**2)
     .  /(2d0*zab2(k1,k2,k3,k4)*zb(k3,k2)*zb(k4,k3))

      A1phiAQggmpmmR=A1phiAQggmpmmR+sum

      return
      end


      double complex function A1phiAQggmpmmF(k1,k2,k3,k4,za,zb)
c--- This is an implementation of Eq. (5.13) in
c---  S.~Badger, John.~M.~Campbell, R.~Keith Ellis and Ciaran Williams
c---  "Analytic results for the one-loop NMHV H-qbar-q-g-g amplitude."
c---   preprint DESY 09-180, FERMILAB-PUB-09-505-T, IPPP/09/86
c---   arXiv: 0910.4481 [hep-ph]
      implicit none
      include 'constants.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'scale.f'
      include 'epinv.f'
      integer k1,k2,k3,k4
      double complex A0phiAQggmpmm,lnrat
      double complex l12,zab2,BGRL1,BGRL2hat,BGRL3hat
      double precision s12,s13,s14,s23,s24,s123,s124
      zab2(k1,k2,k3,k4)=za(k1,k2)*zb(k2,k4)+za(k1,k3)*zb(k3,k4)

      s12=s(k1,k2)
      s13=s(k1,k3)
      s14=s(k1,k4)
      s23=s(k2,k3)
      s24=s(k2,k4)
      s123=s12+s13+s23
      s124=s12+s14+s24

      l12=lnrat(musq,-s12)

      A1phiAQggmpmmF=A0phiAQggmpmm(k1,k2,k3,k4,za,zb)
     . *(-2d0/3d0*(epinv+l12)-10d0/9d0)
     & +2d0/3d0*za(k1,k3)*zab2(k4,k1,k3,k2)**2
     & /(za(k1,k2)*zb(k1,k2)*zb(k2,k3)*s123)*lnrat(-s123,-s12)
     & -2d0/3d0*(s24-s124)*zab2(k3,k1,k4,k2)**2
     & /(za(k1,k2)*zb(k1,k4)*zb(k2,k4)*zb(k1,k2)*s124)*lnrat(-s124,-s12)
     & -2d0/3d0*za(k1,k3)*zab2(k4,k1,k3,k2)**2
     & /(za(k1,k2)*zb(k2,k3)*zb(k1,k2))*BGRL1(s123,s12)
     & +2d0/3d0*za(k1,k4)*zab2(k3,k1,k4,k2)**2
     & /(za(k1,k2)*zb(k2,k4)*zb(k1,k2))*BGRL1(s124,s12)
     & +za(k1,k3)*za(k3,k4)*zab2(k4,k1,k3,k2)/3d0*BGRL2hat(s123,s12)
     & +za(k1,k4)*za(k3,k4)*zab2(k3,k1,k4,k2)/3d0*BGRL2hat(s124,s12)
     & +2d0/3d0*za(k1,k3)**2*za(k3,k4)*zb(k1,k2)*zab2(k4,k1,k2,k3)
     & *BGRL3hat(s123,s12)
     & +2d0/3d0*za(k1,k4)**2*za(k3,k4)*zb(k1,k2)*zab2(k3,k1,k2,k4)
     & *BGRL3hat(s124,s12)
      A1phiAQggmpmmF=A1phiAQggmpmmF
     & -za(k1,k3)*za(k3,k4)*zab2(k4,k1,k3,k2)
     & /(6d0*za(k1,k2)*zb(k1,k2)*s123)
     & -za(k1,k4)*za(k3,k4)*zab2(k3,k1,k4,k2)
     & /(6d0*za(k1,k2)*zb(k1,k2)*s124)
     & +(
     &  -za(k1,k3)*za(k1,k4)*zb(k1,k2)*zb(k3,k4))
     & /(3d0*za(k1,k2)*zb(k1,k2)*zb(k3,k4)**2)

      return
      end
