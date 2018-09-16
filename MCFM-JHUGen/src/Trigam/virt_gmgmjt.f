************************************************************************
*    Author: J. Campbell, January 2014                                 *
*                                                                      *
*    Virtual leading color amplitude for the gamma-gamma-jet process   *
*    q(i1)^- + qb(i2)^+ + gamma(i3)^+ + gamma(i4)^+ + jet(i5)^-        *
*                                                                      *
*    Adapted from the original routines of C. Williams, March 2013     *
*                                                                      *
************************************************************************
      function virt_gmgmjt_gammaMHV(i1,i2,i3,i4,i5,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: virt_gmgmjt_gammaMHV

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'scale.f'
      include 'epinv.f'
      integer:: i1,i2,i3,i4,i5
      complex(dp):: l13,l23,l25
      complex(dp):: amp_2gam1g
      complex(dp):: Vpole,Alo,zab,zab2,Lsm1,L0,L1,lnrat

c--- statement functions
      zab(i1,i2,i3)=+za(i1,i2)*zb(i2,i3)
      zab2(i1,i2,i3,i4)=zab(i1,i2,i4)+zab(i1,i3,i4)

      l13=lnrat(musq,-s(i1,i3))
      l23=lnrat(musq,-s(i2,i3))
      l25=lnrat(musq,-s(i2,i5))
      Alo=amp_2gam1g(i1,i2,i5,i4,i3,za,zb)

      Vpole=(epinv**2+epinv*l13+0.5_dp*l13**2)
     &     +(epinv**2+epinv*l23+0.5_dp*l23**2)
     &     +3._dp/2._dp*(epinv+l25+2._dp)

      virt_gmgmjt_gammaMHV=
     & +Vpole*Alo

     & -za(i1,i2)**3*za(i4,i5)**2
     &  /(za(i1,i3)*za(i1,i4)*za(i2,i3)*za(i2,i4)**3)
     &  *Lsm1(-s(i4,i5),-s(i1,i3),-s(i2,i5),-s(i1,i3))
     & +za(i1,i5)**2/(za(i1,i3)*za(i3,i4)*za(i2,i4))
     &  *Lsm1(-s(i2,i4),-s(i1,i5),-s(i2,i3),-s(i1,i5))
     & -za(i1,i2)*za(i1,i5)**2/(za(i1,i3)*za(i1,i4)*za(i2,i3)*za(i2,i4))
     &  *Lsm1(-s(i4,i5),-s(i2,i3),-s(i1,i5),-s(i2,i3))
     & -(za(i1,i3)*za(i4,i5))**2
     &  /(za(i1,i4)*za(i2,i3)*za(i3,i4)**3)
     &  *Lsm1(-s(i1,i4),-s(i2,i5),-s(i1,i3),-s(i2,i5))
     & -za(i1,i2)*za(i1,i5)**2/(za(i1,i3)*za(i1,i4)*za(i2,i3)*za(i2,i4))
     &  *Lsm1(-s(i1,i3),-s(i4,i5),-s(i2,i3),-s(i4,i5))

     & +za(i1,i2)*za(i2,i5)**2*zb(i3,i2)/(za(i2,i3)*za(i2,i4)**2)
     &  *L0(-s(i1,i3),-s(i4,i5))/s(i4,i5)
     & -(za(i1,i3)*za(i2,i4)+za(i1,i2)*za(i3,i4))*za(i4,i5)**2*zb(i4,i3)
     &  /(za(i2,i4)*za(i3,i4))**2
     &  *L0(-s(i1,i3),-s(i2,i5))/s(i2,i5)
     & -za(i1,i3)*(za(i4,i5)*zb(i4,i3))**2/(2._dp*za(i2,i4)*za(i3,i4))
     &  *L1(-s(i1,i3),-s(i2,i5))/s(i2,i5)**2
     & +(za(i1,i5)*za(i3,i4)-za(i1,i3)*za(i4,i5))
     &  *za(i3,i5)*zb(i4,i3)/(za(i2,i3)*za(i3,i4)**2)
     &  *L0(-s(i1,i4),-s(i2,i5))/s(i2,i5)
     & -za(i1,i4)*(za(i3,i5)*zb(i4,i3))**2
     &  /(2._dp*za(i2,i3)*za(i3,i4))
     &  *L1(-s(i1,i4),-s(i2,i5))/s(i2,i5)**2
     & +za(i1,i2)*za(i2,i5)*za(i1,i5)/(za(i1,i3)*za(i2,i3)*za(i2,i4)**2)
     & *lnrat(-s(i4,i5),-s(i2,i5))

     & +za(i3,i4)*zb(i3,i4)**2
     &  /(2._dp*za(i2,i3)*za(i2,i4)*zb(i1,i5)*zb(i2,i5))

c      write(6,*) 'virt_gmgmjt_gammaMHV,virt_gmgmjt_gammaMHV/Alo',
c     & virt_gmgmjt_gammaMHV,virt_gmgmjt_gammaMHV/Alo

c      write(6,*) 'log(-s13)',
c     & -za(i1,i2)*za(i2,i5)**2*zb(i3,i2)/(za(i2,i3)*za(i2,i4)**2)
c     &  /(1._dp-s(i1,i3)/s(i4,i5))/s(i4,i5)
c     & +(za(i1,i3)*za(i2,i4)+za(i1,i2)*za(i3,i4))*za(i4,i5)**2*zb(i4,i3)
c     &  /(za(i2,i4)*za(i3,i4))**2
c     &  /(1._dp-s(i1,i3)/s(i2,i5))/s(i2,i5)
c     & +za(i1,i3)*(za(i4,i5)*zb(i4,i3))**2/(2._dp*za(i2,i4)*za(i3,i4))
c     &  /(1._dp-s(i1,i3)/s(i2,i5))**2/s(i2,i5)**2
c      write(6,*) 'log(-s45)',
c     & -za(i1,i2)*za(i2,i5)**2*zb(i3,i2)/(za(i2,i3)*za(i2,i4)**2)
c     &  /(1._dp-s(i1,i3)/s(i4,i5))/s(i4,i5)*(-1._dp)
c     & -za(i1,i2)*za(i2,i5)*za(i1,i5)/(za(i1,i3)*za(i2,i3)*za(i2,i4)**2)
c      write(6,*) 'log(-s14)',
c     & +(-za(i1,i5)*za(i3,i5)*zb(i4,i3)/(2._dp*za(i2,i3)*za(i3,i4))
c     &   +za(i1,i3)*za(i3,i5)*za(i4,i5)*zb(i4,i3)
c     &   /(za(i2,i3)*za(i3,i4)**2))
c     &  /(1._dp-s(i1,i4)/s(i2,i5))/s(i2,i5)
c     & +za(i1,i3)*za(i2,i5)*za(i3,i5)*zb(i3,i2)*zb(i4,i3)
c     &  /(2._dp*za(i2,i3)*za(i3,i4))
c     &  /(1._dp-s(i1,i4)/s(i2,i5))**2/s(i2,i5)**2

      return
      end


      function virt_gmgmjt_gluonMHV(i1,i2,i3,i4,i5,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: virt_gmgmjt_gluonMHV

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'scale.f'
      include 'epinv.f'
      integer:: i1,i2,i3,i4,i5
      complex(dp):: l13,l23
      complex(dp):: amp_2gam1g
      complex(dp):: Vpole,Alo,zab,zab2,Lsm1,L0,L1,lnrat,L1norat

c--- statement functions
      zab(i1,i2,i3)=+za(i1,i2)*zb(i2,i3)
      zab2(i1,i2,i3,i4)=zab(i1,i2,i4)+zab(i1,i3,i4)

      l13=lnrat(musq,-s(i1,i3))
      l23=lnrat(musq,-s(i2,i3))
      Alo=amp_2gam1g(i1,i2,i3,i4,i5,za,zb)

      Vpole=(epinv**2+epinv*l13+0.5_dp*l13**2)
     &     +(epinv**2+epinv*l23+0.5_dp*l23**2)
     &     +3._dp/2._dp*(epinv+l23+2._dp)

      virt_gmgmjt_gluonMHV=
     & +Vpole*Alo

     & +za(i1,i3)**2/(za(i1,i4)*za(i2,i5)*za(i4,i5))
     &  *Lsm1(-s(i2,i5),-s(i1,i4),-s(i2,i3),-s(i1,i4))
     & -za(i1,i3)**2/(za(i1,i5)*za(i2,i4)*za(i4,i5))
     &  *Lsm1(-s(i2,i4),-s(i1,i5),-s(i2,i3),-s(i1,i5))
     & -za(i1,i3)**2/(za(i1,i5)*za(i2,i4)*za(i4,i5))
     &  *Lsm1(-s(i1,i5),-s(i2,i4),-s(i1,i3),-s(i2,i4))
     & +za(i1,i3)**2/(za(i1,i4)*za(i2,i5)*za(i4,i5))
     &  *Lsm1(-s(i1,i4),-s(i2,i5),-s(i1,i3),-s(i2,i5))

     & -za(i1,i4)*(za(i3,i5)*zb(i5,i4))**2/(2._dp*za(i2,i5)*za(i4,i5))
     &  *L1(-s(i1,i4),-s(i2,i3))/s(i2,i3)**2
     & -za(i1,i3)*za(i3,i5)*zb(i5,i4)/(za(i2,i5)*za(i4,i5))
     &  *L0(-s(i1,i4),-s(i2,i3))/s(i2,i3)
     & +za(i1,i5)*(za(i3,i4)*zb(i5,i4))**2/(2._dp*za(i2,i4)*za(i4,i5))
     &  *L1(-s(i1,i5),-s(i2,i3))/s(i2,i3)**2
     & -za(i1,i3)*za(i3,i4)*zb(i5,i4)/(za(i2,i4)*za(i4,i5))
     &  *L0(-s(i1,i5),-s(i2,i3))/s(i2,i3)

     & -zb(i4,i5)/(2._dp*zb(i1,i3)*zb(i2,i3))
     &  *(zb(i2,i4)/za(i2,i5)-zb(i2,i5)/za(i2,i4))
     & +za(i1,i3)*zb(i4,i5)/(2._dp*za(i2,i3)*zb(i2,i3)*za(i4,i5))
     &  *(za(i3,i4)/za(i2,i4)+za(i3,i5)/za(i2,i5))


!      write(6,*) 'virt_gmgmjt_gluonMHV,virt_gmgmjt_gluonMHV/Alo',
!     & virt_gmgmjt_gluonMHV,virt_gmgmjt_gluonMHV/Alo

      return
      end

!----rational amplitude not needed for NLO, needed for N^3LO gg
      function virt_gmgmjt_nfallp(i1,i2,i3,i4,i5,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: virt_gmgmjt_nfallp
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'scale.f'
      include 'epinv.f'
      integer:: i1,i2,i3,i4,i5
      complex(dp):: Alo,zab,zab2,Lsm1,L0,L1,L1norat

c--- statement functions
      zab(i1,i2,i3)=+za(i1,i2)*zb(i2,i3)
      zab2(i1,i2,i3,i4)=zab(i1,i2,i4)+zab(i1,i3,i4)

      virt_gmgmjt_nfallp= (2*(za(i2,i3)*za(i4,i5)*zb(i4,i1)*zb(i5,i3)-
     -      za(i2,i4)*za(i3,i5)*zb(i3,i1)*zb(i5,i4)))/
     -  (za(i1,i2)*za(i3,i4)*za(i3,i5)*za(i4,i5)*zb(i2,i1))


      return
      end


      function virt_gmgmjt_nfgammaMHV(i1,i2,i3,i4,i5,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: virt_gmgmjt_nfgammaMHV
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'scale.f'
      include 'epinv.f'
      integer:: i1,i2,i3,i4,i5
      complex(dp):: Alo,zab,zab2,Lsm1,L0,L1,L1norat

c--- statement functions
      zab(i1,i2,i3)=+za(i1,i2)*zb(i2,i3)
      zab2(i1,i2,i3,i4)=zab(i1,i2,i4)+zab(i1,i3,i4)

      virt_gmgmjt_nfgammaMHV=

     & -2._dp*((za(i1,i4)*za(i3,i5))**2+(za(i1,i3)*za(i4,i5))**2)
     &  /(za(i1,i2)*za(i3,i4)**4)
     &  *Lsm1(-s(i3,i5),-s(i1,i2),-s(i4,i5),-s(i1,i2))

     & -za(i1,i4)*za(i3,i5)*zb(i4,i3)/(za(i1,i2)*za(i3,i4)**3)
     &  *(2._dp*za(i1,i4)*za(i3,i5)+4._dp*za(i1,i3)*za(i4,i5))
     &  *L0(-s(i1,i2),-s(i3,i5))/s(i3,i5)
     & -2._dp*(za(i1,i4)*za(i3,i5))**2*za(i4,i5)*zb(i4,i3)*zb(i5,i4)
     &  /(za(i1,i2)*za(i3,i4)**3)
     &  *L1(-s(i1,i2),-s(i3,i5))/s(i3,i5)**2
     & -za(i1,i3)*za(i4,i5)*zb(i3,i4)/(za(i1,i2)*za(i4,i3)**3)
     &  *(2._dp*za(i1,i3)*za(i4,i5)+4._dp*za(i1,i4)*za(i3,i5))
     &  *L0(-s(i1,i2),-s(i4,i5))/s(i4,i5)
     & -2._dp*(za(i1,i3)*za(i4,i5))**2*za(i3,i5)*zb(i3,i4)*zb(i5,i3)
     &  /(za(i1,i2)*za(i4,i3)**3)
     &  *L1(-s(i1,i2),-s(i4,i5))/s(i4,i5)**2

     & -2._dp*za(i3,i5)*za(i4,i5)*zb(i2,i5)**2*zb(i3,i4)
     &  /(za(i3,i4)**3*zb(i1,i2)*zb(i3,i5)*zb(i4,i5))

c      write(6,*) 'virt_gmgmjt_nfgammaMHV',virt_gmgmjt_nfgammaMHV

      return
      end


      function L1norat(x,y)
      implicit none
      include 'types.f'
      complex(dp):: L1norat
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      real(dp):: x,y,denom
      complex(dp):: L0
      denom=one-x/y
      L1norat=(L0(x,y)+cone*czip)/cplx1(denom)
      return
      end

