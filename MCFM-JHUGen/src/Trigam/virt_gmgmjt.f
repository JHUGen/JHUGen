************************************************************************
*    Author: J. Campbell, January 2014                                 *
*                                                                      *
*    Virtual leading color amplitude for the gamma-gamma-jet process   *
*    q(i1)^- + qb(i2)^+ + gamma(i3)^+ + gamma(i4)^+ + jet(i5)^-        *
*                                                                      *
*    Adapted from the original routines of C. Williams, March 2013     *
*                                                                      *
************************************************************************
      double complex function virt_gmgmjt_gammaMHV(i1,i2,i3,i4,i5,za,zb) 
      implicit none 
      include 'constants.f'
      include 'zprods_decl.f'
      include 'sprods_com.f' 
      include 'scale.f'
      include 'epinv.f'
      integer i1,i2,i3,i4,i5
      double complex l13,l23,l25
      double complex amp_2gam1g
      double complex Vpole,Alo,zab,zab2,Lsm1,L0,L1,lnrat

c--- statement functions
      zab(i1,i2,i3)=+za(i1,i2)*zb(i2,i3)
      zab2(i1,i2,i3,i4)=zab(i1,i2,i4)+zab(i1,i3,i4)

      l13=lnrat(musq,-s(i1,i3))
      l23=lnrat(musq,-s(i2,i3))
      l25=lnrat(musq,-s(i2,i5))
      Alo=amp_2gam1g(i1,i2,i5,i4,i3,za,zb) 

      Vpole=(epinv**2+epinv*l13+0.5d0*l13**2)
     &     +(epinv**2+epinv*l23+0.5d0*l23**2)
     &     +3d0/2d0*(epinv+l25+2d0)
   
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
     & -za(i1,i3)*(za(i4,i5)*zb(i4,i3))**2/(2d0*za(i2,i4)*za(i3,i4))
     &  *L1(-s(i1,i3),-s(i2,i5))/s(i2,i5)**2
     & +(za(i1,i5)*za(i3,i4)-za(i1,i3)*za(i4,i5))
     &  *za(i3,i5)*zb(i4,i3)/(za(i2,i3)*za(i3,i4)**2)
     &  *L0(-s(i1,i4),-s(i2,i5))/s(i2,i5)
     & -za(i1,i4)*(za(i3,i5)*zb(i4,i3))**2
     &  /(2d0*za(i2,i3)*za(i3,i4))
     &  *L1(-s(i1,i4),-s(i2,i5))/s(i2,i5)**2
     & +za(i1,i2)*za(i2,i5)*za(i1,i5)/(za(i1,i3)*za(i2,i3)*za(i2,i4)**2)
     & *lnrat(-s(i4,i5),-s(i2,i5))
     
     & +za(i3,i4)*zb(i3,i4)**2
     &  /(2d0*za(i2,i3)*za(i2,i4)*zb(i1,i5)*zb(i2,i5))
          
c      write(6,*) 'virt_gmgmjt_gammaMHV,virt_gmgmjt_gammaMHV/Alo',
c     & virt_gmgmjt_gammaMHV,virt_gmgmjt_gammaMHV/Alo
     
c      write(6,*) 'log(-s13)',
c     & -za(i1,i2)*za(i2,i5)**2*zb(i3,i2)/(za(i2,i3)*za(i2,i4)**2)
c     &  /(1d0-s(i1,i3)/s(i4,i5))/s(i4,i5)
c     & +(za(i1,i3)*za(i2,i4)+za(i1,i2)*za(i3,i4))*za(i4,i5)**2*zb(i4,i3)
c     &  /(za(i2,i4)*za(i3,i4))**2
c     &  /(1d0-s(i1,i3)/s(i2,i5))/s(i2,i5)
c     & +za(i1,i3)*(za(i4,i5)*zb(i4,i3))**2/(2d0*za(i2,i4)*za(i3,i4))
c     &  /(1d0-s(i1,i3)/s(i2,i5))**2/s(i2,i5)**2
c      write(6,*) 'log(-s45)',
c     & -za(i1,i2)*za(i2,i5)**2*zb(i3,i2)/(za(i2,i3)*za(i2,i4)**2)
c     &  /(1d0-s(i1,i3)/s(i4,i5))/s(i4,i5)*(-1d0)
c     & -za(i1,i2)*za(i2,i5)*za(i1,i5)/(za(i1,i3)*za(i2,i3)*za(i2,i4)**2)
c      write(6,*) 'log(-s14)',
c     & +(-za(i1,i5)*za(i3,i5)*zb(i4,i3)/(2d0*za(i2,i3)*za(i3,i4))
c     &   +za(i1,i3)*za(i3,i5)*za(i4,i5)*zb(i4,i3)
c     &   /(za(i2,i3)*za(i3,i4)**2))
c     &  /(1d0-s(i1,i4)/s(i2,i5))/s(i2,i5)
c     & +za(i1,i3)*za(i2,i5)*za(i3,i5)*zb(i3,i2)*zb(i4,i3)
c     &  /(2d0*za(i2,i3)*za(i3,i4))
c     &  /(1d0-s(i1,i4)/s(i2,i5))**2/s(i2,i5)**2
          
      return 
      end
      
      
      double complex function virt_gmgmjt_gluonMHV(i1,i2,i3,i4,i5,za,zb) 
      implicit none 
      include 'constants.f'
      include 'zprods_decl.f'
      include 'sprods_com.f' 
      include 'scale.f'
      include 'epinv.f'
      integer i1,i2,i3,i4,i5
      double complex l13,l23
      double complex amp_2gam1g
      double complex Vpole,Alo,zab,zab2,Lsm1,L0,L1,lnrat,L1norat

c--- statement functions
      zab(i1,i2,i3)=+za(i1,i2)*zb(i2,i3)
      zab2(i1,i2,i3,i4)=zab(i1,i2,i4)+zab(i1,i3,i4)

      l13=lnrat(musq,-s(i1,i3))
      l23=lnrat(musq,-s(i2,i3))
      Alo=amp_2gam1g(i1,i2,i3,i4,i5,za,zb) 

      Vpole=(epinv**2+epinv*l13+0.5d0*l13**2)
     &     +(epinv**2+epinv*l23+0.5d0*l23**2)
     &     +3d0/2d0*(epinv+l23+2d0)
   
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
     
     & -za(i1,i4)*(za(i3,i5)*zb(i5,i4))**2/(2d0*za(i2,i5)*za(i4,i5))
     &  *L1(-s(i1,i4),-s(i2,i3))/s(i2,i3)**2
     & -za(i1,i3)*za(i3,i5)*zb(i5,i4)/(za(i2,i5)*za(i4,i5))
     &  *L0(-s(i1,i4),-s(i2,i3))/s(i2,i3)
     & +za(i1,i5)*(za(i3,i4)*zb(i5,i4))**2/(2d0*za(i2,i4)*za(i4,i5))
     &  *L1(-s(i1,i5),-s(i2,i3))/s(i2,i3)**2
     & -za(i1,i3)*za(i3,i4)*zb(i5,i4)/(za(i2,i4)*za(i4,i5))
     &  *L0(-s(i1,i5),-s(i2,i3))/s(i2,i3)
      
     & -zb(i4,i5)/(2d0*zb(i1,i3)*zb(i2,i3))
     &  *(zb(i2,i4)/za(i2,i5)-zb(i2,i5)/za(i2,i4))
     & +za(i1,i3)*zb(i4,i5)/(2d0*za(i2,i3)*zb(i2,i3)*za(i4,i5))
     &  *(za(i3,i4)/za(i2,i4)+za(i3,i5)/za(i2,i5))
      
          
c      write(6,*) 'virt_gmgmjt_gluonMHV,virt_gmgmjt_gluonMHV/Alo',
c     & virt_gmgmjt_gluonMHV,virt_gmgmjt_gluonMHV/Alo

      return 
      end
      
      
      double complex function virt_gmgmjt_nfgammaMHV(i1,i2,i3,i4,i5,
     & za,zb)
      implicit none 
      include 'constants.f'
      include 'zprods_decl.f'
      include 'sprods_com.f' 
      include 'scale.f'
      include 'epinv.f'
      integer i1,i2,i3,i4,i5
      double complex Alo,zab,zab2,Lsm1,L0,L1,L1norat

c--- statement functions
      zab(i1,i2,i3)=+za(i1,i2)*zb(i2,i3)
      zab2(i1,i2,i3,i4)=zab(i1,i2,i4)+zab(i1,i3,i4)
   
      virt_gmgmjt_nfgammaMHV=

     & -2d0*((za(i1,i4)*za(i3,i5))**2+(za(i1,i3)*za(i4,i5))**2)
     &  /(za(i1,i2)*za(i3,i4)**4)
     &  *Lsm1(-s(i3,i5),-s(i1,i2),-s(i4,i5),-s(i1,i2))
     
     & -za(i1,i4)*za(i3,i5)*zb(i4,i3)/(za(i1,i2)*za(i3,i4)**3)
     &  *(2d0*za(i1,i4)*za(i3,i5)+4d0*za(i1,i3)*za(i4,i5))
     &  *L0(-s(i1,i2),-s(i3,i5))/s(i3,i5)
     & -2d0*(za(i1,i4)*za(i3,i5))**2*za(i4,i5)*zb(i4,i3)*zb(i5,i4)
     &  /(za(i1,i2)*za(i3,i4)**3)
     &  *L1(-s(i1,i2),-s(i3,i5))/s(i3,i5)**2
     & -za(i1,i3)*za(i4,i5)*zb(i3,i4)/(za(i1,i2)*za(i4,i3)**3)
     &  *(2d0*za(i1,i3)*za(i4,i5)+4d0*za(i1,i4)*za(i3,i5))
     &  *L0(-s(i1,i2),-s(i4,i5))/s(i4,i5)
     & -2d0*(za(i1,i3)*za(i4,i5))**2*za(i3,i5)*zb(i3,i4)*zb(i5,i3)
     &  /(za(i1,i2)*za(i4,i3)**3)
     &  *L1(-s(i1,i2),-s(i4,i5))/s(i4,i5)**2

     & -2d0*za(i3,i5)*za(i4,i5)*zb(i2,i5)**2*zb(i3,i4)
     &  /(za(i3,i4)**3*zb(i1,i2)*zb(i3,i5)*zb(i4,i5))
               
c      write(6,*) 'virt_gmgmjt_nfgammaMHV',virt_gmgmjt_nfgammaMHV

      return 
      end
      
      
      double complex function L1norat(x,y)
      implicit none
      include 'constants.f'
      double precision x,y,denom
      double complex L0
      denom=one-x/y
      L1norat=(L0(x,y)+cone*czip)/dcmplx(denom)
      return
      end

