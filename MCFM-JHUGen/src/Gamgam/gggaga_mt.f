!=====this routine calculates the msq for
!=====gg=>gaga including the effect of the top quark mass
!===== CW Jan 2016
      subroutine gggaga_mt(p,msqgg)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'constants.f'
      include 'nf.f'
      include 'masses.f'
      include 'zprods_decl.f'
      include 'ewcharge.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'debug.f'
      include 'scale.f'
      include 'sprods_com.f'
      real(dp):: p(mxpart,4),msqgg
      complex(dp):: amp_massless(2,2,2,2),amp_mass(2,2,2,2)
      complex(dp):: amp_tot(2,2,2,2)
      integer i,h1,h2,h3,h4
      real(dp):: fac_top,fac_light,fac,Qsum,statfac
      real(dp):: test,msqgggaga
      parameter(statfac=0.5_dp)

      amp_massless(:,:,:,:)=czip
      amp_mass(:,:,:,:)=czip
      msqgg=0._dp
     
      Qsum=Q(1)**2+Q(2)**2+Q(3)**2+Q(4)**2+Q(5)**2
    
      fac_light=Qsum*2._dp*ason2pi*esq
      fac_top=Q(2)**2*2._dp*ason2pi*esq
    
      
      fac=avegg*V*statfac
      
      call spinoru(4,p,za,zb)
      debug=.false.
!      if(debug) then
!         scale=1._dp
!         musq=1._dp
!!=======KC point
!         p(1,4)=-3.0000000000000000_dp
!         p(1,1)= 2.1213203435596424_dp
!         p(1,2)= 1.0606601717798212_dp
!         p(1,3)=1.8371173070873839_dp
!         p(2,4)=-4.8000000000000007_dp
!         p(2,1)=-3.3471659869806940_dp
!         p(2,2)=-1.7645881924740974_dp
!         p(2,3)=-2.9534231607713375_dp
!         p(3,4)=3.7905768352654885_dp
!         p(3,1)=-1.3965283077293904_dp
!         p(3,2)= 3.5239440162638527_dp  
!         p(3,3)= 0.0000000000000000_dp
!         p(4,4)=4.0094231647345122_dp 
!         p(4,1)=2.6223739511504425_dp   
!         p(4,2)= -2.8200159955695767_dp      
!         p(4,3)= 1.1163058536839536_dp
!         call spinorz(4,p,za,zb)
!      endif

      call fill_amp_gggaga_mass(1,2,3,4,za,zb,amp_mass,mt**2)
      call fill_amp_gggaga_mass(1,2,3,4,za,zb,amp_massless,zip)
      amp_tot(:,:,:,:)=fac_light*amp_massless(:,:,:,:)
     &     +fac_top*amp_mass(:,:,:,:)

      
      msqgg=zip
      do h1=1,2
         do h2=1,2
            do h3=1,2
               do h4=1,2
                  msqgg=msqgg+abs(amp_tot(h1,h2,h3,h4))**2 
               enddo
            enddo
         enddo
      enddo


      msqgg=msqgg*fac

      return
      end

      subroutine fill_qcdloop_gggaga(ss,tt,uu,mt2)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'scale.f'
      real(dp):: ss,tt,uu,mt2
      complex(dp):: qlI4,qlI3,qlI2
      complex(dp) D0(3),C0(3),B0(3)
      common/basis_int_gggaga/D0,C0,B0
!$omp threadprivate(/basis_int_gggaga/) 

!     notation:
!     D(1) = Ds13s12
!     D(2) = Ds14s12
!     D(3) = Ds14s13

      D0(1)=qlI4(zip,zip,zip,zip,tt,ss,mt2,mt2,mt2,mt2,musq,0)
      D0(2)=qlI4(zip,zip,zip,zip,uu,ss,mt2,mt2,mt2,mt2,musq,0)
      D0(3)=qlI4(zip,zip,zip,zip,tt,uu,mt2,mt2,mt2,mt2,musq,0)
      
!     notation:
!     C(1) = Cs13
!     C(2) = Cs12
!     C(3) = Cs14

      C0(1)=qlI3(tt,zip,zip,mt2,mt2,mt2,musq,0)
      C0(2)=qlI3(zip,zip,ss,mt2,mt2,mt2,musq,0)
      C0(3)=qlI3(uu,zip,zip,mt2,mt2,mt2,musq,0)

!     notation:
!     B(1) = Cs13
!     B(2) = Cs12
!     B(3) = Cs14

      B0(1)=qlI2(tt,mt2,mt2,musq,0)
      B0(2)=qlI2(ss,mt2,mt2,musq,0)
      B0(3)=qlI2(uu,mt2,mt2,musq,0)

      
      return
      end
      
      subroutine fill_amp_gggaga_massless(i1,i2,i3,i4,za,zb,amp)
      implicit none 
      include 'types.f'
      include 'mxpart.f'
      include 'constants.f'
      include 'zprods_decl.f'
      include 'first.f'
      complex(dp):: D0(3),C0(3),B0(3)
      common/basis_int_gggaga/D0,C0,B0
!$omp threadprivate(/basis_int_gggaga/) 
      
      integer i1,i2,i3,i4
      complex(dp):: amp(2,2,2,2)
      complex(dp) Ls0,qlI4,qlI3
      complex(dp):: kcp

      kcp=im
      
      if(first) then
         first=.false.
         call qlinit()
      endif
      
      amp(2,2,2,2) = zb(i1,i2)*zb(i3,i4)
     &     /(za(i1,i2)*za(i3,i4))



      
      return
      end

                      
      subroutine fill_amp_gggaga_mass(i1,i2,i3,i4,za,zb,amp,mt2)
      implicit none 
      include 'types.f'
      include 'mxpart.f'
       include 'constants.f'
      include 'zprods_decl.f'
      include 'first.f' 
      integer i1,i2,i3,i4
      complex(dp):: amp(2,2,2,2),qlI4,qlI3
      complex(dp):: kcp
      real(dp):: mt2,ss,tt,uu
      complex(dp):: D0(3),C0(3),B0(3)
      common/basis_int_gggaga/D0,C0,B0
!$omp threadprivate(/basis_int_gggaga/) 
      complex(dp):: gggaga_pppp_mt,gggaga_pmpp_mt
      complex(dp):: gggaga_mmpp_mt,gggaga_mpmp_mt
      integer h1,h2,h3,h4
     
      kcp=im
!      mt2=(0.4255266775_dp)**2

      if(first) then
         call qlinit()
         first=.false.
      endif
      amp(:,:,:,:)=czip
      ss=real(za(i1,i2)*zb(i2,i1))
      tt=real(za(i1,i3)*zb(i3,i1))
      uu=real(za(i1,i4)*zb(i4,i1))
      
      call fill_qcdloop_gggaga(ss,tt,uu,mt2)
      
!====== ++++      
      amp(2,2,2,2)=gggaga_pppp_mt(i1,i2,i3,i4,za,zb,mt2)
!======= +-++ and friends
      amp(2,1,2,2)=gggaga_pmpp_mt(i1,i2,i3,i4,za,zb,mt2,.false.)
      amp(1,2,2,2)=gggaga_pmpp_mt(i2,i1,i3,i4,za,zb,mt2,.true.)
      amp(2,2,1,2)=gggaga_pmpp_mt(i4,i3,i2,i1,za,zb,mt2,.false.)
      amp(2,2,2,1)=gggaga_pmpp_mt(i3,i4,i1,i2,za,zb,mt2,.false.)
      
! need to calc      amp(2,2,1,2)=gggaga_pmpp_mt(i4,i3,i1,i2,za,zb,mt2)

!======= --++ and friends 
      amp(1,1,2,2)=gggaga_mmpp_mt(i1,i2,i3,i4,za,zb,mt2)
      amp(1,2,1,2)=gggaga_mpmp_mt(i1,i2,i3,i4,za,zb,mt2,.false.)
      amp(1,2,2,1)=gggaga_mpmp_mt(i1,i2,i4,i3,za,zb,mt2,.true.)
      
      
!=======amplitudes obtained by conjg 
      amp(1,1,1,1) =gggaga_pppp_mt(i1,i2,i3,i4,zb,za,mt2)
      amp(1,2,1,1) =gggaga_pmpp_mt(i1,i2,i3,i4,zb,za,mt2,.false.)
      amp(2,1,1,1) =gggaga_pmpp_mt(i2,i1,i3,i4,zb,za,mt2,.true.)
      amp(2,2,1,1) =gggaga_mmpp_mt(i1,i2,i3,i4,zb,za,mt2)
      amp(2,1,2,1) =gggaga_mpmp_mt(i1,i2,i3,i4,zb,za,mt2,.false.)
      amp(2,1,1,2) =gggaga_mpmp_mt(i1,i2,i4,i3,zb,za,mt2,.true.)
      amp(1,1,2,1) =gggaga_pmpp_mt(i4,i3,i2,i1,zb,za,mt2,.false.)
      amp(1,1,1,2) =gggaga_pmpp_mt(i3,i4,i1,i2,zb,za,mt2,.false.)
    
     
 !     write(6,*) '******** massive test ********'
 !     do h1=1,2
 !        do h2=1,2
 !           do h3=1,2
 !              do h4=1,2
 !                
 !      write(6,*) h1,h2,h3,h4,amp(h1,h2,h3,h4)*kcp,abs(amp(h1,h2,h3,h4))
 !     enddo
 !     enddo
 !     enddo
 !     enddo
 !     write(6,*) '*******************************'

!      stop
      
      
      return
      end

      function gggaga_pppp_mt(i1,i2,i3,i4,za,zb,mt2)
      implicit none
      include 'types.f'
      complex(dp) :: gggaga_pppp_mt
      integer i1,i2,i3,i4
      real(dp) ::mt2
 
      include 'mxpart.f'
      include 'constants.f'
      include 'zprods_decl.f'
      complex(dp):: D0(3),C0(3),B0(3)
      common/basis_int_gggaga/D0,C0,B0
!$omp threadprivate(/basis_int_gggaga/) 

      gggaga_pppp_mt=
     &     -mt2**2*(2*zb(i1,i2)*zb(i3,i4))/(za(i1,i2)*za(i3,i4))*(
     &   D0(1)+D0(2)+D0(3))
     &    +zb(i1,i2)*zb(i3,i4)/(za(i1,i2)*za(i3,i4))

      return
      end

!========box coeffs

      function gggaga_mmpp_mt(i1,i2,i3,i4,za,zb,mt2)
      implicit none
      include 'types.f'
      complex(dp) :: gggaga_mmpp_mt
      integer i1,i2,i3,i4
      real(dp) ::mt2
      include 'mxpart.f'
      include 'constants.f'
      include 'zprods_decl.f'
      complex(dp):: D0(3),C0(3),B0(3)
      common/basis_int_gggaga/D0,C0,B0
!$omp threadprivate(/basis_int_gggaga/) 
      complex(dp):: d1coeff,d2coeff,d3coeff
      complex(dp):: c1coeff,c2coeff,c3coeff
      complex(dp):: b1coeff,b2coeff,b3coeff,rat

!========box coeffs
      d1coeff= (-2*mt2**2*za(i1,i2)**2*zb(i3,i2)*zb(i4,i1))/
     -   (za(i1,i4)*za(i2,i3)*zb(i2,i1)**2) + 
     -  (mt2*za(i1,i2)**2*zb(i4,i1)*zb(i4,i3))/
     -   (za(i2,i3)*zb(i2,i1))
      d2coeff=(-2*mt2**2*za(i1,i2)**2*zb(i3,i2)*zb(i4,i1))/
     -   (za(i1,i4)*za(i2,i3)*zb(i2,i1)**2) + 
     -  (mt2*za(i1,i2)**2*zb(i4,i1)*zb(i4,i3))/
     -   (za(i2,i3)*zb(i2,i1))
      d3coeff=(zb(i4,i1)*(-8*mt2*za(i1,i2)*za(i1,i4)*za(i2,i4)*
     -       zb(i2,i1)*zb(i3,i2)*zb(i4,i1)*zb(i4,i2) + 
     -      za(i1,i4)**2*za(i2,i4)*zb(i4,i1)*zb(i4,i2)*
     -       (za(i1,i4)*zb(i3,i2)*zb(i4,i1)**2 - 
     -         za(i2,i4)*zb(i3,i1)*zb(i4,i2)**2) + 
     -      2*mt2*za(i1,i2)**2*zb(i2,i1)**2*
     -       (-2*mt2*zb(i3,i2) + za(i1,i4)*zb(i2,i1)*zb(i4,i3))))
     -   /(2.*za(i1,i4)*za(i2,i3)*zb(i2,i1)**4)
      
!=======tri coeffs
      c1coeff=(-4*mt2*za(i1,i2)*zb(i3,i1)*zb(i4,i2))/
     -   (za(i3,i4)*zb(i2,i1)**2) + 
     -  (za(i1,i3)*zb(i3,i1)*
     -     (-(za(i2,i3)*zb(i3,i2)**2*zb(i4,i1)) + 
     -       za(i1,i3)*zb(i3,i1)**2*zb(i4,i2)))/
     -     (za(i3,i4)*zb(i2,i1)**3)
      c2coeff=czip
      c3coeff= (4*mt2*za(i1,i2)*zb(i3,i2)*zb(i4,i1))/
     -   (za(i3,i4)*zb(i2,i1)**2) - 
     -  (za(i1,i4)*zb(i4,i1)*
     -     (za(i1,i4)*zb(i3,i2)*zb(i4,i1)**2 - 
     -       za(i2,i4)*zb(i3,i1)*zb(i4,i2)**2))/
     -   (za(i3,i4)*zb(i2,i1)**3)
!=======bubs

      b1coeff=-(za(i2,i4)**2*(za(i1,i4)*zb(i4,i1) - 
     -      za(i2,i4)*zb(i4,i2))*zb(i4,i3))/
     -     (za(i3,i4)**3*zb(i3,i1)**2)
      b2coeff=czip
      b3coeff=-b1coeff
!====== rat

      rat= -(za(i1,i2)**2*zb(i3,i2)*zb(i4,i1))/
     -     (za(i1,i4)*za(i2,i3)*zb(i2,i1)**2)
      
      gggaga_mmpp_mt=d1coeff*D0(1)+d2coeff*D0(2)+d3coeff*D0(3)
     &     +c1coeff*C0(1)+c2coeff*C0(2)+c3coeff*C0(3)
     &     +b1coeff*B0(1)+b2coeff*B0(2)+b3coeff*B0(3)
     &     +rat
      
      return
      end

      function gggaga_mpmp_mt(i1,i2,i3,i4,za,zb,mt2,swapi)
      implicit none
      include 'types.f'
      complex(dp) :: gggaga_mpmp_mt
      integer i1,i2,i3,i4
      real(dp) ::mt2
      include 'mxpart.f'
      include 'constants.f'
      include 'zprods_decl.f'
      complex(dp):: D0(3),C0(3),B0(3)
      common/basis_int_gggaga/D0,C0,B0
!$omp threadprivate(/basis_int_gggaga/) 
      complex(dp):: d1coeff,d2coeff,d3coeff
      complex(dp):: c1coeff,c2coeff,c3coeff
      complex(dp):: b1coeff,b2coeff,b3coeff,rat,temp
      logical swapi
      if(swapi) then
         temp=D0(1)
         D0(1)=D0(2)
         D0(2)=temp
         temp=C0(1)
         C0(1)=C0(3)
         C0(3)=temp
         temp=B0(1)
         B0(1)=B0(3)
         B0(3)=temp
      endif
         
      
      
!========box coeffs
      d1coeff= -0.5_dp*((2*mt2*za(i1,i3)*za(i1,i4)*zb(i4,i2)**2)/
     -   (za(i2,i4)*zb(i3,i2)) + 
     -  (4*mt2**2*za(i1,i4)*za(i2,i3)*zb(i4,i2)**2)/
     -   (za(i2,i4)**2*zb(i3,i2)*zb(i4,i1)))
      d2coeff=-0.5_dp*((za(i1,i2)*za(i2,i3)*zb(i2,i1)*zb(i3,i2)*
     -     (za(i2,i3)*zb(i3,i2)**2*zb(i4,i1) + 
     -       za(i1,i2)*zb(i2,i1)**2*zb(i4,i3)) + 
     -    2*mt2*za(i1,i3)*zb(i3,i1)*
     -     (2*mt2*zb(i3,i1)*zb(i4,i2) + 
     -       za(i1,i2)*zb(i2,i1)**2*zb(i4,i3) + 
     -       za(i2,i3)*zb(i3,i2)*
     -        (zb(i3,i2)*zb(i4,i1) - 2*zb(i2,i1)*zb(i4,i3))))/
     -  (za(i2,i4)*zb(i3,i1)**3))
      d3coeff=-0.5_dp*((2*mt2*za(i1,i3)*za(i1,i4)*zb(i4,i2)**2)/
     -   (za(i2,i4)*zb(i3,i2)) + 
     -  (4*mt2**2*za(i1,i4)*za(i2,i3)*zb(i4,i2)**2)/
     -   (za(i2,i4)**2*zb(i3,i2)*zb(i4,i1)))
!=======tri coeffs
      c1coeff=czip
      c2coeff=(4*mt2*za(i1,i2)*za(i1,i3)*zb(i2,i1)*zb(i3,i2)*
     -     zb(i4,i1))/(za(i1,i4)*za(i2,i3)*zb(i3,i1)**3) + 
     -  (za(i1,i2)*zb(i2,i1)*zb(i3,i2)*
     -     (-(za(i1,i4)*zb(i3,i2)*zb(i4,i1)**2) - 
     -       za(i1,i2)*zb(i2,i1)**2*zb(i4,i3)))/
     -     (za(i1,i4)*zb(i3,i1)**4)
      
      c3coeff= (4*mt2*za(i1,i3)*za(i1,i4)*zb(i2,i1)*zb(i4,i1)*
     -     zb(i4,i3))/(za(i1,i2)*za(i3,i4)*zb(i3,i1)**3) - 
     -  (za(i1,i4)*zb(i4,i1)*zb(i4,i3)*
     -     (za(i1,i4)*zb(i3,i2)*zb(i4,i1)**2 + 
     -       za(i1,i2)*zb(i2,i1)**2*zb(i4,i3)))/
     -   (za(i1,i2)*zb(i3,i1)**4)
!=======bubs

      b1coeff=czip
      
      b2coeff=(za(i1,i4)*zb(i4,i2)*
     -    (2*za(i1,i2)*za(i3,i4)*zb(i3,i2) + 
     -      za(i1,i4)*za(i2,i4)*zb(i4,i2)))/
     -  (za(i2,i4)**3*zb(i3,i2)**2)
      b3coeff=-b2coeff
!====== rat

      
      rat=(2*za(i1,i2)*za(i1,i3)*zb(i2,i1)*zb(i3,i2)*zb(i4,i1))/
     -   (za(i1,i4)*za(i2,i3)*zb(i3,i1)**3) + 
     -  (za(i1,i4)*za(i2,i3)*zb(i4,i2)**2)/
     -   (za(i2,i4)**2*zb(i3,i2)*zb(i4,i1)) + 
     -  (2*za(i1,i3)*za(i1,i4)*zb(i2,i1)*zb(i4,i1)*zb(i4,i3))/
     -     (za(i1,i2)*za(i3,i4)*zb(i3,i1)**3)

      
 !     write(6,*) d1coeff*im
 !     write(6,*) d2coeff*im
 !     write(6,*) d3coeff*im
 !     
 !     write(6,*) c1coeff*im
 !     write(6,*) c2coeff*im
 !!     write(6,*) c3coeff*im
!      
!      write(6,*) b1coeff*im
!      write(6,*) b2coeff*im
!      write(6,*) b3coeff*im
!      
!      write(6,*) rat*im
!      pause
      gggaga_mpmp_mt=d1coeff*D0(1)+d2coeff*D0(2)+d3coeff*D0(3)
     &     +c1coeff*C0(1)+c2coeff*C0(2)+c3coeff*C0(3)
     &     +b1coeff*B0(1)+b2coeff*B0(2)+b3coeff*B0(3)
     &     +rat

      if(swapi) then
         temp=D0(1)
         D0(1)=D0(2)
         D0(2)=temp
         temp=C0(1)
         C0(1)=C0(3)
         C0(3)=temp
         temp=B0(1)
         B0(1)=B0(3)
         B0(3)=temp
      endif
      
      return
      end
      
      function gggaga_pmpp_mt(i1,i2,i3,i4,za,zb,mt2,swapi)
      implicit none
      include 'types.f'
      complex(dp) :: gggaga_pmpp_mt
      integer i1,i2,i3,i4
      real(dp) ::mt2
      include 'mxpart.f'
      include 'constants.f'
      include 'zprods_decl.f'
      complex(dp):: D0(3),C0(3),B0(3)
      common/basis_int_gggaga/D0,C0,B0
!$omp threadprivate(/basis_int_gggaga/) 
      logical swapi
      complex temp
      
      if(swapi) then
         temp=D0(1)
         D0(1)=D0(2)
         D0(2)=temp
         temp=C0(1)
         C0(1)=C0(3)
         C0(3)=temp
         temp=B0(1)
         B0(1)=B0(3)
         B0(3)=temp
      endif

      gggaga_pmpp_mt=
     &  D0(1)*((mt2*za(i1,i2)**2*zb(i3,i1)**2)/za(i1,i4)**2 - 
     &  (2*mt2**2*zb(i3,i1)*zb(i4,i1)*zb(i4,i3))/
     &   (za(i1,i4)*zb(i2,i1)*zb(i4,i2)))
     &   +D0(2)*( (mt2*za(i1,i2)**2*zb(i4,i1)**2)/za(i1,i3)**2 + 
     -  (2*mt2**2*zb(i3,i1)*zb(i4,i1)*zb(i4,i3))/
     -   (za(i1,i3)*zb(i2,i1)*zb(i3,i2)))
     & +D0(3)*((2*mt2**2*za(i1,i2)*za(i2,i3)*zb(i3,i1)**2*zb(i4,i2))/
     -   (za(i1,i3)**2*za(i2,i4)*zb(i2,i1)*zb(i3,i2)) - 
     -  (mt2*za(i2,i3)*zb(i3,i1)**2*zb(i4,i1)*zb(i4,i2))/
     -   (za(i1,i3)*zb(i2,i1)**2))
     & +C0(1)*((-2*mt2*zb(i3,i1)*(za(i1,i2)*za(i1,i3)*za(i3,i4)*
     &       zb(i3,i1) - za(i1,i4)**2*za(i2,i3)*zb(i4,i1)))/
     &  (za(i1,i3)*za(i1,i4)**2*za(i3,i4)*zb(i2,i1)))
     & +C0(2)*((2*mt2*za(i1,i2)**3*zb(i2,i1))/
     &   (za(i1,i3)**2*za(i1,i4)**2) - 
     &  (2*mt2*za(i1,i2)*zb(i3,i1)*zb(i4,i1))/
     &   (za(i1,i3)*za(i1,i4)*zb(i2,i1)))
     &  +C0(3)*((-2*mt2*(za(i1,i2)**3*za(i3,i4)*zb(i2,i1)**2 + 
     &      za(i1,i2)**2*za(i1,i3)*za(i3,i4)*zb(i2,i1)*
     &       zb(i3,i1) + za(i1,i3)*za(i1,i4)**2*za(i2,i3)*
     &       zb(i3,i1)*zb(i4,i1)))/
     &  (za(i1,i3)**2*za(i1,i4)**2*za(i3,i4)*zb(i2,i1)))      
     &     -((za(i1,i2)*za(i2,i3)*zb(i3,i1)**2*zb(i4,i2))/
     &     (za(i1,i3)**2*za(i2,i4)*zb(i2,i1)*zb(i3,i2)))

      if(swapi) then
         temp=D0(1)
         D0(1)=D0(2)
         D0(2)=temp
         temp=C0(1)
         C0(1)=C0(3)
         C0(3)=temp
         temp=B0(1)
         B0(1)=B0(3)
         B0(3)=temp
      endif
      return
      end
