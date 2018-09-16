
      function qqbgWH_topamp(i1,i2,i3,i4,i5,za,zb,mt2) 
!===== C.Williams July 2015 
!====== amplitude for q(i1)^-+qb(i2)^+ + ell^-(i3)+ell^+(i4)+g(i5) + Higgs 
!===== where the Higgs is radiated from the top loop 
      implicit none 
      include 'types.f' 
      complex(dp):: qqbgWH_topamp
      include 'constants.f'
      include 'qcdcouple.f' 
      include 'ewcouple.f' 
      include 'mxpart.f'
      include 'scale.f'
      include 'masses.f'
      include 'zprods_decl.f' 
      include 'sprods_com.f' 
      integer:: i1,i2,i3,i4,i5
!====== basis integrals B0mt(1) = Bub(s12345,mt) B0mt(2)=Bub(s1234,mt) 
      complex(dp)::ffDDHK_ql
      real(dp)::s12345,s1234,s134,s234,s34,xd
      complex(dp):: zab2,zab3
      complex(dp):: helpart,intfunc,fac,prop34
      complex(dp):: tri,bmH,bmgH
      complex(dp):: qlI3,qlI2
      real(dp):: mt2
      logical useeft_wh
      common/useeft_wh/useeft_wh
!$omp threadprivate(/useeft_wh/) 

 !====== statement functions 
      include 'cplx.h'
 
      zab2(i1,i2,i3,i4)=za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)
      zab3(i1,i2,i3,i5,i4)=za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)
     &     +za(i1,i5)*zb(i5,i4)

!====== useful invariants 
      s34=s(i3,i4)
      s134=s(i1,i3)+s(i3,i4)+s(i1,i4) 
      s234=s(i2,i3)+s(i2,i4)+s(i3,i4) 
      s1234=s134+s234+s(i1,i2)-s(i3,i4)
!====== s12345 is the Higgs mass 
      s12345=s(i1,i5)+s(i2,i5)+s(i3,i5)+s(i4,i5)+s1234

      prop34=cplx1(s34)/cplx2(s34-wmass**2,wmass*wwidth)
      
      if(useeft_wh) then 
      fac=-im*sqrt(gsq)*gwsq*prop34*as/(3_dp*pi*sqrt(vevsq))/rt2
      else
      fac=im*sqrt(gsq)*ason4pi*gwsq*prop34*mt2/wmass*gw/rt2/2_dp
      endif

!=== Spin independent function this is just the form factor 
!==== and can be written in terms of existing MCFM functions  
      
 !     intfunc=-ason4pi*sqrt(esq/xw)/wmass
 !    &     *ffDDHK_ql(s12345,s1234,mt2)*s1234*(four*xn*mt2/s1234)

!===== helicity dependent part (includes 1/s1234 from gluon propagator) 

      helpart=((s134*zab2(i3,i2,i4,i5)*zab3(i1,i2,i3,i4,i5)*zb(i4,i2)
     &     +s234*za(i1,i3)*zb(i5,i2)
     & *(za(i1,i2)*zb(i4,i1)*zb(i5,i2)-za(i2,i3)*zb(i4,i3)*zb(i5,i2)
     &  +s134*zb(i5,i4))))/(s1234*s134*s234*s34)

      xd=one/(s12345-s1234)
  
      if(useeft_wh) then 
      intfunc=1_dp
      else
         tri=qlI3(s12345,s1234,0_dp,mt2,mt2,mt2,musq,0) 
         bmgH=qlI2(s1234,mt2,mt2,musq,0) 
         bmH=qlI2(s12345,mt2,mt2,musq,0) 
         intfunc=8_dp*((0.5_dp*(one-4_dp*mt2*xd)*tri
     &        + s1234*xd**2*(bmgH-bmH)-xd))
      endif
!      write(6,*) 'tri coefficient ',0.5_dp*(one-4_dp*mt2*xd)*helpart*8
!      write(6,*) 'bub(mH) coefficient ',-s1234*xd**2*helpart*8_dp
!      write(6,*) 'bub(mHg) coefficient ',s1234*xd**2*helpart*8_dp
!      write(6,*) 'rational piece ',-xd*helpart*8_dp
!      write(6,*) 'total ',intfunc*helpart
!      stop

!===== total       
      qqbgWH_topamp=helpart*intfunc*fac

      return 
      end

