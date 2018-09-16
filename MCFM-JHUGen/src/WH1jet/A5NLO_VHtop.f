!===== C.Williams July 15, routine to construct interference between 
!===== top-loop mediated diagrams and tree level in VH +jet production
      subroutine A5NLO_VHtop(i1,i2,i3,i4,i5,za,zb,A5LO,A5NLOt) 
      implicit none 
      include 'types.f'
      include 'mxpart.f' 
      include 'masses.f' 
      include 'scale.f' 
      include 'couple.f' 
      include 'zprods_decl.f' 
!-==== this routine is really just a wrapper       
      integer:: i1,i2,i3,i4,i5
      complex(dp):: A5LO,A5NLOt
      complex(dp)::qqbgWH_topamp,Alo5_VHt
      external qqbgWH_topamp,Alo5_VHt
      real(dp):: mt2,mrun2,massfrun
      integer,parameter:: nloop=1
      logical useeft_wh
      common/useeft_wh/useeft_wh
!$omp threadprivate(/useeft_wh/) 

!      useeft_wh=.true. ! for comparison with vh@nnlo
      useeft_wh=.false. ! for best prediction
     
      mt2=mt**2
!      if(useeft_wh) then 
!      mrun2=mt**2
!      else
!      mrun2=massfrun(mt,scale,amz,nloop)**2 
!      endif

      A5LO=Alo5_VHt(i1,i2,i3,i4,i5,za,zb)
      A5NLOt=qqbgWH_topamp(i1,i2,i3,i4,i5,za,zb,mt2) 
!===== add in bottom loops too
      if(useeft_wh) return
      mt2=mb**2
!      mrun2=massfrun(mb,scale,amz,nloop)**2
      A5NLOt=A5NLOt+qqbgWH_topamp(i1,i2,i3,i4,i5,za,zb,mt2) 

      return 
      end

      function Alo5_VHt(i1,i2,i3,i4,i5,za,zb) 
      implicit none
      include 'types.f' 
      complex(dp):: Alo5_VHt
!======= this is the tree-level amplitude again, defined for a positive
!======= helicity gluon (i5), only reason this is used is to make factors
!+====== extracted for interference more convienent compared to those 
!======= defined in the _v.f routine 

!===== i1 is the negative helicity quark, i3 is the negative helicity lepton
      include 'constants.f' 
      include 'mxpart.f'
      include 'ewcouple.f' 
      include 'zprods_decl.f' 
      include 'sprods_com.f' 
      include 'masses.f' 
      include 'qcdcouple.f' 
      integer:: i1,i2,i3,i4,i5
      complex(dp):: prop125,prop34
      complex(dp):: fac,amp
      complex(dp)::zab2
      real(dp):: s34,s125
!=====begin statement function
      include 'cplx.h'
!=====end statement function

      zab2(i1,i2,i3,i4)=za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)

      s125=s(i1,i2)+s(i2,i5)+s(i1,i5)
      s34=s(i3,i4)

      prop125=cplx1(s125)/cplx2(s125-wmass**2,wmass*wwidth)
      prop34=cplx1(s34)/cplx2(s34-wmass**2,wmass*wwidth)

      
!====== overall pre-factor 
      fac=prop34*prop125*sqrt(gsq*gwsq)*gwsq*wmass*rt2*im  

!====== amplitude 
      amp=-((za(i1,i3)*zab2(i1,i2,i5,i4))
     &     /(s125*s34*za(i1,i5)*za(i2,i5)))

      Alo5_VHt=amp*fac
      return 
      end


      
