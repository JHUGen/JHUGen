      subroutine hgamgamdecay(p,j,k,hdecay)
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'msbarmasses.f'
      include 'ewcouple.f'
      include 'ewcharge.f'
      include 'couple.f'
      include 'part.f'
      include 'first.f'
      integer j,k
      double complex Iw,Iq,Ftriangle
      double precision p(mxpart,4),prefac,mh,hdecay
      double precision x_t,x_b,x_w,x,mt_eff,mb_eff,massfrun
      save mt_eff,mb_eff
!$threadprivate(mt_eff,mb_eff)
C---statement functions
      Iq(x)=dcmplx(4d0*x)*(ctwo+dcmplx(4d0*x-1d0)*Ftriangle(x))
      Iw(x)=-ctwo*(dcmplx(6d0*x+1d0)
     & +dcmplx(6d0*x*(2d0*x-1d0))*Ftriangle(x))
C---end statement functions
      mh=sqrt((p(j,4)+p(k,4))**2
     &       -(p(j,1)+p(k,1))**2
     &       -(p(j,2)+p(k,2))**2
     &       -(p(j,3)+p(k,3))**2)

      if (first) then
c--- run mt to appropriate scale
        if (part .eq. 'lord') then
          mb_eff=massfrun(mb_msbar,hmass,amz,1)
          mt_eff=massfrun(mt_msbar,hmass,amz,1)
        else
          mb_eff=massfrun(mb_msbar,hmass,amz,2)
          mt_eff=massfrun(mt_msbar,hmass,amz,2)
        endif
        first=.false.
      endif

C---Total width multiplied by a factor of 16*pi*mh 
C---to get matrix element squared.
c---maybe it would be better to add esq at a higher scale.
      prefac=(esq/(4d0*pi))**2*Gf*mh**4/(8d0*rt2*pi**2)
      x_b=(mb_eff/mh)**2
      x_t=(mt_eff/mh)**2
      x_w=(wmass/mh)**2
      hdecay=prefac
     & *abs(xn*(Q(1)**2*Iq(x_b)+Q(2)**2*Iq(x_t))+Iw(x_w))**2
      return
      end
