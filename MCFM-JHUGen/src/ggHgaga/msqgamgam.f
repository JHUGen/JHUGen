      function msqgamgam(mh)
      implicit none
      include 'types.f'
      real(dp):: msqgamgam
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'msbarmasses.f'
      include 'ewcouple.f'
      include 'ewcharge.f'
      include 'couple.f'
      include 'kpart.f'
      include 'first.f'
      complex(dp):: Iw,Iq,Ftriangle
      real(dp):: prefac,mh
      real(dp):: x_t,x_b,x_w,x,mt_eff,mb_eff,massfrun
      save mt_eff,mb_eff
!$omp threadprivate(mt_eff,mb_eff)

C---statement functions
      Iq(x)=four*x*(ctwo+(four*x-one)*Ftriangle(x))
      Iw(x)=-ctwo*(cplx1(six*x+one)
     & +(six*x*(two*x-one))*Ftriangle(x))
C---end statement functions


      if (first) then
c--- run mt to appropriate scale
        if (kpart==klord) then
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
      prefac=(esq/(four*pi))**2*Gf*mh**4/(eight*rt2*pi**2)
      x_b=(mb_eff/mh)**2
      x_t=(mt_eff/mh)**2
      x_w=(wmass/mh)**2
      msqgamgam=prefac
     & *abs(xn*(Q(1)**2*Iq(x_b)+Q(2)**2*Iq(x_t))+Iw(x_w))**2

      return
      end

