      subroutine wconstruct(etvec,plept,lepton,etmissing)
      implicit none
      include 'types.f'
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      real(dp):: etvec(2),plept(4),pw(4),etmissing(4),
     & pte,ptn,ptwsq,a,b,discr,plnp,plnm,Ee,ple,mtsq,Enp,Enm
      real(dp):: yrap1,yrap2
      integer:: j
      character*2 lepton
      
      do j=1,2
        pw(j)=plept(j)+etvec(j)
        etmissing(j)=etvec(j)
      enddo
 
      pte=sqrt(plept(1)**2+plept(2)**2)
      ptn=sqrt(etvec(1)**2+etvec(2)**2)
      ptwsq=(pw(1)**2+pw(2)**2)
      Ee=plept(4)
      ple=plept(3)
      mtsq=(pte+ptn)**2-ptwsq

c  + (Mw^2 - MT^2) * Ee^2 * ( 4*pte*ptn + Mw^2 - Mt^2 )
c 
c - pln
c  * ple*( 2*pte*ptn + Mw^2 - MT^2 )
c 
c + pln^2
c  * ( pte^2 ) + 0.

c      write(6,*) 'wmass**2,mtsq',wmass**2,mtsq
      if (wmass**2 > mtsq) then
        discr=Ee*sqrt((wmass**2-mtsq)*(wmass**2-mtsq+4._dp*pte*ptn))
        a=pte**2
        b=-(2._dp*pte*ptn+wmass**2-mtsq)*ple
c      c=-0.25_dp*(wmass**2-mtsq+2*pte*ptn)**2
c     # +ptn**2*Ee**2
c      write(6,*) 'discr',sqrt(b**2-4*a*c),discr
      else
        discr=0._dp
        a=pte**2
        b=-2*pte*ptn*ple
      endif

c      wsq=pw(4)**2-pw(1)**2-pw(2)**2-pw(3)**2
c      write(6,*) 'cal:mw',wsq,plntrue
      plnp=(-b+discr)/(two*a)
      plnm=(-b-discr)/(two*a)

      Enm=sqrt(plnm**2+ptn**2)
      Enp=sqrt(plnp**2+ptn**2)
      yrap1=(Enm+plnm)/(Enm-plnm)
      yrap2=(Enp+plnp)/(Enp-plnp)
      if (yrap1 < 1.e-13_dp) then
        yrap1=100._dp
      else
        yrap1=0.5_dp*log(yrap1)
      endif
      if (yrap2 < 1.e-13_dp) then
        yrap2=100._dp
      else
        yrap2=0.5_dp*log(yrap2)
      endif
      
      if     ((lepton=='ea') .or. (lepton=='ma')) then
c--- choose the solution with the larger rapidity
        if (yrap1 > yrap2) then
          etmissing(4)=Enm    
          etmissing(3)=plnm   
        else
          etmissing(4)=Enp    
          etmissing(3)=plnp   
        endif
      elseif ((lepton=='el') .or. (lepton=='ml')) then
c--- choose the solution with the smaller rapidity
        if (yrap1 < yrap2) then
          etmissing(4)=Enm    
          etmissing(3)=plnm   
        else
          etmissing(4)=Enp    
          etmissing(3)=plnp   
        endif
      else
        write(6,*) 'Problem: not a lepton in wconstruct!'
        stop
      endif
                    
      return
      end
      

