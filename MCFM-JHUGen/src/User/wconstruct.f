      subroutine wconstruct(etvec,plept,lepton,etmissing)
      implicit none
      include 'constants.f'
      include 'masses.f'
      double precision etvec(2),plept(4),pw(4),etmissing(4),
     . pte,ptn,ptwsq,a,b,discr,plnp,plnm,Ee,ple,mtsq,Enp,Enm
      double precision yrap1,yrap2
      integer j
      character*2 lepton
      
      do j=1,2
        pw(j)=plept(j)+etvec(j)
        etmissing(j)=etvec(j)
      enddo
 
      pte=dsqrt(plept(1)**2+plept(2)**2)
      ptn=dsqrt(etvec(1)**2+etvec(2)**2)
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
      if (wmass**2 .gt. mtsq) then
        discr=Ee*dsqrt((wmass**2-mtsq)*(wmass**2-mtsq+4d0*pte*ptn))
        a=pte**2
        b=-(2d0*pte*ptn+wmass**2-mtsq)*ple
c      c=-0.25d0*(wmass**2-mtsq+2*pte*ptn)**2
c     # +ptn**2*Ee**2
c      write(6,*) 'discr',sqrt(b**2-4*a*c),discr
      else
        discr=0d0
        a=pte**2
        b=-2*pte*ptn*ple
      endif

c      wsq=pw(4)**2-pw(1)**2-pw(2)**2-pw(3)**2
c      write(6,*) 'cal:mw',wsq,plntrue
      plnp=(-b+discr)/(two*a)
      plnm=(-b-discr)/(two*a)

      Enm=dsqrt(plnm**2+ptn**2)
      Enp=dsqrt(plnp**2+ptn**2)
      yrap1=(Enm+plnm)/(Enm-plnm)
      yrap2=(Enp+plnp)/(Enp-plnp)
      if (yrap1 .lt. 1d-13) then
        yrap1=100d0
      else
        yrap1=0.5d0*dlog(yrap1)
      endif
      if (yrap2 .lt. 1d-13) then
        yrap2=100d0
      else
        yrap2=0.5d0*dlog(yrap2)
      endif
      
      if     ((lepton.eq.'ea') .or. (lepton.eq.'ma')) then
c--- choose the solution with the larger rapidity
        if (yrap1 .gt. yrap2) then
          etmissing(4)=Enm    
          etmissing(3)=plnm   
        else
          etmissing(4)=Enp    
          etmissing(3)=plnp   
        endif
      elseif ((lepton.eq.'el') .or. (lepton.eq.'ml')) then
c--- choose the solution with the smaller rapidity
        if (yrap1 .lt. yrap2) then
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
      

