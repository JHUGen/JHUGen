      function C0DDHK(s12,s34,msq)
      implicit none
      include 'types.f'
      include 'cplx.h'
      complex(dp):: C0DDHK

C-----Author: R.K.Ellis July 2012
C-----The scalar triangle function taken from
C-----Djouadi,Driesen,Hollik,Kraft arXiv:hep-ph/9701342v1, Eq(5)
      real(dp):: s12,s34,msq,tauZ,tauH
      complex(dp):: FFF
c      complex(dp):: qlC0DDHK
      tauZ=4._dp*msq/s34
      tauH=4._dp*msq/s12
      C0DDHK=-2._dp*tauZ*tauH/(tauZ-tauH)*(FFF(tauZ)-FFF(tauH))
C-----changed sign in order to agree with normal QCDLoop definition of scalar integral
      C0DDHK=-C0DDHK/(4._dp*msq)

c      write(6,*) C0DDHK
c      C0DDHK=qlC0DDHK(s12,s34,msq)
c      write(6,*) C0DDHK
c      pause
      return
      end


      function C2DDHK(s12,s34,msq)
      implicit none
      include 'types.f'
      include 'cplx.h'
      complex(dp):: C2DDHK

C-----Author: R.K.Ellis July 2012
C-----The combination of triangle functions taken from
C-----Djouadi,Driesen,Hollik,Kraft arXiv:hep-ph/9701342v1, Eq(5)
      real(dp):: s12,s34,msq,tauZ,tauH
      complex(dp):: FFF,GGG
c      complex(dp):: qlC2DDHK
      tauZ=4._dp*msq/s34
      tauH=4._dp*msq/s12
      C2DDHK=cplx1(tauZ*tauH/(2._dp*(tauZ-tauH)))
     & +tauZ*tauH**2/(2._dp*(tauZ-tauH)**2)
     & *(tauZ*(FFF(tauZ)-FFF(tauH))+2._dp*(GGG(tauZ)-GGG(tauH)))
C-----changed sign in order to agree with QCDLoop result
      C2DDHK=-C2DDHK/(4._dp*msq)

c      write(6,*) C2DDHK
c      C2DDHK=qlC2DDHK(s12,s34,msq)
c      write(6,*) C2DDHK
c      pause

      return
      end

      function FFF(tau)
      implicit none
      include 'types.f'
      include 'cplx.h'
      complex(dp):: FFF

C-----Author: R.K.Ellis July 2012
C-----The basic triangle function taken from
C-----Djouadi,Driesen,Hollik,Kraft arXiv:hep-ph/9701342v1, Eq(6)
      real(dp):: xlog,tau,rttauinv,rtomtau,pi
      pi=2._dp*asin(1._dp)
      if (tau >= 1._dp) then
      rttauinv=1._dp/sqrt(tau)
      FFF=cplx1(asin(rttauinv)**2)
      else
      rtomtau=sqrt(1._dp-tau)
      xlog=log((1._dp+rtomtau)/(1._dp-rtomtau))
      FFF=-0.25_dp*cplx2(xlog,-pi)**2
      endif
      return
      end

      function GGG(tau)
      implicit none
      include 'types.f'
      include 'cplx.h'
      complex(dp):: GGG

C-----Author: R.K.Ellis July 2012
C-----The other needed function taken from
C-----Djouadi,Driesen,Hollik,Kraft arXiv:hep-ph/9701342v1, Eq(7)
      real(dp):: tau,rttauinv,rtomtau,pi,xlog
      pi=2._dp*asin(1._dp)
      if (tau >= 1._dp) then
      rtomtau=sqrt(tau-1._dp)
      rttauinv=1._dp/sqrt(tau)
      GGG=cplx1(rtomtau*asin(rttauinv))
      else
      rtomtau=sqrt(1._dp-tau)
      xlog=log((1._dp+rtomtau)/(1._dp-rtomtau))
      GGG=rtomtau/2._dp*cplx2(xlog,-pi)
      endif
      return
      end


