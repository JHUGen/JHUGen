      double complex function C0DDHK(s12,s34,msq)
      implicit none
C-----Author: R.K.Ellis July 2012
C-----The scalar triangle function taken from
C-----Djouadi,Driesen,Hollik,Kraft arXiv:hep-ph/9701342v1, Eq(5)
      double precision s12,s34,msq,tauZ,tauH
      double complex FFF
c      double complex qlC0DDHK
      tauZ=4d0*msq/s34
      tauH=4d0*msq/s12
      C0DDHK=-2d0*tauZ*tauH/(tauZ-tauH)*(FFF(tauZ)-FFF(tauH))
C-----changed sign in order to agree with normal QCDLoop definition of scalar integral
      C0DDHK=-C0DDHK/(4d0*msq)

c      write(6,*) C0DDHK
c      C0DDHK=qlC0DDHK(s12,s34,msq)
c      write(6,*) C0DDHK
c      pause
      return
      end


      double complex function C2DDHK(s12,s34,msq)
      implicit none
C-----Author: R.K.Ellis July 2012
C-----The combination of triangle functions taken from
C-----Djouadi,Driesen,Hollik,Kraft arXiv:hep-ph/9701342v1, Eq(5)
      double precision s12,s34,msq,tauZ,tauH
      double complex FFF,GGG
c      double complex qlC2DDHK
      tauZ=4d0*msq/s34
      tauH=4d0*msq/s12
      C2DDHK=dcmplx(tauZ*tauH/(2d0*(tauZ-tauH)))
     & +tauZ*tauH**2/(2d0*(tauZ-tauH)**2)
     & *(tauZ*(FFF(tauZ)-FFF(tauH))+2d0*(GGG(tauZ)-GGG(tauH)))
C-----changed sign in order to agree with QCDLoop result
      C2DDHK=-C2DDHK/(4d0*msq)

c      write(6,*) C2DDHK
c      C2DDHK=qlC2DDHK(s12,s34,msq)
c      write(6,*) C2DDHK
c      pause

      return
      end

      double complex function FFF(tau)
      implicit none
C-----Author: R.K.Ellis July 2012
C-----The basic triangle function taken from
C-----Djouadi,Driesen,Hollik,Kraft arXiv:hep-ph/9701342v1, Eq(6)
      double precision xlog,tau,rttauinv,rtomtau,pi
      pi=2d0*asin(1d0)
      if (tau .ge. 1d0) then
      rttauinv=1d0/sqrt(tau)
      FFF=dcmplx(asin(rttauinv)**2)
      else
      rtomtau=sqrt(1d0-tau)
      xlog=log((1d0+rtomtau)/(1d0-rtomtau))
      FFF=-0.25d0*dcmplx(xlog,-pi)**2
      endif
      return
      end

      double complex function GGG(tau)
      implicit none
C-----Author: R.K.Ellis July 2012
C-----The other needed function taken from
C-----Djouadi,Driesen,Hollik,Kraft arXiv:hep-ph/9701342v1, Eq(7)
      double precision tau,rttauinv,rtomtau,pi,xlog
      pi=2d0*asin(1d0)
      if (tau .ge. 1d0) then
      rtomtau=sqrt(tau-1d0)
      rttauinv=1d0/sqrt(tau)
      GGG=dcmplx(rtomtau*asin(rttauinv))
      else
      rtomtau=sqrt(1d0-tau)
      xlog=log((1d0+rtomtau)/(1d0-rtomtau))
      GGG=rtomtau/2d0*dcmplx(xlog,-pi)
      endif
      return
      end


