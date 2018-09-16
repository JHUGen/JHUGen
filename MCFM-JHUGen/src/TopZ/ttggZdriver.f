      subroutine ttggZdriver(p,etatb,etat,ampsq1,ampsq2,ampsq0)
      implicit none
      include 'types.f'
C     Calculate the matrix elements squared for the process
C      g(-p1)+g(-p2)-->l(p3)+a(p4)+t(p5)+t~(p6)
C     eta5 is the vector chosen to make the  t~ massleses
C     eta6 is the vector chosen to make the  t massleses

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zcouple.f'
      include 'ewcharge.f'
      include 'masses.f'
      real(dp):: p(mxpart,4),s34,ampsq1,ampsq2,ampsq0
      integer:: h1,h2,h3,h4,h56,etatb,etat
      complex(dp):: ampAB(2,2,2,2,2),ampBA(2,2,2,2,2),prop,
     & ampQED(2,2,2,2,2),
     & XABL(2,2,2,2,2),XABR(2,2,2,2,2),
     & XBAL(2,2,2,2,2),XBAR(2,2,2,2,2)
      call ttggZamp(p,6,5,1,2,3,4,etatb,etat,XABL,XABR)
      call ttggZamp(p,6,5,2,1,3,4,etatb,etat,XBAL,XBAR)

C     In the call s34 is the mass^2 of the lepton pair
      s34=2d0*(p(3,4)*p(4,4)-p(3,1)*p(4,1)-p(3,2)*p(4,2)-p(3,3)*p(4,3))
      prop=s34/cplx2(s34-zmass**2,zmass*zwidth)

      do h1=1,2
      do h2=1,2
      do h3=1,2
      do h4=1,2
      h56=1
      ampAB(h1,h2,h3,h4,h56)=
     &  (Q(2)*q1+L(2)*l1*prop)*XABL(h1,h2,h3,h4,h56)
     & +(Q(2)*q1+R(2)*l1*prop)*XABR(h1,h2,h3,h4,h56)
      ampBA(h1,h2,h3,h4,h56)=
     &  (Q(2)*q1+L(2)*l1*prop)*XBAL(h1,h2,h4,h3,h56)
     & +(Q(2)*q1+R(2)*l1*prop)*XBAR(h1,h2,h4,h3,h56)
      h56=2
      ampAB(h1,h2,h3,h4,h56)=
     &  (Q(2)*q1+L(2)*r1*prop)*XABL(h1,h2,h3,h4,h56)
     & +(Q(2)*q1+R(2)*r1*prop)*XABR(h1,h2,h3,h4,h56)
      ampBA(h1,h2,h3,h4,h56)=
     &  (Q(2)*q1+L(2)*r1*prop)*XBAL(h1,h2,h4,h3,h56)
     & +(Q(2)*q1+R(2)*r1*prop)*XBAR(h1,h2,h4,h3,h56)
      enddo
      enddo
      enddo
      enddo

      ampsq0=0d0
      ampsq1=0d0
      ampsq2=0d0
      do h1=1,2
      do h2=1,2
      do h3=1,2
      do h4=1,2
      do h56=1,2
      ampQED(h1,h2,h3,h4,h56)=
     & ampAB(h1,h2,h3,h4,h56)+ampBA(h1,h2,h3,h4,h56)
      ampsq1=ampsq1
     &  +real(ampAB(h1,h2,h3,h4,h56)*conjg(ampAB(h1,h2,h3,h4,h56)))
      ampsq2=ampsq2
     &  +real(ampBA(h1,h2,h3,h4,h56)*conjg(ampBA(h1,h2,h3,h4,h56)))
      ampsq0=ampsq0
     & -real(ampQED(h1,h2,h3,h4,h56)*conjg(ampQED(h1,h2,h3,h4,h56)))
     & /xn**2
      enddo
      enddo
      enddo
      enddo
      enddo

      return
      end
