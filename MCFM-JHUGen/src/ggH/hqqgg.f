      subroutine hqqgg(p1,p2,p3,p4,ampsq,ampsq_ab,ampsq_ba,ampsq_sym)
      implicit none
      include 'types.f'

C     Taken from Kauffman,Desai,Risal
C     PRD 55 1997 (4009)
C     and checked with hep-ph/9903330
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_com.f'
      include 'sprods_com.f'
      integer:: p1,p2,p3,p4,j1,j2,j3
      complex(dp):: ab(2,2,2),ba(2,2,2),abppp,abppm,bappm
c      complex(dp):: abmmm,abmmp,bammp
      real(dp):: ampsq,ampsq_ab,ampsq_ba,ampsq_sym
C====statement functions
C  Eq 25
      abppp(p1,p2,p3,p4)=
     & +(za(p2,p1)*zb(p1,p3)+za(p2,p4)*zb(p4,p3))**2*zb(p1,p4)
     & /((s(p1,p2)+s(p1,p4)+s(p2,p4))*za(p2,p4))
     & *(1._dp/s(p1,p2)+1._dp/s(p1,p4))
     & -(za(p2,p1)*zb(p1,p4)+za(p2,p3)*zb(p3,p4))**2
     & *zb(p1,p3)/((s(p1,p2)+s(p1,p3)+s(p2,p3))*s(p1,p2)*za(p2,p3))
     & +(+za(p2,p3)*zb(p3,p1)+za(p2,p4)*zb(p4,p1))**2
     & /(zb(p1,p2)*za(p2,p3)*za(p2,p4)*za(p3,p4))

c      abmmm(p1,p2,p3,p4)=
c     & +(zb(p2,p1)*za(p1,p3)+zb(p2,p4)*za(p4,p3))**2
c     & *za(p1,p4)/((s(p1,p2)+s(p1,p4)+s(p2,p4))
c     & *zb(p2,p4))*(1._dp/s(p1,p2)+1._dp/s(p1,p4))
c     & -(+zb(p2,p1)*za(p1,p4)+zb(p2,p3)*za(p3,p4))**2
c     & *za(p1,p3)/((s(p1,p2)+s(p1,p3)+s(p2,p3))*s(p1,p2)*zb(p2,p3))
c     & +(+zb(p2,p3)*za(p3,p1)+zb(p2,p4)*za(p4,p1))**2
c     & /(za(p1,p2)*zb(p2,p3)*zb(p2,p4)*zb(p3,p4))
C--Eq 26
      abppm(p1,p2,p3,p4)=
     & -za(p2,p4)**3/(za(p1,p2)*za(p2,p3)*za(p3,p4))
     & +zb(p1,p3)**3/(zb(p1,p2)*zb(p1,p4)*zb(p3,p4))
c      abmmp(p1,p2,p3,p4)=
c     & -zb(p2,p4)**3/(zb(p1,p2)*zb(p2,p3)*zb(p3,p4))
c     & +za(p1,p3)**3/(za(p1,p2)*za(p1,p4)*za(p3,p4))

C--Eq 27
      bappm(p1,p2,p3,p4)=
     & -zb(p1,p3)**2*zb(p2,p3)/(zb(p1,p2)*zb(p2,p4)*zb(p3,p4))
     & +za(p1,p4)*za(p2,p4)**2/(za(p1,p2)*za(p1,p3)*za(p3,p4))

c      bammp(p1,p2,p3,p4)=
c     & -za(p1,p3)**2*za(p2,p3)/(za(p1,p2)*za(p2,p4)*za(p3,p4))
c     & +zb(p1,p4)*zb(p2,p4)**2/(zb(p1,p2)*zb(p1,p3)*zb(p3,p4))
C====end statement functions

C====It has been checked that taking the complex conjugate
C====only gets the answer different by an overall (and hence irrelevant)
C====phase (in some crossings).


CORIG      ab(2,2,2)=abppp(p1,p2,p3,p4)
CORIGc      ab(1,1,1)=abmmm(p1,p2,p3,p4)
CORIG      ab(1,1,1)=conjg(ab(2,2,2))
CORIG
CORIG
CORIG      ba(2,2,2)=abppp(p1,p2,p4,p3)
CORIGc      ba(1,1,1)=abmmm(p1,p2,p4,p3)
CORIG      ba(1,1,1)=conjg(ba(2,2,2))
CORIG
CORIG      ab(1,2,2)=abppp(p2,p1,p3,p4)
CORIGc      ab(2,1,1)=abmmm(p2,p1,p3,p4)
CORIG      ab(2,1,1)=conjg(ab(1,2,2))
CORIG
CORIG      ba(1,2,2)=abppp(p2,p1,p4,p3)
CORIGc      ba(2,1,1)=abmmm(p2,p1,p4,p3)
CORIG      ba(2,1,1)=conjg(ba(1,2,2))
CORIG
CORIG      ab(2,2,1)=abppm(p1,p2,p3,p4)
CORIGc      ab(1,1,2)=abmmp(p1,p2,p3,p4)
CORIG      ab(1,1,2)=conjg(ab(2,2,1))
CORIG
CORIG      ba(2,1,2)=abppm(p1,p2,p4,p3)
CORIGc      ba(1,2,1)=abmmp(p1,p2,p4,p3)
CORIG      ba(1,2,1)=conjg(ba(2,1,2))
CORIG
CORIG      ba(2,2,1)=bappm(p1,p2,p3,p4)
CORIGc      ba(1,1,2)=bammp(p1,p2,p3,p4)
CORIG      ba(1,1,2)=conjg(ba(2,2,1))
CORIG
CORIG      ab(2,1,2)=bappm(p1,p2,p4,p3)
CORIGc      ab(1,2,1)=bammp(p1,p2,p4,p3)
CORIG      ab(1,2,1)=conjg(ab(2,1,2))
CORIG
CORIG



      ba(2,2,2)=abppp(p1,p2,p3,p4)
      ba(1,1,1)=conjg(ba(2,2,2))

      ab(2,2,2)=abppp(p1,p2,p4,p3)
      ab(1,1,1)=conjg(ab(2,2,2))

      ab(1,2,2)=abppp(p2,p1,p3,p4)
      ab(2,1,1)=conjg(ab(1,2,2))

      ba(1,2,2)=abppp(p2,p1,p4,p3)
      ba(2,1,1)=conjg(ba(1,2,2))

      ba(2,2,1)=abppm(p1,p2,p3,p4)
      ba(1,1,2)=conjg(ba(2,2,1))

      ab(2,1,2)=abppm(p1,p2,p4,p3)
      ab(1,2,1)=conjg(ab(2,1,2))

      ab(2,2,1)=bappm(p1,p2,p3,p4)
      ab(1,1,2)=conjg(ab(2,2,1))

      ba(2,1,2)=bappm(p1,p2,p4,p3)
      ba(1,2,1)=conjg(ba(2,1,2))


c--- calculate the matrix element as the sum of two colour orderings,
c--- plus a colour-suppressed QE.e-_dplike piece which is symmetric
c--- in the ordering of the two gluons
      ampsq_ab=0._dp
      ampsq_ba=0._dp
      ampsq_sym=0._dp
      do j1=1,2
      do j2=1,2
      do j3=1,2
C         write(*,*) 'j1,j2,j3',j1,j2,j3
C         write(*,*) cf*xn**2/2._dp*abs(ab(j1,j2,j3))**2,
C     &      cf*xn**2/2._dp*abs(ba(j1,j2,j3))**2,
C     &        -cf/2._dp*abs(ab(j1,j2,j3)+ba(j1,j2,j3))**2
      ampsq_ab=ampsq_ab
     &  +cf*xn**2/2._dp*abs(ab(j1,j2,j3))**2
      ampsq_ba=ampsq_ba
     &  +cf*xn**2/2._dp*abs(ba(j1,j2,j3))**2
      ampsq_sym=ampsq_sym
     &  -cf/2._dp*abs(ab(j1,j2,j3)+ba(j1,j2,j3))**2
      enddo
      enddo
      enddo

      ampsq=ampsq_ab+ampsq_ba+ampsq_sym

      return
      end
