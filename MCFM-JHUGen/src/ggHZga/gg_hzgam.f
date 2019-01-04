      subroutine gg_hzgam(p,msq)
      implicit none
C----Author: R.K.Ellis July 2012
c----Lowest order matrix element for H production
c----in the heavy quark (mt=Infinity) limit.
C----averaged over initial colours and spins
c    g(-p1)+g(-p2)-->H --> Z/gamma^*--(l(p3)+a(p4)) + gamma(p5) 
c---
      include 'constants.f'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'sprods_com.f'
      double precision msq(-nf:nf,-nf:nf),p(mxpart,4)
      double precision decay,gg,Asq,HZgamMSQ

c---set msq=0 to initialize
      msq(:,:)=0d0

      call dotem(5,p,s)
      decay=HZgamMSQ(3,4,5)
     & /((s(1,2)-hmass**2)**2+(hmass*hwidth)**2)
      
      Asq=(as/(3d0*pi))**2/vevsq
      gg=0.5d0*V*Asq*s(1,2)**2

c---calculate propagators
      msq(0,0)=avegg*gg*decay

      return
      end


