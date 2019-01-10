      subroutine gg_zz_int_freenorm(p,hcoupl,msq)
      implicit none
c--- Author: F
c--- For now, work in the approximation of BDK, i.e. that we
c--- retain only 1/mt^2 for top quark loops
c--- Box contributions are then complete (terms of order 1/mt^4 discarded)
c--- No triangle contributions

      include 'constants.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'zcouple.f'
      include 'masses.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      integer h1,h2,h34,h56,up,dn
      double precision p(mxpart,4),msq(fn:nf,fn:nf),msqgg,fac
      double precision cvec(2),cax(2),
     & pttwo,ptWsafetycut_massless,cl1(2),cl2(2)
      double complex Avec(2,2,2,2),prop34,prop56
      double complex a64v
      parameter(up=2,dn=1)
      double precision mfsq,tau,tauinv,rt
      double complex fachiggs,amphiggs,f
      double complex Atot(2,2,2,2),Ahiggs(2,2,2,2),Acont(2,2,2,2)
      double complex e3De4(2,2)
      double complex hcoupl

c--- omit massless loops for pt(W) < "ptWsafetycut_massless" (for num. stability)
c      ptWsafetycut_massless=0.05d0
c--- omit massless loops for pt(W) < "ptWsafetycut_massless" (for num. stability)
      ptWsafetycut_massless=2d0

c--- set up spinor products      
      call spinoru(6,p,za,zb)

c--------------------------------------------------------------------------------------
c--------------------------------------------------------------------------------------
c--- CONTINUUM AMPLITUDE BELOW
c--------------------------------------------------------------------------------------
     
c--- fill amplitudes
c--- labels are: helicity of gluons, lepton 3 and lepton 5

      if (pttwo(3,4,p) .lt. ptWsafetycut_massless) then
c--- ensure numerical stability: set massless loops to zero
c--- for pt(W) < "ptWsafetycut_massless" GeV
        do h1=1,2
        do h2=1,2
        do h34=1,2
        do h56=1,2
          Avec(h1,h2,h34,h56)=czip
        enddo
        enddo
        enddo
        enddo
      else
      Avec(2,2,1,1)=a64v('q+qb-g-g-',3,4,1,2,6,5,zb,za)*(-im)
      Avec(2,1,1,1)=a64v('q+qb-g-g+',3,4,1,2,6,5,zb,za)*(-im)
      Avec(1,2,1,1)=a64v('q+qb-g+g-',3,4,1,2,6,5,zb,za)*(-im)
      Avec(1,1,1,1)=a64v('q+qb-g+g+',3,4,1,2,6,5,zb,za)*(-im)

      Avec(2,2,2,1)=a64v('q+qb-g+g+',3,4,1,2,5,6,za,zb)*(-im)
      Avec(2,1,2,1)=a64v('q+qb-g+g-',3,4,1,2,5,6,za,zb)*(-im)
      Avec(1,2,2,1)=a64v('q+qb-g-g+',3,4,1,2,5,6,za,zb)*(-im)
      Avec(1,1,2,1)=a64v('q+qb-g-g-',3,4,1,2,5,6,za,zb)*(-im)

      Avec(2,2,1,2)=a64v('q+qb-g-g-',3,4,1,2,5,6,zb,za)*(-im)
      Avec(2,1,1,2)=a64v('q+qb-g-g+',3,4,1,2,5,6,zb,za)*(-im)
      Avec(1,2,1,2)=a64v('q+qb-g+g-',3,4,1,2,5,6,zb,za)*(-im)
      Avec(1,1,1,2)=a64v('q+qb-g+g+',3,4,1,2,5,6,zb,za)*(-im)

      Avec(2,2,2,2)=a64v('q+qb-g+g+',3,4,1,2,6,5,za,zb)*(-im)
      Avec(2,1,2,2)=a64v('q+qb-g+g-',3,4,1,2,6,5,za,zb)*(-im)
      Avec(1,2,2,2)=a64v('q+qb-g-g+',3,4,1,2,6,5,za,zb)*(-im)
      Avec(1,1,2,2)=a64v('q+qb-g-g-',3,4,1,2,6,5,za,zb)*(-im)
      endif


c--- propagator factors
      prop34=s(3,4)/dcmplx(s(3,4)-zmass**2,zmass*zwidth)
      prop56=s(5,6)/dcmplx(s(5,6)-zmass**2,zmass*zwidth)

c--- left and right-handed lepton couplings as an array
c--- cl1 associated with Z(3+4), cl2 associated with Z(5+6)
      cl1(1)=l1 
      cl1(2)=r1
      cl2(1)=l2
      cl2(2)=r2

c--- vector and axial couplings as an array for up/down quarks
      cvec(up)=half*(L(up)+R(up))
      cvec(dn)=half*(L(dn)+R(dn))
      cax(up)=half*(L(up)-R(up))
      cax(dn)=half*(L(dn)-R(dn))
      
c--- assume 5 massless flavors in the loop, (2 up, 3 down)
      do h1=1,2
      do h2=1,2
      do h34=1,2
      do h56=1,2
         Acont(h1,h2,h34,h56) = Avec(h1,h2,h34,h56)
c--- vector couplings in the loop
     & *(2d0*(Qu*q1+cvec(up)*cl1(h34)*prop34)
     &      *(Qu*q1+cvec(up)*cl2(h56)*prop56)
     &  +3d0*(Qd*q1+cvec(dn)*cl1(h34)*prop34)
     &      *(Qd*q1+cvec(dn)*cl2(h56)*prop56)
c--- axial couplings in the loop
     &  +2d0*(cax(up)*cl1(h34)*prop34)*(cax(up)*cl2(h56)*prop56)
     &  +3d0*(cax(dn)*cl1(h34)*prop34)*(cax(dn)*cl2(h56)*prop56)
     & ) * (2d0*4d0*esq*gsq/(16d0*pisq)*esq)
      enddo
      enddo
      enddo
      enddo

c--------------------------------------------------------------------------------------
c--------------------------------------------------------------------------------------
c--- HIGGS AMPLITUDES BELOW
c--------------------------------------------------------------------------------------

c---  Z decay amplitud
c---  labels are: helicity of lepton 3 and lepton 5
      e3De4(1,1) = l1*l2*2d0*za(3,5)*zb(6,4)/(s(3,4)*s(5,6))
      e3De4(2,2) = r1*r2*2d0*zb(3,5)*za(6,4)/(s(3,4)*s(5,6))
      e3De4(1,2) = l1*r2*2d0*za(3,6)*zb(5,4)/(s(3,4)*s(5,6))
      e3De4(2,1) = r1*l2*2d0*zb(3,6)*za(5,4)/(s(3,4)*s(5,6))

      
c--- fill amplitudes with contributions of Higgs: top loop
      mfsq=mt**2
      tau=s(1,2)/(4d0*mfsq)
      tauinv=1d0/tau
      fachiggs=cone/dcmplx(s(1,2)-hmass**2,hmass*hwidth)

      fachiggs=fachiggs*(2d0*gwsq*gsq/(16d0*pisq)*gwsq/2d0)
      fachiggs=fachiggs*(2d0*zmass*xw/dsqrt(one-xw)/wmass)
      fachiggs=fachiggs*prop34*prop56

!====== SEYMOUR ISA APPROX 
!      num_c=dcmplx(hmass**2/s(1,2)) 
!      fachiggs=num_c/dcmplx(s(1,2)-hmass**2,s(1,2)*hwidth/hmass)

      if (tau .le. 1d0) then
         f=dcmplx(dasin(sqrt(tau))**2)
      elseif (tau .gt. 1d0) then
         rt=sqrt(1d0-tauinv)
         f=-0.25d0*(dcmplx(log((1d0+rt)/(1d0-rt)))-im*pi)**2
      else
         f=czip
      endif
      amphiggs=mfsq*(cone+(cone-dcmplx(tauinv))*f)*im

c--- CROSS-CHECK: point-like limit
c---      amphiggs=im * s(1,2)/6d0
c---

c--   labels are: helicity of gluons, lepton 3 and lepton 5
      do h34=1,2
      do h56=1,2
        Ahiggs(1,1,h34,h56)=fachiggs*
     .        amphiggs*za(1,2)/zb(2,1)*e3De4(h34,h56)
        Ahiggs(1,2,h34,h56)=czip
        Ahiggs(2,1,h34,h56)=czip
        Ahiggs(2,2,h34,h56)=fachiggs*
     .       amphiggs*zb(1,2)/za(2,1)*e3De4(h34,h56)
      enddo
      enddo

c--- fill amplitudes with contributions of Higgs: bottom loop
      mfsq=mb**2
      tau=s(1,2)/(4d0*mfsq)
      tauinv=1d0/tau

      if (tau .le. 1d0) then
         f=dcmplx(dasin(sqrt(tau))**2)
      elseif (tau .gt. 1d0) then
         rt=sqrt(1d0-tauinv)
         f=-0.25d0*(dcmplx(log((1d0+rt)/(1d0-rt)))-im*pi)**2
      else
         f=czip
      endif
      amphiggs=mfsq*(cone+(cone-dcmplx(tauinv))*f)*im

      do h34=1,2
      do h56=1,2
      Ahiggs(1,1,h34,h56)=Ahiggs(1,1,h34,h56)+
     . fachiggs*amphiggs*za(1,2)/zb(2,1)*e3De4(h34,h56)
      Ahiggs(2,2,h34,h56)=Ahiggs(2,2,h34,h56)+
     . fachiggs*amphiggs*zb(1,2)/za(2,1)*e3De4(h34,h56)
      enddo
      enddo

c--------------------------------------------------------------------------------------
c--------------------------------------------------------------------------------------
c--- TOTAL AMPLITUDE BELOW
c--------------------------------------------------------------------------------------
      
      do h1=1,2
      do h2=1,2
      do h34=1,2
      do h56=1,2
         Atot(h1,h2,h34,h56)=Acont(h1,h2,h34,h56)
     .        +hcoupl*Ahiggs(h1,h2,h34,h56)
      enddo
      enddo
      enddo
      enddo
      
      msqgg=0d0
      do h1=1,2
      do h2=1,2
      do h34=1,2
      do h56=1,2
      msqgg=msqgg+cdabs(Atot(h1,h2,h34,h56))**2
c Uncomment to remove background
c     . -cdabs(Acont(h1,h2,h34,h56))**2
      enddo
      enddo
      enddo
      enddo

      
c--- overall factor from diagrams
      fac=avegg*V

      msqgg=msqgg*fac

      msq(0,0)=msqgg
      
      return
      end
      
