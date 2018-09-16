      subroutine gg_ww_int(p,msq)
      implicit none
      include 'types.f'

c--- Author: R. K. Ellis, May 2011
c--- For now, work in the approximation of two massless isodoublets
c--- Box contributions are then complete
c--- Triangle (vector) pieces always vanish
c--- Triangle (axial) pieces cancel for massless isodoublets

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'masses.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      include 'Higgsint.f'
      include 'scale.f'
      include 'noglue.f'
      include 'anom_higgs.f'
      include 'kprocess.f'
      include 'docheck.f'
      include 'first.f'
      integer:: h1,h2,nu,i,j,k,om,del1,del2,k12h,k34h,k56h11,k34h11,e
      real(dp):: p(mxpart,4),msq(fn:nf,fn:nf),msqgg,fac
      real(dp):: mfsq,tau,tauinv,rt,pttwo,rescale
      complex(dp):: Avec(2,2),Ahiggs(2,2),Agen3(2,2),Atot(2,2),
     & faccont,fachiggs,amphiggs,f,e3De4,sum(2,2,-2:0)

      real(dp):: phi,muk,rho,ssig,csig,theta,
     & p1true(4),p2true(4),p3true(4),p4true(4),p5true(4),p6true(4),
     & dot,s12,s34,s56,dot1256,afac,bfac,gden,delta,
     & dot1234,dot3456,pmax,ptWsafetycut_massive,ptWsafetycut_massless
      complex(dp):: fvs,fvf,box(2,2,-2:0),triang(2,2,-2:0),
     & bub(2,2,-2:0)
      character*9 pp,pm
      logical:: includegens1and2,includegen3,
     & caseggWW4l,caseHWWHpI,caseHWWint
      parameter(pp='q+qb-g+g+',pm='q+qb-g+g-')
      parameter(del1=7,del2=8)
      parameter(k12h=9,k34h=10,k56h11=11,k34h11=12)
      save caseggWW4l,caseHWWHpI,caseHWWint
!$omp threadprivate(caseggWW4l,caseHWWHpI,caseHWWint)

c--- set this to true to include generations 1 and 2 of (light) quarks
      includegens1and2=.true.
c--- set this to true to include 3rd generation of quarks (t,b)
      includegen3=.true.

c--- omit massive loops for pt(W) < "ptWsafetycut_massive" (for num. stability)
      ptWsafetycut_massive=2._dp

c--- omit massless loops for pt(W) < "ptWsafetycut_massless" (for num. stability)
      ptWsafetycut_massless=0.05_dp

      if (first) then
        write(6,*)'****************************************************'
        write(6,*)'*                                                  *'
        if (includegens1and2) then
        write(6,*)'*  gg->WW box loop includes gens. 1 and 2          *'
        else
        write(6,*)'*  gg->WW box loop does not include gens. 1 and 2  *'
        endif
        if (includegen3) then
        write(6,*)'*  gg->WW box loop includes 3rd generation         *'
        else
        write(6,*)'*  gg->WW box loop does not include 3rd gen.       *'
        endif
        write(6,*)'*                                                  *'
        write(6,*)'****************************************************'
        first=.false.
        if (includegen3) call qlinit
        caseggWW4l=.false.
        caseHWWHpI=.false.
        caseHWWint=.false.
        if ((kcase==kggWW4l)
     & .or. (kcase==kWWqqbr)) caseggWW4l=.true.
        if (kcase==kHWWHpi) caseHWWHpI=.true.
        if (kcase==kHWWint) caseHWWint=.true.
      endif

c--- if neither contribution is included print warning message and stop
      if ((includegens1and2 .eqv. .false.) .and.
     &    (includegen3      .eqv. .false.)) then
         write(6,*) 'Box loop is set to zero, please edit gg_WW_int.f'
         stop
      endif
c--- if noglue print warning message and stop
      if (noglue) then
         write(6,*) 'Please set noglue .false. in input file'
         stop
      endif

c--- logical:: variable "docheck"
c---   .true .  --> print out coefficients of integrals at special point
c---   .false.  --> run as normal
      docheck=.false.

c--- set flag to signal calculation of Higgs related contributions only,
c---  e.g. entering interference
      if     (caseggWW4l) then
        Higgsint=.false.
      elseif ((caseHWWHpi) .or. (caseHWWint)) then
        Higgsint=.true.
      else
        write(6,*) 'Unexpected case in gg_WW_int: ',kcase
        stop
      endif

      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0._dp
      enddo
      enddo

c--- rescale
      pmax=abs(p(1,4))
      do i=2,6
      if (abs(p(i,4)) > pmax) pmax=abs(p(i,4))
      enddo
      pmax=1._dp
      do i=1,6
      do j=1,4
      p(i,j)=p(i,j)/pmax
      enddo
      enddo
      mt=mt/pmax
      musq=musq/pmax**2

c--- These lines set up the point for the numerical check
       if (docheck) then
         Higgsint=.false.
         includegen3=.true.
         include 'kinpoint.f'
         mt=0.5_dp
         hmass=5.6_dp
         hwidth=3.21_dp
         zmass=0.53342_dp
         zwidth=0.3414523523_dp
         do nu=1,4
         om=nu-1
         if (nu==1) om=4
         p(1,om)=p1true(nu)
         p(2,om)=p2true(nu)
         p(3,om)=p3true(nu)
         p(4,om)=p4true(nu)
         p(5,om)=p5true(nu)
         p(6,om)=p6true(nu)
         enddo
      endif

c--- if computing 3rd generation, set up extra flattened vectors
      if (includegen3) then
C--- define flattened vectors (k12 and k56)
      s12=2._dp*dot(p,1,2)
      s56=2._dp*dot(p,5,6)
      s34=2._dp*dot(p,3,4)
      dot1256=0.5_dp*(s34-s12-s56)
      delta=dot1256**2-s12*s56
      gden=dot1256+sqrt(delta)
      afac=s12/gden
      bfac=s56/gden

      do nu=1,4
      p(del2,nu)=one/(one-afac*bfac)
     & *(p(1,nu)+p(2,nu)-afac*(p(5,nu)+p(6,nu)))
      p(del1,nu)=one/(one-afac*bfac)
     & *(p(5,nu)+p(6,nu)-bfac*(p(1,nu)+p(2,nu)))
      enddo

C--- define flattened vectors (k12 and k34)
      dot1234=0.5_dp*(s56-s12-s34)
      delta=dot1234**2-s12*s34
      gden=dot1234+sqrt(delta)
      afac=s12/gden
      bfac=s34/gden

      do nu=1,4
      p(k12h,nu)=one/(one-afac*bfac)
     & *(p(1,nu)+p(2,nu)-afac*(p(3,nu)+p(4,nu)))
      p(k34h,nu)=one/(one-afac*bfac)
     & *(p(3,nu)+p(4,nu)-bfac*(p(1,nu)+p(2,nu)))
      enddo

C--- define flattened vectors (k56 and k34)
      dot3456=0.5_dp*(s12-s56-s34)
      delta=dot3456**2-s56*s34
      gden=dot3456+sqrt(delta)
      afac=s56/gden
      bfac=s34/gden

      do nu=1,4
      p(k56h11,nu)=one/(one-afac*bfac)
     & *(p(5,nu)+p(6,nu)-afac*(p(3,nu)+p(4,nu)))
      p(k34h11,nu)=one/(one-afac*bfac)
     & *(p(3,nu)+p(4,nu)-bfac*(p(5,nu)+p(6,nu)))
      enddo

      endif
c--- end of 3rd generation initialization


c--- set up spinor products (including for flat vectors, for 3rd gen)
      if (docheck) then
        call spinorz(12,p,za,zb)        ! Use spinorz for numerical check
      else
        if (includegen3) then
          call spinoru(12,p,za,zb)
        else
          call spinoru(6,p,za,zb)
        endif
      endif

c--- fill amplitudes with contributions of continuum W pairs
c--- note: 1 = minus, 2=plus

c      Avec(2,2)=a64v(mm,3,4,1,2,6,5,zb,za)*(-half*im)
c      Avec(2,1)=a64v(mp,3,4,1,2,6,5,zb,za)*(-half*im)
c      Avec(1,2)=a64v(pm,3,4,1,2,6,5,zb,za)*(-half*im)
c      Avec(1,1)=a64v(pp,3,4,1,2,6,5,zb,za)*(-half*im)

c      write(6,*) 'Avec(1,1)',Avec(1,1)*2._dp
c      write(6,*) 'Avec(1,2)',Avec(1,2)*2._dp
c      write(6,*) 'Avec(2,1)',Avec(2,1)*2._dp
c      write(6,*) 'Avec(2,2)',Avec(2,2)*2._dp
c      write(6,*)

C A simpler way of doing the same thing

c--- fill amplitudes used for generations 1 and 2
      if (includegens1and2) then
      Avec(2,2)=im*(fvs(pp,4,3,1,2,5,6,za,zb)+fvf(pp,4,3,1,2,5,6,za,zb))
      Avec(2,1)=im*(fvs(pm,4,3,1,2,5,6,za,zb)+fvf(pm,4,3,1,2,5,6,za,zb))
      Avec(1,2)=im*(fvs(pm,3,4,1,2,6,5,zb,za)+fvf(pm,3,4,1,2,6,5,zb,za))
      Avec(1,1)=im*(fvs(pp,3,4,1,2,6,5,zb,za)+fvf(pp,3,4,1,2,6,5,zb,za))
      else
        do h1=1,2
        do h2=1,2
        Avec(h1,h2)=czip
        enddo
        enddo
      endif

c      write(6,*) 'Avec(2,1)',Avec(1,1)
c      write(6,*) 'Avec(2,2)',Avec(1,2)
c      write(6,*) 'Avec(2,1)',Avec(2,1)
c      write(6,*) 'Avec(2,2)',Avec(2,2)

c--- factor two for complete massless isodoublets
      faccont=ctwo

c--- fill amplitudes used for 3rd generation
      if (includegen3) then
        do h1=1,2
        do h2=1,2
        do e=-2,0
        box(h1,h2,e)=czip
        triang(h1,h2,e)=czip
        bub(h1,h2,e)=czip
        enddo
        enddo
        enddo

c--- compute integrals and their coefficients
        if (docheck) then
          call massivebox(1,2,3,4,5,6,za,zb,box)
          call massivetri(1,2,3,4,5,6,za,zb,triang)
          call massivebub(1,2,3,4,5,6,za,zb,bub)
        else
          call massivebox6(1,2,3,4,5,6,za,zb,box)
          call massivetri6(1,2,3,4,5,6,za,zb,triang)
          call massivebub(1,2,3,4,5,6,za,zb,bub)
        endif

c--- This contribution is finite so we only retain "0" piece
        e=0
        do h1=1,2
        do h2=1,2
          sum(h1,h2,e)=box(h1,h2,e)+triang(h1,h2,e)+bub(h1,h2,e)
          Agen3(h1,h2)=sum(h1,h2,0)
        enddo
        enddo
      else
        do h1=1,2
        do h2=1,2
        Agen3(h1,h2)=czip
        enddo
        enddo
      endif

c--- rescale back
      do i=1,6
      do j=1,4
      p(i,j)=p(i,j)*pmax
      enddo
      enddo
      mt=mt*pmax
      musq=musq*pmax**2

c--- fill amplitudes with contributions of Higgs: top loop
      mfsq=mt**2
      tau=s(1,2)/(4._dp*mfsq)
      tauinv=1._dp/tau
      fachiggs=cone/cplx2(s(1,2)-hmass**2,hmass*hwidth)
c--- shat-dependent width
c      fachiggs=cone/cplx2(s(1,2)-hmass**2,hwidth*s(1,2)/hmass)
c--- piece proportional to hmass*hwidth
c      fachiggs=im*imag(fachiggs)
c--- piece proportional to (s(1,2)-hmass**2)
c      fachiggs=real(fachiggs)

      if (tau <= 1._dp) then
         f=cplx1(asin(sqrt(tau))**2)
      elseif (tau > 1._dp) then
         rt=sqrt(1._dp-tauinv)
         f=-0.25_dp*(cplx1(log((1._dp+rt)/(1._dp-rt)))-im*pi)**2
      else
         f=czip
      endif
      e3De4=2._dp*za(3,5)*zb(6,4)/(s(3,4)*s(5,6))
      amphiggs=mfsq*(cone+(cone-cplx1(tauinv))*f)*im*e3De4
      Ahiggs(1,1)=fachiggs*amphiggs*za(1,2)/zb(2,1)
      Ahiggs(1,2)=czip
      Ahiggs(2,1)=czip
      Ahiggs(2,2)=fachiggs*amphiggs*zb(1,2)/za(2,1)

c--- fill amplitudes with contributions of Higgs: bottom loop
      mfsq=mb**2
      tau=s(1,2)/(4._dp*mfsq)
      tauinv=1._dp/tau

      if (tau <= 1._dp) then
         f=cplx1(asin(sqrt(tau))**2)
      elseif (tau > 1._dp) then
         rt=sqrt(1._dp-tauinv)
         f=-0.25_dp*(cplx1(log((1._dp+rt)/(1._dp-rt)))-im*pi)**2
      else
         f=czip
      endif
      amphiggs=mfsq*(cone+(cone-cplx1(tauinv))*f)*im*e3De4

      Ahiggs(1,1)=Ahiggs(1,1)+fachiggs*amphiggs*za(1,2)/zb(2,1)
      Ahiggs(2,2)=Ahiggs(2,2)+fachiggs*amphiggs*zb(1,2)/za(2,1)

      if (docheck) then
c--- numerical check includes top loops only
        Ahiggs(1,1)=Ahiggs(1,1)-fachiggs*amphiggs*za(1,2)/zb(2,1)
        Ahiggs(2,2)=Ahiggs(2,2)-fachiggs*amphiggs*zb(1,2)/za(2,1)
       do h1=1,2
        do h2=1,2
        write(6,*) 'h1,h2, massless',h1,h2,faccont*Avec(h1,h2)
        write(6,*) 'h1,h2, massive ',h1,h2,Agen3(h1,h2)
        write(6,*) 'h1,h2, Higgs   ',h1,h2,Ahiggs(h1,h2)
        write(6,*) 'h1,h2, ---- SUM',h1,h2,
     &               faccont*Avec(h1,h2)+Agen3(h1,h2)+Ahiggs(h1,h2)
        enddo
        enddo
c        pause
      endif

c--- ensure numerical stability: set massive loops to zero
c--- for pt(W) < "ptWsafetycut_massive" GeV
      if (pttwo(3,4,p) < ptWsafetycut_massive) then
        do h1=1,2
        do h2=1,2
          Agen3(h1,h2)=czip
        enddo
        enddo
      endif

c--- ensure numerical stability: set massless loops to zero
c--- for pt(W) < "ptWsafetycut_massless" GeV
      if (pttwo(3,4,p) < ptWsafetycut_massless) then
        do h1=1,2
        do h2=1,2
          Avec(h1,h2)=czip
        enddo
        enddo
      endif

c--- Rescale for width study
      if((keep_smhiggs_norm).and.(anom_higgs)) then
         rescale=chi_higgs**2
         Ahiggs(:,:)=Ahiggs(:,:)*rescale
      endif

      msqgg=0._dp
      do h1=1,2
      do h2=1,2
      Atot(h1,h2)=faccont*Avec(h1,h2)+Agen3(h1,h2)+Ahiggs(h1,h2)

      if     (caseggWW4l) then
c--- This accumulates total contributions
        msqgg=msqgg+abs(Atot(h1,h2))**2
      elseif (caseHWWHpi) then
c--- This only accumulates contributions containing the Higgs diagram,
c---  i.e. the Higgs diagrams squared and the interference
        msqgg=msqgg+abs(Atot(h1,h2))**2
     &             -abs(faccont*Avec(h1,h2)+Agen3(h1,h2))**2
      elseif (caseHWWint) then
c--- This only accumulates the interference
        msqgg=msqgg+abs(Atot(h1,h2))**2
     &             -abs(faccont*Avec(h1,h2)+Agen3(h1,h2))**2
     &             -abs(Ahiggs(h1,h2))**2
      else
        write(6,*) 'Unexpected case in gg_WW_int: ',kcase
        stop
      endif

      enddo
      enddo

c--- overall factor from diagrams
      fac=avegg*V*(2._dp*gwsq*gsq/(16._dp*pisq)*gwsq/2._dp)**2
     & *s(3,4)**2/((s(3,4)-wmass**2)**2+(wwidth*wmass)**2)
     & *s(5,6)**2/((s(5,6)-wmass**2)**2+(wwidth*wmass)**2)
      fac=fac/pmax**4

      msq(0,0)=msqgg*fac
      return
      end


