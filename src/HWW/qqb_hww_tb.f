      subroutine qqb_hww_tb(p,msq)
      implicit none
c--- Author: J. M. Campbell, June 2011
c--- Implementation of gg -> H -> WW -> leptons,
c--- including both top and bottom quark loops
      include 'constants.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'masses.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      include 'anom_higgs.f' 
      integer h1,h2,j,k
      double precision p(mxpart,4),msq(fn:nf,fn:nf),msqgg,fac
      double precision mfsq,tau,tauinv,rt,rescale
      double complex Ahiggs(2,2),fachiggs,amphiggs,f,e3De4
!      double complex num_c

      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      enddo
      enddo

      call spinoru(6,p,za,zb)
      
c--- fill amplitudes with contributions of Higgs: top loop
      mfsq=mt**2
      tau=s(1,2)/(4d0*mfsq)
      tauinv=1d0/tau
      fachiggs=cone/dcmplx(s(1,2)-hmass**2,hmass*hwidth)

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
      e3De4=2d0*za(3,5)*zb(6,4)/(s(3,4)*s(5,6))
      amphiggs=mfsq*(cone+(cone-dcmplx(tauinv))*f)*im*e3De4
      Ahiggs(1,1)=fachiggs*amphiggs*za(1,2)/zb(2,1)
      Ahiggs(1,2)=czip
      Ahiggs(2,1)=czip
      Ahiggs(2,2)=fachiggs*amphiggs*zb(1,2)/za(2,1)

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
      amphiggs=mfsq*(cone+(cone-dcmplx(tauinv))*f)*im*e3De4

      Ahiggs(1,1)=Ahiggs(1,1)+fachiggs*amphiggs*za(1,2)/zb(2,1)
      Ahiggs(2,2)=Ahiggs(2,2)+fachiggs*amphiggs*zb(1,2)/za(2,1)

c--- Rescale for width study
      if((keep_smhiggs_norm).and.(anom_higgs)) then 
         rescale=chi_higgs**2 
         Ahiggs(:,:)=Ahiggs(:,:)*rescale
      endif

      msqgg=0d0
      do h1=1,2
      do h2=1,2
      msqgg=msqgg+cdabs(Ahiggs(h1,h2))**2
      enddo
      enddo

c--- overall factor from diagrams
      fac=avegg*V*(2d0*gwsq*gsq/(16d0*pisq)*gwsq/2d0)**2
     & *s(3,4)**2/((s(3,4)-wmass**2)**2+(wwidth*wmass)**2)
     & *s(5,6)**2/((s(5,6)-wmass**2)**2+(wwidth*wmass)**2)
      
      msq(0,0)=msqgg*fac

      return
      end
      
      
