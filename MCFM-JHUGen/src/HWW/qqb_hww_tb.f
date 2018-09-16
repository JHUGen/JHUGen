      subroutine qqb_hww_tb(p,msq)
      implicit none
      include 'types.f'

c--- Author: J. M. Campbell, June 2011
c--- Implementation of gg -> H -> WW -> leptons,
c--- including both top and bottom quark loops
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'masses.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      include 'anom_higgs.f'
      integer:: h1,h2,j,k
      real(dp):: p(mxpart,4),msq(fn:nf,fn:nf),msqgg,fac
      real(dp):: mfsq,tau,tauinv,rt,rescale
      complex(dp):: Ahiggs(2,2),fachiggs,amphiggs,f,e3De4
!      complex(dp):: num_c

      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0._dp
      enddo
      enddo

      call spinoru(6,p,za,zb)

c--- fill amplitudes with contributions of Higgs: top loop
      mfsq=mt**2
      tau=s(1,2)/(4._dp*mfsq)
      tauinv=1._dp/tau
      fachiggs=cone/cplx2(s(1,2)-hmass**2,hmass*hwidth)

!====== SEYMOUR ISA APPROX
!      num_c=cplx1(hmass**2/s(1,2))
!      fachiggs=num_c/cplx2(s(1,2)-hmass**2,s(1,2)*hwidth/hmass)

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

c--- Rescale for width study
      if((keep_smhiggs_norm).and.(anom_higgs)) then
         rescale=chi_higgs**2
         Ahiggs(:,:)=Ahiggs(:,:)*rescale
      endif

      msqgg=0._dp
      do h1=1,2
      do h2=1,2
      msqgg=msqgg+abs(Ahiggs(h1,h2))**2
      enddo
      enddo

c--- overall factor from diagrams
      fac=avegg*V*(2._dp*gwsq*gsq/(16._dp*pisq)*gwsq/2._dp)**2
     & *s(3,4)**2/((s(3,4)-wmass**2)**2+(wwidth*wmass)**2)
     & *s(5,6)**2/((s(5,6)-wmass**2)**2+(wwidth*wmass)**2)

      msq(0,0)=msqgg*fac

      return
      end


