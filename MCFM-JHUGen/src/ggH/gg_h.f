      subroutine gg_h(p,msq)
      implicit none
      include 'types.f'

c----Lowest order matrix element for H production
c----in the heavy quark (mt=Infinity) limit.
C----averaged over initial colours and spins
c     g(-p1)+g(-p2)-->H -->  b(p3)+b(p4))
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'hdecaymode.f'
      integer:: j,k
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4),s,s12
      real(dp):: hdecay,gg,Asq,hprod,msqhgamgam
c      real(dp):: tn,ftn
      s(j,k)=two*(p(j,4)*p(k,4)-p(j,1)*p(k,1)
     &           -p(j,2)*p(k,2)-p(j,3)*p(k,3))

c--- set msq=0 to initialize
      msq(:,:)=zip

      s12=s(1,2)

C   Deal with Higgs decay
      if (hdecaymode == 'tlta') then
          call htautaudecay(p,3,4,hdecay)
      elseif (hdecaymode == 'bqba') then
          call hbbdecay(p,3,4,hdecay)
      elseif (hdecaymode == 'gaga') then
          hdecay=msqhgamgam(s12)
      else
      write(6,*) 'Unimplemented process in gg_hgg_v'
      stop
      endif
      hdecay=hdecay/((s12-hmass**2)**2+(hmass*hwidth)**2)

      hprod=one
C--- The following lines may be uncommented in order to implement
c--- the correct dependence on the top quark mass
c      tn = 4._dp*(mt/hmass)**2
c      if (tn<(1.0_dp)) then
c         ftn = 0.5_dp*((log((1._dp+sqrt(1._dp-tn))
c     &        /(1._dp-sqrt(1._dp-tn))))**2-pisq)
c      else
c         ftn = -2._dp*(asin(1.0_dp/sqrt(tn)))**2
c      endif
c      hprod=((3._dp*tn/4._dp)*(2._dp+(tn-1._dp)*ftn))**2
C---- end of mass corrections in production
c      write(6,*) 'hprod',hprod

      Asq=(as/(three*pi))**2/vevsq
      gg=half*Asq*V*s12**2
      msq(0,0)=avegg*gg*hdecay*hprod

      return
      end
