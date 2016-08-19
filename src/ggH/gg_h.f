      subroutine gg_h(p,msq)
      implicit none
c----Lowest order matrix element for H production
c----in the heavy quark (mt=Infinity) limit.
C----averaged over initial colours and spins
c     g(-p1)+g(-p2)-->H -->  b(p3)+b(p4))
      include 'constants.f'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'hdecaymode.f'
      integer j,k
      double precision msq(-nf:nf,-nf:nf),p(mxpart,4),s,s12
      double precision hdecay,gg,Asq,hprod,msqgamgam
c      double precision tn,ftn
      s(j,k)=2d0*(p(j,4)*p(k,4)-p(j,1)*p(k,1)
     .           -p(j,2)*p(k,2)-p(j,3)*p(k,3))

c--- set msq=0 to initialize
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      enddo
      enddo
      s12=s(1,2)
      
C   Deal with Higgs decay
      if (hdecaymode == 'tlta') then
          call htautaudecay(p,3,4,hdecay)
      elseif (hdecaymode == 'bqba') then
          call hbbdecay(p,3,4,hdecay)
      elseif (hdecaymode == 'gaga') then
          hdecay=msqgamgam(hmass)
      else
      write(6,*) 'Unimplemented process in gg_hgg_v'
      stop
      endif
      hdecay=hdecay/((s12-hmass**2)**2+(hmass*hwidth)**2)

      hprod=1d0
C--- The following lines may be uncommented in order to implement
c--- the correct dependence on the top quark mass      
c      tn = 4d0*(mt/hmass)**2
c      if (tn.lt.(1.0d0)) then
c         ftn = 0.5d0*((dlog((1d0+dsqrt(1d0-tn))
c     &        /(1d0-dsqrt(1d0-tn))))**2-pisq)
c      else
c         ftn = -2d0*(dasin(1.0d0/dsqrt(tn)))**2
c      endif
c      hprod=((3d0*tn/4d0)*(2d0+(tn-1d0)*ftn))**2
C---- end of mass corrections in production
c      write(6,*) 'hprod',hprod

      Asq=(as/(3d0*pi))**2/vevsq
      gg=0.5d0*Asq*V*s12**2
      msq(0,0)=avegg*gg*hdecay*hprod
            
      return
      end
