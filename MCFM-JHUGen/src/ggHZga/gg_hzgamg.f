      subroutine gg_hzgamg(p,msq)
      implicit none
c----Author: R.K. Ellis July 2012
c----matrix element for H production and decay
c----in the heavy quark (mt=Infinity) limit.
C----averaged over initial colours and spins
c    g(-p1)+g(-p2)-->H --> Z/gamma^*--(l(p3)+a(p4)) + gamma(p5) +g(p6) 
      include 'constants.f'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'sprods_com.f'
      integer j,k,iglue
      double precision msq(-nf:nf,-nf:nf),p(mxpart,4)
      double precision ss,tt,uu,shsq,hdecay,Asq,fac,HZgamMSQ
      parameter(iglue=6)

      msq(:,:)=0d0

C   Deal with Higgs decay to Zgam
      call dotem(iglue,p,s)

c--   calculate propagators
      ss=s(1,2)
      tt=s(1,iglue)
      uu=s(2,iglue)
      shsq=ss+tt+uu
      
      hdecay=HZgamMSQ(3,4,5)
     & /((shsq-hmass**2)**2+(hmass*hwidth)**2)


      Asq=(as/(3d0*pi))**2/vevsq
      fac=Asq*gsq*hdecay

      msq(0,0)=
     . avegg*fac*V*xn*(shsq**4+ss**4+tt**4+uu**4)/(ss*tt*uu) 
      msq(1,-1)=+aveqq*fac*V/2d0*(tt**2+uu**2)/ss
      msq(0,+1)=-aveqg*fac*V/2d0*(ss**2+tt**2)/uu
      msq(+1,0)=-aveqg*fac*V/2d0*(ss**2+uu**2)/tt

c--set msq=0 to initialize
      do j=-nf,nf
      do k=-nf,nf
      if ((k .eq. -j) .and. (j .ne. 0)) then
      msq(j,k)=msq(1,-1)
      elseif ((j .eq. 0) .and. (k .ne. 0)) then
      msq(j,k)=msq(0,1)
      elseif ((j .ne. 0) .and. (k .eq. 0)) then
      msq(j,k)=msq(1,0)
      endif
      enddo
      enddo
      
      return
      end
