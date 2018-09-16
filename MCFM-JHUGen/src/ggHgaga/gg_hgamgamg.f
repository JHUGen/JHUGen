      subroutine gg_hgamgamg(p,msq)
      implicit none
      include 'types.f'
      
c---- matrix element for H production
c----in the heavy quark (mt=Infinity) limit.
C----averaged over initial colours and spins
c     g(-p1)+g(-p2)-->H(-->  b(p3)+b~(p4))+g(p5)
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      integer:: j,k,iglue
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4)
      real(dp):: ss,tt,uu,mhsq,hdecay,s(mxpart,mxpart)
      real(dp):: Asq,fac,msqgamgam
      parameter(iglue=5)

      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0._dp
      enddo
      enddo

      call dotem(iglue,p,s)

C   Deal with Higgs decay to gamgam
      hdecay=msqgamgam(hmass)/((s(3,4)-hmass**2)**2+(hmass*hwidth)**2)

c--   calculate propagators
      ss=s(1,2)
      tt=s(1,iglue)
      uu=s(2,iglue)
      mhsq=ss+tt+uu

      Asq=(as/(3._dp*pi))**2/vevsq
      fac=Asq*gsq*hdecay

      msq(0,0)=
     & avegg*fac*V*xn*(mhsq**4+ss**4+tt**4+uu**4)/(ss*tt*uu) 
      msq(1,-1)=+aveqq*fac*V/2._dp*(tt**2+uu**2)/ss
      msq(0,+1)=-aveqg*fac*V/2._dp*(ss**2+tt**2)/uu
      msq(+1,0)=-aveqg*fac*V/2._dp*(ss**2+uu**2)/tt

c--set msq=0 to initialize
      do j=-nf,nf
      do k=-nf,nf
      if ((k == -j) .and. (j .ne. 0)) then
      msq(j,k)=msq(1,-1)
      elseif ((j == 0) .and. (k .ne. 0)) then
      msq(j,k)=msq(0,1)
      elseif ((j .ne. 0) .and. (k == 0)) then
      msq(j,k)=msq(1,0)
      endif
      enddo
      enddo
      
      return
      end
