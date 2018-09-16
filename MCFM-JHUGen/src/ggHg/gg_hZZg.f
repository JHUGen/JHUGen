      subroutine gg_hZZg(p,msq)
      implicit none
      include 'types.f'

C-----Author R.K. Ellis December 2009
c---- matrix element for H production
c----in the heavy quark (mt=Infinity) limit.
C----averaged over initial colours and spins
C---- f(p1)+f(p2) --> H(-->Z^+(e^-(p3),e^+(p4))Z(mu^-(p5),mu^+(p6)))+f(p7)'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'zcouple.f'
      integer:: j,k,iglue
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4)
      real(dp):: ss,tt,uu,mhsq,hdecay,s(mxpart,mxpart)
      real(dp):: s3456,Asq,fac
      parameter(iglue=7)

c--set msq=0 to initialize
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0._dp
      enddo
      enddo

      call dotem(iglue,p,s)

      s3456=s(3,4)+s(3,5)+s(3,6)+s(4,5)+s(4,6)+s(5,6)
      hdecay=gwsq**3*zmass**2*4._dp*xw**2/(one-xw)*
     & ( ((l1*l2)**2+(r1*r2)**2)*s(3,5)*s(4,6)
     &  +((r1*l2)**2+(r2*l1)**2)*s(3,6)*s(4,5))
      hdecay=hdecay/((s(3,4)-zmass**2)**2+(zmass*zwidth)**2)
      hdecay=hdecay/((s(5,6)-zmass**2)**2+(zmass*zwidth)**2)
      hdecay=hdecay/((s3456-hmass**2)**2+(hmass*hwidth)**2)

      Asq=(as/(3._dp*pi))**2/vevsq
      fac=Asq*gsq*hdecay

c--   calculate propagators
      ss=s(1,2)
      tt=s(1,iglue)
      uu=s(2,iglue)
      mhsq=ss+tt+uu

      msq(0,0)=
     & avegg*fac*V*xn*(mhsq**4+ss**4+tt**4+uu**4)/(ss*tt*uu)
      msq(1,-1)=+aveqq*fac*V/2._dp*(tt**2+uu**2)/ss
      msq(0,+1)=-aveqg*fac*V/2._dp*(ss**2+tt**2)/uu
      msq(+1,0)=-aveqg*fac*V/2._dp*(ss**2+uu**2)/tt

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
