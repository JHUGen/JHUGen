      subroutine gg_hg_gvec(p,n,in,msq)
      implicit none
      include 'types.f'

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'hdecaymode.f'
C  in is the label of the momentum contracted with n
      integer:: j,k,in
      real(dp):: msq(-nf:nf,-nf:nf)
      real(dp):: n(4),p(mxpart,4),hdecay,s34,fac,
     & qqghn,ggghn,p1p2(-1:1,-1:1),msqhgamgam

      msq(:,:)=zip

      s34=(p(3,4)+p(4,4))**2
     &   -(p(3,1)+p(4,1))**2-(p(3,2)+p(4,2))**2-(p(3,3)+p(4,3))**2

C   Deal with Higgs decay
      if (hdecaymode == 'tlta') then
          call htautaudecay(p,3,4,hdecay)
      elseif (hdecaymode == 'bqba') then
          call hbbdecay(p,3,4,hdecay)
      elseif (hdecaymode == 'gaga') then
           hdecay=msqhgamgam(s34)
      else
      write(6,*) 'Unimplemented process in gg_hgg_gvec'
      stop
      endif
      hdecay=hdecay/((s34-hmass**2)**2+(hmass*hwidth)**2)

      fac=hdecay

      p1p2(:,:)=zip

      if (in == 1) then
      p1p2(0,-1)=-aveqg*fac*qqghn(2,5,1,p,n)
      p1p2(0,+1)=-aveqg*fac*qqghn(2,5,1,p,n)
      p1p2(0,0)=+avegg*fac*ggghn(5,2,1,p,n)
      elseif (in == 2) then
      p1p2(+1,0)=-aveqg*fac*qqghn(1,5,2,p,n)
      p1p2(-1,0)=-aveqg*fac*qqghn(5,1,2,p,n)
      p1p2(0,0)=+avegg*fac*ggghn(1,5,2,p,n)
      elseif (in == 5) then
      p1p2(1,-1)=+aveqq*fac*qqghn(1,2,5,p,n)
      p1p2(-1,1)=+aveqq*fac*qqghn(2,1,5,p,n)
      p1p2(0,0)=+avegg*fac*ggghn(1,2,5,p,n)
      endif

      do j=-nf,nf
      do k=-nf,nf
      if     ((j > 0) .and. (k == -j)) then
          msq(j,k)=p1p2(1,-1)
      elseif ((j < 0) .and. (k == -j)) then
          msq(j,k)=p1p2(-1,1)
      elseif ((j == 0) .and. (k == 0)) then
          msq(j,k)=p1p2(0,0)
      elseif ((j > 0) .and. (k == 0)) then
          msq(j,k)=p1p2(+1,0)
      elseif ((j < 0) .and. (k == 0)) then
          msq(j,k)=p1p2(-1,0)
      elseif ((j == 0) .and. (k > 0)) then
          msq(j,k)=p1p2(0,+1)
      elseif ((j == 0) .and. (k < 0)) then
          msq(j,k)=p1p2(0,-1)
      endif
      enddo
      enddo

      return
      end

      function qqghn(j1,j2,j5,p,n)
      implicit none
      include 'types.f'
      real(dp):: qqghn

C---calculates the amplitude squared for the process
c   q(p1)+qbar(p2) --> H((p3+p4)+g(p5)
c   contracted with the vector n(mu)
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      integer:: j1,j2,j5
      real(dp):: Asq,p(mxpart,4),n(4),nDn,nDp1,nDp2,dot,s,t,u

      nDp1=n(4)*p(j1,4)-n(3)*p(j1,3)-n(2)*p(j1,2)-n(1)*p(j1,1)
      nDp2=n(4)*p(j2,4)-n(3)*p(j2,3)-n(2)*p(j2,2)-n(1)*p(j2,1)
      nDn=n(4)**2-n(3)**2-n(2)**2-n(1)**2

      call checkndotp(p,n,j5)

      Asq=(as/(three*pi))**2/vevsq

      s=two*Dot(p,j1,j2)
      t=two*Dot(p,j1,j5)
      u=two*Dot(p,j2,j5)

c--- RKE Schoonship answer, ncalc.out
c     Id,bit=g^2*A^2*V/2*(2*(nDp1*u-nDp2*t)^2/s^2
c           + 0.5*nDn*(u+t)^2/s)

      qqghn=-Asq*gsq*V/two*(two*(nDp1*u-nDp2*t)**2/s**2
     &                    +0.5_dp*nDn*(u+t)**2/s)

      return
      end

      function ggghn(j1,j2,j5,p,n)
      implicit none
      include 'types.f'
      real(dp):: ggghn

C---calculates the amplitude squared for the process
c   g(p1)+g(p2) --> H((p3+p4)+g(p5)
c   contracted with the vector n(mu)
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      integer:: j1,j2,j5
      real(dp):: Asq,p(mxpart,4),n(4),nDn,nDp1,nDp2,dot,s,t,u,sh

      nDp1=n(4)*p(j1,4)-n(3)*p(j1,3)-n(2)*p(j1,2)-n(1)*p(j1,1)
      nDp2=n(4)*p(j2,4)-n(3)*p(j2,3)-n(2)*p(j2,2)-n(1)*p(j2,1)
      nDn=n(4)**2-n(3)**2-n(2)**2-n(1)**2

      call checkndotp(p,n,j5)

      Asq=(as/(three*pi))**2/vevsq

      s=two*Dot(p,j1,j2)
      t=two*Dot(p,j1,j5)
      u=two*Dot(p,j2,j5)
      sh=s+t+u

c--- JMC answer, gggH.frm
c -f(b,a,c)^2/s12/s13/s23*(
c -(n.n)/2*(s12^4+s13^4+s23^4+mHsq^4-2*(s13^2*s23^2+s12^2*mHsq^2))
c +2*(p1.n*s23-p2.n*s13)^2/s13/s23*(s13^2*s23^2+s12^2*mHsq^2)/s12);

      ggghn=Asq*gsq*V*xn*(
     & -nDn/two*(s**4+t**4+u**4+sh**4-two*(t**2*u**2+s**2*sh**2))
     & +two*(nDp1*u-nDp2*t)**2/t/u*(t**2*u**2+s**2*sh**2)/s)/s/t/u

      return
      end

