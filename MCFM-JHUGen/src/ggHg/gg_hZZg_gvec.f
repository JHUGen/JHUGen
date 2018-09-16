      subroutine gg_hZZg_gvec(p,n,in,msq)
      implicit none
      include 'types.f'

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'ewcouple.f'
      include 'zcouple.f'
C  in is the label of the momentum contracted with n
      integer:: j,k,in,iglue
      real(dp):: msq(-nf:nf,-nf:nf),s34,s35,s36,s45,s46,s56
      real(dp):: n(4),p(mxpart,4),dot,hdecay,s3456,fac,
     & qqghn,ggghn,p1p2(-1:1,-1:1)
      parameter(iglue=7)
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0._dp
      enddo
      enddo

C   Deal with Higgs decay
      s34=2._dp*dot(p,3,4)
      s35=2._dp*dot(p,3,5)
      s36=2._dp*dot(p,3,6)
      s45=2._dp*dot(p,4,5)
      s46=2._dp*dot(p,4,6)
      s56=2._dp*dot(p,5,6)
      s3456=s34+s35+s36+s45+s46+s56
      hdecay=gwsq**3*zmass**2*4._dp*xw**2/(one-xw)*
     & ( ((l1*l2)**2+(r1*r2)**2)*s35*s46
     &  +((r1*l2)**2+(r2*l1)**2)*s36*s45)
      hdecay=hdecay/((s34-zmass**2)**2+(zmass*zwidth)**2)
      hdecay=hdecay/((s56-zmass**2)**2+(zmass*zwidth)**2)
      hdecay=hdecay/((s3456-hmass**2)**2+(hmass*hwidth)**2)

      fac=hdecay

      do j=-1,+1
      do k=-1,+1
      p1p2(j,k)=0._dp
      enddo
      enddo

      if (in == 1) then
      p1p2(0,-1)=-aveqg*fac*qqghn(2,iglue,1,p,n)
      p1p2(0,+1)=-aveqg*fac*qqghn(2,iglue,1,p,n)
      p1p2(0,0)=+avegg*fac*ggghn(iglue,2,1,p,n)
      elseif (in == 2) then
      p1p2(+1,0)=-aveqg*fac*qqghn(1,iglue,2,p,n)
      p1p2(-1,0)=-aveqg*fac*qqghn(iglue,1,2,p,n)
      p1p2(0,0)=+avegg*fac*ggghn(1,iglue,2,p,n)
      elseif (in == 7) then
      p1p2(1,-1)=+aveqq*fac*qqghn(1,2,iglue,p,n)
      p1p2(-1,1)=+aveqq*fac*qqghn(2,1,iglue,p,n)
      p1p2(0,0)=+avegg*fac*ggghn(1,2,iglue,p,n)
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
          msq(j,k)=
     &    p1p2(+1,0)
      elseif ((j < 0) .and. (k == 0)) then
          msq(j,k)=
     &    p1p2(-1,0)
      elseif ((j == 0) .and. (k > 0)) then
          msq(j,k)=
     &    p1p2(0,+1)
      elseif ((j == 0) .and. (k < 0)) then
          msq(j,k)=
     &    p1p2(0,-1)
      endif
      enddo
      enddo

      return
      end

