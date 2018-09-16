      subroutine gg_hgagag_gvec(p,n,in,msq)
      implicit none
      include 'types.f'
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
C  in is the label of the momentum contracted with n
      integer:: j,k,in,iglue
      real(dp):: msq(-nf:nf,-nf:nf),s34
      real(dp):: n(4),p(mxpart,4),dot,hdecay,fac,
     & qqghn,ggghn,p1p2(-1:1,-1:1),msqgamgam
      parameter(iglue=5)
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0._dp
      enddo
      enddo

C   Deal with Higgs decay
      s34=2._dp*dot(p,3,4)
      hdecay=msqgamgam(hmass)/((s34-hmass**2)**2+(hmass*hwidth)**2)

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
      elseif (in == 5) then     
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

