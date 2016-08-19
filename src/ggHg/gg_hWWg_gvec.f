      subroutine gg_hWWg_gvec(p,n,in,msq)
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'ewcouple.f'
C  in is the label of the momentum contracted with n
      integer j,k,in,iglue
      double precision msq(-nf:nf,-nf:nf),s34,s35,s36,s45,s46,s56
      double precision n(4),p(mxpart,4),dot,hdecay,s3456,fac,
     . qqghn,ggghn,p1p2(-1:1,-1:1)
      parameter(iglue=7)
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      enddo
      enddo

C   Deal with Higgs decay
      s34=2d0*dot(p,3,4)
      s35=2d0*dot(p,3,5)
      s36=2d0*dot(p,3,6)
      s45=2d0*dot(p,4,5)
      s46=2d0*dot(p,4,6)
      s56=2d0*dot(p,5,6)
      s3456=s34+s35+s36+s45+s46+s56
      hdecay=gwsq**3*wmass**2*s35*s46
      hdecay=hdecay/(((s3456-hmass**2)**2+(hmass*hwidth)**2)
     .   *((s34-wmass**2)**2+(wmass*wwidth)**2)
     .   *((s56-wmass**2)**2+(wmass*wwidth)**2))

      fac=hdecay

      do j=-1,+1
      do k=-1,+1
      p1p2(j,k)=0d0
      enddo
      enddo

      if (in .eq. 1) then
      p1p2(0,-1)=-aveqg*fac*qqghn(2,iglue,1,p,n)
      p1p2(0,+1)=-aveqg*fac*qqghn(2,iglue,1,p,n)
      p1p2(0,0)=+avegg*fac*ggghn(iglue,2,1,p,n)
      elseif (in .eq. 2) then
      p1p2(+1,0)=-aveqg*fac*qqghn(1,iglue,2,p,n)
      p1p2(-1,0)=-aveqg*fac*qqghn(iglue,1,2,p,n)
      p1p2(0,0)=+avegg*fac*ggghn(1,iglue,2,p,n)
      elseif (in .eq. 7) then     
      p1p2(1,-1)=+aveqq*fac*qqghn(1,2,iglue,p,n)
      p1p2(-1,1)=+aveqq*fac*qqghn(2,1,iglue,p,n)
      p1p2(0,0)=+avegg*fac*ggghn(1,2,iglue,p,n)
      endif

      do j=-nf,nf
      do k=-nf,nf
      if     ((j .gt. 0) .and. (k .eq. -j)) then
          msq(j,k)=p1p2(1,-1)
      elseif ((j .lt. 0) .and. (k .eq. -j)) then
          msq(j,k)=p1p2(-1,1)
      elseif ((j .eq. 0) .and. (k .eq. 0)) then
          msq(j,k)=p1p2(0,0)
      elseif ((j .gt. 0) .and. (k .eq. 0)) then
          msq(j,k)=
     &    p1p2(+1,0)
      elseif ((j .lt. 0) .and. (k .eq. 0)) then
          msq(j,k)=
     &    p1p2(-1,0)
      elseif ((j .eq. 0) .and. (k .gt. 0)) then
          msq(j,k)=
     &    p1p2(0,+1)
      elseif ((j .eq. 0) .and. (k .lt. 0)) then
          msq(j,k)=
     &    p1p2(0,-1)
      endif
      enddo
      enddo
 
      return
      end

