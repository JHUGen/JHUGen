      subroutine genff(nperms,p,wt,msq)
      implicit none
      include 'types.f'
c--- final-final subtraction.
c--- nperms is an argument between 1 and 2
c--- this routine calculates the jacobian associated with calculating
c--- a given final state, by contracting a final-final dipole

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'qcdcouple.f'
      include 'debug.f'
      integer:: i4,i5,j,k,nperms
      real(dp):: p(mxpart,4),z,dot,q(mxpart,4),
     & msq(-nf:nf,-nf:nf),wt4,wt5_4
      real(dp):: facq,wt,s3i4,y,omy,jacbit,wt0
      parameter(wt0=1._dp/eight/pisq)
      integer,parameter:: j4(2)=(/4,5/)
      integer,parameter:: j5(2)=(/5,4/)

      i4=j4(nperms)
      i5=j5(nperms)

c      write(6,*) 'genff:nperms',nperms
c      write(6,*) 'genff:i4',i4
c      write(6,*) 'genff:i5',i5

c first of calculate the variable's which one started with
c---nb all incoming
      s3i4=two*dot(p,3,i4)
      y=dot(p,3,i4)/(dot(p,3,i4)+dot(p,3,i5)+dot(p,i4,i5))
      z=dot(p,3,i5)/(dot(p,3,i5)+dot(p,i4,i5))
      omy=one-y

      do j=1,4
      q(1,j)=p(1,j)
      q(2,j)=p(2,j)
      q(3,j)=p(3,j)
      q(6,j)=p(6,j)
      q(7,j)=p(7,j)
      q(i4,j)=p(i4,j)+p(3,j)-p(i5,j)*y/omy
      q(i5,j)=p(i5,j)/omy
      enddo

      jacbit=four*sqrt(y)/(0.5_dp/sqrt(z)+0.5_dp/sqrt(1._dp-z))
      wt5_4=wt0*omy*dot(q,4,5)*jacbit

      call wt4gen(q,wt4)
c---calculate total weight
      wt=wt4*wt5_4

      if (debug) then
      write(6,*) 'wt5_4 in genff',wt5_4
      write(6,*) 'wt4 in genff',wt4
      write(6,*) 'wt in genff',wt
      endif


      call qqb_wbb(q,msq)


      facq=-gsq/s3i4*(two/(one-z*(one-y))-one-z)

      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0._dp
      if     ((j > 0) .and. (k < 0)) then
      msq(j,k)=facq*msq(j,k)
      elseif ((j < 0) .and. (k > 0)) then
      msq(j,k)=facq*msq(j,k)
      endif
      enddo
      enddo

      return
      end

