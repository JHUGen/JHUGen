      subroutine genif(nperms,p,wt,msq)
      implicit none
      include 'types.f'
c----initial-final subtraction.
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'qcdcouple.f'
      include 'xmin.f'
      include 'debug.f'
      logical:: justjac

      integer:: i1,i2,i3,i4,i5,j,k,nperms
      real(dp):: p(mxpart,4),u,z,dot,q(mxpart,4),
     & msq(-nf:nf,-nf:nf),omxmin
      real(dp):: facq,si1i3,si1i4,si3i4,x,omx,wt0,
     & wt5_4,wt4,wt,jacbit
      parameter(wt0=one/eight/pisq)
      common/justjac/justjac
      integer,parameter:: j1(4)=(/1,1,2,2/)
      integer,parameter:: j2(4)=(/2,2,1,1/)
      integer,parameter:: j3(4)=(/3,3,3,3/)
      integer,parameter:: j4(4)=(/4,5,4,5/)
      integer,parameter:: j5(4)=(/5,4,5,4/)

      i1=j1(nperms)
      i2=j2(nperms)
      i3=j3(nperms)
      i4=j4(nperms)
      i5=j5(nperms)

c first of calculate the variable's which one started with
c---nb all incoming
      si1i3=two*dot(p,i1,i3)
      si1i4=two*dot(p,i1,i4)
      si3i4=two*dot(p,i3,i4)
c----note momenta all incoming
      x=one+si3i4/(si1i3+si1i4)
      u=si1i3/(si1i3+si1i4)
      z=one-u
      omx=one-x
      omxmin=one-xmin
c -xx(i1)
c---at this stage the p are momenta including radiation

      do j=1,4
      q(i1,j)=x*p(i1,j)
      q(i2,j)=p(i2,j)
      q(i3,j)=p(i3,j)
      q(i4,j)=p(i4,j)+p(i3,j)+omx*p(i1,j)
      q(i5,j)=p(i5,j)
      q(6,j)=p(6,j)
      q(7,j)=p(7,j)
      enddo

      jacbit=four*sqrt(omx)/(half/sqrt(z)+half/sqrt(u))
      wt5_4=-wt0*dot(q,i1,i4)/x**2*omxmin*jacbit
      call wt4gen(q,wt4)

c---calculate total weight
      wt=wt5_4*wt4



      if (debug) then
      write(6,*) 'x in genif',x
      write(6,*) 'omxmin in genif',omxmin
      write(6,*) 'i1 in genif',i1
      write(6,*) 'i2 in genif',i2
      write(6,*) 'i3 in genif',i3
      write(6,*) 'i4 in genif',i4
      write(6,*) 'i5 in genif',i5
      write(6,*) 'wt5_4 in genif',wt5_4
      write(6,*) 'wt4 in genif',wt4
      write(6,*) 'wt in genif',wt
      endif

      if (justjac) return


      call qqb_wbb(q,msq)


c----initial wrt to final
      facq=-gsq/(x*si1i3)*(two/(one-x+u)-one-x)
     &     -gsq/(x*si3i4)*(two/(one-x+u)-one-z)

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
