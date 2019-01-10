      subroutine gen3from2(q,z,rtalpha,phit,p,jac,*)
c----jac is the total wt of the whole business (2) and (3 from 2)
c----q are the input momenta
c----p are the output momenta
      implicit none
      include 'constants.f'
      include 'debug.f'
      integer nmin,nmax,j,iseed,k
      double precision p(mxpart,4),q(mxpart,4),z,rtalpha,phit,
     . wt3_2,msq(-nf:nf,-nf:nf)
      double precision sum(0:8),wtc(8),apweight(8),jac,ran0,myran
      common/apwt/apweight
      common/nmin/nmin
      common/nmax/nmax
      data iseed/1768/
      integer,parameter:: i1(8)=(/1,2,4,5,1,1,2,2/)
      integer,parameter:: i2(8)=(/2,1,5,4,4,5,4,5/)


      do j=1,7
      do k=1,4
      p(j,k)=q(j,k)
      enddo
      enddo
 
      sum(nmin-1)=0d0

      do j=nmin,nmax
      apweight(j)=1d0/dfloat(nmax-nmin+1)
      sum(j)=sum(j-1)+apweight(j)
      if (debug) then
      write(6,*) 'j',j 
      write(6,*) 'apweight(j)',apweight(j) 
      write(6,*) 'sum(j)',sum(j) 
      endif
      enddo
  
      myran=ran0(iseed)

      do j=nmin,nmax
      if ((myran .gt. sum(j-1)) .and. (myran .lt. sum(j))) then
c---genrad is a switchyard routine routing to genrii,genrif,genrff
c---genrad modifies the vector p to provide new ones
      call genrad(p,i1(j),i2(j),5,z,rtalpha,phit,wt3_2,*999)
c---although genrad returns wt3_2 we shall not use it 
c---in this step we have generated the new p's (only one set)
c---only one option is pursued in this do-loop 
      endif      
      enddo
      if (debug) then 
      write(6,*) 'wt3_2 in gen3from2.f',wt3_2
      call writeout(p)
      endif

c---Sum over channels
c---Initialize jac
      jac=0d0
      do j=nmin,nmax
         if ((j .eq. 1) .or. (j .eq. 2)) 
     .      call genii(j,p,wtc(j),msq)
         if ((j .eq. 3) .or. (j .eq. 4)) 
     .      call genff(j-2,p,wtc(j),msq)
         if ((j .eq. 5) .or. (j .eq. 6) .or. (j .eq. 7) .or. (j .eq. 8))
     .      call genif(j-4,p,wtc(j),msq)
c        jac=jac+apweight(j)/wtc(j)
        jac=jac+1d0/wtc(j)
      enddo
      jac=1d0/jac

      if (debug) write(6,*) 
      if (debug) write(6,*) 'this is the result of reconstruction'
      if (debug) write(6,*) 'jac in gen3from2',jac
c      if (debug) pause
      return 

 999  jac=0d0
      return 1
      end

