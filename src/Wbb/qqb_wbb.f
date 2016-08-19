      subroutine qqb_wbb(p,msq)
c---  Matrix elements squared 
c     q(-p1)+qb(-p2) --> nu(p3)+e^+(p4)+b(p5)+bb(p6)
c---  averaged(summed) over initial(final) colours and spins
      implicit none
      include 'constants.f'
      include 'ckm.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      include 'masses.f'
      include 'first.f'
      integer j,k
      double precision p(mxpart,4),msq(-nf:nf,-nf:nf),msqwbb
      double precision qqb,qbq
      
      if (first) then
       write(6,*)
       write(6,*) '****************** Process info ********************'
       write(6,*) '*                                                  *'
       write(6,*) '* mb=0 for this process, although cuts are applied *'
       write(6,*) '* to simulate the effect of the b-mass:            *'
       write(6,*) '*                                                  *'
       write(6,99) ' *                pt(b) > ',dsqrt(mbsq),
     .  '                *'
       write(6,99) ' *                m(bb) > ',two*dsqrt(mbsq),
     .  '                *'
       write(6,*) '****************************************************'
       first=.false.
      endif
      
C---Initialize to zero
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      enddo
      enddo

C---Fill spinor products
      call spinoru(6,p,za,zb)

c ensure that we have a hard process
      if (
     .      (s(5,6) .lt. four*mbsq) 
     . .or. (s(1,5)*s(2,5)/s(1,2) .lt. mbsq) 
     . .or. (s(1,6)*s(2,6)/s(1,2) .lt. mbsq) ) return
      
C--calculate matrix element squared
      qqb=msqwbb(1,2,5,6)
      qbq=msqwbb(2,1,5,6)

      do j=-(nf-1),(nf-1)
      do k=-(nf-1),(nf-1)
      if     ((j .gt. 0) .and. (k .lt. 0)) then
               msq(j,k)=Vsq(j,k)*qqb
      elseif ((j .lt. 0) .and. (k .gt. 0)) then
               msq(j,k)=Vsq(j,k)*qbq
      endif
      enddo
      enddo

      return
   
   99 format(a26,f6.3,a21)   
      
      end







