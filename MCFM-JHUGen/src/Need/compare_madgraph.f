      subroutine compare_madgraph(p,mcfmroutine,madroutine)
      implicit none
      
      include 'constants.f'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'nwz.f'
      double precision 
     . p(mxpart,4),msq(-nf:nf,-nf:nf),msqc(-nf:nf,-nf:nf)
      integer i,j,k
      double precision aemmz,ee
      real*8 P1(0:3),P2(0:3),P3(0:3),P4(0:3),P5(0:3),P6(0:3),P7(0:3)                 
      common/madmom/p1,p2,p3,p4,p5,p6,p7
      external mcfmroutine,madroutine

c--- implement the momentum exchange      
      do i=1,4
        if (i.lt.4) then
          j=i
        else
          j=0
        endif 
        p1(j)=-p(1,i)
        p2(j)=-p(2,i)
        p3(j)=p(3,i)
        p4(j)=p(4,i)
        p5(j)=p(5,i)
        p6(j)=p(6,i)
        p7(j)=p(7,i)
      enddo

c--- make call to MCFM routines
      aemmz=1d0/128d0  
      esq=fourpi*aemmz    
      ee=dsqrt(esq)
      gwsq=fourpi*aemmz/xw
      gw=sqrt(gwsq)
      gsq=1d0
      call ckmfill(nwz)
      

      call initialize

      call mcfmroutine(p,msq)
      call madroutine(p,msqc)


      write(*,*) '        MCFM             MADG               ratio'

      
c      do j=-nf,nf
c      do k=-nf,nf
c      write(6,98) j,k,msq(j,k),msqc(j,k)
c      enddo
c      enddo

c      pause


      do j=-nf,nf
      do k=-nf,nf
c      do j=-4,4
c      do k=-4,4
c      if ((j.lt.0 .and. k.lt.0) .or. (j.gt.0 .and. k.gt.0)) then
c      if (abs(msq(j,k)/msqc(j,k)-1d0) .gt.1d-6)
      if ((msq(j,k).gt.0d0) .or. (msqc(j,k).gt.0d0) ) then
      write(6,99) j,k,msq(j,k),msqc(j,k),msq(j,k)/msqc(j,k)
      endif
      enddo
      enddo
      pause


   98 format(2i3,g17.6,g17.6)
   99 format(2i3,2e17.9,f21.12)

      
      return
      end
