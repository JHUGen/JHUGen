      subroutine compare_madgraph(p,mcfmroutine,madroutine)
      implicit none
      include 'types.f'
      
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'nwz.f'
      real(dp):: 
     & p(mxpart,4),msq(-nf:nf,-nf:nf),msqc(-nf:nf,-nf:nf),
     & aemmz,ee,
     & P1(0:3),P2(0:3),P3(0:3),P4(0:3),P5(0:3),P6(0:3),P7(0:3)                 
      integer:: i,j,k
      common/madmom/p1,p2,p3,p4,p5,p6,p7
      external mcfmroutine,madroutine

c--- implement the momentum exchange      
      do i=1,4
        if (i<4) then
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
      aemmz=1._dp/128._dp  
      esq=fourpi*aemmz    
      ee=sqrt(esq)
      gwsq=fourpi*aemmz/xw
      gw=sqrt(gwsq)
      gsq=1._dp
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
c      if ((j<0 .and. k<0) .or. (j>0 .and. k>0)) then
c      if (abs(msq(j,k)/msqc(j,k)-1._dp) >1.e-6_dp)
      if ((msq(j,k)>0._dp) .or. (msqc(j,k)>0._dp) ) then
      write(6,99) j,k,msq(j,k),msqc(j,k),msq(j,k)/msqc(j,k)
      endif
      enddo
      enddo
      pause


   98 format(2i3,g17.6,g17.6)
   99 format(2i3,2e17.9,f21.12)

      
      return
      end
