
      subroutine qqb_zaa_frag(p,msq)
****************************************************************
*     Matrix element Fragmentation contribution for            *
*     f(-p1)+fbar(-p2) -> Z(l(p3)+lb(p4))+gam(p5) + gam(p6)    *
*     averaged over initial colours and spins                  *
****************************************************************
      implicit none
      include 'constants.f'
      include 'ewcharge.f'
      include 'frag.f'
      integer j,k,i
      double precision msq(-nf:nf,-nf:nf),p(mxpart,4)
      double precision msq0(-nf:nf,-nf:nf)
      double precision D(0:5),fsq
      common/D/D
!$omp threadprivate(/D/)

c-----initialize MSQ
      do j=-nf,nf
      do k=-nf,nf
         msq(j,k) =0d0
         msq0(j,k)=0d0
      enddo
      enddo
c-----set fragmentation scale
      fsq=frag_scale**2
c-----generate array D(j) corresponding to MCFM notation 0=gluon 1=down 2=up ....
      do i=0,5
         D(i)=0d0
         if     (fragset .eq. 'BFGset_I') then
            call get_frag(z_frag,fsq,1,i,D(i))   
         elseif (fragset .eq. 'BFGsetII') then  
            call get_frag(z_frag,fsq,2,i,D(i))   
         elseif (fragset .eq. 'GdRG__LO') then 
            call GGdR_frag(z_frag,i,D(i),0)
         else
            write(6,*) 'Unrecognized fragmentation set name: ',fragset
            stop        
         endif
      enddo
c-----call qqb_zaj LO matelem
      call qqb_zaj(p,msq0)
c-----fill msq_frag
      do j=-nf,nf
      do k=-nf,nf
c-----qqb      
      if ((j .gt. 0) .and. (k .lt. 0)) then
          if (j .eq. -k) msq(j,k)=Q(j)**2*msq0(j,k)*D(0)
c-----qbq      
      elseif ((j .lt. 0) .and. (k .gt. 0)) then
          if (j .eq. -k) msq(j,k)=Q(k)**2*msq0(j,k)*D(0)
c-----qg
      elseif ((j .gt. 0) .and. (k .eq. 0)) then
            msq(j,k)=Q(j)**2*msq0(j,k)*D(abs(j))
c-----qbg      
      elseif ((j .lt. 0) .and. (k .eq. 0)) then
            msq(j,k)=Q(j)**2*msq0(j,k)*D(abs(j))
c-----gq
      elseif ((j .eq. 0) .and. (k .gt. 0)) then
            msq(j,k)=Q(k)**2*msq0(j,k)*D(abs(k))
c-----gqb      
      elseif ((j .eq. 0) .and. (k .lt. 0)) then
            msq(j,k)=Q(k)**2*msq0(j,k)*D(abs(k))
      endif

      enddo
      enddo
c-----done
      return
      end


      
