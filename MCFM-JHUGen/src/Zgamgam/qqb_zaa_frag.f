
      subroutine qqb_zaa_frag(p,msq)
      implicit none
      include 'types.f'
****************************************************************
*     Matrix element Fragmentation contribution for            *
*     f(-p1)+fbar(-p2) -> Z(l(p3)+lb(p4))+gam(p5) + gam(p6)    *
*     averaged over initial colours and spins                  *
****************************************************************
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ewcharge.f'
      include 'frag.f'
      integer:: j,k,i
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4)
      real(dp):: msq0(-nf:nf,-nf:nf)
      real(dp):: D(0:5),fsq
      common/D/D
!$omp threadprivate(/D/)

c-----initialize MSQ
      do j=-nf,nf
      do k=-nf,nf
         msq(j,k) =0._dp
         msq0(j,k)=0._dp
      enddo
      enddo
c-----set fragmentation scale
      fsq=frag_scale**2
c-----generate array D(j) corresponding to MCFM notation 0=gluon 1=down 2=up ....
      do i=0,5
         D(i)=0._dp
         if     (fragset == 'BFGset_I') then
            call get_frag(z_frag,fsq,1,i,D(i))   
         elseif (fragset == 'BFGsetII') then  
            call get_frag(z_frag,fsq,2,i,D(i))   
         elseif (fragset == 'GdRG__LO') then 
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
      if ((j > 0) .and. (k < 0)) then
          if (j == -k) msq(j,k)=Q(j)**2*msq0(j,k)*D(0)
c-----qbq      
      elseif ((j < 0) .and. (k > 0)) then
          if (j == -k) msq(j,k)=Q(k)**2*msq0(j,k)*D(0)
c-----qg
      elseif ((j > 0) .and. (k == 0)) then
            msq(j,k)=Q(j)**2*msq0(j,k)*D(abs(j))
c-----qbg      
      elseif ((j < 0) .and. (k == 0)) then
            msq(j,k)=Q(j)**2*msq0(j,k)*D(abs(j))
c-----gq
      elseif ((j == 0) .and. (k > 0)) then
            msq(j,k)=Q(k)**2*msq0(j,k)*D(abs(k))
c-----gqb      
      elseif ((j == 0) .and. (k < 0)) then
            msq(j,k)=Q(k)**2*msq0(j,k)*D(abs(k))
      endif

      enddo
      enddo
c-----done
      return
      end


      
