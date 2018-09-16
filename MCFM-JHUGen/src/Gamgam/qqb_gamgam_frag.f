      subroutine qqb_gamgam_frag(p,msq)
C=====
C-----Matrix element for f(-p1)+f(-p2)->gamma(p3)+g(p4)
      implicit none
      include 'constants.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'ewcharge.f'
      include 'sprods_com.f'
      include 'frag.f'
      
      integer j,k,i
      double precision msq(-nf:nf,-nf:nf),p(mxpart,4),fac,
     .  qa,aq,qg,gq,ag,ga
      double precision D(0:5),fsq
      common/D/D
!$omp threadprivate(/D/)
      
      fsq=frag_scale**2
c---- Generate array D(j) corresponding to MCFM notation 0=gluon 1=down 2=up ....
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



      call dotem(3,p,s)
      fac=4d0*V*gsq*esq

      qa=fac*aveqq*(s(1,3)/s(2,3)+s(2,3)/s(1,3))
      aq=qa
      qg=-fac*aveqg*(s(1,3)/s(1,2)+s(1,2)/s(1,3))
      ag=qg
      gq=-fac*aveqg*(s(1,2)/s(2,3)+s(2,3)/s(1,2))
      ga=gq

      do j=-nf,nf
      do k=-nf,nf
c--set msq=0 to initalize
      msq(j,k)=0d0
C--qa      
      if ((j .gt. 0) .and. (k .lt. 0)) then
          if (j .eq. -k) msq(j,k)=Q(j)**2*qa*D(0)
C--aq      
      elseif ((j .lt. 0) .and. (k .gt. 0)) then
          if (j .eq. -k) msq(j,k)=Q(k)**2*aq*D(0)
C--qg
      elseif ((j .gt. 0) .and. (k .eq. 0)) then
            msq(j,k)=Q(j)**2*qg*D(abs(j))
C--ag      
      elseif ((j .lt. 0) .and. (k .eq. 0)) then
            msq(j,k)=Q(j)**2*ag*D(abs(j))
C--gq
      elseif ((j .eq. 0) .and. (k .gt. 0)) then
            msq(j,k)=Q(k)**2*gq*D(abs(k))
C--ga      
      elseif ((j .eq. 0) .and. (k .lt. 0)) then
            msq(j,k)=Q(k)**2*ga*D(abs(k))
      endif

      enddo
      enddo


      return
      end


      
