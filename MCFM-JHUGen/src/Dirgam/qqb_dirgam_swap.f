      subroutine qqb_dirgam_swap(p,msq)
C=====
C-----Matrix element for f(-p1)+f(-p2)->g(p3)+gamma(p4)
!--- Need for gamgam dipoles
      implicit none
      include 'constants.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'ewcharge.f'
      include 'sprods_com.f'
      integer j,k
      double precision msq(-nf:nf,-nf:nf),p(mxpart,4),fac,
     .  qa,aq,qg,gq,ag,ga
    

      call dotem(4,p,s)
      fac=4d0*V*gsq*esq

      qa=fac*aveqq*(s(1,4)/s(2,4)+s(2,4)/s(1,4))
      aq=qa
      qg=-fac*aveqg*(s(1,4)/s(1,2)+s(1,2)/s(1,4))
      ag=qg
      gq=-fac*aveqg*(s(1,2)/s(2,4)+s(2,4)/s(1,2))
      ga=gq

      do j=-nf,nf
      do k=-nf,nf
c--set msq=0 to initalize
      msq(j,k)=0d0
C--qa      
      if ((j .gt. 0) .and. (k .lt. 0)) then
          if (j .eq. -k) msq(j,k)=Q(j)**2*qa
C--aq      
      elseif ((j .lt. 0) .and. (k .gt. 0)) then
          if (j .eq. -k) msq(j,k)=Q(k)**2*aq
C--qg
      elseif ((j .gt. 0) .and. (k .eq. 0)) then
            msq(j,k)=Q(j)**2*qg
C--ag      
      elseif ((j .lt. 0) .and. (k .eq. 0)) then
            msq(j,k)=Q(j)**2*ag
C--gq
      elseif ((j .eq. 0) .and. (k .gt. 0)) then
            msq(j,k)=Q(k)**2*gq
C--ga      
      elseif ((j .eq. 0) .and. (k .lt. 0)) then
            msq(j,k)=Q(k)**2*ga
      endif

      enddo
      enddo


!--- restore p4 and p3
!      do j=1,4
!         ptemp(j)=p(3,j)
!         p(3,j)=p(4,j)
!         p(4,j)=ptemp(j)
!      enddo


     

      return
      end


      
