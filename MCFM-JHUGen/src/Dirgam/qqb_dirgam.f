      subroutine qqb_dirgam(p,msq)
      implicit none
      include 'types.f'
C=====
C-----Matrix element for f(-p1)+f(-p2)->gamma(p3)+g(p4)
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'ewcharge.f'
      include 'sprods_com.f'
      integer:: j,k
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4),fac,
     &  qa,aq,qg,gq,ag,ga

      call dotem(3,p,s)
      fac=4._dp*V*gsq*esq

      qa=fac*aveqq*(s(1,3)/s(2,3)+s(2,3)/s(1,3))
      aq=qa
      qg=-fac*aveqg*(s(1,3)/s(1,2)+s(1,2)/s(1,3))
      ag=qg
      gq=-fac*aveqg*(s(1,2)/s(2,3)+s(2,3)/s(1,2))
      ga=gq

      do j=-nf,nf
      do k=-nf,nf
c--set msq=0 to initalize
      msq(j,k)=0._dp
C--qa      
      if ((j > 0) .and. (k < 0)) then
          if (j == -k) msq(j,k)=Q(j)**2*qa
C--aq      
      elseif ((j < 0) .and. (k > 0)) then
          if (j == -k) msq(j,k)=Q(k)**2*aq
C--qg
      elseif ((j > 0) .and. (k == 0)) then
            msq(j,k)=Q(j)**2*qg
C--ag      
      elseif ((j < 0) .and. (k == 0)) then
            msq(j,k)=Q(j)**2*ag
C--gq
      elseif ((j == 0) .and. (k > 0)) then
            msq(j,k)=Q(k)**2*gq
C--ga      
      elseif ((j == 0) .and. (k < 0)) then
            msq(j,k)=Q(k)**2*ga
      endif

      enddo
      enddo


      return
      end


      
