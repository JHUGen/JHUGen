      subroutine qqb_gamgam(p,msq)
      implicit none
      include 'types.f'
C=====
C-----Matrix element for f(-p1)+f(-p2)->gamma(p3)+gamma(p4)
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'qcdcouple.f'
      include 'ewcharge.f'
      include 'ewcouple.f'
      include 'sprods_com.f'
      include 'noglue.f'
      integer:: j,k
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4),fac,
     &  qa,aq,gg,facgg,msqgggaga,Qsum
      real(dp),parameter::statfac=0.5_dp

      do j=-nf,nf
      do k=-nf,nf
c--set msq=0 to initalize
      msq(j,k)=0._dp
      enddo
      enddo

      call dotem(3,p,s)

      fac=8._dp*xn*esq**2*statfac

      Qsum=+Q(1)**2+Q(2)**2+Q(3)**2+Q(4)**2+Q(5)**2 
      facgg=4._dp*esq*gsq/(16._dp*pisq)*Qsum

      gg=zip
! JC removing gg contribution
!      if (omitgg) then
!        gg=0._dp
!      else
!        gg=avegg*V*facgg**2*msqgggaga(s(1,2),s(1,3),s(2,3))*statfac
!      endif
      
      qa=fac*aveqq*(s(1,3)/s(2,3)+s(2,3)/s(1,3))
      aq=qa
      
      do j=-nf,nf
      k=-j
C--qa      
      if (j > 0) then
        msq(j,k)=Q(j)**4*qa
C--aq      
      elseif (k > 0) then
        msq(j,k)=Q(k)**4*aq
C--gg
      elseif (j == 0) then
        msq(j,k)=gg
      endif
 
      enddo

      return
      end
