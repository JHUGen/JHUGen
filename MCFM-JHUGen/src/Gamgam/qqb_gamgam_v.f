      subroutine qqb_gamgam_v(p,msq)
      implicit none
************************************************************************
*     Authors: R.K. Ellis and John M. Campbell                         *
*     December, 2010.                                                  *
************************************************************************
*                   .                                                  *
*     Matrix element for gamma + gamma production,                     *
*     averaged over initial colours and spins                          *
*                   .                                                  *
*     q(-p1)+qbar(-p2) --> gamma(p3)+gamma(p4)                         *
*                   .                                                  *
************************************************************************
      include 'constants.f'
      include 'ewcouple.f'
      include 'ewcharge.f'
      include 'qcdcouple.f'
      include 'sprods_com.f'
      include 'scheme.f'
      integer j,k
      double precision qa,aq,gg
      double precision msq(-nf:nf,-nf:nf),p(mxpart,4),fac,facgg,Qsum,
     & qagamgam,statfac,virtgamgam
      parameter(statfac=0.5d0)

      scheme='tH-V'

      fac=8d0*xn*esq**2*ason2pi*statfac

      call dotem(4,p,s)
      qa=+qagamgam(3,4,2)*fac*aveqq
      aq=qa

c--- initialize gg 2-loop matrix elements
      Qsum=+Q(1)**2+Q(2)**2+Q(3)**2+Q(4)**2+Q(5)**2 
      facgg=4d0*esq*gsq/(16d0*pisq)*Qsum
      gg=avegg*V*facgg**2*statfac*virtgamgam(s(1,2),s(1,3),s(2,3))

      do j=-nf,nf
      do k=-nf,nf
c--set msq=0 to initalize
      msq(j,k)=0d0
C--qa      
      if ((j .gt. 0) .and. (k .lt. 0)) then
          if (j .eq. -k) msq(j,k)=Q(j)**4*qa
C--aq      
      elseif ((j .lt. 0) .and. (k .gt. 0)) then
          if (j .eq. -k) msq(j,k)=Q(k)**4*aq
C--gg      
      elseif ((j .eq. 0) .and. (k .eq. 0)) then
          msq(j,k)=gg
      endif

      enddo
      enddo

      return
      end

      double precision function qagamgam(i1,i2,i3)
      implicit none
c----Matrix element for gamma + gamma production
c----in order alpha_s
c---
c---  0 -> gamma(p1) + gamma(p2) + q(p3) + qb(p4)
c---
      include 'constants.f'
      include 'sprods_com.f'
      include 'epinv.f'
      include 'scale.f'
      include 'scheme.f'
      integer i1,i2,i3
      double complex lnrat,l12,l13,l23
      double precision s12,s13,s23,T0,deltar
      s12=s(i1,i2)
      s13=s(i1,i3)
      s23=s(i2,i3)
      T0=s13/s23+s23/s13
      l12=lnrat(-s12,musq)
      l13=lnrat(-s13,musq)
      l23=lnrat(-s23,musq)

      if     (scheme .eq. 'tH-V') then
        deltar=0d0
      elseif (scheme .eq. 'dred') then
        deltar=1d0
      else
        write(6,*) 'Invalid scheme in qqb_gamgam_v.f'
        stop
      endif
      
CId,anscdr=cf*((-2*epinv**2-epinv*(3-2*[ln(-s12)])-7-[ln(-s12)]**2)*T0
C    +[ln(-s13)]*(s23-2*s12)/s13
C    +[ln(-s23)]*(s13-2*s12)/s23
C    +(s12**2+s23**2)/s13/s23*(([ln(-s12)]-[ln(-s23)])**2+pisq)
c    +(s12**2+s13**2)/s13/s23*(([ln(-s12)]-[ln(-s13)])**2+pisq)
C    -4*[ln(-s12)]);


      qagamgam=cf*((-2d0*epinv**2-epinv*(3d0-2d0*dble(l12))
     . -7d0+deltar-dble(l12**2))*T0
     . +dble(l13)*(s23-2d0*s12)/s13
     . +dble(l23)*(s13-2d0*s12)/s23
     . +(s12**2+s23**2)/s13/s23*(dble((l12-l23)**2)+pisq)
     . +(s12**2+s13**2)/s13/s23*(dble((l12-l13)**2)+pisq)
     . -4d0*dble(l12))

      return 
      end


