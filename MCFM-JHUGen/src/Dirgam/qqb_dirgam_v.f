      subroutine qqb_dirgam_v(p,msq)
      implicit none
      include 'types.f'
      
c----Matrix element for gamma + jet production
c----in order alpha_s^2
C----averaged over initial colours and spins
c     q(-p1)+qbar(-p2)-->gamma(p3)+f(p4)
c---
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ewcouple.f'
      include 'ewcharge.f'
      include 'qcdcouple.f'
      include 'sprods_com.f'
      include 'scheme.f'
      real(dp):: qa,aq,qg,ag,gq,ga
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4),fac,qaggam
      integer:: j,k
      fac=4._dp*V*gsq*esq*ason2pi

      scheme='tH-V'

      call dotem(4,p,s)
      qa=+qaggam(1,2,4)*fac*aveqq
c      aq=+qaggam(2,1,4)*fac*aveqq
      aq=qa
      qg=-qaggam(1,4,2)*fac*aveqg
c      ag=-qaggam(4,1,2)*fac*aveqg
      ag=qg

      gq=-qaggam(2,4,1)*fac*aveqg
c      ga=-qaggam(4,2,1)*fac*aveqg
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
            msq(j,k)=Q(-k)**2*ga
      endif

      enddo
      enddo

      return
      end

      function qaggam(i1,i2,i3)
      implicit none
      include 'types.f'
      real(dp):: qaggam
      
c----Matrix element for gamma + jet production
c----in order alpha_s^2
C----averaged over initial colours and spins
c     q(-p1)+qbar(-p2)-->g(p3)+gamma(p4)
c---
C---- Matrix element derived from ERT
C---- %\cite{Ellis:1980wv}
C---- \bibitem{Ellis:1980wv}
C---- R.~K.~Ellis, D.~A.~Ross and A.~E.~Terrano,
C---- %``The Perturbative Calculation Of Jet Structure In E+ E- Annihilation,''
C---- Nucl.\ Phys.\ B {\bf 178}, 421 (1981).
C---- %%CITATION = NUPHA,B178,421;%%

C---- Matrix element also checked from Eq. (2.28)
C---- %\cite{Aurenche:1986ff}
C---- \bibitem{Aurenche:1986ff}
C---- P.~Aurenche, R.~Baier, A.~Douiri, M.~Fontannaz and D.~Schiff,
C---- %``Scheme Invariant Higher Order QCD Predictions For Large P(T)
C---- Photoproduction Reactions,''
C---- Nucl.\ Phys.\ B {\bf 286}, 553 (1987).
C---- %%CITATION = NUPHA,B286,553;%%

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'sprods_com.f'
      include 'epinv.f'
      include 'scale.f'
      include 'scheme.f'
      integer:: i1,i2,i3
      complex(dp):: lnrat,l12,l13,l23
      real(dp):: s12,s13,s23,T0,deltar

      s12=s(i1,i2)
      s13=s(i1,i3)
      s23=s(i2,i3)
      T0=s13/s23+s23/s13
      l12=lnrat(-s12,musq)
      l13=lnrat(-s13,musq)
      l23=lnrat(-s23,musq)

      if     (scheme == 'tH-V') then
        deltar=0._dp
      elseif (scheme == 'dred') then
        deltar=1._dp
      else
        write(6,*) 'Invalid scheme in qqb_gamgam_v.f'
        stop
      endif
      
      qaggam=-(2._dp*xn-1._dp/xn)*(epinv**2-pisq/six)*T0
     & -epinv*(3._dp*cf-2._dp*cf*real(l12)+xn*real(l12-l13-l23)
     &  +11._dp*xn/6._dp-nf/3._dp)*T0
     & +((xn-cf)*pisq/3._dp-7._dp*cf+deltar*(cf+xn/6._dp)-cf*real(l12**2)
     & +xn*(0.5_dp*real(l12**2)-real(l13*l23)))*T0
     & +real(l13)*(cf*(s23-2._dp*s12)/s13-xn)
     & +real(l23)*(cf*(s13-2._dp*s12)/s23-xn)
     & -((s12**2+s23**2)/s13/s23*(real((l12-l23)**2)+pisq)
     &  +(s12**2+s13**2)/s13/s23*(real((l12-l13)**2)+pisq)
     &  -4._dp*real(l12))/(2._dp*xn)
      return 
      end



c      write(6,*) 'original,T0',T0
c      write(6,*) 'i1,i2,i3,qaggam',i1,i2,i3,qaggam
c      if ((i1==1) .and. (i2==2) .and.(i3==4)) then
c      ss=s(i1,i2)
c      tt=s(i1,i3)
c      uu=s(i2,i3)
c      lnss=log(+ss/musq)
c      lntt=log(-tt/musq)
c      lnuu=log(-uu/musq)
c      T0=uu/tt+tt/uu
     
c      write(6,*) 'ss,tt,uu',ss,tt,uu
c      write(6,*) 'lnss,l12',lnss,l12
c      write(6,*) 'lntt,l13',lntt,l13
c      write(6,*) 'lnuu,l23',lnuu,l23
c      write(6,*) 'T0',T0
c      check= -(2._dp*xn-1._dp/xn)*epinv**2*T0
c     & -epinv*(3._dp*cf+11._dp/6._dp*xn+lnss/xn-xn*(lntt+lnuu)-nf/3._dp)*T0
c     & -7._dp/2._dp*tt*uu**(-1)*xn + 7._dp/2._dp*tt*uu**(-1)*xn**(-1)
c     & -7._dp/2._dp*tt**(-1)*uu*xn + 7._dp/2._dp*tt**(-1)*uu*xn**(-1)
c     .+lnss*(2._dp*xn**(-1))
c     .+lnss**2
c     .*(-tt*uu**(-1)*xn**(-1)-tt**(-1)*uu*xn**(-1)-2._dp*xn**(-1))
c     .+lntt
c     .*(3._dp/2._dp*tt**(-1)*uu*xn-3._dp/2._dp*tt**(-1)*uu*xn**(-1)-xn**(-1))
c     .+lntt*lnss
c     .*(2._dp*tt*uu**(-1)*xn**(-1)+tt**(-1)*uu*xn**(-1)+2._dp*xn**(-1))
c     .+lntt*lnuu
c     .*(-tt*uu**(-1)*xn-tt**(-1)*uu*xn)
c     .+lntt**2
c     .*(-tt*uu**(-1)*xn**(-1)-1._dp/2._dp*tt**(-1)*uu*xn**(-1)-xn**(-1))
c     .+lnuu
c     .*(3._dp/2._dp*tt*uu**(-1)*xn-3._dp/2._dp*tt*uu**(-1)*xn**(-1)-xn**(-1))
c     .+lnuu*lnss
c     .*(tt*uu**(-1)*xn**(-1)+2._dp*tt**(-1)*uu*xn**(-1)+2._dp*xn**(-1))
c     .+lnuu**2
c     .*(-1._dp/2._dp*tt*uu**(-1)*xn**(-1)-tt**(-1)*uu*xn**(-1)-xn**(-1))
c     .+pi**2
c     .*(1._dp/6._dp*tt*uu**(-1)*xn-1._dp/3._dp*tt*uu**(-1)*xn**(-1)
c     .+1._dp/6._dp*tt**(-1)*uu*xn-1._dp/3._dp*tt**(-1)*uu*xn**(-1))
c      write(6,*) 'check,124',check
c      elseif ((i1==1) .and. (i2==4) .and.(i3==2)) then
c      ss=s(i1,i3)
c      tt=s(i1,i2)
c      uu=s(i2,i3)

c      T0=-(uu/ss+ss/uu)
c      lnss=log(ss/musq)
c      lntt=log(-tt/musq)
c      lnuu=log(-uu/musq)
c      write(6,*) 'ss,tt,uu',ss,tt,uu
c      write(6,*) 'lnss,l13',lnss,l13
c      write(6,*) 'lntt,l12',lntt,l12
c      write(6,*) 'lnuu,l23',lnuu,l23
c      write(6,*) 'T0',T0

c      check=
c     &  -(11._dp*xn-2._dp*nf)/6._dp*epinv*T0
c     & +(-(2._dp*xn-1._dp/xn)*epinv**2
c     & +epinv*(-lntt/xn+xn*lnss+xn*lnuu-3._dp*cf))*T0
c     & +7._dp/2._dp*ss*uu**(-1)*xn-7._dp/2._dp*ss*uu**(-1)*xn**(-1)
c     & +7._dp/2._dp*ss**(-1)*uu*xn-7._dp/2._dp*ss**(-1)*uu*xn**(-1)
c     & +lntt
c     & *(-2._dp*xn**(-1))
c     & +lntt*lnuu
c     & *(-ss*uu**(-1)*xn**(-1)-2*ss**(-1)*uu*xn**(-1)-2*xn**(-1))
c     & +lntt**2
c     & *(ss*uu**(-1)*xn**(-1)+ss**(-1)*uu*xn**(-1)+2*xn**(-1))
c     & +lnuu
c     & *(-3._dp/2._dp*ss*uu**(-1)*xn+3._dp/2._dp*ss*uu**(-1)*xn**(-1)+xn**(-1))
c     & +lnuu**2
c     & *(1._dp/2._dp*ss*uu**(-1)*xn**(-1)+ss**(-1)*uu*xn**(-1)+xn**(-1))
c     & +lnss
c     & *(-3._dp/2._dp*ss**(-1)*uu*xn+3._dp/2._dp*ss**(-1)*uu*xn**(-1)+xn**(-1))
c     & +lnss*lntt
c     & *(-2*ss*uu**(-1)*xn**(-1)-ss**(-1)*uu*xn**(-1)-2*xn**(-1))
c     & +lnss*lnuu
c     & *(ss*uu**(-1)*xn+ss**(-1)*uu*xn)
c     & +lnss**2
c     & *(ss*uu**(-1)*xn**(-1)+1._dp/2._dp*ss**(-1)*uu*xn**(-1)+xn**(-1))
c     & +pi**2
c     & *(-1._dp/6._dp*ss*uu**(-1)*xn+1._dp/3._dp*ss*uu**(-1)*xn**(-1)
c     & +1._dp/3._dp*uu*ss**(-1)*xn**(-1)-1._dp/6._dp*ss**(-1)*uu*xn
c     & +xn**(-1)*(ss-tt)/2._dp/ss)
c      write(6,*) 'check,142',check

c      endif

c      pause


