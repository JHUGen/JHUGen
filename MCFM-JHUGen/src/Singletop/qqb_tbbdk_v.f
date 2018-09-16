      subroutine qqb_tbbdk_v(p,msqv)
      implicit none
      include 'types.f'

************************************************************************
*     Author: R.K. Ellis                                               *
*     February, 2012.                                                  *
*     calculate the virtual matrix elements for the process            *
*                                                                      *
*     u(-p1) +dbar(-p2)=t(nu(p3)+e+(p4)+b(p5))+bbar(p6)                *
*     or                                                               *
*     d(-p1) +ubar(-p2)=t~(e-(p3)+nu~(p4)+bb(p5))+b(p6)                *
*                                                                      *
*     Top (antitop) is kept strictly on-shell                          *
*     although all spin correlations are retained.                     *
*                                                                      *
*     NOTE: this routine is a replacement for qqb_tbb_v.f that         *
*           includes the effect of the b-quark mass.                   *
*                                                                      *
************************************************************************
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'ckm.f'
      include 'nwz.f'
      include 'masses.f'
      include 'scheme.f'
      real(dp):: p(mxpart,4),msqv(-nf:nf,-nf:nf),
     & qqb,qbq,fac,facho
      integer:: j,k,hb,hc,j1
      complex(dp)::  prop,mtop(2,2),manti(2,2),
     & mtot12(2,2),mtot12v(2,2),mtot21(2,2),mtot21v(2,2),
     & m12(2,2),m21(2,2),m12v(2,2),m21v(2,2)

      scheme='dred'
C----set all elements to zero
      do j=-nf,nf
      do k=-nf,nf
      msqv(j,k)=0._dp
      enddo
      enddo

      if (nwz == +1) then

      call schantoponshellv(1,2,p,m12,m12v)
      call schantoponshellv(2,1,p,m21,m21v)
      call tdecay(p,3,4,5,mtop)

      do hb=1,2
      do hc=1,2
      mtot12(hb,hc)=czip
      mtot12v(hb,hc)=czip
      mtot21(hb,hc)=czip
      mtot21v(hb,hc)=czip
      do j1=1,2
      mtot12(hb,hc)=mtot12(hb,hc)+mtop(hb,j1)*m12(j1,hc)
      mtot21(hb,hc)=mtot21(hb,hc)+mtop(hb,j1)*m21(j1,hc)
      mtot12v(hb,hc)=mtot12v(hb,hc)+mtop(hb,j1)*m12v(j1,hc)
      mtot21v(hb,hc)=mtot21v(hb,hc)+mtop(hb,j1)*m21v(j1,hc)
      enddo
      enddo
      enddo

      elseif (nwz == -1) then

      call schanatoponshellv(1,2,p,m12,m12v)
      call schanatoponshellv(2,1,p,m21,m21v)
      call adecay(p,3,4,5,manti)

      do hb=1,2
      do hc=1,2
      mtot12(hb,hc)=czip
      mtot12v(hb,hc)=czip
      mtot21(hb,hc)=czip
      mtot21v(hb,hc)=czip
      do j1=1,2
      mtot12(hb,hc)=mtot12(hb,hc)+m12(hb,j1)*manti(j1,hc)
      mtot21(hb,hc)=mtot21(hb,hc)+m21(hb,j1)*manti(j1,hc)
      mtot12v(hb,hc)=mtot12v(hb,hc)+m12v(hb,j1)*manti(j1,hc)
      mtot21v(hb,hc)=mtot21v(hb,hc)+m21v(hb,j1)*manti(j1,hc)
      enddo
      enddo
      enddo

      endif

      prop=cplx2(zip,mt*twidth)
      fac=xn**2*aveqq*gw**8/abs(prop)**2
C  Factor of two for interference effectively included by multiplying only
C  by ason2pi
      facho=fac*ason2pi*cf
      qqb=0._dp
      qbq=0._dp
      do hb=1,2
      do hc=1,2
      qqb=qqb+facho*real(conjg(mtot12(hb,hc))*mtot12v(hb,hc))
      qbq=qbq+facho*real(conjg(mtot21(hb,hc))*mtot21v(hb,hc))
      enddo
      enddo


      do j=-nf,nf
      do k=-nf,nf
      msqv(j,k)=0._dp
      if     ((j > 0) .and. (k < 0)) then
      msqv(j,k)=Vsq(j,k)*qqb
      elseif ((j < 0) .and. (k > 0)) then
      msqv(j,k)=Vsq(j,k)*qbq
      endif
      enddo
      enddo

      return
      end

