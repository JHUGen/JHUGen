      subroutine qqb_QQbdkBSY(p,msq)
      implicit none
      include 'types.f'

************************************************************************
*     Author: R.K. Ellis                                               *
*     October, 2011.                                                   *
*     calculate the element squared and subtraction terms              *
*     for the process                                                  *
*                                                                      *
*     q(-p1) +qbar(-p2)=bbar(p6)+e-(p7)+nubar(p8)+nu(p3)+e+(p4)+b(p5)  *
*                                                                      *
*     Top is kept strictly on-shell although all spin correlations     *
*     are retained.                                                    *
*                                                                      *
************************************************************************
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'masses.f'
      include 'kprocess.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      include 'zabprods_decl.f'
      include 'msq_cs.f'
      include 'etadef.f'
      include 'qdef.f'
      integer:: j,k,nu,cs,j1,j3
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4),q(mxpart,4)
      real(dp):: fac,qqb,c1,c4,ss,s34,s35,s45,s67,s68,s78
      complex(dp)::  amp(2),prop,loab(2,2),loba(2,2),loqed(2,2)
      complex(dp)::  BSYA0qqppmp,BSYA0ggpppp,BSYA0ggppmp
      external BSYA0qqppmp,BSYA0ggpppp,BSYA0ggppmp
      ss(j,k)=2d0
     & *(p(j,4)*p(k,4)-p(j,1)*p(k,1)-p(j,2)*p(k,2)-p(j,3)*p(k,3))

      s34=ss(3,4)
      s35=ss(3,5)
      s45=ss(4,5)
      s67=ss(6,7)
      s68=ss(6,8)
      s78=ss(7,8)

C----set all elements to zero
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      do cs=0,2
      msq_cs(cs,j,k)=0d0
      enddo
      enddo
      enddo

      prop=cplx2(s34-wmass**2,wmass*wwidth)
     &    *cplx2(s78-wmass**2,wmass*wwidth)
     &    *cplx2(zip,mt*twidth)**2
      fac=V*gwsq**4*gsq**2/abs(prop)**2*s35*s68
C--include factor for hadronic decays
      if ((kcase==ktt_bbh) .or. (kcase==ktt_hdk)) fac=2d0*xn*fac

C-----make top and topb massless wrt e+(4) and e-(7) momentum
      c1=mt**2/(s34+s45)
      c4=mt**2/(s67+s78)

      do nu=1,4
      q(1,nu)=p(3,nu)+p(4,nu)+p(5,nu)-c1*p(4,nu)
      q(2,nu)=p(1,nu)
      q(3,nu)=p(2,nu)
      q(4,nu)=p(6,nu)+p(7,nu)+p(8,nu)-c4*p(7,nu)
      q(e1,nu)=p(4,nu)
      q(e4,nu)=p(7,nu)
      enddo

      call spinoru(6,q,za,zb)
      call spinorextend(za,zb)
      do j1=1,8
      do j3=1,8
      zab(j1,q1,j3)=+za(j1,1)*zb(1,j3)+c1*za(j1,e1)*zb(e1,j3)
      zba(j3,q1,j1)=zab(j1,q1,j3)
      enddo
      enddo

C---currently s(1,2) and s(1,3) are given as s(1f,2) and s(1f,3)
C---but we want the full s(1,2) and s(1,3)
C---hence restore them

      s(1,2)=real(zab(2,q1,2))
      s(1,3)=real(zab(3,q1,3))
      s(2,1)=s(1,2)
      s(3,1)=s(1,3)

      call qqbgen(e1,2,3,e4,za,zb,zab,zba,
     & BSYA0qqppmp,amp)
      qqb=fac*aveqq*(abs(amp(1))**2+abs(amp(2))**2)

      call gluegen(e1,2,3,e4,za,zb,zab,zba,
     & BSYA0ggpppp,BSYA0ggppmp,loab,loba,.true.)

      do j=1,2
      do k=1,2
      enddo
      enddo
      do j=1,2
      do k=1,2
      enddo
      enddo

c--- note that filling of colour structures here looks unnatural:
c---    1 <--> loba , 2 <--> loab
c---  but this does correspond to the filling in qqb_QQb.f
      do j=1,2
      do k=1,2
      loqed(j,k)=loab(j,k)+loba(j,k)
      msq_cs(1,0,0)=msq_cs(1,0,0)+fac*avegg*xn*abs(loba(j,k))**2
      msq_cs(2,0,0)=msq_cs(2,0,0)+fac*avegg*xn*abs(loab(j,k))**2
      msq_cs(0,0,0)=msq_cs(0,0,0)-fac*avegg/xn*abs(loqed(j,k))**2
      enddo
      enddo

C---fill qb-q, gg and q-qb elements
      do j=-nf,nf
      if ((j < 0) .or. (j > 0)) then
          msq(j,-j)=qqb
C Division of quark into color structures is arbitrary
          msq_cs(1,j,-j)=qqb/3d0
          msq_cs(2,j,-j)=qqb/3d0
          msq_cs(0,j,-j)=qqb/3d0
      elseif (j == 0) then
          msq(0,0)=msq_cs(1,0,0)+msq_cs(2,0,0)+msq_cs(0,0,0)
      endif
      enddo

      return
      end
