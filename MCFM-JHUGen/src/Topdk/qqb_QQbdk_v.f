      subroutine qqb_QQbdk_v(p,msq)
      implicit none
      include 'types.f'

************************************************************************
*     Author: R.K. Ellis                                               *
*     Ocotber 2011.                                                    *
*     Virtual matrix elements squared for the process                  *
*                                                                      *
*     q(-p1) +qbar(-p2)=bbar(p6)+e-(p7)+nubar(p8)+nu(p3)+e+(p4)+b(p5)  *
************************************************************************
c--- Computed using the helicity amplitudes from:Badger,Sattler,Yundin
c--- arXiv:1101.5947[hep-ph]
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'masses.f'
      include 'zprods_com.f'
      include 'sprods_com.f'
      include 'zabprods_decl.f'
      include 'scheme.f'
      include 'etadef.f'
      include 'qdef.f'
      include 'epinv.f'
      include 'scale.f'
      include 'kprocess.f'
      include 'plabel.f'
      integer:: j,k,nu,j1,j3
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4),q(mxpart,4)
      real(dp):: s34,s35,s45,s67,s68,s78,fac,ren
      real(dp):: resqqb0,resqqb1,resqbq1,resgg0,resgg1,ss,c1,c4
      complex(dp):: prop,amp0(2),amp1(2),loab(2,2),loba(2,2),
     & ab43(2,2),ba43(2,2),ab41(2,2),ba41(2,2),loqed(2,2)
      complex(dp):: BSYA0qqppmp,BSYA0ggpppp,BSYA0ggppmp,
     & A41qqppmp,A41ggpppp,A41ggppmp,A43ggpppp,A43ggppmp
      real(dp):: nlf,nhf
      parameter(nlf=5._dp,nhf=1._dp)
      external BSYA0qqppmp,BSYA0ggpppp,BSYA0ggppmp,
     & A41qqppmp,A41ggpppp,A41ggppmp,A43ggpppp,A43ggppmp
C---statement function
      ss(j,k)=2._dp
     & *(p(j,4)*p(k,4)-p(j,1)*p(k,1)-p(j,2)*p(k,2)-p(j,3)*p(k,3))

      scheme='dred'

      s34=ss(3,4)
      s35=ss(3,5)
      s45=ss(4,5)
      s67=ss(6,7)
      s68=ss(6,8)
      s78=ss(7,8)
C----set all elements to zero
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0._dp
      enddo
      enddo

      prop=cplx2(s34-wmass**2,wmass*wwidth)
     &    *cplx2(s78-wmass**2,wmass*wwidth)
     &    *cplx2(zip,mt*twidth)**2
      fac=V*gwsq**4*gsq**2/abs(prop)**2*s35*s68
C----Include factor for hadronic decays of W
      if (plabel(3) == 'pp') fac=2._dp*xn*fac
      if (plabel(7) == 'pp') fac=2._dp*xn*fac

C-----setup top and anti-top momenta and calculate integral functions
      do nu=1,4
      q(1,nu)=p(3,nu)+p(4,nu)+p(5,nu)
      q(2,nu)=p(1,nu)
      q(3,nu)=p(2,nu)
      q(4,nu)=p(6,nu)+p(7,nu)+p(8,nu)
      enddo
      call integralfill(q)


C-----make top and topb massless wrt e+(4) and e-(7) momentum
      c1=mt**2/(s34+s45)
      c4=mt**2/(s67+s78)

      do nu=1,4
      q(1,nu)=q(1,nu)-c1*p(4,nu)
      q(2,nu)=p(1,nu)
      q(3,nu)=p(2,nu)
      q(4,nu)=q(4,nu)-c4*p(7,nu)
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

      s(1,2)=real(zab(2,q1,2))
      s(1,3)=real(zab(3,q1,3))
      s(2,1)=s(1,2)
      s(3,1)=s(1,3)

      call qqbgen(e1,2,3,e4,za,zb,zab,zba,BSYA0qqppmp,amp0)
      call qqbgen(e1,2,3,e4,za,zb,zab,zba,A41qqppmp,amp1)

      resqqb0=fac*(abs(amp0(1))**2+abs(amp0(2))**2)
      resqqb1=xn*fac
     & *real(conjg(amp0(1))*amp1(1)+conjg(amp0(2))*amp1(2))

      call qqbgen(e1,3,2,e4,za,zb,zab,zba,BSYA0qqppmp,amp0)
      call qqbgen(e1,3,2,e4,za,zb,zab,zba,A41qqppmp,amp1)
      resqbq1=xn*fac
     & *real(conjg(amp0(1))*amp1(1)+conjg(amp0(2))*amp1(2))

      call gluegen(e1,2,3,e4,za,zb,zab,zba,
     & BSYA0ggpppp,BSYA0ggppmp,loab,loba,.true.)
      call gluegen(e1,2,3,e4,za,zb,zab,zba,
     & A41ggpppp,A41ggppmp,ab41,ba41,.true.)
      call gluegen(e1,2,3,e4,za,zb,zab,zba,
     & A43ggpppp,A43ggppmp,ab43,ba43,.false.)


c--- note that filling of colour structures here looks unnatural:
c---    1 <--> loba , 2 <--> loab
c---  but this does correspond to the filling in qqb_QQb.f
      resgg0=0._dp
      resgg1=0._dp
      do j=1,2
      do k=1,2
      loqed(j,k)=loab(j,k)+loba(j,k)

      resgg0=resgg0
     & +xn*fac*(abs(loab(j,k))**2+abs(loba(j,k))**2
     &           -abs(loqed(j,k))**2/xn**2)

      resgg1=resgg1+xn**2*fac*real(
     & +ab41(j,k)*conjg(loab(j,k))+ba41(j,k)*conjg(loba(j,k))
     & +(ab43(j,k)-ab41(j,k)-ba41(j,k))*conjg(loqed(j,k))/xn**2)

      enddo
      enddo

c--- renormalization
      ren=(
     & +2._dp*((nlf/xn-11._dp/6._dp*xn)*epinv+xn/6._dp)
     & +2._dp*nhf/xn*(epinv+log(musq/mt**2))
     & -Cf*(3._dp*(epinv+log(musq/mt**2))+5._dp))
      resqqb1=resqqb1+ren*resqqb0
      resqbq1=resqbq1+ren*resqqb0

      ren=(
     & +2._dp*((nlf/3._dp-11._dp/6._dp*xn)*epinv+xn/6._dp)
     & -(xn**2-1._dp)/2._dp/xn*(3._dp*(epinv+log(musq/mt**2))+5._dp))
      resgg1=resgg1+ren*resgg0

C----set all elements to zero
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0._dp
      enddo
      enddo

C---fill qb-q, gg and q-qb elements
      do j=-nf,nf
      if (j < 0) then
          msq(j,-j)=ason2pi*aveqq*resqbq1
c          msq(j,-j)=aveqq*resqqb0
      elseif (j == 0) then
          msq(j,j)=ason2pi*avegg*resgg1
c          msq(j,j)=avegg*resgg0
      elseif (j > 0) then
          msq(j,-j)=ason2pi*aveqq*resqqb1
c          msq(j,-j)=aveqq*resqqb0
      endif
      enddo

      return
      end
