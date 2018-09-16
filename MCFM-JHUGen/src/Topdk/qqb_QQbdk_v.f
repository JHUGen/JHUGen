      subroutine qqb_QQbdk_v(p,msq)
      implicit none
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
      include 'process.f'
      include 'plabel.f'
      integer j,k,nu,j1,j3
      double precision msq(-nf:nf,-nf:nf),p(mxpart,4),q(mxpart,4)
      double precision s34,s35,s45,s67,s68,s78,fac,ren
      double precision resqqb0,resqqb1,resqbq1,resgg0,resgg1,ss,c1,c4
      double complex prop,amp0(2),amp1(2),loab(2,2),loba(2,2),
     & ab43(2,2),ba43(2,2),ab41(2,2),ba41(2,2),loqed(2,2)
      double complex BSYA0qqppmp,BSYA0ggpppp,BSYA0ggppmp,
     & A41qqppmp,A41ggpppp,A41ggppmp,A43ggpppp,A43ggppmp
      double precision nlf,nhf
      parameter(nlf=5d0,nhf=1d0)
      external BSYA0qqppmp,BSYA0ggpppp,BSYA0ggppmp,
     & A41qqppmp,A41ggpppp,A41ggppmp,A43ggpppp,A43ggppmp
C---statement function
      ss(j,k)=2d0
     . *(p(j,4)*p(k,4)-p(j,1)*p(k,1)-p(j,2)*p(k,2)-p(j,3)*p(k,3))

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
      msq(j,k)=0d0
      enddo
      enddo

      prop=dcmplx(s34-wmass**2,wmass*wwidth)
     .    *dcmplx(s78-wmass**2,wmass*wwidth)
     .    *dcmplx(zip,mt*twidth)**2
      fac=V*gwsq**4*gsq**2/abs(prop)**2*s35*s68
C----Include factor for hadronic decays of W
      if (plabel(3) .eq. 'pp') fac=2d0*xn*fac
      if (plabel(7) .eq. 'pp') fac=2d0*xn*fac

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

      s(1,2)=dble(zab(2,q1,2))
      s(1,3)=dble(zab(3,q1,3))
      s(2,1)=s(1,2)
      s(3,1)=s(1,3)

      call qqbgen(e1,2,3,e4,za,zb,zab,zba,BSYA0qqppmp,amp0)
      call qqbgen(e1,2,3,e4,za,zb,zab,zba,A41qqppmp,amp1)

      resqqb0=fac*(cdabs(amp0(1))**2+cdabs(amp0(2))**2)
      resqqb1=xn*fac
     & *dble(dconjg(amp0(1))*amp1(1)+dconjg(amp0(2))*amp1(2))

      call qqbgen(e1,3,2,e4,za,zb,zab,zba,BSYA0qqppmp,amp0)
      call qqbgen(e1,3,2,e4,za,zb,zab,zba,A41qqppmp,amp1)
      resqbq1=xn*fac
     & *dble(dconjg(amp0(1))*amp1(1)+dconjg(amp0(2))*amp1(2))

      call gluegen(e1,2,3,e4,za,zb,zab,zba,
     & BSYA0ggpppp,BSYA0ggppmp,loab,loba,.true.)
      call gluegen(e1,2,3,e4,za,zb,zab,zba,
     & A41ggpppp,A41ggppmp,ab41,ba41,.true.)
      call gluegen(e1,2,3,e4,za,zb,zab,zba,
     & A43ggpppp,A43ggppmp,ab43,ba43,.false.)


c--- note that filling of colour structures here looks unnatural:
c---    1 <--> loba , 2 <--> loab
c---  but this does correspond to the filling in qqb_QQb.f
      resgg0=0d0
      resgg1=0d0
      do j=1,2
      do k=1,2
      loqed(j,k)=loab(j,k)+loba(j,k)

      resgg0=resgg0
     & +xn*fac*(abs(loab(j,k))**2+abs(loba(j,k))**2
     &           -abs(loqed(j,k))**2/xn**2)

      resgg1=resgg1+xn**2*fac*dreal(
     & +ab41(j,k)*dconjg(loab(j,k))+ba41(j,k)*dconjg(loba(j,k))
     & +(ab43(j,k)-ab41(j,k)-ba41(j,k))*dconjg(loqed(j,k))/xn**2)

      enddo
      enddo

c--- renormalization
      ren=(
     & +2d0*((nlf/xn-11d0/6d0*xn)*epinv+xn/6d0)
     & +2d0*nhf/xn*(epinv+log(musq/mt**2))
     & -Cf*(3d0*(epinv+log(musq/mt**2))+5d0))
      resqqb1=resqqb1+ren*resqqb0
      resqbq1=resqbq1+ren*resqqb0
     
      ren=(
     & +2d0*((nlf/3d0-11d0/6d0*xn)*epinv+xn/6d0)
     & -(xn**2-1d0)/2d0/xn*(3d0*(epinv+log(musq/mt**2))+5d0))
      resgg1=resgg1+ren*resgg0
      
C----set all elements to zero
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      enddo
      enddo

C---fill qb-q, gg and q-qb elements
      do j=-nf,nf
      if (j .lt. 0) then
          msq(j,-j)=ason2pi*aveqq*resqbq1
c          msq(j,-j)=aveqq*resqqb0
      elseif (j .eq. 0) then
          msq(j,j)=ason2pi*avegg*resgg1
c          msq(j,j)=avegg*resgg0
      elseif (j .gt. 0) then
          msq(j,-j)=ason2pi*aveqq*resqqb1
c          msq(j,-j)=aveqq*resqqb0
      endif
      enddo
      
      return
      end
