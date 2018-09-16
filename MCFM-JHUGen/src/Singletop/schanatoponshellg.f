      subroutine schanatoponshellg(q1,q2,q7,p,iswitch,mi,mf)
      implicit none
      include 'types.f'
      
************************************************************************
*     Author: R.K. Ellis, January 2012                                 *
*                                                                      *
*     d(-p1)+u~(-p2)--> t~(ee,nb,pc)+pb(p6)+g(p7)                      *
*                                                                      *
*     amplitudes with radiation in initial state mi, final state mf    *
*      (returned separately since they do not interfere)               *
*                                                                      *
*     keeping polarization information for t (also nmaed "a"         *
*     iswitch= 0 for no gluon emission                                 *
*     iswitch=+1 for gluon emission in top decay                       *
************************************************************************
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'zprods_decl.f'
      real(dp):: p(mxpart,4),q(mxpart,4),
     & ala,alb,dot,s12,s127,mt2,twopbDg,twopaDg
      complex(dp):: mi(2,2,2),mf(2,2,2),iza,izb,cprop,cprop7
      integer:: p1,p2,a,ee,b,g,si,i,j,k,iswitch,q1,q2,q7
      parameter(p1=1,p2=2,ee=3,a=4,b=5,g=6)
C-----matrix element for p1+p2(g) -> t+c where both t and b are on shell
C-----t rendered massless wrt ee, and p rendered massless wrt p2
c--- statement functions
      iza(i,j)=cone/za(i,j)
      izb(i,j)=cone/zb(i,j)
c--- end statement functions
c      write(6,*) 'iswitch',iswitch
c      call writeout(p)
c      pause
C---zero all arrays
      do i=1,2
      do j=1,2
      do k=1,2
      mi(i,j,k)=czip
      mf(i,j,k)=czip
      enddo
      enddo
      enddo
      do si=1,4
      q(p1,si)=p(q1,si)
      q(p2,si)=p(q2,si)
      if (iswitch == 0) then
      q(a,si)=p(3,si)+p(4,si)+p(5,si)
      elseif (iswitch == 1) then
      q(a,si)=p(3,si)+p(4,si)+p(5,si)+p(7,si)
      endif
      q(ee,si)=p(3,si)
      q(b,si)=p(6,si)
      q(g,si)=p(q7,si)
      enddo
      mt2=mt**2
      s12=2._dp*dot(q,p1,p2)
      s127=s12+2._dp*dot(q,p1,g)+2._dp*dot(q,p2,g)
      twopaDg=2._dp*dot(q,a,g)
      twopbDg=2._dp*dot(q,b,g)
      cprop=cplx2(s12-wmass**2,wmass*wwidth)
      cprop7=cplx2(s127-wmass**2,wmass*wwidth)
C---- now render "a" massless wrt to vector ee
C---- now render "pb" massless wrt to vector p1
      ala=mt2/(2._dp*dot(q,a,ee))
      alb=mb**2/(2._dp*dot(q,b,p2))
      do si=1,4
      q(a,si)=q(a,si)-ala*q(ee,si)
      q(b,si)=q(b,si)-alb*q(p2,si)
      enddo
      call spinoru(6,q,za,zb)
c      call writeout(q)
c      pause
      
C----order of indices is polb,polg,pola
      mi(1,1,2)= + cprop7**(-1) * ( za(b,g)*zb(a,p1)*izb(g,p2) - za(b,
     &    p2)*zb(a,p1)*zb(p1,p2)*izb(g,p1)*izb(g,p2) )

      mi(1,2,2)= + cprop7**(-1) * ( za(b,p2)*za(p1,p2)*zb(a,p1)*iza(g,
     &    p1)*iza(g,p2) + za(b,p2)*zb(a,g)*iza(g,p1) )

      mi(2,1,2)= + cprop7**(-1) * ( za(g,p2)*zb(a,p1)*iza(b,p2)*izb(g,
     &    p2)*mb )

      mi(2,2,2)= czip 

      
      mf(1,1,2)= + cprop**(-1) * ( za(b,g)*za(b,p2)*zb(b,p2)*zb(a,p1)*
     &    izb(g,p2)*twopbDg**(-1) + za(b,g)*za(g,p2)*zb(a,p1)*
     &    twopbDg**(-1) - za(b,p2)*za(a,g)*zb(a,p1)*zb(a,p2)*izb(g,p2)*
     &    twopaDg**(-1) + za(b,p2)*za(g,ee)*zb(a,p2)*zb(ee,p1)*izb(g,p2
     &    )*ala*twopaDg**(-1) - za(b,p2)*za(g,ee)*zb(p1,p2)*iza(a,ee)*
     &    izb(g,p2)*mt2*twopaDg**(-1) )

      mf(1,2,2)= + cprop**(-1) * (  - za(b,p1)*za(b,p2)*zb(b,g)*zb(a,p1
     &    )*iza(g,p1)*twopbDg**(-1) + za(b,p2)*za(a,p1)*zb(a,g)*zb(a,p1
     &    )*iza(g,p1)*twopaDg**(-1) + za(b,p2)*za(ee,p1)*zb(a,g)*zb(ee,
     &    p1)*iza(g,p1)*ala*twopaDg**(-1) + za(b,p2)*za(ee,p1)*zb(g,p1)
     &    *iza(a,ee)*iza(g,p1)*mt2*twopaDg**(-1) + za(b,p2)*zb(a,g)*zb(
     &    g,p1)*twopaDg**(-1) + za(p1,p2)*zb(a,p1)*zb(g,p2)*iza(g,p1)*
     &    izb(b,p2)*twopbDg**(-1)*mb**2 )

      mf(2,1,2)= + cprop**(-1) * ( za(g,p2)**2*zb(a,p1)*iza(b,p2)*
     &    twopbDg**(-1)*mb )

      mf(2,2,2)= czip 

      
      mi(1,1,1)= + cprop7**(-1) * (  - za(b,g)*zb(ee,p1)*izb(a,ee)*izb(
     &    g,p2)*mt + za(b,p2)*zb(ee,p1)*zb(p1,p2)*izb(a,ee)*izb(g,p1)*
     &    izb(g,p2)*mt )

      mi(1,2,1)= + cprop7**(-1) * (  - za(b,p2)*za(p1,p2)*zb(ee,p1)*
     &    iza(g,p1)*iza(g,p2)*izb(a,ee)*mt + za(b,p2)*zb(g,ee)*iza(g,p1
     &    )*izb(a,ee)*mt )

      mi(2,1,1)= + cprop7**(-1) * (  - za(g,p2)*zb(ee,p1)*iza(b,p2)*
     &    izb(a,ee)*izb(g,p2)*mt*mb )

      mi(2,2,1)= czip 

      
      mf(1,1,1)= + cprop**(-1) * (  - za(b,g)*za(b,p2)*zb(b,p2)*zb(ee,
     &    p1)*izb(a,ee)*izb(g,p2)*mt*twopbDg**(-1) - za(b,g)*za(g,p2)*
     &    zb(ee,p1)*izb(a,ee)*mt*twopbDg**(-1) + za(b,p2)*za(a,g)*zb(a,
     &    p1)*zb(ee,p2)*izb(a,ee)*izb(g,p2)*mt*twopaDg**(-1) - za(b,p2)
     &    *za(a,g)*zb(p1,p2)*izb(g,p2)*mt*twopaDg**(-1) - za(b,p2)*za(g
     &    ,ee)*zb(ee,p1)*zb(ee,p2)*izb(a,ee)*izb(g,p2)*ala*mt*
     &    twopaDg**(-1) )

      mf(1,2,1)= + cprop**(-1) * ( za(b,p1)*za(b,p2)*zb(b,g)*zb(ee,p1)*
     &    iza(g,p1)*izb(a,ee)*mt*twopbDg**(-1) + za(b,p2)*za(a,p1)*zb(a
     &    ,p1)*zb(g,ee)*iza(g,p1)*izb(a,ee)*mt*twopaDg**(-1) - za(b,p2)
     &    *za(a,p1)*zb(g,p1)*iza(g,p1)*mt*twopaDg**(-1) + za(b,p2)*za(
     &    ee,p1)*zb(g,ee)*zb(ee,p1)*iza(g,p1)*izb(a,ee)*ala*mt*
     &    twopaDg**(-1) + za(b,p2)*zb(g,ee)*zb(g,p1)*izb(a,ee)*mt*
     &    twopaDg**(-1) - za(p1,p2)*zb(g,p2)*zb(ee,p1)*iza(g,p1)*izb(b,
     &    p2)*izb(a,ee)*mt*twopbDg**(-1)*mb**2 )

      mf(2,1,1)= + cprop**(-1) * (  - za(g,p2)**2*zb(ee,p1)*iza(b,p2)*
     &    izb(a,ee)*mt*twopbDg**(-1)*mb )

      mf(2,2,1)= czip 

      return
      end
