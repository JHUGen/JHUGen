      subroutine schanatoponshellv(q1,q2,p,m,mv)
      implicit none
      include 'types.f'
      
************************************************************************
*     Author: R.K. Ellis, February 2012                                *
*                                                                      *
*     d(-p1)+u~(-p2)--> t~(e-,nb,pc)+pb(p6)                            *
*                                                                      *
*     keeping polarization information for t                           *
************************************************************************
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'zprods_decl.f'
      include 'includect.f'
      real(dp):: p(mxpart,4),q(mxpart,4),
     & ala,alb,dot,s12,mt2,ct
      complex(dp):: m(2,2),mv(2,2),iza,izb,wprop
      complex(dp):: X0L,C0L,C0R,C1L,C1R,C1Lon2,C1Ron2
      integer:: p1,p2,a,ee,b,si,i,j,q1,q2
      logical:: oldincludect
      parameter(p1=1,p2=2,ee=3,a=4,b=5)
C-----matrix element for .e+_dpu~ -> t~+b where both t~ and b are on shell
C-----t rendered massless wrt ee, and b rendered massless wrt p2
c--- statement functions
      iza(i,j)=cone/za(i,j)
      izb(i,j)=cone/zb(i,j)
c--- end statement functions

c--- corrections in production do not need CT to be included     
      oldincludect=includect
      includect=.false.

C---zero all arrays
      do i=1,2
      do j=1,2
      m(i,j)=czip
      mv(i,j)=czip
      enddo
      enddo
      do si=1,4
      q(p1,si)=p(q1,si)
      q(p2,si)=p(q2,si)
      q(a,si)=p(3,si)+p(4,si)+p(5,si)
      q(ee,si)=p(3,si)
      q(b,si)=p(6,si)
      enddo
      mt2=mt**2
      s12=2._dp*dot(q,p1,p2)
      call coefsdkmass(s12,0._dp,0._dp,ct,X0L,c0R,c1L,C1R)
      call coefsdkmass(s12,mt,mb,ct,C0L,c0R,c1L,C1R)
      C1Lon2=0.5_dp*c1L
      C1Ron2=0.5_dp*c1R
      wprop=cplx2(s12-wmass**2,wmass*wwidth)
     
C---- now render "a" massless wrt to vector e
C---- now render "pb" massless wrt to vector p2
      ala=mt2/(2._dp*dot(q,a,ee))
      alb=mb**2/(2._dp*dot(q,b,p2))
      do si=1,4
      q(a,si)=q(a,si)-ala*q(ee,si)
      q(b,si)=q(b,si)-alb*q(p2,si)
      enddo
      call spinoru(5,q,za,zb)
      
C----order of indices is polb,pola
      m(1,2)= - za(b,p2)*zb(a,p1)*wprop**(-1)

      m(1,1)=za(b,p2)*zb(ee,p1)*izb(a,ee)*wprop**(-1)*mt

      m(2,2)=czip

      m(2,1)=czip

      mv(1,2)=za(b,ee)*za(b,p2)*zb(b,p1)*iza(a,ee)*c1Ron2*wprop**(-1)
     &  - za(b,p2)*zb(b,p1)*zb(a,p2)*izb(b,p2)*c1Lon2*wprop**(-1)*
     & mt**(-1)*mb - za(b,p2)*zb(a,p1)*x0L*wprop**(-1) - za(b,p2)*zb(a,
     & p1)*c0L*wprop**(-1) + za(ee,p2)*zb(p1,p2)*iza(a,ee)*izb(b,p2)*
     & c0R*wprop**(-1)*mt*mb

      mv(1,1)= - za(b,a)*za(b,p2)*zb(b,p1)*c1Ron2*wprop**(-1)*mt**(-1)
     &  + za(b,p2)*zb(b,p1)*zb(ee,p2)*izb(b,p2)*izb(a,ee)*c1Lon2*
     & wprop**(-1)*mb + za(b,p2)*zb(ee,p1)*izb(a,ee)*x0L*wprop**(-1)*mt
     &  + za(b,p2)*zb(ee,p1)*izb(a,ee)*c0L*wprop**(-1)*mt - za(a,p2)*
     & zb(p1,p2)*izb(b,p2)*c0R*wprop**(-1)*mb

      mv(2,2)= - za(b,p2)*zb(b,a)*zb(b,p1)*c1Lon2*wprop**(-1)*mt**(-1)
     &  + za(ee,p2)*zb(b,p1)*iza(a,ee)*c1Ron2*wprop**(-1)*mb + za(ee,p2
     & )*zb(b,p1)*iza(a,ee)*c0R*wprop**(-1)*mt

      mv(2,1)=za(b,p2)*zb(b,ee)*zb(b,p1)*izb(a,ee)*c1Lon2*wprop**(-1)
     &  - za(a,p2)*zb(b,p1)*c1Ron2*wprop**(-1)*mt**(-1)*mb - za(a,p2)*
     & zb(b,p1)*c0R*wprop**(-1)

c--- restore original value of includect
      includect=oldincludect

      return
      end
