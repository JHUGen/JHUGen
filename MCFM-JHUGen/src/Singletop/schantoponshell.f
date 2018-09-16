      subroutine schantoponshell(q1,q2,p,iswitch,m)
      implicit none
      include 'types.f'
      
************************************************************************
*     Author: R.K. Ellis, January 2012                                 *
*                                                                      *
*     u(-p1)+d~(-p2)--> t(nu,eb,pb)+pc~(p6)                            *
*                                                                      *
*     keeping polarization information for t                           *
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
     & alt,alc,dot,s12,mt2
      complex(dp):: m(2,2),iza,cprop
      integer:: p1,p2,t,eb,c,si,i,j,iswitch,q1,q2
      parameter(p1=1,p2=2,t=3,eb=4,c=5)
C-----matrix element for .e+_dpu~ -> t+b~ where both t and b~ are on shell
C-----t rendered massless wrt eb, and b rendered massless wrt p1
c--- statement functions
      iza(i,j)=cone/za(i,j)
c--- end statement functions
     
C---zero all arrays
      do i=1,2
      do j=1,2
      m(i,j)=czip
      enddo
      enddo
      do si=1,4
      q(p1,si)=p(q1,si)
      q(p2,si)=p(q2,si)
      if (iswitch == 0) then
      q(t,si)=p(3,si)+p(4,si)+p(5,si)
      elseif (iswitch == 1) then
      q(t,si)=p(3,si)+p(4,si)+p(5,si)+p(7,si)
      endif
      q(eb,si)=p(4,si)
      q(c,si)=p(6,si)
      enddo
      mt2=mt**2
      s12=2._dp*dot(q,p1,p2)
      cprop=cplx2(s12-wmass**2,wmass*wwidth)
     
C---- now render "t" massless wrt to vector eb
C---- now render "pc" massless wrt to vector p1
      alt=mt2/(2._dp*dot(q,t,eb))
      alc=mb**2/(2._dp*dot(q,c,p1))
      do si=1,4
      q(t,si)=q(t,si)-alt*q(eb,si)
      q(c,si)=q(c,si)-alc*q(p1,si)
      enddo
      call spinoru(5,q,za,zb)
      
C----order of indices is polt,polc
      m(1,2)= - za(t,p2)*zb(c,p1)*cprop**(-1)

      m(1,1)=czip

      
      m(2,2)=za(eb,p2)*zb(c,p1)*iza(t,eb)*cprop**(-1)*mt

      m(2,1)=czip

      return
      end
