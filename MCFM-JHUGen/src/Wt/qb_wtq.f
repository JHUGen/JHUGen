      subroutine qb_wtq(mq,qwidth,p,ix,is,ie,in,jn,je,jb,iy,msq)
      implicit none
      include 'types.f'
************************************************************************
*     Author: Francesco Tramontano                                     *
*     February, 2005.                                                  *
*     Real correction to W+t, radiation in production                  *
*                                                                      *
*    Matrix element squared and averaged over initial colours and spins*
*--- W+t production (nwz=-1)                                           *
*     q(-ix) + b(-is) --> W + t(pneb) + q(iy)                          *
*                         |   |                                        *
*                         |   --> nu(jn) + e^+(jn) + b(jb)             *
*                         |                                            *
*                         --> e^-(ie) + nubar(in)                      *
*--- W+tbar production (nwz=+1)                                        *
*     q(-p1) + qbar(-p2) --> W + t(p567) + f(p8)                       *
*                            |   |                                     *
*                            |   --> e^-(p5) + nubar(p6) + b(p7)       *
*                            |                                         *
*                            --> nu(p3) + e^+(p4)                      *
************************************************************************
c---- helicities: 1=minus 2=plus

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'masses.f'
      integer:: ix,is,ie,in,jn,je,jb,iy
      real(dp):: p(mxpart,4),fac,prop,dot,msq,mq,qwidth
      complex(dp):: amp(2)
      fac=aveqq*gsq**2*gw**8*xn*cf*half
      prop=(two*dot(p,3,4)-wmass**2)**2+(wmass*wwidth)**2
      prop=prop*((two*dot(p,5,6)-wmass**2)**2+(wmass*wwidth)**2)
      prop=prop*((two*(dot(p,5,6)+dot(p,5,7)+dot(p,6,7))-mq**2)**2
     & +(mq*qwidth)**2)
      call amps_4quark(mq,qwidth,p,ix,is,ie,in,jn,je,jb,iy,amp)
      msq=abs(amp(1))**2+abs(amp(2))**2
      msq=fac*msq/prop
      return
      end


      subroutine amps_4quark(mq,qwidth,p,ix,is,ie,in,jn,je,jb,iy,amp)
      implicit none
      include 'types.f'

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      integer:: j,ix,is,ie,in,jn,je,jb,iy
      real(dp):: p(mxpart,4),dot,xDy,txy2,sxy2,t(4),q(mxpart,4),
     & mq,qwidth
      complex(dp):: amp(2),zab(mxpart,mxpart),zba(mxpart,mxpart)
      do j=1,4
      t(j)=p(jn,j)+p(je,j)+p(jb,j)
      q(1,j)=t(j)
      q(2,j)=p(ix,j)
      q(3,j)=p(iy,j)
      enddo
      xDy=dot(p,ix,iy)
      sxy2=two*(dot(p,is,ix)+dot(p,is,iy)+dot(p,ix,iy))
      txy2=two*(dot(q,1,2)+dot(q,1,3)+dot(q,2,3))
      call spinoru(8,p,za,zb)
      call spinork(8,p,zab,zba,t)
      amp(1)=  + zab(iy,je)*xDy**(-1) * ( za(jb,jn)*za(iy,ie
     &    )*zb(in,is)*zb(ix,iy) )
      amp(1) = amp(1) + zab(ie,je)*sxy2**(-1)*xDy**(-1)*txy2 * (
     &     - za(jb,jn
     &    )*za(ix,iy)*zb(in,ix)*zb(ix,is) - za(jb,jn)*za(iy,iy)*zb(in,
     &    iy)*zb(ix,is) - za(jb,jn)*za(is,iy)*zb(in,is)*zb(ix,is) )
      amp(1) = amp(1) + zab(ie,je)*zab(iy,ix)*xDy**(-1) * (
     &    za(jb,jn)*zb(in,is) )

      amp(2)=  + zab(ix,je)*xDy**(-1) * ( za(jb,jn)*za(ix,ie
     &    )*zb(in,is)*zb(iy,ix) )
      amp(2) = amp(2) + zab(ie,je)*sxy2**(-1)*xDy**(-1)*txy2 * (
     &     - za(jb,jn
     &    )*za(ix,ix)*zb(in,ix)*zb(iy,is) - za(jb,jn)*za(iy,ix)*zb(in,
     &    iy)*zb(iy,is) - za(jb,jn)*za(is,ix)*zb(in,is)*zb(iy,is) )
      amp(2) = amp(2) + zab(ie,je)*zab(ix,iy)*xDy**(-1) * (
     &    za(jb,jn)*zb(in,is) )

c--- Use the "overall scheme" to include the top width when necessary
      if (txy2+mq**2 > zero) then
        txy2=sqrt(txy2**2+(mq*qwidth)**2)
      endif
      amp(1)=amp(1)/txy2
      amp(2)=amp(2)/txy2

      return
      end
