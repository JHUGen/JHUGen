      subroutine qqb_w_twdk_gvec(p,n,in,msq)
      implicit none
      include 'types.f'

************************************************************************
*     Authors: John Campbell & Francesco Tramontano                    *
*     February, 2005.                                                  *
*    Matrix element squared and averaged over initial colours and spins*
*     with gluon index contracted with vector n                        *
*         ip emitter                                                   *
*         kp spectator                                                 *
*         in label of gluon which is contracted with n                 *
*--- W+t production (nwz=-1)                                           *
*     q(-p1) + qbar(-p2) --> W + t(p567)                               *
*                            |   |                                     *
*                            |   --> nu(p5) + e^+(p6) + b(p7)          *
*                            |                                         *
*                            --> e^-(p3) + nubar(p4)                   *
*--- W+tbar production (nwz=+1)                                        *
*     q(-p1) + qbar(-p2) --> W + t(p567)                               *
*                            |   |                                     *
*                            |   --> e^-(p5) + nubar(p6) + b(p7)       *
*                            |                                         *
*                            --> nu(p3) + e^+(p4)                      *
************************************************************************
c---ip emitter
c---kp spectator
c---in label of gluon which is contracted with n
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'nwz.f'
      integer:: j,k,in,i3,i4,i5,i6,iq
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4),n(4),
     & p1p2(-1:1,-1:1),wtgvecn

c--- set up lepton variables depending on whether it's t or tbar
      if     (nwz == -1) then
        i3=3
        i4=4
        i5=5
        i6=6
        iq=+1 ! quark initial state
      elseif (nwz == +1) then
        i3=4
        i4=3
        i5=6
        i6=5
        iq=-1 ! antiquark initial state
      else
        write(6,*) 'Error in qqb_w_twdk_gvec, nwz is not +1 or -1 :',nwz
        stop
      endif

      do j=-1,+1
      do k=-1,+1
      p1p2(j,k)=zero
      enddo
      enddo

      if (in == 1) then
      p1p2(0,iq)=aveqg*wtgvecn(mt,twidth,1,2,i3,i4,i5,i6,7,p,n)
      elseif (in == 2) then
      p1p2(iq,0)=aveqg*wtgvecn(mt,twidth,2,1,i3,i4,i5,i6,7,p,n)
      endif

      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=zero
      enddo
      enddo

      do j=-nf,nf,nf
      do k=-nf,nf,nf
      if     ((j > 0) .and. (k == 0)) then
          msq(j,k)=p1p2(+1,0)
      elseif ((j < 0) .and. (k == 0)) then
          msq(j,k)=p1p2(-1,0)
      elseif ((j == 0) .and. (k > 0)) then
          msq(j,k)=p1p2(0,+1)
      elseif ((j == 0) .and. (k < 0)) then
          msq(j,k)=p1p2(0,-1)
      endif

      enddo
      enddo

      return
      end

      function wtgvecn(mq,qwidth,ig,is,ie,in,jn,je,jb,p,vec)
      implicit none
      include 'types.f'
      real(dp):: wtgvecn

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'masses.f'
      include 'zprods_com.f'
      integer:: ig,is,ie,in,je,jn,jb
      real(dp):: p(mxpart,4),vec(4),prop,fac,mq,qwidth
      complex(dp):: amp
      complex(dp):: zab(mxpart,mxpart),zba(mxpart,mxpart)
      common/zabprods/zab,zba
!$omp threadprivate(/zabprods/)

      call checkndotp(p,vec,ig)

      call spinoru(7,p,za,zb)
      call spinork(7,p,zab,zba,vec)
      call ampsn(mq,p,ig,is,ie,in,jn,je,jb,amp)

      prop=(real(za(ie,in)*zb(in,ie))-wmass**2)**2+(wmass*wwidth)**2
      prop=prop*((real(za(jn,je)*zb(je,jn))-wmass**2)**2
     & +(wmass*wwidth)**2)
      prop=prop*(
     .(real(za(jn,je)*zb(je,jn)+za(jn,jb)*zb(jb,jn)+za(je,jb)*zb(jb,je))
     & -mq**2)**2+(mq*qwidth)**2)

      fac=xn*cf*gsq*gw**8

      wtgvecn=fac*abs(amp)**2/prop

      return
      end


      subroutine ampsn(mq,p,ig,is,ie,in,jn,je,jb,amp)
      implicit none
      include 'types.f'

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_com.f'
      complex(dp):: amp
      real(dp):: p(mxpart,4),dot,taugt,taugs,mq
      integer:: ig,is,ie,in,je,jn,jb
      complex(dp):: zab(mxpart,mxpart),zba(mxpart,mxpart)
      common/zabprods/zab,zba
!$omp threadprivate(/zabprods/)

      taugs=two*dot(p,ig,is)
      taugt=two*(dot(p,ig,jn)+dot(p,ig,je)+dot(p,ig,jb))
      amp= za(ig,ie)*za(jn,jb)*zb(is,in)*zb(je,jn)*zab(jn,ig)*
     & taugt**(-1) + za(ig,ie)*za(jn,jb)*zb(is,in)*zb(je,jb)*zab(jb,ig)
     & *taugt**(-1) - za(ie,je)*za(jn,jb)*zb(is,in)*zb(je,jn)*zab(jn,je
     & )*taugt**(-1) - za(ie,je)*za(jn,jb)*zb(is,in)*zb(je,jb)*zab(jb,
     & je)*taugt**(-1) + za(ie,jn)*za(jn,jb)*zb(ig,in)*zb(je,jn)*zab(ig
     & ,is)*taugs**(-1) + za(ie,jn)*za(jn,jb)*zb(is,in)*zb(je,jn)*zab(
     & is,is)*taugs**(-1) - za(ie,jn)*za(jn,jb)*zb(is,in)*zb(je,jn)*
     & zab(jn,jn)*taugt**(-1) - za(ie,jn)*za(jn,jb)*zb(is,in)*zb(je,jb)
     & *zab(jb,jn)*taugt**(-1) + za(ie,jb)*za(jn,jb)*zb(ig,in)*zb(je,jb
     & )*zab(ig,is)*taugs**(-1) - za(ie,jb)*za(jn,jb)*zb(is,in)*zb(je,
     & jn)*zab(jn,jb)*taugt**(-1) + za(ie,jb)*za(jn,jb)*zb(is,in)*zb(je
     & ,jb)*zab(is,is)*taugs**(-1) - za(ie,jb)*za(jn,jb)*zb(is,in)*zb(
     & je,jb)*zab(jb,jb)*taugt**(-1) + za(jn,jb)*zb(is,in)*zab(ie,je)*
     & mq**2*taugt**(-1)

      return
      end
