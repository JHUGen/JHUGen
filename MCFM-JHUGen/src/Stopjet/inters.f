c--- similar to "inter", but set up for s-channel single top + jet;
c---  in particular:
c---     corrections to the light quark line may be removed
c---     u,g1,d,g2 are already fixed for this crossing
c---     average is fixed for "qqbar" rather than "qg"
      subroutine inters(pp,me_lc,me_slc)
      implicit none
      include 'types.f'

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'nwz.f'
      include 'stopscales.f'
      real(dp):: me_lc,me_slc
      complex(dp):: gs1(2,2,2,2),gs2(2,2,2,2),gs3(2,2,2,2)
      real(dp):: pp(mxpart,4),dot,
     & g1g2b,g1g2t,g2b,g2t,g1b,g1t,ud,udg2,mq,ma,gsq_H,gsq_L
      integer:: i,j,k,l
      integer:: u,g1,g2,d,t,b,b1,t1,b2,t2

c color matrix:
c      integer:: CF(3,3)
c      DATA (CF(i,1  ),i=1  ,3  ) /     16,    0,    0/
c      DATA (CF(i,2  ),i=1  ,3  ) /      0,   16,   -2/
c      DATA (CF(i,3  ),i=1  ,3  ) /      0,   -2,   16/

c particle identifiers:-  NB: changed by JC
      u=1
      g1=5
      d=2
      g2=6
c--- set mass of quark and antiquark according to nwz
      if (nwz == +1) then
        mq=mt
        ma=mb
        t=3
        b=4
      else
        mq=mb
        ma=mt
        t=4
        b=3
      endif

      t1=7
      b1=8
      t2=9
      b2=10

c project spin of the t and b on g1 (this defines t1 and b1,
c  used in the call to the reals_qq)
      g1b=dot(pp,g1,b)
      g1t=dot(pp,g1,t)
      do i=1,4
         pp(t1,i)=pp(t,i)-pp(g1,i)*mq**2/2._dp/g1t
         pp(t2,i)=pp(g1,i)
         pp(b1,i)=pp(b,i)-pp(g1,i)*ma**2/2._dp/g1b
         pp(b2,i)=pp(g1,i)
      enddo

c Get the spinor products (and invariants):
      call spinoru(8,pp,za,zb)

c inner products for the 'common terms' below:
      ud    = s( u, d)/2._dp
      udg2  = s( u, d)/2._dp+s( u,g2)/2._dp+s(d,g2)/2._dp
      g2b   = s(g2,b1)/2._dp+s(g1,g2)*ma**2/2._dp/s(g1,b1)
      g2t   = s(g2,t1)/2._dp+s(g1,g2)*mq**2/2._dp/s(g1,t1)
      g1g2b = s(g1,g2)/2._dp+s(g1,b1)/2._dp+g2b
      g1g2t = s(g1,g2)/2._dp+s(g1,t1)/2._dp+g2t


c gluon attached to w current
c      call reals1(u,g1,d,t1,b1,g2,mq,ma,za,zb,gs1)

c color structure 1
      call reals2(u,g1,d,t1,b1,g2,mq,ma,za,zb,gs2)

c color structure 2
      call reals3(u,g1,d,t1,b1,g2,mq,ma,za,zb,gs3)

c--- note that gs1 corresponds to corrections to the W current and thus
c--- should receive a factor of gsq_L, all others should be gsq_H
      gsq_L=fourpi*as_L
      gsq_H=fourpi*as_H

      me_lc=0._dp
      me_slc=0._dp
      do i=1,2
         do j=1,2
            do k=1,2
               do l=1,2
c color diagonal terms:
c                  me=me+gsq_L*gs1(i,j,k,l)*conjg(gs1(i,j,k,l))*16._dp
c     &       /(g1t*(wmass**2-2._dp*udg2))**2
c--- 1/12/08: JC replaced "16._dp" by (xn^2*Cf)*[(2._dp*Cf-Ca)+Ca]/2._dp
c---                  and "-2._dp" by (xn^2*Cf)*(2._dp*Cf-Ca)/2._dp
                me_lc=me_lc
     &       +gsq_H*gs2(i,j,k,l)*conjg(gs2(i,j,k,l))*(xn**2*cf)*Ca/2._dp
     &       /(g1g2b*g1g2t*g1t*g2b*(wmass**2 - 2._dp*ud))**2
                me_lc=me_lc
     &       +gsq_H*gs3(i,j,k,l)*conjg(gs3(i,j,k,l))*(xn**2*cf)*Ca/2._dp
     &       /(g1g2b*g1g2t*g1b*g2t*(wmass**2 - 2._dp*ud))**2

                me_slc=me_slc
     &       +gsq_H*gs2(i,j,k,l)*conjg(gs2(i,j,k,l))
     &                     *(xn**2*cf)*(2._dp*cf-ca)/2._dp
     &       /(g1g2b*g1g2t*g1t*g2b*(wmass**2 - 2._dp*ud))**2
                me_slc=me_slc
     &       +gsq_H*gs3(i,j,k,l)*conjg(gs3(i,j,k,l))
     &                     *(xn**2*cf)*(2._dp*cf-ca)/2._dp
     &       /(g1g2b*g1g2t*g1b*g2t*(wmass**2 - 2._dp*ud))**2
c color suppressed terms:
                me_slc=me_slc
     &       +gsq_H*gs2(i,j,k,l)*conjg(gs3(i,j,k,l))
     &                     *(xn**2*cf)*(2._dp*cf-ca)/2._dp
     &       /(g1g2b**2*g1g2t**2*g1t*g2b*g1b*g2t*(wmass**2 - 2._dp*ud)**2)
                me_slc=me_slc
     &       +gsq_H*gs3(i,j,k,l)*conjg(gs2(i,j,k,l))
     &                     *(xn**2*cf)*(2._dp*cf-ca)/2._dp
     &       /(g1g2b**2*g1g2t**2*g1t*g2b*g1b*g2t*(wmass**2 - 2._dp*ud)**2)

c--- LC and SLC explicit
c  leading color:
c                  me=me+gsq_L*gs1(i,j,k,l)*conjg(gs1(i,j,k,l))*16._dp
c     &       /(g1t*(wmass**2-2._dp*udg2))**2
c                  me=me+gsq_H*gs2(i,j,k,l)*conjg(gs2(i,j,k,l))*18._dp
c     &       /(g1g2b*g1g2t*g1t*g2b*(wmass**2 - 2._dp*ud))**2
c                  me=me+gsq_H*gs3(i,j,k,l)*conjg(gs3(i,j,k,l))*18._dp
c     &       /(g1g2b*g1g2t*g1b*g2t*(wmass**2 - 2._dp*ud))**2
c sub-leading color:
c                  me=me-gsq_H*gs2(i,j,k,l)*conjg(gs3(i,j,k,l))*2._dp
c     &       /(g1g2b**2*g1g2t**2*g1t*g2b*g1b*g2t*(wmass**2 - 2._dp*ud)**2)
c                  me=me-gsq_H*gs3(i,j,k,l)*conjg(gs2(i,j,k,l))*2._dp
c     &       /(g1g2b**2*g1g2t**2*g1t*g2b*g1b*g2t*(wmass**2 - 2._dp*ud)**2)
c                  me=me-gsq_H*gs2(i,j,k,l)*conjg(gs2(i,j,k,l))*2._dp
c     &       /(g1g2b*g1g2t*g1t*g2b*(wmass**2 - 2._dp*ud))**2
c                  me=me-gsq_H*gs3(i,j,k,l)*conjg(gs3(i,j,k,l))*2._dp
c     &       /(g1g2b*g1g2t*g1b*g2t*(wmass**2 - 2._dp*ud))**2
               enddo
            enddo
         enddo
      enddo

c Multiply by coupling constants:
c      me=me*GG(1)**4*GWF(1)**4    !MG/ME value
c      me=me*gsq**2*gwsq**2/4._dp    !MCFM

c--- only one factor of the coupling is included here now, for the
c--- heavy quark line at LO, the other is incorporated above
      me_lc=me_lc*gsq_H*gwsq**2/4._dp    !MCFM
      me_slc=me_slc*gsq_H*gwsq**2/4._dp    !MCFM

c 'Average' over incoming helicities and color, with an extra factor
c of two for identical gluons in the final state
      me_lc=me_lc/36._dp/2._dp
      me_slc=me_slc/36._dp/2._dp

      return
      end

