      subroutine inter(pp,u,g1,t,b,d,g2,mq,ma,me)
c--- Wrapper for the quark-gluon initiated reals for t-channel
c--- single top with massive b-quark.
c--- By R. Frederix, July 16, 2008.
      implicit none
      include 'constants.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'stopscales.f'
      double precision me
      double complex gs1(2,2,2,2),gs2(2,2,2,2),gs3(2,2,2,2)
      double precision pp(mxpart,4),dot,
     .  g1g2b,g1g2t,g2b,g2t,g1b,g1t,ud,udg2,mq,ma,gsq_L,gsq_H
      integer i,j,k,l
      integer u,g1,g2,d,t,b,b1,t1,proj

c color matrix:
c      integer CF(3,3)
c      DATA (CF(i,1  ),i=1  ,3  ) /     16,    0,    0/    
c      DATA (CF(i,2  ),i=1  ,3  ) /      0,   16,   -2/    
c      DATA (CF(i,3  ),i=1  ,3  ) /      0,   -2,   16/    


      t1=7
      b1=8

c momentum to project the spin of top and bottom quarks
      proj=g1  ! hard-coded in reals2 and reals3 to g1


c project spin of the t and b on proj (this defines t1 and b1,
c used in the call to the reals1, reals2 and reals3)
      g1b=dot(pp,proj,b)
      g1t=dot(pp,proj,t)
      do i=1,4
         pp(t1,i)=pp(t,i)-pp(proj,i)*mq**2/2d0/g1t
         pp(b1,i)=pp(b,i)-pp(proj,i)*ma**2/2d0/g1b
      enddo

c Get the spinor products (and invariants):
      call spinoru(8,pp,za,zb)


c inner products for the 'common terms' below:
      ud    = s( u, d)/2d0
      udg2  = s( u, d)/2d0+s( u,g2)/2d0+s(d,g2)/2d0
      g2b   = s(g2,b1)/2d0+s(g1,g2)*ma**2/2d0/s(g1,b1)
      g2t   = s(g2,t1)/2d0+s(g1,g2)*mq**2/2d0/s(g1,t1)
      g1g2b = s(g1,g2)/2d0+s(g1,b1)/2d0+g2b
      g1g2t = s(g1,g2)/2d0+s(g1,t1)/2d0+g2t


c gluon attached to w current
      call reals1(u,g1,d,t1,b1,g2,proj,mq,ma,za,zb,gs1)

c color structure 1
      call reals2(u,g1,d,t1,b1,g2,mq,ma,za,zb,gs2)

c color structure 2
      call reals3(u,g1,d,t1,b1,g2,mq,ma,za,zb,gs3)

c--- note that gs1 corresponds to corrections to the W current and thus
c--- should receive a factor of gsq_L, all others should be gsq_H
      gsq_L=fourpi*as_L
      gsq_H=fourpi*as_H

      
      me=0d0
      do i=1,2
         do j=1,2
            do k=1,2
               do l=1,2
c color diagonal terms:
                  me=me+gsq_L*cdabs(gs1(i,j,k,l))**2*16d0
     &       /((wmass**2-2d0*udg2))**2
                  me=me+gsq_H*cdabs(gs2(i,j,k,l))**2*16d0
     &       /(g1g2b*g1g2t*g1t*g2b*(wmass**2 - 2d0*ud))**2
                  me=me+gsq_H*cdabs(gs3(i,j,k,l))**2*16d0
     &       /(g1g2b*g1g2t*g1b*g2t*(wmass**2 - 2d0*ud))**2
c color suppressed terms:
                 me=me-gsq_H*dble(gs2(i,j,k,l)*dconjg(gs3(i,j,k,l)))*2d0
     &       /(g1g2b**2*g1g2t**2*g1t*g2b*g1b*g2t*(wmass**2 - 2d0*ud)**2)
                 me=me-gsq_H*dble(gs3(i,j,k,l)*dconjg(gs2(i,j,k,l)))*2d0
     &       /(g1g2b**2*g1g2t**2*g1t*g2b*g1b*g2t*(wmass**2 - 2d0*ud)**2)
               enddo
            enddo
         enddo
      enddo

c Multiply by coupling constants:
c      me=me*GG(1)**4*GWF(1)**4    !MG/ME value
c      me=me*gsq**2*gwsq**2/4d0    !MCFM

c--- only one factor of the coupling is included here now, for the
c--- heavy quark line at LO, the other is incorporated above
      me=me*gsq_H*gwsq**2/4d0    !MCFM


c 'Average' over incoming helicities and color:
      me=me/96d0

      return
      end

