      subroutine interdk(pp,u,g1,i3,i4,d,g2,mq,ma,me)
c--- Wrapper for the quark-gluon initiated reals for t-channel
c--- single top with massive b-quark.
c--- By R. Frederix, July 16, 2008.
c---    added decay: J. Campbell, May 2011 
      implicit none
      include 'constants.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'nwz.f'
      include 'stopscales.f'
      double precision me
      double complex gs1(2,2,2,2),gs2(2,2,2,2),gs3(2,2,2,2),
     & gs1dk(2,2,2),gs2dk(2,2,2),gs3dk(2,2,2)
      double precision pp(mxpart,4),dot,s345,fac,
     .  g1g2b,g1g2t,g2b,g2t,g1b,g1t,ud,udg2,mq,ma,gsq_L,gsq_H
      integer i,j,k
      integer u,g1,g2,d,b1,t1,proj,i3,i4

c color matrix:
c      integer CF(3,3)
c      DATA (CF(i,1  ),i=1  ,3  ) /     16,    0,    0/    
c      DATA (CF(i,2  ),i=1  ,3  ) /      0,   16,   -2/    
c      DATA (CF(i,3  ),i=1  ,3  ) /      0,   -2,   16/    

      if (nwz .eq. +1) then
        t1=9
        b1=10
      else
        t1=10
        b1=9
      endif
      
c momentum to project the spin of top and bottom quarks
      proj=g1  ! hard-coded in reals2 and reals3 to g1


c project spin of the t and b on proj (this defines t1 and b1,
c used in the call to the reals1, reals2 and reals3)
      g1b=dot(pp,proj,6)
      g1t=dot(pp,proj,3)+dot(pp,proj,4)+dot(pp,proj,5)
      do i=1,4
         pp(t1,i)=pp(3,i)+pp(4,i)+pp(5,i)-pp(proj,i)*mt**2/2d0/g1t
         pp(b1,i)=pp(6,i)-pp(proj,i)*mb**2/2d0/g1b
      enddo

c Get the spinor products (and invariants):
      call spinoru(10,pp,za,zb)


c inner products for the 'common terms' below:
      ud    = s( u, d)/2d0
      udg2  = s( u, d)/2d0+s( u,g2)/2d0+s(d,g2)/2d0
      g2b   = s(g2,10)/2d0+s(g1,g2)*ma**2/2d0/s(g1,10)
      g2t   = s(g2, 9)/2d0+s(g1,g2)*mq**2/2d0/s(g1, 9)
      g1g2b = s(g1,g2)/2d0+s(g1,10)/2d0+g2b
      g1g2t = s(g1,g2)/2d0+s(g1, 9)/2d0+g2t
c--- redefine g1b and g1t to handle anti-top case
      g1b   = s(g1,10)/2d0
      g1t   = s(g1, 9)/2d0


c gluon attached to w current
      call reals1(u,g1,d,9,10,g2,proj,mq,ma,za,zb,gs1)

c color structure 1
      call reals2(u,g1,d,9,10,g2,mq,ma,za,zb,gs2)

c color structure 2
      call reals3(u,g1,d,9,10,g2,mq,ma,za,zb,gs3)

c--- note that gs1 corresponds to corrections to the W current and thus
c--- should receive a factor of gsq_L, all others should be gsq_H
      gsq_L=fourpi*as_L
      gsq_H=fourpi*as_H

      s345=(pp(3,4)+pp(4,4)+pp(5,4))**2
     &    -(pp(3,1)+pp(4,1)+pp(5,1))**2
     &    -(pp(3,2)+pp(4,2)+pp(5,2))**2
     &    -(pp(3,3)+pp(4,3)+pp(5,3))**2
     
c--- now dress up with appropriate factors to include the top quark decay
      fac=gwsq*dsqrt(2d0*dot(pp,i3,5))
     &    /dsqrt((2d0*dot(pp,3,4)-wmass**2)**2+(wmass*wwidth)**2)
     &    /dsqrt((s345-mt**2)**2+(mt*twidth)**2)
     
      do i=1,2
      do j=1,2
      do k=1,2
      if     (nwz .eq. +1) then
        gs1dk(i,j,k)=fac*(
     &   +gs1(i,2,j,k)*zb(i4,t1)
     &   +gs1(i,1,j,k)*mt*zb(i4,proj)/zb(t1,proj))
        gs2dk(i,j,k)=fac*(
     &   +gs2(i,2,j,k)*zb(i4,t1)
     &   +gs2(i,1,j,k)*mt*zb(i4,proj)/zb(t1,proj))
        gs3dk(i,j,k)=fac*(
     &   +gs3(i,2,j,k)*zb(i4,t1)
     &   +gs3(i,1,j,k)*mt*zb(i4,proj)/zb(t1,proj))
      elseif (nwz .eq. -1) then
        gs1dk(i,j,k)=fac*(
     &   +gs1(i,j,1,k)*za(i4,t1)
     &   +gs1(i,j,2,k)*mt*za(i4,proj)/za(t1,proj))
        gs2dk(i,j,k)=fac*(
     &   +gs2(i,j,1,k)*za(i4,t1)
     &   +gs2(i,j,2,k)*mt*za(i4,proj)/za(t1,proj))
        gs3dk(i,j,k)=fac*(
     &   +gs3(i,j,1,k)*za(i4,t1)
     &   +gs3(i,j,2,k)*mt*za(i4,proj)/za(t1,proj))
      else
        write(6,*) 'nwz must be +1 or -1 in interdk'
        stop
      endif
      enddo
      enddo
      enddo
      
      
      me=0d0
      do i=1,2
         do j=1,2
            do k=1,2
c color diagonal terms:
                  me=me+gsq_L*cdabs(gs1dk(i,j,k))**2*16d0
     &       /((wmass**2-2d0*udg2))**2
                  me=me+gsq_H*cdabs(gs2dk(i,j,k))**2*16d0
     &       /(g1g2b*g1g2t*g1t*g2b*(wmass**2 - 2d0*ud))**2
                  me=me+gsq_H*cdabs(gs3dk(i,j,k))**2*16d0
     &       /(g1g2b*g1g2t*g1b*g2t*(wmass**2 - 2d0*ud))**2
c color suppressed terms:
                 me=me-gsq_H*dble(gs2dk(i,j,k)*dconjg(gs3dk(i,j,k)))*2d0
     &       /(g1g2b**2*g1g2t**2*g1t*g2b*g1b*g2t*(wmass**2 - 2d0*ud)**2)
                 me=me-gsq_H*dble(gs3dk(i,j,k)*dconjg(gs2dk(i,j,k)))*2d0
     &       /(g1g2b**2*g1g2t**2*g1t*g2b*g1b*g2t*(wmass**2 - 2d0*ud)**2)
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

c--- ensure array pp is not contaminated upon return
      do k=1,4
      pp(t1,k)=zip
      pp(b1,k)=zip
      enddo
      
      return
      end

