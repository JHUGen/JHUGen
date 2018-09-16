      subroutine interdk_gg(pp,me12,me21,intf)
      implicit none
      include 'types.f'
c--- Wrapper for the gluon-gluon initiated reals for t-channel
c--- single top with massive b-quark. Includes both gluon
c--- permutations.
c--- By R. Frederix, July 16, 2008.
c---    added decay: J. Campbell, May 2011

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
      real(dp):: me12,me21,intf,pp(mxpart,4),
     &   dot,projb,projt,wprop12,wprop21,mq,ma,gsq_L,gsq_H,s345,fac
      complex(dp):: gs1(2,2,2,2),gs2(2,2,2,2),gs1dk(2,2,2),
     & gs2dk(2,2,2)
      integer:: i,j,k
      integer:: u,g1,g2,d,b1,t1,proj,i3,i4

c color matrix:
c      integer:: CF(2,2)
c      DATA (CF(i,1  ),i=1  ,2  ) /     16,    2/
c      DATA (CF(i,2  ),i=1  ,2  ) /      2,   16/



c particle identifiers:
      u=8
      d=7

c set mass of quark and antiquark according to nwz
      if (nwz == +1) then
        mq=mt
        ma=mb
        i3=3
        i4=4
        t1=9
        b1=10
      else
        mq=mb
        ma=mt
        i3=4
        i4=3
        t1=10
        b1=9
      endif

      proj=1 ! momentum to project the spin of top and bottom quarks

c project spin of the t and b on proj (this defines t1 and b1,
c used in the call to the reals1)
      projb=dot(pp,proj,6)
      projt=dot(pp,proj,3)+dot(pp,proj,4)+dot(pp,proj,5)
      do i=1,4
         pp(t1,i)=pp(3,i)+pp(4,i)+pp(5,i)-pp(proj,i)*mt**2/2._dp/projt
         pp(b1,i)=pp(6,i)-pp(proj,i)*mb**2/2._dp/projb
      enddo

c-------------------------
c First gluon permutation:
c-------------------------
      g1=1  !g1 is gluon attached to massive fermion line
      g2=2  !g2 is gluon attached to massless fermion line


c Get the spinor products (and invariants):
      call spinoru(10,pp,za,zb)

c W boson propagator:
      wprop12= wmass**2 - (s( u, d)+s( u,g2)+s(d,g2))

c gluon attached to w current
      call reals1(u,g1,d,9,10,g2,proj,mq,ma,za,zb,gs1)




c--------------------------
c Second gluon permutation:
c--------------------------
      g1=2  !g1 is gluon attached to massive fermion line
      g2=1  !g2 is gluon attached to massless fermion line

c Get the spinor products (and invariants):
      call spinoru(10,pp,za,zb)

c W boson propagator:
      wprop21= wmass**2 - (s( u, d)+s( u,g2)+s(d,g2))

c gluon attached to w current
      call reals1(u,g1,d,9,10,g2,proj,mq,ma,za,zb,gs2)


      s345=(pp(3,4)+pp(4,4)+pp(5,4))**2
     &    -(pp(3,1)+pp(4,1)+pp(5,1))**2
     &    -(pp(3,2)+pp(4,2)+pp(5,2))**2
     &    -(pp(3,3)+pp(4,3)+pp(5,3))**2

c--- now dress up with appropriate factors to include the top quark decay
      fac=gwsq*sqrt(2._dp*dot(pp,i3,5))
     &    /sqrt((2._dp*dot(pp,3,4)-wmass**2)**2+(wmass*wwidth)**2)
     &    /sqrt((s345-mt**2)**2+(mt*twidth)**2)

      do i=1,2
      do j=1,2
      do k=1,2
      if     (nwz == +1) then
        gs1dk(i,j,k)=fac*(
     &   +gs1(i,2,j,k)*zb(i4,t1)
     &   +gs1(i,1,j,k)*mt*zb(i4,proj)/zb(t1,proj))
        gs2dk(i,j,k)=fac*(
     &   +gs2(i,2,j,k)*zb(i4,t1)
     &   +gs2(i,1,j,k)*mt*zb(i4,proj)/zb(t1,proj))
      elseif (nwz == -1) then
        gs1dk(i,j,k)=fac*(
     &   +gs1(i,j,1,k)*za(i4,t1)
     &   +gs1(i,j,2,k)*mt*za(i4,proj)/za(t1,proj))
        gs2dk(i,j,k)=fac*(
     &   +gs2(i,j,1,k)*za(i4,t1)
     &   +gs2(i,j,2,k)*mt*za(i4,proj)/za(t1,proj))
      else
        write(6,*) 'nwz must be +1 or -1 in interdk_gg'
        stop
      endif
      enddo
      enddo
      enddo


c----------------------------------------
c Square and sum all helicity amplitudes:
c----------------------------------------
      me12=0._dp
      me21=0._dp
      intf=0._dp
      do i=1,2
        do j=1,2
          do k=1,2
c color diagonal terms:
             me12=me12+
     &            abs(gs1dk(i,j,k))**2/wprop12**2*16._dp
             me21=me21+
     &            abs(gs2dk(j,k,i))**2/wprop21**2*16._dp
c off-diagonal terms (also interchange i <--> l):
             intf=intf+4._dp*real(
     &            gs1dk(i,j,k)*conjg(gs2dk(j,k,i)))/wprop12/wprop21
          enddo
        enddo
      enddo

c--- note that these are all corrections to the W current and thus
c--- should receive a factor of gsq_L
      gsq_L=fourpi*as_L
      gsq_H=fourpi*as_H

c Multiply by coupling constants:
c      me=me*GG(1)**4*GWF(1)**4    !MG/ME value
      me12=me12*gsq_L*gsq_H*gwsq**2/4._dp    !MCFM
      me21=me21*gsq_L*gsq_H*gwsq**2/4._dp    !MCFM
      intf=intf*gsq_L*gsq_H*gwsq**2/4._dp    !MCFM


c 'Average' over incoming helicities and color:
      me12=me12/256._dp
      me21=me21/256._dp
      intf=intf/256._dp

c--- ensure array pp is not contaminated upon return
      do k=1,4
      pp(t1,k)=zip
      pp(b1,k)=zip
      enddo

      return
      end

