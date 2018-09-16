      subroutine inter_gg(pp,me12,me21,intf)
      implicit none
      include 'types.f'
c--- Wrapper for the gluon-gluon initiated reals for t-channel
c--- single top with massive b-quark. Includes both gluon
c--- permutations.
c--- By R. Frederix, July 16, 2008.


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
     &   dot,projb,projt,wprop12,wprop21,mq,ma,gsq_L,gsq_H
      complex(dp):: gs1(2,2,2,2),gs2(2,2,2,2)
      integer:: i,j,k,l
      integer:: u,g1,g2,d,t,b,b1,t1,proj

c color matrix:
c      integer:: CF(2,2)
c      DATA (CF(i,1  ),i=1  ,2  ) /     16,    2/
c      DATA (CF(i,2  ),i=1  ,2  ) /      2,   16/



c particle identifiers:
      u=6
      d=5

c set mass of quark and antiquark according to nwz
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
      proj=1 ! momentum to project the spin of top and bottom quarks

c project spin of the t and b on proj (this defines t1 and b1,
c used in the call to the reals1)
      projb=dot(pp,proj,b)
      projt=dot(pp,proj,t)
      do i=1,4
         pp(t1,i)=pp(t,i)-pp(proj,i)*mq**2/2._dp/projt
         pp(b1,i)=pp(b,i)-pp(proj,i)*ma**2/2._dp/projb
      enddo

c-------------------------
c First gluon permutation:
c-------------------------
      g1=1  !g1 is gluon attached to massive fermion line
      g2=2  !g2 is gluon attached to massless fermion line


c Get the spinor products (and invariants):
      call spinoru(8,pp,za,zb)

c W boson propagator:
      wprop12= wmass**2 - (s( u, d)+s( u,g2)+s(d,g2))

c gluon attached to w current
      call reals1(u,g1,d,t1,b1,g2,proj,mq,ma,za,zb,gs1)




c--------------------------
c Second gluon permutation:
c--------------------------
      g1=2  !g1 is gluon attached to massive fermion line
      g2=1  !g2 is gluon attached to massless fermion line

c Get the spinor products (and invariants):
      call spinoru(8,pp,za,zb)

c W boson propagator:
      wprop21= wmass**2 - (s( u, d)+s( u,g2)+s(d,g2))

c gluon attached to w current
      call reals1(u,g1,d,t1,b1,g2,proj,mq,ma,za,zb,gs2)



c----------------------------------------
c Square and sum all helicity amplitudes:
c----------------------------------------
      me12=0._dp
      me21=0._dp
      intf=0._dp
      do i=1,2
        do j=1,2
          do k=1,2
            do l=1,2
c color diagonal terms:
               me12=me12+
     &              abs(gs1(i,j,k,l))**2/wprop12**2*16._dp
               me21=me21+
     &              abs(gs2(l,j,k,i))**2/wprop21**2*16._dp
c off-diagonal terms (also interchange i <--> l):
               intf=intf+4._dp*real(
     &              gs1(i,j,k,l)*conjg(gs2(l,j,k,i)))/wprop12/wprop21
            enddo
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
      return
      end

