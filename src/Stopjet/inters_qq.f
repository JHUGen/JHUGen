c--- similar to "inter", but set up for s-channel single top + jet;
c---  in particular:
c---     u,q1,d,q2 are already fixed for this crossing
      subroutine inters_qq(pp,me)
      implicit none
      include 'constants.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'nwz.f'
      include 'stopscales.f'
      double precision me
      double complex gs(2,2,2)
      double precision pp(mxpart,4),dot,bDq2,tDq2,bDg,tDg,ud,q1b,q1t,
     . mq,ma,gsq_H
      integer i,j,k,l
      integer u,q1,q2,d,t,b,b1,t1,b2,t2

c color matrix:
c      DATA (CF(i,1  ),i=1  ,1  ) /     6    /                

c particle identifiers:
      u=1
      q1=5
      d=2
      q2=6
c--- set mass of quark and antiquark according to nwz
      if (nwz .eq. +1) then
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

c project spin of the t and b on g1 (this defines t1 and b1, used in the call to the reals_qq)
      q1b=dot(pp,q1,b)
      q1t=dot(pp,q1,t)
      do i=1,4
         pp(t1,i)=pp(t,i)-pp(q1,i)*mq**2/2d0/q1t
         pp(t2,i)=pp(q1,i)
         pp(b1,i)=pp(b,i)-pp(q1,i)*ma**2/2d0/q1b
         pp(b2,i)=pp(q1,i)
      enddo

c Get the spinor products (and invariants):
      call spinoru(8,pp,za,zb)


c inner products for the 'common terms' below:
      ud   = s( u, d)/2d0
      bDq2 = s(q2,b1)/2d0+s(q1,q2)*ma**2/2d0/s(q1,b1)
      tDq2 = s(q2,t1)/2d0+s(q1,q2)*mq**2/2d0/s(q1,t1)
      bDg  = s(q1,q2)/2d0+s(q1,b1)/2d0+bDq2
      tDg  = s(q1,q2)/2d0+s(q1,t1)/2d0+tDq2

    
c Only one color structure
      call reals_qq(u,q1,d,t1,b1,q2,q1,mq,ma,za,zb,gs)
      
      me=0d0
      do i=1,2
         do j=1,2
            do k=1,2
c color diagonal terms:
c--- 1/12/08: JC replaced "6d0" by (xn^2*Cf)/2d0
               me=me+gs(i,j,k)*dconjg(gs(i,j,k))*(xn**2*cf)/2d0
     &   /(bDg*tDg*(wmass**2-2d0*ud))**2
            enddo
         enddo
      enddo

c--- note that these are all corrections to the heavy quark line and
c--- thus should receive a factor of gsq_H    
      gsq_H=fourpi*as_H
      
c Multiply by coupling constants:
c      me=me*GG(1)**4*GWF(1)**4    !MG/ME value
      me=me*gsq_H**2*gwsq**2/4d0    !MCFM

c 'Average' over incoming helicities and color:
      me=me/36d0

      return
      end

