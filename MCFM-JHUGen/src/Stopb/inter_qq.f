      subroutine inter_qq(pp,u,q1,t,b,d,q2,mq,ma, me,intf)
c--- Wrapper for the quark-quark initiated reals for t-channel
c--- single top with massive b-quark.
c--- By R. Frederix, July 17, 2008.
      implicit none
      include 'constants.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'ewcouple.f'
      include 'masses.f'
      include 'stopscales.f'
      double precision me(4),intf(4),mq,ma,
     .     wprop_a,wprop_b,wprop_c,wprop_d
      double complex gs_a(2,2,2),gs_b(2,2,2),gs_c(2,2,2),gs_d(2,2,2)
      double precision pp(mxpart,4),dot,bDg,tDg,ud,q1b,q1t,gsq_H
      integer i,j,k
      integer u,q1,q2,d,t,b,b1,t1,proj

c color matrix:
c      DATA (CF(i,1  ),i=1  ,2  ) /     6,   2    /                
c      DATA (CF(i,2  ),i=1  ,2  ) /     2,   6    /                

      t1=7
      b1=8

c momentum to project the spin of top and bottom quarks
      proj=1

c project spin of the t and b on proj (this defines t1 and b1,
c used in the call to the reals_qq)
      q1b=dot(pp,proj,b)
      q1t=dot(pp,proj,t)
      do i=1,4
         pp(t1,i)=pp(t,i)-pp(proj,i)*mq**2/2d0/q1t
         pp(b1,i)=pp(b,i)-pp(proj,i)*ma**2/2d0/q1b
      enddo

c Get the spinor products (and invariants):
      call spinoru(8,pp,za,zb)

c------------------------
c First configuration
c------------------------
c Nothing to permute

c W propagator and common factors:
      ud   = s( u, d)/2d0
      bDg=(s(q1,q2)+s(b,q1)+s(b,q2))/2d0
      tDg=(s(q1,q2)+s(t,q1)+s(t,q2))/2d0
      wprop_a = bDg*tDg*(wmass**2-2d0*ud)

c Only one structure
      call reals_qq(u,q1,d,t1,b1,q2,proj,mq,ma,za,zb,gs_a)


c------------------------
c Second configuration
c------------------------
c interchange final state d <--> q2

c W propagator and common factors:
      ud   = s( u, q2)/2d0
      bDg=(s(q1,d)+s(b,q1)+s(b,d))/2d0
      tDg=(s(q1,d)+s(t,q1)+s(t,d))/2d0
      wprop_b = bDg*tDg*(wmass**2-2d0*ud)
    
c Only one structure
      call reals_qq(u,q1,q2,t1,b1,d,proj,mq,ma,za,zb,gs_b)


c------------------------
c Third configuration
c------------------------
c interchange final state d <--> q2
c and initial state u <--> q1

c W propagator and common factors:
      ud   = s( q1, q2)/2d0
      bDg=(s(u,d)+s(b,u)+s(b,d))/2d0
      tDg=(s(u,d)+s(t,u)+s(t,d))/2d0
      wprop_c = bDg*tDg*(wmass**2-2d0*ud)

c Only one structure
      call reals_qq(q1,u,q2,t1,b1,d,proj,mq,ma,za,zb,gs_c)


c------------------------
c Fourth configuration
c------------------------
c interchange u <--> q1

c W propagator and common factors:
      ud   = s( q1, d)/2d0
      bDg=(s(u,q2)+s(b,u)+s(b,q2))/2d0
      tDg=(s(u,q2)+s(t,u)+s(t,q2))/2d0
      wprop_d = bDg*tDg*(wmass**2-2d0*ud)

c Only one structure
      call reals_qq(q1,u,d,t1,b1,q2,proj,mq,ma,za,zb,gs_d)


c----------------------------------------
c Square and sum all helicity amplitudes:
c----------------------------------------
      do i=1,4
         me(i)=0d0
         intf(i)=0d0
      enddo

      do k=1,2
         do j=1,2
            do i=1,2
c Diagonal color structure:
               me(1)=me(1)+
     &              cdabs(gs_a(i,j,k))**2*6d0/wprop_a**2
               me(2)=me(2)+
     &              cdabs(gs_b(i,j,k))**2*6d0/wprop_b**2
               me(3)=me(3)+
     &              cdabs(gs_c(i,j,k))**2*6d0/wprop_c**2
               me(4)=me(4)+
     &              cdabs(gs_d(i,j,k))**2*6d0/wprop_d**2
            enddo
c Off-diagonal color structure (Only add '1' helicity):
            intf(1)=intf(1)-dreal(
     &           gs_a(1,j,k)*dconjg(gs_b(1,j,k))*4d0)/wprop_a/wprop_b
            intf(2)=intf(2)-dreal(
     &           gs_b(1,j,k)*dconjg(gs_c(1,j,k))*4d0)/wprop_b/wprop_c
            intf(3)=intf(3)-dreal(
     &           gs_c(1,j,k)*dconjg(gs_d(1,j,k))*4d0)/wprop_c/wprop_d
            intf(4)=intf(4)-dreal(
     &           gs_d(1,j,k)*dconjg(gs_a(1,j,k))*4d0)/wprop_d/wprop_a
         enddo
      enddo

c--- note that these are all corrections to the heavy quark line and
c--- thus should receive a factor of gsq_H
      gsq_H=fourpi*as_H

c Multiply by coupling constants:
      do i=1,4
         me(i)=me(i)*gsq_H**2*gwsq**2/4d0
         intf(i)=intf(i)*gsq_H**2*gwsq**2/4d0
      enddo

c 'Average' over incoming helicities and color:
      do i=1,4
         me(i)=me(i)/36d0
         intf(i)=intf(i)/36d0
      enddo

      return
      end

