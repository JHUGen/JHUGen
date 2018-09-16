      subroutine interdk_qq(pp,u,q1,i3,i4,d,q2,mq,ma, me,intf)
      implicit none
      include 'types.f'
c--- Wrapper for the quark-quark initiated reals for t-channel
c--- single top with massive b-quark.
c--- By R. Frederix, July 17, 2008.

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'ewcouple.f'
      include 'masses.f'
      include 'nwz.f'
      include 'stopscales.f'
      real(dp):: me(4),intf(4),mq,ma,
     &     wprop_a,wprop_b,wprop_c,wprop_d,s345,fac
      complex(dp):: gs_a(2,2,2),gs_b(2,2,2),gs_c(2,2,2),gs_d(2,2,2)
      complex(dp):: gsdk_a(2,2),gsdk_b(2,2),gsdk_c(2,2),gsdk_d(2,2)
      real(dp):: pp(mxpart,4),dot,bDg,tDg,ud,q1b,q1t,gsq_H
      integer:: i,j,k
      integer:: u,q1,q2,d,b1,t1,proj,i3,i4

c color matrix:
c      DATA (CF(i,1  ),i=1  ,2  ) /     6,   2    /
c      DATA (CF(i,2  ),i=1  ,2  ) /     2,   6    /

      if (nwz == +1) then
        t1=9
        b1=10
      else
        t1=10
        b1=9
      endif

c momentum to project the spin of top and bottom quarks
      proj=1

c project spin of the t and b on proj (this defines t1 and b1,
c used in the call to the reals_qq)
      q1b=dot(pp,proj,6)
      q1t=dot(pp,proj,3)+dot(pp,proj,4)+dot(pp,proj,5)
      do i=1,4
         pp(t1,i)=pp(3,i)+pp(4,i)+pp(5,i)-pp(proj,i)*mt**2/2._dp/q1t
         pp(b1,i)=pp(6,i)-pp(proj,i)*mb**2/2._dp/q1b
      enddo

c Get the spinor products (and invariants):
      call spinoru(10,pp,za,zb)

c------------------------
c First configuration
c------------------------
c Nothing to permute

c W propagator and common factors:
      ud   = s( u, d)/2._dp
      bDg=(s(q1,q2)+s(6,q1)+s(6,q2))/2._dp
      tDg=(s(q1,q2)+s(3,q1)+s(3,q2)+s(4,q1)+s(4,q2)+s(5,q1)+s(5,q2))/2._dp
      wprop_a = bDg*tDg*(wmass**2-2._dp*ud)

c Only one structure
      call reals_qq(u,q1,d,9,10,q2,proj,mq,ma,za,zb,gs_a)


c------------------------
c Second configuration
c------------------------
c interchange final state d <--> q2

c W propagator and common factors:
      ud   = s( u, q2)/2._dp
      bDg=(s(q1,d)+s(6,q1)+s(6,d))/2._dp
      tDg=(s(q1,d)+s(3,q1)+s(3,d)+s(4,q1)+s(4,d)+s(5,q1)+s(5,d))/2._dp
      wprop_b = bDg*tDg*(wmass**2-2._dp*ud)

c Only one structure
      call reals_qq(u,q1,q2,9,10,d,proj,mq,ma,za,zb,gs_b)


c------------------------
c Third configuration
c------------------------
c interchange final state d <--> q2
c and initial state u <--> q1

c W propagator and common factors:
      ud   = s( q1, q2)/2._dp
      bDg=(s(u,d)+s(6,u)+s(6,d))/2._dp
      tDg=(s(u,d)+s(3,u)+s(3,d)+s(4,u)+s(4,d)+s(5,u)+s(5,d))/2._dp
      wprop_c = bDg*tDg*(wmass**2-2._dp*ud)

c Only one structure
      call reals_qq(q1,u,q2,9,10,d,proj,mq,ma,za,zb,gs_c)


c------------------------
c Fourth configuration
c------------------------
c interchange u <--> q1

c W propagator and common factors:
      ud   = s( q1, d)/2._dp
      bDg=(s(u,q2)+s(6,u)+s(6,q2))/2._dp
      tDg=(s(u,q2)+s(3,u)+s(3,q2)+s(4,u)+s(4,q2)+s(5,u)+s(5,q2))/2._dp
      wprop_d = bDg*tDg*(wmass**2-2._dp*ud)

c Only one structure
      call reals_qq(q1,u,d,9,10,q2,proj,mq,ma,za,zb,gs_d)


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
      if     (nwz == +1) then
        gsdk_a(i,j)=fac*(
     &   +gs_a(i,2,j)*zb(i4,t1)
     &   +gs_a(i,1,j)*mt*zb(i4,proj)/zb(t1,proj))
        gsdk_b(i,j)=fac*(
     &   +gs_b(i,2,j)*zb(i4,t1)
     &   +gs_b(i,1,j)*mt*zb(i4,proj)/zb(t1,proj))
        gsdk_c(i,j)=fac*(
     &   +gs_c(i,2,j)*zb(i4,t1)
     &   +gs_c(i,1,j)*mt*zb(i4,proj)/zb(t1,proj))
        gsdk_d(i,j)=fac*(
     &   +gs_d(i,2,j)*zb(i4,t1)
     &   +gs_d(i,1,j)*mt*zb(i4,proj)/zb(t1,proj))
      elseif (nwz == -1) then
        gsdk_a(i,j)=fac*(
     &   +gs_a(i,j,1)*za(i4,t1)
     &   +gs_a(i,j,2)*mt*za(i4,proj)/za(t1,proj))
        gsdk_b(i,j)=fac*(
     &   +gs_b(i,j,1)*za(i4,t1)
     &   +gs_b(i,j,2)*mt*za(i4,proj)/za(t1,proj))
        gsdk_c(i,j)=fac*(
     &   +gs_c(i,j,1)*za(i4,t1)
     &   +gs_c(i,j,2)*mt*za(i4,proj)/za(t1,proj))
        gsdk_d(i,j)=fac*(
     &   +gs_d(i,j,1)*za(i4,t1)
     &   +gs_d(i,j,2)*mt*za(i4,proj)/za(t1,proj))
      else
        write(6,*) 'nwz must be +1 or -1 in interdk_qq'
        stop
      endif

      enddo
      enddo

c----------------------------------------
c Square and sum all helicity amplitudes:
c----------------------------------------
      do i=1,4
         me(i)=0._dp
         intf(i)=0._dp
      enddo

         do j=1,2
            do i=1,2
c Diagonal color structure:
               me(1)=me(1)+
     &              abs(gsdk_a(i,j))**2*6._dp/wprop_a**2
               me(2)=me(2)+
     &              abs(gsdk_b(i,j))**2*6._dp/wprop_b**2
               me(3)=me(3)+
     &              abs(gsdk_c(i,j))**2*6._dp/wprop_c**2
               me(4)=me(4)+
     &              abs(gsdk_d(i,j))**2*6._dp/wprop_d**2
            enddo
c Off-diagonal color structure (Only add '1' helicity):
            intf(1)=intf(1)-real(
     &           gsdk_a(1,j)*conjg(gsdk_b(1,j))*4._dp)/wprop_a/wprop_b
            intf(2)=intf(2)-real(
     &           gsdk_b(1,j)*conjg(gsdk_c(1,j))*4._dp)/wprop_b/wprop_c
            intf(3)=intf(3)-real(
     &           gsdk_c(1,j)*conjg(gsdk_d(1,j))*4._dp)/wprop_c/wprop_d
            intf(4)=intf(4)-real(
     &           gsdk_d(1,j)*conjg(gsdk_a(1,j))*4._dp)/wprop_d/wprop_a
         enddo

c--- note that these are all corrections to the heavy quark line and
c--- thus should receive a factor of gsq_H
      gsq_H=fourpi*as_H

c Multiply by coupling constants:
      do i=1,4
         me(i)=me(i)*gsq_H**2*gwsq**2/4._dp
         intf(i)=intf(i)*gsq_H**2*gwsq**2/4._dp
      enddo

c 'Average' over incoming helicities and color:
      do i=1,4
         me(i)=me(i)/36._dp
         intf(i)=intf(i)/36._dp
      enddo

c--- ensure array pp is not contaminated upon return
      do k=1,4
      pp(t1,k)=zip
      pp(b1,k)=zip
      enddo

      return
      end

