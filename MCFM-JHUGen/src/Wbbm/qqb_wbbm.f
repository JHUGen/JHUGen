      subroutine qqb_wbbm(p,msq)
      implicit none
      include 'types.f'
c--- Note: optimization in the "sumsq" function by RKE, April 2009.
c---  Matrix elements squared
c     q(-p1)+qb(-p2) --> nu(p3)+e^+(p4)+b(p5)+bb(p6)
c---  averaged(summed) over initial(final) colours and spins

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ckm.f'
      include 'masses.f'
      include 'sprods_com.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'heavyflav.f'
      integer:: j,k
      real(dp):: p(mxpart,4),msq(-nf:nf,-nf:nf)
      real(dp):: qqb,qbq,sumsq
      real(dp):: faclo,mQsq
      faclo=V*gsq**2*gwsq**2*aveqq

C--- Set up the correct mass, according to 'flav'
      if     (flav == 6) then
        mQsq=mt**2
      elseif (flav == 5) then
        mQsq=mb**2
      elseif (flav == 4) then
        mQsq=mc**2
      else
        write(6,*) 'Wrong flavour in qqb_wbbm.f: flav=',flav
        stop
      endif

C---Initialize whole array to zero
      msq(:,:)=zero

C---Fill dot-products
      call dotem(6,p,s)

      qqb=faclo*sumsq(1,2,3,4,6,5,mQsq)
      qbq=faclo*sumsq(2,1,3,4,6,5,mQsq)

      do j=-(flav-1),(flav-1)
      do k=-(flav-1),(flav-1)
      if     ((j > 0) .and. (k < 0)) then
               msq(j,k)=Vsq(j,k)*qqb
      elseif ((j < 0) .and. (k > 0)) then
               msq(j,k)=Vsq(j,k)*qbq
      endif
      enddo
      enddo

      return
      end

      function sumsq(p1,p2,p3,p4,p5,p6,mQsq)
      implicit none
      include 'types.f'
      real(dp):: sumsq

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'sprods_com.f'
      include 'masses.f'
      integer:: p1,p2,p3,p4,p5,p6
      real(dp):: s56,s134,s234,prop,mQsq
      s134=s(p1,p3)+s(p1,p4)+s(p3,p4)
      s234=s(p2,p3)+s(p2,p4)+s(p3,p4)
      s56=s134+s234+s(p1,p2)-s(p3,p4)
c---overall factor of 4 removed
      prop=s56**2*((s(p3,p4)-wmass**2)**2+(wmass*wwidth)**2)

      sumsq=one/(s134*s234)
     & *(+s(p1,p4)*s(p2,p4)
     & *(s(p2,p6)*s(p3,p5)+s(p2,p5)*s(p3,p6)+two*s(p3,p5)*s(p3,p6))
     & +s(p1,p3)*s(p2,p3)
     & *(s(p1,p6)*s(p4,p5)+s(p1,p5)*s(p4,p6)+two*s(p4,p5)*s(p4,p6))
     & -s(p2,p3)*s(p3,p4)
     & *(s(p1,p6)*s(p4,p5)+s(p1,p5)*s(p4,p6)+two*s(p1,p5)*s(p1,p6))
     & -s(p1,p4)*s(p3,p4)
     & *(s(p2,p6)*s(p3,p5)+s(p2,p5)*s(p3,p6)+two*s(p2,p5)*s(p2,p6)))

      sumsq=sumsq+((s(p1,p2)*s(p3,p4)-s(p1,p3)*s(p2,p4))
     &   *(s(p3,p6)*s(p4,p5)+s(p3,p5)*s(p4,p6)
     &    +s(p2,p6)*s(p4,p5)+s(p2,p5)*s(p4,p6)
     &    +s(p1,p6)*s(p3,p5)+s(p1,p6)*s(p2,p5)
     &    +s(p1,p5)*s(p3,p6)+s(p1,p5)*s(p2,p6))
     &  -s(p1,p4)*s(p2,p3)*(+s(p3,p6)*s(p4,p5)+s(p3,p5)*s(p4,p6)
     &  +s(p1,p6)*s(p2,p5)+s(p1,p5)*s(p2,p6)-two*s(p3,p4)*s(p5,p6))
     & )/(s134*s234)

      sumsq = sumsq + two*s134**(-1)*s234**(-1)*mQsq * (
     &     +(s(p1,p2)*s(p3,p4)-s(p1,p3)*s(p2,p4))
     *     *(s(p1,p2)+s(p3,p4)+s(p1,p3)+s(p2,p4))
     &     +s(p1,p4)*s(p2,p3)*(s(p1,p3)-s(p1,p2)+s(p2,p4)+s(p3,p4)))

      sumsq = sumsq+s(p1,p4)/s134**2 * (
     &     +(s(p1,p3)+s(p3,p4))
     &     *(s(p2,p5)*(s(p1,p6)+s(p4,p6))+s(p2,p6)*(s(p1,p5)+s(p4,p5)))
     &     - s(p1,p4)*(s(p2,p5)*s(p3,p6)+s(p2,p6)*s(p3,p5))
     &     +two*mQsq*(
     &     (s(p1,p2)+s(p2,p4))*(s(p1,p3)+s(p3,p4))-s(p1,p4)*s(p2,p3)))
      sumsq = sumsq + s(p2,p3)/s234**2 * (
     &     +(s(p2,p4)+s(p3,p4))
     &     *(s(p1,p5)*(s(p2,p6)+s(p3,p6))+s(p1,p6)*(s(p2,p5)+s(p3,p5)))
     &     -s(p2,p3)*(s(p1,p5)*s(p4,p6)+s(p1,p6)*s(p4,p5))
     &     +two*mQsq*(
     &     (s(p1,p2)+s(p1,p3))*(s(p2,p4)+s(p3,p4))-s(p1,p4)*s(p2,p3)))

      sumsq=sumsq/prop
      return
      end
