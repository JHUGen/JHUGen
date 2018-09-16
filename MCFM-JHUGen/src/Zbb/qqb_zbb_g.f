      subroutine qqb_zbb_g(p,msq)
      implicit none
      include 'types.f'
************************************************************************
*     Author: R.K. Ellis                                               *
*     July, 1998.                                                      *
************************************************************************
c---Matrix element squared averaged over initial colors and spins
c     p(-p1)+p(-p2) -->  g*  + Z + p(p7)
c                        |     |
c                        |     --> e-(p3)+e^+(p4)
c                        |
c                         ---> bb(p5)+b(p6)

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      include 'zcouple.f'
      include 'ewcouple.f'
      include 'ewcharge.f'
      include 'noglue.f'
      include 'heavyflav.f'
      include 'nflav.f'
      include 'first.f'
      integer:: j,k,hq,Qh,hg,lh
      real(dp):: P(mxpart,4),msq(-nf:nf,-nf:nf),mmsq(2,2)
      real(dp):: fac,LRq(2),LRb(2),lr1(2),scalesq
      real(dp):: msq_qqb(2,2,2,2,4),msq_qbq(2,2,2,2,4),
     &                 msq_qg(2,2,2,2,4),msq_qbg(2,2,2,2,4),
     &                 msq_gqb(2,2,2,2,4),msq_gq(2,2,2,2,4)
      complex(dp):: prop,czq,czb
      save scalesq

      if (first) then
       if     (flav == 5) then
         scalesq=mbsq
       elseif (flav == 4) then
         scalesq=mcsq
       else
         write(6,*) 'Invalid flav in qqb_zbb_g.f, flav=',flav
       endif
       first=.false.
      endif

      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0._dp
      enddo
      enddo

C---call spinor routine and load common block twopij
      call spinoru(7,p,za,zb)
      prop=s(3,4)/cplx2((s(3,4)-zmass**2),zmass*zwidth)

C so our notation is
C     0 ---> q(p6)+g(p1)+g(p2)+g(p7)+qbar(p5)+a(p4)+l(p3)
C whilst xzqqggg notation is
C     0 ---> q(p1)+g(p2)+g(p3)+g(p4)+qbar(p5)+a(p6)+l(p7)
      if (gqonly) then
        do j=1,2
        do k=1,2
          mmsq(j,k)=0._dp
        enddo
        enddo
      else
        call xzqqggg(6,1,2,7,5,4,3,mmsq)
      endif

      if (
     &      (s(5,6) < four*scalesq)
     & .or. (s(1,5)*s(2,5)/s(1,2) < scalesq)
     & .or. (s(1,6)*s(2,6)/s(1,2) < scalesq) ) return

c--- note that (2,1) and (4,3) are switched due to crossing from NT
      if (gqonly) then
        do hq=1,2
        do Qh=1,2
        do hg=1,2
        do lh=1,2
        do j=1,4
          msq_qqb(hq,Qh,hg,lh,j)=0._dp
          msq_qbq(hq,Qh,hg,lh,j)=0._dp
        enddo
        enddo
        enddo
        enddo
        enddo
      else
        call msq_qqQQg(2,1,6,5,7,4,3,msq_qqb)
        call msq_qqQQg(1,2,6,5,7,4,3,msq_qbq)
      endif
      call msq_qqQQg(7,1,6,5,2,4,3,msq_qg)
      call msq_qqQQg(2,7,6,5,1,4,3,msq_gqb)
      call msq_qqQQg(7,2,6,5,1,4,3,msq_gq)
      call msq_qqQQg(1,7,6,5,2,4,3,msq_qbg)

c--- note the factor of 4._dp*xw**2 relative to wbb
      fac=4._dp*gsq**3*esq**2
c--- extra factor of 2**3=8 to compensate for Ta normalization
      fac=fac*8._dp

      LRb(1)=L(flav)
      LRb(2)=R(flav)

      do j=-nflav,nflav
      do k=-nflav,nflav
      if( j .ne. 0 .and. k .ne. 0 .and. j .ne. -k) goto 19

      msq(j,k)=0._dp

      if     ((j == 0) .and. (k == 0)) then
C---L(flav),R(flav) for coupling to down quark
      msq(j,k)=abs(Q(flav)*q1+prop*L(flav)*l1)**2*mmsq(1,1)
     &        +abs(Q(flav)*q1+prop*R(flav)*l1)**2*mmsq(2,1)
     &        +abs(Q(flav)*q1+prop*L(flav)*r1)**2*mmsq(1,2)
     &        +abs(Q(flav)*q1+prop*R(flav)*r1)**2*mmsq(2,2)
      elseif ((j > 0) .and. (k < 0)) then
        LRq(1)=L(j)
        LRq(2)=R(j)
        lr1(1)=l1
        lr1(2)=r1
        do hq=1,2
        do Qh=1,2
        do hg=1,2
        do lh=1,2
c--- couplings of Z to q (12) and Z to b (56)
c--- the 2nd of these is conjugated for convenience in interference
          czq=Q(j)*q1+prop*LRq(hq)*lr1(lh)
          czb=Q(flav)*q1+conjg(prop)*LRb(Qh)*lr1(lh)
          msq(j,k)=msq(j,k)+aveqq*fac*(
     &      abs(czq)**2*msq_qqb(hq,Qh,hg,lh,1)
     &     +abs(czb)**2*msq_qqb(hq,Qh,hg,lh,2)
     &     +real(czq*czb*cplx2(msq_qqb(hq,Qh,hg,lh,3),
     &                           msq_qqb(hq,Qh,hg,lh,4))))
        enddo
        enddo
        enddo
        enddo
      elseif ((j < 0) .and. (k > 0)) then
        LRq(1)=L(k)
        LRq(2)=R(k)
        lr1(1)=l1
        lr1(2)=r1
        do hq=1,2
        do Qh=1,2
        do hg=1,2
        do lh=1,2
          czq=Q(k)*q1+prop*LRq(hq)*lr1(lh)
          czb=Q(flav)*q1+conjg(prop)*LRb(Qh)*lr1(lh)
          msq(j,k)=msq(j,k)+aveqq*fac*(
     &      abs(czq)**2*msq_qbq(hq,Qh,hg,lh,1)
     &     +abs(czb)**2*msq_qbq(hq,Qh,hg,lh,2)
     &     +real(czq*czb*cplx2(msq_qbq(hq,Qh,hg,lh,3),
     &                           msq_qbq(hq,Qh,hg,lh,4))))
        enddo
        enddo
        enddo
        enddo
      elseif ((j > 0) .and. (k == 0)) then
        LRq(1)=L(j)
        LRq(2)=R(j)
        lr1(1)=l1
        lr1(2)=r1
        do hq=1,2
        do Qh=1,2
        do hg=1,2
        do lh=1,2
          czq=Q(j)*q1+prop*LRq(hq)*lr1(lh)
          czb=Q(flav)*q1+conjg(prop)*LRb(Qh)*lr1(lh)
          msq(j,k)=msq(j,k)+aveqg*fac*(
     &      abs(czq)**2*msq_qg(hq,Qh,hg,lh,1)
     &     +abs(czb)**2*msq_qg(hq,Qh,hg,lh,2)
     &     +real(czq*czb*cplx2(msq_qg(hq,Qh,hg,lh,3),
     &                           msq_qg(hq,Qh,hg,lh,4))))
        enddo
        enddo
        enddo
        enddo
      elseif ((j < 0) .and. (k == 0)) then
        LRq(1)=L(-j)
        LRq(2)=R(-j)
        lr1(1)=l1
        lr1(2)=r1
        do hq=1,2
        do Qh=1,2
        do hg=1,2
        do lh=1,2
          czq=Q(-j)*q1+prop*LRq(hq)*lr1(lh)
          czb=Q(flav)*q1+conjg(prop)*LRb(Qh)*lr1(lh)
          msq(j,k)=msq(j,k)+aveqg*fac*(
     &      abs(czq)**2*msq_qbg(hq,Qh,hg,lh,1)
     &     +abs(czb)**2*msq_qbg(hq,Qh,hg,lh,2)
     &     +real(czq*czb*cplx2(msq_qbg(hq,Qh,hg,lh,3),
     &                           msq_qbg(hq,Qh,hg,lh,4))))
        enddo
        enddo
        enddo
        enddo
      elseif ((j == 0) .and. (k > 0)) then
        LRq(1)=L(k)
        LRq(2)=R(k)
        lr1(1)=l1
        lr1(2)=r1
        do hq=1,2
        do Qh=1,2
        do hg=1,2
        do lh=1,2
          czq=Q(k)*q1+prop*LRq(hq)*lr1(lh)
          czb=Q(flav)*q1+conjg(prop)*LRb(Qh)*lr1(lh)
          msq(j,k)=msq(j,k)+aveqg*fac*(
     &      abs(czq)**2*msq_gq(hq,Qh,hg,lh,1)
     &     +abs(czb)**2*msq_gq(hq,Qh,hg,lh,2)
     &     +real(czq*czb*cplx2(msq_gq(hq,Qh,hg,lh,3),
     &                           msq_gq(hq,Qh,hg,lh,4))))
        enddo
        enddo
        enddo
        enddo
      elseif ((j == 0) .and. (k < 0)) then
        LRq(1)=L(-k)
        LRq(2)=R(-k)
        lr1(1)=l1
        lr1(2)=r1
        do hq=1,2
        do Qh=1,2
        do hg=1,2
        do lh=1,2
          czq=Q(-k)*q1+prop*LRq(hq)*lr1(lh)
          czb=Q(flav)*q1+conjg(prop)*LRb(Qh)*lr1(lh)
          msq(j,k)=msq(j,k)+aveqg*fac*(
     &      abs(czq)**2*msq_gqb(hq,Qh,hg,lh,1)
     &     +abs(czb)**2*msq_gqb(hq,Qh,hg,lh,2)
     &     +real(czq*czb*cplx2(msq_gqb(hq,Qh,hg,lh,3),
     &                           msq_gqb(hq,Qh,hg,lh,4))))
        enddo
        enddo
        enddo
        enddo
      endif

   19 continue
      enddo
      enddo
      return
      end


