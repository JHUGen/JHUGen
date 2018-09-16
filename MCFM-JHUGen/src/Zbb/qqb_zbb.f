      subroutine qqb_zbb(p,msq)
      implicit none
      include 'types.f'
c---  Matrix elements squared
c     q(-p1)+qb(-p2) --> bbar(p5)+b(p6)+e^-(p3)+e^+(p4)
c---  averaged(summed) over initial(final) colours and spins

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      include 'zcouple.f'
      include 'ewcharge.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'msq_cs.f'
      include 'mmsq_cs.f'
      include 'heavyflav.f'
      include 'nflav.f'
      include 'first.f'
      integer:: j,k,nu,ics,j1,j2,j3
      real(dp):: p(mxpart,4),msq(-nf:nf,-nf:nf),mmsq(2,2),
     & pswap(mxpart,4),faclo,scalesq
      complex(dp):: tamp,prop
c      complex(dp):: qqb5,qbq5,qqb6,qbq6,qqb7,qbq7,qqb8,qbq8
      complex(dp):: qqb_a(2,2,2),qqb_b(2,2,2)
      complex(dp):: qbq_a(2,2,2),qbq_b(2,2,2)
      character qflav
      integer,parameter::swap(2)=(/2,1/)
      save scalesq

      if (first) then
       if     (flav == 5) then
         scalesq=mbsq
         qflav='b'
       elseif (flav == 4) then
         scalesq=mcsq
         qflav='c'
       else
         write(6,*) 'Invalid flav in qqb_zbb.f, flav=',flav
       endif
       write(6,*)
       write(6,*) '****************** Process info ********************'
       write(6,*) '*                                                  *'
       write(6,*) '* m'//qflav//
     &  '=0 for this process, although cuts are applied *'
       write(6,*) '* to simulate the effect of the '//qflav//
     &  ' mass:            *'
       write(6,*) '*                                                  *'
       write(6,99) ' *                pt('//qflav//
     &  ') > ',sqrt(scalesq),
     &  '                *'
       write(6,99) ' *                m('//qflav//qflav//
     &  ') > ',two*sqrt(scalesq),'                *'
       write(6,*) '****************************************************'
       first=.false.
      endif

c--initialize to zero
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0._dp
      enddo
      enddo

c ---Call the two gluon process which is defined in xzqqgg
C ---in the notation
C     0 ---> q(p1)+g(p2)+g(p3)+qbar(p4)+a(p5)+  l(p6)
C ---compared with ours which is:-
c     0 ---> b(p6)+g(p1)+g(p2)+bb(p5)+e^+(p4)+e^-(p3)

      do nu=1,4
      pswap(1,nu)=p(6,nu)
      pswap(2,nu)=p(1,nu)
      pswap(3,nu)=p(2,nu)
      pswap(4,nu)=p(5,nu)
      pswap(5,nu)=p(4,nu)
      pswap(6,nu)=p(3,nu)
      enddo
      call spinoru(6,pswap,za,zb)
      call xzqqgg(mmsq)

C---Fill spinor products
      call spinoru(6,p,za,zb)
      prop=s(3,4)/cplx2((s(3,4)-zmass**2),zmass*zwidth)

c ensure that we have a hard process
      if (  (s(5,6) < four*scalesq)
     & .or. (s(1,5)*s(2,5)/s(1,2) < scalesq)
     & .or. (s(1,6)*s(2,6)/s(1,2) < scalesq) ) return

c--- qqb
      call ampqqb_qqb(1,2,6,5,qqb_a,qqb_b)
C Instead of calling ampqqb_qqb(2,1,5,6,qbq_a,qbq_b)
c--- qbq from symmetries
      do j1=1,2
      do j2=1,2
      do j3=1,2
      qbq_a(j1,j2,j3)=-qqb_a(swap(j1),j2,j3)
      qbq_b(j1,j2,j3)=+qqb_b(swap(j1),j2,j3)
      enddo
      enddo
      enddo

      faclo=4._dp*V*gsq**2*esq**2*aveqq

      do j=-nflav,nflav
      k=-j
          if ((j == 0) .and. (k == 0)) then
            msq(j,k)=
     &      +abs(Q(flav)*q1+L(flav)*l1*prop)**2*mmsq(1,1)
     &      +abs(Q(flav)*q1+R(flav)*l1*prop)**2*mmsq(2,1)
     &      +abs(Q(flav)*q1+L(flav)*r1*prop)**2*mmsq(1,2)
     &      +abs(Q(flav)*q1+R(flav)*r1*prop)**2*mmsq(2,2)
            do ics=0,2
            msq_cs(ics,j,k)=
     &      +abs(Q(flav)*q1+L(flav)*l1*prop)**2*mmsq_cs(ics,1,1)
     &      +abs(Q(flav)*q1+R(flav)*l1*prop)**2*mmsq_cs(ics,2,1)
     &      +abs(Q(flav)*q1+L(flav)*r1*prop)**2*mmsq_cs(ics,1,2)
     &      +abs(Q(flav)*q1+R(flav)*r1*prop)**2*mmsq_cs(ics,2,2)
            enddo
          elseif ((j > 0) .and. (k < 0)) then
            tamp=(Q(j)*q1+L(j)*l1*prop)*qqb_a(1,1,1)
     &          +(Q(flav)*q1+L(flav)*l1*prop)*qqb_b(1,1,1)
            msq(j,k)=faclo*abs(tamp)**2
            tamp=(Q(j)*q1+L(j)*l1*prop)*qqb_a(1,2,1)
     &          +(Q(flav)*q1+R(flav)*l1*prop)*qqb_b(1,2,1)
            msq(j,k)=msq(j,k)+faclo*abs(tamp)**2
            tamp=(Q(j)*q1+L(j)*r1*prop)*qqb_a(1,1,2)
     &          +(Q(flav)*q1+L(flav)*r1*prop)*qqb_b(1,1,2)
            msq(j,k)=msq(j,k)+faclo*abs(tamp)**2
            tamp=(Q(j)*q1+L(j)*r1*prop)*qqb_a(1,2,2)
     &          +(Q(flav)*q1+R(flav)*r1*prop)*qqb_b(1,2,2)
            msq(j,k)=msq(j,k)+faclo*abs(tamp)**2
            tamp=(Q(j)*q1+R(j)*l1*prop)*qqb_a(2,1,1)
     &          +(Q(flav)*q1+L(flav)*l1*prop)*qqb_b(2,1,1)
            msq(j,k)=msq(j,k)+faclo*abs(tamp)**2
            tamp=(Q(j)*q1+R(j)*l1*prop)*qqb_a(2,2,1)
     &          +(Q(flav)*q1+R(flav)*l1*prop)*qqb_b(2,2,1)
            msq(j,k)=msq(j,k)+faclo*abs(tamp)**2
            tamp=(Q(j)*q1+R(j)*r1*prop)*qqb_a(2,1,2)
     &          +(Q(flav)*q1+L(flav)*r1*prop)*qqb_b(2,1,2)
            msq(j,k)=msq(j,k)+faclo*abs(tamp)**2
            tamp=(Q(j)*q1+R(j)*r1*prop)*qqb_a(2,2,2)
     &          +(Q(flav)*q1+R(flav)*r1*prop)*qqb_b(2,2,2)
            msq(j,k)=msq(j,k)+faclo*abs(tamp)**2
          elseif ((j < 0) .and. (k > 0)) then
            tamp=(Q(k)*q1+L(k)*l1*prop)*qbq_a(1,1,1)
     &          +(Q(flav)*q1+L(flav)*l1*prop)*qbq_b(1,1,1)
            msq(j,k)=faclo*abs(tamp)**2
            tamp=(Q(k)*q1+L(k)*l1*prop)*qbq_a(1,2,1)
     &          +(Q(flav)*q1+R(flav)*l1*prop)*qbq_b(1,2,1)
            msq(j,k)=msq(j,k)+faclo*abs(tamp)**2
            tamp=(Q(k)*q1+L(k)*r1*prop)*qbq_a(1,1,2)
     &          +(Q(flav)*q1+L(flav)*r1*prop)*qbq_b(1,1,2)
            msq(j,k)=msq(j,k)+faclo*abs(tamp)**2
            tamp=(Q(k)*q1+L(k)*r1*prop)*qbq_a(1,2,2)
     &          +(Q(flav)*q1+R(flav)*r1*prop)*qbq_b(1,2,2)
            msq(j,k)=msq(j,k)+faclo*abs(tamp)**2
            tamp=(Q(k)*q1+R(k)*l1*prop)*qbq_a(2,1,1)
     &          +(Q(flav)*q1+L(flav)*l1*prop)*qbq_b(2,1,1)
            msq(j,k)=msq(j,k)+faclo*abs(tamp)**2
            tamp=(Q(k)*q1+R(k)*l1*prop)*qbq_a(2,2,1)
     &          +(Q(flav)*q1+R(flav)*l1*prop)*qbq_b(2,2,1)
            msq(j,k)=msq(j,k)+faclo*abs(tamp)**2
            tamp=(Q(k)*q1+R(k)*r1*prop)*qbq_a(2,1,2)
     &          +(Q(flav)*q1+L(flav)*r1*prop)*qbq_b(2,1,2)
            msq(j,k)=msq(j,k)+faclo*abs(tamp)**2
            tamp=(Q(k)*q1+R(k)*r1*prop)*qbq_a(2,2,2)
     &          +(Q(flav)*q1+R(flav)*r1*prop)*qbq_b(2,2,2)
            msq(j,k)=msq(j,k)+faclo*abs(tamp)**2
          endif
      enddo

      return

   99 format(a26,f6.3,a21)

      end





