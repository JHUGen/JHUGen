      subroutine qqb_zbbm(p,msq)
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
      include 'heavyflav.f'
      include 'nflav.f'
      integer:: j,k,nu,h1,h3,h5,h6,iflav
      real(dp):: p(mxpart,4),p12(mxpart,4),p21(mxpart,4),
     & msq(-nf:nf,-nf:nf),sumleptL,sumleptR,
     & facqq,facgg,p1Dp(5:6),p2Dp(5:6),lr(2),mQ
      complex(dp):: tamp,prop,coupL,coupR
      complex(dp):: qqb_a(2,2,2,2,2),qqb_b(2,2,2,2,2)
      complex(dp):: qbq_a(2,2,2,2,2),qbq_b(2,2,2,2,2)
      integer,parameter::swap(4)=(/2,1,3,4/)

      if     (flav == 6) then
        mQ=mt
        iflav=2
      elseif (flav == 5) then
        mQ=mb
        iflav=1
      else
        write(6,*) 'Invalid flav in qqb_zbbmas.f, flav=',flav
      endif

c--initialize to zero
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0._dp
      enddo
      enddo

      do j=5,6
      p1Dp(j)=p(1,4)*p(j,4)-p(1,1)*p(j,1)-p(1,2)*p(j,2)-p(1,3)*p(j,3)
      p2Dp(j)=p(2,4)*p(j,4)-p(2,1)*p(j,1)-p(2,2)*p(j,2)-p(2,3)*p(j,3)
      enddo

c---define modified (zero-mass) vectors
      do j=1,6
      do nu=1,4
      if (j<5) then
         p12(j,nu)=p(j,nu)
         p21(j,nu)=p(swap(j),nu)
      elseif (j==5) then
         p12(j,nu)=p(j,nu)-0.5_dp*mQ**2*p(2,nu)/p2Dp(5)
         p21(j,nu)=p(j,nu)-0.5_dp*mQ**2*p(1,nu)/p1Dp(5)
      elseif (j==6) then
         p12(j,nu)=p(j,nu)-0.5_dp*mQ**2*p(1,nu)/p1Dp(6)
         p21(j,nu)=p(j,nu)-0.5_dp*mQ**2*p(2,nu)/p2Dp(6)
      endif
      enddo
      enddo

C---Fill spinor products
      call dotem(6,p,s)
      prop=s(3,4)/cplx2((s(3,4)-zmass**2),zmass*zwidth)
      facgg=4._dp*V*gsq**2*esq**2*xn*avegg
      facqq=4._dp*V*gsq**2*esq**2*aveqq

      call spinoru(6,p12,za,zb)
      call mamps(mQ,1,2,3,4,5,6,qqb_a,qqb_b)

      coupL=Q(iflav)*q1+L(iflav)*l1*prop
      coupR=Q(iflav)*q1+R(iflav)*l1*prop
      call gamps0(mQ,1,2,3,4,5,6,sumleptL,coupL,coupR)
      coupL=Q(iflav)*q1+L(iflav)*r1*prop
      coupR=Q(iflav)*q1+R(iflav)*r1*prop
      call gamps0(mQ,1,2,4,3,5,6,sumleptR,coupL,coupR)

      call spinoru(6,p21,za,zb)
      call mamps(mQ,1,2,3,4,5,6,qbq_a,qbq_b)

      lr(1)=l1
      lr(2)=r1

      do j=-nflav,nflav
      k=-j
          if ((j == 0) .and. (k == 0)) then
            msq(j,k)=facgg*(sumleptL+sumleptR)
          elseif ((j > 0) .and. (k < 0)) then
            do h1=1,2
            do h3=1,2
            do h5=1,2
            do h6=1,2
            tamp=
     &          +(Q(j)*q1+L(j)*lr(h3)*prop)*qqb_a(1,h1,h3,h5,h6)
     &          +(Q(j)*q1+R(j)*lr(h3)*prop)*qqb_a(2,h1,h3,h5,h6)
     &          +(Q(iflav)*q1+L(iflav)*lr(h3)*prop)*qqb_b(1,h1,h3,h5,h6)
     &          +(Q(iflav)*q1+R(iflav)*lr(h3)*prop)*qqb_b(2,h1,h3,h5,h6)
            msq(j,k)=msq(j,k)+facqq*abs(tamp)**2
            enddo
            enddo
            enddo
            enddo


          elseif ((j < 0) .and. (k > 0)) then
            do h1=1,2
            do h3=1,2
            do h5=1,2
            do h6=1,2
            tamp=
     &          +(Q(k)*q1+L(k)*lr(h3)*prop)*qbq_a(1,h1,h3,h5,h6)
     &          +(Q(k)*q1+R(k)*lr(h3)*prop)*qbq_a(2,h1,h3,h5,h6)
     &          +(Q(iflav)*q1+L(iflav)*lr(h3)*prop)*qbq_b(1,h1,h3,h5,h6)
     &          +(Q(iflav)*q1+R(iflav)*lr(h3)*prop)*qbq_b(2,h1,h3,h5,h6)
            msq(j,k)=msq(j,k)+facqq*abs(tamp)**2
            enddo
            enddo
            enddo
            enddo


          endif
      enddo
      return
      end


