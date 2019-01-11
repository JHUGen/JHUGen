      subroutine qqb_zccm(p,msq)
c---  Matrix elements squared
c     q(-p1)+qb(-p2) --> cbar(p5)+c(p6)+e^-(p3)+e^+(p4)
c---  averaged(summed) over initial(final) colours and spins
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      include 'zcouple.f'
      include 'ewcharge.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      integer j,k,nu,h1,h3,h5,h6
      double precision p(mxpart,4),p12(mxpart,4),p21(mxpart,4),
     . msq(-nf:nf,-nf:nf),sumleptL,sumleptR,
     . facqq,facgg,p1Dp(5:6),p2Dp(5:6),lr(2)
      double complex tamp,prop,coupL,coupR
      double complex qqb_a(2,2,2,2,2),qqb_b(2,2,2,2,2)
      double complex qbq_a(2,2,2,2,2),qbq_b(2,2,2,2,2)
      integer,parameter::swap(4)=(/2,1,3,4/)

c--initialize to zero
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      enddo
      enddo

      do j=5,6
      p1Dp(j)=p(1,4)*p(j,4)-p(1,1)*p(j,1)-p(1,2)*p(j,2)-p(1,3)*p(j,3)
      p2Dp(j)=p(2,4)*p(j,4)-p(2,1)*p(j,1)-p(2,2)*p(j,2)-p(2,3)*p(j,3)
      enddo

c---define modified (zero-mass) vectors
      do j=1,6
      do nu=1,4
      if (j.lt.5) then
         p12(j,nu)=p(j,nu)
         p21(j,nu)=p(swap(j),nu)
      elseif (j.eq.5) then
         p12(j,nu)=p(j,nu)-0.5d0*mc**2*p(2,nu)/p2Dp(5)
         p21(j,nu)=p(j,nu)-0.5d0*mc**2*p(1,nu)/p1Dp(5)
      elseif (j.eq.6) then
         p12(j,nu)=p(j,nu)-0.5d0*mc**2*p(1,nu)/p1Dp(6)
         p21(j,nu)=p(j,nu)-0.5d0*mc**2*p(2,nu)/p2Dp(6)
      endif
      enddo
      enddo

C---Fill spinor products
      call dotem(6,p,s)
      prop=s(3,4)/dcmplx((s(3,4)-zmass**2),zmass*zwidth)
      facgg=4d0*V*gsq**2*esq**2*xn*avegg
      facqq=4d0*V*gsq**2*esq**2*aveqq

      call spinoru(6,p12,za,zb)
      call mamps(mc,1,2,3,4,5,6,qqb_a,qqb_b)

      coupL=Q(2)*q1+L(2)*l1*prop
      coupR=Q(2)*q1+R(2)*l1*prop
      call gamps0(mc,1,2,3,4,5,6,sumleptL,coupL,coupR)
      coupL=Q(2)*q1+L(2)*r1*prop
      coupR=Q(2)*q1+R(2)*r1*prop
      call gamps0(mc,1,2,4,3,5,6,sumleptR,coupL,coupR)

      call spinoru(6,p21,za,zb)
      call mamps(mc,1,2,3,4,5,6,qbq_a,qbq_b)

      lr(1)=l1
      lr(2)=r1

      do j=-(nf-1),(nf-1)
      k=-j
          if ((j .eq. 0) .and. (k .eq. 0)) then
            msq(j,k)=facgg*(sumleptL+sumleptR)
          elseif ((j .gt. 0) .and. (k .lt. 0)) then
            do h1=1,2
            do h3=1,2
            do h5=1,2
            do h6=1,2
            tamp=
     .          +(Q(j)*q1+L(j)*lr(h3)*prop)*qqb_a(1,h1,h3,h5,h6)
     .          +(Q(j)*q1+R(j)*lr(h3)*prop)*qqb_a(2,h1,h3,h5,h6)
     .          +(Q(2)*q1+L(2)*lr(h3)*prop)*qqb_b(1,h1,h3,h5,h6)
     .          +(Q(2)*q1+R(2)*lr(h3)*prop)*qqb_b(2,h1,h3,h5,h6)
            msq(j,k)=msq(j,k)+facqq*abs(tamp)**2
            enddo
            enddo
            enddo
            enddo


          elseif ((j .lt. 0) .and. (k .gt. 0)) then
            do h1=1,2
            do h3=1,2
            do h5=1,2
            do h6=1,2
            tamp=
     .          +(Q(k)*q1+L(k)*lr(h3)*prop)*qbq_a(1,h1,h3,h5,h6)
     .          +(Q(k)*q1+R(k)*lr(h3)*prop)*qbq_a(2,h1,h3,h5,h6)
     .          +(Q(2)*q1+L(2)*lr(h3)*prop)*qbq_b(1,h1,h3,h5,h6)
     .          +(Q(2)*q1+R(2)*lr(h3)*prop)*qbq_b(2,h1,h3,h5,h6)
            msq(j,k)=msq(j,k)+facqq*abs(tamp)**2
            enddo
            enddo
            enddo
            enddo


          endif
      enddo
      return
      end


