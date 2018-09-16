      subroutine qqb_w2jet_soft(P,msq)
      implicit none
      include 'types.f'
c---Matrix element SUBTRACTION squared averaged over initial colors and spins
c     q(-p1)+qbar(-p2) -->  b(p5)+bb(p6) + Z +g(p7)
c                                          |
c                                          --> mu(p3)+nubar(p4)
c                            
c   positively charged W only
c--all momenta incoming
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'qcdcouple.f'
      include 'masses.f'
      include 'flags.f'
      include 'msq_cs.f'
      include 'jetlabel.f'
      include 'nqcdjets.f'
      integer:: i,j,k,n,qq
      integer:: perm(5:7,5:7)
      real(dp):: P(mxpart,4),msq(-nf:nf,-nf:nf),pjet(mxpart,4),
     & msq0(-nf:nf,-nf:nf),s(mxpart,mxpart),Rcut
      real(dp):: eik56,eik15,eik16,eik26,eik25,eik12
      real(dp):: eik17_2,eik17_5,eik17_6,eik27_5,eik27_6
      real(dp):: eik27_1,eik57_1,eik67_1,eik57_2,eik67_2
      real(dp):: eik57_6,eik67_5,eik65_1,eik76_5,eik65_2
      real(dp):: eik15_6,eik56_1,eik56_7,eik25_6,eik56_2
      real(dp):: eik1a_b(6),eikba_1(6),eikab_c(6),
     &                 eikcb_a(6),eikbc_2(6),eik2c_b(6),
     &                 eik1b_2(6),eik2b_1(6),eikf,eikc,
     &                 msq_ac(6,0:2,-nf:nf,-nf:nf),p_ac(mxpart,4)
      common/rcut/rcut
      real(dp):: mqq(0:2,-nf:nf,-nf:nf)
      real(dp):: 
     & sub17_2(4),sub27_1(4),sub57_6(4),sub67_5(4),
     & sub17_5(4),sub57_1(4),sub27_5(4),sub57_2(4),
     & sub17_6(4),sub67_1(4),sub27_6(4),sub67_2(4)
      parameter(qq=1)
      common/mqq/mqq
      integer,parameter:: a(6)=(/5,5,7,6,6,7/),
     & b(6)=(/6,7,5,7,5,6/),c(6)=(/7,6,6,5,7,5/),
     & perm(5:7,5:7)=reshape((/0,4,6,2,0,3,1,5,0/),(/3,3/))
!$omp threadprivate/mqq/)
      
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0._dp
      enddo
      enddo

      if (Gflag) then
      do j=1,4
      do k=1,4
      p_ac(j,k)=p(j,k)
      enddo
      enddo
      call dotem(7,p,s)

      do n=1,6
      eik1a_b(n)=+s(1,b(n))/(-s(1,a(n))+s(a(n),b(n)))/s(1,a(n))*two*gsq
      eikba_1(n)=-s(b(n),1)/(s(b(n),a(n))-s(a(n),1))
     &             /s(b(n),a(n))*two*gsq
      eikab_c(n)=+s(a(n),c(n))/(s(a(n),b(n))+s(b(n),c(n)))
     &             /s(a(n),b(n))*two*gsq
      eikcb_a(n)=+s(c(n),a(n))/(s(c(n),b(n))+s(b(n),a(n)))
     &             /s(c(n),b(n))*two*gsq
      eikbc_2(n)=-s(b(n),2)/(s(b(n),c(n))-s(c(n),2))
     &             /s(b(n),c(n))*two*gsq
      eik2c_b(n)=+s(2,b(n))/(-s(2,c(n))+s(c(n),b(n)))/s(2,c(n))*two*gsq
      eik1b_2(n)=+s(1,2)/(s(1,b(n))+s(2,b(n)))/s(1,b(n))*two*gsq
      eik2b_1(n)=+s(1,2)/(s(1,b(n))+s(2,b(n)))/s(2,b(n))*two*gsq
      do k=1,4
        p_ac(5,k)=p(a(n),k)
        p_ac(6,k)=p(c(n),k)
        p_ac(7,k)=0._dp
      enddo
      call qqb_w2jet(p_ac,msq0)
      call genclust2(p_ac,rcut,pjet,0)
      do i=0,2
      do j=-nf,nf
      do k=-nf,nf
        if (jets .ne. nqcdjets) then
          msq_ac(n,i,j,k)=0._dp
        else
          msq_ac(n,i,j,k)=msq_cs(i,j,k)
        endif
      enddo
      enddo
      enddo
      
      enddo
      
      eik17_2=+s(1,2)/(s(1,7)+s(2,7))/s(1,7)*two*gsq
      eik27_1=+s(1,2)/(s(1,7)+s(2,7))/s(2,7)*two*gsq

      eik57_6=+s(5,6)/(s(5,7)+s(6,7))/s(5,7)*two*gsq
      eik67_5=+s(5,6)/(s(5,7)+s(6,7))/s(6,7)*two*gsq

      eik17_6=+s(1,6)/(-s(1,7)+s(6,7))/s(1,7)*two*gsq
      eik67_1=-s(1,6)/(-s(1,7)+s(6,7))/s(6,7)*two*gsq

      eik27_5=+s(2,5)/(-s(2,7)+s(5,7))/s(2,7)*two*gsq
      eik57_2=-s(2,5)/(-s(2,7)+s(5,7))/s(5,7)*two*gsq

      eik27_6=+s(2,6)/(-s(2,7)+s(6,7))/s(2,7)*two*gsq
      eik67_2=-s(2,6)/(-s(2,7)+s(6,7))/s(6,7)*two*gsq

      eik17_5=+s(1,5)/(-s(1,7)+s(5,7))/s(1,7)*two*gsq
      eik57_1=-s(1,5)/(-s(1,7)+s(5,7))/s(5,7)*two*gsq

      eik57_6=+s(6,5)/(s(6,7)+s(5,7))/s(5,7)*two*gsq
      eik67_5=+s(5,6)/(s(5,7)+s(6,7))/s(6,7)*two*gsq

c---new
      eik15_6=+s(1,6)/(-s(1,5)+s(6,5))/s(1,5)*two*gsq

      eik56_1=-s(1,5)/(-s(1,6)+s(5,6))/s(5,6)*two*gsq
      eik65_1=-s(1,6)/(-s(1,5)+s(5,6))/s(5,6)*two*gsq
      eik65_2=-s(2,6)/(-s(2,5)+s(5,6))/s(5,6)*two*gsq

      eik56_7=+s(7,5)/(+s(7,6)+s(5,6))/s(5,6)*two*gsq
      eik76_5=+s(7,5)/(+s(5,6)+s(7,6))/s(7,6)*two*gsq

      eik25_6=+s(2,6)/(-s(2,5)+s(6,5))/s(2,5)*two*gsq
      eik56_2=-s(2,5)/(-s(2,6)+s(5,6))/s(5,6)*two*gsq

      call qqb_w2jet(p,msq0)

      do j=-nf,nf
      do k=-nf,nf
          if ((j > 0) .and. (k<0)) then
          msq(j,k)=0._dp
          do n=1,6
          msq(j,k)=msq(j,k)+(eik1a_b(n)+eikba_1(n))
     &               *msq_ac(perm(b(n),c(n)),1,j,k)*xn/3._dp
     &                     +(eikab_c(n)+eikcb_a(n))
     &               *msq_ac(perm(a(n),c(n)),1,j,k)*xn/3._dp
     &                     +(eikbc_2(n)+eik2c_b(n))
     &               *msq_ac(perm(a(n),b(n)),1,j,k)*xn/3._dp                   
          enddo
c--- done with leading colour
          elseif ((j < 0) .and. (k>0)) then
          msq(j,k)=+(eik57_2+eik27_5)*msq_cs(1,j,k)*xn
     &             +(eik67_1+eik17_6)*msq_cs(1,j,k)*xn
     &             +(eik57_6+eik67_5)*msq_cs(1,j,k)*xn
     &             +(eik67_2+eik27_6)*msq_cs(2,j,k)*xn
     &             +(eik57_1+eik17_5)*msq_cs(2,j,k)*xn
     &             +(eik57_6+eik67_5)*msq_cs(2,j,k)*xn
          msq(j,k)=+(eik25_6+eik65_2)*msq_cs(1,j,k)*xn/3._dp
     &             +(eik56_7+eik76_5)*msq_cs(1,j,k)*xn/3._dp
     &             +(eik67_1+eik17_6)*msq_cs(1,j,k)*xn/3._dp
c--- done with leading colour
           elseif ((j == 0) .and. (k == 0)) then
          msq(j,k)=+(eik17_5+eik57_1)*msq_cs(1,j,k)*xn
     &             +(eik27_6+eik67_2)*msq_cs(1,j,k)*xn
     &             +(eik17_2+eik27_1)*msq_cs(1,j,k)*xn
     &             +(eik27_5+eik57_2)*msq_cs(2,j,k)*xn
     &             +(eik17_6+eik67_1)*msq_cs(2,j,k)*xn
     &             +(eik17_2+eik27_1)*msq_cs(2,j,k)*xn
c--- done with leading colour
c     &             -(eik57_6+eik67_5)*msq_cs(1,j,k)/xn
c     &             -(eik57_6+eik67_5)*msq_cs(2,j,k)/xn
c     &             +(eik57_1+eik17_5)*msq_cs(0,j,k)*xn
c     &             +(eik67_2+eik27_6)*msq_cs(0,j,k)*xn
c     &             +(eik67_1+eik17_6)*msq_cs(0,j,k)*xn
c     &             +(eik57_2+eik27_5)*msq_cs(0,j,k)*xn
c--- done with one gluon QE.e-_dplike (1/N**2 colour-suppressed)
c     &             -(eik57_6+eik67_5)*msq_cs(0,j,k)*(xn+one/xn)

c             write(*,*) 'soft ',msq(j,k)
c             write(*,*) 'soft 1',eik17_2*msq_cs(1,j,k)*xn
c             write(*,*) 'soft 2',eik17_5*msq_cs(1,j,k)*xn
c             write(*,*) 'soft 3',eik57_1*msq_cs(1,j,k)*xn
c--- done with all gluons QE.e-_dplike
             elseif ((j > 0) .and. (k==0)) then
        msq(j,k)=0._dp
c-- (a,b,c) = (5,6,7) and (5,7,6)
        do i=1,2
        eikc=eikf(s,b(i),c(i),5)*msq_ac(perm(a(i),b(i)),1,j,k)*xn
        msq(j,k)=msq(j,k)+eikc
        eikc=eikf(s,5,c(i),b(i))*msq_ac(perm(a(i),b(i)),1,j,k)*xn
        msq(j,k)=msq(j,k)+eikc
        eikc=(eikf(s,c(i),b(i),2)+eikf(s,2,b(i),c(i)))
     &       *msq_ac(perm(a(i),c(i)),1,j,k)*xn
        msq(j,k)=msq(j,k)+eikc
        enddo
        
c-- (a,b,c) = (7,5,6) and (6,5,7)
        do i=3,5,2
        eikc=(eikf(s,5,c(i),2)+eikf(s,2,c(i),5))
     &       *msq_ac(perm(b(i),a(i)),2,j,k)*xn
        msq(j,k)=msq(j,k)+eikc
        eikc=eikf(s,1,a(i),2)*msq_ac(perm(b(i),c(i)),1,j,k)*xn
        msq(j,k)=msq(j,k)+eikc
        eikc=eikf(s,2,a(i),1)*msq_ac(perm(b(i),c(i)),1,j,k)*xn
        msq(j,k)=msq(j,k)+eikc
       enddo
        
c-- (a,b,c) = (6,7,5) and (7,6,5)
        do i=4,6,2
        eikc=(eikf(s,a(i),b(i),2)+eikf(s,2,b(i),a(i)))
     &       *msq_ac(perm(c(i),a(i)),2,j,k)*xn
        msq(j,k)=msq(j,k)+eikc
        eikc=(eikf(s,1,a(i),b(i))+eikf(s,b(i),a(i),1))
     &       *msq_ac(perm(c(i),b(i)),2,j,k)*xn
        msq(j,k)=msq(j,k)+eikc
        enddo
        
            elseif ((j == 0) .and. (k>0)) then
        msq(j,k)=0._dp
c-- (a,b,c) = (5,6,7) and (5,7,6)
        do i=1,2
        eikc=eikf(s,b(i),c(i),5)*msq_ac(perm(a(i),b(i)),1,j,k)*xn
        msq(j,k)=msq(j,k)+eikc
        eikc=eikf(s,5,c(i),b(i))*msq_ac(perm(a(i),b(i)),1,j,k)*xn
        msq(j,k)=msq(j,k)+eikc
        eikc=(eikf(s,c(i),b(i),1)+eikf(s,1,b(i),c(i)))
     &       *msq_ac(perm(a(i),c(i)),1,j,k)*xn
        msq(j,k)=msq(j,k)+eikc
        enddo
        
c-- (a,b,c) = (7,5,6) and (6,5,7)
        do i=3,5,2
        eikc=(eikf(s,5,c(i),1)+eikf(s,1,c(i),5))
     &       *msq_ac(perm(b(i),a(i)),2,j,k)*xn
        msq(j,k)=msq(j,k)+eikc
        eikc=eikf(s,2,a(i),1)*msq_ac(perm(b(i),c(i)),1,j,k)*xn
        msq(j,k)=msq(j,k)+eikc
        eikc=eikf(s,1,a(i),2)*msq_ac(perm(b(i),c(i)),1,j,k)*xn
        msq(j,k)=msq(j,k)+eikc
       enddo
        
c-- (a,b,c) = (6,7,5) and (7,6,5)
        do i=4,6,2
        eikc=(eikf(s,a(i),b(i),1)+eikf(s,1,b(i),a(i)))
     &       *msq_ac(perm(c(i),a(i)),2,j,k)*xn
        msq(j,k)=msq(j,k)+eikc
        eikc=(eikf(s,2,a(i),b(i))+eikf(s,b(i),a(i),2))
     &       *msq_ac(perm(c(i),b(i)),2,j,k)*xn
        msq(j,k)=msq(j,k)+eikc
        enddo

            elseif ((j == 0) .and. (k<0)) then
        msq(j,k)=0._dp
c-- (a,b,c) = (5,6,7) and (5,7,6)
        do i=1,2
        eikc=(eikf(s,b(i),c(i),2)+eikf(s,2,c(i),b(i)))
     &       *msq_ac(perm(a(i),b(i)),1,j,k)*xn
        msq(j,k)=msq(j,k)+eikc
        eikc=(eikf(s,c(i),b(i),1)+eikf(s,1,b(i),c(i)))
     &       *msq_ac(perm(a(i),c(i)),1,j,k)*xn
        msq(j,k)=msq(j,k)+eikc
        enddo
        
c-- (a,b,c) = (7,5,6) and (6,5,7)
        do i=3,5,2
        eikc=(eikf(s,5,a(i),1)+eikf(s,1,a(i),5))
     &       *msq_ac(perm(b(i),c(i)),1,j,k)*xn
        msq(j,k)=msq(j,k)+eikc
        eikc=eikf(s,1,c(i),2)*msq_ac(perm(b(i),a(i)),2,j,k)*xn
        msq(j,k)=msq(j,k)+eikc
        eikc=eikf(s,2,c(i),1)*msq_ac(perm(b(i),a(i)),2,j,k)*xn
        msq(j,k)=msq(j,k)+eikc
       enddo
        
c-- (a,b,c) = (6,7,5) and (7,6,5)
        do i=4,6,2
        eikc=eikf(s,5,a(i),b(i))*msq_ac(perm(c(i),b(i)),2,j,k)*xn
        msq(j,k)=msq(j,k)+eikc
        eikc=eikf(s,b(i),a(i),5)*msq_ac(perm(c(i),b(i)),2,j,k)*xn
        msq(j,k)=msq(j,k)+eikc
        eikc=(eikf(s,a(i),b(i),1)+eikf(s,1,b(i),a(i)))
     &       *msq_ac(perm(c(i),a(i)),2,j,k)*xn
        msq(j,k)=msq(j,k)+eikc
        enddo

            elseif ((j < 0) .and. (k==0)) then
        msq(j,k)=0._dp
c-- (a,b,c) = (5,6,7) and (5,7,6)
        do i=1,2
        eikc=(eikf(s,b(i),c(i),1)+eikf(s,1,c(i),b(i)))
     &       *msq_ac(perm(a(i),b(i)),1,j,k)*xn
        msq(j,k)=msq(j,k)+eikc
        eikc=(eikf(s,c(i),b(i),2)+eikf(s,2,b(i),c(i)))
     &       *msq_ac(perm(a(i),c(i)),1,j,k)*xn
        msq(j,k)=msq(j,k)+eikc
        enddo
        
c-- (a,b,c) = (7,5,6) and (6,5,7)
        do i=3,5,2
        eikc=(eikf(s,5,a(i),2)+eikf(s,2,a(i),5))
     &       *msq_ac(perm(b(i),c(i)),1,j,k)*xn
        msq(j,k)=msq(j,k)+eikc
        eikc=eikf(s,1,c(i),2)*msq_ac(perm(b(i),a(i)),2,j,k)*xn
        msq(j,k)=msq(j,k)+eikc
        eikc=eikf(s,2,c(i),1)*msq_ac(perm(b(i),a(i)),2,j,k)*xn
        msq(j,k)=msq(j,k)+eikc
       enddo
        
c-- (a,b,c) = (6,7,5) and (7,6,5)
        do i=4,6,2
        eikc=eikf(s,5,a(i),b(i))*msq_ac(perm(c(i),b(i)),2,j,k)*xn
        msq(j,k)=msq(j,k)+eikc
        eikc=eikf(s,b(i),a(i),5)*msq_ac(perm(c(i),b(i)),2,j,k)*xn
        msq(j,k)=msq(j,k)+eikc
        eikc=(eikf(s,a(i),b(i),2)+eikf(s,2,b(i),a(i)))
     &       *msq_ac(perm(c(i),a(i)),2,j,k)*xn
        msq(j,k)=msq(j,k)+eikc
        enddo
 
        endif
      if (msq(j,k) .ne. 0._dp) write(*,*) 'soft ',j,k,msq(j,k)           
      enddo
      enddo
      endif
      
      if (Qflag) then
      call dotem(7,p,s)
      eik56=two*gsq*s(5,6)/(s(5,7)*s(6,7))
      eik15=two*gsq*s(1,5)/(s(1,7)*s(5,7))
      eik16=two*gsq*s(1,6)/(s(1,7)*s(6,7))
      eik26=two*gsq*s(2,6)/(s(2,7)*s(6,7))
      eik25=two*gsq*s(2,5)/(s(2,7)*s(5,7))
      eik12=two*gsq*s(1,2)/(s(1,7)*s(2,7))

      sub57_6(qq)=+s(5,6)/(s(5,7)+s(6,7))/s(5,7)*two*gsq
      sub67_5(qq)=+s(5,6)/(s(5,7)+s(6,7))/s(6,7)*two*gsq

      sub17_5(qq)=+s(1,5)/(s(5,7)-s(1,7))/s(1,7)*two*gsq
      sub57_1(qq)=-s(1,5)/(s(5,7)-s(1,7))/s(5,7)*two*gsq

      sub17_6(qq)=+s(1,6)/(s(6,7)-s(1,7))/s(1,7)*two*gsq
      sub67_1(qq)=-s(1,6)/(s(6,7)-s(1,7))/s(6,7)*two*gsq

      sub27_6(qq)=+s(2,6)/(s(6,7)-s(2,7))/s(2,7)*two*gsq
      sub67_2(qq)=-s(2,6)/(s(6,7)-s(2,7))/s(6,7)*two*gsq

      sub27_5(qq)=+s(2,5)/(s(5,7)-s(2,7))/s(2,7)*two*gsq
      sub57_2(qq)=-s(2,5)/(s(5,7)-s(2,7))/s(5,7)*two*gsq

      sub17_2(qq)=+s(1,2)/(s(1,7)+s(2,7))/s(1,7)*two*gsq
      sub27_1(qq)=+s(1,2)/(s(1,7)+s(2,7))/s(2,7)*two*gsq


      call qqb_w2jet(p,msq0)
 

      
      do j=-nf,nf
      do k=-nf,nf
            if ((j > 0) .and. (k>0)) then
c---QQ
       msq(j,k)=msq(j,k)
     & +sub57_6(qq)
     & *((xn+1._dp/xn)*mqq(0,j,k)+two*(mqq(1,j,k)+mqq(2,j,k))/xn)
     & +sub67_5(qq)
     & *((xn+1._dp/xn)*mqq(0,j,k)+two*(mqq(1,j,k)+mqq(2,j,k))/xn)
     & +sub17_5(qq)
     & *((xn-two/xn)*mqq(2,j,k)-mqq(0,j,k)/xn-mqq(1,j,k)/xn)
     & +sub57_1(qq)
     & *((xn-two/xn)*mqq(2,j,k)-mqq(0,j,k)/xn-mqq(1,j,k)/xn)
     & +sub17_6(qq)
     & *((xn-two/xn)*mqq(1,j,k)-mqq(0,j,k)/xn-mqq(2,j,k)/xn)
     & +sub67_1(qq)
     & *((xn-two/xn)*mqq(1,j,k)-mqq(0,j,k)/xn-mqq(2,j,k)/xn)
     & +sub27_6(qq)
     & *((xn-two/xn)*mqq(2,j,k)-mqq(0,j,k)/xn-mqq(1,j,k)/xn)
     & +sub67_2(qq)
     & *((xn-two/xn)*mqq(2,j,k)-mqq(0,j,k)/xn-mqq(1,j,k)/xn)
     & +sub27_5(qq)
     & *((xn-two/xn)*mqq(1,j,k)-mqq(0,j,k)/xn-mqq(2,j,k)/xn)
     & +sub57_2(qq)
     & *((xn-two/xn)*mqq(1,j,k)-mqq(0,j,k)/xn-mqq(2,j,k)/xn)
     & +sub17_2(qq)
     & *((xn+1._dp/xn)*mqq(0,j,k)+two*(mqq(1,j,k)+mqq(2,j,k))/xn)
     & +sub27_1(qq)
     & *((xn+1._dp/xn)*mqq(0,j,k)+two*(mqq(1,j,k)+mqq(2,j,k))/xn)

            elseif ((j < 0) .and. (k<0)) then
c---QbarQbar
       msq(j,k)=msq(j,k)
     & +sub57_6(qq)
     & *((xn+1._dp/xn)*mqq(0,j,k)+two*(mqq(1,j,k)+mqq(2,j,k))/xn)
     & +sub67_5(qq)
     & *((xn+1._dp/xn)*mqq(0,j,k)+two*(mqq(1,j,k)+mqq(2,j,k))/xn)
     & +sub17_5(qq)
     & *((xn-two/xn)*mqq(2,j,k)-mqq(0,j,k)/xn-mqq(1,j,k)/xn)
     & +sub57_1(qq)
     & *((xn-two/xn)*mqq(2,j,k)-mqq(0,j,k)/xn-mqq(1,j,k)/xn)
     & +sub17_6(qq)
     & *((xn-two/xn)*mqq(1,j,k)-mqq(0,j,k)/xn-mqq(2,j,k)/xn)
     & +sub67_1(qq)
     & *((xn-two/xn)*mqq(1,j,k)-mqq(0,j,k)/xn-mqq(2,j,k)/xn)
     & +sub27_6(qq)
     & *((xn-two/xn)*mqq(2,j,k)-mqq(0,j,k)/xn-mqq(1,j,k)/xn)
     & +sub67_2(qq)
     & *((xn-two/xn)*mqq(2,j,k)-mqq(0,j,k)/xn-mqq(1,j,k)/xn)
     & +sub27_5(qq)
     & *((xn-two/xn)*mqq(1,j,k)-mqq(0,j,k)/xn-mqq(2,j,k)/xn)
     & +sub57_2(qq)
     & *((xn-two/xn)*mqq(1,j,k)-mqq(0,j,k)/xn-mqq(2,j,k)/xn)
     & +sub17_2(qq)
     & *((xn+1._dp/xn)*mqq(0,j,k)+two*(mqq(1,j,k)+mqq(2,j,k))/xn)
     & +sub27_1(qq)
     & *((xn+1._dp/xn)*mqq(0,j,k)+two*(mqq(1,j,k)+mqq(2,j,k))/xn)
            elseif ((j > 0) .and. (k<0)) then
c---QQbar
       msq(j,k)=msq(j,k)
     & +sub57_6(qq)
     & *((xn-two/xn)*mqq(2,j,k)-mqq(0,j,k)/xn-mqq(1,j,k)/xn)
     & +sub67_5(qq)
     & *((xn-two/xn)*mqq(2,j,k)-mqq(0,j,k)/xn-mqq(1,j,k)/xn)
     & +sub17_5(qq)
     & *((xn-two/xn)*mqq(1,j,k)-mqq(0,j,k)/xn-mqq(2,j,k)/xn)
     & +sub57_1(qq)
     & *((xn-two/xn)*mqq(1,j,k)-mqq(0,j,k)/xn-mqq(2,j,k)/xn)
     & +sub17_6(qq)
     & *((xn+1._dp/xn)*mqq(0,j,k)+two*(mqq(1,j,k)+mqq(2,j,k))/xn)
     & +sub67_1(qq)
     & *((xn+1._dp/xn)*mqq(0,j,k)+two*(mqq(1,j,k)+mqq(2,j,k))/xn)
     & +sub27_6(qq)
     & *((xn-two/xn)*mqq(1,j,k)-mqq(0,j,k)/xn-mqq(2,j,k)/xn)
     & +sub67_2(qq)
     & *((xn-two/xn)*mqq(1,j,k)-mqq(0,j,k)/xn-mqq(2,j,k)/xn)
     & +sub27_5(qq)
     & *((xn+1._dp/xn)*mqq(0,j,k)+two*(mqq(1,j,k)+mqq(2,j,k))/xn)
     & +sub57_2(qq)
     & *((xn+1._dp/xn)*mqq(0,j,k)+two*(mqq(1,j,k)+mqq(2,j,k))/xn)
     & +sub17_2(qq)
     & *((xn-two/xn)*mqq(2,j,k)-mqq(0,j,k)/xn-mqq(1,j,k)/xn)
     & +sub27_1(qq)
     & *((xn-two/xn)*mqq(2,j,k)-mqq(0,j,k)/xn-mqq(1,j,k)/xn)
            elseif ((j < 0) .and. (k>0)) then
c---QbarQ
       msq(j,k)=msq(j,k)
     & +sub57_6(qq)
     & *((xn-two/xn)*mqq(2,j,k)-mqq(0,j,k)/xn-mqq(1,j,k)/xn)
     & +sub67_5(qq)
     & *((xn-two/xn)*mqq(2,j,k)-mqq(0,j,k)/xn-mqq(1,j,k)/xn)
     & +sub17_5(qq)
     & *((xn+1._dp/xn)*mqq(0,j,k)+two*(mqq(1,j,k)+mqq(2,j,k))/xn)
     & +sub57_1(qq)
     & *((xn+1._dp/xn)*mqq(0,j,k)+two*(mqq(1,j,k)+mqq(2,j,k))/xn)
     & +sub17_6(qq)
     & *((xn-two/xn)*mqq(1,j,k)-mqq(0,j,k)/xn-mqq(2,j,k)/xn)
     & +sub67_1(qq)
     & *((xn-two/xn)*mqq(1,j,k)-mqq(0,j,k)/xn-mqq(2,j,k)/xn)
     & +sub27_6(qq)
     & *((xn+1._dp/xn)*mqq(0,j,k)+two*(mqq(1,j,k)+mqq(2,j,k))/xn)
     & +sub67_2(qq)
     & *((xn+1._dp/xn)*mqq(0,j,k)+two*(mqq(1,j,k)+mqq(2,j,k))/xn)
     & +sub27_5(qq)
     & *((xn-two/xn)*mqq(1,j,k)-mqq(0,j,k)/xn-mqq(2,j,k)/xn)
     & +sub57_2(qq)
     & *((xn-two/xn)*mqq(1,j,k)-mqq(0,j,k)/xn-mqq(2,j,k)/xn)
     & +sub17_2(qq)
     & *((xn-two/xn)*mqq(2,j,k)-mqq(0,j,k)/xn-mqq(1,j,k)/xn)
     & +sub27_1(qq)
     & *((xn-two/xn)*mqq(2,j,k)-mqq(0,j,k)/xn-mqq(1,j,k)/xn)

            elseif ((j > 0) .and. (k==0)) then
               msq(j,k)=msq(j,k)
            elseif ((j < 0) .and. (k==0)) then
               msq(j,k)=msq(j,k)
            elseif ((j == 0) .and. (k>0)) then
               msq(j,k)=msq(j,k)
            elseif ((j == 0) .and. (k<0)) then
               msq(j,k)=msq(j,k)
            endif

      enddo
      enddo
      endif
      return
      end


      
      function eikf(s,a,b,c)
      implicit none
      include 'types.f'
      real(dp):: eikf
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'qcdcouple.f'
      integer:: a,b,c
      real(dp):: sa,sb,sc
      real(dp):: s(mxpart,mxpart)

      sa=1._dp
      sb=1._dp
      sc=1._dp
      if (a <= 2) sa=-1._dp
      if (b <= 2) sb=-1._dp
      if (c <= 2) sc=-1._dp
            
      eikf=sa*sc*s(a,c)/(sa*sb*s(a,b)+sb*sc*s(b,c))
     &                 /(sa*sb*s(a,b))*two*gsq
            
       return
       end     
            
            
