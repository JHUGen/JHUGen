      subroutine qqb_wbb_soft(P,msq)
      implicit none
      include 'types.f'
c---Matrix element SUBTRACTION squared averaged over initial colors and spins
c     q(-p1)+qbar(-p2) -->  b(p5)+bb(p6) + W +g(p7)
c                                          |
c                                          --> nu(p3)+e^+(p4)
c                            
c   positively charged W only
c--all momenta incoming
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'qcdcouple.f'
      include 'masses.f'
      integer:: j,k

      real(dp):: P(mxpart,4),msq(-nf:nf,-nf:nf),
     & msq0(-nf:nf,-nf:nf),s(mxpart,mxpart),facqqb,facqbq
      real(dp):: eik17_2,eik17_5,eik17_6,eik27_5,eik27_6,eik57_6
      real(dp):: eik27_1,eik57_1,eik67_1,eik57_2,eik67_2,eik67_5
      
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0._dp
      enddo
      enddo
      call dotem(7,p,s)
c + p1Dp7^-1*p2Dp7^-1
c  * ( - n^-1*xlo^2*p1Dp2 )
c + p7Dp5^-1*p7Dp6^-1
c  * ( - n^-1*xlo^2*p5Dp6 )
c + p1Dp7^-1*p7Dp6^-1
c  * ( 2*n^-1*xlo^2*p1Dp6 )
c + p2Dp7^-1*p7Dp5^-1
c  * ( 2*n^-1*xlo^2*p2Dp5 )
c + p1Dp7^-1*p7Dp5^-1
c  * ( n*xlo^2*p1Dp5 - 2*n^-1*xlo^2*p1Dp5 )
c + p2Dp7^-1*p7Dp6^-1
c  * ( n*xlo^2*p2Dp6 - 2*n^-1*xlo^2*p2Dp6 )

c      eik12=s(1,2)/s(1,7)/s(2,7)
c      eik56=s(5,6)/s(5,7)/s(6,7)

c      eik15=s(1,5)/s(1,7)/s(5,7)
c      eik16=s(1,6)/s(1,7)/s(6,7)

c      eik25=s(2,5)/s(2,7)/s(5,7)
c      eik26=s(2,6)/s(2,7)/s(6,7)


c      facqqb=((+two*(eik17_6+eik67_1+eik27_5+eik57_2
c     & -eik27_6-eik67_2-eik17_5-eik57_1)
c     & -eik17_2-eik27_1-eik57_6-eik67_5
c     & )/xn
c     & +xn*(eik17_5+eik57_1+eik27_6+eik67_2))


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




c      facqbq=((+two*(eik27_6+eik67_2+eik17_5+eik57_1
c     & -eik17_6-eik67_1-eik27_5-eik57_2)
c     & -eik27_1-eik17_2-eik57_6-eik67_5
c     & )/xn
c     & +xn*(eik27_5+eik57_2+eik17_6+eik67_1))

      facqqb=(eik27_6+eik67_2+eik17_5+eik57_1)*(xn-two/xn)
     &      +(eik17_6+eik67_1+eik27_5+eik57_2)*(two/xn)
     &      +(eik27_1+eik17_2+eik57_6+eik67_5)*(-one/xn)
      facqbq=(eik17_6+eik67_1+eik27_5+eik57_2)*(xn-two/xn)
     &      +(eik27_6+eik67_2+eik17_5+eik57_1)*(two/xn)
     &      +(eik17_2+eik27_1+eik57_6+eik67_5)*(-one/xn)


c??????? exchange ???/

      call qqb_wbb(p,msq0)

c      write(6,*) 'p(1,4)',p(1,4)
c      write(6,*) 'p(2,4)',p(2,4)
c      write(6,*) 'msq0(2,-1)',msq0(2,-1),'msq0(-1,2)',msq0(-1,2)
c      write(6,*) 
      
      do j=2,2
      do k=-1,-1
            if ((j > 0) .and. (k<0)) then
               msq(j,k)=facqqb*msq0(j,k)
      
       write(6,*) '-eik17_2*msq0(j,k)/xn',-eik17_2*msq0(j,k)/xn
       write(6,*) '-eik27_1*msq0(j,k)/xn',-eik27_1*msq0(j,k)/xn

       write(6,*) '-eik57_6*msq0(j,k)/xn',-eik57_6*msq0(j,k)/xn
       write(6,*) '-eik67_5*msq0(j,k)/xn',-eik67_5*msq0(j,k)/xn

       write(6,*) 'two/xn*(eik17_5+eik57_1)*msq0(j,k)',
     & +(xn-two/xn)*(eik17_5+eik57_1)*msq0(j,k)
       write(6,*) 'two/xn*(eik17_6+eik67_1)*msq0(j,k)',
     & two/xn*(eik17_6+eik67_1)*msq0(j,k)
       write(6,*) '(eik27_6+eik67_2)*msq0(j,k)',
     & +(xn-two/xn)*(eik27_6+eik67_2)*msq0(j,k)
       write(6,*) 'two/xn*(eik27_5+eik57_2)*msq0(j,k)',
     & +two/xn*(eik27_5+eik57_2)*msq0(j,k)


            elseif ((j < 0) .and. (k>0)) then
               msq(j,k)=facqbq*msq0(j,k)
            endif
      enddo
      enddo
      return
      end

 
