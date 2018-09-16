      double complex function amp_qqggg(i1,qh,i2,h2,i3,h3,i4,h4,i5,lh,
     . j6,j7)
C--Results taken from Nagy and Trocsanyi-hep-ph/9806317
C--Appendix A
      implicit none
      include 'constants.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      integer i1,i2,i3,i4,i5,i6,i7,j,k,qh,h2,h3,h4,lh,g1,g2,g3,g4,lg,
     . j6,j7
      double complex t2a,t3a,t3b
      double complex xa(mxpart,mxpart),xb(mxpart,mxpart)
      double precision t123,t234,t345,t567,t167

      t2a(i1,i2,i3,i4)=xa(i1,i2)*xb(i2,i4)+xa(i1,i3)*xb(i3,i4)
c      t2b(i1,i2,i3,i4)=xb(i1,i2)*xa(i2,i4)+xb(i1,i3)*xa(i3,i4)
      t3a(i1,i2,i3,i4,i5,i6)=xa(i1,i2)*xb(i2,i4)*xa(i4,i6)
     .                      +xa(i1,i2)*xb(i2,i5)*xa(i5,i6)
     .                      +xa(i1,i3)*xb(i3,i4)*xa(i4,i6)
     .                      +xa(i1,i3)*xb(i3,i5)*xa(i5,i6)
      t3b(i1,i2,i3,i4,i5,i6)=xb(i1,i2)*xa(i2,i4)*xb(i4,i6)
     .                      +xb(i1,i2)*xa(i2,i5)*xb(i5,i6)
     .                      +xb(i1,i3)*xa(i3,i4)*xb(i4,i6)
     .                      +xb(i1,i3)*xa(i3,i5)*xb(i5,i6)

C---A(qh,h2,h3,h4,lh)
C---h=1 LH
C---h=2 RH

      g1=qh
      g2=h2
      g3=h3
      g4=h4
      lg=lh
      i6=j6
      i7=j7
      if (g1*lg.eq.2)then
      i7=j6
      i6=j7
      lg=3-lg  
      endif 

      if (g1.eq.2) then
      do j=1,mxpart
      do k=1,mxpart
      xa(j,k)=za(j,k)
      xb(j,k)=zb(j,k)
      enddo
      enddo
      elseif (g1.eq.1) then
      g1=2
      g2=3-g2
      g3=3-g3
      g4=3-g4
      lg=3-lg
      do j=1,mxpart
      do k=1,mxpart
      xa(j,k)=zb(k,j)
      xb(j,k)=za(k,j)
      enddo
      enddo

      endif
      t123=s(i1,i2)+s(i2,i3)+s(i3,i1)
      t167=s(i1,i6)+s(i6,i7)+s(i7,i1)
      t567=s(i5,i6)+s(i6,i7)+s(i7,i5)
      t345=s(i3,i4)+s(i4,i5)+s(i5,i3)
      t234=s(i2,i3)+s(i3,i4)+s(i4,i2)

C--A43 (2,2,2,2,2)
      if (
     . (g1.eq.2).and.(g2.eq.2).and.(g3.eq.2).and.(g4.eq.2).and.(lg.eq.2)
     . ) then 
      amp_qqggg =
     .-xa(i6,i5)**2*xb(i6,i7)/(xa(i1,i2)*xa(i2,i3)*xa(i3,i4)*xa(i4,i5))
C--A49
C        A(2,1,1,1,2) =
      elseif 
     .((g1.eq.2).and.(g2.eq.1).and.(g3.eq.1).and.(g4.eq.1).and.(lg.eq.2)
     . ) then
      amp_qqggg =
     .+xb(i1,i7)**2*xa(i6,i7)/(xb(i1,i2)*xb(i2,i3)*xb(i3,i4)*xb(i4,i5))
c--- extra minus sign after CP symmetry
C--A44
C      A(2,2,2,1,2) =
      elseif 
     .((g1.eq.2).and.(g2.eq.2).and.(g3.eq.2).and.(g4.eq.1).and.(lg.eq.2)
     . ) then
      amp_qqggg =
     . +xa(i6,i5)*t2a(i4,i5,i6,i7)/(xa(i2,i3)*xa(i3,i4)*xb(i3,i4)*t567)
     .*(xb(i2,i3)*t2a(i4,i2,i3,i1)/t234+t2a(i4,i1,i2,i3)/xa(i1,i2)) 
     . +t2a(i4,i5,i6,i7)*t2a(i6,i5,i4,i3)
     . /(xa(i1,i2)*xa(i2,i3)*xa(i3,i4)*xb(i3,i4)*xb(i4,i5))
     . -t2a(i4,i1,i2,i7)*t2a(i6,i5,i4,i3)*xa(i4,i5)*xb(i5,i3)
     . /(xa(i1,i2)*xa(i2,i4)*xb(i4,i5)*s(i3,i4)*t345) 
     . -xb(i1,i7)*t2a(i6,i1,i7,i2)*xa(i4,i5)**2*xb(i5,i3)**2
     . /(xb(i4,i5)*xa(i4,i2)*s(i3,i4)*t345*t167) 
     . +(xb(i1,i7)*xa(i4,i5)*xb(i5,i3))
     . /(xa(i2,i3)*xb(i3,i4)*xb(i4,i5)*t167)
     . *(t2a(i6,i1,i7,i2)/xa(i3,i4)+t2a(i6,i1,i7,i3)/xa(i2,i4))
     . -xb(i1,i7)*xa(i4,i5)*xb(i2,i3)/(xa(i2,i3)*xb(i3,i4)*t234*t167)
     . *(t2a(i6,i1,i7,i2)*xa(i2,i4)/xa(i3,i4)+t2a(i6,i1,i7,i3)) 


C--A45
C      A(2,2,1,2,2)=
      elseif 
     .((g1.eq.2).and.(g2.eq.2).and.(g3.eq.1).and.(g4.eq.2).and.(lg.eq.2)
     . ) then
      amp_qqggg =
     . +xa(i3,i1)*xb(i1,i2)*t2a(i3,i1,i2,i4)*t2a(i3,i5,i6,i7)*xa(i6,i5)
     . /(xa(i1,i2)*xa(i3,i4)*s(i2,i3)*t123*t567) 
     . +xa(i3,i1)*xb(i1,i2)*t2a(i3,i1,i2,i7)*xa(i6,i5)*xa(i3,i5)
     . /(xa(i1,i2)*xa(i3,i4)*xa(i4,i5)*s(i2,i3)*t123)
     . -t2a(i3,i1,i2,i7)*t2a(i6,i1,i7,i2)*xa(i3,i5)**2
     . /(xa(i1,i2)*xa(i2,i3)*xa(i3,i4)*xa(i4,i5)*xb(i2,i3)*t345) 
     . +t2a(i3,i1,i2,i7)*t2a(i6,i5,i3,i4)*xa(i3,i5)*xb(i4,i2)
     . /(xb(i2,i3)*xa(i1,i2)*xa(i2,i3)*s(i3,i4)*t345)
     . -(xb(i4,i2)**2*xb(i1,i7)*t3a(i6,i1,i7,i2,i4,i3)*xa(i3,i5))
     . /(s(i2,i3)*s(i3,i4)*t234*t167) 
     . +(xb(i4,i2)*t2a(i3,i5,i6,i7)*xa(i6,i5))/(s(i2,i3)*s(i3,i4)*t567)
     . *((xb(i4,i2)*t2a(i3,i2,i4,i1))/t234-t2a(i3,i1,i2,i4)/xa(i1,i2))
     . +xb(i1,i7)*t2a(i6,i1,i7,i2)*xa(i3,i5)**2/(s(i2,i3)*t345*t167)
     . *(t2a(i3,i5,i4,i2)/(xa(i3,i4)*xa(i4,i5))
     . +xb(i4,i2)*xb(i5,i4)/s(i3,i4))
 
C---A46
C      A(2,1,2,2,2) =
      elseif 
     .((g1.eq.2).and.(g2.eq.1).and.(g3.eq.2).and.(g4.eq.2).and.(lg.eq.2)
     . ) then
      amp_qqggg =
     . +xb(i4,i3)**2*t2a(i2,i3,i4,i1)*t2a(i2,i5,i6,i7)*xa(i6,i5)
     . /(s(i2,i3)*s(i3,i4)*t234*t567)
     . +xb(i1,i3)*t2a(i2,i3,i4,i1)*t2a(i2,i5,i6,i7)*xa(i6,i5)
     . /(xb(i1,i2)*xa(i3,i4)*xa(i4,i2)*s(i2,i3)*t567)
     . -xb(i1,i3)**2*t2a(i2,i1,i3,i4)*t2a(i2,i5,i6,i7)*xa(i6,i5)
     . /(xb(i1,i2)*xa(i2,i4)*s(i2,i3)*t123*t567)
     . -xb(i1,i3)**2*t2a(i2,i1,i3,i7)*xa(i6,i5)*xa(i2,i5)
     . /(xb(i1,i2)*xa(i2,i4)*xa(i4,i5)*s(i2,i3)*t123)
     . -xb(i1,i7)*xa(i2,i5)*xb(i4,i3)**2*t3a(i6,i1,i7,i3,i4,i2)
     . /(s(i2,i3)*s(i3,i4)*t234*t167) 
     . +xb(i1,i3)*xb(i1,i7)*xa(i2,i5)
     . /(xb(i1,i2)*s(i2,i3)*xa(i2,i4)*t345)
     . *(t2a(i6,i5,i4,i3)*(xa(i2,i5)/xa(i4,i5)-xa(i3,i2)/xa(i3,i4))
     . - t2a(i6,i5,i3,i4)*(xa(i4,i2))/xa(i3,i4)) 
     . +xb(i1,i7)*t2a(i6,i1,i7,i3)*xa(i2,i5)
     . /(xa(i2,i4)*s(i2,i3)*t345*t167)
     . *(t2a(i2,i5,i4,i3)
     . *((xa(i2,i5))/(xa(i4,i5))-(xa(i3,i2))/(xa(i3,i4)))
     . -t2a(i2,i5,i3,i4)*xa(i4,i2)/xa(i3,i4))
 
C--A47
C      A(2,2,1,1,2) =
      elseif 
     .((g1.eq.2).and.(g2.eq.2).and.(g3.eq.1).and.(g4.eq.1).and.(lg.eq.2)
     . ) then
      amp_qqggg =
     . +xb(i1,i2)*t3b(i2,i3,i4,i5,i6,i7)*xa(i6,i5)/(s(i2,i3)*t567)
     . *(xa(i4,i3)**2/(s(i3,i4)*t234)
     . -xa(i3,i1)/(xb(i3,i4)*xb(i4,i2)*xa(i1,i2)))
     . -xa(i3,i1)**2*xb(i1,i2)**2*t2a(i4,i5,i6,i7)*xa(i6,i5)
     . /(xa(i1,i2)*xb(i2,i4)*s(i2,i3)*t123*t567)
     . +xa(i3,i1)*xb(i1,i2)*t2a(i3,i1,i2,i7)*t2a(i6,i5,i4,i2)
     . /(xa(i1,i2)*xb(i2,i4)*xb(i4,i5)*s(i2,i3)*t123) 
     . +t2a(i3,i1,i2,i7)*t2a(i6,i1,i7,i2)
     . /(xa(i1,i2)*xb(i3,i4)*xb(i4,i5)*s(i2,i3))
     . +xb(i1,i7)*t2a(i6,i1,i7,i2)*t2a(i3,i5,i4,i2)
     . /(xb(i3,i4)*xb(i4,i5)*s(i2,i3)*t167)
     . -xb(i1,i7)*t2a(i6,i1,i7,i2)*t2a(i5,i3,i4,i2)*xa(i4,i3)**2
     . /(s(i2,i3)*s(i3,i4)*t234*t167)

C--A48
C      A(2,1,1,2,2) =
      elseif 
     .((g1.eq.2).and.(g2.eq.1).and.(g3.eq.1).and.(g4.eq.2).and.(lg.eq.2)
     . ) then
      amp_qqggg =
     . +xa(i2,i3)**2*xb(i1,i4)*t3b(i4,i2,i3,i5,i6,i7)*xa(i6,i5)
     . /(s(i2,i3)*s(i3,i4)*t234*t567)
     . -xb(i1,i4)*t2a(i3,i5,i6,i7)*xa(i6,i5)
     . /(xb(i4,i2)*s(i3,i4)*t123*t567)
     . *((xb(i4,i2)*t2a(i2,i1,i3,i4)+xb(i4,i3)*t2a(i3,i1,i2,i4))
     . /xb(i2,i3) 
     . -xb(i1,i4)*t2a(i3,i1,i2,i4)/xb(i1,i2)) 
     . -xb(i1,i4)*xa(i6,i5)*xa(i3,i5)
     . /(xb(i4,i2)*xa(i4,i5)*s(i3,i4)*t123)
     . *((xb(i4,i2)*t2a(i2,i1,i3,i7)+xb(i4,i3)*t2a(i3,i1,i2,i7))
     . /xb(i2,i3)
     .  -xb(i1,i4)*t2a(i3,i1,i2,i7)/xb(i1,i2)) 
     . +xb(i1,i4)*xb(i1,i7)*t2a(i6,i5,i3,i4)*xa(i3,i5)**2
     . /(xb(i1,i2)*xb(i2,i4)*xa(i4,i5)*s(i3,i4)*t345)
     . -xb(i1,i7)*t2a(i6,i1,i7,i4)*t2a(i2,i5,i3,i4)*xa(i3,i5)**2
     . /(xb(i4,i2)*xa(i4,i5)*s(i3,i4)*t345*t167)
     . -xb(i1,i7)*t2a(i6,i1,i7,i4)*t2a(i5,i2,i3,i4)*xa(i3,i5)
     . /(xa(i4,i5)*xb(i4,i2)*xb(i2,i3)*s(i3,i4)*t167)
     . -xb(i1,i7)*t2a(i6,i1,i7,i4)*t2a(i5,i2,i3,i4)*xa(i2,i3)**2
     . /(s(i2,i3)*s(i3,i4)*t234*t167)
 

      elseif 
     .((g1.eq.2).and.(g2.eq.1).and.(g3.eq.2).and.(g4.eq.1).and.(lg.eq.2)
     . ) then

C--A49
C----Constructed from A45 with P operation (above) and 6 and 7 exchanged
C      A(2,1,2,1,2)=
      amp_qqggg =
     . -xb(i5,i3)*xa(i4,i5)*t2a(i2,i4,i5,i3)*t2a(i6,i1,i7,i3)*xb(i1,i7)
     . /(xb(i4,i5)*xb(i2,i3)*s(i3,i4)*t345*t167)      
     . -xb(i5,i3)*xa(i4,i5)*t2a(i6,i4,i5,i3)*xb(i1,i7)*xb(i1,i3)
     . /(xb(i4,i5)*xb(i2,i3)*xb(i1,i2)*s(i3,i4)*t345)     
     . +t2a(i6,i4,i5,i3)*t2a(i4,i5,i6,i7)*xb(i1,i3)**2
     . /(xb(i4,i5)*xb(i3,i4)*xb(i2,i3)*xb(i1,i2)*xa(i3,i4)*t123)     
     . -t2a(i6,i4,i5,i3)*t2a(i2,i1,i3,i7)*xb(i1,i3)*xa(i4,i2)
     . /(xa(i3,i4)*xb(i4,i5)*xb(i3,i4)*s(i2,i3)*t123)     
     . +(xa(i4,i2)**2*xa(i6,i5)*t3b(i3,i2,i4,i5,i6,i7)*xb(i1,i3))
     . /(s(i3,i4)*s(i2,i3)*t234*t567)      
     . -(xa(i4,i2)*t2a(i6,i1,i7,i3)*xb(i1,i7))/(s(i3,i4)*s(i2,i3)*t167)
     . *((xa(i4,i2)*t2a(i5,i2,i4,i3))/t234-t2a(i2,i4,i5,i3)/xb(i4,i5))     
     . -xa(i6,i5)*t2a(i4,i5,i6,i7)*xb(i1,i3)**2/(s(i3,i4)*t123*t567)
     . *(t2a(i4,i1,i2,i3)/(xb(i2,i3)*xb(i1,i2))
     . +xa(i4,i2)*xa(i2,i1)/s(i2,i3))
c--- extra minus sign after CP symmetry
 
      endif

c--- divide out by the photon propagator that we put back in later

      amp_qqggg = amp_qqggg/s(j6,j7)

      return 
      end

