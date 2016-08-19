      subroutine h4qg(i1,i2,i3,i4,i5,q4ghsq,q4idghsq)
C-----Author R.K.Ellis
C-----November 2004
C-----Matrix element squared for 
C     q(p2)+Q(p4) --> q(p1)+Q(p3)+g(p5)+h
C-----Summed over colors and spins
      implicit none 
      include 'constants.f'
      include 'zprods_com.f'
      integer i1,i2,i3,i4,i5,hq1,hq2,hg,j
      double precision q4ghsq,q4idghsq,sign
      double complex ma(4,2,2,2),mb(4,2,2,2),
     .    q4ghppp1,q4ghppp3,q4ghpmp1,q4ghpmp3

C---Amplitudes are labelled by their 
C---1) Color structure m(1),m(2),m(2),m(3),m(4)
C---2) Helicity of incoming line 2
C---3) Helicity of incoming line 4
C---4) Helicity of gluon
C---
C---
C---  M=g^3/sqrt(2)
C---  *(+T(C5,i1,i4)*delta(i3,i2)*m(1)
C---    +T(C5,i3,i2)*delta(i1,i4)*m(2)
C---    +1/xn*T(C5,i1,i2)*delta(i3,i4)*m(3)
C---    +1/xn*T(C5,i3,i4)*delta(i1,i2)*m(4))
C---
C---
       

      ma(1,2,2,2)=+q4ghppp1(i1,i2,i3,i4,i5,za,zb)
      ma(3,2,2,2)=+q4ghppp3(i1,i2,i3,i4,i5,za,zb)
      ma(4,2,2,2)=+q4ghppp3(i3,i4,i1,i2,i5,za,zb)
      ma(2,2,2,2)=-ma(4,2,2,2)-ma(3,2,2,2)-ma(1,2,2,2)

      ma(1,1,1,2)=-q4ghppp1(i4,i3,i2,i1,i5,za,zb)
      ma(3,1,1,2)=-q4ghppp3(i2,i1,i4,i3,i5,za,zb)
      ma(4,1,1,2)=-q4ghppp3(i4,i3,i2,i1,i5,za,zb)
      ma(2,1,1,2)=-ma(4,1,1,2)-ma(3,1,1,2)-ma(1,1,1,2)
      
      ma(2,2,1,1)=-q4ghpmp1(i3,i4,i1,i2,i5,zb,za)
      ma(3,2,1,1)=+q4ghpmp3(i2,i1,i4,i3,i5,zb,za)
      ma(4,2,1,1)=-q4ghpmp3(i3,i4,i1,i2,i5,zb,za)
      ma(1,2,1,1)=-ma(2,2,1,1)-ma(3,2,1,1)-ma(4,2,1,1)

      ma(1,1,2,1)=-q4ghpmp1(i1,i2,i3,i4,i5,zb,za)
      ma(3,1,2,1)=-q4ghpmp3(i1,i2,i3,i4,i5,zb,za)
      ma(4,1,2,1)=+q4ghpmp3(i4,i3,i2,i1,i5,zb,za)
      ma(2,1,2,1)=-ma(3,1,2,1)-ma(4,1,2,1)-ma(1,1,2,1)

c      ma(1,1,1,1)=-q4ghppp1(i1,i2,i3,i4,i5,zb,za)
c      ma(3,1,1,1)=-q4ghppp3(i1,i2,i3,i4,i5,zb,za)
c      ma(4,1,1,1)=-q4ghppp3(i3,i4,i1,i2,i5,zb,za)
c      ma(2,1,1,1)=-ma(3,1,1,1)-ma(4,1,1,1)-ma(1,1,1,1)

c      ma(1,2,2,1)=+q4ghppp1(i4,i3,i2,i1,i5,zb,za)
c      ma(3,2,2,1)=+q4ghppp3(i2,i1,i4,i3,i5,zb,za)
c      ma(4,2,2,1)=+q4ghppp3(i4,i3,i2,i1,i5,zb,za)
c      ma(2,2,2,1)=-ma(1,2,2,1)-ma(3,2,2,1)-ma(4,2,2,1)

c      ma(2,1,2,2)=+q4ghpmp1(i3,i4,i1,i2,i5,za,zb)
c      ma(3,1,2,2)=-q4ghpmp3(i2,i1,i4,i3,i5,za,zb)
c      ma(4,1,2,2)=+q4ghpmp3(i3,i4,i1,i2,i5,za,zb)
c      ma(1,1,2,2)=-ma(3,1,2,2)-ma(4,1,2,2)-ma(2,1,2,2)

c      ma(1,2,1,2)=+q4ghpmp1(i1,i2,i3,i4,i5,za,zb)
c      ma(3,2,1,2)=+q4ghpmp3(i1,i2,i3,i4,i5,za,zb)
c      ma(4,2,1,2)=-q4ghpmp3(i4,i3,i2,i1,i5,za,zb)
c      ma(2,2,1,2)=-ma(4,2,1,2)-ma(3,2,1,2)-ma(1,2,1,2)

c--- additional sign for the crossed amplitudes with initial gluon
      if (i5 .ne. 5) then
        sign=-1d0
      else
        sign=+1d0
      endif
      
c--- fastest to obtain remaining amplitudes by symmetry
      do j=1,4
      ma(j,1,1,1)=sign*dconjg(ma(j,2,2,2))
      ma(j,2,2,1)=sign*dconjg(ma(j,1,1,2))
      ma(j,1,2,2)=sign*dconjg(ma(j,2,1,1))
      ma(j,2,1,2)=sign*dconjg(ma(j,1,2,1))
      enddo
      

C now do same piece of work for identical quarks
C reverse i2 <--> i4
C exchange helicities of 2 and 4
C multiply with minus sign
      mb(1,2,2,2)=-q4ghppp1(i1,i4,i3,i2,i5,za,zb)
      mb(3,2,2,2)=-q4ghppp3(i1,i4,i3,i2,i5,za,zb)
      mb(4,2,2,2)=-q4ghppp3(i3,i2,i1,i4,i5,za,zb)
      mb(2,2,2,2)=-mb(4,2,2,2)-mb(3,2,2,2)-mb(1,2,2,2)

      mb(1,1,1,2)=+q4ghppp1(i2,i3,i4,i1,i5,za,zb)
      mb(3,1,1,2)=+q4ghppp3(i4,i1,i2,i3,i5,za,zb)
      mb(4,1,1,2)=+q4ghppp3(i2,i3,i4,i1,i5,za,zb)
      mb(2,1,1,2)=-mb(4,1,1,2)-mb(3,1,1,2)-mb(1,1,1,2)
      
      mb(1,2,1,1)=+q4ghpmp1(i1,i4,i3,i2,i5,zb,za)
      mb(3,2,1,1)=+q4ghpmp3(i1,i4,i3,i2,i5,zb,za)
      mb(4,2,1,1)=-q4ghpmp3(i2,i3,i4,i1,i5,zb,za)
      mb(2,2,1,1)=-mb(3,2,1,1)-mb(4,2,1,1)-mb(1,2,1,1)

      mb(2,1,2,1)=+q4ghpmp1(i3,i2,i1,i4,i5,zb,za)
      mb(3,1,2,1)=-q4ghpmp3(i4,i1,i2,i3,i5,zb,za)
      mb(4,1,2,1)=+q4ghpmp3(i3,i2,i1,i4,i5,zb,za)
      mb(1,1,2,1)=-mb(2,1,2,1)-mb(3,1,2,1)-mb(4,1,2,1)

c      mb(1,1,1,1)=+q4ghppp1(i1,i4,i3,i2,i5,zb,za)
c      mb(3,1,1,1)=+q4ghppp3(i1,i4,i3,i2,i5,zb,za)
c      mb(4,1,1,1)=+q4ghppp3(i3,i2,i1,i4,i5,zb,za)
c      mb(2,1,1,1)=-mb(3,1,1,1)-mb(4,1,1,1)-mb(1,1,1,1)

c      mb(1,2,2,1)=-q4ghppp1(i2,i3,i4,i1,i5,zb,za)
c      mb(3,2,2,1)=-q4ghppp3(i4,i1,i2,i3,i5,zb,za)
c      mb(4,2,2,1)=-q4ghppp3(i2,i3,i4,i1,i5,zb,za)
c      mb(2,2,2,1)=-mb(1,2,2,1)-mb(3,2,2,1)-mb(4,2,2,1)

c      mb(1,1,2,2)=-q4ghpmp1(i1,i4,i3,i2,i5,za,zb)
c      mb(3,1,2,2)=-q4ghpmp3(i1,i4,i3,i2,i5,za,zb)
c      mb(4,1,2,2)=+q4ghpmp3(i2,i3,i4,i1,i5,za,zb)
c      mb(2,1,2,2)=-mb(4,1,2,2)-mb(3,1,2,2)-mb(1,1,2,2)

c      mb(2,2,1,2)=-q4ghpmp1(i3,i2,i1,i4,i5,za,zb)
c      mb(3,2,1,2)=+q4ghpmp3(i4,i1,i2,i3,i5,za,zb)
c      mb(4,2,1,2)=-q4ghpmp3(i3,i2,i1,i4,i5,za,zb)
c      mb(1,2,1,2)=-mb(3,2,1,2)-mb(4,2,1,2)-mb(2,2,1,2)

c--- fastest to obtain remaining amplitudes by symmetry
      do j=1,4
      mb(j,1,1,1)=sign*dconjg(mb(j,2,2,2))
      mb(j,2,2,1)=sign*dconjg(mb(j,1,1,2))
      mb(j,1,2,2)=sign*dconjg(mb(j,2,1,1))
      mb(j,2,1,2)=sign*dconjg(mb(j,1,2,1))
      enddo
      

c      write(6,*) i1,i2,i3,i4,i5
      q4ghsq=0d0
      q4idghsq=0d0
      do hq1=1,2
      do hq2=1,2
      do hg=1,2 
c      write(6,'(a8,3i3,2f21.12)') 'mb 1 ; ',hq1,hq2,hg,mb(1,hq1,hq2,hg)
c     . -sign*dconjg(mb(1,3-hq1,3-hq2,3-hg))
c      write(6,'(a8,3i3,2f21.12)') 'mb 2 ; ',hq1,hq2,hg,mb(2,hq1,hq2,hg)
c     . -sign*dconjg(mb(2,3-hq1,3-hq2,3-hg))
c      write(6,'(a8,3i3,2f21.12)') 'mb 3 ; ',hq1,hq2,hg,mb(3,hq1,hq2,hg)
c     . -sign*dconjg(mb(3,3-hq1,3-hq2,3-hg))
c      write(6,'(a8,3i3,2f21.12)') 'mb 4 ; ',hq1,hq2,hg,mb(4,hq1,hq2,hg)
c     . -sign*dconjg(mb(4,3-hq1,3-hq2,3-hg))

      q4ghsq=q4ghsq+xn**2
     . *(cdabs(ma(1,hq1,hq2,hg)**2)+cdabs(ma(2,hq1,hq2,hg)**2))
     .  +cdabs(ma(3,hq1,hq2,hg)**2)+cdabs(ma(4,hq1,hq2,hg)**2)
     . +2d0*dble(Dconjg(+ma(1,hq1,hq2,hg)+ma(2,hq1,hq2,hg))
     .                *(+ma(3,hq1,hq2,hg)+ma(4,hq1,hq2,hg)))
      q4idghsq=q4idghsq+xn**2
     . *(+cdabs(ma(1,hq1,hq2,hg)**2)+cdabs(ma(2,hq1,hq2,hg)**2)
     .   +cdabs(mb(1,hq1,hq2,hg)**2)+cdabs(mb(2,hq1,hq2,hg)**2))
     .   +cdabs(ma(3,hq1,hq2,hg)**2)+cdabs(ma(4,hq1,hq2,hg)**2)
     .   +cdabs(mb(3,hq1,hq2,hg)**2)+cdabs(mb(4,hq1,hq2,hg)**2)
     . +2d0*dble(
     . +Dconjg(+ma(1,hq1,hq2,hg)+ma(2,hq1,hq2,hg))
     .       *(+ma(3,hq1,hq2,hg)+ma(4,hq1,hq2,hg))
     . +Dconjg(+mb(1,hq1,hq2,hg)+mb(2,hq1,hq2,hg))
     .       *(+mb(3,hq1,hq2,hg)+mb(4,hq1,hq2,hg)))

C--- Interference if the helicities are the same
c--- New version, JMC on 8/19/05; the -ve sign has already been
c--- applied in the amplitudes, so the signs here should be switched.
c--- See the text before Eq.(B24) in DFM
      if (hq1.eq.hq2) then
      q4idghsq=q4idghsq
     . +2d0*xn*dble(
     .    +Dconjg(ma(1,hq1,hq2,hg))
     .     *(mb(1,hq1,hq2,hg)+mb(2,hq1,hq2,hg)+mb(3,hq1,hq2,hg))
     .    +Dconjg(mb(1,hq1,hq2,hg))*ma(3,hq1,hq2,hg)
     .    +Dconjg(ma(2,hq1,hq2,hg))
     .     *(mb(1,hq1,hq2,hg)+mb(2,hq1,hq2,hg)+mb(4,hq1,hq2,hg))
     .    +Dconjg(mb(2,hq1,hq2,hg))*ma(4,hq1,hq2,hg))

     . +2d0/xn*dble(
     . +Dconjg(ma(3,hq1,hq2,hg)+ma(4,hq1,hq2,hg))
     .       *(mb(3,hq1,hq2,hg)+mb(4,hq1,hq2,hg)))
      endif

      enddo
      enddo
      enddo

c      pause
      
      q4ghsq=0.5d0*Cf*q4ghsq
      q4idghsq=0.5d0*Cf*q4idghsq

      return
      end
