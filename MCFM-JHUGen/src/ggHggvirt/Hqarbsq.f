      double precision function Hqarbsq(i1,i2,i3,i4)
      implicit none
C----Matrix element squared for the process
C----q(p1)+r(p3) -> q(p2)+r(p4)+Higgs
C----summed over incoming/outgoing colors and spins.
C    with a factor of gsq*Asq removed.

c       +V*((s13*s24-s23*s14)**2+s12**2*s34**2)/s34**2/s12**2
c       +1/2*V*((s13-s24)**2+(s14-s23)**2)/s34/s12
c       -1/2*e*V*(s13+s14+s23+s24)**2/s34/s12;

      include 'constants.f'
      include 'sprods_com.f'
      integer i1,i2,i3,i4
      double precision s12,s13,s14,s23,s24,s34
      s12=s(i1,i2)
      s13=s(i1,i3)
      s14=s(i1,i4)
      s23=s(i2,i3)
      s24=s(i2,i4)
      s34=s(i3,i4)
C GZ 
!      write(*,*) 'in Hqarbsq'
!       s12 = 1.000000000000000d0
!       s13= -0.2861267766967376d0
!       s14= -0.3784302994364646d0
!       s23= -0.6215181025806800d0
!       s24= -0.3421507931851938d0
!       s34=  0.6319744100498450d0
 
      Hqarbsq=V*(((s13*s24-s23*s14)**2+s12**2*s34**2)/s34**2/s12**2
     . +0.5d0*((s13-s24)**2+(s14-s23)**2)/s34/s12)
c      write(*,*) 'In Hqarbsq',Hqarbsq
      return
      end
