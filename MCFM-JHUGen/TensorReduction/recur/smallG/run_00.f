      subroutine run_00(k,l,DetGr,f,Gtwiddle,Gtt,
     . Shat1,Shat2,S00,N0)
      implicit none
C---Fixes D00 according to 5.42 Denner-Dittmaier
C---knowing D0 with corrections of order Delta_3 Dii
      include 'TRconstants.f'
      include 'pvDnames.f'
      include 'pvDv.f'
      include 'Darraydef.f'
      include 'Darrays.f'
      integer ep,N0,k,l,n,m
      double precision DetGr,Gtwiddle(3,3),Gtt(3,3,3,3),
     . f(3)
      double complex S00(-2:0),Shat1(3,-2:0),Shat2(3,3,-2:0),bit

      do ep=-2,0
      bit=czip
      do n=1,3
      do m=1,3
      bit=bit
     . +Gtt(k,n,l,m)*(f(n)*Shat1(m,ep)-f(n)*f(m)*Dv(dd0+N0,ep))
      enddo
      enddo
      Dv(dd00+N0,ep)=-(
     . +DetGr*Dv(dii(z2(k,l))+N0,ep)
     . -Gtwiddle(k,l)*S00(ep)
     . -Gtwiddle(1,l)*Shat2(1,k,ep)
     . -Gtwiddle(2,l)*Shat2(2,k,ep)
     . -Gtwiddle(3,l)*Shat2(3,k,ep)
     . +Gtwiddle(k,l)*Shat2(1,1,ep)
     . +Gtwiddle(k,l)*Shat2(2,2,ep)
     . +Gtwiddle(k,l)*Shat2(3,3,ep)
     . +bit)
     . /(4d0*Gtwiddle(k,l))
      enddo


      return
      end
  



