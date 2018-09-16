      subroutine run_00ii(k,l,i1,i2,DetGr,f,Gtwiddle,Gtt,
     . Shat3,Shat4,S00ii,Shat4zz,N0)
      implicit none
C--   Fixes D00ii using 5.47
C     knowing D0000,Dii,D00i with correction of order delta*Diiii
      include 'TRconstants.f' 
      include 'pvDnames.f' 
      include 'pvDv.f' 
      include 'Darraydef.f' 
      include 'Darrays.f' 
      integer ep,N0,k,l,n,m,i1,i2,np
      parameter(np=3)
      double precision DetGr,Gtwiddle(np,np),Gtt(np,np,np,np),
     . f(np)
      double complex S00ii(z2max,-2:0),Shat4zz(np,z1max,-2:0),
     . Shat3(np,z2max,-2:0),Shat4(np,z3max,-2:0),bit,pole

       
      do ep=-2,0
      bit=czip
      do n=1,np
      do m=1,np
      bit=bit
     . +Gtt(k,n,l,m)*(f(n)*Shat3(m,z2(i1,i2),ep)
     . +2d0*(delta(n,i1)*Shat4zz(m,i2,ep)+delta(n,i2)*Shat4zz(m,i1,ep))
     . -f(n)*f(m)*Dv(dii(z2(i1,i2))+N0,ep)
     . -2d0*(f(n)*delta(m,i1)+f(m)*delta(n,i1))*Dv(dzzi(i2)+N0,ep)
     . -2d0*(f(n)*delta(m,i2)+f(m)*delta(n,i2))*Dv(dzzi(i1)+N0,ep)
     . -4d0*(delta(n,i1)*delta(m,i2)+delta(n,i2)*delta(m,i1))
     . *Dv(dd0000+N0,ep))
      enddo
      enddo
      pole=czip
      if (ep .gt. -2)
     .  pole=-4d0*Gtwiddle(k,l)*Dv(dzzii(z2(i1,i2))+N0,ep-1)
     
      Dv(dzzii(z2(i1,i2))+N0,ep)=
     . -(pole
     . +DetGr*Dv(diiii(z4(k,l,i1,i2))+N0,ep)
     . -Gtwiddle(k,l)*S00ii(z2(i1,i2),ep)
     . -Gtwiddle(1,l)*Shat4(1,z3(k,i1,i2),ep)
     . -Gtwiddle(2,l)*Shat4(2,z3(k,i1,i2),ep)
     . -Gtwiddle(3,l)*Shat4(3,z3(k,i1,i2),ep)
     . +Gtwiddle(k,l)
     . *(Shat4(1,z3(1,i1,i2),ep)
     .  +Shat4(2,z3(2,i1,i2),ep)
     .  +Shat4(3,z3(3,i1,i2),ep))
     . +bit)/(12d0*Gtwiddle(k,l))

      enddo


      return
      end
  



