      subroutine run_00iii(k,l,i1,i2,i3,DetGr,f,Gtwiddle,Gtt,
     . Shat4,Shat5,S00iii,Shat5zz,N0)
      implicit none
C----fixes D00iii
c     known Diii,D00ii,D0000i with corrections of order Delta Diiiii
      include 'TRconstants.f'  
      include 'pvDnames.f'  
      include 'pvDv.f'  
      include 'Darraydef.f'  
      include 'Darrays.f'  
      integer ep,N0,k,l,n,m,i1,i2,i3,np
      parameter(np=3)
      double precision DetGr,Gtwiddle(np,np),Gtt(np,np,np,np),
     . f(np)
      double complex S00iii(z3max,-2:0),Shat5zz(np,z2max,-2:0),
     . Shat4(np,z3max,-2:0),Shat5(np,z4max,-2:0),bit,pole

       
      do ep=-2,0
      bit=czip
      do n=1,np
      do m=1,np
      bit=bit
     . +Gtt(k,n,l,m)*(f(n)*Shat4(m,z3(i1,i2,i3),ep)
     . +2*(delta(n,i1)*Shat5zz(m,z2(i2,i3),ep)
     .    +delta(n,i2)*Shat5zz(m,z2(i1,i3),ep)
     .    +delta(n,i3)*Shat5zz(m,z2(i1,i2),ep))
     . -f(n)*f(m)*Dv(diii(z3(i1,i2,i3))+N0,ep)
     . -2*(f(n)*delta(m,i1)+f(m)*delta(n,i1))*Dv(dzzii(z2(i2,i3))+N0,ep)
     . -2*(f(n)*delta(m,i2)+f(m)*delta(n,i2))*Dv(dzzii(z2(i1,i3))+N0,ep)
     . -2*(f(n)*delta(m,i3)+f(m)*delta(n,i3))*Dv(dzzii(z2(i1,i2))+N0,ep)
     . -4*(delta(n,i1)*delta(m,i2)+delta(n,i2)*delta(m,i1))
     . *Dv(dzzzzi(i3)+N0,ep)
     . -4*(delta(n,i2)*delta(m,i3)+delta(n,i3)*delta(m,i2))
     . *Dv(dzzzzi(i1)+N0,ep)
     . -4*(delta(n,i3)*delta(m,i1)+delta(n,i1)*delta(m,i3))
     . *Dv(dzzzzi(i2)+N0,ep))
      enddo
      enddo
      pole=czip
      if (ep .gt. -2) 
     . pole=-4*Gtwiddle(k,l)*Dv(dzziii(z3(i1,i2,i3))+N0,ep-1)

      Dv(dzziii(z3(i1,i2,i3))+N0,ep)=
     . -(pole
     . +DetGr*Dv(diiiii(z5(k,l,i1,i2,i3))+N0,ep)
     . -Gtwiddle(k,l)*S00iii(z3(i1,i2,i3),ep)
     . -Gtwiddle(1,l)*Shat5(1,z4(k,i1,i2,i3),ep)
     . -Gtwiddle(2,l)*Shat5(2,z4(k,i1,i2,i3),ep)
     . -Gtwiddle(3,l)*Shat5(3,z4(k,i1,i2,i3),ep)
     . +Gtwiddle(k,l)
     . *(Shat5(1,z4(1,i1,i2,i3),ep)
     .  +Shat5(2,z4(2,i1,i2,i3),ep)
     .  +Shat5(3,z4(3,i1,i2,i3),ep))
     . +bit)/(16d0*Gtwiddle(k,l))

      enddo


      return
      end
  



