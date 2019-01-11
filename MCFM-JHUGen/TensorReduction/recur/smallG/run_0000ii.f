      subroutine run_0000ii(k,l,i1,i2,DetGr,f,Gtwiddle,Gtt, 
     . Shat5zz,Shat6zzzz,S0000ii,Shat6zz,N0) 
      implicit none 
C------Fixes D0000ii 
C-----knowing D000000 with correction of order Delta D0000ii
      include 'TRconstants.f'   
      include 'pvDnames.f'   
      include 'pvDv.f'   
      include 'Darraydef.f'   
      include 'Darrays.f'   
      integer ep,N0,k,l,n,m,i1,i2,np 
      parameter(np=3) 
      double precision DetGr,Gtwiddle(np,np),Gtt(np,np,np,np), 
     . f(np) 
      double complex Shat5zz(np,z2max,-2:0),Shat6zzzz(np,z1max,-2:0), 
     . S0000ii(z2max,-2:0),Shat6zz(np,z3max,-2:0),bit,pole
        
      do ep=-2,0
      bit=czip 
      do n=1,np 
      do m=1,np 
      bit=bit 
     . +Gtt(k,n,l,m)*(f(n)*Shat5zz(m,z2(i1,i2),ep) 
     . +2*(delta(n,i1)*Shat6zzzz(m,i2,ep)
     .    +delta(n,i2)*Shat6zzzz(m,i1,ep)) 
     . -f(n)*f(m)*Dv(dzzii(z2(i1,i2))+N0,ep) 
     . -2*(f(n)*delta(m,i1)+f(m)*delta(n,i1))*Dv(dzzzzi(i2)+N0,ep)
     . -2*(f(n)*delta(m,i2)+f(m)*delta(n,i2))*Dv(dzzzzi(i1)+N0,ep)
     . -4*(delta(n,i1)*delta(m,i2)+delta(m,i1)*delta(n,i2))
     . *Dv(dd000000+N0,ep))
      enddo 
      enddo 

      pole=czip 
      if (ep .gt. -2)  
     . pole=-4*Gtwiddle(k,l)*Dv(dzzzzii(z2(i1,i2))+N0,ep-1) 
      Dv(dzzzzii(z2(i1,i2))+N0,ep)=
     . -(pole 
     . +DetGr*Dv(dzziiii(z4(k,l,i1,i2))+N0,ep) 
     . -Gtwiddle(k,l)*S0000ii(z2(i1,i2),ep) 
     . -Gtwiddle(1,l)*Shat6zz(1,z3(k,i1,i2),ep) 
     . -Gtwiddle(2,l)*Shat6zz(2,z3(k,i1,i2),ep) 
     . -Gtwiddle(3,l)*Shat6zz(3,z3(k,i1,i2),ep) 
     . +Gtwiddle(k,l) 
     . *(Shat6zz(1,z3(1,i1,i2),ep) 
     .  +Shat6zz(2,z3(2,i1,i2),ep) 
     .  +Shat6zz(3,z3(3,i1,i2),ep)) 
     . +bit)/(16d0*Gtwiddle(k,l))
      enddo 
      return
      end
