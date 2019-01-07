      subroutine runC_0000i(k,l,i1,DetGr,f,Gtwiddle,Gtt, 
     . Shat4zz,Shat5zzzz,S0000i,Shat5zz,N0) 
      implicit none 
C     Fixes C0000i
C     known C00i,C0000  corrections of order Delta C00iii 
      include 'TRconstants.f'   
      include 'pvCnames.f'   
      include 'pvCv.f'   
      include 'Carraydef.f'   
      include 'Carrays.f'   
      integer ep,N0,k,l,n,m,i1,np 
      parameter(np=2) 
      double precision DetGr,Gtwiddle(np,np),Gtt(np,np,np,np), 
     . f(np) 
      double complex Shat4zz(np,z1max,-2:0),Shat5zzzz(np,-2:0), 
     . S0000i(np,-2:0),Shat5zz(np,z2max,-2:0),bit,pole
        
      do ep=-2,0
      bit=czip 
      do n=1,np 
      do m=1,np 
      bit=bit 
     . +Gtt(k,n,l,m)*(f(n)*Shat4zz(m,i1,ep) 
     . +2d0*(delta(n,i1)*Shat5zzzz(m,ep)) 
     . -f(n)*f(m)*Cv(czzi(i1)+N0,ep) 
     . -2d0*(f(n)*delta(m,i1)+f(m)*delta(n,i1))*Cv(cc0000+N0,ep))
      enddo 
      enddo 
      pole=czip 
      if (ep .gt. -2)  
     . pole=-4d0*Gtwiddle(k,l)*Cv(czzzzi(i1)+N0,ep-1) 
      Cv(czzzzi(i1)+N0,ep)=
     . -(pole 
     . +DetGr*Cv(czziii(z3(k,l,i1))+N0,ep) 
     . -Gtwiddle(k,l)*S0000i(i1,ep) 
     . -Gtwiddle(1,l)*Shat5zz(1,z2(k,i1),ep) 
     . -Gtwiddle(2,l)*Shat5zz(2,z2(k,i1),ep) 
     . +Gtwiddle(k,l) 
     . *(Shat5zz(1,z2(1,i1),ep) 
     .  +Shat5zz(2,z2(2,i1),ep)) 
     . +bit)/(14d0*Gtwiddle(k,l))

c      if (ep .eq. -1) then 
c      write(6,*) 'C0000i ',ep,k,l,i1,
c     . pole 
c      write(6,*) 'C0000i ',ep,k,l,i1,
c     . +DetGr*Cv(czziii(z3(k,l,i1))+N0,ep) 
c      write(6,*) 'C0000i ',ep,k,l,i1,
c     . -Gtwiddle(k,l)*S0000i(i1,ep) 
c      write(6,*) 'C0000i ',ep,k,l,i1,
c     . -Gtwiddle(1,l)*Shat5zz(1,z2(k,i1),ep) 
c      write(6,*) 'C0000i ',ep,k,l,i1,
c     . -Gtwiddle(2,l)*Shat5zz(2,z2(k,i1),ep) 
c      write(6,*) 'C0000i ',ep,k,l,i1,
c     . +Gtwiddle(k,l) 
c     . *(Shat5zz(1,z2(1,i1),ep) 
c     .  +Shat5zz(2,z2(2,i1),ep)) 
c      write(6,*) 'C0000i ',ep,k,l,i1,
c     . +bit
c      endif
      
      enddo 
      return
      end
