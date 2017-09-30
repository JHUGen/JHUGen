      subroutine runC_00iiii(k,l,i1,i2,i3,i4,DetGr,f,Gtwiddle,Gtt,
     . Shat5,Shat6,S00iiii,Shat6zz,N0)
C------Fixes C00iiii 
C-----knowing C00iii,C0000ii with correction of order Delta Ciiiiii
      implicit none
      include 'TRconstants.f'  
      include 'pvCnames.f'  
      include 'pvCv.f'  
      include 'Carraydef.f'  
      include 'Carrays.f'  
      integer ep,N0,k,l,n,m,i1,i2,i3,i4,np
      parameter(np=2)
      double precision DetGr,Gtwiddle(np,np),Gtt(np,np,np,np),f(np)
      double complex S00iiii(z4max,-2:0),Shat6zz(np,z3max,-2:0),
     . Shat5(np,z4max,-2:0),Shat6(np,z5max,-2:0),bit,pole

       
      do ep=-2,0
      bit=czip
      do n=1,np
      do m=1,np
      bit=bit
     . +Gtt(k,n,l,m)*(f(n)*Shat5(m,z4(i1,i2,i3,i4),ep)
     . +2*(delta(n,i1)*Shat6zz(m,z3(i2,i3,i4),ep)
     .    +delta(n,i2)*Shat6zz(m,z3(i1,i3,i4),ep)
     .    +delta(n,i3)*Shat6zz(m,z3(i1,i2,i4),ep)
     .    +delta(n,i4)*Shat6zz(m,z3(i1,i2,i3),ep))
     . -f(n)*f(m)*Cv(ciiii(z4(i1,i2,i3,i4))+N0,ep)
     . -2*(f(n)*delta(m,i1)+f(m)*delta(n,i1))
     . *Cv(czziii(z3(i2,i3,i4))+N0,ep)
     . -2*(f(n)*delta(m,i2)+f(m)*delta(n,i2))
     . *Cv(czziii(z3(i1,i3,i4))+N0,ep)
     . -2*(f(n)*delta(m,i3)+f(m)*delta(n,i3))
     . *Cv(czziii(z3(i1,i2,i4))+N0,ep)
     . -2*(f(n)*delta(m,i4)+f(m)*delta(n,i4))
     . *Cv(czziii(z3(i1,i2,i3))+N0,ep)

     . -4*(delta(n,i1)*delta(m,i2)+delta(n,i2)*delta(m,i1))
     . *Cv(czzzzii(z2(i3,i4))+N0,ep)
     . -4*(delta(n,i1)*delta(m,i3)+delta(n,i3)*delta(m,i1))
     . *Cv(czzzzii(z2(i2,i4))+N0,ep)
     . -4*(delta(n,i1)*delta(m,i4)+delta(n,i4)*delta(m,i1))
     . *Cv(czzzzii(z2(i2,i3))+N0,ep)
     . -4*(delta(n,i2)*delta(m,i3)+delta(n,i3)*delta(m,i2))
     . *Cv(czzzzii(z2(i1,i4))+N0,ep)
     . -4*(delta(n,i2)*delta(m,i4)+delta(n,i4)*delta(m,i2))
     . *Cv(czzzzii(z2(i1,i3))+N0,ep)
     . -4*(delta(n,i3)*delta(m,i4)+delta(n,i4)*delta(m,i3))
     . *Cv(czzzzii(z2(i1,i2))+N0,ep))
      enddo
      enddo
      pole=czip
      if (ep .gt. -2) 
     . pole=-4*Gtwiddle(k,l)*Cv(czziiii(z4(i1,i2,i3,i4))+N0,ep-1)
      Cv(czziiii(z4(i1,i2,i3,i4))+N0,ep)=
     . -(pole
     . +DetGr*Cv(ciiiiii(z6(k,l,i1,i2,i3,i4))+N0,ep)
     . -Gtwiddle(k,l)*S00iiii(z4(i1,i2,i3,i4),ep)
     . -Gtwiddle(1,l)*Shat6(1,z5(k,i1,i2,i3,i4),ep)
     . -Gtwiddle(2,l)*Shat6(2,z5(k,i1,i2,i3,i4),ep)
     . +Gtwiddle(k,l)
     . *(Shat6(1,z5(1,i1,i2,i3,i4),ep)
     .  +Shat6(2,z5(2,i1,i2,i3,i4),ep))
     . +bit)/(22d0*Gtwiddle(k,l))

      enddo


      return
      end
  



