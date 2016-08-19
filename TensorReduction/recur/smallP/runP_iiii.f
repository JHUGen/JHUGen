      subroutine runP_iiii(k,i1,i2,i3,i4,f,Gr,Shat5,N0)
      implicit none
      include 'pvDnames.f'
      include 'pvDv.f'
      include 'Darraydef.f'
      include 'Darrays.f'
      integer ep,N0,k,i1,i2,i3,i4,np
      parameter(np=3)
      double precision f(np),Gr(np,np)
      double complex Shat5(np,z4max,-2:0)
       
      do ep=-2,0
      Dv(diiii(z4(i1,i2,i3,i4))+N0,ep)=
     . (Shat5(k,z4(i1,i2,i3,i4),ep)
     . -2d0*delta(k,i1)*Dv(dzziii(z3(i2,i3,i4))+N0,ep)
     . -2d0*delta(k,i2)*Dv(dzziii(z3(i1,i3,i4))+N0,ep)
     . -2d0*delta(k,i3)*Dv(dzziii(z3(i1,i2,i4))+N0,ep)
     . -2d0*delta(k,i4)*Dv(dzziii(z3(i1,i2,i3))+N0,ep)
     . -Gr(k,1)*Dv(diiiii(z5(1,i1,i2,i3,i4))+N0,ep) 
     . -Gr(k,2)*Dv(diiiii(z5(2,i1,i2,i3,i4))+N0,ep)
     . -Gr(k,3)*Dv(diiiii(z5(3,i1,i2,i3,i4))+N0,ep))/f(k) 
      enddo
      
      return
      end
