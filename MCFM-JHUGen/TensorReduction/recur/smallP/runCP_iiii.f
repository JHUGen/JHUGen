      subroutine runCP_iiii(k,i1,i2,i3,i4,f,Gr,Shat5,N0)
      implicit none
      include 'pvCnames.f'
      include 'pvCv.f'
      include 'Carraydef.f'
      include 'Carrays.f'
      integer ep,N0,k,i1,i2,i3,i4,np
      parameter(np=2)
      double precision f(np),Gr(np,np)
      double complex Shat5(np,z4max,-2:0)
       
      do ep=-2,0
      Cv(ciiii(z4(i1,i2,i3,i4))+N0,ep)=
     . (Shat5(k,z4(i1,i2,i3,i4),ep)
     . -2d0*delta(k,i1)*Cv(czziii(z3(i2,i3,i4))+N0,ep)
     . -2d0*delta(k,i2)*Cv(czziii(z3(i1,i3,i4))+N0,ep)
     . -2d0*delta(k,i3)*Cv(czziii(z3(i1,i2,i4))+N0,ep)
     . -2d0*delta(k,i4)*Cv(czziii(z3(i1,i2,i3))+N0,ep)
     . -Gr(k,1)*Cv(ciiiii(z5(1,i1,i2,i3,i4))+N0,ep) 
     . -Gr(k,2)*Cv(ciiiii(z5(2,i1,i2,i3,i4))+N0,ep))/f(k) 
      enddo
      
      return
      end
