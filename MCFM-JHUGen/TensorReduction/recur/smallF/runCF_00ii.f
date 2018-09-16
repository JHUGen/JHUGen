      subroutine runCF_00ii(i1,i2,f,Gr,Shat4,N0)
C---  Expression for rearrangement of Eq. 5.70
C---  Calculates C00ii
C---  Small terms of order f(i)*Cijk,Gr(i,j)*Cijkl
      implicit none
      include 'pvCnames.f'
      include 'pvCv.f'
      include 'Carraydef.f'
      include 'Carrays.f'
      include 'pvweenumber.f' 
      integer ep,N0,i1,i2,np
      parameter(np=2)
      double precision f(np),Gr(np,np),den
      double complex Shat4(np,z3max,-2:0)
       
      do ep=-2,0
      if (i1 .eq. i2) then
        den=6d0
      else
        den=4d0
      endif      
      Cv(czzii(z2(i1,i2))+N0,ep)=
     . (Shat4(i1,z3(i1,i1,i2),ep)
     . -f(i1)*Cv(ciii(z3(i1,i1,i2))+N0,ep)
     . -Gr(i1,1)*Cv(ciiii(z4(1,i1,i1,i2))+N0,ep) 
     . -Gr(i1,2)*Cv(ciiii(z4(2,i1,i1,i2))+N0,ep))/den
     
      enddo
      
      return
      end
