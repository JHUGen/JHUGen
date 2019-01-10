      subroutine runC_i(j,i1,DetGr,Xtwiddle0,Gtwiddle,Shat2,N0)
      implicit none
      include 'TRconstants.f'  
      include 'pvCnames.f'  
      include 'pvCv.f'  
      include 'Carraydef.f'  
      include 'Carrays.f'  
      integer ep,N0,j,i1,n,np
      parameter(np=2)
      double precision DetGr,Xtwiddle0(np),Gtwiddle(np,np)
      double complex Shat2(np,np,-2:0),Shat2s(np,np,-2:0),bit
       
      do ep=-2,0
      bit=czip
      do n=1,np
      Shat2s(n,i1,ep)=Shat2(n,i1,ep)-2d0*delta(n,i1)*Cv(N0+cc00,ep)
      bit=bit+Gtwiddle(j,n)*Shat2s(n,i1,ep)
      enddo
      Cv(ci(i1)+N0,ep)=
     . -(bit-DetGr*Cv(cii(z2(j,i1))+N0,ep))/Xtwiddle0(j)

      enddo
      
      return
      end
