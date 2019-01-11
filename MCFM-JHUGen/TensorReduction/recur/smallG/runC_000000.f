      subroutine runC_000000(k,l,DetGr,f,Gtwiddle,Gtt,
     . Shat5zzzz,Shat6zzzz,S000000,N0)
      implicit none
C------Fixes C000000 
C-----knowing C0000 with correction of order Delta C0000ii
      include 'TRconstants.f'  
      include 'pvCnames.f'  
      include 'pvCv.f'  
      include 'Carraydef.f'  
      include 'Carrays.f'  
      integer ep,N0,k,l,n,m,np
      parameter(np=2)
      double precision DetGr,Gtwiddle(np,np),Gtt(np,np,np,np),f(np)
      double complex Shat5zzzz(np,-2:0),S000000(-2:0),
     . Shat6zzzz(np,z1max,-2:0),bit,pole

       
      do ep=-2,0
      bit=czip
      do n=1,np
      do m=1,np
      bit=bit
     . +Gtt(k,n,l,m)*(f(n)*Shat5zzzz(m,ep)-f(n)*f(m)*Cv(cc0000+N0,ep))
      enddo
      enddo
      pole=czip
      if (ep .gt. -2) pole=-4d0*Gtwiddle(k,l)*Cv(cc000000+N0,ep-1)

      Cv(cc000000+N0,ep)=
     . -(pole
     . +DetGr*Cv(czzzzii(z2(k,l))+N0,ep)
     . -Gtwiddle(k,l)*S000000(ep)
     . -Gtwiddle(1,l)*Shat6zzzz(1,k,ep)
     . -Gtwiddle(2,l)*Shat6zzzz(2,k,ep)
     . +Gtwiddle(k,l)
     . *(Shat6zzzz(1,1,ep)+Shat6zzzz(2,2,ep))
     . +bit)/(14d0*Gtwiddle(k,l))

      enddo


      return
      end
  



