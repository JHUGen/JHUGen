       subroutine runC_0000(k,l,DetGr,f,Gtwiddle,Gtt,
     . Shat3zz,Shat4zz,S0000,N0)
      implicit none
      include 'TRconstants.f'  
      include 'pvCnames.f'  
      include 'pvCv.f'  
      include 'Carraydef.f'  
      include 'Carrays.f'  
      integer ep,N0,k,l,n,m,np
      parameter(np=2)
      double precision DetGr,Gtwiddle(np,np),Gtt(np,np,np,np),f(np)
      double complex Shat3zz(np,-2:0),S0000(-2:0),
     . Shat4zz(np,z1max,-2:0),bit,pole

       
      do ep=-2,0
      bit=czip
      do n=1,np
      do m=1,np
      bit=bit
     . +Gtt(k,n,l,m)*(f(n)*Shat3zz(m,ep)-f(n)*f(m)*Cv(cc00+N0,ep))
      enddo
      enddo
      pole=czip
      if (ep .gt. -2) pole=-4d0*Gtwiddle(k,l)*Cv(cc0000+N0,ep-1)

      Cv(cc0000+N0,ep)=-(
     . +pole
     . +DetGr*Cv(czzii(z2(k,l))+N0,ep)
     . -Gtwiddle(k,l)*S0000(ep)
     . -Gtwiddle(1,l)*Shat4zz(1,k,ep)
     . -Gtwiddle(2,l)*Shat4zz(2,k,ep)
     . +Gtwiddle(k,l)*(Shat4zz(1,1,ep)+Shat4zz(2,2,ep))
     . +bit)/(10d0*Gtwiddle(k,l))

c      if (ep .eq. -1) then
c      write(6,*) 'C0000 ',k,l,
c     . +pole/(10d0*Gtwiddle(k,l))
c      write(6,*) 'C0000 ',k,l,
c     . +DetGr*Cv(czzii(z2(k,l))+N0,ep)/(10d0*Gtwiddle(k,l))
c      write(6,*) 'C0000 ',k,l,
c     . -Gtwiddle(k,l)*S0000(ep)/(10d0*Gtwiddle(k,l))
c      write(6,*) 'C0000 ',k,l,
c     . -Gtwiddle(1,l)*Shat4zz(1,k,ep)/(10d0*Gtwiddle(k,l))
c      write(6,*) 'C0000 ',k,l,
c     . -Gtwiddle(2,l)*Shat4zz(2,k,ep)/(10d0*Gtwiddle(k,l))
c      write(6,*) 'C0000 ',k,l,
c     . +Gtwiddle(k,l)*(Shat4zz(1,1,ep)+Shat4zz(2,2,ep))
c     .  /(10d0*Gtwiddle(k,l))
c      write(6,*) 'C0000 ',k,l,
c     . +bit/(10d0*Gtwiddle(k,l))
c      endif
      enddo

      return
      end
  



