      subroutine runC_iiii(j,i1,i2,i3,i4,DetGr,Xtwiddle0,Gtwiddle,
     . Shat5,N0)
      implicit none
      include 'pvCnames.f'  
      include 'pvCv.f'  
      include 'Carraydef.f'  
      include 'Carrays.f'  
      integer ep,N0,j,i1,i2,i3,i4,n,np
      parameter(np=2)
      double precision DetGr,Xtwiddle0(np),Gtwiddle(np,np)
      double complex Shat5(np,z4max,-2:0),Shat5s(np,z4max,-2:0)
       

      do ep=-2,0
      do n=1,np
      Shat5s(n,z4(i1,i2,i3,i4),ep)=Shat5(n,z4(i1,i2,i3,i4),ep)
     .   -2d0*(delta(n,i1)*Cv(N0+czziii(z3(i2,i3,i4)),ep)
     .        +delta(n,i2)*Cv(N0+czziii(z3(i1,i3,i4)),ep)
     .        +delta(n,i3)*Cv(N0+czziii(z3(i1,i2,i4)),ep)
     .        +delta(n,i4)*Cv(N0+czziii(z3(i1,i2,i3)),ep))
c      if (ep.eq.-1) then
c      write(6,*) 'Shat5s ',n,i1,i2,i3,i4,
c     .         Shat5(n,z4(i1,i2,i3,i4),-1)
c      write(6,*) 'Shat5s ',n,i1,i2,i3,i4,
c     .         -2d0*delta(n,i1)*Cv(N0+czziii(z3(i2,i3,i4)),-1)
c      write(6,*) 'Shat5s ',n,i1,i2,i3,i4,
c     .        -2d0*delta(n,i2)*Cv(N0+czziii(z3(i1,i3,i4)),-1)
c      write(6,*) 'Shat5s ',n,i1,i2,i3,i4,
c     .        -2d0*delta(n,i3)*Cv(N0+czziii(z3(i1,i2,i4)),-1)
c      write(6,*) 'Shat5s ',n,i1,i2,i3,i4,
c     .        -2d0*delta(n,i4)*Cv(N0+czziii(z3(i1,i2,i3)),-1)
c      endif
      enddo 

      Cv(ciiii(z4(i1,i2,i3,i4))+N0,ep)=-(
     . +Gtwiddle(j,1)*Shat5s(1,z4(i1,i2,i3,i4),ep)
     . +Gtwiddle(j,2)*Shat5s(2,z4(i1,i2,i3,i4),ep)
     . -DetGr*Cv(ciiiii(z5(j,i1,i2,i3,i4))+N0,ep))/Xtwiddle0(j)
     
c      write(6,*) 'Ciiii ',i1,i2,i3,i4,ep,
c     . +Gtwiddle(j,1)*Shat5s(1,z4(i1,i2,i3,i4),ep)/Xtwiddle0(j)
c      write(6,*) 'Ciiii ',i1,i2,i3,i4,ep,
c     . +Gtwiddle(j,2)*Shat5s(1,z4(i1,i2,i3,i4),ep)/Xtwiddle0(j)
c      write(6,*) 'Ciiii ',i1,i2,i3,i4,ep,
c     . -DetGr*Cv(ciiiii(z5(j,i1,i2,i3,i4))+N0,ep)/Xtwiddle0(j)
     
      enddo

c      pause

      return
      end
