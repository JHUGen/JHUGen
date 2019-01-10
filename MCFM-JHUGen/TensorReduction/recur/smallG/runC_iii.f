      subroutine runC_iii(j,i1,i2,i3,DetGr,Xtwiddle0,Gtwiddle,
     . Shat4,N0)
      implicit none
      include 'pvCnames.f'  
      include 'pvCv.f'  
      include 'Carraydef.f'  
      include 'Carrays.f'  
      integer ep,N0,j,i1,i2,i3,n,np
      parameter(np=2)
      double precision DetGr,Xtwiddle0(np),Gtwiddle(np,np)
      double complex Shat4(np,z3max,-2:0),Shat4s(np,z3max,-2:0)
       
      do ep=-2,0
      do n=1,np
      Shat4s(n,z3(i1,i2,i3),ep)=Shat4(n,z3(i1,i2,i3),ep)
     .   -2d0*(delta(n,i1)*Cv(N0+czzii(z2(i2,i3)),ep)
     .        +delta(n,i2)*Cv(N0+czzii(z2(i1,i3)),ep)
     .        +delta(n,i3)*Cv(N0+czzii(z2(i1,i2)),ep))
c      if ((ep.eq.-1) .and. (z3(i1,i2,i3) .eq. z3(2,2,2))) then
c      write(6,'(a8,4i3,2(f20.15))') 'Shat4s ',n,i1,i2,i3,
c     .         Shat4(n,z3(i1,i2,i3),ep)
c      write(6,'(a8,4i3,2(f20.15))') 'Shat4s ',n,i1,i2,i3,
c     .        -2d0*delta(n,i1)*Cv(N0+czzii(z2(i2,i3)),ep)
c      write(6,'(a8,4i3,2(f20.15))') 'Shat4s ',n,i1,i2,i3,
c     .        -2d0*delta(n,i2)*Cv(N0+czzii(z2(i1,i3)),ep)
c      write(6,'(a8,4i3,2(f20.15))') 'Shat4s ',n,i1,i2,i3,
c     .        -2d0*delta(n,i3)*Cv(N0+czzii(z2(i1,i2)),ep)
c      endif
      enddo 
      Cv(ciii(z3(i1,i2,i3))+N0,ep)=-(
     . +Gtwiddle(j,1)*Shat4s(1,z3(i1,i2,i3),ep)
     . +Gtwiddle(j,2)*Shat4s(2,z3(i1,i2,i3),ep)
     . -DetGr*Cv(ciiii(z4(j,i1,i2,i3))+N0,ep))/Xtwiddle0(j)

c      if (ep .eq. -1) then
c      write(6,*) 'Ciii ',j,i1,i2,i3,
c     . +Gtwiddle(j,1)*Shat4s(1,z3(i1,i2,i3),ep)/Xtwiddle0(j)
c      write(6,*) 'Ciii ',j,i1,i2,i3,
c     . +Gtwiddle(j,2),Shat4s(2,z3(i1,i2,i3),ep),Xtwiddle0(j)
c      write(6,*) 'Ciii ',j,i1,i2,i3,
c     . -DetGr*Cv(ciiii(z4(j,i1,i2,i3))+N0,ep)/Xtwiddle0(j)
c      endif

      enddo

c      pause

      return
      end
