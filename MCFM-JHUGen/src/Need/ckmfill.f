      subroutine ckmfill(nwz)
      implicit none
      include 'types.f'

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ckm.f'
      include 'ckm1.f'
      include 'ewcouple.f'
      include 'zcouple.f'
      integer:: nwz,j,k
      real(dp):: Vud,Vus,Vub,Vcd,Vcs,Vcb,rtxw
      common/cabib/Vud,Vus,Vub,
     &             Vcd,Vcs,Vcb

c---- initialize Vsq
      do j=-nf,nf
      do k=-nf,nf
      Vsq(j,k)=0._dp
      enddo
      enddo
c set all other couplings (for W+2 jets) to zero too
      flsq=0._dp
      frsq=0._dp
      fl=0._dp
      fr=0._dp

      do j=-nf,nf
      do k=-nf,nf
      gl(j,k)=0._dp
      gr(j,k)=0._dp
      glsq(j,k)=0._dp
      grsq(j,k)=0._dp
      enddo
      enddo

C case Z0
      if (nwz == 0) then
      Vsq(1,-1)=1._dp
      Vsq(2,-2)=1._dp
      Vsq(3,-3)=1._dp
      Vsq(4,-4)=1._dp
      Vsq(5,-5)=1._dp
      rtxw=sqrt(xw)
      do j=1,5
        gl(+j,-j)=l(j)*gw*rtxw
        gr(+j,-j)=r(j)*gw*rtxw
      enddo
      fl=le*gw*rtxw
      fr=re*gw*rtxw
      do j=1,nf
        glsq(j,-j)=gl(j,-j)**2
        grsq(j,-j)=gr(j,-j)**2
      enddo
      flsq=fl**2
      frsq=fr**2
      return

C case W+
      elseif (nwz == 1) then

      Vsq(2,-1)=Vud**2
      Vsq(2,-3)=Vus**2
      Vsq(2,-5)=Vub**2
      Vsq(4,-1)=Vcd**2
      Vsq(4,-3)=Vcs**2
      Vsq(4,-5)=Vcb**2

      glsq(2,-1)=gwsq/2._dp
      glsq(2,-3)=gwsq/2._dp
      glsq(2,-5)=gwsq/2._dp
      glsq(4,-1)=gwsq/2._dp
      glsq(4,-3)=gwsq/2._dp
      glsq(4,-5)=gwsq/2._dp
      flsq=gwsq/2._dp
      fl=sqrt(flsq)


C case W-
      elseif (nwz == -1) then

      Vsq(1,-2)=Vud**2
      Vsq(3,-2)=Vus**2
      Vsq(5,-2)=Vub**2
      Vsq(1,-4)=Vcd**2
      Vsq(3,-4)=Vcs**2
      Vsq(5,-4)=Vcb**2

      glsq(1,-2)=gwsq/2._dp
      glsq(3,-2)=gwsq/2._dp
      glsq(5,-2)=gwsq/2._dp
      glsq(1,-4)=gwsq/2._dp
      glsq(3,-4)=gwsq/2._dp
      glsq(5,-4)=gwsq/2._dp
      flsq=gwsq/2._dp
      fl=sqrt(flsq)

C case (W+ + W-)
      elseif (nwz == 2) then

C case W+
      Vsq(2,-1)=Vud**2
      Vsq(2,-3)=Vus**2
      Vsq(2,-5)=Vub**2
      Vsq(4,-1)=Vcd**2
      Vsq(4,-3)=Vcs**2
      Vsq(4,-5)=Vcb**2
      glsq(2,-1)=gwsq/2._dp
      glsq(2,-3)=gwsq/2._dp
      glsq(2,-5)=gwsq/2._dp
      glsq(4,-1)=gwsq/2._dp
      glsq(4,-3)=gwsq/2._dp
      glsq(4,-5)=gwsq/2._dp
C case W-
      Vsq(1,-2)=Vud**2
      Vsq(3,-2)=Vus**2
      Vsq(5,-2)=Vub**2
      Vsq(1,-4)=Vcd**2
      Vsq(3,-4)=Vcs**2
      Vsq(5,-4)=Vcb**2
      glsq(1,-2)=gwsq/2._dp
      glsq(3,-2)=gwsq/2._dp
      glsq(5,-2)=gwsq/2._dp
      glsq(1,-4)=gwsq/2._dp
      glsq(3,-4)=gwsq/2._dp
      glsq(5,-4)=gwsq/2._dp

      flsq=gwsq/2._dp
      fl=sqrt(flsq)

      endif

      do j=-nf,nf
      do k=-nf,nf
      Vsq(j,k)=Vsq(k,j)
      glsq(j,k)=glsq(k,j)
      grsq(j,k)=grsq(k,j)
      VV(j,k)=sqrt(Vsq(j,k))
      gl(j,k)=sqrt(glsq(j,k))
      gr(j,k)=sqrt(grsq(j,k))
      enddo
      enddo
      do j=1,5
      Vsum(+j)=Vsq(+j,-1)+Vsq(+j,-2)+Vsq(+j,-3)+Vsq(+j,-4)+Vsq(+j,-5)
      Vsum(-j)=Vsq(-j,+1)+Vsq(-j,+2)+Vsq(-j,+3)+Vsq(-j,+4)+Vsq(-j,+5)
      enddo
      Vsum(0)=0
      return
      end
