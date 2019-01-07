      double complex function zbab(j1,p1,p2,j2)
      implicit none
      include 'constants.f'
      include 'momwbbm.f'
      integer j1,j2

      double precision p1(4),p2(4),E,px,py,pz
      double complex Up(4),Vp(4),zp1(4),zp2(4),spstrng2,rtpp

C---Calculate spinor 1
      E=mom(j1,4)      
      px=+mom(j1,3)      
      py=-mom(j1,2)      
      pz=+mom(j1,1)      
      rtpp=sqrt(dcmplx(E+pz))
      Up(1)=czip
      Up(2)=czip
      Up(3)=rtpp
      Up(4)=dcmplx(px,-py)/rtpp

C---Calculate spinor 2
      E=mom(j2,4)      
      px=+mom(j2,3)      
      py=-mom(j2,2)      
      pz=+mom(j2,1)      
      rtpp=sqrt(dcmplx(E+pz))
      Vp(1)=czip
      Vp(2)=czip
      Vp(3)=dcmplx(px,-py)/rtpp
      Vp(4)=-rtpp

      zp1(1)=dcmplx(p1(1))
      zp1(2)=dcmplx(p1(2))
      zp1(3)=dcmplx(p1(3))
      zp1(4)=dcmplx(p1(4))

      zp2(1)=dcmplx(p2(1))
      zp2(2)=dcmplx(p2(2))
      zp2(3)=dcmplx(p2(3))
      zp2(4)=dcmplx(p2(4))

      zbab=spstrng2(Up,zp1,zp2,Vp)
      return
      end
