      double complex function zabab(j1,p1,p2,p3,j2)
      implicit none
      include 'constants.f'
      include 'momwbbm.f'
      integer j1,j2

      double precision p1(4),p2(4),p3(4),E,px,py,pz
      double complex Um(4),Vp(4),zp1(4),zp2(4),zp3(4),spstrng3,rtpp

C---Calculate spinor 1
      E=mom(j1,4)      
      px=+mom(j1,3)      
      py=-mom(j1,2)      
      pz=+mom(j1,1)      
      rtpp=sqrt(dcmplx(E+pz))
      Um(1)=dcmplx(px,+py)/rtpp
      Um(2)=-rtpp
      Um(3)=czip
      Um(4)=czip

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

      zp3(1)=dcmplx(p3(1))
      zp3(2)=dcmplx(p3(2))
      zp3(3)=dcmplx(p3(3))
      zp3(4)=dcmplx(p3(4))

      zabab=spstrng3(Um,zp1,zp2,zp3,Vp)
      
      return
      end
