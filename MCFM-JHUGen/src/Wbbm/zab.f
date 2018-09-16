      function zab(j1,p1,j2)
      implicit none
      include 'types.f'
      complex(dp):: zab
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'momwbbm.f'
      integer:: j1,j2

      real(dp):: p1(4),E,px,py,pz
      complex(dp):: Um(4),Vp(4),zp1(4),spstrng1,rtpp

C---Calculate spinor 1
      E=mom(j1,4)      
      px=+mom(j1,3)      
      py=-mom(j1,2)      
      pz=+mom(j1,1)      
      rtpp=sqrt(cplx1(E+pz))
      Um(1)=cplx2(px,+py)/rtpp
      Um(2)=-rtpp
      Um(3)=czip
      Um(4)=czip

C---Calculate spinor 2
      E=mom(j2,4)      
      px=+mom(j2,3)      
      py=-mom(j2,2)      
      pz=+mom(j2,1)      
      rtpp=sqrt(cplx1(E+pz))
      Vp(1)=czip
      Vp(2)=czip
      Vp(3)=cplx2(px,-py)/rtpp
      Vp(4)=-rtpp

      zp1(1)=cplx1(p1(1))
      zp1(2)=cplx1(p1(2))
      zp1(3)=cplx1(p1(3))
      zp1(4)=cplx1(p1(4))

      zab=spstrng1(Um,zp1,Vp)
      return
      end
