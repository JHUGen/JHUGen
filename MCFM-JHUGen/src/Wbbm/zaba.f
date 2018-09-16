      function zaba(j1,p1,p2,j2)
      implicit none
      include 'types.f'
      complex(dp):: zaba
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'momwbbm.f'
      integer:: j1,j2

      real(dp):: p1(4),p2(4),E,px,py,pz
      complex(dp):: Um(4),Vm(4),zp1(4),zp2(4),spstrng2,rtpp

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
      Vm(1)=rtpp
      Vm(2)=cplx2(px,py)/rtpp
      Vm(3)=czip
      Vm(4)=czip

      zp1(1)=cplx1(p1(1))
      zp1(2)=cplx1(p1(2))
      zp1(3)=cplx1(p1(3))
      zp1(4)=cplx1(p1(4))

      zp2(1)=cplx1(p2(1))
      zp2(2)=cplx1(p2(2))
      zp2(3)=cplx1(p2(3))
      zp2(4)=cplx1(p2(4))

      zaba=spstrng2(Um,zp1,zp2,Vm)
      return
      end
