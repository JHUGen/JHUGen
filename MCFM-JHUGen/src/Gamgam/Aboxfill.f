      subroutine Aboxfill(j1,j2,j3,j4,j5,za,zb,Abox)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      complex(dp)::A51ppppp,A51mpppp,A51mmppp,A51mpmpp,
     & Abox(2,2,2,2,2)
      integer:: h(5),j1,j2,j3,j4,j5
      
      Abox(2,2,2,2,2)=A51ppppp(j1,j2,j3,j4,j5,za,zb)
      Abox(1,1,1,1,1)=A51ppppp(j1,j2,j3,j4,j5,zb,za)

      call helfill(1,2,2,2,2,j1,j2,j3,j4,j5,h)
      Abox(h(1),h(2),h(3),h(4),h(5))=A51mpppp(j1,j2,j3,j4,j5,za,zb)      
      call helfill(2,1,2,2,2,j1,j2,j3,j4,j5,h)
      Abox(h(1),h(2),h(3),h(4),h(5))=A51mpppp(j2,j3,j4,j5,j1,za,zb)      
      call helfill(2,2,1,2,2,j1,j2,j3,j4,j5,h)
      Abox(h(1),h(2),h(3),h(4),h(5))=A51mpppp(j3,j4,j5,j1,j2,za,zb)
      call helfill(2,2,2,1,2,j1,j2,j3,j4,j5,h)
      Abox(h(1),h(2),h(3),h(4),h(5))=A51mpppp(j4,j5,j1,j2,j3,za,zb)
      call helfill(2,2,2,2,1,j1,j2,j3,j4,j5,h)
      Abox(h(1),h(2),h(3),h(4),h(5))=A51mpppp(j5,j1,j2,j3,j4,za,zb)

      call helfill(2,1,1,1,1,j1,j2,j3,j4,j5,h)
      Abox(h(1),h(2),h(3),h(4),h(5))=A51mpppp(j1,j2,j3,j4,j5,zb,za)
      call helfill(1,2,1,1,1,j1,j2,j3,j4,j5,h)
      Abox(h(1),h(2),h(3),h(4),h(5))=A51mpppp(j2,j3,j4,j5,j1,zb,za)
      call helfill(1,1,2,1,1,j1,j2,j3,j4,j5,h)
      Abox(h(1),h(2),h(3),h(4),h(5))=A51mpppp(j3,j4,j5,j1,j2,zb,za)
      call helfill(1,1,1,2,1,j1,j2,j3,j4,j5,h)
      Abox(h(1),h(2),h(3),h(4),h(5))=A51mpppp(j4,j5,j1,j2,j3,zb,za)
      call helfill(1,1,1,1,2,j1,j2,j3,j4,j5,h)
      Abox(h(1),h(2),h(3),h(4),h(5))=A51mpppp(j5,j1,j2,j3,j4,zb,za)

      call helfill(1,1,2,2,2,j1,j2,j3,j4,j5,h)
      Abox(h(1),h(2),h(3),h(4),h(5))=A51mmppp(j1,j2,j3,j4,j5,za,zb)
      call helfill(2,1,1,2,2,j1,j2,j3,j4,j5,h)
      Abox(h(1),h(2),h(3),h(4),h(5))=A51mmppp(j2,j3,j4,j5,j1,za,zb)
      call helfill(2,2,1,1,2,j1,j2,j3,j4,j5,h)
      Abox(h(1),h(2),h(3),h(4),h(5))=A51mmppp(j3,j4,j5,j1,j2,za,zb)
      call helfill(2,2,2,1,1,j1,j2,j3,j4,j5,h)
      Abox(h(1),h(2),h(3),h(4),h(5))=A51mmppp(j4,j5,j1,j2,j3,za,zb)
      call helfill(1,2,2,2,1,j1,j2,j3,j4,j5,h)
      Abox(h(1),h(2),h(3),h(4),h(5))=A51mmppp(j5,j1,j2,j3,j4,za,zb)

      call helfill(2,2,1,1,1,j1,j2,j3,j4,j5,h)
      Abox(h(1),h(2),h(3),h(4),h(5))=A51mmppp(j1,j2,j3,j4,j5,zb,za)
      call helfill(1,2,2,1,1,j1,j2,j3,j4,j5,h)
      Abox(h(1),h(2),h(3),h(4),h(5))=A51mmppp(j2,j3,j4,j5,j1,zb,za)
      call helfill(1,1,2,2,1,j1,j2,j3,j4,j5,h)
      Abox(h(1),h(2),h(3),h(4),h(5))=A51mmppp(j3,j4,j5,j1,j2,zb,za)
      call helfill(1,1,1,2,2,j1,j2,j3,j4,j5,h)
      Abox(h(1),h(2),h(3),h(4),h(5))=A51mmppp(j4,j5,j1,j2,j3,zb,za)
      call helfill(2,1,1,1,2,j1,j2,j3,j4,j5,h)
      Abox(h(1),h(2),h(3),h(4),h(5))=A51mmppp(j5,j1,j2,j3,j4,zb,za)

      call helfill(1,2,1,2,2,j1,j2,j3,j4,j5,h)
      Abox(h(1),h(2),h(3),h(4),h(5))=A51mpmpp(j1,j2,j3,j4,j5,za,zb)
      call helfill(2,1,2,2,1,j1,j2,j3,j4,j5,h)
      Abox(h(1),h(2),h(3),h(4),h(5))=A51mpmpp(j5,j1,j2,j3,j4,za,zb)
      call helfill(1,2,2,1,2,j1,j2,j3,j4,j5,h)
      Abox(h(1),h(2),h(3),h(4),h(5))=A51mpmpp(j4,j5,j1,j2,j3,za,zb)
      call helfill(2,2,1,2,1,j1,j2,j3,j4,j5,h)
      Abox(h(1),h(2),h(3),h(4),h(5))=A51mpmpp(j3,j4,j5,j1,j2,za,zb)
      call helfill(2,1,2,1,2,j1,j2,j3,j4,j5,h)
      Abox(h(1),h(2),h(3),h(4),h(5))=A51mpmpp(j2,j3,j4,j5,j1,za,zb)

      call helfill(2,1,2,1,1,j1,j2,j3,j4,j5,h)
      Abox(h(1),h(2),h(3),h(4),h(5))=A51mpmpp(j1,j2,j3,j4,j5,zb,za)
      call helfill(1,2,1,1,2,j1,j2,j3,j4,j5,h)
      Abox(h(1),h(2),h(3),h(4),h(5))=A51mpmpp(j5,j1,j2,j3,j4,zb,za)
      call helfill(2,1,1,2,1,j1,j2,j3,j4,j5,h)
      Abox(h(1),h(2),h(3),h(4),h(5))=A51mpmpp(j4,j5,j1,j2,j3,zb,za)
      call helfill(1,1,2,1,2,j1,j2,j3,j4,j5,h)
      Abox(h(1),h(2),h(3),h(4),h(5))=A51mpmpp(j3,j4,j5,j1,j2,zb,za)
      call helfill(1,2,1,2,1,j1,j2,j3,j4,j5,h)
      Abox(h(1),h(2),h(3),h(4),h(5))=A51mpmpp(j2,j3,j4,j5,j1,zb,za)

      return
      end



      subroutine helfill(h1,h2,h3,h4,h5,j1,j2,j3,j4,j5,harray)
      implicit none
      include 'types.f'
      
      integer:: h1,h2,h3,h4,h5,j1,j2,j3,j4,j5,harray(5)
      
      harray(j1)=h1
      harray(j2)=h2
      harray(j3)=h3
      harray(j4)=h4
      harray(j5)=h5
      
      return
      end
      
