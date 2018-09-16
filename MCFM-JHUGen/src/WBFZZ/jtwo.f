      subroutine jtwo(n2,n3,n4,n5,n6,n1,za,zb,zab,zba,j2,jw2)
      implicit none
      include 'types.f'
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      complex(dp):: zab(mxpart,4,mxpart),zba(mxpart,4,mxpart),
     & j2(4,2,2,2,2),j2_34_56_1(4,2,2,2,2),j2_56_34_1(4,2,2,2,2),
     & jw2(4,2,2,2),jw2_34_56_1(4,2,2,2),jw2_56_34_1(4,2,2,2)
      integer:: n1,n2,n3,n4,n5,n6,jdu,h21,h34,h56
C---The two Z-current divided by (-i)
      call jtwo3456(n2,n3,n4,n5,n6,n1,
     & za,zb,zab,zba,j2_34_56_1,jw2_34_56_1)
      call jtwo3456(n2,n5,n6,n3,n4,n1,
     & za,zb,zab,zba,j2_56_34_1,jw2_56_34_1)

      do jdu=1,2
      do h56=1,2
      do h34=1,2
      jw2(1:4,jdu,h34,h56)=
     &  jw2_34_56_1(1:4,jdu,h34,h56)
     & +jw2_56_34_1(1:4,jdu,h56,h34)
      do h21=1,2
      j2(1:4,jdu,h21,h34,h56)=
     &  j2_34_56_1(1:4,jdu,h21,h34,h56)
     & +j2_56_34_1(1:4,jdu,h21,h56,h34)
      enddo
      enddo
      enddo
      enddo

      return
      end 


