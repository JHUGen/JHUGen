! T. Dennen
! fill arrays of lo and nlo colour-ordered helicity amplitudes
! q(j1) qb(j2) a(j3) a(j4) a(j5) a(j6)
      subroutine aaaa_fill(j1,j2,j3,j4,j5,j6,za,zb,aaaa,aaaa_lo,ord)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      integer:: j1,j2,j3,j4,j5,j6,ord
      complex(dp):: aaaa(2,2,2,2,2), aaaa_lo(2,2,2,2,2)
      complex(dp):: Boxint(195), Triint(90), Bubint(25), Ratint
      complex(dp):: dcoeffs(195), ccoeffs(90), bcoeffs(25), rats
      complex(dp):: aaaa_MHV_lo, aaaa_NMHV_lo

! aaaa and aaaa_lo are dimension 5 arrays, aaaa(hel(j1),hel(j3),...) etc
! because helicity of quark line is conserved
      aaaa(:,:,:,:,:)=czip
      aaaa_lo(:,:,:,:,:)=czip

      if (ord==0) then
            Ratint = -cone*0.5
      else
            Ratint = czip
      endif

      call Boxint_init(j1,j2,j3,j4,j5,j6,Boxint,ord)
      call Triint_init(j1,j2,j3,j4,j5,j6,Triint,ord)
      call Bubint_init(j1,j2,j3,j4,j5,j6,Bubint,ord)

      call aaaa_MHV_init(j1,j2,j3,j4,j5,j6,za,zb,
     & dcoeffs,ccoeffs,bcoeffs,rats)
      aaaa(1,1,2,2,2) = sum(dcoeffs*Boxint)+sum(ccoeffs*Triint)
     &                 +sum(bcoeffs*Bubint)+rats*Ratint
      aaaa_lo(1,1,2,2,2) = aaaa_MHV_lo(j1,j2,j3,j4,j5,j6,za,zb)

      call aaaa_NMHV_init(j1,j2,j3,j4,j5,j6,za,zb,
     & dcoeffs,ccoeffs,bcoeffs,rats)
      aaaa(1,1,1,2,2) = sum(dcoeffs*Boxint)+sum(ccoeffs*Triint)
     &                 +sum(bcoeffs*Bubint)+rats*Ratint
      aaaa_lo(1,1,1,2,2) = aaaa_NMHV_lo(j1,j2,j3,j4,j5,j6,za,zb)

      call aaaa_MHV_init(j1,j2,j3,j4,j5,j6,zb,za,
     & dcoeffs,ccoeffs,bcoeffs,rats)
      aaaa(2,2,1,1,1) = sum(dcoeffs*Boxint)+sum(ccoeffs*Triint)
     &                 +sum(bcoeffs*Bubint)+rats*Ratint
      aaaa_lo(2,2,1,1,1) = aaaa_MHV_lo(j1,j2,j3,j4,j5,j6,zb,za)

      call aaaa_NMHV_init(j1,j2,j3,j4,j5,j6,zb,za,
     & dcoeffs,ccoeffs,bcoeffs,rats)
      aaaa(2,2,2,1,1) = sum(dcoeffs*Boxint)+sum(ccoeffs*Triint)
     &                 +sum(bcoeffs*Bubint)+rats*Ratint
      aaaa_lo(2,2,2,1,1) = aaaa_NMHV_lo(j1,j2,j3,j4,j5,j6,zb,za)


      call Boxint_init(j1,j2,j4,j3,j5,j6,Boxint,ord)
      call Triint_init(j1,j2,j4,j3,j5,j6,Triint,ord)
      call Bubint_init(j1,j2,j4,j3,j5,j6,Bubint,ord)

      call aaaa_MHV_init(j1,j2,j4,j3,j5,j6,za,zb,
     & 
     & dcoeffs,ccoeffs,bcoeffs,rats)
      aaaa(1,2,1,2,2) = sum(dcoeffs*Boxint)+sum(ccoeffs*Triint)
     &                 +sum(bcoeffs*Bubint)+rats*Ratint
      aaaa_lo(1,2,1,2,2) = aaaa_MHV_lo(j1,j2,j4,j3,j5,j6,za,zb)

      call aaaa_MHV_init(j1,j2,j4,j3,j5,j6,zb,za,
     & dcoeffs,ccoeffs,bcoeffs,rats)
      aaaa(2,1,2,1,1) = sum(dcoeffs*Boxint)+sum(ccoeffs*Triint)
     &                 +sum(bcoeffs*Bubint)+rats*Ratint
      aaaa_lo(2,1,2,1,1) = aaaa_MHV_lo(j1,j2,j4,j3,j5,j6,zb,za)


      call Boxint_init(j1,j2,j5,j4,j3,j6,Boxint,ord)
      call Triint_init(j1,j2,j5,j4,j3,j6,Triint,ord)
      call Bubint_init(j1,j2,j5,j4,j3,j6,Bubint,ord)

      call aaaa_MHV_init(j1,j2,j5,j4,j3,j6,za,zb,
     & dcoeffs,ccoeffs,bcoeffs,rats)
      aaaa(1,2,2,1,2) = sum(dcoeffs*Boxint)+sum(ccoeffs*Triint)
     &                 +sum(bcoeffs*Bubint)+rats*Ratint
      aaaa_lo(1,2,2,1,2) = aaaa_MHV_lo(j1,j2,j5,j4,j3,j6,za,zb)

      call aaaa_NMHV_init(j1,j2,j5,j4,j3,j6,za,zb,
     & dcoeffs,ccoeffs,bcoeffs,rats)
      aaaa(1,2,1,1,2) = sum(dcoeffs*Boxint)+sum(ccoeffs*Triint)
     &                 +sum(bcoeffs*Bubint)+rats*Ratint
      aaaa_lo(1,2,1,1,2) = aaaa_NMHV_lo(j1,j2,j5,j4,j3,j6,za,zb)

      call aaaa_MHV_init(j1,j2,j5,j4,j3,j6,zb,za,
     & dcoeffs,ccoeffs,bcoeffs,rats)
      aaaa(2,1,1,2,1) = sum(dcoeffs*Boxint)+sum(ccoeffs*Triint)
     &                 +sum(bcoeffs*Bubint)+rats*Ratint
      aaaa_lo(2,1,1,2,1) = aaaa_MHV_lo(j1,j2,j5,j4,j3,j6,zb,za)

      call aaaa_NMHV_init(j1,j2,j5,j4,j3,j6,zb,za,
     & dcoeffs,ccoeffs,bcoeffs,rats)
      aaaa(2,1,2,2,1) = sum(dcoeffs*Boxint)+sum(ccoeffs*Triint)
     &                 +sum(bcoeffs*Bubint)+rats*Ratint
      aaaa_lo(2,1,2,2,1) = aaaa_NMHV_lo(j1,j2,j5,j4,j3,j6,zb,za)


      call Boxint_init(j1,j2,j6,j4,j5,j3,Boxint,ord)
      call Triint_init(j1,j2,j6,j4,j5,j3,Triint,ord)
      call Bubint_init(j1,j2,j6,j4,j5,j3,Bubint,ord)

      call aaaa_MHV_init(j1,j2,j6,j4,j5,j3,za,zb,
     & dcoeffs,ccoeffs,bcoeffs,rats)
      aaaa(1,2,2,2,1) = sum(dcoeffs*Boxint)+sum(ccoeffs*Triint)
     &                 +sum(bcoeffs*Bubint)+rats*Ratint
      aaaa_lo(1,2,2,2,1) = aaaa_MHV_lo(j1,j2,j6,j4,j5,j3,za,zb)

      call aaaa_NMHV_init(j1,j2,j6,j4,j5,j3,za,zb,
     & dcoeffs,ccoeffs,bcoeffs,rats)
      aaaa(1,2,1,2,1) = sum(dcoeffs*Boxint)+sum(ccoeffs*Triint)
     &                 +sum(bcoeffs*Bubint)+rats*Ratint
      aaaa_lo(1,2,1,2,1) = aaaa_NMHV_lo(j1,j2,j6,j4,j5,j3,za,zb)

      call aaaa_MHV_init(j1,j2,j6,j4,j5,j3,zb,za,
     & dcoeffs,ccoeffs,bcoeffs,rats)
      aaaa(2,1,1,1,2) = sum(dcoeffs*Boxint)+sum(ccoeffs*Triint)
     &                 +sum(bcoeffs*Bubint)+rats*Ratint
      aaaa_lo(2,1,1,1,2) = aaaa_MHV_lo(j1,j2,j6,j4,j5,j3,zb,za)

      call aaaa_NMHV_init(j1,j2,j6,j4,j5,j3,zb,za,
     & dcoeffs,ccoeffs,bcoeffs,rats)
      aaaa(2,1,2,1,2) = sum(dcoeffs*Boxint)+sum(ccoeffs*Triint)
     &                 +sum(bcoeffs*Bubint)+rats*Ratint
      aaaa_lo(2,1,2,1,2) = aaaa_NMHV_lo(j1,j2,j6,j4,j5,j3,zb,za)



!========================

      call Boxint_init(j2,j1,j3,j4,j5,j6,Boxint,ord)
      call Triint_init(j2,j1,j3,j4,j5,j6,Triint,ord)
      call Bubint_init(j2,j1,j3,j4,j5,j6,Bubint,ord)

      call aaaa_MHV_init(j2,j1,j3,j4,j5,j6,za,zb,
     & dcoeffs,ccoeffs,bcoeffs,rats)
      aaaa(2,1,2,2,2) = sum(dcoeffs*Boxint)+sum(ccoeffs*Triint)
     &                 +sum(bcoeffs*Bubint)+rats*Ratint
      aaaa_lo(2,1,2,2,2) = aaaa_MHV_lo(j2,j1,j3,j4,j5,j6,za,zb)

      call aaaa_NMHV_init(j2,j1,j3,j4,j5,j6,za,zb,
     & dcoeffs,ccoeffs,bcoeffs,rats)
      aaaa(2,1,1,2,2) = sum(dcoeffs*Boxint)+sum(ccoeffs*Triint)
     &                 +sum(bcoeffs*Bubint)+rats*Ratint
      aaaa_lo(2,1,1,2,2) = aaaa_NMHV_lo(j2,j1,j3,j4,j5,j6,za,zb)

      call aaaa_MHV_init(j2,j1,j3,j4,j5,j6,zb,za,
     & dcoeffs,ccoeffs,bcoeffs,rats)
      aaaa(1,2,1,1,1) = sum(dcoeffs*Boxint)+sum(ccoeffs*Triint)
     &                 +sum(bcoeffs*Bubint)+rats*Ratint
      aaaa_lo(1,2,1,1,1) = aaaa_MHV_lo(j2,j1,j3,j4,j5,j6,zb,za)

      call aaaa_NMHV_init(j2,j1,j3,j4,j5,j6,zb,za,
     & dcoeffs,ccoeffs,bcoeffs,rats)
      aaaa(1,2,2,1,1) = sum(dcoeffs*Boxint)+sum(ccoeffs*Triint)
     &                 +sum(bcoeffs*Bubint)+rats*Ratint
      aaaa_lo(1,2,2,1,1) = aaaa_NMHV_lo(j2,j1,j3,j4,j5,j6,zb,za)


      call Boxint_init(j2,j1,j4,j3,j5,j6,Boxint,ord)
      call Triint_init(j2,j1,j4,j3,j5,j6,Triint,ord)
      call Bubint_init(j2,j1,j4,j3,j5,j6,Bubint,ord)

      call aaaa_MHV_init(j2,j1,j4,j3,j5,j6,za,zb,
     & dcoeffs,ccoeffs,bcoeffs,rats)
      aaaa(2,2,1,2,2) = sum(dcoeffs*Boxint)+sum(ccoeffs*Triint)
     &                 +sum(bcoeffs*Bubint)+rats*Ratint
      aaaa_lo(2,2,1,2,2) = aaaa_MHV_lo(j2,j1,j4,j3,j5,j6,za,zb)

      call aaaa_MHV_init(j2,j1,j4,j3,j5,j6,zb,za,
     & dcoeffs,ccoeffs,bcoeffs,rats)
      aaaa(1,1,2,1,1) = sum(dcoeffs*Boxint)+sum(ccoeffs*Triint)
     &                 +sum(bcoeffs*Bubint)+rats*Ratint
      aaaa_lo(1,1,2,1,1) = aaaa_MHV_lo(j2,j1,j4,j3,j5,j6,zb,za)


      call Boxint_init(j2,j1,j5,j4,j3,j6,Boxint,ord)
      call Triint_init(j2,j1,j5,j4,j3,j6,Triint,ord)
      call Bubint_init(j2,j1,j5,j4,j3,j6,Bubint,ord)

      call aaaa_MHV_init(j2,j1,j5,j4,j3,j6,za,zb,
     & dcoeffs,ccoeffs,bcoeffs,rats)
      aaaa(2,2,2,1,2) = sum(dcoeffs*Boxint)+sum(ccoeffs*Triint)
     &                 +sum(bcoeffs*Bubint)+rats*Ratint
      aaaa_lo(2,2,2,1,2) = aaaa_MHV_lo(j2,j1,j5,j4,j3,j6,za,zb)

      call aaaa_NMHV_init(j2,j1,j5,j4,j3,j6,za,zb,
     & dcoeffs,ccoeffs,bcoeffs,rats)
      aaaa(2,2,1,1,2) = sum(dcoeffs*Boxint)+sum(ccoeffs*Triint)
     &                 +sum(bcoeffs*Bubint)+rats*Ratint
      aaaa_lo(2,2,1,1,2) = aaaa_NMHV_lo(j2,j1,j5,j4,j3,j6,za,zb)

      call aaaa_MHV_init(j2,j1,j5,j4,j3,j6,zb,za,
     & dcoeffs,ccoeffs,bcoeffs,rats)
      aaaa(1,1,1,2,1) = sum(dcoeffs*Boxint)+sum(ccoeffs*Triint)
     &                 +sum(bcoeffs*Bubint)+rats*Ratint
      aaaa_lo(1,1,1,2,1) = aaaa_MHV_lo(j2,j1,j5,j4,j3,j6,zb,za)

      call aaaa_NMHV_init(j2,j1,j5,j4,j3,j6,zb,za,
     & dcoeffs,ccoeffs,bcoeffs,rats)
      aaaa(1,1,2,2,1) = sum(dcoeffs*Boxint)+sum(ccoeffs*Triint)
     &                 +sum(bcoeffs*Bubint)+rats*Ratint
      aaaa_lo(1,1,2,2,1) = aaaa_NMHV_lo(j2,j1,j5,j4,j3,j6,zb,za)


      call Boxint_init(j2,j1,j6,j4,j5,j3,Boxint,ord)
      call Triint_init(j2,j1,j6,j4,j5,j3,Triint,ord)
      call Bubint_init(j2,j1,j6,j4,j5,j3,Bubint,ord)

      call aaaa_MHV_init(j2,j1,j6,j4,j5,j3,za,zb,
     & dcoeffs,ccoeffs,bcoeffs,rats)
      aaaa(2,2,2,2,1) = sum(dcoeffs*Boxint)+sum(ccoeffs*Triint)
     &                 +sum(bcoeffs*Bubint)+rats*Ratint
      aaaa_lo(2,2,2,2,1) = aaaa_MHV_lo(j2,j1,j6,j4,j5,j3,za,zb)

      call aaaa_NMHV_init(j2,j1,j6,j4,j5,j3,za,zb,
     & dcoeffs,ccoeffs,bcoeffs,rats)
      aaaa(2,2,1,2,1) = sum(dcoeffs*Boxint)+sum(ccoeffs*Triint)
     &                 +sum(bcoeffs*Bubint)+rats*Ratint
      aaaa_lo(2,2,1,2,1) = aaaa_NMHV_lo(j2,j1,j6,j4,j5,j3,za,zb)

      call aaaa_MHV_init(j2,j1,j6,j4,j5,j3,zb,za,
     & dcoeffs,ccoeffs,bcoeffs,rats)
      aaaa(1,1,1,1,2) = sum(dcoeffs*Boxint)+sum(ccoeffs*Triint)
     &                 +sum(bcoeffs*Bubint)+rats*Ratint
      aaaa_lo(1,1,1,1,2) = aaaa_MHV_lo(j2,j1,j6,j4,j5,j3,zb,za)

      call aaaa_NMHV_init(j2,j1,j6,j4,j5,j3,zb,za,
     & dcoeffs,ccoeffs,bcoeffs,rats)
      aaaa(1,1,2,1,2) = sum(dcoeffs*Boxint)+sum(ccoeffs*Triint)
     &                 +sum(bcoeffs*Bubint)+rats*Ratint
      aaaa_lo(1,1,2,1,2) = aaaa_NMHV_lo(j2,j1,j6,j4,j5,j3,zb,za)



      return
      end


      subroutine aaaa_MHV_init(j1,j2,j3,j4,j5,j6,za,zb,
     & dcoeffs,ccoeffs,bcoeffs,rats)
      implicit none
      include 'types.f'
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      integer:: j1,j2,j3,j4,j5,j6
      complex(dp):: dcoeffs(195), ccoeffs(90), bcoeffs(25), rats
      call aaaa_MHV_d_init(j1,j2,j3,j4,j5,j6,za,zb,dcoeffs)
      call aaaa_MHV_c_init(j1,j2,j3,j4,j5,j6,za,zb,ccoeffs)
      call aaaa_MHV_b_init(j1,j2,j3,j4,j5,j6,za,zb,bcoeffs)
      call aaaa_MHV_r_init(j1,j2,j3,j4,j5,j6,za,zb,rats)
      return
      end


      subroutine aaaa_NMHV_init(j1,j2,j3,j4,j5,j6,za,zb,
     & dcoeffs,ccoeffs,bcoeffs,rats)
      implicit none
      include 'types.f'
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      integer:: j1,j2,j3,j4,j5,j6
      complex(dp):: dcoeffs(195), ccoeffs(90), bcoeffs(25), rats
      call aaaa_NMHV_d_init(j1,j2,j3,j4,j5,j6,za,zb,dcoeffs)
      call aaaa_NMHV_c_init(j1,j2,j3,j4,j5,j6,za,zb,ccoeffs)
      call aaaa_NMHV_b_init(j1,j2,j3,j4,j5,j6,za,zb,bcoeffs)
      call aaaa_NMHV_r_init(j1,j2,j3,j4,j5,j6,za,zb,rats)
      return
      end




      function aaaa_MHV_lo(i1,i2,i3,i4,i5,i6,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: aaaa_MHV_lo
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      integer:: i1,i2,i3,i4,i5,i6
      aaaa_MHV_lo =
     &  (2*im*za(i1,i2)**2*za(i1,i3)**2*za(i2,i3))/
     & (za(i1,i4)*za(i1,i5)*za(i1,i6)*za(i3,i2)*za(i4,i2)*za(i5,i2)*
     &   za(i6,i2))
      return
      end


      function aaaa_NMHV_lo(i1,i2,i3,i4,i5,i6,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: aaaa_NMHV_lo
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer:: i1,i2,i3,i4,i5,i6
      complex(dp):: zab2
      real(dp):: t
      zab2(i1,i2,i3,i4) = za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)
      t(i1,i2,i3) = s(i1,i2)+s(i2,i3)+s(i3,i1)
      aaaa_NMHV_lo =
     &  -2*im*((za(i1,i4)*zab2(i3,i2,i5,i6)*zb(i2,i5))/
     &    (t(i1,i4,i6)*za(i1,i6)*za(i2,i5)*zb(i1,i4)*zb(i2,i3)) + 
     &   (za(i1,i3)*zab2(i4,i2,i5,i6)*zb(i2,i5))/
     &    (t(i1,i3,i6)*za(i1,i6)*za(i2,i5)*zb(i1,i3)*zb(i2,i4)) + 
     &   (za(i1,i4)*zab2(i3,i2,i6,i5)*zb(i2,i6))/
     &    (t(i1,i4,i5)*za(i1,i5)*za(i2,i6)*zb(i1,i4)*zb(i2,i3)) + 
     &   (za(i1,i3)*zab2(i4,i2,i6,i5)*zb(i2,i6))/
     &    (t(i1,i3,i5)*za(i1,i5)*za(i2,i6)*zb(i1,i3)*zb(i2,i4)) + 
     &   (s(i1,i2)*t(i1,i3,i4)*zab2(i1,i5,i6,i2))/
     &    (za(i2,i5)*za(i2,i6)*za(i5,i1)*za(i6,i1)*zb(i2,i3)*zb(i2,i4)*
     &      zb(i3,i1)*zb(i4,i1)))
      return
      end

      
