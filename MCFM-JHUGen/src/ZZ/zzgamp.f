      subroutine zzgamp(k1,k2,k3,k4,k5,k6,k7,za,zb,z12,z34,z56)
      implicit none
      include 'types.f'
      
C-----Author:R.K. Ellis, Novemeber 2013
c---  This is the new code for the amplitudes for ZZ+gluon production
c---  (singly-resonant diagrams are included)
C---  Process calculated is 
C---  0 --> q(-p1)+qb(-p2)+e-(p3)+e+(p4)+mu-(p5)+mu(p6)+g(p7)
C---  
C---  Order of indices for z12 is h1,hg,h34,h56
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'srdiags.f'
      complex(dp):: z12(2,2,2,2),z34(2,2,2,2),z56(2,2,2,2)
      complex(dp):: DRP,DRM,SRmm,SRmp,SRpm,SRpp
      integer:: p1,p2,p3,p4,p5,p6,p7,h34,h56
      integer:: k1,k2,k3,k4,k5,k6,k7
C--   order of indices is z12(hq,h34,h56,hg)
      
      p1=k1
      p2=k2
      p7=k7

      do h34=1,2
         if (h34 == 1) then
            p3=k3
            p4=k4
         elseif (h34 == 2) then
            p3=k4
            p4=k3
         endif
         do h56=1,2
            if (h56 == 1) then
               p5=k5
               p6=k6
            elseif (h56 == 2) then
               p5=k6
               p6=k5
            endif
            z12(1,h34,h56,1)=DRM(p1,p2,p3,p4,p5,p6,p7,za,zb)
            z12(1,h34,h56,2)=DRP(p1,p2,p3,p4,p5,p6,p7,za,zb)
C----- ( -sign required when performing (za<-->zb))
            z12(2,3-h34,3-h56,2)=-DRM(p1,p2,p3,p4,p5,p6,p7,zb,za)
            z12(2,3-h34,3-h56,1)=-DRP(p1,p2,p3,p4,p5,p6,p7,zb,za)
         enddo
      enddo

c---  do not calculate single resonant diagrams unnecessarily
      if (srdiags .eqv. .false.) return

         do h56=1,2
            if (h56 == 1) then
               p5=k5
               p6=k6
            elseif (h56 == 2) then
               p5=k6
               p6=k5
            endif

            z56(1,1,h56,1)=SRmm(k1,k2,k3,k4,p5,p6,k7,za,zb)
            z56(1,1,h56,2)=SRmp(k1,k2,k3,k4,p5,p6,k7,za,zb)
            z56(2,1,h56,1)=SRpm(k1,k2,k3,k4,p5,p6,k7,za,zb)
            z56(2,1,h56,2)=SRpp(k1,k2,k3,k4,p5,p6,k7,za,zb)
C----- ( -sign required when performing (za<-->zb))
            z56(2,2,3-h56,2)=-SRmm(k1,k2,k3,k4,p5,p6,k7,zb,za)
            z56(2,2,3-h56,1)=-SRmp(k1,k2,k3,k4,p5,p6,k7,zb,za)
            z56(1,2,3-h56,2)=-SRpm(k1,k2,k3,k4,p5,p6,k7,zb,za)
            z56(1,2,3-h56,1)=-SRpp(k1,k2,k3,k4,p5,p6,k7,zb,za)
         enddo

         do h34=1,2
            if (h34 == 1) then
               p3=k3
               p4=k4
            elseif (h34 == 2) then
               p3=k4
               p4=k3
            endif

            z34(1,h34,1,1)=SRmm(k1,k2,k5,k6,p3,p4,k7,za,zb)
            z34(1,h34,1,2)=SRmp(k1,k2,k5,k6,p3,p4,k7,za,zb)
            z34(2,h34,1,1)=SRpm(k1,k2,k5,k6,p3,p4,k7,za,zb)
            z34(2,h34,1,2)=SRpp(k1,k2,k5,k6,p3,p4,k7,za,zb)
C----- ( -sign required when performing (za<-->zb))
            z34(2,3-h34,2,2)=-SRmm(k1,k2,k5,k6,p3,p4,k7,zb,za)
            z34(2,3-h34,2,1)=-SRmp(k1,k2,k5,k6,p3,p4,k7,zb,za)
            z34(1,3-h34,2,2)=-SRpm(k1,k2,k5,k6,p3,p4,k7,zb,za)
            z34(1,3-h34,2,1)=-SRpp(k1,k2,k5,k6,p3,p4,k7,zb,za)

         enddo
      return
      end



      function DRPa(p1,p2,p3,p4,p5,p6,p7,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: DRPa
       
C---  Unsymmetrized version of the doubly resonant piece
C---  1_q^-,2_q^+,3_e^-,4_eb^+,5_mu^-,6_mb^-,7_g^+
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer:: p1,p2,p3,p4,p5,p6,p7
      real(dp):: t
      complex(dp):: zab2
C----Begin statement functions
      t(p1,p2,p3)=s(p1,p2)+s(p2,p3)+s(p1,p3)
      zab2(p1,p2,p3,p4)=za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)
C----End statement functions

      DRPa=
     & (za(p5,p1)*zab2(p1,p2,p7,p4)*zab2(p3,p1,p5,p6)/t(p1,p5,p6)
     & +za(p2,p7)*za(p5,p1)**2*zb(p2,p4)*zb(p5,p6)*zab2(p3,p2,p4,p7)
     & /(t(p2,p3,p4)*t(p1,p5,p6))
     & )/(za(p2,p7)*za(p7,p1)*s(p3,p4)*s(p5,p6))

      return
      end

      function DRMa(p1,p2,p3,p4,p5,p6,p7,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: DRMa
       
C---  Unsymmetrized version of the doubly resonant piece
C---  1_q^-,2_qb^+,3_e^-,4_eb^+,5_mu^-,6_mb^-,7_g^-
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer:: p1,p2,p3,p4,p5,p6,p7
      real(dp):: t
      complex(dp):: zab2
C----Begin statement functions
      t(p1,p2,p3)=s(p1,p2)+s(p2,p3)+s(p1,p3)
      zab2(p1,p2,p3,p4)=za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)
C----End statement functions

      DRMa=
     & (zb(p2,p6)*zab2(p3,p1,p7,p2)*zab2(p5,p2,p6,p4)/t(p2,p5,p6)
     & +za(p3,p1)*za(p5,p6)*zb(p2,p6)**2*zb(p7,p1)*zab2(p7,p3,p1,p4)
     & /(t(p1,p3,p4)*t(p2,p5,p6)))
     & /(zb(p2,p7)*zb(p7,p1)*s(p3,p4)*s(p5,p6))
      return
      end


      function SRmp(p1,p2,p3,p4,p5,p6,p7,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: SRmp
       
C---  Unsymmetrized version of the singly resonant piece
C---  1_q^-,2_qb^+,3_e^-,4_eb^+,5_mu^-,6_mb^-,7_g^+
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer:: p1,p2,p3,p4,p5,p6,p7
      real(dp):: t,s3456
      complex(dp):: zab2
      zab2(p1,p2,p3,p4)=za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)
      s3456=t(p1,p2,p7)
      SRmp=
     & (za(p1,p3)*zb(p6,p4)/t(p4,p5,p6)
     & *(zab2(p5,p4,p6,p2)*za(p2,p1)+zab2(p5,p4,p6,p7)*za(p7,p1))
     & +za(p5,p3)*zab2(p1,p3,p5,p6)*zab2(p1,p2,p7,p4)/t(p3,p5,p6))
     & /(za(p2,p7)*za(p1,p7)*s3456*s(p5,p6))

      end

      function SRmm(p1,p2,p3,p4,p5,p6,p7,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: SRmm
       
C---  Unsymmetrized version of the singly resonant piece
C---  1_q^-,2_qb^+,3_e^-,4_eb^+,5_mu^-,6_mb^-,7_g^-
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer:: p1,p2,p3,p4,p5,p6,p7
      real(dp):: t,s3456
      complex(dp):: zab2
      zab2(p1,p2,p3,p4)=za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)
      s3456=t(p1,p2,p7)
      SRmm=(za(p3,p5)*zb(p2,p4)/t(p3,p5,p6)
     & *(zb(p2,p1)*zab2(p1,p3,p5,p6)+zb(p2,p7)*zab2(p7,p3,p5,p6))
     &  +zb(p6,p4)*zab2(p3,p1,p7,p2)*zab2(p5,p4,p6,p2)/t(p4,p5,p6))
     & /(zb(p2,p7)*zb(p7,p1)*s3456*s(p5,p6))
      return
      end

      function SRpp(p1,p2,p3,p4,p5,p6,p7,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: SRpp
       
C---  Unsymmetrized version of the singly resonant piece
C---  1_q^+,2_qb^-,3_e^-,4_eb^+,5_mu^-,6_mb^-,7_g^+
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer:: p1,p2,p3,p4,p5,p6,p7
      real(dp):: t,s3456
      complex(dp):: zab2
      zab2(p1,p2,p3,p4)=za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)
      s3456=t(p1,p2,p7)
      SRpp=
     & (-za(p2,p3)*zb(p4,p6)/t(p4,p5,p6)
     & *(zab2(p5,p4,p6,p1)*za(p1,p2)+zab2(p5,p4,p6,p7)*za(p7,p2))
     & -za(p3,p5)*zab2(p2,p1,p7,p4)*zab2(p2,p3,p5,p6)/t(p3,p5,p6))
     & /(za(p7,p1)*za(p2,p7)*s3456*s(p5,p6))
      return
      end

      function SRpm(p1,p2,p3,p4,p5,p6,p7,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: SRpm
       
C---  Unsymmetrized version of the singly resonant piece
C---  1_q^+,2_qb^-,3_e^-,4_eb^+,5_mu^-,6_mb^-,7_g^-
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer:: p1,p2,p3,p4,p5,p6,p7
      real(dp):: t,s3456
      complex(dp):: zab2
      zab2(p1,p2,p3,p4)=za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)
      s3456=t(p1,p2,p7)
      SRpm=
     & (zab2(p3,p2,p7,p1)*zb(p4,p6)*zab2(p5,p4,p6,p1)/t(p4,p5,p6) 
     & -za(p3,p5)*zb(p1,p4)/t(p3,p5,p6)
     & *(zb(p1,p2)*zab2(p2,p3,p5,p6)+zb(p1,p7)*zab2(p7,p3,p5,p6)))
     & /(zb(p7,p1)*zb(p2,p7)*s3456*s(p5,p6))
      return
      end


      function DRP(p1,p2,p3,p4,p5,p6,p7,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: DRP
       
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      integer:: p1,p2,p3,p4,p5,p6,p7
      complex(dp):: DRPa
      DRP=
     & DRPa(p1,p2,p3,p4,p5,p6,p7,za,zb)+DRPa(p1,p2,p5,p6,p3,p4,p7,za,zb)
      return
      end

      function DRM(p1,p2,p3,p4,p5,p6,p7,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: DRM
       
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      integer:: p1,p2,p3,p4,p5,p6,p7
      complex(dp):: DRMa
      DRM=
     & DRMa(p1,p2,p3,p4,p5,p6,p7,za,zb)+DRMa(p1,p2,p5,p6,p3,p4,p7,za,zb)
      return
      end



