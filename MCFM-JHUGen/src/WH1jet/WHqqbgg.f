      function WHqqbgg(p1,p2,p7,p8)
C-----Calculate the process
!      0 --> qb(-p1)+q(p2)+W(p3,p4)H(pX)+g(p7)+g(p8)
!      without using momentum conservation
!     returns matrix element squared
!     we remove the factor of
!     gwsq^3*Wmass^2/((s1278-MW^2)^2+MW^2*GamW^2)/((s34-MW^2)^2+MW^2*GamW^2)
!     color gives

      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'zprods_com.f'
      include 'sprods_com.f'
      integer::p1,p2,p7,p8,j1,j2,j3,j4,j5,h1,h2
      integer,parameter::p3=3,p4=4
      real(dp)::WHqqbgg,s3
      complex(dp)::zba3,zba2,amp78(2,2),amp87(2,2),amqed(2,2)
      s3(j1,j2,j3)=s(j1,j2)+s(j2,j3)+s(j3,j1)
      zba3(j1,j2,j3,j4,j5)=
     & +zb(j1,j2)*za(j2,j5)+zb(j1,j3)*za(j3,j5)+zb(j1,j4)*za(j4,j5)
      zba2(j1,j2,j3,j4)=zb(j1,j2)*za(j2,j4)+zb(j1,j3)*za(j3,j4)

      amp78(1,1)=-zb(p1,p4)*zba3(p1,p2,p7,p8,p3)
     & /(zb(p1,p8)*zb(p2,p7)*zb(p7,p8))
      amp78(2,2)=+za(p2,p3)*zba3(p4,p1,p7,p8,p2)
     & /(za(p1,p8)*za(p2,p7)*za(p7,p8))
      amp78(1,2)=(-zba2(p4,p1,p8,p7)*zba2(p8,p2,p7,p3)
     & /(za(p1,p8)*zb(p2,p7))
     & +za(p1,p7)*za(p2,p3)*zb(p1,p8)*zba3(p4,p1,p7,p8,p7)
     & /(za(p1,p8)*s3(p1,p7,p8))
     & +za(p2,p7)*zb(p1,p4)*zb(p2,p8)*zba2(p8,p2,p7,p3)
     & /(zb(p2,p7)*s3(p2,p7,p8)))/s(p7,p8)
      amp78(2,1)=(-za(p2,p3)*za(p2,p8)*zb(p1,p4)*zb(p1,p7)
     & /(za(p2,p7)*zb(p1,p8))
     & +za(p2,p3)*zb(p1,p7)**2*zba3(p4,p1,p7,p8,p8)
     & /(zb(p1,p8)*s3(p1,p7,p8))
     & +za(p2,p8)**2*zb(p1,p4)*zba2(p7,p2,p8,p3)
     & /(za(p2,p7)*s3(p2,p7,p8)))/(s(p7,p8))


      amp87(1,1)=-zb(p1,p4)*zba3(p1,p2,p8,p7,p3)
     & /(zb(p1,p7)*zb(p2,p8)*zb(p8,p7))
      amp87(2,2)=+za(p2,p3)*zba3(p4,p1,p8,p7,p2)
     & /(za(p1,p7)*za(p2,p8)*za(p8,p7))
      amp87(2,1)=(-zba2(p4,p1,p7,p8)*zba2(p7,p2,p8,p3)
     & /(za(p1,p7)*zb(p2,p8))
     & +za(p1,p8)*za(p2,p3)*zb(p1,p7)*zba3(p4,p1,p8,p7,p8)
     & /(za(p1,p7)*s3(p1,p8,p7))
     & +za(p2,p8)*zb(p1,p4)*zb(p2,p7)*zba2(p7,p2,p8,p3)
     & /(zb(p2,p8)*s3(p2,p8,p7)))/s(p8,p7)
      amp87(1,2)=(-za(p2,p3)*za(p2,p7)*zb(p1,p4)*zb(p1,p8)
     & /(za(p2,p8)*zb(p1,p7))
     & +za(p2,p3)*zb(p1,p8)**2*zba3(p4,p1,p8,p7,p7)
     & /(zb(p1,p7)*s3(p1,p8,p7))
     & +za(p2,p7)**2*zb(p1,p4)*zba2(p8,p2,p7,p3)
     & /(za(p2,p8)*s3(p2,p8,p7)))/(s(p8,p7))
      amqed(:,:)=amp78(:,:)+amp87(:,:)
      WHqqbgg=zip
      do h1=1,2
      do h2=1,2
      WHqqbgg=WHqqbgg
     &  +real(amp78(h1,h2)*Conjg(amp78(h1,h2)),dp)
     &  +real(amp87(h1,h2)*Conjg(amp87(h1,h2)),dp)
     &  -real(amqed(h1,h2)*Conjg(amqed(h1,h2)),dp)/xnsq
      enddo
      enddo
      return
      end


