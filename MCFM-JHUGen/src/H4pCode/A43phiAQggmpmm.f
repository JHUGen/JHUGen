      function A43phiAQggmpmm(j1,j2,j3,j4,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: A43phiAQggmpmm

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      complex(dp):: A43phiAQggmpmm_unsym
      integer:: j1,j2,j3,j4
      A43phiAQggmpmm=
     & +A43phiAQggmpmm_unsym(j1,j2,j3,j4,za,zb)
     & +A43phiAQggmpmm_unsym(j1,j2,j4,j3,za,zb)
      return
      end

      function A43phiAQggmpmm_unsym(j1,j2,j3,j4,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: A43phiAQggmpmm_unsym
c--- This is an implementation of Eq. (5.17) in
c---  S.~Badger, John.~M.~Campbell, R.~Keith Ellis and Ciaran Williams
c---  "Analytic results for the one-loop NMHV H-qbar-q-g-g amplitude."
c---   preprint DESY 09-180, FERMILAB-PUB-09-505-T, IPPP/09/86
c---   arXiv: 0910.4481 [hep-ph]

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'scale.f'
      include 'epinv.f'
      complex(dp):: zab2,V5L,
     & A0phiAQggmpmm,lnrat,sum
      complex(dp):: l34,l12,l24,l13,Lsm1,Lsm1_2mht
      real(dp):: s3,mhsq,s123,s234,s341,s412
      integer:: j1,j2,j3,j4
      zab2(j1,j2,j3,j4)=+za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
      s3(j1,j2,j3)=s(j1,j2)+s(j1,j3)+s(j2,j3)
      s123=s3(j1,j2,j3)
      s234=s3(j2,j3,j4)
      s341=s3(j3,j4,j1)
      s412=s3(j4,j1,j2)
      mhsq=s(j1,j2)+s(j1,j3)+s(j1,j4)+s(j2,j3)+s(j2,j4)+s(j3,j4)
      l12=lnrat(musq,-s(j1,j2))
      l24=lnrat(musq,-s(j2,j4))
      l34=lnrat(musq,-s(j3,j4))
      l13=lnrat(musq,-s(j1,j3))

c--- new representation of poles
      V5L=
     & -epinv**2-epinv*l12-0.5_dp*l12**2
     & -epinv**2-epinv*l34-0.5_dp*l34**2
     & +epinv**2+epinv*l13+0.5_dp*l13**2
     & +epinv**2+epinv*l24+0.5_dp*l24**2

      sum=+A0phiAQggmpmm(j1,j2,j3,j4,za,zb)*V5L

      sum=sum
     & +((zab2(j4,j1,j2,j3)**2*zb(j1,j2)**2/zb(j2,j3)
     & +zb(j2,j3)**2*zab2(j4,j2,j3,j1)**3
     &  /(zb(j1,j2)*zab2(j4,j1,j2,j3)))/zb(j1,j3)**3
     & -mhsq**2*za(j1,j3)**3
     &  /(za(j1,j2)*zab2(j1,j2,j3,j4)*zab2(j3,j1,j2,j4)))
     & *Lsm1(-s(j1,j2),-s123,-s(j2,j3),-s123)/s123

     & +(mhsq**2*za(j1,j4)**2*za(j2,j4)
     &  /(za(j1,j2)*zab2(j2,j1,j4,j3)*zab2(j4,j1,j2,j3))
     & -zab2(j3,j1,j4,j2)**2*zab2(j3,j2,j4,j1)
     &  /(zb(j1,j4)*zab2(j3,j1,j2,j4)*zb(j1,j2)))
     & *Lsm1(-s(j1,j2),-s412,-s(j1,j4),-s412)/s412

     & +(mhsq**2*za(j1,j3)**2/(zab2(j1,j2,j3,j4)*zab2(j2,j1,j3,j4))
     & -zab2(j4,j1,j3,j2)**2/(zb(j1,j3)*zb(j2,j3)))
     & *Lsm1(-s(j1,j3),-s123,-s(j2,j3),-s123)/s123

     & +s341**2
     & /(zb(j1,j3)*zb(j3,j4)*zab2(j2,j1,j3,j4))
     & *(Lsm1(-s(j1,j3),-s341,-s(j1,j4),-s341)
     &  +Lsm1_2mht(s(j2,j3),s341,s(j1,j4),mhsq))

     & +s341**2
     & /(zb(j1,j4)*zb(j3,j4)*zab2(j2,j1,j4,j3))
     & *(Lsm1(-s(j1,j4),-s341,-s(j3,j4),-s341)
     &  +Lsm1_2mht(s(j1,j2),s341,s(j3,j4),mhsq))

     & -zab2(j1,j3,j4,j2)**2
     & /(zb(j2,j3)*zb(j3,j4)*zab2(j1,j2,j3,j4))
     & *(Lsm1(-s(j2,j3),-s234,-s(j3,j4),-s234)
     &  +Lsm1_2mht(s(j1,j2),s234,s(j3,j4),mhsq))

     & +zb(j2,j4)**2*zab2(j1,j2,j4,j3)**2
     & /(zb(j2,j3)*zb(j3,j4)**3*zab2(j1,j2,j3,j4))
     & *Lsm1(-s(j2,j3),-s234,-s(j2,j4),-s234)

     & -zab2(j1,j3,j4,j2)**2
     & /(zab2(j1,j2,j4,j3)*zb(j2,j4)*zb(j3,j4))
     & *Lsm1_2mht(s(j1,j4),s234,s(j2,j3),mhsq)

     & +(-mhsq**2*za(j1,j3)**2*za(j2,j3)
     & /(za(j1,j2)*zab2(j2,j1,j3,j4)*zab2(j3,j1,j2,j4))
     & +zab2(j4,j1,j3,j2)**2/zab2(j4,j1,j2,j3)
     & *zab2(j4,j2,j3,j1)/(zb(j1,j2)*zb(j1,j3)))
     & *Lsm1_2mht(s(j1,j4),s123,s(j2,j3),mhsq)/s123

     & +(mhsq**2*za(j1,j3)**3
     & /(za(j1,j2)*zab2(j3,j1,j2,j4)*zab2(j1,j2,j3,j4))
     & -zab2(j4,j1,j3,j2)**3/(zab2(j4,j1,j2,j3)*zb(j1,j2)*zb(j2,j3)))
     & *Lsm1_2mht(s(j2,j4),s123,s(j1,j3),mhsq)/s123

     & +(-mhsq**2*za(j1,j3)**2
     & /(zab2(j2,j1,j3,j4)*zab2(j1,j2,j3,j4))
     & +zab2(j4,j1,j3,j2)**2/(zb(j1,j3)*zb(j2,j3)))
     & *Lsm1_2mht(s(j3,j4),s123,s(j1,j2),mhsq)/s123

      A43phiAQggmpmm_unsym=sum

      return
      end
