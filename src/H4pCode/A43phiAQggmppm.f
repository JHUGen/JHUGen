      double complex function A43phiAQggmppm(j1,j2,j3,j4,za,zb)
      implicit none
      include 'constants.f'
      include 'zprods_decl.f'
      double complex A43phiAQggmpmp_unsym,A43phiAQggmppm_unsym
      integer j1,j2,j3,j4
      A43phiAQggmppm=
     & +A43phiAQggmppm_unsym(j1,j2,j3,j4,za,zb)
     & +A43phiAQggmpmp_unsym(j1,j2,j4,j3,za,zb)
      return
      end

c--- The expression below corresponds to one half of A4;3:
c---   Aleft(a1,q2,g3+,g4-)+Aright(a1,q2,g3+,g4-)+Aleft(a1,g3+,q2,g4-)

c--- It is obtained by summing the expressions given in:
c---  L.~J.~Dixon and Y.~Sofianatos,
c---  %``Analytic one-loop amplitudes for a Higgs boson plus four partons,''
c---  arXiv:0906.0008 [hep-ph].
      double complex function A43phiAQggmppm_unsym(j1,j2,j3,j4,za,zb)
      implicit none
      include 'constants.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'scale.f'
      include 'epinv.f'
      double complex V5L,V6l,A0phiAQggmppm,lnrat
      double complex l13,l34,l12,l24,Lsm1DS,Lsm1_2me
      double precision s3,mhsq,s123,s234,s341,s412
      integer j1,j2,j3,j4
      s3(j1,j2,j3)=s(j1,j2)+s(j1,j3)+s(j2,j3)

      s123=s3(j1,j2,j3)
      s234=s3(j2,j3,j4)
      s341=s3(j3,j4,j1)
      s412=s3(j4,j1,j2)
      mhsq=s(j1,j2)+s(j1,j3)+s(j1,j4)+s(j2,j3)+s(j2,j4)+s(j3,j4)      
      l12=lnrat(musq,-s(j1,j2))
      l13=lnrat(musq,-s(j1,j3))
      l24=lnrat(musq,-s(j2,j4))
      l34=lnrat(musq,-s(j3,j4))

c--- This is the same function as in the mpmm amplitude
      V5L=
     & -epinv**2-epinv*l12-0.5d0*l12**2
     & -epinv**2-epinv*l34-0.5d0*l34**2
     & +epinv**2+epinv*l13+0.5d0*l13**2
     & +epinv**2+epinv*l24+0.5d0*l24**2

c--- These are additional boxes
      V6L=
     . +Lsm1_2me(s341,s123,s(j1,j3),mhsq) 
     . +Lsm1_2me(s412,s234,s(j2,j4),mhsq) 
     . -Lsm1_2me(s412,s123,s(j1,j2),mhsq) 
     . -Lsm1_2me(s234,s341,s(j3,j4),mhsq) 
     . +Lsm1DS(s(j2,j4),s(j4,j1),s412)
     . +Lsm1DS(s(j1,j3),s(j3,j2),s123)
     . +Lsm1DS(s(j2,j4),s(j2,j3),s234)
     . +Lsm1DS(s(j1,j3),s(j1,j4),s341)
      
      A43phiAQggmppm_unsym=
     .  A0phiAQggmppm(j1,j2,j3,j4,za,zb)
     .  *(V5L+V6L
     .    -Lsm1DS(s(j2,j3),s(j3,j4),s234)
     .    -Lsm1DS(s(j4,j1),s(j1,j2),s3(j4,j1,j2)))

     . +(za(j1,j4)**3/(za(j1,j2)*za(j3,j4)*za(j1,j3))
     .  +za(j1,j4)**2/(za(j2,j3)*za(j1,j3)))
     .  *(Lsm1DS(s(j1,j2),s(j2,j3),s123)
     .   +Lsm1DS(s(j3,j4),s(j4,j1),s3(j3,j4,j1)))
     
      return
      end
