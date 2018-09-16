!---- CW Feb 2015 
!---- amplitude for ++ gluons => Z(i3,i4) + H 


      function ggHZ_pp_box(i1,i2,i3,i4,za,zb,mt2)
!---- this is helicity amplitude for 
!===== g(i1)^+ + g(i2)^+ + (Z=>(l1(i3)^- + l2(i3)^+)) + H 
!===== and the left-handed proejection of the Z is chosen

      implicit none
      include 'types.f' 
      complex(kind=dp)::  ggHZ_pp_box
      include 'constants.f' 
      include 'mxpart.f'
    
      include 'zprods_decl.f' 
      include 'sprods_com.f' 
      include 'scale.f'
      integer:: i1,i2,i3,i4 
      real(kind=dp):: mt2
      real(kind=dp):: mZsq,mHsq,s15,s25,s12,t
      integer:: Nbox,Ntri,i
      parameter(Nbox=3,Ntri=6)

      complex(kind=dp):: qlI4,qlI3,qlI2
      integer:: d25_12,d15_12,d15_25
      parameter(d25_12=1,d15_12=2,d15_25=3)
      integer:: c25_Z,cH_25,c12,c15_H,cZ_15
      parameter(c25_Z=1,cH_25=2,c12=3,c15_H=4,cZ_15=5)
      complex(kind=dp):: D0(Nbox),C0(Ntri)
      complex(kind=dp):: di(Nbox),ci(Ntri)
      common/ggZH_basisint/D0,C0
!$omp threadprivate(/ggZH_basisint/)

     

      t(i1,i2,i3)=s(i1,i2)+s(i1,i3)+s(i2,i3)

      ggHZ_pp_box=czip
      
      
!======= kinematic configurations
      s12=s(i1,i2)
      mZsq=s(i3,i4) 
      s25=t(i1,i3,i4)
      s15=t(i2,i3,i4) 
      mHsq=s(i1,i2)+s(i1,i3)+s(i1,i4)+s(i2,i3)+s(i2,i4)+s(i3,i4)
      
      

      di(:)=czip
      ci(:)=czip

!====== integrals are now filled elsewhere and stored in common block
!----- fill integrals from QCD loop : boxes
!
!      D0(d25_12)=qlI4(mZsq,zip,zip,mHsq,s25,s12,mt2,mt2,mt2,mt2,musq,0) 
!      D0(d15_12)=qlI4(mZsq,zip,zip,mHsq,s15,s12,mt2,mt2,mt2,mt2,musq,0) 
!      D0(d15_25)=qlI4(mHsq,zip,mZsq,zip,s15,s25,mt2,mt2,mt2,mt2,musq,0) 

!----- fill integrals from QCD loop : triangles
!      C0(c25_Z)=qlI3(s25,zip,mZsq,mt2,mt2,mt2,musq,0) 
!      C0(cH_25)=qlI3(mHsq,zip,s25,mt2,mt2,mt2,musq,0) 
!      C0(c12)=qlI3(s12,zip,zip,mt2,mt2,mt2,musq,0) 
!      C0(c15_H)=qlI3(s15,zip,mHsq,mt2,mt2,mt2,musq,0) 
!      C0(cZ_15)=qlI3(mZsq,zip,s15,mt2,mt2,mt2,musq,0) 
      

!------- coefficients of integrals 

      di(d25_12)=
     &(zb(i2,i1)*(2*(za(i2,i3)*zb(i3,i1) + za(i2,i4)*zb(i4,i1))*
     -       (-2*mt2*za(i2,i3) + za(i1,i2)*za(i3,i4)*zb(i4,i1))*
     -       zb(i4,i2)*(za(i1,i3)*zb(i3,i2) + za(i1,i4)*zb(i4,i2))
     -       - za(i1,i3)*(za(i2,i3)*zb(i3,i1) + 
     -         za(i2,i4)*zb(i4,i1))*
     -       (t(i1,i3,i4)*za(i1,i2)*zb(i2,i1)*zb(i4,i2) + 
     -         2*(2*mt2 - za(i1,i2)*zb(i2,i1))*zb(i4,i1)*
     -          (za(i1,i3)*zb(i3,i2) + za(i1,i4)*zb(i4,i2))) - 
     -      t(i1,i3,i4)*za(i1,i2)*
     -       (za(i2,i3)*zb(i2,i1)*zb(i4,i1)*
     -          (za(i1,i3)*zb(i3,i2) + za(i1,i4)*zb(i4,i2)) + 
     -         za(i3,i4)*((za(i2,i3)*zb(i3,i1) + 
     -               za(i2,i4)*zb(i4,i1))*zb(i4,i2)**2 - 
     -            zb(i4,i1)**2*
     -             (za(i1,i3)*zb(i3,i2) + za(i1,i4)*zb(i4,i2))))))/
     -  (2.*s(i3,i4)*za(i1,i2)*
     -    (za(i2,i3)*zb(i3,i1) + za(i2,i4)*zb(i4,i1))*
     -    (za(i1,i3)*zb(i3,i2) + za(i1,i4)*zb(i4,i2)))

      di(d15_12)= (zb(i2,i1)*(2*zb(i4,i1)*
     -       (za(i2,i3)*zb(i3,i1) + za(i2,i4)*zb(i4,i1))*
     -       (za(i1,i3)*zb(i3,i2) + za(i1,i4)*zb(i4,i2))*
     -       (-2*mt2*za(i1,i3) - za(i1,i2)*za(i3,i4)*zb(i4,i2)) - 
     -      za(i2,i3)*(za(i1,i3)*zb(i3,i2) + za(i1,i4)*zb(i4,i2))*
     -       (t(i2,i3,i4)*za(i1,i2)*zb(i2,i1)*zb(i4,i1) + 
     -         2*(2*mt2 - za(i1,i2)*zb(i2,i1))*
     -          (za(i2,i3)*zb(i3,i1) + za(i2,i4)*zb(i4,i1))*
     -          zb(i4,i2)) + 
     -      t(i2,i3,i4)*za(i1,i2)*
     -       (-(za(i1,i3)*zb(i2,i1)*
     -            (za(i2,i3)*zb(i3,i1) + za(i2,i4)*zb(i4,i1))*
     -            zb(i4,i2)) + 
     -         za(i3,i4)*(-((za(i2,i3)*zb(i3,i1) + 
     -                 za(i2,i4)*zb(i4,i1))*zb(i4,i2)**2) + 
     -            zb(i4,i1)**2*
     -             (za(i1,i3)*zb(i3,i2) + za(i1,i4)*zb(i4,i2))))))/
     -  (2.*s(i3,i4)*za(i1,i2)*
     -    (za(i2,i3)*zb(i3,i1) + za(i2,i4)*zb(i4,i1))*
     -    (za(i1,i3)*zb(i3,i2) + za(i1,i4)*zb(i4,i2)))


      di(d15_25)= -((4*mt2*zb(i2,i1)*(za(i1,i3)**2*zb(i3,i2)*zb(i4,i1)*
     -           (za(i2,i3)*zb(i3,i1) + za(i2,i4)*zb(i4,i1)) - 
     -          za(i1,i2)*za(i3,i4)*
     -           (za(i2,i3)*zb(i3,i1) + za(i2,i4)*zb(i4,i1))*
     -           zb(i4,i2)**2 + 
     -          za(i1,i3)*zb(i4,i2)*
     -           (za(i2,i3)**2*zb(i3,i1)*zb(i3,i2) + 
     -             za(i2,i4)*zb(i4,i1)*
     -              (za(i1,i4)*zb(i4,i1) + za(i2,i4)*zb(i4,i2)) + 
     -             za(i2,i3)*
     -              (za(i1,i4)*zb(i3,i1)*zb(i4,i1) + 
     -                za(i2,i4)*
     -                 (2*zb(i3,i2)*zb(i4,i1) + zb(i2,i1)*zb(i4,i3))
     -                ))))/
     -      (za(i1,i2)*(za(i2,i3)*zb(i3,i1) + za(i2,i4)*zb(i4,i1))*
     -        (za(i1,i3)*zb(i3,i2) + za(i1,i4)*zb(i4,i2))) + 
     -     (2*za(i1,i3)*(za(i1,i3)**2*zb(i3,i1)*zb(i3,i2)*zb(i4,i1)*
     -            (za(i2,i3)*zb(i3,i1) + za(i2,i4)*zb(i4,i1))*
     -            (za(i2,i3)*zb(i3,i2) + za(i2,i4)*zb(i4,i2)) + 
     -           za(i1,i3)*(za(i2,i3)**3*zb(i3,i1)*zb(i3,i2)**2*
     -               (zb(i3,i2)*zb(i4,i1) + zb(i2,i1)*zb(i4,i3)) + 
     -              za(i2,i3)**2*zb(i3,i1)*zb(i3,i2)*
     -               (za(i1,i4)*zb(i4,i1) + za(i2,i4)*zb(i4,i2))*
     -               (2*zb(i3,i2)*zb(i4,i1) + zb(i2,i1)*zb(i4,i3))
     -               + za(i2,i4)**2*zb(i4,i1)**2*zb(i4,i2)*
     -               (za(i3,i4)*zb(i3,i2)*zb(i4,i3) + 
     -                 za(i1,i4)*
     -                  (2*zb(i3,i2)*zb(i4,i1) + 
     -                    zb(i2,i1)*zb(i4,i3))) + 
     -              za(i2,i3)*za(i2,i4)*
     -               (za(i2,i4)*zb(i3,i2)*zb(i4,i1)*zb(i4,i2)*
     -                  (zb(i3,i2)*zb(i4,i1) + zb(i2,i1)*zb(i4,i3))
     -                  + za(i1,i4)*zb(i4,i1)*
     -                  (2*zb(i3,i2)*zb(i4,i1) + 
     -                     zb(i2,i1)*zb(i4,i3))**2 - 
     -                 za(i3,i4)*zb(i4,i3)*
     -                  (-(zb(i3,i2)**2*zb(i4,i1)**2) + 
     -                    zb(i2,i1)*zb(i3,i2)*zb(i4,i1)*zb(i4,i3) + 
     -                    zb(i2,i1)**2*zb(i4,i3)**2))) + 
     -           za(i1,i4)*(za(i2,i4)**2*zb(i4,i1)**2*zb(i4,i2)**2*
     -               (za(i1,i4)*zb(i4,i1) + za(i2,i4)*zb(i4,i2)) + 
     -              za(i2,i3)**3*zb(i3,i2)*
     -               (2*zb(i3,i2)**2*zb(i4,i1)**2 + 
     -                 3*zb(i2,i1)*zb(i3,i2)*zb(i4,i1)*zb(i4,i3) + 
     -                 zb(i2,i1)**2*zb(i4,i3)**2) + 
     -              za(i2,i3)*za(i2,i4)*zb(i4,i1)*zb(i4,i2)*
     -               (-(za(i3,i4)*zb(i3,i2)*zb(i4,i1)*zb(i4,i3)) + 
     -                 za(i1,i4)*zb(i4,i1)*
     -                  (2*zb(i3,i2)*zb(i4,i1) + 
     -                    zb(i2,i1)*zb(i4,i3)) + 
     -                 2*za(i2,i4)*zb(i4,i2)*
     -                  (2*zb(i3,i2)*zb(i4,i1) + 
     -                    zb(i2,i1)*zb(i4,i3))) + 
     -              za(i2,i3)**2*
     -               (za(i1,i4)*zb(i3,i2)*zb(i4,i1)**2*
     -                  (zb(i3,i2)*zb(i4,i1) + zb(i2,i1)*zb(i4,i3))
     -                  + za(i3,i4)*zb(i4,i3)*
     -                  (-(zb(i3,i2)**2*zb(i4,i1)**2) + 
     -                    zb(i2,i1)*zb(i3,i2)*zb(i4,i1)*zb(i4,i3) + 
     -                    zb(i2,i1)**2*zb(i4,i3)**2) + 
     -                 za(i2,i4)*zb(i4,i2)*
     -                  (5*zb(i3,i2)**2*zb(i4,i1)**2 + 
     -                    5*zb(i2,i1)*zb(i3,i2)*zb(i4,i1)*
     -                     zb(i4,i3) + zb(i2,i1)**2*zb(i4,i3)**2))))
     -          + za(i1,i2)*
     -         (za(i1,i4)*zb(i4,i2)*
     -            (-(za(i2,i4)*za(i3,i4)*zb(i4,i1)**2*zb(i4,i2)*
     -                 (za(i1,i4)*zb(i4,i1) + za(i2,i4)*zb(i4,i2)))
     -               + za(i2,i3)*zb(i4,i1)*zb(i4,i2)*
     -               (za(i1,i4)*
     -                  (za(i2,i4)*zb(i2,i1) - za(i3,i4)*zb(i3,i1))*
     -                  zb(i4,i1) - 
     -                 2*za(i2,i4)*za(i3,i4)*
     -                  (zb(i3,i2)*zb(i4,i1) + zb(i2,i1)*zb(i4,i3)))
     -                - za(i2,i3)**2*
     -               (-(za(i1,i4)*zb(i2,i1)*zb(i3,i1)*zb(i4,i1)*
     -                    zb(i4,i2)) + 
     -                 za(i3,i4)*
     -                  (zb(i3,i2)*zb(i4,i1) + 
     -                     zb(i2,i1)*zb(i4,i3))**2)) + 
     -           za(i1,i3)**2*zb(i3,i2)*
     -            (za(i2,i3)*zb(i3,i1) + za(i2,i4)*zb(i4,i1))*
     -            (za(i2,i3)*zb(i2,i1)*
     -               (2*zb(i3,i2)*zb(i4,i1) + zb(i2,i1)*zb(i4,i3))
     -               + zb(i4,i1)*
     -               (za(i2,i4)*zb(i2,i1)*zb(i4,i2) - 
     -                 za(i3,i4)*
     -                  (zb(i3,i2)*zb(i4,i1) + 
     -                    2*zb(i2,i1)*zb(i4,i3)))) + 
     -           za(i1,i3)*(za(i2,i4)*zb(i4,i1)**2*zb(i4,i2)*
     -               (za(i1,i4)*
     -                  (za(i2,i4)*zb(i2,i1) - 
     -                    2*za(i3,i4)*zb(i3,i1))*zb(i4,i2) + 
     -                 za(i3,i4)*zb(i3,i2)*
     -                  (za(i2,i4)*zb(i4,i2) - 
     -                    2*za(i3,i4)*zb(i4,i3))) + 
     -              za(i2,i3)**2*
     -               (za(i1,i4)*zb(i2,i1)*zb(i3,i1)*zb(i4,i2)*
     -                  (3*zb(i3,i2)*zb(i4,i1) + 
     -                    zb(i2,i1)*zb(i4,i3)) + 
     -                 za(i3,i4)*
     -                  (zb(i3,i2)**3*zb(i4,i1)**2 - 
     -                    zb(i2,i1)**2*zb(i3,i2)*zb(i4,i3)**2)) + 
     -              za(i2,i3)*
     -               (za(i1,i4)*zb(i4,i1)*zb(i4,i2)*
     -                  (-2*za(i3,i4)*zb(i3,i1)**2*zb(i4,i2) + 
     -                    za(i2,i4)*zb(i2,i1)*
     -                     (3*zb(i3,i2)*zb(i4,i1) + 
     -                       zb(i3,i1)*zb(i4,i2) + 
     -                       zb(i2,i1)*zb(i4,i3))) + 
     -                 2*za(i3,i4)*
     -                  (za(i2,i4)*zb(i3,i2)**2*zb(i4,i1)**2*
     -                     zb(i4,i2) + 
     -                    za(i3,i4)*zb(i4,i3)*
     -                     (-(zb(i3,i2)**2*zb(i4,i1)**2) + 
     -                       zb(i2,i1)*zb(i3,i2)*zb(i4,i1)*
     -                       zb(i4,i3) + zb(i2,i1)**2*zb(i4,i3)**2))
     -                 ))))/
     -      (za(i1,i2)**2*(za(i2,i3)*zb(i3,i1) + 
     -          za(i2,i4)*zb(i4,i1))*
     -        (za(i1,i3)*zb(i3,i2) + za(i1,i4)*zb(i4,i2))))/
     -  (2.*s(i3,i4))


!---- triangle coefficients 
      ci(c25_Z)=(za(i1,i3)*(za(i1,i3)*zb(i3,i1) + za(i1,i4)*zb(i4,i1))*
     -    (za(i1,i3)*zb(i3,i2)*zb(i4,i1)*
     -       (za(i2,i3)*zb(i3,i1) + za(i2,i4)*zb(i4,i1)) + 
     -      zb(i4,i2)*(za(i2,i3)**2*zb(i3,i1)*zb(i3,i2) + 
     -         za(i1,i2)*zb(i2,i1)*
     -          (za(i2,i3)*zb(i3,i1) + za(i2,i4)*zb(i4,i1)) + 
     -         za(i2,i4)*zb(i4,i1)*
     -          (za(i1,i4)*zb(i4,i1) + za(i2,i4)*zb(i4,i2)) + 
     -         za(i2,i3)*(za(i1,i4)*zb(i3,i1)*zb(i4,i1) + 
     -            za(i2,i4)*
     -             (2*zb(i3,i2)*zb(i4,i1) + zb(i2,i1)*zb(i4,i3))))))
     -   /(za(i1,i2)**2*za(i3,i4)*
     -    (za(i2,i3)*zb(i3,i1) + za(i2,i4)*zb(i4,i1))*
     -    (za(i1,i3)*zb(i3,i2) + za(i1,i4)*zb(i4,i2))*zb(i4,i3))

      ci(cZ_15)= (za(i2,i3)*(za(i2,i3)*zb(i3,i2) + za(i2,i4)*zb(i4,i2))*
     -    (za(i2,i3)*zb(i3,i1)*zb(i4,i2)*
     -       (za(i1,i3)*zb(i3,i2) + za(i1,i4)*zb(i4,i2)) + 
     -      zb(i4,i1)*(za(i1,i3)**2*zb(i3,i1)*zb(i3,i2) + 
     -         za(i1,i2)*zb(i2,i1)*
     -          (za(i1,i3)*zb(i3,i2) + za(i1,i4)*zb(i4,i2)) + 
     -         za(i1,i4)*zb(i4,i2)*
     -          (za(i1,i4)*zb(i4,i1) + za(i2,i4)*zb(i4,i2)) + 
     -         za(i1,i3)*(za(i2,i4)*zb(i3,i2)*zb(i4,i2) + 
     -            za(i1,i4)*
     -             (2*zb(i3,i1)*zb(i4,i2) - zb(i2,i1)*zb(i4,i3))))))
     -   /(za(i1,i2)**2*za(i3,i4)*
     -    (za(i2,i3)*zb(i3,i1) + za(i2,i4)*zb(i4,i1))*
     -    (za(i1,i3)*zb(i3,i2) + za(i1,i4)*zb(i4,i2))*zb(i4,i3))

      ci(cH_25)=-((za(i2,i3)*(za(i1,i2)**2*zb(i2,i1)**2*zb(i4,i1) + 
     -        (za(i2,i3)*zb(i3,i2) + za(i2,i4)*zb(i4,i2))*
     -         (za(i1,i3)*zb(i3,i1)*zb(i4,i1) + 
     -           za(i1,i4)*zb(i4,i1)**2 + 
     -           (za(i2,i3)*zb(i3,i1) + za(i2,i4)*zb(i4,i1))*
     -            zb(i4,i2)) + 
     -        za(i1,i2)*zb(i2,i1)*
     -         (za(i1,i3)*zb(i3,i1)*zb(i4,i1) + 
     -           za(i1,i4)*zb(i4,i1)**2 + 
     -           2*za(i2,i3)*zb(i3,i1)*zb(i4,i2) + 
     -           2*za(i2,i4)*zb(i4,i1)*zb(i4,i2) - 
     -           za(i2,i3)*zb(i2,i1)*zb(i4,i3))))/
     -    (s(i3,i4)*za(i1,i2)**2*
     -      (za(i2,i3)*zb(i3,i1) + za(i2,i4)*zb(i4,i1))))

      ci(c15_H)= -((za(i1,i3)*(za(i1,i2)**2*zb(i2,i1)**2*zb(i4,i2) + 
     -        (za(i1,i3)*zb(i3,i1) + za(i1,i4)*zb(i4,i1))*
     -         (za(i2,i3)*zb(i3,i2)*zb(i4,i2) + 
     -           za(i2,i4)*zb(i4,i2)**2 + 
     -           zb(i4,i1)*(za(i1,i3)*zb(i3,i2) + 
     -              za(i1,i4)*zb(i4,i2))) + 
     -        za(i1,i2)*zb(i2,i1)*
     -         (2*za(i1,i3)*zb(i3,i2)*zb(i4,i1) + 
     -           za(i2,i3)*zb(i3,i2)*zb(i4,i2) + 
     -           2*za(i1,i4)*zb(i4,i1)*zb(i4,i2) + 
     -           za(i2,i4)*zb(i4,i2)**2 + 
     -           za(i1,i3)*zb(i2,i1)*zb(i4,i3))))/
     -    (s(i3,i4)*za(i1,i2)**2*
     -      (za(i1,i3)*zb(i3,i2) + za(i1,i4)*zb(i4,i2))))

      ci(c12)=  -((zb(i2,i1)*((zb(i4,i1)*
     -           (-(za(i2,i3)*zb(i2,i1)) + za(i3,i4)*zb(i4,i1)))/
     -         (za(i2,i3)*zb(i3,i1) + za(i2,i4)*zb(i4,i1)) - 
     -        (zb(i4,i2)*(za(i1,i3)*zb(i2,i1) + 
     -             za(i3,i4)*zb(i4,i2)))/
     -         (za(i1,i3)*zb(i3,i2) + za(i1,i4)*zb(i4,i2))))/
     -    s(i3,i4))

!===== debug print out info 
!      write(6,*) 'box 1*im',di(1)*im
!      write(6,*) 'box 2*im',di(2)*im 
!      write(6,*) 'box 3*im',di(3)*im 
      
!      write(6,*) 'triangle (s25,mZ)*im',ci(c25_Z)*im 
!      write(6,*) 'triangle (s25,mH)*im',ci(cH_25)*im 
!      write(6,*) 'triangle (s12,zip)*im',ci(c12)*im 
!      write(6,*) 'triangle (mH,s15)*im',ci(c15_H)*im 
!      write(6,*) 'triangle (mZ,s15)*im',ci(cZ_15)*im 
!====== debug stop

      ggHZ_pp_box=czip
!===== make amplitude 
        do i=1,Nbox 
           ggHZ_pp_box=ggHZ_pp_box+D0(i)*di(i)
        enddo
        do i=1,Ntri 
           ggHZ_pp_box=ggHZ_pp_box+C0(i)*ci(i)
        enddo

!---- debug amp check
!        write(6,*) 'total amp = ',im*ggHZ_pp_box
!---- debug stop
      return 
      end 
