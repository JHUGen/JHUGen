!---- CW Feb 2015 
!---- amplitude for ++ gluons => Z(i3,i4) + H 


      function ggHZ_mp_box(i1,i2,i3,i4,za,zb,mt2)
!---- this is helicity amplitude for 
!===== g(i1)^+ + g(i2)^+ + (Z=>(l1(i3)^- + l2(i3)^+)) + H 
!===== and the left-handed proejection of the Z is chosen
      implicit none
      include 'types.f'
      complex(kind=dp)::  ggHZ_mp_box
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
      integer:: c25_Z,cH_25,c12,c15_H,cZ_15,c12_Z_H
      parameter(c25_Z=1,cH_25=2,c12=3,c15_H=4,cZ_15=5,c12_Z_H=6)
!--- nb I swapped the notation wrt to KC for 3m triangle, so that I 
!---- can unify two rotuines and save calls to QCDLoop. 
      complex(kind=dp):: D0(Nbox),C0(Ntri)
      complex(kind=dp):: di(Nbox),ci(Ntri)

      complex(kind=dp)::  ggHZ_mp_2me_b1,ggHZ_mp_2mh_b1
      external  ggHZ_mp_2me_b1,ggHZ_mp_2mh_b1
      complex(kind=dp)::  ggHZ_mp_tri1,ggHZ_mp_tri2,ggHZ_mp_3mtri
      external  ggHZ_mp_tri1,ggHZ_mp_tri2,ggHZ_mp_3mtri
      common/ggZH_basisint/D0,C0
!$omp threadprivate(/ggZH_basisint/)

   
      t(i1,i2,i3)=s(i1,i2)+s(i1,i3)+s(i2,i3)

      ggHZ_mp_box=czip
      
     
!======= kinematic configurations
      s12=s(i1,i2)
      mZsq=s(i3,i4) 
      s25=t(i1,i3,i4)
      s25=t(i1,i3,i4)
      s15=t(i2,i3,i4) 
      mHsq=s(i1,i2)+s(i1,i3)+s(i1,i4)+s(i2,i3)+s(i2,i4)+s(i3,i4)
      
      

      di(:)=czip
      ci(:)=czip


!----- fill integrals from QCD loop : boxes

!      D0(d25_12)=qlI4(mZsq,zip,zip,mHsq,s25,s12,mt2,mt2,mt2,mt2,musq,0) 
!      D0(d15_12)=qlI4(mZsq,zip,zip,mHsq,s15,s12,mt2,mt2,mt2,mt2,musq,0) 
!      D0(d15_25)=qlI4(mHsq,zip,mZsq,zip,s15,s25,mt2,mt2,mt2,mt2,musq,0) 

!----- fill integrals from QCD loop : triangles
!      C0(c25_Z)=qlI3(s25,zip,mZsq,mt2,mt2,mt2,musq,0) 
!      C0(cH_25)=qlI3(mHsq,zip,s25,mt2,mt2,mt2,musq,0) 
!      C0(c12)=qlI3(s12,zip,zip,mt2,mt2,mt2,musq,0) 
!      C0(c15_H)=qlI3(s15,zip,mHsq,mt2,mt2,mt2,musq,0) 
!      C0(cZ_15)=qlI3(mZsq,zip,s15,mt2,mt2,mt2,musq,0) 
!      C0(c12_Z_H)=qlI3(s12,mZsq,mHsq,mt2,mt2,mt2,musq,0) 
   
!------ box coefficients 
      di(d25_12)= ggHZ_mp_2mh_b1(i1,i2,i3,i4,za,zb,mt2)
      di(d15_12)= -ggHZ_mp_2mh_b1(i2,i1,i4,i3,zb,za,mt2)
      di(d15_25)= ggHZ_mp_2me_b1(i1,i2,i3,i4,za,zb,mt2)
  

!==== triangle coefficients 

      ci(c25_Z)=ggHZ_mp_tri2(i1,i2,i3,i4,za,zb,mt2) 
      ci(cZ_15)=-ggHZ_mp_tri2(i2,i1,i4,i3,zb,za,mt2) 
      ci(cH_25)=-ggHZ_mp_tri1(i1,i2,i3,i4,za,zb,mt2) 
      ci(c15_H)=ggHZ_mp_tri1(i2,i1,i4,i3,zb,za,mt2) 
      ci(c12_Z_H)=ggHZ_mp_3mtri(i1,i2,i3,i4,za,zb)
      ci(c12)= (-(za(i1,i3)*(za(i1,i4)*za(i2,i3) + 
     -         za(i1,i3)*za(i2,i4))*zb(i2,i1)*
     -       (za(i2,i3)*zb(i3,i1) + za(i2,i4)*zb(i4,i1))*
     -       zb(i4,i3)) + 
     -    za(i1,i2)**2*za(i3,i4)*zb(i2,i1)*
     -     (-(za(i2,i3)*zb(i2,i1)) + za(i3,i4)*zb(i4,i1))*
     -     zb(i4,i3) + za(i1,i2)*
     -     (2*za(i1,i3)**2*za(i2,i4)*zb(i2,i1)*zb(i3,i1)*
     -        zb(i4,i1) + 
     -       2*za(i2,i4)**2*za(i3,i4)*zb(i4,i1)*
     -        zb(i4,i2)**2 - 
     -       za(i2,i3)**2*zb(i2,i1)*
     -        (2*za(i2,i4)*zb(i3,i2)*zb(i4,i2) + 
     -          za(i1,i4)*zb(i2,i1)*zb(i4,i3)) + 
     -       za(i1,i3)*zb(i2,i1)*
     -        (2*za(i1,i4)*za(i2,i4)*zb(i4,i1)**2 + 
     -          2*za(i2,i4)**2*zb(i4,i1)*zb(i4,i2) - 
     -          za(i2,i3)*za(i3,i4)*zb(i3,i1)*zb(i4,i3) + 
     -          za(i2,i4)*
     -           (-(za(i2,i3)*zb(i2,i1)) + 
     -             2*za(i3,i4)*zb(i4,i1))*zb(i4,i3)) + 
     -       za(i2,i3)*(za(i1,i4)*zb(i2,i1)*zb(i4,i1)*
     -           (-2*za(i2,i4)*zb(i4,i2) + 
     -             za(i3,i4)*zb(i4,i3)) - 
     -          2*za(i2,i4)*zb(i4,i2)*
     -           (za(i2,i4)*zb(i2,i1)*zb(i4,i2) + 
     -             za(i3,i4)*
     -              (-(zb(i3,i1)*zb(i4,i2)) + 
     -                2*zb(i2,i1)*zb(i4,i3))))))/
     -  (2.*za(i2,i4)*za(i3,i4)*
     -    (za(i2,i3)*zb(i3,i1) + za(i2,i4)*zb(i4,i1))**2*
     -    zb(i4,i3))

!===== debug print out info 
 !     write(6,*) 'box 1*im',di(1)*im
 !     write(6,*) 'box 2*im',di(2)*im 
 !     write(6,*) 'box 3*im',di(3)*im 
      
 !     write(6,*) 'triangle (s25,mZ)*im',ci(c25_Z)*im 
 !     write(6,*) 'triangle (s25,mH)*im',ci(cH_25)*im 
 !     write(6,*) 'triangle (s12,zip)*im',ci(c12)*im 
 !     write(6,*) 'triangle (mH,s15)*im',ci(c15_H)*im 
 !     write(6,*) 'triangle (mZ,s15)*im',ci(cZ_15)*im 
 !     write(6,*) 'triangle (s12,Z,H)*im',ci(c12_Z_H)*im 
 !     stop
!====== debug stop

      ggHZ_mp_box=czip
!===== make amplitude 
        do i=1,Nbox 
           ggHZ_mp_box=ggHZ_mp_box+D0(i)*di(i)
        enddo
        do i=1,Ntri 
           ggHZ_mp_box=ggHZ_mp_box+C0(i)*ci(i)
        enddo

!---- debug amp check
!       write(6,*) 'total amp = ',im*ggHZ_mp_box
!       stop
!---- debug stop
      return 
      end 


      function ggHZ_mp_2mh_b1(i1,i2,i3,i4,za,zb,mt2)
      implicit none
      include 'types.f'
      complex(kind=dp):: ggHZ_mp_2mh_b1
!---- coefficient of 2mh box 
      include 'constants.f'  
      include 'mxpart.f'
      include 'zprods_decl.f' 
      include 'sprods_com.f' 
      integer:: i1,i2,i3,i4 
      real(kind=dp):: mt2,t
      t(i1,i2,i3)=s(i1,i2)+s(i1,i3)+s(i2,i3)

        ggHZ_mp_2mh_b1=(-2*za(i1,i3)*(2*mt2*zb(i4,i1)*
     -        (za(i2,i3)*zb(i3,i1) + za(i2,i4)*zb(i4,i1))*
     -        (za(i1,i3)*zb(i3,i2) + za(i1,i4)*zb(i4,i2)) + 
     -       t(i1,i3,i4)*za(i1,i2)*zb(i2,i1)*
     -        (t(i1,i3,i4)*zb(i4,i1) + 
     -          (za(i2,i3)*zb(i3,i1) + za(i2,i4)*zb(i4,i1))*zb(i4,i2))) + 
     -    zb(i4,i2)*(-(t(i1,i3,i4)*za(i1,i2)*za(i3,i4)*
     -          (za(i2,i3)*zb(i3,i1) + za(i2,i4)*zb(i4,i1))*zb(i4,i2)) + 
     -       2*za(i2,i3)*(t(i1,i3,i4)**2*za(i1,i2)*zb(i2,i1) + 
     -          2*mt2*(za(i2,i3)*zb(i3,i1) + za(i2,i4)*zb(i4,i1))*
     -           (za(i1,i3)*zb(i3,i2) + za(i1,i4)*zb(i4,i2)))) + 
     -    t(i1,i3,i4)*za(i1,i3)**2*zb(i2,i1)*
     -     (za(i2,i3)*zb(i3,i1) + za(i2,i4)*zb(i4,i1))*zb(i4,i3))/
     -  (2.*s(i3,i4)*(za(i2,i3)*zb(i3,i1) + za(i2,i4)*zb(i4,i1))**2)
        
        return 
        end
      
      function ggHZ_mp_2me_b1(i1,i2,i3,i4,za,zb,mt2)
      implicit none
      include 'types.f'
      complex(kind=dp):: ggHZ_mp_2me_b1
!---- coefficient of 2mh box 
      include 'constants.f' 
      include 'mxpart.f'
    
      include 'zprods_decl.f' 
      include 'sprods_com.f' 
      integer:: i1,i2,i3,i4 
      real(kind=dp):: mt2

      ggHZ_mp_2me_b1=
     &(-(za(i1,i2)**3*za(i3,i4)**2*zb(i2,i1)**2*zb(i4,i3)*
     -       (za(i2,i3)**2*zb(i3,i1)**2*zb(i3,i2)*zb(i4,i1)*
     -          zb(i4,i2) + 
     -         2*za(i2,i3)*zb(i3,i1)*zb(i4,i1)*
     -          (za(i2,i4)*zb(i3,i1)*zb(i4,i2)**2 + 
     -            zb(i4,i3)*
     -             (za(i1,i3)*zb(i2,i1)*zb(i3,i1) + 
     -               za(i3,i4)*
     -                (-(zb(i3,i1)*zb(i4,i2)) + 
     -                  zb(i2,i1)*zb(i4,i3)))) + 
     -         za(i2,i4)*zb(i4,i1)**2*
     -          (za(i2,i4)*zb(i4,i2)*
     -             (zb(i3,i1)*zb(i4,i2) + zb(i2,i1)*zb(i4,i3))
     -              - 2*zb(i4,i3)*
     -             (-(za(i1,i3)*zb(i2,i1)*zb(i3,i1)) + 
     -               za(i3,i4)*
     -                (zb(i3,i1)*zb(i4,i2) + 
     -                  zb(i2,i1)*zb(i4,i3)))))) - 
     -    (za(i1,i3)*zb(i3,i1) + za(i1,i4)*zb(i4,i1))*
     -     (za(i2,i3)*zb(i3,i2) + za(i2,i4)*zb(i4,i2))*
     -     (-2*za(i1,i4)**2*za(i2,i3)*zb(i4,i1)**2*
     -        (za(i2,i3)**2*zb(i3,i1)*zb(i3,i2)**2*
     -           zb(i4,i1) + 
     -          za(i2,i4)**2*zb(i3,i1)*zb(i4,i1)*
     -           zb(i4,i2)**2 - 
     -          za(i2,i3)*za(i2,i4)*
     -           (zb(i3,i2)**2*zb(i4,i1)**2 - 
     -             4*zb(i3,i1)*zb(i3,i2)*zb(i4,i1)*
     -              zb(i4,i2) + zb(i3,i1)**2*zb(i4,i2)**2)) + 
     -       za(i1,i3)**2*za(i2,i4)*zb(i3,i1)*
     -        (-2*za(i2,i4)**2*zb(i4,i1)**2*zb(i4,i2)*
     -           (zb(i3,i2)*zb(i4,i1) - 2*zb(i3,i1)*zb(i4,i2))
     -            + za(i2,i3)**2*zb(i3,i1)*
     -           (zb(i3,i2)**2*zb(i4,i1)**2 + 
     -             zb(i3,i1)*zb(i3,i2)*zb(i4,i1)*zb(i4,i2) + 
     -             zb(i2,i1)*zb(i3,i2)*zb(i4,i1)*zb(i4,i3)) + 
     -          za(i2,i3)*za(i2,i4)*zb(i4,i1)*
     -           (-3*zb(i3,i2)**2*zb(i4,i1)**2 + 
     -             zb(i3,i1)*zb(i3,i2)*zb(i4,i1)*zb(i4,i2) + 
     -             zb(i3,i2)*zb(i4,i1)*
     -              (6*zb(i3,i1)*zb(i4,i2) + 
     -                zb(i2,i1)*zb(i4,i3)))) + 
     -       za(i1,i3)*za(i1,i4)*zb(i4,i1)*
     -        (-2*za(i2,i3)**3*zb(i3,i1)**3*zb(i3,i2)*
     -           zb(i4,i2) + 
     -          2*za(i2,i4)**3*zb(i3,i1)*zb(i4,i1)**2*
     -           zb(i4,i2)**2 + 
     -          za(i2,i3)**2*za(i2,i4)*zb(i3,i1)*
     -           (7*zb(i3,i2)**2*zb(i4,i1)**2 + 
     -             zb(i3,i1)*zb(i3,i2)*zb(i4,i1)*zb(i4,i2) + 
     -             zb(i3,i2)*zb(i4,i1)*
     -              (-10*zb(i3,i1)*zb(i4,i2) + 
     -                zb(i2,i1)*zb(i4,i3))) - 
     -          za(i2,i3)*za(i2,i4)**2*zb(i4,i1)*
     -           (zb(i3,i2)**2*zb(i4,i1)**2 + 
     -             zb(i3,i1)*zb(i4,i2)*
     -              (5*zb(i3,i1)*zb(i4,i2) + 
     -                zb(i2,i1)*zb(i4,i3)) - 
     -             zb(i3,i2)*zb(i4,i1)*
     -              (8*zb(i3,i1)*zb(i4,i2) + 
     -                zb(i2,i1)*zb(i4,i3))))) + 
     -    za(i1,i2)*(za(i1,i3)**3*zb(i2,i1)*zb(i3,i1)**2*
     -        (za(i2,i3)**2*zb(i3,i1)**2 - 
     -          za(i2,i4)**2*zb(i4,i1)**2)*
     -        (za(i2,i3)*zb(i3,i2) + za(i2,i4)*zb(i4,i2))*
     -        (-(zb(i3,i2)*zb(i4,i1)) + zb(i3,i1)*zb(i4,i2))
     -        - za(i1,i3)**2*
     -        (za(i2,i3)**3*zb(i3,i1)**3*zb(i3,i2)*
     -           (3*za(i1,i4)*zb(i2,i1)*zb(i4,i1)*
     -              (zb(i3,i2)*zb(i4,i1) - 
     -                zb(i3,i1)*zb(i4,i2)) + 
     -             zb(i4,i2)*
     -              (-2*za(i3,i4)*zb(i3,i1)*zb(i3,i2)*
     -                 zb(i4,i1) + 
     -                za(i2,i4)*zb(i2,i1)*
     -                 (zb(i3,i2)*zb(i4,i1) + 
     -                   zb(i3,i1)*zb(i4,i2)))) + 
     -          za(i2,i3)*za(i2,i4)**2*zb(i3,i1)*zb(i4,i1)*
     -           (za(i3,i4)*zb(i3,i1)*zb(i4,i2)*
     -              (zb(i3,i2)**2*zb(i4,i1)**2 - 
     -                3*zb(i3,i1)*zb(i3,i2)*zb(i4,i1)*
     -                 zb(i4,i2) - 4*zb(i3,i1)**2*zb(i4,i2)**2
     -                ) - 
     -             2*zb(i2,i1)**2*zb(i3,i2)*zb(i4,i1)*
     -              zb(i4,i3)*(-2*mt2 + za(i3,i4)*zb(i4,i3))
     -              + zb(i2,i1)*
     -              (zb(i3,i1)*zb(i3,i2)*zb(i4,i1)*zb(i4,i2)*
     -                 (-8*mt2 + 3*za(i1,i4)*zb(i4,i1) + 
     -                   7*za(i2,i4)*zb(i4,i2) - 
     -                   17*za(i3,i4)*zb(i4,i3)) - 
     -                zb(i3,i2)**2*zb(i4,i1)**2*
     -                 (-8*mt2 + za(i1,i4)*zb(i4,i1) + 
     -                   za(i2,i4)*zb(i4,i2) - 
     -                   4*za(i3,i4)*zb(i4,i3)) + 
     -                2*zb(i3,i1)**2*zb(i4,i2)**2*
     -                 (4*mt2 - za(i1,i4)*zb(i4,i1) + 
     -                   2*za(i3,i4)*zb(i4,i3)))) + 
     -          za(i2,i4)**3*zb(i4,i1)**2*
     -           (za(i3,i4)*zb(i3,i1)**2*zb(i4,i2)**2*
     -              (zb(i3,i2)*zb(i4,i1) - 
     -                3*zb(i3,i1)*zb(i4,i2)) + 
     -             zb(i2,i1)**2*zb(i4,i3)*
     -              (2*mt2*zb(i3,i2)*zb(i4,i1) - 
     -                za(i3,i4)*zb(i3,i1)*zb(i4,i2)*zb(i4,i3))
     -               + zb(i2,i1)*
     -              (-2*mt2*zb(i3,i2)**2*zb(i4,i1)**2 + 
     -                zb(i3,i1)**2*zb(i4,i2)**2*
     -                 (4*mt2 + za(i1,i4)*zb(i4,i1) + 
     -                   za(i2,i4)*zb(i4,i2) - 
     -                   8*za(i3,i4)*zb(i4,i3)) + 
     -                zb(i3,i1)*zb(i3,i2)*zb(i4,i1)*zb(i4,i2)*
     -                 (2*mt2 - za(i1,i4)*zb(i4,i1) + 
     -                   za(i2,i4)*zb(i4,i2) + 
     -                   3*za(i3,i4)*zb(i4,i3)))) + 
     -          za(i2,i3)**2*za(i2,i4)*zb(i3,i1)**2*
     -           (-(za(i3,i4)*zb(i3,i1)**2*zb(i4,i2)**2*
     -                (5*zb(i3,i2)*zb(i4,i1) + 
     -                  zb(i3,i1)*zb(i4,i2))) + 
     -             zb(i2,i1)**2*zb(i4,i3)*
     -              (za(i3,i4)*zb(i3,i1)*zb(i4,i2)*
     -                 zb(i4,i3) + 
     -                2*zb(i3,i2)*zb(i4,i1)*
     -                 (mt2 - za(i3,i4)*zb(i4,i3))) + 
     -             zb(i2,i1)*
     -              (-(zb(i3,i1)**2*zb(i4,i2)**2*
     -                   (-4*mt2 + 3*za(i1,i4)*zb(i4,i1) + 
     -                     za(i2,i4)*zb(i4,i2))) + 
     -                2*zb(i3,i2)**2*zb(i4,i1)**2*
     -                 (5*mt2 + za(i1,i4)*zb(i4,i1) - 
     -                   4*za(i3,i4)*zb(i4,i3)) + 
     -                zb(i3,i1)*zb(i3,i2)*zb(i4,i1)*zb(i4,i2)*
     -                 (-10*mt2 + za(i1,i4)*zb(i4,i1) + 
     -                   7*za(i2,i4)*zb(i4,i2) + 
     -                   5*za(i3,i4)*zb(i4,i3))))) + 
     -       za(i1,i4)**2*zb(i4,i1)*
     -        (2*za(i2,i4)**3*za(i3,i4)*zb(i3,i1)*
     -           zb(i4,i1)**3*zb(i4,i2)**3 + 
     -          za(i2,i3)**4*zb(i2,i1)*zb(i3,i1)**2*zb(i3,i2)*
     -           zb(i4,i2)*
     -           (zb(i3,i2)*zb(i4,i1) + zb(i3,i1)*zb(i4,i2))
     -           + za(i2,i3)**3*zb(i3,i1)*
     -           (2*za(i3,i4)*zb(i3,i1)*zb(i3,i2)**2*
     -              zb(i4,i1)**2*zb(i4,i2) - 
     -             zb(i2,i1)*
     -              (zb(i3,i1)**2*zb(i4,i2)**2*
     -                 (4*mt2 + za(i2,i4)*zb(i4,i2)) - 
     -                zb(i3,i1)*zb(i3,i2)*zb(i4,i1)*zb(i4,i2)*
     -                 (8*mt2 + 7*za(i2,i4)*zb(i4,i2)) + 
     -                6*za(i3,i4)*zb(i3,i2)**2*zb(i4,i1)**2*
     -                 zb(i4,i3))) + 
     -          za(i2,i3)*za(i2,i4)**2*zb(i4,i1)**2*zb(i4,i2)*
     -           (2*za(i3,i4)*zb(i3,i1)*zb(i4,i2)*
     -              (2*zb(i3,i2)*zb(i4,i1) + 
     -                zb(i3,i1)*zb(i4,i2)) + 
     -             zb(i2,i1)*
     -              (zb(i3,i1)*zb(i4,i2)*
     -                 (8*mt2 + za(i2,i4)*zb(i4,i2) - 
     -                   10*za(i3,i4)*zb(i4,i3)) + 
     -                zb(i3,i2)*zb(i4,i1)*
     -                 (-4*mt2 + za(i2,i4)*zb(i4,i2) + 
     -                   4*za(i3,i4)*zb(i4,i3)))) + 
     -          za(i2,i3)**2*za(i2,i4)*zb(i4,i1)*
     -           (2*za(i3,i4)*zb(i3,i1)*zb(i3,i2)*zb(i4,i1)*
     -              zb(i4,i2)*
     -              (zb(i3,i2)*zb(i4,i1) + 
     -                2*zb(i3,i1)*zb(i4,i2)) + 
     -             zb(i2,i1)*
     -              (zb(i3,i1)*zb(i3,i2)*zb(i4,i1)*zb(i4,i2)*
     -                 (4*mt2 + 7*za(i2,i4)*zb(i4,i2) - 
     -                   20*za(i3,i4)*zb(i4,i3)) + 
     -                2*zb(i3,i1)**2*zb(i4,i2)**2*
     -                 (2*mt2 + za(i3,i4)*zb(i4,i3)) + 
     -                zb(i3,i2)**2*zb(i4,i1)**2*
     -                 (-(za(i2,i4)*zb(i4,i2)) + 
     -                   6*za(i3,i4)*zb(i4,i3))))) + 
     -       za(i1,i3)*za(i1,i4)*
     -        (za(i2,i3)**4*zb(i2,i1)*zb(i3,i1)**3*zb(i3,i2)*
     -           zb(i4,i2)*
     -           (zb(i3,i2)*zb(i4,i1) + zb(i3,i1)*zb(i4,i2))
     -           - za(i2,i4)**3*zb(i4,i1)**3*zb(i4,i2)*
     -           (za(i3,i4)*zb(i3,i1)*zb(i4,i2)*
     -              (zb(i3,i2)*zb(i4,i1) - 
     -                5*zb(i3,i1)*zb(i4,i2)) - 
     -             za(i3,i4)*zb(i2,i1)**2*zb(i4,i3)**2 + 
     -             zb(i2,i1)*
     -              (zb(i3,i1)*zb(i4,i2)*
     -                 (8*mt2 + za(i2,i4)*zb(i4,i2) - 
     -                   6*za(i3,i4)*zb(i4,i3)) + 
     -                zb(i3,i2)*zb(i4,i1)*
     -                 (-4*mt2 + za(i2,i4)*zb(i4,i2) + 
     -                   za(i3,i4)*zb(i4,i3)))) + 
     -          za(i2,i3)**3*zb(i3,i1)**2*
     -           (za(i3,i4)*zb(i3,i1)*zb(i3,i2)*zb(i4,i1)*
     -              zb(i4,i2)*
     -              (3*zb(i3,i2)*zb(i4,i1) + 
     -                zb(i3,i1)*zb(i4,i2)) + 
     -             2*mt2*zb(i2,i1)**2*zb(i3,i1)*zb(i4,i2)*
     -              zb(i4,i3) + 
     -             zb(i2,i1)*
     -              (zb(i3,i1)**2*zb(i4,i2)**2*
     -                 (2*mt2 - za(i2,i4)*zb(i4,i2)) + 
     -                zb(i3,i1)*zb(i3,i2)*zb(i4,i1)*zb(i4,i2)*
     -                 (-6*mt2 + 2*za(i1,i4)*zb(i4,i1) + 
     -                   6*za(i2,i4)*zb(i4,i2) - 
     -                   3*za(i3,i4)*zb(i4,i3)) - 
     -                zb(i3,i2)**2*zb(i4,i1)**2*
     -                 (-8*mt2 + 2*za(i1,i4)*zb(i4,i1) + 
     -                   za(i2,i4)*zb(i4,i2) + 
     -                   4*za(i3,i4)*zb(i4,i3)))) + 
     -          za(i2,i3)**2*za(i2,i4)*zb(i3,i1)*zb(i4,i1)*
     -           (za(i3,i4)*zb(i3,i1)*zb(i4,i2)*
     -              (2*zb(i3,i2)**2*zb(i4,i1)**2 + 
     -                9*zb(i3,i1)*zb(i3,i2)*zb(i4,i1)*
     -                 zb(i4,i2) + zb(i3,i1)**2*zb(i4,i2)**2)
     -              + zb(i2,i1)**2*zb(i4,i3)*
     -              (2*za(i3,i4)*zb(i3,i2)*zb(i4,i1)*
     -                 zb(i4,i3) + 
     -                zb(i3,i1)*zb(i4,i2)*
     -                 (4*mt2 - za(i3,i4)*zb(i4,i3))) + 
     -             zb(i2,i1)*
     -              (zb(i3,i1)**2*zb(i4,i2)**2*
     -                 (8*mt2 + 2*za(i1,i4)*zb(i4,i1) + 
     -                   za(i2,i4)*zb(i4,i2)) + 
     -                zb(i3,i2)**2*zb(i4,i1)**2*
     -                 (4*mt2 - 2*za(i1,i4)*zb(i4,i1) - 
     -                   za(i2,i4)*zb(i4,i2) + 
     -                   18*za(i3,i4)*zb(i4,i3)) - 
     -                zb(i3,i1)*zb(i3,i2)*zb(i4,i1)*zb(i4,i2)*
     -                 (8*mt2 + 27*za(i3,i4)*zb(i4,i3)))) + 
     -          za(i2,i3)*za(i2,i4)**2*zb(i4,i1)**2*
     -           (za(i3,i4)*zb(i3,i1)*zb(i4,i2)*
     -              (-(zb(i3,i2)**2*zb(i4,i1)**2) + 
     -                7*zb(i3,i1)*zb(i3,i2)*zb(i4,i1)*
     -                 zb(i4,i2) + 6*zb(i3,i1)**2*zb(i4,i2)**2
     -                ) + 
     -             2*zb(i2,i1)**2*zb(i4,i3)*
     -              (mt2*zb(i3,i1)*zb(i4,i2) + 
     -                za(i3,i4)*zb(i3,i2)*zb(i4,i1)*zb(i4,i3))
     -               + zb(i2,i1)*
     -              (zb(i3,i2)**2*zb(i4,i1)**2*
     -                 (-4*mt2 + za(i2,i4)*zb(i4,i2) - 
     -                   2*za(i3,i4)*zb(i4,i3)) + 
     -                zb(i3,i1)*zb(i3,i2)*zb(i4,i1)*zb(i4,i2)*
     -                 (2*mt2 - 2*za(i1,i4)*zb(i4,i1) - 
     -                   6*za(i2,i4)*zb(i4,i2) + 
     -                   23*za(i3,i4)*zb(i4,i3)) + 
     -                zb(i3,i1)**2*zb(i4,i2)**2*
     -                 (2*za(i1,i4)*zb(i4,i1) + 
     -                   za(i2,i4)*zb(i4,i2) - 
     -                   2*(mt2 + 9*za(i3,i4)*zb(i4,i3)))))))
     -     + za(i1,i2)**2*zb(i2,i1)*
     -     (-(za(i1,i3)**2*zb(i2,i1)*zb(i3,i1)*
     -          (za(i2,i3)*zb(i3,i1) + za(i2,i4)*zb(i4,i1))*
     -          (za(i2,i3)*zb(i3,i1)*
     -             (za(i3,i4)*zb(i3,i1)*zb(i4,i2)*zb(i4,i3) + 
     -               zb(i3,i2)*zb(i4,i1)*
     -                (4*mt2 - 3*za(i3,i4)*zb(i4,i3))) + 
     -            za(i2,i4)*zb(i4,i1)*
     -             (-3*za(i3,i4)*zb(i3,i1)*zb(i4,i2)*
     -                zb(i4,i3) + 
     -               zb(i3,i2)*zb(i4,i1)*
     -                (4*mt2 + za(i3,i4)*zb(i4,i3))))) + 
     -       za(i1,i4)*(za(i2,i4)**2*za(i3,i4)*zb(i4,i1)**3*
     -           zb(i4,i2)*
     -           (zb(i3,i1)*zb(i4,i2)*
     -              (4*mt2 + za(i2,i4)*zb(i4,i2) - 
     -                4*za(i3,i4)*zb(i4,i3)) - 
     -             zb(i2,i1)*zb(i4,i3)*
     -              (-4*mt2 + za(i2,i4)*zb(i4,i2) + 
     -                2*za(i3,i4)*zb(i4,i3))) + 
     -          za(i2,i3)**2*zb(i3,i1)*zb(i4,i1)*
     -           (2*za(i3,i4)*
     -              (2*mt2*zb(i3,i1)**2*zb(i4,i2)**2 - 
     -                2*zb(i3,i1)*
     -                 (2*mt2*zb(i2,i1) + 
     -                   za(i3,i4)*zb(i3,i2)*zb(i4,i1))*
     -                 zb(i4,i2)*zb(i4,i3) + 
     -                3*za(i3,i4)*zb(i2,i1)*zb(i3,i2)*
     -                 zb(i4,i1)*zb(i4,i3)**2) + 
     -             za(i2,i4)*zb(i3,i1)*zb(i4,i2)**2*
     -              (za(i3,i4)*
     -                 (2*zb(i3,i2)*zb(i4,i1) + 
     -                   zb(i3,i1)*zb(i4,i2)) + 
     -                zb(i2,i1)*
     -                 (8*mt2 - 7*za(i3,i4)*zb(i4,i3)))) + 
     -          za(i2,i3)**3*zb(i3,i1)**2*zb(i4,i2)*
     -           (za(i3,i4)*zb(i3,i1)*zb(i3,i2)*zb(i4,i1)*
     -              zb(i4,i2) + 
     -             zb(i2,i1)*
     -              (-2*za(i3,i4)*zb(i3,i2)*zb(i4,i1)*
     -                 zb(i4,i3) + 
     -                zb(i3,i1)*zb(i4,i2)*
     -                 (4*mt2 - za(i3,i4)*zb(i4,i3)))) + 
     -          za(i2,i3)*za(i2,i4)*zb(i4,i1)**2*
     -           (za(i2,i4)*zb(i4,i2)*
     -              (za(i3,i4)*zb(i3,i1)*zb(i4,i2)*
     -                 (zb(i3,i2)*zb(i4,i1) + 
     -                   2*zb(i3,i1)*zb(i4,i2)) + 
     -                zb(i2,i1)*
     -                 (2*za(i3,i4)*zb(i3,i2)*zb(i4,i1)*
     -                    zb(i4,i3) + 
     -                   zb(i3,i1)*zb(i4,i2)*
     -                    (4*mt2 - 7*za(i3,i4)*zb(i4,i3)))) - 
     -             2*za(i3,i4)*
     -              (3*za(i3,i4)*zb(i2,i1)*zb(i3,i2)*
     -                 zb(i4,i1)*zb(i4,i3)**2 + 
     -                2*zb(i3,i1)**2*zb(i4,i2)**2*
     -                 (-2*mt2 + za(i3,i4)*zb(i4,i3)) + 
     -                zb(i3,i1)*zb(i4,i2)*zb(i4,i3)*
     -                 (2*za(i3,i4)*zb(i3,i2)*zb(i4,i1) + 
     -                   zb(i2,i1)*
     -                    (2*mt2 - 5*za(i3,i4)*zb(i4,i3))))))
     -        + za(i1,i3)*
     -        (za(i2,i3)**3*zb(i3,i1)**3*zb(i3,i2)*zb(i4,i2)*
     -           (za(i3,i4)*zb(i3,i1)*zb(i4,i2) + 
     -             zb(i2,i1)*(4*mt2 - za(i3,i4)*zb(i4,i3))) + 
     -          za(i2,i3)**2*zb(i3,i1)**2*
     -           (-2*za(i1,i4)*zb(i2,i1)*zb(i4,i1)*
     -              (-2*za(i3,i4)*zb(i3,i2)*zb(i4,i1)*
     -                 zb(i4,i3) + 
     -                zb(i3,i1)*zb(i4,i2)*
     -                 (2*mt2 + za(i3,i4)*zb(i4,i3))) + 
     -             za(i3,i4)*
     -              (-(za(i3,i4)*zb(i3,i1)**2*zb(i4,i2)**2*
     -                   zb(i4,i3)) + 
     -                2*zb(i2,i1)*zb(i3,i2)*zb(i4,i1)*
     -                 zb(i4,i3)*
     -                 (-4*mt2 + za(i3,i4)*zb(i4,i3)) + 
     -                zb(i3,i1)*zb(i4,i2)*
     -                 (zb(i3,i2)*zb(i4,i1)*
     -                    (4*mt2 - 3*za(i3,i4)*zb(i4,i3)) + 
     -                   zb(i2,i1)*zb(i4,i3)*
     -                    (4*mt2 + za(i3,i4)*zb(i4,i3)))) + 
     -             za(i2,i4)*zb(i4,i2)*
     -              (za(i3,i4)*zb(i3,i1)*zb(i4,i2)*
     -                 (2*zb(i3,i2)*zb(i4,i1) + 
     -                   zb(i3,i1)*zb(i4,i2)) + 
     -                zb(i2,i1)*
     -                 (-2*za(i3,i4)*zb(i3,i1)*zb(i4,i2)*
     -                    zb(i4,i3) + 
     -                   zb(i3,i2)*zb(i4,i1)*
     -                    (8*mt2 + za(i3,i4)*zb(i4,i3))))) + 
     -          za(i2,i4)**2*zb(i4,i1)**2*
     -           (2*za(i1,i4)*zb(i2,i1)*zb(i3,i1)*zb(i4,i1)*
     -              zb(i4,i2)*(-2*mt2 + za(i3,i4)*zb(i4,i3))
     -              + za(i3,i4)*
     -              (zb(i3,i1)**2*zb(i4,i2)**2*
     -                 (za(i2,i4)*zb(i4,i2) - 
     -                   5*za(i3,i4)*zb(i4,i3)) + 
     -                zb(i3,i1)*zb(i4,i2)*
     -                 (2*zb(i2,i1)*zb(i4,i3)*
     -                    (2*mt2 + za(i2,i4)*zb(i4,i2) - 
     -                     3*za(i3,i4)*zb(i4,i3)) + 
     -                   zb(i3,i2)*zb(i4,i1)*
     -                    (4*mt2 + za(i3,i4)*zb(i4,i3))) + 
     -                zb(i2,i1)*zb(i4,i3)*
     -                 (-(za(i3,i4)*zb(i2,i1)*zb(i4,i3)**2) + 
     -                   zb(i3,i2)*zb(i4,i1)*
     -                    (4*mt2 - za(i2,i4)*zb(i4,i2) + 
     -                     za(i3,i4)*zb(i4,i3))))) + 
     -          za(i2,i3)*za(i2,i4)*zb(i3,i1)*zb(i4,i1)*
     -           (4*za(i1,i4)*zb(i2,i1)*zb(i4,i1)*
     -              (-2*mt2*zb(i3,i1)*zb(i4,i2) + 
     -                za(i3,i4)*zb(i3,i2)*zb(i4,i1)*zb(i4,i3))
     -               + za(i2,i4)*zb(i4,i2)*
     -              (za(i3,i4)*zb(i3,i1)*zb(i4,i2)*
     -                 (zb(i3,i2)*zb(i4,i1) + 
     -                   2*zb(i3,i1)*zb(i4,i2)) + 
     -                zb(i2,i1)*zb(i3,i2)*zb(i4,i1)*
     -                 (4*mt2 + za(i3,i4)*zb(i4,i3))) - 
     -             za(i3,i4)*
     -              (6*za(i3,i4)*zb(i3,i1)**2*zb(i4,i2)**2*
     -                 zb(i4,i3) - 
     -                zb(i3,i1)*zb(i4,i2)*
     -                 (2*zb(i3,i2)*zb(i4,i1)*
     -                    (4*mt2 - za(i3,i4)*zb(i4,i3)) + 
     -                   zb(i2,i1)*zb(i4,i3)*
     -                    (8*mt2 + 7*za(i3,i4)*zb(i4,i3))) + 
     -                zb(i2,i1)*zb(i4,i3)*
     -                 (za(i3,i4)*zb(i2,i1)*zb(i4,i3)**2 + 
     -                   zb(i3,i2)*zb(i4,i1)*
     -                    (4*mt2 + 9*za(i3,i4)*zb(i4,i3)))))))
     -    )/(2.*za(i1,i2)**2*za(i3,i4)*zb(i2,i1)**2*zb(i3,i1)*
     -    (za(i2,i3)*zb(i3,i1) + za(i2,i4)*zb(i4,i1))**3*
     -    zb(i4,i3))

      return 
      end

      function ggHZ_mp_tri1(i1,i2,i3,i4,za,zb,mt2)
      implicit none
      include 'types.f'
      complex(kind=dp)::ggHZ_mp_tri1
!---- coefficient of 2mh box 
      include 'constants.f' 
      include 'mxpart.f'
    
      include 'zprods_decl.f' 
      include 'sprods_com.f' 
      integer:: i1,i2,i3,i4 
      real(kind=dp):: mt2,t
      t(i1,i2,i3)=s(i1,i2)+s(i1,i3)+s(i2,i3)

      ggHZ_mp_tri1=(-(za(i1,i2)**3*za(i3,i4)*zb(i2,i1)**2*zb(i4,i1)*
     -       zb(i4,i2)) + 
     -    za(i1,i2)**2*zb(i2,i1)*zb(i4,i1)*
     -     (2*za(i1,i3)**2*zb(i2,i1)*zb(i3,i1) + 
     -       za(i1,i3)*(2*za(i1,i4)*zb(i2,i1)*zb(i4,i1) + 
     -          (3*za(i2,i4)*zb(i2,i1) - 
     -             za(i3,i4)*zb(i3,i1))*zb(i4,i2)) - 
     -       zb(i4,i2)*(3*za(i1,i4)*za(i2,i3)*zb(i2,i1) + 
     -          2*za(i3,i4)*
     -           (za(i2,i3)*zb(i3,i2) + za(i2,i4)*zb(i4,i2)))
     -       ) - za(i1,i4)*za(i2,i3)*za(i3,i4)*zb(i4,i1)*
     -     zb(i4,i2)*(za(i2,i3)*zb(i3,i2) + 
     -       za(i2,i4)*zb(i4,i2))*zb(i4,i3) + 
     -    za(i1,i3)**2*zb(i3,i1)*
     -     (za(i2,i4)**2*zb(i4,i1)*zb(i4,i2)**2 - 
     -       2*za(i2,i3)**2*zb(i2,i1)*zb(i3,i2)*zb(i4,i3) + 
     -       za(i2,i3)*za(i2,i4)*zb(i4,i2)*
     -        (zb(i3,i2)*zb(i4,i1) - 3*zb(i2,i1)*zb(i4,i3)))
     -     + za(i1,i3)*(za(i2,i4)*za(i3,i4)*zb(i4,i1)*
     -        zb(i4,i2)*(za(i2,i3)*zb(i3,i2) + 
     -          za(i2,i4)*zb(i4,i2))*zb(i4,i3) - 
     -       za(i1,i4)*za(i2,i3)*
     -        (za(i2,i4)*zb(i4,i1)*zb(i4,i2)*
     -           (zb(i3,i2)*zb(i4,i1) + 
     -             3*zb(i2,i1)*zb(i4,i3)) + 
     -          za(i2,i3)*
     -           (zb(i3,i2)**2*zb(i4,i1)**2 + 
     -             2*zb(i2,i1)*zb(i3,i2)*zb(i4,i1)*
     -              zb(i4,i3) - zb(i2,i1)**2*zb(i4,i3)**2)))
     -     - za(i1,i2)*(za(i1,i3)**2*zb(i2,i1)*zb(i3,i1)*
     -        (-3*za(i2,i4)*zb(i4,i1)*zb(i4,i2) + 
     -          za(i2,i3)*
     -           (-2*zb(i3,i2)*zb(i4,i1) + 
     -             2*zb(i2,i1)*zb(i4,i3))) + 
     -       zb(i4,i2)*(za(i3,i4)*
     -           (2*za(i2,i3)*zb(i2,i1) + 
     -             za(i3,i4)*zb(i4,i1))*
     -           (za(i2,i3)*zb(i3,i2) + za(i2,i4)*zb(i4,i2))*
     -           zb(i4,i3) + 
     -          2*za(i1,i4)*za(i2,i3)*zb(i2,i1)*
     -           (2*za(i2,i4)*zb(i4,i1)*zb(i4,i2) + 
     -             za(i2,i3)*
     -              (2*zb(i3,i2)*zb(i4,i1) - 
     -                zb(i2,i1)*zb(i4,i3)))) + 
     -       za(i1,i3)*(za(i1,i4)*zb(i2,i1)*zb(i4,i1)*
     -           (-2*za(i2,i4)*zb(i4,i1)*zb(i4,i2) + 
     -             za(i2,i3)*
     -              (-(zb(i3,i2)*zb(i4,i1)) + 
     -                3*zb(i2,i1)*zb(i4,i3))) + 
     -          zb(i4,i2)*
     -           (za(i2,i4)*
     -              (-4*za(i2,i4)*zb(i2,i1) + 
     -                za(i3,i4)*zb(i3,i1))*zb(i4,i1)*
     -              zb(i4,i2) + 
     -             za(i2,i3)*
     -              (za(i3,i4)*zb(i3,i1)*
     -                 (zb(i3,i2)*zb(i4,i1) - 
     -                   zb(i2,i1)*zb(i4,i3)) + 
     -                2*za(i2,i4)*zb(i2,i1)*
     -                 (-2*zb(i3,i2)*zb(i4,i1) + 
     -                   zb(i2,i1)*zb(i4,i3)))))))/
     -  (2.*s(i3,i4)*za(i1,i2)*zb(i2,i1)*
     -    (za(i2,i3)*zb(i3,i1) + za(i2,i4)*zb(i4,i1))**2)
      return 
      end

      function ggHZ_mp_tri2(i1,i2,i3,i4,za,zb,mt2)
      implicit none
      include 'types.f'
      complex(kind=dp):: ggHZ_mp_tri2
!---- coefficient of 2mh box 
      include 'constants.f' 
      include 'mxpart.f'
      include 'zprods_decl.f' 
      include 'sprods_com.f' 
      integer:: i1,i2,i3,i4 
      real(kind=dp):: mt2,t
      t(i1,i2,i3)=s(i1,i2)+s(i1,i3)+s(i2,i3)

      ggHZ_mp_tri2= (zb(i4,i1)*(2*za(i1,i3)**3*zb(i2,i1)*zb(i3,i1)*
     -       (za(i2,i3)*zb(i3,i1) + za(i2,i4)*zb(i4,i1)) + 
     -      2*za(i1,i4)*za(i2,i3)*za(i3,i4)*zb(i4,i1)*
     -       zb(i4,i2)*(za(i1,i2)*zb(i2,i1) + 
     -         za(i2,i3)*zb(i3,i2) + za(i2,i4)*zb(i4,i2)) + 
     -      2*za(i1,i3)*za(i3,i4)*
     -       (za(i1,i2)*zb(i2,i1)*
     -          (-(za(i1,i4)*zb(i4,i1)**2) + 
     -            za(i2,i3)*zb(i3,i1)*zb(i4,i2)) + 
     -         za(i2,i3)*
     -          (zb(i3,i1)*zb(i4,i2)*
     -             (za(i2,i3)*zb(i3,i2) + 
     -               za(i2,i4)*zb(i4,i2)) + 
     -            za(i1,i4)*zb(i4,i1)*
     -             (-(zb(i3,i2)*zb(i4,i1)) + 
     -               zb(i3,i1)*zb(i4,i2)))) + 
     -      za(i1,i3)**2*
     -       (2*za(i1,i4)*zb(i2,i1)*zb(i4,i1)*
     -          (za(i2,i3)*zb(i3,i1) + za(i2,i4)*zb(i4,i1))
     -          + za(i3,i4)*
     -          (-2*za(i1,i2)*zb(i2,i1)*zb(i3,i1)*
     -             zb(i4,i1) + 
     -            za(i2,i4)*zb(i4,i1)*
     -             (zb(i3,i2)*zb(i4,i1) - 
     -               zb(i3,i1)*zb(i4,i2) + 
     -               zb(i2,i1)*zb(i4,i3)) + 
     -            za(i2,i3)*zb(i3,i1)*
     -             (-(zb(i3,i2)*zb(i4,i1)) + 
     -               zb(i3,i1)*zb(i4,i2) + 
     -               zb(i2,i1)*zb(i4,i3))))))/
     -  (2.*za(i2,i3)*za(i3,i4)*zb(i2,i1)*
     -    (za(i2,i3)*zb(i3,i1) + za(i2,i4)*zb(i4,i1))**2*
     -    zb(i4,i3))
      
      return 
      end


       function ggHZ_mp_3mtri(i1,i2,i3,i4,za,zb) 
      implicit none 
      include 'types.f'
      complex(kind=dp):: ggHZ_mp_3mtri
      include 'constants.f' 
      include 'mxpart.f'
    
      include 'zprods_decl.f' 
      include 'sprods_com.f' 
      integer:: i1,i2,i3,i4 
      real(kind=dp):: mt2,t,gam
      real(kind=dp):: P12(4),P34(4)
      integer:: nu

      
      t(i1,i2,i3)=s(i1,i2)+s(i1,i3)+s(i2,i3)
      
!==== part 1      
      ggHZ_mp_3mtri=  (-2*((za(i1,i3)*zb(i3,i1))/2. + 
     -       (za(i2,i3)*zb(i3,i2))/2. + 
     -       (za(i1,i4)*zb(i4,i1))/2. + 
     -       (za(i2,i4)*zb(i4,i2))/2. + 
     -       Sqrt(((za(i1,i3)*zb(i3,i1))/2. + 
     -            (za(i2,i3)*zb(i3,i2))/2. + 
     -            (za(i1,i4)*zb(i4,i1))/2. + 
     -            (za(i2,i4)*zb(i4,i2))/2.)**2 - 
     -         za(i1,i2)*za(i3,i4)*zb(i2,i1)*zb(i4,i3)))
     -      **2*(1 - (za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -          zb(i4,i3))/
     -        ((za(i1,i3)*zb(i3,i1))/2. + 
     -           (za(i2,i3)*zb(i3,i2))/2. + 
     -           (za(i1,i4)*zb(i4,i1))/2. + 
     -           (za(i2,i4)*zb(i4,i2))/2. + 
     -           Sqrt(((za(i1,i3)*zb(i3,i1))/2. + 
     -                (za(i2,i3)*zb(i3,i2))/2. + 
     -                (za(i1,i4)*zb(i4,i1))/2. + 
     -                (za(i2,i4)*zb(i4,i2))/2.)**2 - 
     -             za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -              zb(i4,i3)))**2)**8*
     -    ((((za(i1,i3)*zb(i3,i1))/2. + 
     -           (za(i2,i3)*zb(i3,i2))/2. + 
     -           (za(i1,i4)*zb(i4,i1))/2. + 
     -           (za(i2,i4)*zb(i4,i2))/2. + 
     -           Sqrt(((za(i1,i3)*zb(i3,i1))/2. + 
     -                (za(i2,i3)*zb(i3,i2))/2. + 
     -                (za(i1,i4)*zb(i4,i1))/2. + 
     -                (za(i2,i4)*zb(i4,i2))/2.)**2 - 
     -             za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -              zb(i4,i3)))*
     -         (za(i1,i2)*zb(i2,i1) - 
     -           (za(i1,i2)*zb(i2,i1)*
     -              (za(i1,i3)*zb(i3,i1) + 
     -                za(i1,i4)*zb(i4,i1)))/
     -            ((za(i1,i3)*zb(i3,i1))/2. + 
     -              (za(i2,i3)*zb(i3,i2))/2. + 
     -              (za(i1,i4)*zb(i4,i1))/2. + 
     -              (za(i2,i4)*zb(i4,i2))/2. + 
     -              Sqrt(((za(i1,i3)*zb(i3,i1))/2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2.)**2 - 
     -                za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                 zb(i4,i3))))*
     -         (za(i1,i2)*zb(i2,i1)*
     -            (((za(i2,i3)*zb(i3,i1) + 
     -                   za(i2,i4)*zb(i4,i1))*
     -                 (za(i1,i3)*zb(i3,i2) + 
     -                   za(i1,i4)*zb(i4,i2)))/
     -               (1 - 
     -                  (za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3))/
     -                   ((za(i1,i3)*zb(i3,i1))/2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2. + 
     -                   Sqrt(((za(i1,i3)*zb(i3,i1))/
     -                   2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2.)**2 - 
     -                   za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3)))**2)**2 - 
     -              ((za(i1,i3)*zb(i3,i1) + 
     -                   za(i1,i4)*zb(i4,i1) - 
     -                   (za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3))/
     -                   ((za(i1,i3)*zb(i3,i1))/2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2. + 
     -                   Sqrt(((za(i1,i3)*zb(i3,i1))/
     -                   2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2.)**2 - 
     -                   za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3))))*
     -                 (za(i2,i3)*zb(i3,i2) + 
     -                   za(i2,i4)*zb(i4,i2) - 
     -                   (za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3))/
     -                   ((za(i1,i3)*zb(i3,i1))/2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2. + 
     -                   Sqrt(((za(i1,i3)*zb(i3,i1))/
     -                   2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2.)**2 - 
     -                   za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3)))))/
     -               (1 - 
     -                  (za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3))/
     -                   ((za(i1,i3)*zb(i3,i1))/2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2. + 
     -                   Sqrt(((za(i1,i3)*zb(i3,i1))/
     -                   2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2.)**2 - 
     -                   za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3)))**2)**2)*
     -            ((za(i1,i2)**2*zb(i2,i1)**2*
     -                 (za(i2,i3)*zb(i3,i1) + 
     -                   za(i2,i4)*zb(i4,i1))**2*
     -                 (-(za(i1,i2)*zb(i4,i2)) + 
     -                   (za(i1,i2)*za(i1,i3)*zb(i2,i1)*
     -                   zb(i4,i3))/
     -                   ((za(i1,i3)*zb(i3,i1))/2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2. + 
     -                   Sqrt(((za(i1,i3)*zb(i3,i1))/
     -                   2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2.)**2 - 
     -                   za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3))))*
     -                 (za(i1,i3)*zb(i3,i1) + 
     -                   za(i1,i4)*zb(i4,i1) - 
     -                   (za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3))/
     -                   ((za(i1,i3)*zb(i3,i1))/2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2. + 
     -                   Sqrt(((za(i1,i3)*zb(i3,i1))/
     -                   2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2.)**2 - 
     -                   za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3))))**2*
     -                 (za(i3,i4)*zb(i4,i1) + 
     -                   (za(i2,i3)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3))/
     -                   ((za(i1,i3)*zb(i3,i1))/2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2. + 
     -                   Sqrt(((za(i1,i3)*zb(i3,i1))/
     -                   2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2.)**2 - 
     -                   za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3)))))/
     -               (((za(i1,i3)*zb(i3,i1))/2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2. + 
     -                   Sqrt(((za(i1,i3)*zb(i3,i1))/
     -                   2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2.)**2 - 
     -                   za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3)))**2*
     -                 (1 - 
     -                   (za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3))/
     -                   ((za(i1,i3)*zb(i3,i1))/2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2. + 
     -                   Sqrt(((za(i1,i3)*zb(i3,i1))/
     -                   2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2.)**2 - 
     -                   za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3)))**2)**6) + 
     -              ((za(i2,i3)*zb(i3,i1) + 
     -                   za(i2,i4)*zb(i4,i1))**2*
     -                 (-(za(i2,i3)*zb(i2,i1)) - 
     -                   (za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i1))/
     -                   ((za(i1,i3)*zb(i3,i1))/2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2. + 
     -                   Sqrt(((za(i1,i3)*zb(i3,i1))/
     -                   2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2.)**2 - 
     -                   za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3))))*
     -                 (za(i1,i2)*zb(i2,i1) - 
     -                   (za(i1,i2)*zb(i2,i1)*
     -                   (za(i1,i3)*zb(i3,i1) + 
     -                   za(i1,i4)*zb(i4,i1)))/
     -                   ((za(i1,i3)*zb(i3,i1))/2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2. + 
     -                   Sqrt(((za(i1,i3)*zb(i3,i1))/
     -                   2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2.)**2 - 
     -                   za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3))))**2*
     -                 (-(za(i1,i3)*zb(i4,i3)) + 
     -                   (za(i1,i2)*za(i3,i4)*zb(i4,i2)*
     -                   zb(i4,i3))/
     -                   ((za(i1,i3)*zb(i3,i1))/2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2. + 
     -                   Sqrt(((za(i1,i3)*zb(i3,i1))/
     -                   2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2.)**2 - 
     -                   za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3)))))/
     -               (1 - 
     -                  (za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3))/
     -                   ((za(i1,i3)*zb(i3,i1))/2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2. + 
     -                   Sqrt(((za(i1,i3)*zb(i3,i1))/
     -                   2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2.)**2 - 
     -                   za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3)))**2)**6) - 
     -           (za(i1,i2)*zb(i2,i1)*
     -              (za(i2,i3)*zb(i3,i1) + 
     -                 za(i2,i4)*zb(i4,i1))**2*
     -              (za(i1,i3)*zb(i3,i1) + 
     -                za(i1,i4)*zb(i4,i1) - 
     -                (za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3))/
     -                 ((za(i1,i3)*zb(i3,i1))/2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2. + 
     -                   Sqrt(((za(i1,i3)*zb(i3,i1))/
     -                   2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2.)**2 - 
     -                   za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3))))*
     -              (za(i1,i2)*zb(i4,i2)*
     -                 (((za(i2,i3)*zb(i3,i1) + 
     -                   za(i2,i4)*zb(i4,i1))*
     -                   (za(i1,i3)*zb(i3,i2) + 
     -                   za(i1,i4)*zb(i4,i2))*
     -                   (-(za(i2,i3)*zb(i2,i1)) - 
     -                   (za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i1))/
     -                   ((za(i1,i3)*zb(i3,i1))/2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2. + 
     -                   Sqrt(((za(i1,i3)*zb(i3,i1))/
     -                   2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2.)**2 - 
     -                   za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3))))*
     -                   (za(i1,i2)*zb(i2,i1) - 
     -                   (za(i1,i2)*zb(i2,i1)*
     -                   (za(i1,i3)*zb(i3,i1) + 
     -                   za(i1,i4)*zb(i4,i1)))/
     -                   ((za(i1,i3)*zb(i3,i1))/2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2. + 
     -                   Sqrt(((za(i1,i3)*zb(i3,i1))/
     -                   2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2.)**2 - 
     -                   za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3)))))/
     -                   (1 - 
     -                   (za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3))/
     -                   ((za(i1,i3)*zb(i3,i1))/2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2. + 
     -                   Sqrt(((za(i1,i3)*zb(i3,i1))/
     -                   2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2.)**2 - 
     -                   za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3)))**2)**4 + 
     -                   (za(i1,i2)**2*zb(i2,i1)**2*
     -                   (za(i2,i3)*zb(i3,i1) + 
     -                   za(i2,i4)*zb(i4,i1))*
     -                   (za(i1,i3)*zb(i3,i2) + 
     -                   za(i1,i4)*zb(i4,i2))*
     -                   (za(i1,i3)*zb(i3,i1) + 
     -                   za(i1,i4)*zb(i4,i1) - 
     -                   (za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3))/
     -                   ((za(i1,i3)*zb(i3,i1))/2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2. + 
     -                   Sqrt(((za(i1,i3)*zb(i3,i1))/
     -                   2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2.)**2 - 
     -                   za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3))))*
     -                   (za(i3,i4)*zb(i4,i1) + 
     -                   (za(i2,i3)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3))/
     -                   ((za(i1,i3)*zb(i3,i1))/2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2. + 
     -                   Sqrt(((za(i1,i3)*zb(i3,i1))/
     -                   2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2.)**2 - 
     -                   za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3)))))/
     -                   (((za(i1,i3)*zb(i3,i1))/2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2. + 
     -                   Sqrt(((za(i1,i3)*zb(i3,i1))/
     -                   2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2.)**2 - 
     -                   za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3)))**2*
     -                   (1 - 
     -                   (za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3))/
     -                   ((za(i1,i3)*zb(i3,i1))/2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2. + 
     -                   Sqrt(((za(i1,i3)*zb(i3,i1))/
     -                   2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2.)**2 - 
     -                   za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3)))**2)**4)) + 
     -                za(i1,i3)*zb(i2,i1)*
     -                 (-((za(i1,i2)*zb(i2,i1)*
     -                   (za(i2,i3)*zb(i3,i1) + 
     -                   za(i2,i4)*zb(i4,i1))*
     -                   (-(za(i1,i2)*zb(i4,i2)) + 
     -                   (za(i1,i2)*za(i1,i3)*zb(i2,i1)*
     -                   zb(i4,i3))/
     -                   ((za(i1,i3)*zb(i3,i1))/2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2. + 
     -                   Sqrt(((za(i1,i3)*zb(i3,i1))/
     -                   2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2.)**2 - 
     -                   za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3))))*
     -                   (za(i1,i3)*zb(i3,i1) + 
     -                   za(i1,i4)*zb(i4,i1) - 
     -                   (za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3))/
     -                   ((za(i1,i3)*zb(i3,i1))/2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2. + 
     -                   Sqrt(((za(i1,i3)*zb(i3,i1))/
     -                   2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2.)**2 - 
     -                   za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3))))**2)/
     -                   (((za(i1,i3)*zb(i3,i1))/2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2. + 
     -                   Sqrt(((za(i1,i3)*zb(i3,i1))/
     -                   2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2.)**2 - 
     -                   za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3)))*
     -                   (1 - 
     -                   (za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3))/
     -                   ((za(i1,i3)*zb(i3,i1))/2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2. + 
     -                   Sqrt(((za(i1,i3)*zb(i3,i1))/
     -                   2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2.)**2 - 
     -                   za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3)))**2)**4)) + 
     -                   ((za(i2,i3)*zb(i3,i1) + 
     -                   za(i2,i4)*zb(i4,i1))*
     -                   (za(i1,i2)*zb(i2,i1) - 
     -                   (za(i1,i2)*zb(i2,i1)*
     -                   (za(i1,i3)*zb(i3,i1) + 
     -                   za(i1,i4)*zb(i4,i1)))/
     -                   ((za(i1,i3)*zb(i3,i1))/2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2. + 
     -                   Sqrt(((za(i1,i3)*zb(i3,i1))/
     -                   2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2.)**2 - 
     -                   za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3))))**2*
     -                   (-(za(i1,i3)*zb(i4,i3)) + 
     -                   (za(i1,i2)*za(i3,i4)*zb(i4,i2)*
     -                   zb(i4,i3))/
     -                   ((za(i1,i3)*zb(i3,i1))/2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2. + 
     -                   Sqrt(((za(i1,i3)*zb(i3,i1))/
     -                   2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2.)**2 - 
     -                   za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3)))))/
     -                   (1 - 
     -                   (za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3))/
     -                   ((za(i1,i3)*zb(i3,i1))/2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2. + 
     -                   Sqrt(((za(i1,i3)*zb(i3,i1))/
     -                   2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2.)**2 - 
     -                   za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3)))**2)**4)))/
     -            (1 - 
     -               (za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                  zb(i4,i3))/
     -                ((za(i1,i3)*zb(i3,i1))/2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2. + 
     -                   Sqrt(((za(i1,i3)*zb(i3,i1))/
     -                   2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2.)**2 - 
     -                   za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3)))**2)**3))/
     -       (1 - (za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -            zb(i4,i3))/
     -          ((za(i1,i3)*zb(i3,i1))/2. + 
     -             (za(i2,i3)*zb(i3,i2))/2. + 
     -             (za(i1,i4)*zb(i4,i1))/2. + 
     -             (za(i2,i4)*zb(i4,i2))/2. + 
     -             Sqrt(((za(i1,i3)*zb(i3,i1))/2. + 
     -                  (za(i2,i3)*zb(i3,i2))/2. + 
     -                  (za(i1,i4)*zb(i4,i1))/2. + 
     -                  (za(i2,i4)*zb(i4,i2))/2.)**2 - 
     -               za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                zb(i4,i3)))**2) + 
     -      za(i3,i4)*zb(i4,i3)*
     -       (-((za(i1,i2)**4*zb(i2,i1)**4*
     -              (za(i2,i3)*zb(i3,i1) + 
     -                 za(i2,i4)*zb(i4,i1))**3*
     -              (za(i1,i3)*zb(i3,i2) + 
     -                za(i1,i4)*zb(i4,i2))*
     -              (za(i1,i2)*zb(i2,i1) + 
     -                (za(i1,i3)*zb(i3,i1))/2. + 
     -                (za(i2,i3)*zb(i3,i2))/2. + 
     -                (za(i1,i4)*zb(i4,i1))/2. + 
     -                (za(i2,i4)*zb(i4,i2))/2. + 
     -                Sqrt(((za(i1,i3)*zb(i3,i1))/2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2.)**2 - 
     -                  za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3)))*
     -              (-(za(i1,i2)*zb(i4,i2)) + 
     -                (za(i1,i2)*za(i1,i3)*zb(i2,i1)*
     -                   zb(i4,i3))/
     -                 ((za(i1,i3)*zb(i3,i1))/2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2. + 
     -                   Sqrt(((za(i1,i3)*zb(i3,i1))/
     -                   2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2.)**2 - 
     -                   za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3))))*
     -              (za(i1,i3)*zb(i3,i1) + 
     -                 za(i1,i4)*zb(i4,i1) - 
     -                 (za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3))/
     -                  ((za(i1,i3)*zb(i3,i1))/2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2. + 
     -                   Sqrt(((za(i1,i3)*zb(i3,i1))/
     -                   2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2.)**2 - 
     -                   za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3))))**3*
     -              (za(i3,i4)*zb(i4,i1) + 
     -                (za(i2,i3)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3))/
     -                 ((za(i1,i3)*zb(i3,i1))/2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2. + 
     -                   Sqrt(((za(i1,i3)*zb(i3,i1))/
     -                   2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2.)**2 - 
     -                   za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3)))))/
     -            (((za(i1,i3)*zb(i3,i1))/2. + 
     -                 (za(i2,i3)*zb(i3,i2))/2. + 
     -                 (za(i1,i4)*zb(i4,i1))/2. + 
     -                 (za(i2,i4)*zb(i4,i2))/2. + 
     -                 Sqrt(((za(i1,i3)*zb(i3,i1))/2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2.)**2 - 
     -                   za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3)))**4*
     -              (1 - 
     -                 (za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3))/
     -                  ((za(i1,i3)*zb(i3,i1))/2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2. + 
     -                   Sqrt(((za(i1,i3)*zb(i3,i1))/
     -                   2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2.)**2 - 
     -                   za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3)))**2)**9)) + 
     -         (za(i1,i2)*zb(i2,i1)*
     -            (za(i2,i3)*zb(i3,i1) + 
     -               za(i2,i4)*zb(i4,i1))**3*
     -            (-(za(i2,i3)*zb(i2,i1)) - 
     -              (za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                 zb(i4,i1))/
     -               ((za(i1,i3)*zb(i3,i1))/2. + 
     -                 (za(i2,i3)*zb(i3,i2))/2. + 
     -                 (za(i1,i4)*zb(i4,i1))/2. + 
     -                 (za(i2,i4)*zb(i4,i2))/2. + 
     -                 Sqrt(((za(i1,i3)*zb(i3,i1))/2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2.)**2 - 
     -                   za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3))))*
     -            (za(i1,i2)*zb(i2,i1) - 
     -               (za(i1,i2)*zb(i2,i1)*
     -                  (za(i1,i3)*zb(i3,i1) + 
     -                   za(i1,i4)*zb(i4,i1)))/
     -                ((za(i1,i3)*zb(i3,i1))/2. + 
     -                  (za(i2,i3)*zb(i3,i2))/2. + 
     -                  (za(i1,i4)*zb(i4,i1))/2. + 
     -                  (za(i2,i4)*zb(i4,i2))/2. + 
     -                  Sqrt(((za(i1,i3)*zb(i3,i1))/
     -                  2. + (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2.)**2 - 
     -                   za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3))))**2*
     -            (za(i1,i3)*zb(i3,i1) + 
     -              za(i1,i4)*zb(i4,i1) - 
     -              (za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                 zb(i4,i3))/
     -               ((za(i1,i3)*zb(i3,i1))/2. + 
     -                 (za(i2,i3)*zb(i3,i2))/2. + 
     -                 (za(i1,i4)*zb(i4,i1))/2. + 
     -                 (za(i2,i4)*zb(i4,i2))/2. + 
     -                 Sqrt(((za(i1,i3)*zb(i3,i1))/2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2.)**2 - 
     -                   za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3))))*
     -            ((za(i1,i2)**2*zb(i2,i1)*zb(i4,i2)*
     -                 (za(i1,i3)*zb(i3,i2) + 
     -                   za(i1,i4)*zb(i4,i2)))/
     -               (1 - 
     -                 (za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3))/
     -                  ((za(i1,i3)*zb(i3,i1))/2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2. + 
     -                   Sqrt(((za(i1,i3)*zb(i3,i1))/
     -                   2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2.)**2 - 
     -                   za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3)))**2) - 
     -              (za(i1,i2)*zb(i2,i1)*
     -                 (za(i1,i3)*zb(i3,i2) + 
     -                   za(i1,i4)*zb(i4,i2))*
     -                 (za(i1,i2)*zb(i2,i1) + 
     -                   (za(i1,i3)*zb(i3,i1))/2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2. + 
     -                   Sqrt(((za(i1,i3)*zb(i3,i1))/
     -                   2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2.)**2 - 
     -                   za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3)))*
     -                 (-(za(i1,i3)*zb(i4,i3)) + 
     -                   (za(i1,i2)*za(i3,i4)*zb(i4,i2)*
     -                   zb(i4,i3))/
     -                   ((za(i1,i3)*zb(i3,i1))/2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2. + 
     -                   Sqrt(((za(i1,i3)*zb(i3,i1))/
     -                   2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2.)**2 - 
     -                   za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3)))))/
     -               (((za(i1,i3)*zb(i3,i1))/2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2. + 
     -                   Sqrt(((za(i1,i3)*zb(i3,i1))/
     -                   2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2.)**2 - 
     -                   za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3)))*
     -                 (1 - 
     -                   (za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3))/
     -                   ((za(i1,i3)*zb(i3,i1))/2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2. + 
     -                   Sqrt(((za(i1,i3)*zb(i3,i1))/
     -                   2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2.)**2 - 
     -                   za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3)))**2)**2)))/
     -          (((za(i1,i3)*zb(i3,i1))/2. + 
     -              (za(i2,i3)*zb(i3,i2))/2. + 
     -              (za(i1,i4)*zb(i4,i1))/2. + 
     -              (za(i2,i4)*zb(i4,i2))/2. + 
     -              Sqrt(((za(i1,i3)*zb(i3,i1))/2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2.)**2 - 
     -                za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                 zb(i4,i3)))*
     -            (1 - 
     -               (za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                  zb(i4,i3))/
     -                ((za(i1,i3)*zb(i3,i1))/2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2. + 
     -                   Sqrt(((za(i1,i3)*zb(i3,i1))/
     -                   2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2.)**2 - 
     -                   za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3)))**2)**7) + 
     -         ((za(i2,i3)*zb(i3,i1) + 
     -               za(i2,i4)*zb(i4,i1))**2*
     -            (za(i1,i2)*zb(i2,i1) - 
     -               (za(i1,i2)*zb(i2,i1)*
     -                  (za(i1,i3)*zb(i3,i1) + 
     -                   za(i1,i4)*zb(i4,i1)))/
     -                ((za(i1,i3)*zb(i3,i1))/2. + 
     -                  (za(i2,i3)*zb(i3,i2))/2. + 
     -                  (za(i1,i4)*zb(i4,i1))/2. + 
     -                  (za(i2,i4)*zb(i4,i2))/2. + 
     -                  Sqrt(((za(i1,i3)*zb(i3,i1))/
     -                  2. + (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2.)**2 - 
     -                   za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3))))**3*
     -            (-(za(i1,i3)*zb(i4,i3)) + 
     -              (za(i1,i2)*za(i3,i4)*zb(i4,i2)*
     -                 zb(i4,i3))/
     -               ((za(i1,i3)*zb(i3,i1))/2. + 
     -                 (za(i2,i3)*zb(i3,i2))/2. + 
     -                 (za(i1,i4)*zb(i4,i1))/2. + 
     -                 (za(i2,i4)*zb(i4,i2))/2. + 
     -                 Sqrt(((za(i1,i3)*zb(i3,i1))/2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2.)**2 - 
     -                   za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3))))*
     -            ((za(i1,i2)*zb(i2,i1)*
     -                 (za(i2,i3)*zb(i3,i1) + 
     -                   za(i2,i4)*zb(i4,i1))*
     -                 (za(i1,i3)*zb(i3,i2) + 
     -                   za(i1,i4)*zb(i4,i2))*
     -                 (-(za(i2,i3)*zb(i2,i1)) - 
     -                   (za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i1))/
     -                   ((za(i1,i3)*zb(i3,i1))/2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2. + 
     -                   Sqrt(((za(i1,i3)*zb(i3,i1))/
     -                   2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2.)**2 - 
     -                   za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3)))))/
     -               (1 - 
     -                  (za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3))/
     -                   ((za(i1,i3)*zb(i3,i1))/2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2. + 
     -                   Sqrt(((za(i1,i3)*zb(i3,i1))/
     -                   2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2.)**2 - 
     -                   za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3)))**2)**3 + 
     -              ((za(i1,i3)*zb(i3,i1) + 
     -                   za(i1,i4)*zb(i4,i1) - 
     -                   (za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3))/
     -                   ((za(i1,i3)*zb(i3,i1))/2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2. + 
     -                   Sqrt(((za(i1,i3)*zb(i3,i1))/
     -                   2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2.)**2 - 
     -                   za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3))))*
     -                 (((za(i1,i2)*zb(i2,i1) + 
     -                   (za(i1,i3)*zb(i3,i1))/2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2. + 
     -                   Sqrt(((za(i1,i3)*zb(i3,i1))/
     -                   2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2.)**2 - 
     -                   za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3)))*
     -                   (-(za(i2,i3)*zb(i2,i1)) - 
     -                   (za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i1))/
     -                   ((za(i1,i3)*zb(i3,i1))/2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2. + 
     -                   Sqrt(((za(i1,i3)*zb(i3,i1))/
     -                   2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2.)**2 - 
     -                   za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3))))*
     -                   (za(i1,i2)*zb(i2,i1) - 
     -                   (za(i1,i2)*zb(i2,i1)*
     -                   (za(i2,i3)*zb(i3,i2) + 
     -                   za(i2,i4)*zb(i4,i2)))/
     -                   ((za(i1,i3)*zb(i3,i1))/2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2. + 
     -                   Sqrt(((za(i1,i3)*zb(i3,i1))/
     -                   2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2.)**2 - 
     -                   za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3)))))/
     -                   (1 - 
     -                   (za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3))/
     -                   ((za(i1,i3)*zb(i3,i1))/2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2. + 
     -                   Sqrt(((za(i1,i3)*zb(i3,i1))/
     -                   2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2.)**2 - 
     -                   za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3)))**2)**2 - 
     -                   za(i1,i2)*zb(i2,i1)*
     -                   (-((za(i1,i2)*za(i1,i3)*
     -                   zb(i2,i1)**2*
     -                   (za(i2,i3)*zb(i3,i1) + 
     -                   za(i2,i4)*zb(i4,i1)))/
     -                   (((za(i1,i3)*zb(i3,i1))/2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2. + 
     -                   Sqrt(((za(i1,i3)*zb(i3,i1))/
     -                   2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2.)**2 - 
     -                   za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3)))*
     -                   (1 - 
     -                   (za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3))/
     -                   ((za(i1,i3)*zb(i3,i1))/2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2. + 
     -                   Sqrt(((za(i1,i3)*zb(i3,i1))/
     -                   2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2.)**2 - 
     -                   za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3)))**2))) + 
     -                   ((-(za(i2,i3)*zb(i2,i1)) - 
     -                   (za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i1))/
     -                   ((za(i1,i3)*zb(i3,i1))/2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2. + 
     -                   Sqrt(((za(i1,i3)*zb(i3,i1))/
     -                   2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2.)**2 - 
     -                   za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3))))*
     -                   (za(i2,i3)*zb(i3,i2) + 
     -                   za(i2,i4)*zb(i4,i2) - 
     -                   (za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3))/
     -                   ((za(i1,i3)*zb(i3,i1))/2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2. + 
     -                   Sqrt(((za(i1,i3)*zb(i3,i1))/
     -                   2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2.)**2 - 
     -                   za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3)))))/
     -                   (1 - 
     -                   (za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3))/
     -                   ((za(i1,i3)*zb(i3,i1))/2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2. + 
     -                   Sqrt(((za(i1,i3)*zb(i3,i1))/
     -                   2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2.)**2 - 
     -                   za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3)))**2)**2)))/
     -               (1 - 
     -                 (za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3))/
     -                  ((za(i1,i3)*zb(i3,i1))/2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2. + 
     -                   Sqrt(((za(i1,i3)*zb(i3,i1))/
     -                   2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2.)**2 - 
     -                   za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3)))**2)))/
     -          (1 - (za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                zb(i4,i3))/
     -              ((za(i1,i3)*zb(i3,i1))/2. + 
     -                 (za(i2,i3)*zb(i3,i2))/2. + 
     -                 (za(i1,i4)*zb(i4,i1))/2. + 
     -                 (za(i2,i4)*zb(i4,i2))/2. + 
     -                 Sqrt(((za(i1,i3)*zb(i3,i1))/2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2.)**2 - 
     -                   za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3)))**2)**6 + 
     -         (za(i1,i2)**2*zb(i2,i1)**2*
     -            (za(i2,i3)*zb(i3,i1) + 
     -               za(i2,i4)*zb(i4,i1))**2*
     -            (za(i1,i2)*zb(i2,i1) - 
     -              (za(i1,i2)*zb(i2,i1)*
     -                 (za(i1,i3)*zb(i3,i1) + 
     -                   za(i1,i4)*zb(i4,i1)))/
     -               ((za(i1,i3)*zb(i3,i1))/2. + 
     -                 (za(i2,i3)*zb(i3,i2))/2. + 
     -                 (za(i1,i4)*zb(i4,i1))/2. + 
     -                 (za(i2,i4)*zb(i4,i2))/2. + 
     -                 Sqrt(((za(i1,i3)*zb(i3,i1))/2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2.)**2 - 
     -                   za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3))))*
     -            (za(i1,i3)*zb(i3,i1) + 
     -               za(i1,i4)*zb(i4,i1) - 
     -               (za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                  zb(i4,i3))/
     -                ((za(i1,i3)*zb(i3,i1))/2. + 
     -                  (za(i2,i3)*zb(i3,i2))/2. + 
     -                  (za(i1,i4)*zb(i4,i1))/2. + 
     -                  (za(i2,i4)*zb(i4,i2))/2. + 
     -                  Sqrt(((za(i1,i3)*zb(i3,i1))/
     -                  2. + (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2.)**2 - 
     -                   za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3))))**2*
     -            ((za(i1,i2)**3*zb(i2,i1)**2*
     -                 (za(i2,i3)*zb(i3,i1) + 
     -                   za(i2,i4)*zb(i4,i1))*zb(i4,i2)*
     -                 (za(i1,i3)*zb(i3,i2) + 
     -                   za(i1,i4)*zb(i4,i2))*
     -                 (za(i3,i4)*zb(i4,i1) + 
     -                   (za(i2,i3)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3))/
     -                   ((za(i1,i3)*zb(i3,i1))/2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2. + 
     -                   Sqrt(((za(i1,i3)*zb(i3,i1))/
     -                   2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2.)**2 - 
     -                   za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3)))))/
     -               (((za(i1,i3)*zb(i3,i1))/2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2. + 
     -                   Sqrt(((za(i1,i3)*zb(i3,i1))/
     -                   2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2.)**2 - 
     -                   za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3)))*
     -                 (1 - 
     -                   (za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3))/
     -                   ((za(i1,i3)*zb(i3,i1))/2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2. + 
     -                   Sqrt(((za(i1,i3)*zb(i3,i1))/
     -                   2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2.)**2 - 
     -                   za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3)))**2)**3) + 
     -              ((-(za(i1,i2)*zb(i4,i2)) + 
     -                   (za(i1,i2)*za(i1,i3)*zb(i2,i1)*
     -                   zb(i4,i3))/
     -                   ((za(i1,i3)*zb(i3,i1))/2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2. + 
     -                   Sqrt(((za(i1,i3)*zb(i3,i1))/
     -                   2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2.)**2 - 
     -                   za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3))))*
     -                 ((za(i1,i2)*zb(i2,i1)*
     -                   (za(i2,i3)*zb(i3,i1) + 
     -                   za(i2,i4)*zb(i4,i1))*
     -                   (za(i1,i3)*zb(i3,i2) + 
     -                   za(i1,i4)*zb(i4,i2))*
     -                   (za(i3,i4)*zb(i4,i1) + 
     -                   (za(i2,i3)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3))/
     -                   ((za(i1,i3)*zb(i3,i1))/2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2. + 
     -                   Sqrt(((za(i1,i3)*zb(i3,i1))/
     -                   2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2.)**2 - 
     -                   za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3)))))/
     -                   (1 - 
     -                   (za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3))/
     -                   ((za(i1,i3)*zb(i3,i1))/2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2. + 
     -                   Sqrt(((za(i1,i3)*zb(i3,i1))/
     -                   2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2.)**2 - 
     -                   za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3)))**2)**3 + 
     -                   ((za(i1,i3)*zb(i3,i1) + 
     -                   za(i1,i4)*zb(i4,i1) - 
     -                   (za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3))/
     -                   ((za(i1,i3)*zb(i3,i1))/2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2. + 
     -                   Sqrt(((za(i1,i3)*zb(i3,i1))/
     -                   2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2.)**2 - 
     -                   za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3))))*
     -                   (((za(i1,i2)*zb(i2,i1) + 
     -                   (za(i1,i3)*zb(i3,i1))/2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2. + 
     -                   Sqrt(((za(i1,i3)*zb(i3,i1))/
     -                   2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2.)**2 - 
     -                   za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3)))*
     -                   (za(i1,i2)*zb(i2,i1) - 
     -                   (za(i1,i2)*zb(i2,i1)*
     -                   (za(i2,i3)*zb(i3,i2) + 
     -                   za(i2,i4)*zb(i4,i2)))/
     -                   ((za(i1,i3)*zb(i3,i1))/2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2. + 
     -                   Sqrt(((za(i1,i3)*zb(i3,i1))/
     -                   2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2.)**2 - 
     -                   za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3))))*
     -                   (za(i3,i4)*zb(i4,i1) + 
     -                   (za(i2,i3)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3))/
     -                   ((za(i1,i3)*zb(i3,i1))/2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2. + 
     -                   Sqrt(((za(i1,i3)*zb(i3,i1))/
     -                   2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2.)**2 - 
     -                   za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3)))))/
     -                   (1 - 
     -                   (za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3))/
     -                   ((za(i1,i3)*zb(i3,i1))/2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2. + 
     -                   Sqrt(((za(i1,i3)*zb(i3,i1))/
     -                   2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2.)**2 - 
     -                   za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3)))**2)**2 - 
     -                   za(i1,i2)*zb(i2,i1)*
     -                   ((za(i1,i3)*zb(i2,i1)*
     -                   (za(i2,i3)*zb(i3,i1) + 
     -                   za(i2,i4)*zb(i4,i1)))/
     -                   (1 - 
     -                   (za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3))/
     -                   ((za(i1,i3)*zb(i3,i1))/2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2. + 
     -                   Sqrt(((za(i1,i3)*zb(i3,i1))/
     -                   2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2.)**2 - 
     -                   za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3)))**2) + 
     -                   ((za(i2,i3)*zb(i3,i2) + 
     -                   za(i2,i4)*zb(i4,i2) - 
     -                   (za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3))/
     -                   ((za(i1,i3)*zb(i3,i1))/2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2. + 
     -                   Sqrt(((za(i1,i3)*zb(i3,i1))/
     -                   2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2.)**2 - 
     -                   za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3))))*
     -                   (za(i3,i4)*zb(i4,i1) + 
     -                   (za(i2,i3)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3))/
     -                   ((za(i1,i3)*zb(i3,i1))/2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2. + 
     -                   Sqrt(((za(i1,i3)*zb(i3,i1))/
     -                   2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2.)**2 - 
     -                   za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3)))))/
     -                   (1 - 
     -                   (za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3))/
     -                   ((za(i1,i3)*zb(i3,i1))/2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2. + 
     -                   Sqrt(((za(i1,i3)*zb(i3,i1))/
     -                   2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2.)**2 - 
     -                   za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3)))**2)**2)))/
     -                   (1 - 
     -                   (za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3))/
     -                   ((za(i1,i3)*zb(i3,i1))/2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2. + 
     -                   Sqrt(((za(i1,i3)*zb(i3,i1))/
     -                   2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2.)**2 - 
     -                   za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3)))**2)))/
     -               (1 - 
     -                 (za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3))/
     -                  ((za(i1,i3)*zb(i3,i1))/2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2. + 
     -                   Sqrt(((za(i1,i3)*zb(i3,i1))/
     -                   2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2.)**2 - 
     -                   za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3)))**2)))/
     -          (((za(i1,i3)*zb(i3,i1))/2. + 
     -               (za(i2,i3)*zb(i3,i2))/2. + 
     -               (za(i1,i4)*zb(i4,i1))/2. + 
     -               (za(i2,i4)*zb(i4,i2))/2. + 
     -               Sqrt(((za(i1,i3)*zb(i3,i1))/2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2.)**2 - 
     -                 za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                  zb(i4,i3)))**2*
     -            (1 - 
     -               (za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                  zb(i4,i3))/
     -                ((za(i1,i3)*zb(i3,i1))/2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2. + 
     -                   Sqrt(((za(i1,i3)*zb(i3,i1))/
     -                   2. + 
     -                   (za(i2,i3)*zb(i3,i2))/2. + 
     -                   (za(i1,i4)*zb(i4,i1))/2. + 
     -                   (za(i2,i4)*zb(i4,i2))/2.)**2 - 
     -                   za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -                   zb(i4,i3)))**2)**5))))/
     -  (za(i1,i2)**2*zb(i2,i1)**2*
     -    (za(i2,i3)*zb(i3,i1) + za(i2,i4)*zb(i4,i1))**
     -     4*(za(i1,i2)*zb(i2,i1) - 
     -       (za(i1,i2)*zb(i2,i1)*
     -          (za(i1,i3)*zb(i3,i1) + 
     -            za(i1,i4)*zb(i4,i1)))/
     -        ((za(i1,i3)*zb(i3,i1))/2. + 
     -          (za(i2,i3)*zb(i3,i2))/2. + 
     -          (za(i1,i4)*zb(i4,i1))/2. + 
     -          (za(i2,i4)*zb(i4,i2))/2. + 
     -          Sqrt(((za(i1,i3)*zb(i3,i1))/2. + 
     -               (za(i2,i3)*zb(i3,i2))/2. + 
     -               (za(i1,i4)*zb(i4,i1))/2. + 
     -               (za(i2,i4)*zb(i4,i2))/2.)**2 - 
     -            za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -             zb(i4,i3))))**2*
     -    (za(i1,i3)*zb(i3,i1) + za(i1,i4)*zb(i4,i1) - 
     -       (za(i1,i2)*za(i3,i4)*zb(i2,i1)*zb(i4,i3))/
     -        ((za(i1,i3)*zb(i3,i1))/2. + 
     -          (za(i2,i3)*zb(i3,i2))/2. + 
     -          (za(i1,i4)*zb(i4,i1))/2. + 
     -          (za(i2,i4)*zb(i4,i2))/2. + 
     -          Sqrt(((za(i1,i3)*zb(i3,i1))/2. + 
     -               (za(i2,i3)*zb(i3,i2))/2. + 
     -               (za(i1,i4)*zb(i4,i1))/2. + 
     -               (za(i2,i4)*zb(i4,i2))/2.)**2 - 
     -            za(i1,i2)*za(i3,i4)*zb(i2,i1)*
     -             zb(i4,i3))))**2*
     -    (-(za(i1,i2)*za(i3,i4)*zb(i2,i1)*zb(i4,i3)) + 
     -      ((za(i1,i3)*zb(i3,i1))/2. + 
     -         (za(i2,i3)*zb(i3,i2))/2. + 
     -         (za(i1,i4)*zb(i4,i1))/2. + 
     -         (za(i2,i4)*zb(i4,i2))/2. + 
     -         Sqrt(((za(i1,i3)*zb(i3,i1))/2. + 
     -              (za(i2,i3)*zb(i3,i2))/2. + 
     -              (za(i1,i4)*zb(i4,i1))/2. + 
     -              (za(i2,i4)*zb(i4,i2))/2.)**2 - 
     -           za(i1,i2)*za(i3,i4)*zb(i2,i1)*zb(i4,i3)
     -           ))**2))

!+==== part 2 
!=====normalization 
      ggHZ_mp_3mtri=-ggHZ_mp_3mtri/(two*s(i3,i4))

      return 
      end
