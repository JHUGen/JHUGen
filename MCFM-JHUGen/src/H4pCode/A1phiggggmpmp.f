c--- Results taken from:
c---   E.~W.~N.~Glover, P.~Mastrolia and C.~Williams,
c---   %``One-loop phigggg-MHV amplitudes using the unitarity bootstrap: the general
c---   %helicity case,''
c---   JHEP {\bf 0808}, 017 (2008)
c---   [arXiv:0804.4149 [hep-ph]].

c--- with the sign of rational terms reversed

c---  Note that a factor of c_\Gamma is missing in Eq. (5.1)

      function A1phiggggmpmp(j1,j2,j3,j4,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: A1phiggggmpmp

C     arXiv:0804.4149v3, Eq.(5.15)
      integer:: j1,j2,j3,j4
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      complex(dp):: C4mpmp,Rhat4mpmp,CR4mpmp
      A1phiggggmpmp=C4mpmp(j1,j2,j3,j4,za,zb)
     &            +CR4mpmp(j1,j2,j3,j4,za,zb)
     &          +Rhat4mpmp(j1,j2,j3,j4,za,zb)

      return
      end

      function C4mpmp(j1,j2,j3,j4,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: C4mpmp

C     arXiv:0804.4149v3, Eq.(5.1) (factor of C_\Gamma removed)
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'sprods_com.f'
      include 'scprods_com.f'
      include 'zprods_decl.f'
      integer:: j1,j2,j3,j4,i
      complex(dp):: C4mpmpsub,A0phiggggmpmp,F31m,F42me,F41m,sum
      integer,parameter::ii(7)=(/1,2,3,4,1,2,3/)

c--- set up 's-comma' products
      sc(1,1)=zip
      sc(1,2)=s(j1,j2)
      sc(1,3)=s(j1,j2)+s(j1,j3)+s(j2,j3)
      sc(1,4)=s(j1,j2)+s(j1,j3)+s(j1,j4)+s(j2,j3)+s(j2,j4)+s(j3,j4)
      sc(2,1)=sc(1,4) !s(j1,j2)+s(j1,j3)+s(j1,j4)+s(j2,j3)+s(j2,j4)+s(j3,j4)
      sc(2,2)=zip
      sc(2,3)=s(j2,j3)
      sc(2,4)=s(j2,j3)+s(j2,j4)+s(j3,j4)
      sc(3,1)=s(j3,j4)+s(j3,j1)+s(j4,j1)
      sc(3,2)=sc(1,4) !s(j1,j2)+s(j1,j3)+s(j1,j4)+s(j2,j3)+s(j2,j4)+s(j3,j4)
      sc(3,3)=zip
      sc(3,4)=s(j3,j4)
      sc(4,1)=s(j4,j1)
      sc(4,2)=s(j4,j1)+s(j4,j2)+s(j1,j2)
      sc(4,3)=sc(1,4) !s(j1,j2)+s(j1,j3)+s(j1,j4)+s(j2,j3)+s(j2,j4)+s(j3,j4)
      sc(4,4)=zip

      sum=czip
      do i=1,4
      sum=sum
     & -0.5_dp*F42me(sc(ii(i),ii(i+3)),sc(ii(i+1),ii(i+2))
     &             ,sc(ii(i+1),ii(i+3)),sc(ii(i),ii(i+2)))
     & -0.5_dp*F41m(sc(ii(i),ii(i+2)),sc(ii(i),ii(i+1)),
     &                               sc(ii(i+1),ii(i+2)))
     & +F31m(sc(ii(i),ii(i+2)))-F31m(sc(ii(i),ii(i+3)))
      enddo
      sum=sum
     & +C4mpmpsub(j1,j2,j3,j4,za,zb)
     & +C4mpmpsub(j1,j4,j3,j2,za,zb)
     & +C4mpmpsub(j3,j2,j1,j4,za,zb)
     & +C4mpmpsub(j3,j4,j1,j2,za,zb)
      C4mpmp=sum*A0phiggggmpmp(j1,j2,j3,j4,za,zb)
      return
      end

      function C4mpmpsub(j1,j2,j3,j4,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: C4mpmpsub

C     arXiv:0804.4149v3, Eq.(5.1) (factor of C_\Gamma removed)
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      include 'nflav.f'
      integer:: j1,j2,j3,j4
      complex(dp):: trm,trm3241,trm3421,BGRL3,BGRL2,BGRL1,F41mF
      real(dp):: s234
C     c.f. arXiv:0804.4149v3 Eq.(3.38)
      trm(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j3)*za(j3,j4)*zb(j4,j1)
      s234=s(j2,j3)+s(j2,j4)+s(j3,j4)
      trm3241=trm(j3,j2,j4,j1)
      trm3421=trm(j3,j4,j2,j1)
c--- MODIFIED: added an overall factor of (-1) here
      C4mpmpsub=+4._dp*(1._dp-real(nflav,dp)/(4._dp*xn))
     & *(trm3241*trm3421/(2._dp*(s(j2,j4)*s(j1,j3))**2)
     & *F41mF(s234,s(j2,j3),s(j3,j4))
     & -trm3241*trm3421/(s(j2,j4)*s(j1,j3)**2)
     & *BGRL1(s(j2,j3),s234))
     & +2._dp*(1._dp-real(nflav,dp)/xn)
     & *(-0.5_dp*(trm3241*trm3421/(s(j2,j4)*s(j1,j3))**2)**2
     & *F41mF(s234,s(j2,j3),s(j3,j4))
     & -trm3241*trm3421/s(j1,j3)**4
     & *(+trm3241**2/(3._dp*s(j2,j4))*BGRL3(s(j2,j3),s234)
     & +trm3421*trm3241/(2._dp*s(j2,j4)**2)*BGRL2(s(j2,j3),s234)
     & -trm3421*trm3241/(s(j2,j4)**3)*BGRL1(s(j2,j3),s234)))

      return
      end



      function Rhat4mpmp(j1,j2,j3,j4,za,zb)

C     arXiv:0804.4149v3, Eq.(5.17)
c---                    (factor of C_\Gamma = 1/(16 pi^2) removed)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'nflav.f'
      complex(dp)::Rhat4mpmp
      integer:: j1,j2,j3,j4
      complex(dp):: A0phiggggmpmp
      real(dp):: Np,s3
      s3(j1,j2,j3)=s(j1,j2)+s(j2,j3)+s(j3,j1)
      Np=2._dp*(1._dp-real(nflav,dp)/xn)

c--- Note: implicit use of A0(A,...)=-i*[A0phi(...)-A0phidagger(...)]
c---       and the fact that A0phidagger=(c.c. of A0phi)
      Rhat4mpmp=
     & -2._dp*(-im)*(A0phiggggmpmp(j1,j2,j3,j4,za,zb)
     &            -A0phiggggmpmp(j2,j3,j4,j1,zb,za))
     & +Np/12._dp*zb(j2,j4)**4/(zb(j1,j2)*zb(j2,j3)*zb(j3,j4)*zb(j4,j1))
     & *(-s(j2,j3)*s(j3,j4)/(s(j2,j4)*s3(j4,j1,j2))
     &   +3._dp*s(j2,j3)*s(j3,j4)/s(j2,j4)**2
     &   -s(j1,j2)*s(j4,j1)/(s(j2,j4)*s3(j2,j3,j4))
     &   +3._dp*s(j1,j2)*s(j4,j1)/s(j2,j4)**2)
      return
      end


      function CR4mpmp(j1,j2,j3,j4,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: CR4mpmp

C     arXiv:0804.4149v3, Eq.(5.2)
c---                    (factor of C_\Gamma = 1/(16 pi^2) removed)
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      integer:: j1,j2,j3,j4
      complex(dp):: CR4mpmpsub
      CR4mpmp=
     & +CR4mpmpsub(j1,j2,j3,j4,za,zb)+CR4mpmpsub(j1,j4,j3,j2,za,zb)
     & +CR4mpmpsub(j3,j2,j1,j4,za,zb)+CR4mpmpsub(j3,j4,j1,j2,za,zb)
      return
      end


      function CR4mpmpsub(j1,j2,j3,j4,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: CR4mpmpsub

C     arXiv:0804.4149v3, Eq.(5.2)
c---                    (factor of C_\Gamma = 1/(16 pi^2) removed)
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'nflav.f'
      integer:: j1,j2,j3,j4
      real(dp):: s3,Np
      complex(dp):: zab3
      zab3(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j3)*za(j3,j4)
      s3(j1,j2,j3)=s(j1,j2)+s(j2,j3)+s(j3,j1)
      Np=2._dp*(1._dp-real(nflav,dp)/xn)
      CR4mpmpsub=Np/(2._dp*za(j1,j2)*za(j2,j3)*za(j3,j4)*za(j4,j1))
     & *(-zab3(j3,j2,j4,j1)**3*za(j3,j4)*za(j2,j1)
     & /(3._dp*za(j4,j2)*(s3(j2,j3,j4)-s(j2,j3))**2)
     & -(zab3(j3,j2,j4,j1)*za(j3,j4)*za(j2,j1)/za(j4,j2))**2
     & /(2._dp*(s3(j2,j3,j4)-s(j2,j3))))
     & *(1._dp/s(j2,j3)+1._dp/s3(j2,j3,j4))
c      write(6,*) 'CR4mpmpsub',CR4mpmpsub
      return
      end
