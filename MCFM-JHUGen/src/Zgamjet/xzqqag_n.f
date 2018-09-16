      subroutine zajn_a60h(j1,j2,j3,j4,j5,j6,p,n,za,zb,zab,zba,a6nh)
      implicit none
      include 'types.f'
****************************************************************
*  Color ordered amplitudes for:
*  0 -> q(p1) + qb(p2) + l(p3) + lb(p4) + gam(p5) + glu(p6)
*  where line 6 is contracted with the vector n(mu)
****************************************************************
c-----the order of momentum in the argument
c-----(q,qb,l,lb,ph,g,p,n,za,zb,zab,zba,a6nh)
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      integer:: j1,j2,j3,j4,j5,j6
      real(dp):: nDp5,p(mxpart,4),n(4),isgn
      complex(dp):: zab(mxpart,mxpart),zba(mxpart,mxpart)
      complex(dp):: qcdabn(2,2,2),qcdban(2,2,2),a6nql(2,2,2)
      complex(dp):: a6nh(2,8)
      integer:: i,j,k
c-----n.p5
      nDp5=n(4)*p(j5,4)-n(3)*p(j5,3)-n(2)*p(j5,2)-n(1)*p(j5,1)
c-----check that n.p6=0
      call checkndotp(p,n,j6)
c-----call qqbzgg_n
c-----argument of qcdabn or qcdban:
c-----qcdxxn(pg,pq,lp)
c-----pg:gluon pol,pq:quark pol,lp,lepton pol
c-----1=(-);2=(+)
      call subqcdn(j2,j1,j3,j4,j5,j6,nDp5,za,zb,zab,zba,qcdabn,qcdban)
      call a6zqqbagnQL(j1,j2,j3,j4,j5,j6,za,zb,zab,zba,a6nql)
c-----helicity index
c-----h(1)=(+++)
c-----h(2)=(-++)
c-----h(3)=(++-)
c-----h(4)=(-+-)
c-----h(5)=(+--)
c-----h(6)=(---)
c-----h(7)=(+-+)
c-----h(8)=(--+)
c-----Q-contributions
      a6nh(1,1) = qcdabn(2,2,2) + qcdban(2,2,2)
      a6nh(1,2) = qcdabn(1,2,2) + qcdban(1,2,2)
      a6nh(1,3) = qcdabn(2,2,1) + qcdban(2,2,1)
      a6nh(1,4) = qcdabn(1,2,1) + qcdban(1,2,1)
      a6nh(1,5) = qcdabn(2,1,1) + qcdban(2,1,1)
      a6nh(1,6) = qcdabn(1,1,1) + qcdban(1,1,1)
      a6nh(1,7) = qcdabn(2,1,2) + qcdban(2,1,2)
      a6nh(1,8) = qcdabn(1,1,2) + qcdban(1,1,2)
c-----QL contributions
      a6nh(2,1) = two*a6nql(2,2,2)*isgn(j6)
      a6nh(2,2) = two*a6nql(1,2,2)*isgn(j6) 
      a6nh(2,3) = two*a6nql(2,2,1) 
      a6nh(2,4) = two*a6nql(1,2,1) 
      a6nh(2,5) = two*a6nql(2,1,1)*isgn(j6) 
      a6nh(2,6) = two*a6nql(1,1,1)*isgn(j6)
      a6nh(2,7) = two*a6nql(2,1,2) 
      a6nh(2,8) = two*a6nql(1,1,2) 
c-----donehere  
      return
      end



      subroutine a6zqqbagnQL(j1,j2,j3,j4,j5,j6,za,zb,zab,zba,a6nql)
      implicit none
      include 'types.f'
****************************************************************
*  Helicity amplitudes for:
*  0 -> q(p1) + qb(p2) + l(p3) + lb(p4) + gam(p5) + glu(p6)
*  where the photon is emitted from lepton line, and,
*  line 6 is contracted with the vector n(mu)
****************************************************************
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      integer:: j1,j2,j3,j4,j5,j6
      complex(dp):: zab(mxpart,mxpart),zba(mxpart,mxpart)
c-----1st argument is the photon polarization
c-----2nd argument is the quark
c-----3nd argument is the lepton line
c-----1 is left handed
c-----2 is right handed
      complex(dp):: a6nql(2,2,2)
      complex(dp)::ab41,ab46,ab23,ab63,ab31,ab36,ab24,ab64,
     & DENaa,DENbb
      real(dp):: It345,Is16,Is26
      integer:: pg,pq,pl
c-----initialize
      do pg=1,2
      do pq=1,2
      do pl=1,2
         a6nql(pg,pq,pl)=czip
      enddo
      enddo
      enddo
c-----shorthand notation
c-----kinematic invariants
      It345 = 1._dp/(s(j3,j4)+s(j3,j5)+s(j4,j5))
      Is16  = 1._dp/s(j1,j6)
      Is26  = 1._dp/s(j2,j6)
c-----spinor products
      ab41 = za(j4,j3)*zb(j3,j1)+za(j4,j5)*zb(j5,j1)
      ab46 = za(j4,j3)*zb(j3,j6)+za(j4,j5)*zb(j5,j6)
      ab23 = za(j2,j4)*zb(j4,j3)+za(j2,j5)*zb(j5,j3)
      ab63 = za(j6,j4)*zb(j4,j3)+za(j6,j5)*zb(j5,j3)
      ab31 = za(j3,j4)*zb(j4,j1)+za(j3,j5)*zb(j5,j1)
      ab36 = za(j3,j4)*zb(j4,j6)+za(j3,j5)*zb(j5,j6)
      ab24 = za(j2,j3)*zb(j3,j4)+za(j2,j5)*zb(j5,j4)
      ab64 = za(j6,j3)*zb(j3,j4)+za(j6,j5)*zb(j5,j4)
      DENaa = 1._dp/za(j3,j5)/za(j4,j5)
      DENbb = 1._dp/zb(j3,j5)/zb(j4,j5)
c-----Helicity amplitudes
c-----h(1)=(+++)=(2,2,2)
      a6nql(2,2,2) = 
     .(-Is16*za(j2,j4)*(ab41*zab(j1,j1)+ab46*zab(j6,j1))
     & +Is26*ab41*(zab(j2,j2)*za(j2,j4)+zab(j2,j6)*za(j6,j4)) 
     .)*DENaa*It345
c-----h(2)=(-++)=(1,2,2)
      a6nql(1,2,2) =
     .(+Is16*ab23*(zb(j3,j1)*zab(j1,j1)+zb(j3,j6)*zab(j6,j1))
     & -Is26*zb(j3,j1)*(zab(j2,j2)*ab23+zab(j2,j6)*ab63) 
     .)*DENbb*It345
c-----h(3)=(++-)=(2,2,1)
      a6nql(2,2,1) =
     .(+Is16*za(j2,j3)*(ab31*zab(j1,j1)+ab36*zab(j6,j1))
     & -Is26*ab31*(zab(j2,j2)*za(j2,j3)+zab(j2,j6)*za(j6,j3)) 
     .)*DENaa*It345
c-----h(4)=(-+-)=(1,2,1)
      a6nql(1,2,1) =
     .(-Is16*ab24*(zb(j4,j1)*zab(j1,j1)+zb(j4,j6)*zab(j6,j1))
     & +Is26*zb(j4,j1)*(zab(j2,j2)*ab24+zab(j2,j6)*ab64) 
     .)*DENbb*It345
c-----h(5)=(+--)=(2,1,1)
      a6nql(2,1,1)=conjg(a6nql(1,2,2))
c-----h(6)=(---)=(1,1,1)
      a6nql(1,1,1)=conjg(a6nql(2,2,2))
c-----h(7)=(+-+)=(2,1,2)
      a6nql(2,1,2)=conjg(a6nql(1,2,1))
c-----h(8)=(--+)=(1,1,2)
      a6nql(1,1,2)=conjg(a6nql(2,2,1))
c-----done
      return
      end         

      function isgn(j1)
      implicit none
      include 'types.f'
      real(dp):: isgn
      
      integer:: j1
      if (j1<=2) then
         isgn=-1._dp
      else
         isgn=1._dp
      endif
      return
      end

