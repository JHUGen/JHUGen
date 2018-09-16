      subroutine ubtzdamp_dk(p_dk,j1,j2,j3,j4,j8,amp_dk)
      implicit none
      include 'types.f'
C     Matrix element - not squared - to allow inclusion
C     of decay matrix elements
c     u(-p1)+b(p2)->Z(p3,p4)+t(-> nu(p5) e+(p6) b(p7))+d(p6)
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'zprods_com.f'
      real(dp):: p(mxpart,4),p_dk(mxpart,4),q(mxpart,4),dot
      real(dp):: sud,sbt,stz,suz,sbz,szd,sz,suda,ptDpep
      complex(dp):: amp_dk(2,2)
      integer:: nu,j1,j2,j3,j4,j8,p1,p2,ph1,ph2,k3,e3,p5,eta
      parameter(p1=1,p2=2,ph1=3,ph2=4,k3=5,e3=6,p5=7)
      eta=6
      p(1,:)=p_dk(j1,:)
      p(2,:)=p_dk(j2,:)
      p(3,:)=p_dk(j3,:)
      p(4,:)=p_dk(j4,:)
      p(5,:)=p_dk(5,:)+p_dk(6,:)+p_dk(7,:)
      p(6,:)=p_dk(j8,:)
      ptDpep=p(5,4)*p_dk(eta,4)-p(5,1)*p_dk(eta,1)-p(5,2)*p_dk(eta,2)-
     & p(5,3)*p_dk(eta,3)
      do nu=1,4
      q(p1,nu)=p(1,nu)
      q(p2,nu)=p(2,nu)
      q(ph1,nu)=p(3,nu)
      q(ph2,nu)=p(4,nu)
      q(k3,nu)=p(5,nu)-mt**2*p_dk(eta,nu)/(2d0*ptDpep)
      q(e3,nu)=mt**2*p_dk(eta,nu)/(2d0*ptDpep)
      q(p5,nu)=p(6,nu)
      enddo
      sud=2d0*dot(p,1,6)
      sbt=2d0*dot(p,2,5)+mt**2
      stz=2d0*dot(p,1,2)+2d0*dot(p,2,6)+2d0*dot(p,1,6)
      suz=2d0*dot(p,1,3)+2d0*dot(p,1,4)+2d0*dot(p,3,4)
      sbz=2d0*dot(p,2,3)+2d0*dot(p,2,4)+2d0*dot(p,3,4)
      szd=2d0*dot(p,6,3)+2d0*dot(p,6,4)+2d0*dot(p,3,4)
      sz=2d0*dot(p,3,4)
      suda=2d0*dot(p,1,6)+2d0*dot(p,1,4)+2d0*dot(p,6,4)

c--- top production (j3=3,j4=4) as-is
      if (j3==3) then
      call spinoru(7,q,za,zb)
      else
c--- antitop production (j3=4,j4=3) requires c.c.
      call spinoru(7,q,zb,za)
      endif
c      call writeout(q)
c      pause
      call calc_ubtzdamp_dk(sud,sbt,stz,suz,sbz,szd,sz,suda,
     &  za,zb,amp_dk)
      return
      end
      
      subroutine calc_ubtzdamp_dk(sud,sbt,stz,suz,sbz,szd,
     & sz,suda,za,zb,amp)
      implicit none
      include 'types.f'
c--- First label of amp is for the lepton helicity, second is top spin
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'ewcouple.f'
      include 'zcouple.f'
      include 'zprods_decl.f'
      integer:: j1,j2,p1,p2,k3,e3,p5,pl,pa,j
      parameter(p1=1,p2=2,k3=5,e3=6,p5=7)
C      integer:: ph1,ph2
C      parameter(ph1=3,ph2=4)
      real(dp):: sud,sbt,stz,suz,sbz,szd,sz,suda,mw
      complex(dp):: facuLl,facuRl,facdLl,nrp,nrm,cprop
      complex(dp):: prW,prt,prq,iza,izb,amp(2,2),iprZ
      prW(sz)=cone/cplx2(sz-wmass**2,zip)
      prt(sz)=cone/cplx2(sz-mt**2,zip)
      prq(sz)=cone/cplx1(sz)
      iza(j1,j2)=cone/za(j1,j2)
      izb(j1,j2)=cone/zb(j1,j2)
C   
      mw=wmass

c--- Implementation of Baur-Zeppenfeld treatment of Z width
	     cprop=cplx1(1d0/sqrt((sz-zmass**2)**2+(zmass*zwidth)**2))
	     iprZ=cplx1(sz-zmass**2)

      cprop=cprop/cplx2(zip,mt*twidth)

      do j=1,2
      if (j == 1) then
        pl=3
        pa=4
	    facuLl=cplx1(Qu*q1)*iprZ+cplx1(L(2)*le)*sz
	    facuRl=cplx1(Qu*q1)*iprZ+cplx1(R(2)*le)*sz
	    facdLl=cplx1(Qd*q1)*iprZ+cplx1(L(1)*le)*sz
      else
        pl=4
        pa=3
	    facuLl=cplx1(Qu*q1)*iprZ+cplx1(L(2)*re)*sz
	    facuRl=cplx1(Qu*q1)*iprZ+cplx1(R(2)*re)*sz
	    facdLl=cplx1(Qd*q1)*iprZ+cplx1(L(1)*re)*sz
      endif

      amp(j,1)= + prW(sud)*prW(sbt)*sz**(-1)*facdLl * ( za(p1,k3)*za(p5
     &    ,pl)*zb(p1,p2)*zb(p1,pa) - za(p1,pl)*za(p5,k3)*zb(p1,p2)*zb(
     &    p1,pa) + za(p2,p5)*za(k3,pl)*zb(p1,p2)*zb(p2,pa) - za(p5,e3)*
     &    za(k3,pl)*zb(p1,e3)*zb(p2,pa) - za(p5,k3)*za(p5,pl)*zb(p1,p2)
     &    *zb(p5,pa) - za(p5,k3)*za(p5,pl)*zb(p1,pa)*zb(p2,p5) - za(p5,
     &    k3)*za(k3,pl)*zb(p1,k3)*zb(p2,pa) - za(p5,pl)*za(e3,k3)*zb(p1
     &    ,pa)*zb(p2,e3) )
      amp(j,1) = amp(j,1) + prW(sud)*prW(sbt)*sz**(-1)*facuLl * (  - 
     &    za(p1,k3)*za(p5,pl)*zb(p1,p2)*zb(p1,pa) + za(p1,pl)*za(p5,k3)
     &    *zb(p1,p2)*zb(p1,pa) - za(p2,p5)*za(k3,pl)*zb(p1,p2)*zb(p2,pa
     &    ) + za(p5,e3)*za(k3,pl)*zb(p1,e3)*zb(p2,pa) + za(p5,k3)*za(p5
     &    ,pl)*zb(p1,p2)*zb(p5,pa) + za(p5,k3)*za(p5,pl)*zb(p1,pa)*zb(
     &    p2,p5) + za(p5,k3)*za(k3,pl)*zb(p1,k3)*zb(p2,pa) + za(p5,pl)*
     &    za(e3,k3)*zb(p1,pa)*zb(p2,e3) )
      amp(j,1) = amp(j,1) + prW(sud)*prW(sbt)*izb(e3,k3)*mt**2*mw**(-2)
     & *facdLl * (  - 1.D0/2.D0*za(p5,pl)*zb(p1,pa)*zb(p2,e3) )
      amp(j,1) = amp(j,1) + prW(sud)*prW(sbt)*izb(e3,k3)*mt**2*mw**(-2)
     & *facuLl * ( 1.D0/2.D0*za(p5,pl)*zb(p1,pa)*zb(p2,e3) )
      amp(j,1) = amp(j,1) + prW(sud)*prt(stz)*sz**(-1)*facuLl * (  - 
     &    za(p1,p5)*za(k3,pl)*zb(p1,p2)*zb(p1,pa) - za(p2,p5)*za(k3,pl)
     &    *zb(p1,p2)*zb(p2,pa) )
      amp(j,1) = amp(j,1) + prW(sud)*prt(stz)*izb(e3,k3)*mt**2*sz**(-1)
     & *facuRl * ( za(p5,pl)*zb(p1,p2)*zb(e3,pa) )
      amp(j,1) = amp(j,1) + prW(sud)*prq(sbz)*sz**(-1)*facdLl * ( za(p2
     &    ,pl)*za(p5,k3)*zb(p1,p2)*zb(p2,pa) - za(p5,k3)*za(pl,pa)*zb(
     &    p1,pa)*zb(p2,pa) )
      amp(j,1) = amp(j,1) + prW(sbt)*prq(suz)*sz**(-1)*facuLl * ( za(p5
     &    ,k3)*za(p5,pl)*zb(p1,pa)*zb(p2,p5) + za(p5,k3)*za(e3,pl)*zb(
     &    p1,pa)*zb(p2,e3) + za(p5,k3)*za(k3,pl)*zb(p1,pa)*zb(p2,k3) )
      amp(j,1) = amp(j,1) + prW(sbt)*prq(szd)*sz**(-1)*facdLl * ( za(p1
     &    ,k3)*za(p5,pl)*zb(p1,p2)*zb(p1,pa) + za(p2,k3)*za(p5,pl)*zb(
     &    p1,p2)*zb(p2,pa) + za(p5,pl)*za(e3,k3)*zb(p1,p2)*zb(e3,pa) )

      amp(j,2)= + prW(sud)*prW(sbt)*mt*mw**(-2)*facdLl * (  - 1.D0/2.D0
     &    *za(p5,pl)*zb(p1,pa)*zb(p2,k3) )
      amp(j,2) = amp(j,2) + prW(sud)*prW(sbt)*mt*mw**(-2)*facuLl * ( 1.D
     &    0/2.D0*za(p5,pl)*zb(p1,pa)*zb(p2,k3) )
      amp(j,2) = amp(j,2) + prW(sud)*prW(sbt)*mt*sz**(-1)*facdLl * ( 
     &    za(p5,pl)*zb(p1,pa)*zb(p2,k3) )
      amp(j,2) = amp(j,2) + prW(sud)*prW(sbt)*mt*sz**(-1)*facuLl * ( 
     &     - za(p5,pl)*zb(p1,pa)*zb(p2,k3) )
      amp(j,2) = amp(j,2) + prW(sud)*prW(sbt)*iza(e3,k3)*mt*sz**(-1)*
     & facdLl * ( za(p1,e3)*za(p5,pl)*zb(p1,p2)*zb(p1,pa) - za(p1,pl)*
     &    za(p5,e3)*zb(p1,p2)*zb(p1,pa) + za(p2,p5)*za(e3,pl)*zb(p1,p2)
     &    *zb(p2,pa) - za(p5,e3)*za(p5,pl)*zb(p1,p2)*zb(p5,pa) - za(p5,
     &    e3)*za(p5,pl)*zb(p1,pa)*zb(p2,p5) - za(p5,e3)*za(e3,pl)*zb(p1
     &    ,e3)*zb(p2,pa) - za(p5,k3)*za(e3,pl)*zb(p1,k3)*zb(p2,pa) )
      amp(j,2) = amp(j,2) + prW(sud)*prW(sbt)*iza(e3,k3)*mt*sz**(-1)*
     & facuLl * (  - za(p1,e3)*za(p5,pl)*zb(p1,p2)*zb(p1,pa) + za(p1,pl
     &    )*za(p5,e3)*zb(p1,p2)*zb(p1,pa) - za(p2,p5)*za(e3,pl)*zb(p1,
     &    p2)*zb(p2,pa) + za(p5,e3)*za(p5,pl)*zb(p1,p2)*zb(p5,pa) + za(
     &    p5,e3)*za(p5,pl)*zb(p1,pa)*zb(p2,p5) + za(p5,e3)*za(e3,pl)*
     &    zb(p1,e3)*zb(p2,pa) + za(p5,k3)*za(e3,pl)*zb(p1,k3)*zb(p2,pa)
     &     )
      amp(j,2) = amp(j,2) + prW(sud)*prt(stz)*mt*sz**(-1)*facuRl * ( 
     &    za(p5,pl)*zb(p1,p2)*zb(k3,pa) )
      amp(j,2) = amp(j,2) + prW(sud)*prt(stz)*iza(e3,k3)*mt*sz**(-1)*
     & facuLl * (  - za(p1,p5)*za(e3,pl)*zb(p1,p2)*zb(p1,pa) - za(p2,p5
     &    )*za(e3,pl)*zb(p1,p2)*zb(p2,pa) )
      amp(j,2) = amp(j,2) + prW(sud)*iza(e3,k3)*prq(sbz)*mt*sz**(-1)*
     & facdLl * ( za(p2,pl)*za(p5,e3)*zb(p1,p2)*zb(p2,pa) - za(p5,e3)*
     &    za(pl,pa)*zb(p1,pa)*zb(p2,pa) )
      amp(j,2) = amp(j,2) + prW(sbt)*iza(e3,k3)*prq(suz)*mt*sz**(-1)*
     & facuLl * ( za(p5,e3)*za(p5,pl)*zb(p1,pa)*zb(p2,p5) + za(p5,e3)*
     &    za(e3,pl)*zb(p1,pa)*zb(p2,e3) + za(p5,e3)*za(k3,pl)*zb(p1,pa)
     &    *zb(p2,k3) )
      amp(j,2) = amp(j,2) + prW(sbt)*iza(e3,k3)*prq(szd)*mt*sz**(-1)*
     & facdLl * ( za(p1,e3)*za(p5,pl)*zb(p1,p2)*zb(p1,pa) + za(p2,e3)*
     &    za(p5,pl)*zb(p1,p2)*zb(p2,pa) )
      amp(j,2) = amp(j,2) + prW(sbt)*prq(szd)*mt*sz**(-1)*facdLl * ( 
     &     - za(p5,pl)*zb(p1,p2)*zb(k3,pa) )


c--- Non-resonant diagrams for left-handed lepton helicity only
      if (j. eq. 1) then
      nrm= + prW(sud)*prW(sbt)*suda**(-1)*q1*xw**(-1) * (  - 1.D0/2.D0*
     &    za(p1,p5)*za(k3,pl)*zb(p1,p2)*zb(p1,pa)*iprZ - 1.D0/2.D0*za(
     &    p5,pa)*za(k3,pl)*zb(p1,pa)*zb(p2,pa)*iprZ )
      nrm = nrm + prW(sud)*prW(sbt)*izb(e3,k3)*mt**2*mw**(-2)*q1*
     & xw**(-1) * (  - 1.D0/4.D0*za(p5,pl)*zb(p1,pa)*zb(p2,e3)*iprZ )

      nrp= + prW(sud)*prW(sbt)*mt*mw**(-2)*q1*xw**(-1) * (  - 1.D0/4.D0
     &    *za(p5,pl)*zb(p1,pa)*zb(p2,k3)*iprZ )
      nrp = nrp + prW(sud)*prW(sbt)*iza(e3,k3)*mt*suda**(-1)*q1*
     & xw**(-1) * (  - 1.D0/2.D0*za(p1,p5)*za(e3,pl)*zb(p1,p2)*zb(p1,pa
     &    )*iprZ - 1.D0/2.D0*za(p5,pa)*za(e3,pl)*zb(p1,pa)*zb(p2,pa)*
     &    iprZ )

      amp(j,1)=amp(j,1)+nrm
      amp(j,2)=amp(j,2)+nrp
      endif

      amp(j,1)=amp(j,1)*cprop
      amp(j,2)=amp(j,2)*cprop

      enddo


      return
      end
