      subroutine getggZZamps(p,dolight,dobottom,dotop,
     & Mloop_uptype,Mloop_dntype,Mloop_bquark,Mloop_tquark)
      implicit none
      include 'types.f'
c--- Returns a series of arrays representing the dressed amp[itudes
c--- for the process gg->ZZ; there are:
c---        Mloop_uptype(h1,h2,h34,h56)   single massless up-type quark
c---        Mloop_dntype(h1,h2,h34,h56)   single massless down-type quark
c---        Mloop_bquark(h1,h2,h34,h56)   top quark mass=mt
c---        Mloop_tquark(h1,h2,h34,h56)   bottom quark mass=mb
c---
c--- The overall factor on the amplitude is:
c---
c---      4._dp*esq*gsq/(16._dp*pisq)*esq * delta(a,b)
c---
c--- apportioned as follows:
c--- 8*e^2*gsq from gg->ZZ, extra esq from Z decays, 1/(16*pisq) from
c--- QCDLoop normalization and 1/2 from color factor
c---
c--- The series of flags control whether or not the amplitudes should
c--- be computed:
c---        dolight  -> Mloop_uptype, Mloop_dntype
c---        dobottom -> Mloop_bquark
c---        dotop    -> Mloop_tquark
c---
c--- An array that is identically zero is returned if the result is
c--- expected to be unreliable, namely pt(Z)<ptZsafetycut set below
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'zcouple.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      include 'scale.f'
      include 'docheck.f'
      include 'qlfirst.f'
      include 'first.f'
      logical:: dolight,dobottom,dotop,ggZZuse6d
      integer:: h1,h2,h34,h56,up,dn,om,nu
      real(dp):: p(mxpart,4),cvec(2),cax(2),cl1(2),cl2(2),
     & ptZsafetycut_massless,ptZsafetycut_massive,ptZ,pttwo
      real(dp):: phi,muk,rho,ssig,csig,theta,
     & p1true(4),p2true(4),p3true(4),p4true(4),p5true(4),p6true(4)
      complex(dp):: Avec(2,2,2,2),a64v,prop34,prop56,
     & AmtLL(2,2,2,2),AmtLR(2,2,2,2),AmbLL(2,2,2,2),AmbLR(2,2,2,2),
     & Mloop_uptype(2,2,2,2),Mloop_dntype(2,2,2,2),
     & Mloop_bquark(2,2,2,2),Mloop_tquark(2,2,2,2),
     & Amb_vec,Amb_ax,Amt_vec,Amt_ax
c     & AmtLL_new(2,2,2,2),AmtLR_new(2,2,2,2)
      common/ggZZuse6d/ggZZuse6d
      parameter(up=2,dn=1)
!$omp threadprivate(/ggZZuse6d/)

      ggZZuse6d=.true.    ! TRUE -> use 6d instead of 4d boxes (recommended)
c      ggZZuse6d=.false.  ! FALSE -> use 4d boxes, poorer numerical stability

c--- omit u,d,s,c loops for pt(Z) < "ptZsafetycut_massless" (for num. stability)
      ptZsafetycut_massless=0.1_dp
c--- omit t,b quark loops for pt(Z) < "ptZsafetycut_massive"  (for num. stability)
      ptZsafetycut_massive=0.1_dp

      if (first) then
        write(6,*)'****************************************************'
        write(6,*)'*                                                  *'
        if (dolight) then
        write(6,*)'*  gg->ZZ box loop includes gens. 1 and 2          *'
        else
        write(6,*)'*  gg->ZZ box loop does not include gens. 1 and 2  *'
        endif
        if (dobottom) then
        write(6,*)'*  gg->ZZ box loop includes bottom quark           *'
        else
        write(6,*)'*  gg->ZZ box loop does not include bottom quark   *'
        endif
        if (dotop) then
        write(6,*)'*  gg->ZZ box loop includes top quark              *'
        else
        write(6,*)'*  gg->ZZ box loop does not include top quark      *'
        endif
        write(6,*)'*                                                  *'
        write(6,54) ptZsafetycut_massless,'(gens. 1,2)'
        write(6,54) ptZsafetycut_massive, '(b,t loops)'
        write(6,*)'*                                                  *'
        write(6,*)'****************************************************'
        first=.false. 
        if (qlfirst) then
          qlfirst=.false.
          call qlinit
        endif
   54   format(' *  Numer. stability: pt(Z) >',f6.3,' GeV ',a11,' *')
      endif

c--- set up spinor products      
      call spinoru(6,p,za,zb)

c--- if this variable (passed via common block) is true then
c--- set the kinematics to a special point for numerical check
      if (docheck) then
        include 'kinpoint.f'
        mu=1._dp
        musq=mu**2
        mt=0.4255266775_dp
        do nu=1,4
        om=nu-1
        if (nu==1) om=4
        p(1,om)=p1true(nu)
        p(2,om)=p2true(nu)
        p(3,om)=p3true(nu)
        p(4,om)=p4true(nu)
        p(5,om)=p5true(nu)
        p(6,om)=p6true(nu)
        enddo
        do nu=1,6
        write(6,'(a4,i2,4f18.12)') 'p_',nu,
     &   p(nu,4),p(nu,1),p(nu,2),p(nu,3)
        enddo
        call spinorz(6,p,za,zb)
      endif

c--- fill amplitudes
c--- labels are: helicity of gluons, lepton 3 and lepton 5

c--- compute amplitudes with massless quarks
      if (dolight) then
        Avec(2,2,1,1)=a64v('q+qb-g-g-',3,4,1,2,6,5,zb,za)*(-im)
        Avec(2,1,1,1)=a64v('q+qb-g-g+',3,4,1,2,6,5,zb,za)*(-im)
        Avec(1,2,1,1)=a64v('q+qb-g+g-',3,4,1,2,6,5,zb,za)*(-im)
        Avec(1,1,1,1)=a64v('q+qb-g+g+',3,4,1,2,6,5,zb,za)*(-im)

        Avec(2,2,2,1)=a64v('q+qb-g+g+',3,4,1,2,5,6,za,zb)*(-im)
        Avec(2,1,2,1)=a64v('q+qb-g+g-',3,4,1,2,5,6,za,zb)*(-im)
        Avec(1,2,2,1)=a64v('q+qb-g-g+',3,4,1,2,5,6,za,zb)*(-im)
        Avec(1,1,2,1)=a64v('q+qb-g-g-',3,4,1,2,5,6,za,zb)*(-im)

        Avec(2,2,1,2)=a64v('q+qb-g-g-',3,4,1,2,5,6,zb,za)*(-im)
        Avec(2,1,1,2)=a64v('q+qb-g-g+',3,4,1,2,5,6,zb,za)*(-im)
        Avec(1,2,1,2)=a64v('q+qb-g+g-',3,4,1,2,5,6,zb,za)*(-im)
        Avec(1,1,1,2)=a64v('q+qb-g+g+',3,4,1,2,5,6,zb,za)*(-im)

        Avec(2,2,2,2)=a64v('q+qb-g+g+',3,4,1,2,6,5,za,zb)*(-im)
        Avec(2,1,2,2)=a64v('q+qb-g+g-',3,4,1,2,6,5,za,zb)*(-im)
        Avec(1,2,2,2)=a64v('q+qb-g-g+',3,4,1,2,6,5,za,zb)*(-im)
        Avec(1,1,2,2)=a64v('q+qb-g-g-',3,4,1,2,6,5,za,zb)*(-im)
c--- a64v amplitudes have overall color factor delta(A,B);
c--- put in factor of two to extract delta(A,B)/2 as in massive case
        Avec(:,:,:,:)=two*Avec(:,:,:,:)
      endif
      
c--- compute amplitudes with a massive bottom quark internal loop
      if (dobottom) then
        if (ggZZuse6d) then ! scalar integrals including 6-d box
          call ggZZmassamp_new(p,za,zb,mb,AmbLL,AmbLR)
        else                ! all scalar integrals in 4-d
          call ggZZmassamp(p,za,zb,mb,AmbLL,AmbLR)
        endif
      endif
      
c--- compute amplitudes with a massive top quark internal loop
      if (dotop) then
        if (docheck) write(6,*) '>>>>>> TOP AMPLITUDE <<<<<'
        if (ggZZuse6d) then ! scalar integrals including 6-d box
          call ggZZmassamp_new(p,za,zb,mt,AmtLL,AmtLR)
        else                ! all scalar integrals in 4-d
          call ggZZmassamp(p,za,zb,mt,AmtLL,AmtLR)
        endif
        if (docheck) then
          call ggZZwritetable
          write(6,*) 'mt=',mt
          write(6,*)
          write(6,'(a15,SP,2e20.10)') 'AmtLR(2,2,1,1)',AmtLR(2,2,1,1)
          write(6,'(a15,SP,2e20.10)') 'AmtLR(1,2,1,1)',AmtLR(1,2,1,1)
          write(6,*)
          write(6,'(a15,SP,2e20.10)') 'AmtLL(2,2,1,1)',AmtLL(2,2,1,1)
          write(6,'(a15,SP,2e20.10)') 'AmtLL(1,2,1,1)',AmtLL(1,2,1,1)
          stop
        endif
      endif      

c--- propagator factors
      prop34=s(3,4)/cplx2(s(3,4)-zmass**2,zmass*zwidth)
      prop56=s(5,6)/cplx2(s(5,6)-zmass**2,zmass*zwidth)

c--- left and right-handed lepton couplings as an array
c--- cl1 associated with Z(3+4), cl2 associated with Z(5+6)
      cl1(1)=l1 
      cl1(2)=r1
      cl2(1)=l2
      cl2(2)=r2

c--- vector and axial couplings as an array for up/down quarks
      cvec(up)=half*(L(up)+R(up))
      cvec(dn)=half*(L(dn)+R(dn))
      cax(up)=half*(L(up)-R(up))
      cax(dn)=half*(L(dn)-R(dn))

c--- dress vector and axial amplitudes with appropriate couplings
      do h1=1,2
      do h2=1,2
      do h34=1,2
      do h56=1,2

      if (dolight) then
c--- internal loops of massless up-type quarks      
        Mloop_uptype(h1,h2,h34,h56)=Avec(h1,h2,h34,h56)*(
     & (Qu*q1+cvec(up)*cl1(h34)*prop34)*(Qu*q2+cvec(up)*cl2(h56)*prop56)
     & +(cax(up)*cl1(h34)*prop34)*(cax(up)*cl2(h56)*prop56))      
c--- internal loops of massless down-type quarks      
        Mloop_dntype(h1,h2,h34,h56)=Avec(h1,h2,h34,h56)*(
     & (Qd*q1+cvec(dn)*cl1(h34)*prop34)*(Qd*q2+cvec(dn)*cl2(h56)*prop56)
     & +(cax(dn)*cl1(h34)*prop34)*(cax(dn)*cl2(h56)*prop56))
      else
        Mloop_uptype(h1,h2,h34,h56)=czip
        Mloop_dntype(h1,h2,h34,h56)=czip
      endif
      
      if (dobottom) then
c--- implementation of b-quark loop in terms of vector and axial couplings
        Amb_vec=two*(AmbLL(h1,h2,h34,h56)+AmbLR(h1,h2,h34,h56))
        Amb_ax =two*(AmbLL(h1,h2,h34,h56)-AmbLR(h1,h2,h34,h56))
        Mloop_bquark(h1,h2,h34,h56)=im*(
     &   Amb_vec*
     & (Qd*q1+cvec(dn)*cl1(h34)*prop34)*(Qd*q2+cvec(dn)*cl2(h56)*prop56)
     &+Amb_ax*(cax(dn)*cl1(h34)*prop34)*(cax(dn)*cl2(h56)*prop56))
      else
        Mloop_bquark(h1,h2,h34,h56)=czip
      endif
      
      if (dotop) then
c--- implementation of t-quark loop in terms of vector and axial couplings
      Amt_vec=two*(AmtLL(h1,h2,h34,h56)+AmtLR(h1,h2,h34,h56))
      Amt_ax =two*(AmtLL(h1,h2,h34,h56)-AmtLR(h1,h2,h34,h56))
      Mloop_tquark(h1,h2,h34,h56)=im*(
     & Amt_vec*
     & (Qu*q1+cvec(up)*cl1(h34)*prop34)*(Qu*q2+cvec(up)*cl2(h56)*prop56)
     &+Amt_ax*(cax(up)*cl1(h34)*prop34)*(cax(up)*cl2(h56)*prop56))
      else
        Mloop_tquark(h1,h2,h34,h56)=czip
      endif
      
      enddo
      enddo
      enddo
      enddo

c--- remove contributions if in unstable regions (pt(Z) small)          
      ptZ=pttwo(3,4,p)
      if (ptZ < ptZsafetycut_massless) then
        Mloop_uptype(:,:,:,:)=czip
        Mloop_dntype(:,:,:,:)=czip
        Mloop_bquark(:,:,:,:)=czip
      endif  
      if (ptZ < ptZsafetycut_massive) then
        Mloop_tquark(:,:,:,:)=czip
      endif  
            
      return
      end
      
      
