************************************************************************
*     Author: J. M. Campbell                                           *
*     August, 1999                                                     *
*     Calculates the nj-jet subtraction term corresponding to dipole   *
*     nd with momentum p and dipole kinematics (ip,jp) wrt kp          *
*     Automatically chooses dipole kind                                *
*     Returns the dipoles in sub,subv and matrix elements in msq,msqv  *
*     nd labels the dipole configurations                              *
*     ip labels the emitter parton                                     *
*     jp labels the emitted parton                                     *
*     kp labels the spectator parton                                   *
*     subr_born is the subroutine which call the born process          *
*     subr_corr is the subroutine which call the born process dotted   *
*     with vec for an emitted gluon only                               *
************************************************************************

      subroutine dips(nd,p,ip,jp,kp,sub,subv,msq,msqv,
     . subr_born,subr_corr)
      implicit none
      include 'constants.f'
      include 'qcdcouple.f'
      include 'qqgg.f'
      include 'ptilde.f'
      include 'alfacut.f'
      include 'process.f'
      include 'dynamicscale.f'
      include 'initialscales.f'
      include 'dipolescale.f'
      include 'facscale.f'
      include 'incldip.f'
      double precision p(mxpart,4),ptrans(mxpart,4),sub(4),subv,vecsq
      double precision x,omx,z,omz,y,omy,u,omu,sij,sik,sjk,dot,vec(4)
      double precision msq(-nf:nf,-nf:nf),msqv(-nf:nf,-nf:nf),vtilde
      integer nd,ip,jp,kp,nu,j,k,ipt
c--      logical includedipole
      external subr_born,subr_corr
            
C---Initialize the dipoles to zero
      do j=1,4
      sub(j)=0d0
      enddo
      subv=0d0
      call zeromsq(msq,msqv)
      incldip(nd)=.true.
      
      sij=two*dot(p,ip,jp)
      sik=two*dot(p,ip,kp)
      sjk=two*dot(p,jp,kp)

***********************************************************************
*************************** INITIAL-INITIAL ***************************
***********************************************************************
      if ((ip .le. 2) .and. (kp .le. 2)) then
        omx=-(sij+sjk)/sik
        x=one-omx
        vtilde=sij/sik

C---Modification so that only close to singular subtracted
        if (-vtilde .gt. aii) then
           incldip(nd)=.false.
           return
        endif
        
        call transform(p,ptrans,x,ip,jp,kp)
        call storeptilde(nd,ptrans)
        
c-- Check to see if this dipole will be included
c        incldip(nd)=includedipole(nd,ptrans)
C--if not return
c        if (incldip(nd) .eqv. .false.) return
        
        vecsq=-sij*sjk/sik
        do nu=1,4
          vec(nu)=(p(jp,nu)-vtilde*p(kp,nu))/dsqrt(-vecsq)
        enddo

c--- if using a dynamic scale, set that scale with dipole kinematics      
      if (dynamicscale) then
        call scaleset(initscale,initfacscale,ptrans)
        dipscale(nd)=facscale
      endif

c--- for "gamgam" process, gvec contribution is not written in
c--- a canonical way; instead, pass emitted vector directly as "vec"
        if (case .eq. 'gamgam') then
          do nu=1,4
            vec(nu)=p(jp,nu)
          enddo
      endif
      
        call subr_born(ptrans,msq)
        call subr_corr(ptrans,vec,ip,msqv)

        sub(qq)=-gsq/x/sij*(two/omx-one-x)
        sub(gq)=-gsq/sij
        sub(qg)=-gsq/x/sij*(one-two*x*omx)
        sub(gg)=-2d0*gsq/x/sij*(x/omx+x*omx)
        subv   =-4d0*gsq/x/sij*omx/x

***********************************************************************
*************************** INITIAL-FINAL *****************************
***********************************************************************
      elseif ((ip .le. 2) .and. (kp .gt. 2)) then
        u=sij/(sij+sik)

        omx=-sjk/(sij+sik)
        x=one-omx
        omu=sik/(sij+sik)
C---npart is the number of particles in the final state
C---transform the momenta so that only the first npart+1 are filled
        call transform(p,ptrans,x,ip,jp,kp)
        call storeptilde(nd,ptrans)

c-- Check to see if this dipole will be included
c        incldip(nd)=includedipole(nd,ptrans)
C-- if not return
c        if (incldip(nd) .eqv. .false.) return

c--- if using a dynamic scale, set that scale with dipole kinematics      
      if (dynamicscale) then
        call scaleset(initscale,initfacscale,ptrans)
        dipscale(nd)=facscale
      endif
      
c--- Calculate the matrix element now because it might be needed
c--- in the final-initial segment, regardless of whether or not the
c--- alfa cut fails here
        call subr_born(ptrans,msq)
C---Modification so that only close to singular subtracted
C---Do not set incldip because initial-final can fail 
C---but final initial needs still to be tested

        if ((case .eq. 'H_1jet') .and. (ip.eq.1) .and. (jp.eq.6)) then
c--- do nothing
        else
          if (u .gt. aif) return
        endif
        
        do nu=1,4
           vec(nu)=(p(jp,nu)/u-p(kp,nu)/omu)/dsqrt(sjk)
        enddo
        
        call subr_corr(ptrans,vec,ip,msqv)        
        sub(qq)=-gsq/x/sij*(two/(omx+u)-one-x)
        sub(gq)=-gsq/sij
        sub(qg)=-gsq/x/sij*(one-two*x*omx)
        sub(gg)=-2d0*gsq/x/sij*(one/(omx+u)-one+x*omx)
        subv   =-4d0*gsq/x/sij*(omx/x*u*(one-u))

***********************************************************************
*************************** FINAL-INITIAL *****************************
***********************************************************************
      elseif ((ip .gt. 2) .and. (kp .le. 2)) then
c-- Check to see if this dipole will be included - should have been
c-- already determined at this point in the initial-final phase
c        if (incldip(nd) .eqv. .false.) return
        
c--- note, here we assume that msq kinematics are already taken care of
c--- for msq, although msqv must be recalculated each time
        omx=-sij/(sjk+sik)
C---Modification so that only close to singular subtracted
        if (omx .gt. afi) return
        
        x=one-omx
        z=sik/(sik+sjk)
        omz=sjk/(sik+sjk)
        do nu=1,4
          vec(nu)=z*p(ip,nu)-omz*p(jp,nu)
        enddo
C---call msqv again because vec has changed
        do j=1,mxpart
        do k=1,4
          ptrans(j,k)=ptilde(nd,j,k)
        enddo
        enddo

c--- if using a dynamic scale, set that scale with dipole kinematics      
      if (dynamicscale) then
        call scaleset(initscale,initfacscale,ptrans)
        dipscale(nd)=facscale
      endif      

c--- do something special if we're doing W+2,Z+2jet (jp .ne. 7)
        if ((jp .ne. 7) .and. (case .ne. 'HWWjet')
     .      .and. (case .ne. 'HZZjet')
     .      .and. (case .ne. 'qq_HWW')
     .      .and. (case .ne. 'qq_HZZ')) then
          if (ip .lt. 7) then
C ie for cases 56_i,65_i
            ipt=5
          else
C ie for cases 76_i,75_i
            ipt=6
          endif
        else
C ie for cases 57_i,67_i
          ipt=ip
        endif
      
c--- do something special for HWW/HZZ+2 jets
        if ((case .eq. 'HWW2jt') .or. (case .eq. 'HZZ2jt')) then
          if (jp .ne. 9) then
            if (ip .lt. 9) then
C ie for cases 78_i,87_i
              ipt=7
            else
C ie for cases 98_i,97_i
              ipt=8
            endif
          else
C ie for cases 79_i,89_i
            ipt=ip
          endif
      endif      
        
c--- do something special for direct photon production
        if (case .eq. 'dirgam') ipt=4
      
        call subr_corr(ptrans,vec,ipt,msqv)
                                
        sub(qq)=+gsq/x/sij*(two/(omz+omx)-one-z)
        sub(gq)=+gsq/x/sij
        sub(gg)=+2d0*gsq/x/sij*(one/(omz+omx)+one/(z+omx)-two) 
        subv   =+4d0*gsq/x/sij/sij

***********************************************************************
**************************** FINAL-FINAL ******************************
***********************************************************************
      elseif ((ip .gt. 2) .and. (kp .gt. 2)) then
c------Eq-(5.2)    
        y=sij/(sij+sjk+sik)

C---Modification so that only close to singular subtracted
        if (y .gt. aff) then
          incldip(nd)=.false.
          return
        endif

        z=sik/(sjk+sik)
        omz=one-z
        omy=one-y
C---calculate the ptrans-momenta 

        call transform(p,ptrans,y,ip,jp,kp)
        call storeptilde(nd,ptrans)

c-- Check to see if this dipole will be included
c        incldip(nd)=includedipole(nd,ptrans)
c        if (incldip(nd) .eqv. .false.) return
        
        do nu=1,4
          vec(nu)=z*p(ip,nu)-omz*p(jp,nu)
        enddo

c--- if using a dynamic scale, set that scale with dipole kinematics      
        if (dynamicscale) then
         call scaleset(initscale,initfacscale,ptrans)
         dipscale(nd)=facscale
        endif
      
        call subr_born(ptrans,msq)
        if (case .eq. 'epem3j') then
          ipt=5
        else
          if (ip .lt. kp) then
            ipt=5
          else
            ipt=6
          endif
        endif
               
c--- do something special for HWW/HZZ+2 jets
        if ((case .eq. 'HWW2jt') .or. (case .eq. 'HZZ2jt')) then
          if (ip .lt. kp) then
            ipt=7
          else
            ipt=8
          endif
      endif      

        call subr_corr(ptrans,vec,ipt,msqv)
         
        sub(qq)=gsq/sij*(two/(one-z*omy)-one-z)
        sub(gq)=gsq/sij
        sub(gg)=gsq/sij*(two/(one-z*omy)+two/(one-omz*omy)-four)
        subv   =+4d0*gsq/sij/sij

      endif
      
      return
      end
      
