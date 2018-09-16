************************************************************************
*     Author: J. M. Campbell                                           *
*     August, 2001                                                     *
*                                                                      *
*     Replica of dipolesub.f, except for the fact that extra matrix    *
*     element arrays are called in the Born term                       *
*                                                                      *
*     Differs from dipolesubx.f by the fact that msqx and mvxg         *
*     arrays with 5 and 4 dimensions have been replaced by smaller     *
*     2- and 1-dimensional arrays                                       *
*                                                                      *
*     Calculates the nj-jet subtraction term corresponding to dipole   *
*     nd with momentum p and dipole kinematics (ip,jp) wrt kp          *
*     Automatically chooses dipole kind                                *
*     Returns the dipoles in sub,subv and matrix elements in msq,msqv  *
*     nd labels the dipole configurations                              *
*     ip labels the emitter parton                                     *
*     jp labels the emitted parton                                     *
*     kp labels the spectator parton                                   *
*     msq - lowest order matrix elements at rescaled momentum, msq(j,k)*
*     msqv -  lowest order matrix elements at rescaled momentum        *
*      with emitter contracted with appropriate vector, msqv(j,k)      *
*     subr_born is the subroutine which call the born process          *
*     subr_corr is the subroutine which call the born process dotted   *
*      with vec for an emitted gluon only                              *
*     mqq - 4-quark contribution to lowest order matrix elements sqd.  *
*     msqx - lowest order matrix elements with 4 indices, msqx(j,k,l,m)*
*            Sum_{l,m} msqx(j,k,l,m) = msq(j,k)                        *
*     mg - 2-quark contribution to lowest order matrix elements sqd,   *
*           separated by colours                                       *
*     mvg - 2-quark contribution to lowest order matrix elements sqd,  *
*           separated by colours, contracted with appropriate vector   *
*     mvxg - lowest order matrix elements with 4 indices and           *
*        contracted with appropriate vector, msqvx(j,k,l,m)            *
*        Sum_{l,m} msqvx(j,k,l,m) = msqv(j,k)                          *
************************************************************************

      subroutine dipsx_new(nd,p,ip,jp,kp,sub,subv,msq,msqv,
     & subr_born,subr_corr,mqq,msqx,mg,mvg,mvxg)
      implicit none
      include 'types.f'
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'qcdcouple.f'
      include 'qqgg.f'
      include 'ptilde.f'
      include 'alfacut.f'
      include 'dynamicscale.f'
      include 'initialscales.f'
      include 'dipolescale.f'
      include 'facscale.f'
      include 'incldip.f'
      include 'ppmax.f'
      real(dp):: p(mxpart,4),ptrans(mxpart,4),sub(4),subv,vecsq
      real(dp):: x,omx,z,omz,y,omy,u,omu,sij,sik,sjk,dot,vec(4)
      real(dp):: msq(-nf:nf,-nf:nf),msqv(-nf:nf,-nf:nf),vtilde
      real(dp):: mqq(0:2,-nf:nf,-nf:nf)
      real(dp):: msqx(0:2,ppmax)
      real(dp):: mg(0:2,-nf:nf,-nf:nf)
      real(dp):: mvg(0:2,-nf:nf,-nf:nf)
      real(dp):: mvxg(ppmax)
      integer:: nd,ip,jp,kp,nu,j,k
c--      logical:: includedipole
      external subr_born,subr_corr
      
C---Initialize the dipoles to zero
      do j=1,4
      sub(j)=0._dp
      enddo
      subv=0._dp
      call zeromsq(msq,msqv)
      mqq=0._dp
      msqx=0._dp
      mg=0._dp
      mvg=0._dp
      mvxg=0._dp
      incldip(nd)=.true.

      sij=two*dot(p,ip,jp)
      sik=two*dot(p,ip,kp)
      sjk=two*dot(p,jp,kp)

***********************************************************************
*************************** INITIAL-INITIAL ***************************
***********************************************************************
      if ((ip <= 2) .and. (kp <= 2)) then
        omx=-(sij+sjk)/sik
        x=one-omx
        
        vtilde=sij/sik

C---Modification so that only close to singular subtracted
        if (-vtilde > aii) then
           incldip(nd)=.false.
           return
        endif
        
        call transform(p,ptrans,x,ip,jp,kp)
        call storeptilde(nd,ptrans)
      
        vecsq=-sij*sjk/sik
        do nu=1,4
          vec(nu)=p(jp,nu)-vtilde*p(kp,nu)
        enddo

c--- if using a dynamic scale, set that scale with dipole kinematics      
      if (dynamicscale) then
        call scaleset(initscale,initfacscale,ptrans)
        dipscale(nd)=facscale
      endif
      
        call subr_born(ptrans,msq,mqq,msqx,mg)
        call subr_corr(ptrans,vec,ip,msqv,mvg,mvxg)

        sub(qq)=-gsq/x/sij*(two/omx-one-x)
        sub(gq)=-gsq/sij
        sub(qg)=-gsq/x/sij*(one-two*x*omx)
        sub(gg)=-2._dp*gsq/x/sij*(x/omx+x*omx)
        subv   =+4._dp*gsq/x/sij*omx/x/vecsq

***********************************************************************
*************************** INITIAL-FINAL *****************************
***********************************************************************
      elseif ((ip <= 2) .and. (kp > 2)) then
        u=sij/(sij+sik)

        omx=-sjk/(sij+sik)
        x=one-omx
        omu=sik/(sij+sik)
C---npart is the number of particles in the final state
C---transform the momenta so that only the first npart+1 are filled
        call transform(p,ptrans,x,ip,jp,kp)
        call storeptilde(nd,ptrans)

c--- if using a dynamic scale, set that scale with dipole kinematics      
      if (dynamicscale) then
        call scaleset(initscale,initfacscale,ptrans)
        dipscale(nd)=facscale
      endif
      
c--- Calculate the matrix element now because it might be needed
c--- in the final-initial segment, regardless of whether or not the
c--- alfa cut fails here
        call subr_born(ptrans,msq,mqq,msqx,mg)
C---Modification so that only close to singular subtracted
C---Do not set incldip because initial-final can fail 
C---but final initial needs still to be tested
        if (u > aif) return

        do nu=1,4
           vec(nu)=p(jp,nu)/u-p(kp,nu)/omu
        enddo

        call subr_corr(ptrans,vec,ip,msqv,mvg,mvxg)        
        sub(qq)=-gsq/x/sij*(two/(omx+u)-one-x)
        sub(gq)=-gsq/sij
        sub(qg)=-gsq/x/sij*(one-two*x*omx)
        sub(gg)=-2._dp*gsq/x/sij*(one/(omx+u)-one+x*omx)
        subv   =-4._dp*gsq/x/sij*(omx/x*u*(one-u)/sjk)

***********************************************************************
*************************** FINAL-INITIAL *****************************
***********************************************************************
      elseif ((ip > 2) .and. (kp <= 2)) then
c--- note, here we assume that msq kinematics are already taken care of
c--- for msq, although msqv must be recalculated each time
        omx=-sij/(sjk+sik)
C---Modification so that only close to singular subtracted
        if (omx > afi) return

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
        if (jp .ne.7) then
          if (ip < 7) then
C ie for cases 56_i,65_i
          call subr_corr(ptrans,vec,5,msqv,mvg,mvxg)
          else
C ie for cases 76_i,75_i
          call subr_corr(ptrans,vec,6,msqv,mvg,mvxg)
          endif
        else
C ie for cases 57_i,67_i
          call subr_corr(ptrans,vec,ip,msqv,mvg,mvxg)
        endif
                
        sub(qq)=+gsq/x/sij*(two/(omz+omx)-one-z)
        sub(gq)=+gsq/x/sij
        sub(gg)=+2._dp*gsq/x/sij*(one/(omz+omx)+one/(z+omx)-two) 
        subv   =+4._dp*gsq/x/sij/sij


***********************************************************************
**************************** FINAL-FINAL ******************************
***********************************************************************
      elseif ((ip > 2) .and. (kp > 2)) then
c------Eq-(5.2)    
       y=sij/(sij+sjk+sik)

C---Modification so that only close to singular subtracted
        if (y > aff) then
          incldip(nd)=.false.
          return
        endif

       z=sik/(sjk+sik)
       omz=one-z
       omy=one-y
C---calculate the ptrans-momenta 

       call transform(p,ptrans,y,ip,jp,kp)
       call storeptilde(nd,ptrans)

       do nu=1,4
         vec(nu)=z*p(ip,nu)-omz*p(jp,nu)
       enddo

c--- if using a dynamic scale, set that scale with dipole kinematics      
      if (dynamicscale) then
        call scaleset(initscale,initfacscale,ptrans)
        dipscale(nd)=facscale
      endif
      
       call subr_born(ptrans,msq,mqq,msqx,mg)
       if (ip < kp) then
         call subr_corr(ptrans,vec,5,msqv,mvg,mvxg)
       else
         call subr_corr(ptrans,vec,6,msqv,mvg,mvxg)
       endif
              
       sub(qq)=gsq/sij*(two/(one-z*omy)-one-z)
       sub(gq)=gsq/sij
       sub(gg)=gsq/sij*(two/(one-z*omy)+two/(one-omz*omy)-four)
       subv   =+4._dp*gsq/sij/sij

      endif
      
      return
      end
      
