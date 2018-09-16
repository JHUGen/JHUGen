************************************************************************ 
*     This subroutine calculates dipoles with an identified            *
*     Final state photon                                               *
*     C. Williams Dec 2010                                             * 
*     Returns the dipoles in sub,subv and matrix elements in msq,msqv  *
*     nd labels the dipole configurations                              *
*     ip labels the emitter PHOTON                                     *
*     jp labels the emitted parton                                     *
*     kp labels the spectator parton                                   *
*     subr_born is the subroutine which call the born process          *
*     subr_corr is the subroutine which call the born process dotted   *
*     with vec for an emitted gluon only                               *
************************************************************************

c--- Extension to dipolesfrag that also returns msqx,
c--- 4 dimensional array indexed by initial and final parton labels

      subroutine dipsfragx(nd,p,ip,jp,kp,sub,msq,msqx,subr_born) 
      implicit none
      include 'types.f'
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ewcouple.f'
      include 'ptilde.f'
      include 'dynamicscale.f'
      include 'initialscales.f' 
      include 'dipolescale.f'
      include 'facscale.f'
      include 'betacut.f'
      include 'lastphot.f'
      include 'incldip.f'
      real(dp):: p(mxpart,4),ptrans(mxpart,4),sub
      real(dp):: p_phys(mxpart,4),tmp
      integer:: ipt
      real(dp):: z,omz,sij,sik,sjk,dot,u
      real(dp):: msq(-nf:nf,-nf:nf)
      real(dp):: mdum1(0:2,fn:nf,fn:nf) ! mqq
      real(dp):: msqx(0:2,-nf:nf,-nf:nf,-nf:nf,-nf:nf)
      real(dp):: mdum2(0:2,-nf:nf,-nf:nf) !msqx_cs
      integer:: nd,ip,jp,kp,j,k
      logical:: check_nv,phot_pass
      external subr_born

      z=0._dp
      omz=1._dp
      u=0._dp
      sub=0._dp
      do j=-nf,nf
         do k=-nf,nf
            msq(j,k)=0._dp
         enddo
      enddo
      
      incldip(nd)=.true.

      sij=two*dot(p,ip,jp)
      sik=two*dot(p,ip,kp) 
      sjk=two*dot(p,jp,kp)

     
******************************************************************************* 
************************ INITIAL - INITIAL ************************************
*******************************************************************************

**** I === I not implemented (need photon PDFS rather than frags) need rapidity cuts to remove collinear sing 
      if ((ip <= 2) .and. (kp <= 2)) return  

******************************************************************************* 
************************ INITIAL - FINAL **************************************
*******************************************************************************

**** I === F not implemented (need photon PDFS rather than frags) need rapidity cuts to remove collinear sing
      if ((ip <= 2) .and. (kp > 2)) return  

******************************************************************************* 
************************ FINAL - INITIAL  *************************************
*******************************************************************************

      
      if ((ip > 2) .and. (kp <= 2)) then
         
     
         if(check_nv(p,ip,jp,kp).eqv..false.) then 
            incldip(nd)=.false.
            return 
         endif
         call transformfrag(p,ptrans,z,ip,jp,kp)
         ipt=ip
         if (ip < lastphot) then
             do j=1,4
             tmp=ptrans(ip,j)
             ptrans(ip,j)=ptrans(lastphot,j)
             ptrans(lastphot,j)=tmp
             enddo
             ipt=lastphot
         endif

        if (dynamicscale) then
c--- rescale momentum (c.f. code below) in order to obtain physical
c--- momentum for setting dynamic scale
           p_phys(:,:)=ptrans(:,:)
           p_phys(lastphot,:)=z*ptrans(lastphot,:)
           call scaleset(initscale,initfacscale,p_phys)
           dipscale(nd)=facscale
        endif
        
c-- here, ptrans contains the right momentum for evaluating
c-- the LO matrix element with the last photon replaced by a jet 
         call subr_born(ptrans,msq,mdum1,msqx,mdum2)
         omz=one-z
         sub=two*(esq/sij)*((one+omz**2)/z)
         
c--- now rescale momentum of last photon entry to represent
c--- the observed photon momentum
c--- NB: momentum will no longer be conserved in ptrans  
         ptrans(lastphot,:)=z*ptrans(lastphot,:)
         
         call storeptilde(nd,ptrans) 
c--- store z for use in isolation routine         
         call store_zdip(nd,z)

         
        
         
******************************************************************************* 
************************ FINAL - FINAL ****************************************
*******************************************************************************
         
      elseif ((ip > 2) .and. (kp > 2)) then 
      

        write(6,*) 'Final-final fragmentation dipole not implemented.'
        stop
 
         z=(sij+sik)/(sik+sjk+sij) 
         omz=one-z
         
         u=sij/(sij+sik) 
       
         if(u>bff) then
            incldip(nd)=.false.
            return 
         endif
            
c--- Calculate the ptrans-momenta 
         call transformfrag(p,ptrans,z,ip,jp,kp)
         call storeptilde(nd,ptrans) 
         call store_zdip(nd,z)

c----  Check that photon will still pass cuts 
         if (phot_pass(ptrans,ip,z) .eqv. .false.) return 
c--- if using a dynamic scale, set that scale with dipole kinematics        
        if (dynamicscale) then
           call rescale_z_dip(ptrans,nd,ip)
           call scaleset(initscale,initfacscale,ptrans)
           call return_z_dip(ptrans,nd,ip)
           dipscale(nd)=facscale
        endif


         call subr_born(ptrans,msq,mdum1,msqx,mdum2)
      
         sub=two*(esq/sij)*((one+omz**2)/z)
      
      endif

c      write(6,*) 'dipole ',nd,' momenta'
c      call writeout(ptrans)

      return 
      end

