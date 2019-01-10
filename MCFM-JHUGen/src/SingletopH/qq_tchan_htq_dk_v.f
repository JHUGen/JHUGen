      subroutine qq_tchan_htq_dk_v(q,msq)
c--- Virtual contribution averaged over initial colors and spins
c     u(-p1)+b(p2)->h(p3,p4)+t(nu(p5)+e+(p6)+b(p7))+d(p6)
      implicit none
      include 'constants.f'
      include 'epinv.f'
      include 'scheme.f'
      include 'alpha1.f'
      include 'scale.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'hdecaymode.f'
      include 'masses.f'
      include 'swapxz.f'
      include 'nwz.f'
      include 'decl_kininv.f'
      include 'TRbadpoint.f'
      include 'TRtensorcontrol.f'
      include 'tensorinfo.f'
      include 'anomHiggs.f'
      include 'first.f'
      integer nu,icross,u_b,b_u,db_b,b_db
      integer origipolesfailed,origitotal,origibadpoint
      double precision p(mxpart,4),msq(-nf:nf,-nf:nf),fac,
     & virt(4),q(mxpart,4),p_dk(mxpart,4),hdecay,msqgamgam
      double complex vtot(2,-2:0),vtotep(2),
     & xlower(-2:0),xmiddle(-2:0),xscalar(-2:0),
     & vlower(2,-2:0),vmiddle(2,-2:0),vscalar(2,-2:0)
      double complex lotot(2)
      integer j5,h5
      logical failed
      parameter(u_b=1,b_u=2,db_b=3,b_db=4)
      integer,parameter:: i1(4)=(/1,2,8,8/)
      integer,parameter:: i2(4)=(/2,1,2,1/)
      integer,parameter:: i8(4)=(/8,8,1,2/)
! RR added
      double complex mdecaymb(2,2), mdecay
            
      scheme='dred'

      msq(:,:)=zip
c--- define alpha1
      alpha1=1d0
c--- Note that we must be in an EW scheme where mW=mZ*cosW
c--- in order for the alpha=0 and alpha=1 results to be identical
c--- (and for poles to cancel properly)
      
C   Deal with Higgs decay
      s34=(q(3,4)+q(4,4))**2-(q(3,1)+q(4,1))**2
     &   -(q(3,2)+q(4,2))**2-(q(3,3)+q(4,3))**2
      if     (hdecaymode == 'tlta') then
          call htautaudecay(q,3,4,hdecay)
      elseif (hdecaymode == 'bqba') then
          call hbbdecay(q,3,4,hdecay)
      elseif (hdecaymode == 'gaga') then
          hdecay=msqgamgam(hmass)
      else
      write(6,*) 'Unimplemented process'
      stop
      endif
      hdecay=hdecay/((s34-hmass**2)**2+(hmass*hwidth)**2)


c--- fill strings of gamma matrices on first pass
      if (first) then
        swapxz=.true.
        call fillgam
        call TRsetmaxindex(2,2,0) ! Only need up to Cij and Dij
        call qlinit
        call pvArraysetup
        itotal=0
        ibadpoint=0
        ipolesfailed=0
      endif

c--- setup for pvQCDLoop
      call pvsetmudim(scale)     
      call TRsettensorcontrol(1)

      origitotal=itotal
      origibadpoint=ibadpoint
      origipolesfailed=ipolesfailed

   66 continue
      call pvclearcache

c--- loop over all required crossings of amplitude
      do icross=1,4

      do nu=1,4
      p_dk(1,nu)=q(i1(icross),nu)
      p_dk(2,nu)=q(i2(icross),nu)
      p_dk(3,nu)=q(3,nu)
      p_dk(4,nu)=q(4,nu)
      p_dk(5,nu)=q(5,nu)
      p_dk(6,nu)=q(6,nu)
      p_dk(7,nu)=q(7,nu)
      p_dk(8,nu)=q(i8(icross),nu)

      p(1,nu)=p_dk(1,nu)
      p(2,nu)=p_dk(2,nu)
      p(3,nu)=p_dk(3,nu)
      p(4,nu)=p_dk(4,nu)
      p(5,nu)=p_dk(5,nu)+p_dk(6,nu)+p_dk(7,nu)
      p(6,nu)=p_dk(8,nu)
      enddo
c--- setup the kinematics for the various crossings

      call kininv(p)
      
      if ( nwz .eq. +1) then
      call tdecay(p_dk,5,6,7,mdecaymb)
      mdecay=mdecaymb(1,1)
      elseif (nwz .eq. -1) then
      call adecay(p_dk,5,6,7,mdecaymb)
      mdecay=mdecaymb(2,2)
      endif
      
c--- fill leading order amplitudes
      if     (nwz .eq. +1) then
        call ubthdamp_dk(p,1,2,3,4,6,p_dk(6,:),
     &        lotot)
      elseif (nwz .eq. -1) then
        call ubthdamp_dk(p,1,2,4,3,6,p_dk(6,:),
     &        lotot)
      else
       write(6,*) 'Problem with nwz in qq_tchan_htq_v.f: nwz=',nwz
       stop
      endif
      
c--- setup currents
      do j5=1,1
c----- top spin h5
      h5=2*j5-3

      if (nwz .eq. +1) then
        call stringsh_dk(p,dcmplx(p_dk(6,:)),h5,.false.)
      else
        call stringsh_dk(p,dcmplx(p_dk(6,:)),h5,.true.)
      endif
      call scalarh(p,first,xscalar)
      call middleh(p,first,xmiddle)      
      call lowerh(p,first,xlower)


c --   include anomalous Yukawa and EW couplings for the Higgs
      xlower=cttH*xlower
      xscalar=cWWH*xscalar
      xmiddle=cWWH*xmiddle

      vmiddle(j5,:)=xmiddle(:)
      vscalar(j5,:)=xscalar(:)
      vlower (j5,:)=xlower (:)

      enddo
      
      itotal=itotal+1
      if (pvbadpoint) then
c        write(6,*) 'badpoint set'
        ibadpoint=ibadpoint+1
        goto 99
      endif

c--- total virtual amplitudes
      vtot(:,:)=vlower(:,:)+vmiddle(:,:)+vscalar(:,:)

c -- now include propagator of onshell top
      vtot(:,:)=vtot(:,:)/dcmplx(zip,mt*twidth)
      lotot(:)=lotot(:)/dcmplx(zip,mt*twidth)

c--- block of code for checking poles
      call dopolesh_dk(p,'total',lotot,vtot,first,failed)
      if (failed) ipolesfailed=ipolesfailed+1
c        if (mod(itotal,10000) .eq. 0) then
c          write(6,*) itotal,': bad points for PV reduction ',
c     &     dfloat(ibadpoint)/dfloat(itotal)*100d0,'%'
c          write(6,*) itotal,': poles failing check ',
c     &     dfloat(ipolesfailed)/dfloat(itotal)*100d0,'%'
c          call flush(6)
c        endif
c--- this block of code rejects PS point if pole check fails
      if (failed) then
         pvbadpoint=.true.
         goto 99 
      endif
c--- this block of code rejects PS point if pole check fails

      if (first) first=.false.


c--- wave function renormalization
      vtot(:,-1)=vtot(:,-1)-3d0/2d0*lotot(:)    
      vtot(:, 0)=vtot(:, 0)
     & -(3d0*dlog(musq/mt**2)+5d0)/2d0*lotot(:)

c--- couplings are applied here
      fac=Cf*ason2pi*gwsq**5*xn**2*hdecay
c -- decay ME included here
      vtot(:,:)=vtot(:,:)*mdecay
      lotot(:)=lotot(:)*mdecay
          
      vtotep(:)=vtot(:,-2)*epinv**2+vtot(:,-1)*epinv+vtot(:,0)

      virt(icross)=aveqq*fac*dble(
     & +vtotep(1)*dconjg(lotot(1))+vtotep(2)*dconjg(lotot(2)))
           
      enddo    ! end of loop over crossings
      
c--- fill matrix elements
c--- note: variable names are not correct for t~ (nwz=-1)
      if (nwz .eq. +1) then
        msq(+2,5)=virt(u_b)
        msq(+4,5)=virt(u_b)
        msq(+5,+2)=virt(b_u)
        msq(+5,+4)=virt(b_u)
      
        msq(-1,+5)=virt(db_b)
        msq(-3,+5)=virt(db_b)
        msq(+5,-1)=virt(b_db)
        msq(+5,-3)=virt(b_db)
      else
        msq(+1,-5)=virt(u_b)
        msq(+3,-5)=virt(u_b)
        msq(-5,+1)=virt(b_u)
        msq(-5,+3)=virt(b_u)
      
        msq(-2,-5)=virt(db_b)
        msq(-4,-5)=virt(db_b)
        msq(-5,-2)=virt(b_db)
        msq(-5,-4)=virt(b_db)
      endif
      

   99 continue

c--- This block repeats calculation using PV if requested
      if (  (TRtensorcontrol .eq. 1) .and. (doovred .eqv. .true.)
     &.and. (pvbadpoint) ) then
        itotal=origitotal
        ibadpoint=origibadpoint
        ipolesfailed=origipolesfailed
        doovred=.false.
        dopvred=.true.
        goto 66
      endif

      return
      end
      
      



      subroutine dopolesh_dk(p,desc,lo,virt,first,failed)
      include 'constants.f'
      include 'scale.f'
      include 'masses.f'
      double precision p(mxpart,4),dot
      character*5 desc
      double complex lo(2),virt(2,-2:0)
      double complex spnum
      logical first,failed
      double precision tol
      
      failed=.false.
      tol=1d-5
      
      do j5=1,1
      
      spnum=virt(j5,-1)
      
      virt(j5,-2)=-lo(j5)*3d0
      virt(j5,-1)=lo(j5)*(-11d0/2d0
     & -2d0*dlog(musq/(-2d0*dot(p,1,6)))
     & -2d0*dlog(musq/(-2d0*dot(p,2,5)))
     & +dlog(musq/mt**2)
     & +3d0/2d0) ! Last term is from w.f. renorm
           
      if ( abs(spnum/virt(j5,-1)-1d0) .lt. tol ) then
        continue
      else
        failed=.true.
      endif

      enddo
      
      if (first) write(6,*) 'Poles computed correctly in: ',desc,
     & ' with tolerance ', tol
     
      
      return
      end
