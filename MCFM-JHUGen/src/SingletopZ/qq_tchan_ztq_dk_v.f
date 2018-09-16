      subroutine qq_tchan_ztq_dk_v(q,msq)
      implicit none
      include 'types.f'
c--- Virtual contribution averaged over initial colors and spins
c     u(-q1)+b(-q2)->e-(q3)+e+(q4)+t(->nu(q5) l+(q6) b(q7)) +d(q8)
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'epinv.f'
      include 'scheme.f'
      include 'alpha1.f'
      include 'scale.f'
      include 'ewcouple.f'
      include 'nwz.f'
      include 'qcdcouple.f'
      include 'masses.f'
      include 'swapxz.f'
      include 'TRbadpoint.f'
      include 'TRtensorcontrol.f'
      include 'tensorinfo.f'
      include 'first.f'
      integer:: nu,icross,u_b,b_u,db_b,b_db
      integer:: origitotal,origibadpoint,origipolesfailed
      real(dp):: p(mxpart,4),msq(-nf:nf,-nf:nf),fac,
     & virt(4),q(mxpart,4)
      real(dp):: p_dk(mxpart,4),ps(mxpart,4),ptDpep
      complex(dp):: xupper(-2:0),xlower(-2:0), vupper(2,-2:0),
     & vlower(2,-2:0),vmiddle(2,-2:0),vextra(2,-2:0),vscalar(2,-2:0),
     & vtot(2,-2:0),vtotep(2),upper_tri(2,-2:0), lower_tri(2,-2:0),
     &  lotot(2),lotot1(2,2)
      integer:: j3,h3,h5
      logical:: failed
      parameter(u_b=1,b_u=2,db_b=3,b_db=4)
      integer,parameter:: i1(4)=(/1,2,8,8/)
      integer,parameter:: i2(4)=(/2,1,2,1/)
      integer,parameter:: i8(4)=(/8,8,1,2/)
      integer:: p3,p4,k5,e5,p6,eta
      parameter(p3=3,p4=4,k5=5,e5=6,p6=7,eta=6)
      complex(dp):: mdecaymb(2,2),mdecay

      scheme='dred'

      msq(:,:)=zip
      
      if (nwz == +1) then
      call tdecay(q,5,6,7,mdecaymb)
      mdecay=mdecaymb(1,1)
      elseif (nwz == -1) then
      call adecay(q,5,6,7,mdecaymb)
      mdecay=mdecaymb(2,2)
      endif

c--- gauge parameter: alpha1=0->unitary gauge, alpha1=1->Feynman gauge        
c--- Note that we must be in an EW scheme where mW=mZ*cosW                   
c--- in order for the alpha=0 and alpha=1 results to be identical             
c--- (and for poles to cancel properly)                                       
      alpha1=one
                      
c--- fill strings of gamma matrices on first pass
      if (first) then
        swapxz=.true.
        call fillgam
        call TRsetmaxindex(2,3,0) ! Only need up to Cij and Dijk
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
c -- momenta after decay
      p_dk(1,nu)=q(i1(icross),nu)
      p_dk(2,nu)=q(i2(icross),nu)
      p_dk(3,nu)=q(3,nu)
      p_dk(4,nu)=q(4,nu)
      if (nwz == -1) then
         p_dk(3,nu)=q(4,nu)
         p_dk(4,nu)=q(3,nu)
      endif
      p_dk(5,nu)=q(5,nu)
      p_dk(6,nu)=q(6,nu)
      p_dk(7,nu)=q(7,nu)
      p_dk(8,nu)=q(i8(icross),nu)
      
c -- momenta prior to decay
      p(1,nu)=p_dk(1,nu)
      p(2,nu)=p_dk(2,nu)
      p(3,nu)=p_dk(3,nu)
      p(4,nu)=p_dk(4,nu)
      p(5,nu)=p_dk(5,nu)+p_dk(6,nu)+p_dk(7,nu)
      p(6,nu)=p_dk(8,nu)
      enddo

c -- split pt into massless components
      ptDpep=p(5,4)*p_dk(eta,4)-p(5,1)*p_dk(eta,1)-p(5,2)*p_dk(eta,2)-
     & p(5,3)*p_dk(eta,3)
      ps(1,:)=p(1,:)
      ps(2,:)=p(2,:)
      ps(p3,:)=p(3,:)
      ps(p4,:)=p(4,:)
      ps(k5,:)=p(5,:)-mt**2*p_dk(eta,:)/(two*ptDpep)
      ps(e5,:)=mt**2*p_dk(eta,:)/(two*ptDpep)
      ps(p6,:)=p(6,:)

      call kininv(p)

c--- fill leading order amplitudes
c--- no need to swap momenta as this has happened above already ??
      call ubtzdamp_dk(p_dk,1,2,3,4,8,lotot1)
      if (nwz == 1) then
      lotot(:)=lotot1(:,1)
      elseif (nwz == -1) then
      lotot(:)=-conjg(lotot1(:,1))
      endif


      call scalardk(ps,vscalar)
      call extradk(ps,vextra)
      call middledk(ps,vmiddle)

      do j3=1,2
c----- lepton helicity h3, top spin is always -1
      h3=2*j3-3
      h5=-1
      
      if (nwz == +1) then
      call strings_dk(p,p_dk,h3,h5,.false.)
      elseif (nwz == -1) then
      call strings_dk(p,p_dk,h3,h5,.true.)
      endif

c--- Upper and lower: numerical for box contribution only
      call upperdk_partbox(p,h3,xupper,first)

      call lowerdk_partbox(p,h3,xlower,first)

      vupper(j3,:)= xupper(:)
      vlower(j3,:)= xlower(:)
      enddo
      
      itotal=itotal+1
      if (pvbadpoint) then
c        write(6,*) 'badpoint set'
      ibadpoint=ibadpoint+1
      goto 99
      endif
c--- Upper and lower: analytic for remaining contributions
      
      call upperdk_parttri(ps,upper_tri,first)
      call lowerdk_parttri(ps,lower_tri,first)

      vupper(:,:)=vupper(:,:)+upper_tri(:,:)
      vlower(:,:)=vlower(:,:)+lower_tri(:,:)
      
      vtot(:,:)=vupper(:,:)+vlower(:,:)
     & +vmiddle(:,:)+vscalar(:,:)+vextra(:,:)

      call dopolesdk(p,'total',lotot,vtot,first,failed)

      if (failed) ipolesfailed=ipolesfailed+1
c        if (mod(itotal,100000) == 0) then
c           write(6,*) itotal,': bad points for PV reduction ',
c     &          real(ibadpoint,dp)/real(itotal,dp)*100._dp,'%'
c           write(6,*) itotal,': poles failing check ',
c     &          real(ipolesfailed,dp)/real(itotal,dp)*100._dp,'%'
c           call flush(6)
c        endif
c---  this block of code rejects PS point if pole check fails
        if (failed) then
           pvbadpoint=.true.
           goto 99 
        endif        
        if (first) first=.false.

c---  wave function renormalization
        vtot(:,-1)=vtot(:,-1)-three*lotot(:)    
        vtot(:, 0)=vtot(:, 0)
     &       -(three*log(musq/mt**2)+five)*lotot(:)
      
c--- couplings are applied here
      fac=two*Cf*ason2pi*esq**2*gwsq**4*xn**2

      vtotep(:)=vtot(:,-2)*epinv**2+vtot(:,-1)*epinv+vtot(:,0)
                  
c -- include LO decay matrix elements      
      vtotep(:)=vtotep(:)*mdecay
      lotot(:)=lotot(:)*mdecay
     

      virt(icross)=aveqq*fac*real(
     & +vtotep(1)*conjg(lotot(1))+vtotep(2)*conjg(lotot(2)))

      if (failed) virt(icross)=zip
     
      enddo      ! end of loop over crossings
            
c--- fill matrix elements
      if (nwz == +1) then
      msq(+2,5)=virt(u_b)
      msq(+4,5)=virt(u_b)
      msq(+5,+2)=virt(b_u)
      msq(+5,+4)=virt(b_u)
      msq(-1,+5)=virt(db_b)
      msq(-3,+5)=virt(db_b)
      msq(+5,-1)=virt(b_db)
      msq(+5,-3)=virt(b_db)
      elseif ( nwz == -1) then
      msq(-2,-5)=virt(u_b)
      msq(-4,-5)=virt(u_b)
      msq(-5,-2)=virt(b_u)
      msq(-5,-4)=virt(b_u)
      msq(+1,-5)=virt(db_b)
      msq(+3,-5)=virt(db_b)
      msq(-5,+1)=virt(b_db)
      msq(-5,+3)=virt(b_db)
      endif

   99 continue

c--- This block repeats calculation using PV if requested
      if (  (TRtensorcontrol == 1) .and. (doovred .eqv. .true.)
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
      
      




      subroutine dopolesdk(p,desc,lo,virt,first,failed)
      implicit none
      include 'types.f'
c -- same as dopoles, but no label for top spin
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'scale.f'
      include 'masses.f'
      integer::j3
      real(dp):: p(mxpart,4),dot
      character*5 desc
      complex(dp):: lo(2),virt(2,-2:0)
      complex(dp):: dpnum,spnum
      logical:: first,failed
      real(dp):: tol
      
      failed=.false.
      tol=1e-5_dp
      
      do j3=1,2
         dpnum=virt(j3,-2)
         spnum=virt(j3,-1)

         virt(j3,-2)=-two*lo(j3)*three
         virt(j3,-1)=two*lo(j3)*(-eleven/two
     &        -two*log(musq/(-two*dot(p,1,6)))
     &        -two*log(musq/(-two*dot(p,2,5)))
     &        +log(musq/mt**2)
     &        +three/two)         ! Last term is from w.f. renorm
     
         if ( abs(spnum/virt(j3,-1)-one) < tol ) then
            continue
         else
            failed=.true.
            if (first) then
              write(6,*) 'Problem with poles in ', desc
              write(6,*) 'double pole', virt(j3,-2)/lo(j3),
     &           dpnum/virt(j3,-2)
              write(6,*) 'single pole', virt(j3,-1)/lo(j3),
     &           spnum/lo(j3),  spnum/virt(j3,-1)
c     pause
            endif
         endif
      enddo
            
      if ((first) .and. (failed .eqv. .false.)) write(6,*)
     & 'Poles computed correctly in: ',desc,' with tolerance ', tol
      
      return
      end
