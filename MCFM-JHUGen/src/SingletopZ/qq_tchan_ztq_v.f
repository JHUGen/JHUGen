      subroutine qq_tchan_ztq_v(q,msq)
      implicit none
      include 'types.f'
c--- Virtual contribution averaged over initial colors and spins
c     u(-p1)+b(p2)->e-(p3)+e+(p4)+t(p5)+d(p6)
      
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
     & virt(4),q(mxpart,4),pswap(mxpart,4)
      complex(dp):: xupper(-2:0),xlower(-2:0),
     & vupper(2,2,-2:0),vlower(2,2,-2:0),vmiddle(2,2,-2:0),
     & vextra(2,2,-2:0),vscalar(2,2,-2:0),vtot(2,2,-2:0),vtotep(2,2),
     & upper_tri(2,2,-2:0), lower_tri(2,2,-2:0)
      complex(dp):: lotot(2,2)
      integer:: j3,j5,h3,h5
      logical:: failed
      parameter(u_b=1,b_u=2,db_b=3,b_db=4)
      integer, parameter:: i1(4)=(/1,2,6,6/)
      integer, parameter:: i2(4)=(/2,1,2,1/)
      integer, parameter:: i6(4)=(/6,6,1,2/)
       
      scheme='dred'

      msq(:,:)=zip

c--- gauge parameter: alpha1=0->unitary gauge, alpha1=1->Feynman gauge 
c--- Note that we must be in an EW scheme where mW=mZ*cosW
c--- in order for the alpha=0 and alpha=1 results to be identical
c--- (and for poles to cancel properly)
      alpha1=1d0

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
      p(1,nu)=q(i1(icross),nu)
      p(2,nu)=q(i2(icross),nu)
      p(3,nu)=q(3,nu)
      p(4,nu)=q(4,nu)
      p(5,nu)=q(5,nu)
      p(6,nu)=q(i6(icross),nu)
      enddo
      call kininv(p)
      
c --- swap momenta p3 <-> p4 for tbar
      pswap=p
      if (nwz == -1) then
         pswap(3,:)=p(4,:)
         pswap(4,:)=p(3,:)
      endif
      
      if (nwz == 1) then
         call ubtzdamp(p,1,2,3,4,6,lotot)
      elseif (nwz == -1) then
         call ubtzdamp(p,1,2,4,3,6,lotot)
      endif

c--- setup currents
      do j5=1,2
      do j3=1,2
c----- lepton helicity h3, top spin h5
      h3=2*j3-3
      h5=2*j5-3

      if (nwz == +1) then
      call strings(p,h3,h5,.false.)
      elseif (nwz == -1) then
      call strings(p,h3,h5,.true.)
      endif
         
c--- Upper and lower: numerical for box contribution only
      call upper_partbox(p,h3,xupper,first)
      call lower_partbox(p,h3,xlower,first)
      
      vupper(j3,j5,:)= xupper(:)
      vlower(j3,j5,:)= xlower(:)

      enddo
      enddo
      
      itotal=itotal+1
      if (pvbadpoint) then
c        write(6,*) 'badpoint set'
        ibadpoint=ibadpoint+1
        goto 99
      endif

c --- Upper and lower: analytic for remaining contributions
      call upper_parttri(pswap,upper_tri,first)
      vupper(:,:,:)=vupper(:,:,:)+upper_tri(:,:,:)
      call lower_parttri(pswap,lower_tri,first)
      vlower(:,:,:)=vlower(:,:,:)+lower_tri(:,:,:)

c --- Scalar, middle, extra: analytic only
      call scalar(pswap,vscalar)
      call middle(pswap,vmiddle)
      call extra(pswap,vextra)

c--- total virtual amplitudes
      vtot(:,:,:)=vupper(:,:,:)+vlower(:,:,:)
     & +vmiddle(:,:,:)+vscalar(:,:,:)+vextra(:,:,:)

c -- Pole check
      call dopoles(p,'total',lotot,vtot,first,failed)

      if (failed) ipolesfailed=ipolesfailed+1
c      if (mod(itotal,100000) == 0) then
c         write(6,*) itotal,': bad points for PV reduction ',
c     &        real(ibadpoint,dp)/real(itotal,dp)*100d0,'%'
c         write(6,*) itotal,': poles failing check ',
c     &        real(ipolesfailed,dp)/real(itotal,dp)*100d0,'%'
c         call flush(6)
c      endif
c---  this block of code rejects PS point if pole check fails
      if (failed) then
         pvbadpoint=.true.
         goto 99 
      endif
      if (first) first=.false.

c--- wave function renormalization
      vtot(:,:,-1)=vtot(:,:,-1)-3d0*lotot(:,:)    
      vtot(:,:, 0)=vtot(:,:, 0)
     &     -(3d0*log(musq/mt**2)+5d0)*lotot(:,:)

c--- couplings are applied here
      fac=2d0*Cf*ason2pi*esq**2*gwsq**2*xn**2
      
      vtotep(:,:)=vtot(:,:,-2)*epinv**2+vtot(:,:,-1)*epinv+vtot(:,:,0)
                      
      virt(icross)=aveqq*fac*real(
     & +vtotep(1,1)*conjg(lotot(1,1))+vtotep(1,2)*conjg(lotot(1,2))
     & +vtotep(2,1)*conjg(lotot(2,1))+vtotep(2,2)*conjg(lotot(2,2)))

      if (failed) virt(icross)=0d0

      enddo        ! end of loop over crossings

c--- fill matrix elements
      if ( nwz == +1) then
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
      
      


      subroutine dopoles(p,desc,lo,virt,first,failed)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'scale.f'
      include 'masses.f'
      real(dp):: p(mxpart,4),dot
      integer::j3,j5
      character*5 desc
      complex(dp):: lo(2,2),virt(2,2,-2:0)
      complex(dp):: spnum
      logical:: first,failed
      real(dp):: tol
      
      failed=.false.
      tol=1d-5
      
      do j3=1,2
         do j5=1,2
      
            spnum=virt(j3,j5,-1)
      
            virt(j3,j5,-2)=-2d0*lo(j3,j5)*3d0
            virt(j3,j5,-1)=2d0*lo(j3,j5)*(-11d0/2d0
     &           -2d0*log(musq/(-2d0*dot(p,1,6)))
     &           -2d0*log(musq/(-2d0*dot(p,2,5)))
     &           +log(musq/mt**2)
     &           +3d0/2d0)      ! Last term is from w.f. renorm
     
            if ( abs(spnum/virt(j3,j5,-1)-1d0) < tol ) then
               continue
            else
               failed=.true.
            endif

         enddo
      enddo
      
      if (first) write(6,*) 'Poles computed correctly in: ',desc,
     & ' with tolerance ', tol
      
      return
      end
