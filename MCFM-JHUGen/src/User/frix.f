!---------- Generic Frixione Routine, will take eps and delta_0 from input.DAT
!--------   n will be har.e-_dpcoded but can be changed below.
!--------- C. Williams July 2011


!----- p -momentum array
!----- passed, should be obvious!
!----- j - photon identification in p i.e. p(j,i) = photon(i)
!----- isub whether we are working with a dipole or not

!===== C. Williams July 2015
!===== extended to allow for multiple partons in the cone which can
!==== happen at NNLO,
      subroutine frix(p,passed,j,isub)
      implicit none
      include 'types.f'

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'frag.f'
      include 'npart.f'
      include 'first.f'
      include 'mpicommon.f'
      real(dp):: p(mxpart,4),R,pref,ret_ET
      integer:: j,isub,k
      logical:: passed,is_hadronic,in_cone_n
      integer:: i
      real(dp):: vsmall
      parameter(vsmall=1.e-10_dp)
      real(dp):: ET_had,pt
      integer:: itag

      passed=.true.

      if(first) then
         first=.false.
!----- check for non-zero parameters, if zero exit with warning
         if((epsilon_h<vsmall).or.(cone_ang<vsmall)) then
!$omp master
         if (rank == 0) then
       write(6,*)
       write(6,*)'************** Frixione Isolation    ***************'
       write(6,*)'*   Read zero parameters, not isolating            *'
       write(6,*) '* Warning, this may be unsafe in general *'
       write(6,99)'*  eps_phot = ',epsilon_h,' delta_0 = ',cone_ang, '*'
       write(6,97)'*  n = ',n_pow,'                                   *'
       write(6,*)'****************************************************'
         endif
!$omp end master
         return
         endif


!$omp master
       if (rank == 0) then
       write(6,*)
       write(6,*)'************** Frixione Isolation    ***************'
       write(6,*)'*                                                  *'
       write(6,99)'*  eps_phot = ',epsilon_h,', delta_0 = ',cone_ang,  '*'
       write(6,97)'*  n = ',n_pow,'                                   *'
       write(6,*)'****************************************************'
       endif
!$omp end master
      endif

 99   format (1x,a14,f5.3,a12,f5.3,a16)
! 97   format (1x,a7,i1,a44)
 97   format (1x,a7,f5.2,a40)

!----- Cycle over final state particles, if it is hadronic and inside
!----- photon cone then check its energy passes frixione requirement
!----- else fail and return

      pref=ret_ET(p,j)*epsilon_h/((1._dp-cos(cone_ang))**n_pow)
!      pref=epsilon_h/((1._dp-cos(cone_ang))**n_pow)
      ET_had=0_dp

!===== this section is altered now to allow for an additional parton
!===== < R_ij inside the cone too, which can happen at NNLO
      itag=0
      do i=3,2+npart-isub
!===== reset ET_had for each initial hadron
         ET_had =0_dp
!======first thing, find a hadron inside isolation cone
         if(is_hadronic(i).and.(R(p,i,j)<cone_ang)) then
            ET_had=ET_had+ret_ET(p,i)
!========= now we need to sum over the additional hadronic momenta looking for
!========= friends inside the cone
            do k=3,2+npart-isub
!========= is k nearer the photon then i?
               if(is_hadronic(k).and.(k.ne.i)) then
               if(R(p,k,j).lt.R(p,i,j)) then
                  ET_had=ET_had+ret_ET(p,k)
c                  itag=2
               endif
               endif
!====== include this hadron in total energy
            enddo
!========= now check total energy is isolated against
            passed=in_cone_n(R(p,i,j),ET_had,pref,n_pow)

            if(itag==2) then
               call writeout(p)
               write(6,*) 'photon = ',j
               write(6,*) 'pt(photon) = ',pt(j,p)
               write(6,*) 'R(5,j) = ',R(p,5,j)
               if(isub==0) write(6,*) 'R(6,j)= ',R(p,6,j)
               write(6,*) 'PT(5) = ',pt(5,p)
               write(6,*) 'PT(6) = ',pt(6,p)
               write(6,*) 'RET_ET(5)  = ',ret_ET(p,5)
               if(isub==0) write(6,*) 'RET_ET(6)  = ', ret_ET(p,6)

               write(6,*) 'total ET_had in cone',ET_had
               write(6,*) 'iso condition',pref*(1._dp-cos(R(p,i,j)))
               write(6,*) 'passed = ',passed
               pause
            endif
            if(passed.eqv..false.) return
         endif
      enddo

      return
      end




      function ret_ET(p,j)
      implicit none
      include 'types.f'
      real(dp):: ret_ET

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      real(dp):: p(mxpart,4)
      integer:: j
      real(dp):: ptsq

      ptsq=p(j,1)**2+p(j,2)**2
      ret_ET=p(j,4)*sqrt(ptsq)/(sqrt(ptsq+p(j,3)**2))

      if(ptsq==0) ret_ET=0_dp
      return
      end

      function in_cone_n(Rij,Ejet,pref,n)
       implicit none
      include 'types.f'
      logical:: in_cone_n

      real(dp):: Rij,Ejet,pref
!     integer:: n
      real(dp):: n

      if ( Ejet < pref*(1._dp-cos(Rij))**n) then
         in_cone_n = .true.
      else
         in_cone_n = .false.
      endif

      return
      end













