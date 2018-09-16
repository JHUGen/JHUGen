!====== virtual matrix element for Wgamjet
!====== C.W Nov 2013 
!+====== Currently a wrapper for testing routines
      subroutine qqb_wgamjet_v(p,msq)
      implicit none
      include 'types.f'
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h' 
      include 'zprods_decl.f' 
      include 'scale.f' 
      real(dp):: p(mxpart,4),msq(-nf:nf,-nf:nf) 
      real(dp):: phi,muk,rho,ssig,csig,theta,
     & p1true(4),p2true(4),p3true(4),p4true(4),p5true(4),p6true(4)  
      real(dp):: mu
      logical:: do_check 
      common/do_check_wgamj/do_check
      integer:: nu,om
      complex(dp):: wgamjet_vamp_q2lc_pp,test
      complex(dp):: wgamjet_vamp_q2slc_pp

      do_check=.true.
      call writeout(p)
      if(do_check) then 
!====== include KC
         include 'kinpoint.f'
         mu=1._dp
!===== sadly my mathematica kinpoint has a funny translation 
!===== compared to above!
         musq=mu**2
         do nu=1,4
            om=nu-1
            if (nu==1) om=4
            p(1,om)=p1true(nu)
            p(2,om)=p4true(nu)
            p(3,om)=p3true(nu)
            p(4,om)=p2true(nu)
            p(5,om)=p5true(nu)
            p(6,om)=p6true(nu)
         enddo
         do nu=1,6
            write(6,'(a4,i2,4f18.12)') 'p_',nu,
     &           p(nu,4),p(nu,1),p(nu,2),p(nu,3)
         enddo
         call spinorz(6,p,za,zb)
      endif

     
!======= check LC Q2 amplitude 
      test=wgamjet_vamp_q2lc_pp(1,2,3,4,5,6,za,zb)
      test=wgamjet_vamp_q2slc_pp(1,2,3,4,5,6,za,zb)
      pause
      return 
      end
