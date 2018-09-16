      subroutine scaleset_m34(p,mu0)
      implicit none
      include 'types.f'
c--- subroutine to calculate dynamic scale equal to
c--- invariant mass of particles 3 and 4
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'kprocess.f'
      real(dp):: p(mxpart,4),mu0

      if((kcase==kW_only) .or.
     &   (kcase==kZ_only) .or.
     &   (kcase==kW_1jet) .or.
     &   (kcase==kW_2jet) .or.
     &   (kcase==kW_3jet) .or.
     &   (kcase==kZ_1jet) .or.
     &   (kcase==kZ_2jet) .or.
     &   (kcase==kZ_3jet) .or.
     &   (kcase==kWbbbar) .or.
     &   (kcase==kWbbmas) .or.
     &   (kcase==kZbbbar) .or.
     &   (kcase==kgamgam) .or.
     &   (kcase==kgg2gam) .or.
     &   (kcase==kgmgmjt) .or.
     &   (kcase==kgagajj) .or.
     &   (kcase==kggfus0) .or.
     &   (kcase==kggfus1) .or.
     &   (kcase==kggfus2) .or.
     &   (kcase==kggfus3) .or.
     &   (kcase==kgagajj) .or.
     &   (kcase==kqg_tbq) .or.
     &   (kcase==ktt_tot) .or.
     &   (kcase==kbb_tot) .or.
     &   (kcase==kcc_tot) .or.
     &   (kcase==khttjet) .or.
     &   (kcase==kHigaga) .or.
     &   (kcase==kHgagaj) .or.
     &   (kcase==kqq_Hqq) .or.
     &   (kcase==kdm_jet) .or.
     &   (kcase==kdm_gam) .or.
     &   (kcase==kdm2jet) .or.
     &   (kcase==kdm_gaj) .or.
     &   (kcase==kdmgamj) .or.
     &   (kcase==kqq_Hgg) .or.
     &   (kcase==ktwojet)) then
        mu0=(p(3,4)+p(4,4))**2-(p(3,1)+p(4,1))**2
     &     -(p(3,2)+p(4,2))**2-(p(3,3)+p(4,3))**2       
        mu0=sqrt(abs(mu0))
      else
        write(6,*) 'dynamicscale m(34) not supported for this process.',kcase
        stop
      endif
      
      return
      end
      
