      subroutine read_jetcuts(read_ptmin,read_etamin,read_etamax)
      implicit none
      include 'clustering.f'
      include 'jetcuts.f'
      integer nargs
      character*72 jetcutsfile
      double precision read_ptmin,read_etamin,read_etamax
      double precision ptmin,etamax,ptmin_tev,etamax_tev,
     . ptmin_lhc,etamax_lhc
      logical useTevcuts,useLHCcuts
      
      read_ptmin=ptjetmin
      read_etamin=etajetmin
      read_etamax=etajetmax

      return
      end
      
      
      
