      subroutine pvsetmudim(mu)
      implicit none
      include 'TRscale.f'
      double precision mu
      scale=mu
      musq=scale**2
      return
      end
