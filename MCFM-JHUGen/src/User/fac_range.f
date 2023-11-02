************************************************************************
*     Range for the dynamic factorization scale "rand"                 *
*      low must be <= high                                             *
************************************************************************
      block data fac_range
      implicit none
      include 'facscale_range.f'
      data facscale_low/10d0/
      data facscale_high/40d0/
      end
************************************************************************
