      include 'mxdim.f'
      integer:: ndim,itmx,nprn
      integer(kind=8) ncall
      real(dp):: xl(mxdim),xu(mxdim),acc
      common/bveg1int/ndim,itmx,nprn
      common/bveg1int8/ncall
      common/bveg1real/xl,xu,acc
