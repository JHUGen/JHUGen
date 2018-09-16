        function rn(idummy)
        implicit none
        include 'types.f'
        real(dp)::rn,ran
        integer::ij,kl,idummy
        integer,save::init=1
        common/seed/ij,kl

        if (init.eq.1) then
        init=0
        if (ij .eq. 0) then
        ij=1802
        kl=9373
        endif
        call rmarin(ij,kl)
        end if
*
  10    call ranmar(ran)
        if (ran.lt.1.e-16_dp) goto 10
        rn=ran
*
        end
*
      subroutine ranmar(rvec)
*     -----------------
* universal random number generator proposed by marsaglia and zaman
* in report fsu-scri-87-50
* in this version rvec is a double precision variable.
      implicit none
      include 'types.f'
      real(dp)::uni,rvec
      real(dp)::ranu(97),ranc,rancd,rancm
      common/ raset1 / ranu,ranc,rancd,rancm
      integer:: iranmr,jranmr
      common/ raset2 / iranmr,jranmr
      save /raset1/,/raset2/
      uni = ranu(iranmr) - ranu(jranmr)
      if(uni .lt. 0d0) uni = uni + 1d0
      ranu(iranmr) = uni
      iranmr = iranmr - 1
      jranmr = jranmr - 1
      if(iranmr .eq. 0) iranmr = 97
      if(jranmr .eq. 0) jranmr = 97
      ranc = ranc - rancd
      if(ranc .lt. 0d0) ranc = ranc + rancm
      uni = uni - ranc
      if(uni .lt. 0d0) uni = uni + 1d0
      rvec = uni
      end

      subroutine rmarin(ij,kl)
*     -----------------
* initializing routine for ranmar, must be called before generating
* any pseudorandom numbers with ranmar. the input values should be in
* the ranges 0<=ij<=31328 ; 0<=kl<=30081
      implicit none
      include 'types.f'
      integer:: i,j,k,l,ii,jj,ij,kl,m
      real(dp):: s,t
      real(dp)::ranu(97),ranc,rancd,rancm
      common/ raset1 / ranu,ranc,rancd,rancm
      integer:: iranmr,jranmr
      common/ raset2 / iranmr,jranmr
      save /raset1/,/raset2/
* this shows correspondence between the simplified input seeds ij, kl
* and the original marsaglia-zaman seeds i,j,k,l.
* to get the standard values in the marsaglia-zaman paper (i=12,j=34
* k=56,l=78) put ij=1802, kl=9373
      i = mod( ij/177 , 177 ) + 2
      j = mod( ij     , 177 ) + 2
      k = mod( kl/169 , 178 ) + 1
      l = mod( kl     , 169 )
      do 300 ii = 1 , 97
        s =  0d0
        t = .5d0
        do 200 jj = 1 , 24
          m = mod( mod(i*j,179)*k , 179 )
          i = j
          j = k
          k = m
          l = mod( 53*l+1 , 169 )
          if(mod(l*m,64) .ge. 32) s = s + t
          t = .5d0*t
  200   continue
        ranu(ii) = s
  300 continue
      ranc  =   362436d0 / 16777216d0
      rancd =  7654321d0 / 16777216d0
      rancm = 16777213d0 / 16777216d0
      iranmr = 97
      jranmr = 33
      end
