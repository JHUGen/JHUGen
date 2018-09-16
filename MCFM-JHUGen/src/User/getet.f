      function getet(E,px,py,pz)
      implicit none
      include 'types.f'
      real(dp):: getet
c--- given (E,px,py,pz) for a four-vector, calculates the corresponding
c--- Et or Pt, depending on the parameter that is set in mdata.f

      real(dp):: E,px,py,pz,etsq
      include 'useet.f'

      if (useEt) then
c--- this is the formula for Et
        etsq=px**2+py**2
        getet=sqrt(etsq)*E/sqrt(etsq+pz**2)
      else
c--- this is the formula for pt
        getet=sqrt(px**2+py**2)
      endif

      return
      end

