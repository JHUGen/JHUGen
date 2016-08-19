*
* ADJUSTED TO REMOVE ALL 1/eps**2 PIECES IF EPINV2 IN EPINV2.F = 0d0
* ADJUSTED TO REMOVE ALL 1/eps**2 PIECES IF EPINV2 IN EPINV2.F = 0d0
* ADJUSTED TO REMOVE ALL 1/eps**2 PIECES IF EPINV2 IN EPINV2.F = 0d0
*
      double complex function vvg(st,j1,j2,j3,j4,j5,j6)
************************************************************************
*     Author: R.K. Ellis                                               *
*     July, 1998.                                                      *
************************************************************************
      include 'constants.f'
      include 'sprods_com.f'
      include 'epinv.f'
      include 'scale.f'
      include 'epinv2.f'
      integer j1,j2,j3,j4,j5,j6
      double complex Lnrat,xl12,xl34,xl23,xl56,Vcc,Vsc
      character*9 st


      xl12=Lnrat(musq,-s(j1,j2))
      xl34=Lnrat(musq,-s(j3,j4))
      xl23=Lnrat(musq,-s(j2,j3))
      xl56=Lnrat(musq,-s(j5,j6))

      if(st.eq.'q+g-g+qb-') then
      Vcc=-four
     . -(epinv*epinv2+xl12*epinv+half*xl12**2)
     . -(epinv*epinv2+xl23*epinv+half*xl23**2)
     . -(epinv*epinv2+xl34*epinv+half*xl34**2)-two*(epinv+xl56)
      Vsc=half*(one+epinv+xl56)

      elseif(st.eq.'q+g+g-qb-') then
      Vcc=-four
     .   -(epinv*epinv2+xl12*epinv+half*xl12**2)
     .   -(epinv*epinv2+xl23*epinv+half*xl23**2)
     .   -(epinv*epinv2+xl34*epinv+half*xl34**2)-two*(epinv+xl56)
      Vsc=half*(one+epinv+xl56)

      elseif(st.eq.'q+g+g+qb-') then
      Vcc=-four
     .  -(epinv*epinv2+xl12*epinv+half*xl12**2)
     .  -(epinv*epinv2+xl23*epinv+half*xl23**2)
     .  -(epinv*epinv2+xl34*epinv+half*xl34**2)-two*(epinv+xl56)
      Vsc=half*(one+epinv+xl56)

      elseif(st.eq.'q+g+qb-g-') then
      Vcc=-four
     . -(+(epinv*epinv2+xl12*epinv+half*xl12**2)
     .   +(epinv*epinv2+xl23*epinv+half*xl23**2))-two*(epinv+xl56)
      Vsc=half*(one+epinv+xl56)
 

      elseif(st.eq.'q+g+qb-g+') then
      Vcc=-four
     . -(epinv*epinv2+xl12*epinv+half*xl12**2)
     . -(epinv*epinv2+xl23*epinv+half*xl23**2)-two*(epinv+xl56)
      Vsc=half*(one+epinv+xl56)

      elseif(st.eq.'q+qb-g-g+') then 
      Vcc=-four-(epinv*epinv2+xl12*epinv+half*xl12**2)-two*(epinv+xl56)
      Vsc=half*(one+epinv+xl56)

      elseif(st.eq.'q+qb-g+g-') then
      Vcc=-four-(epinv*epinv2+xl12*epinv+half*xl12**2)-two*(epinv+xl56)
      Vsc=half*(one+epinv+xl56)

      elseif(st.eq.'q+qb-g+g+') then
      Vcc=-four-(epinv*epinv2+xl12*epinv+half*xl12**2)-two*(epinv+xl56)
      Vsc=half*(one+epinv+xl56)
      else 
      write(6,*) 'unimplemented st',st
      stop
      endif

      vvg=Vcc+Vsc
      end
