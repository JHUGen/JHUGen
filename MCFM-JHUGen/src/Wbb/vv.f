      double complex function vv(st,j1,j2,j3,j4,j5,j6)
************************************************************************
*     Author: R.K. Ellis                                               *
*     July, 1998.                                                      *
************************************************************************
      include 'constants.f'
      include 'sprods_com.f'
      include 'scale.f'
      include 'epinv.f'
      integer j1,j2,j3,j4,j5,j6
      double complex Lnrat,xl12,xl34,xl23,xl56
      character*2 st


      xl12=Lnrat(musq,-s(j1,j2))
      xl34=Lnrat(musq,-s(j3,j4))

      if ((st .eq. 'pp') .or. (st .eq. 'pm')) then

c   Vpm =Vpp= 
c     - 2*epinv^2   + epinv* ( 2/3 - xlog12 - xlog34 )
c     + 10/9 - 1/2*xlog12^2 - 1/2*xlog34^2 + 13/6*xlog23 - 3/2*xlog56
c   Vsl = 
c    -2*epinv^2 + epinv* ( - 3 - xlog12 - xlog34 )
c    -15/2 - 3/2*xlog12 - 1/2*xlog12^2 - 3/2*xlog34 - 1/2*xlog34^2

      xl23=Lnrat(musq,-s(j2,j3))
      xl56=Lnrat(musq,-s(j5,j6))
      vv =-epinv*(two*epinv-two/three+xl12+xl34)
     . +10d0/9d0-0.5d0*(xl12**2+xl34**2)+13d0/6d0*xl23-1.5d0*xl56


      elseif (st .eq. 'sl') then
      vv=
     & -epinv*(two*epinv+three+xl12+xl34)
     & -7.5d0-1.5d0*(xl12+xl34)-0.5d0*(xl12**2+xl34**2)


      endif

      end
 
 
 

