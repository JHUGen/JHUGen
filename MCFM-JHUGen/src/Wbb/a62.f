      double complex function a62(st,j1,j2,j3,j4,j5,j6,za,zb) 
      implicit none
************************************************************************
*     Author: R.K. Ellis                                               *
*     July, 1998.                                                      *
*     Updated 8/16/01: UV subtraction included                         *
************************************************************************
*     implementation of Eqs. (2.7) and (2.8) of BDKW hep-ph/9610370
*     with ns=0
*     character string st can take the value pp or pm
      include 'constants.f'
      include 'zprods_decl.f'
      integer j1,j2,j3,j4,j5,j6
      character*2 st
      double complex a6,aa6sf,aa6tp,aa6uv

c----Includes ultraviolet subtraction aa6uv (with +ve sign since
c--- it appears in interference terms,e.g. 1234 x 3214)
      call a6routine(st,j1,j2,j3,j4,j5,j6,za,zb,aa6sf,aa6tp,aa6uv) 
      a62=a6(st,j1,j2,j3,j4,j5,j6,za,zb)/xnsq
     & +(aa6sf*dble(nf)-aa6tp)/xn+aa6uv

      if (st .eq. 'pp') then
      a62=a62+a6('pm',j1,j3,j2,j4,j5,j6,za,zb)*(one+one/xnsq)
     &       -a6('sl',j2,j3,j1,j4,j5,j6,za,zb)/xnsq
      elseif (st .eq. 'pm') then
      a62=a62+a6('pp',j1,j3,j2,j4,j5,j6,za,zb)*(one+one/xnsq)
     &       +a6('sl',j3,j2,j1,j4,j5,j6,za,zb)/xnsq
      else
      write(6,*) 'Unimplemented st in a62',st 
      stop
      endif
      return
      end

