      double complex function a61(st,j1,j2,j3,j4,j5,j6,za,zb) 
************************************************************************
*     Author: R.K. Ellis                                               *
*     July, 1998.                                                      *
************************************************************************
*     implementation of Eqs. (2.7) and (2.8) of BDKW hep-ph/9610370
*     with ns=0
*     character string st can take the value pp or pm
*     a61(pp,j1,j2,j3,j4,j5,j6,za,zb) corresponds to 
*     q(j1,+)+Q(j3,-)+e(j6,+)+q~(j4)+Q~(j2)+e~(j5)
*     a61(pm,j1,j2,j3,j4,j5,j6,za,zb) corresponds to 
*     q(j1,+)+Q(j3,+)+e(j6,+)+q~(j4)+Q~(j2)+e~(j5)
      implicit none
      include 'constants.f'
      include 'zprods_decl.f'
      integer j1,j2,j3,j4,j5,j6
      character*2 st
      double complex a6,aa6sf,aa6tp,aa6uv

c----Includes ultraviolet subtraction aa6uv
      call a6routine(st,j1,j2,j3,j4,j5,j6,za,zb,aa6sf,aa6tp,aa6uv) 
      a61=(one-two/xnsq)*a6(st,j1,j2,j3,j4,j5,j6,za,zb)
     & -(dble(nf)*aa6sf-aa6tp)/xn-aa6uv

      if (st .eq. 'pp') then
      a61=a61+(-two*a6('pm',j1,j3,j2,j4,j5,j6,za,zb)
     &             +a6('sl',j2,j3,j1,j4,j5,j6,za,zb))/xnsq
      elseif (st .eq. 'pm') then
      a61=a61+(-two*a6('pp',j1,j3,j2,j4,j5,j6,za,zb)
     &             -a6('sl',j3,j2,j1,j4,j5,j6,za,zb))/xnsq
      else
      write(6,*) 'Unimplemented st in a61',st 
      stop
      endif
      return
      end

