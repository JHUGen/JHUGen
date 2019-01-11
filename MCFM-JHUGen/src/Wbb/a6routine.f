      subroutine a6routine(st,j1,j2,j3,j4,j5,j6,za,zb,a6sf,a6tp,a6uv) 
************************************************************************
*     Author: R.K. Ellis                                               *
*     July, 1998.                                                      *
************************************************************************
*     a6sf is the sum of a6s and a6f 
*     it is only the sum which is needed.
      implicit none
      include 'constants.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'scale.f'
      include 'masses.f'
      include 'epinv.f'
      integer j1,j2,j3,j4,j5,j6
      double complex atree,virtsf,virtuv,virttp,Lnrat,a6sf,a6tp,a6uv,
     & tree
      character*2 st 
      logical msbar
      common/msbar/msbar
      
      if (st .eq. 'sl') then
      write(6,*) 'error in a6routine',st
      stop
      endif
      tree=atree(st,j1,j2,j3,j4,j5,j6,za,zb)
      if (mt .ne. 0d0) then
      virttp=-2d0/15d0*s(j2,j3)/mt**2
      else
      stop 'mt=0 in a6routine'
      endif
      virtsf=two/three*epinv
     . +two/three*Lnrat(musq,-s(j2,j3))+10d0/9d0
      virtuv=(epinv*(11d0-two/xn*dble(nf))-one)/three
      if (msbar) virtuv=virtuv+two*CF/xn
c---virtuv is the infinite and finite renormalization 
c---to get us to MS bar system
c--the term commented is associated with the number of legs
c  and in this program is taken care of in factorization procedure.
      a6sf=tree*virtsf 
      a6tp=tree*virttp 
      a6uv=tree*virtuv 
      return
      end

 
