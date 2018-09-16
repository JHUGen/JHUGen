      subroutine gluegen(e1,p2,p3,e4,za,zb,zab,zba,pppp,ppmp,lo23,lo32,
     & calculateboth)
      implicit none
      include 'types.f'

C-----Authors: John Campbell and Keith Ellis, November 2011
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'zabprods_decl.f'
      complex(dp):: lo23(2,2),lo32(2,2),pppp,ppmp
      integer:: e1,p2,p3,e4,e1p,e4p
      logical:: calculateboth
C---we are calculating the amplitudes mxyp where x=1,2 and y=1,2
C---from our two base amplitudes which are pppp,ppmp

c--- perform the swaps 5 <-> 7 and 6 <-> 8
      e1p=12-e1
      e4p=14-e4

c--- JC: added minus signs here based on form of tree level amplitudes
c---     (complex conjugation = (-1) * [za <-> zb])

c-----mmmp
      lo23(1,1)=-pppp(e1,p2,p3,e4p,zb,za,zba,zab)
c-----mppp
      lo23(2,2)=+pppp(e1p,p2,p3,e4,za,zb,zab,zba)
c-----mpmp
      lo23(2,1)=+ppmp(e1p,p2,p3,e4,za,zb,zab,zba)
c-----mmpp
      lo23(1,2)=-ppmp(e1,p2,p3,e4p,zb,za,zba,zab)

c--- alternate return if calculating A43
      if (calculateboth .eqv. .false.) then
        lo32(1,1)=lo23(1,1)
        lo32(1,2)=lo23(1,2)
        lo32(2,1)=lo23(2,1)
        lo32(2,2)=lo23(2,2)
        return
      endif

c-----mmmp
      lo32(1,1)=-pppp(e1,p3,p2,e4p,zb,za,zba,zab)
c-----mppp
      lo32(2,2)=+pppp(e1p,p3,p2,e4,za,zb,zab,zba)
c-----mpmp
      lo32(1,2)=+ppmp(e1p,p3,p2,e4,za,zb,zab,zba)
c-----mmpp
      lo32(2,1)=-ppmp(e1,p3,p2,e4p,zb,za,zba,zab)

      return
      end
