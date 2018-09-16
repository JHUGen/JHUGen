      subroutine ggZZmassamp_new(p,za,zb,mt,AmassLL,AmassLR)
      implicit none
      include 'types.f'
      
c--- Author: J. Campbell, September 2013
c---
c--- Given momentum p and spinor products za,zb, calculate contribution
c--- of gg->ZZ through a loop of fermions with mass "mt".
c--- Returns amplitudes AmassLL, AmassLR indexed by
c--- helicities (h1=gluon,h2=gluon,h34=Z(34),h56=Z(p56))
c--- with LL and LR couplings in loop
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'docheck.f'
      include 'ggZZcomputemp.f'
      real(dp):: mt,p(mxpart,4)
      complex(dp):: bub(2,2,2,2,-2:0),box(2,2,2,2,-2:0),
     & tri(2,2,2,2,-2:0),drat(2,2,2,2,3),totrat(2,2,2,2),
     & A(6),Xmp(2,2),Xpp(2,2),Xpm(2,2),Xmm(2,2),
     & AmassLL(2,2,2,2),AmassLR(2,2,2,2)
      
c--- initialize all scalar integrals
      call ZZintegraleval(p,mt)   

c--- First calculate ALR
      call Acalc(1,2,3,4,5,6,mt,A)
      call LRcalc(1,2,3,4,5,6,za,zb,A,xpp,xmp,xpm,xmm)
      AmassLR(1,2,:,:)=Xmp(:,:)
      AmassLR(2,2,:,:)=Xpp(:,:)
      AmassLR(1,1,:,:)=Xmm(:,:)
      AmassLR(2,1,:,:)=Xpm(:,:)

c--- Second calculate ALL            
      computemp=.false.  ! flag to not compute (-,+) box and triangle coeffs
c--- boxes      
      call ZZmassivebox(1,2,3,4,5,6,za,zb,mt,box,drat)
c--- bubbles
      call ZZmassivebub(1,2,3,4,5,6,za,zb,mt,bub,totrat)
c--- triangles -  note: calculated last to allow easy calculation
c--- of three-mass triangle from rational and drat (passed in)
      call ZZmassivetri(1,2,3,4,5,6,za,zb,mt,totrat,drat,tri)
      AmassLL(:,:,:,:)=box(:,:,:,:,0)+tri(:,:,:,:,0)+bub(:,:,:,:,0)

c--- new implementation using 6-d boxes for (-,+) and (+,-) amplitudes
      call ZZmassiveboxtri(1,2,3,4,5,6,za,zb,mt,totrat,box,tri)
      AmassLL(1,2,:,:)=box(1,2,:,:,0)+tri(1,2,:,:,0)+bub(1,2,:,:,0)
      AmassLL(2,1,:,:)=box(2,1,:,:,0)+tri(2,1,:,:,0)+bub(2,1,:,:,0)

      if (docheck) then
c-- check final result for LL
        call ggZZparsecheck('LL',-1,-1,AmassLL(1,1,:,:),'0+1,res/lo')
        call ggZZparsecheck('LL',-1,+1,AmassLL(1,2,:,:),'0+1,res/lo')
        call ggZZparsecheck('LL',+1,-1,AmassLL(2,1,:,:),'0+1,res/lo')
        call ggZZparsecheck('LL',+1,+1,AmassLL(2,2,:,:),'0+1,res/lo')
c-- check that final result for RR is same as for LL
        call ggZZparsecheck('RR',-1,-1,AmassLL(1,1,:,:),'0+1,res/lo')
        call ggZZparsecheck('RR',-1,+1,AmassLL(1,2,:,:),'0+1,res/lo')
        call ggZZparsecheck('RR',+1,-1,AmassLL(2,1,:,:),'0+1,res/lo')
        call ggZZparsecheck('RR',+1,+1,AmassLL(2,2,:,:),'0+1,res/lo')
c-- check final result for LR
        call ggZZparsecheck('LR',-1,-1,AmassLR(1,1,:,:),'0+1,res/lo')
        call ggZZparsecheck('LR',-1,+1,AmassLR(1,2,:,:),'0+1,res/lo')
        call ggZZparsecheck('LR',+1,-1,AmassLR(2,1,:,:),'0+1,res/lo')
        call ggZZparsecheck('LR',+1,+1,AmassLR(2,2,:,:),'0+1,res/lo')
c-- check that final result for RL is same as for LR
        call ggZZparsecheck('RL',-1,-1,AmassLR(1,1,:,:),'0+1,res/lo')
        call ggZZparsecheck('RL',-1,+1,AmassLR(1,2,:,:),'0+1,res/lo')
        call ggZZparsecheck('RL',+1,-1,AmassLR(2,1,:,:),'0+1,res/lo')
        call ggZZparsecheck('RL',+1,+1,AmassLR(2,2,:,:),'0+1,res/lo')
      endif

      return
      end
      
