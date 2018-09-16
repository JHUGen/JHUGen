      subroutine qq_tchan_ztq(p,msq)
      implicit none
      include 'types.f'
c---Matrix element squared averaged over initial colors and spins
c     u(-p1)+b(p2)->e-(p3)+e+(p4)+t(p5)+d(p6)
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'alpha1.f'
      include 'nwz.f'
      real(dp):: p(mxpart,4),msq(-nf:nf,-nf:nf),
     & b_u,u_b,db_b,b_db,d_bb,ub_bb,bb_d,bb_ub,ubtzdsq

      msq(:,:)=0d0

      alpha1=0d0

      if (nwz == 1) then
        u_b=ubtzdsq(p,1,2,3,4,6)
        b_u=ubtzdsq(p,2,1,3,4,6)
        db_b=ubtzdsq(p,6,2,3,4,1)
        b_db=ubtzdsq(p,6,1,3,4,2)

        msq(-1,+5)=db_b
        msq(-3,+5)=db_b
        msq(+2,+5)=u_b
        msq(+4,+5)=u_b
        msq(+5,-1)=b_db
        msq(+5,-3)=b_db
        msq(+5,+2)=b_u
        msq(+5,+4)=b_u
      elseif(nwz == -1) then
        ub_bb=ubtzdsq(p,1,2,4,3,6)
        bb_ub=ubtzdsq(p,2,1,4,3,6)
        d_bb=ubtzdsq(p,6,2,4,3,1)
        bb_d=ubtzdsq(p,6,1,4,3,2)

        msq(+1,-5)=d_bb
        msq(+3,-5)=d_bb
        msq(-2,-5)=ub_bb
        msq(-4,-5)=ub_bb

        msq(-5,+1)=bb_d
        msq(-5,+3)=bb_d
        msq(-5,-2)=bb_ub
        msq(-5,-4)=bb_ub      
      endif
      
      return
      end
