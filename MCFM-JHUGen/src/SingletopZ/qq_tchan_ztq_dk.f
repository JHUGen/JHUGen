      subroutine qq_tchan_ztq_dk(p,msq)
c---Matrix element squared averaged over initial colors and spins
c     u(-p1)+b(p2)->e-(p3)+e+(p4)+t(nu(p5)+mu(p6)+b(p7))+d(p8)
      implicit none
      include 'constants.f'
      include 'alpha1.f'
      include 'nwz.f'
      include 'ewcouple.f'
      integer j
      double precision p(mxpart,4),msq(-nf:nf,-nf:nf),
     &     b_u,u_b,db_b,b_db,d_bb,ub_bb,bb_d,bb_ub
      double complex mdecaymb(2,2),mdecay(2)
      double complex ub_amp(2,2), bu_amp(2,2),dbb_amp(2,2),bdb_amp(2,2),
     & bbd_amp(2,2),ubbb_amp(2,2),bbub_amp(2,2)


c --  Deal with top decay                                                                                                                    
c --  tdecay assumes a massive b, hence 4 polarizations -- we only need one                                                                       mb=0d0
      if (nwz .eq. 1) then
      call tdecay(p,5,6,7,mdecaymb)
      mdecay(1:2)=mdecaymb(1,1:2)
      elseif (nwz .eq. -1) then
      call adecay(p,5,6,7,mdecaymb)
      mdecay(1)=mdecaymb(2,2)
      mdecay(2)=mdecaymb(2,1)
      endif

      msq(:,:)=0d0
      alpha1=0d0

      if (nwz .eq. 1) then

      call ubtzdamp_dk(p,1,2,3,4,8,ub_amp)
      call ubtzdamp_dk(p,2,1,3,4,8,bu_amp)
      call ubtzdamp_dk(p,8,2,3,4,1,dbb_amp)
      call ubtzdamp_dk(p,8,1,3,4,2,bdb_amp)

c --  multiply by decay matrix elements
      do j =1,2
         ub_amp(:,j)=ub_amp(:,j)*mdecay(j)
         bu_amp(:,j)=bu_amp(:,j)*mdecay(j)
         dbb_amp(:,j)=dbb_amp(:,j)*mdecay(j)
         bdb_amp(:,j)=bdb_amp(:,j)*mdecay(j)
      enddo
      
      u_b=cdabs(ub_amp(1,1))**2+cdabs(ub_amp(1,2))**2
     &    + cdabs(ub_amp(2,1))**2+cdabs(ub_amp(2,2))**2
      b_u=cdabs(bu_amp(1,1))**2+cdabs(bu_amp(1,2))**2
     &    + cdabs(bu_amp(2,1))**2+cdabs(bu_amp(2,2))**2
      db_b=cdabs(dbb_amp(1,1))**2+cdabs(dbb_amp(1,2))**2
     &    + cdabs(dbb_amp(2,1))**2+cdabs(dbb_amp(2,2))**2
      b_db=cdabs(bdb_amp(1,1))**2+cdabs(bdb_amp(1,2))**2
     &    + cdabs(bdb_amp(2,1))**2+cdabs(bdb_amp(2,2))**2

      u_b=4d0*xn**2*aveqq*gwsq**4*esq**2*u_b
      b_u=4d0*xn**2*aveqq*gwsq**4*esq**2*b_u
      db_b=4d0*xn**2*aveqq*gwsq**4*esq**2*db_b
      b_db=4d0*xn**2*aveqq*gwsq**4*esq**2*b_db

      msq(-1,+5)=db_b
      msq(-3,+5)=db_b
      msq(+2,+5)=u_b
      msq(+4,+5)=u_b
      msq(+5,-1)=b_db
      msq(+5,-3)=b_db
      msq(+5,+2)=b_u
      msq(+5,+4)=b_u

      elseif(nwz .eq. -1) then

      call ubtzdamp_dk(p,1,2,4,3,8,ubbb_amp)
      call ubtzdamp_dk(p,2,1,4,3,8,bbub_amp)
      call ubtzdamp_dk(p,8,2,4,3,1,dbb_amp)
      call ubtzdamp_dk(p,8,1,4,3,2,bbd_amp)

c --  multiply by decay matrix elements
      do j =1,2
         dbb_amp(:,j)=dbb_amp(:,j)*mdecay(j)
         bbd_amp(:,j)=bbd_amp(:,j)*mdecay(j)
         ubbb_amp(:,j)=ubbb_amp(:,j)*mdecay(j)
         bbub_amp(:,j)=bbub_amp(:,j)*mdecay(j)
      enddo
      
      d_bb=cdabs(dbb_amp(1,1))**2+cdabs(dbb_amp(1,2))**2
     &    + cdabs(dbb_amp(2,1))**2+cdabs(dbb_amp(2,2))**2
      bb_d=cdabs(bbd_amp(1,1))**2+cdabs(bbd_amp(1,2))**2
     &    + cdabs(bbd_amp(2,1))**2+cdabs(bbd_amp(2,2))**2
      ub_bb=cdabs(ubbb_amp(1,1))**2+cdabs(ubbb_amp(1,2))**2
     &    + cdabs(ubbb_amp(2,1))**2+cdabs(ubbb_amp(2,2))**2
      bb_ub=cdabs(bbub_amp(1,1))**2+cdabs(bbub_amp(1,2))**2
     &    + cdabs(bbub_amp(2,1))**2+cdabs(bbub_amp(2,2))**2

      d_bb= 4d0*xn**2*aveqq*gwsq**4*esq**2*d_bb
      bb_d= 4d0*xn**2*aveqq*gwsq**4*esq**2*bb_d
      ub_bb=4d0*xn**2*aveqq*gwsq**4*esq**2*ub_bb
      bb_ub=4d0*xn**2*aveqq*gwsq**4*esq**2*bb_ub

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
