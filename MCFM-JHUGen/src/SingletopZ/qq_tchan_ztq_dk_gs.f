      subroutine qq_tchan_ztq_dk_gs(p,msq)
c---Matrix element SUBTRACTION squared averaged over initial colors and spins
c     u(-p1)+b(p2)->e-(p3)+e+(p4)+t(nu(p5)+e(p6)+b(p7))+d(p8)+g(p9)
      implicit none 
      include 'constants.f'
      include 'ptilde.f'
      include 'nwz.f'
      include 'qqgg.f'
      include 'breit.f'
      include 'masses.f'
      integer j,k,ib,nd
      double precision p(mxpart,4),msq(maxd,-nf:nf,-nf:nf)
      double precision 
     . msq17_2(-nf:nf,-nf:nf),msq27_1(-nf:nf,-nf:nf),
     . msq16_2(-nf:nf,-nf:nf),msq26_1(-nf:nf,-nf:nf),
     . msq17_6(-nf:nf,-nf:nf),msq27_5(-nf:nf,-nf:nf),
     . msq67_1(-nf:nf,-nf:nf),msq57_2(-nf:nf,-nf:nf),
     . msq27_6(-nf:nf,-nf:nf),msq17_5(-nf:nf,-nf:nf),
     . msq67_2(-nf:nf,-nf:nf),msq57_1(-nf:nf,-nf:nf),
     . sub17_2(4),sub27_1(4),sub16_2(4),sub26_1(4),
     . sub17_6(4),sub27_5(4),sub67_1(4),sub57_2(4),
     . sub17_5(4),sub27_6(4),sub67_2(4),sub57_1(4),
     . dummyv(-nf:nf,-nf:nf),dsubv
      external qq_tchan_ztq_dk,donothing_gvec

      ndmax=12

c---- calculate both initial-initial dipoles
c---- note that we do not require the gg dipoles, so the v-type
c---- entries are left as dummies

      qqproc=.true.      

c--- mass2 not used elsewhere, so value does not need to be saved
      
c--- NOTE: to pass to the "decay" (dk) version of this file, the labels
c---       of the dipoles have not been changed.
c---       For the calls to dips_mass, the momenta
c---       to which they refer will be transformed according to:
c---       5 -> (5+6+7), 6 -> 8, 7 -> 9

      mass2=0d0      
      call dips_mass(1,p,1,7,6,sub17_6,dsubv,msq17_6,dummyv,
     . qq_tchan_ztq_dk,donothing_gvec)
      call dips_mass(2,p,6,7,1,sub67_1,dsubv,msq67_1,dummyv,
     . qq_tchan_ztq_dk,donothing_gvec)

      mass2=mt
      call dips_mass(3,p,2,7,5,sub27_5,dsubv,msq27_5,dummyv,
     . qq_tchan_ztq_dk,donothing_gvec)
      call dips_mass(4,p,5,7,2,sub57_2,dsubv,msq57_2,dummyv,
     . qq_tchan_ztq_dk,donothing_gvec)

      mass2=0d0
      call dips_mass(5,p,2,7,6,sub27_6,dsubv,msq27_6,dummyv,
     . qq_tchan_ztq_dk,donothing_gvec)
      call dips_mass(6,p,6,7,2,sub67_2,dsubv,msq67_2,dummyv,
     . qq_tchan_ztq_dk,donothing_gvec)
     
      mass2=mt
      call dips_mass(7,p,1,7,5,sub17_5,dsubv,msq17_5,dummyv,
     . qq_tchan_ztq_dk,donothing_gvec)
      call dips_mass(8,p,5,7,1,sub57_1,dsubv,msq57_1,dummyv,
     . qq_tchan_ztq_dk,donothing_gvec)

      mass2=0d0
      call dips_mass(9,p,1,7,2,sub17_2,dsubv,msq17_2,dummyv,
     . qq_tchan_ztq_dk,donothing_gvec)
      call dips_mass(10,p,2,7,1,sub27_1,dsubv,msq27_1,dummyv,
     . qq_tchan_ztq_dk,donothing_gvec)
      
      call dips_mass(11,p,1,6,2,sub16_2,dsubv,msq16_2,dummyv,
     . qq_tchan_ztq_dk,donothing_gvec)
      call dips_mass(12,p,2,6,1,sub26_1,dsubv,msq26_1,dummyv,
     . qq_tchan_ztq_dk,donothing_gvec)

c--- for nwz=+1, initial state is b, for nwz=-1 it is b~
      ib=5*nwz

      do j=-nf,nf
      do k=-nf,nf

      do nd=1,ndmax
      msq(nd,j,k)=0d0
      enddo

      if  ((j .eq. 0) .and. (k .eq. 0)) then
         goto 21
      elseif ((j .eq. ib) .and. (k .ne. 0)) then
c--- checked 17,27,67 singularities, soft singularity looks ok      
        msq(5,j,k)=2d0*cf*sub27_6(qq)*msq27_6(j,k)
        msq(6,j,k)=2d0*cf*sub67_2(qq)*msq67_2(j,k)
        msq(7,j,k)=2d0*cf*sub17_5(qq)*msq17_5(j,k)
        msq(8,j,k)=2d0*cf*sub57_1(qq)*msq57_1(j,k)
      elseif ((j .ne. 0) .and. (k .eq. ib)) then
c--- checked 17,27,67 singularities, soft singularity looks ok      
        msq(1,j,k)=2d0*cf*sub17_6(qq)*msq17_6(j,k)
        msq(2,j,k)=2d0*cf*sub67_1(qq)*msq67_1(j,k)
        msq(3,j,k)=2d0*cf*sub27_5(qq)*msq27_5(j,k)
        msq(4,j,k)=2d0*cf*sub57_2(qq)*msq57_2(j,k)
c--- checked
      elseif ((j .eq. ib) .and. (k .eq. 0)) then
         if (nwz .eq. +1) then
         msq(10,j,k)=2d0*tr*sub27_1(qg)*(msq27_1(j,+1)+msq27_1(j,+2)
     &   +msq27_1(j,+3)+msq27_1(j,+4)+msq27_1(j,+5))
         msq(12,j,k)=2d0*tr*sub26_1(qg)*(msq26_1(j,-1)+msq26_1(j,-2)
     &   +msq26_1(j,-3)+msq26_1(j,-4)+msq26_1(j,-5))
         else
         msq(10,j,k)=2d0*tr*sub27_1(qg)*(msq27_1(j,-1)+msq27_1(j,-2)
     &   +msq27_1(j,-3)+msq27_1(j,-4)+msq27_1(j,-5))
         msq(12,j,k)=2d0*tr*sub26_1(qg)*(msq26_1(j,+1)+msq26_1(j,+2)
     &   +msq26_1(j,+3)+msq26_1(j,+4)+msq26_1(j,+5))
         endif
      elseif ((j .eq. 0) .and. (k .eq. ib)) then
         if (nwz .eq. +1) then
         msq(9,j,k)=2d0*tr*sub17_2(qg)*(msq17_2(+1,k)+msq17_2(+2,k)
     &   +msq17_2(+3,k)+msq17_2(+4,k)+msq17_2(+5,k))
         msq(11,j,k)=2d0*tr*sub16_2(qg)*(msq16_2(-1,k)+msq16_2(-2,k)
     &   +msq16_2(-3,k)+msq16_2(-4,k)+msq16_2(-5,k))
         else
         msq(9,j,k)=2d0*tr*sub17_2(qg)*(msq17_2(-1,k)+msq17_2(-2,k)
     &   +msq17_2(-3,k)+msq17_2(-4,k)+msq17_2(-5,k))
         msq(11,j,k)=2d0*tr*sub16_2(qg)*(msq16_2(+1,k)+msq16_2(+2,k)
     &   +msq16_2(+3,k)+msq16_2(+4,k)+msq16_2(+5,k))
         endif
      elseif ((j .ne. 0) .and. (k .eq. 0)) then
         msq(10,j,k)=2d0*tr*sub27_1(qg)*msq27_1(j,ib)
      elseif ((j .eq. 0) .and. (k .ne. 0)) then
         msq(9,j,k)=2d0*tr*sub17_2(qg)*msq17_2(ib,k)
      endif
      
 21   continue

      enddo
      enddo

      return      
      end


