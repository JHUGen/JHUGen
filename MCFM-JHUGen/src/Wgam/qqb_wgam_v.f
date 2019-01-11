      subroutine qqb_wgam_v(p,msqv)
      implicit none
C----Author R.K.Ellis August 2002
C====Virtual corrections to
c     q(-p1)+qbar(-p2)-->e^-(p3)+nu(p4)+gamma(p5)
      include 'constants.f'
      include 'scheme.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'ckm.f'
      include 'zprods_decl.f'
      include 'nwz.f'
      integer j,k
      double precision msqv(-nf:nf,-nf:nf),p(mxpart,4),qbq,qqb
      double precision fac
      double complex agamtree,agamvirt

      call spinoru(5,p,za,zb)
    
     

      fac=ason2pi*cf*aveqq*2d0*xn*gwsq**2*esq
    
      scheme='dred'
      if (nwz .eq. -1) then
C ie ub-d
         
      qbq=+fac*2d0*dble(
     .Dconjg(agamtree(1,2,3,4,5,za,zb,+1))*agamvirt(1,2,3,4,5,za,zb,+1))
     .    +fac*2d0*dble(
     .Dconjg(agamtree(1,2,3,4,5,za,zb,-1))*agamvirt(1,2,3,4,5,za,zb,-1))
C ie d-ub
      qqb=+fac*2d0*dble(
     .Dconjg(agamtree(2,1,3,4,5,za,zb,+1))*agamvirt(2,1,3,4,5,za,zb,+1))
     .    +fac*2d0*dble(
     .Dconjg(agamtree(2,1,3,4,5,za,zb,-1))*agamvirt(2,1,3,4,5,za,zb,-1))

      elseif (nwz .eq. +1) then 
C ie db-u
          
      qbq=+fac*2d0*dble(
     .Dconjg(agamtree(2,1,4,3,5,zb,za,+1))*agamvirt(2,1,4,3,5,zb,za,+1))
     .    +fac*2d0*dble(
     .Dconjg(agamtree(2,1,4,3,5,zb,za,-1))*agamvirt(2,1,4,3,5,zb,za,-1))
C ie u-db
   
      qqb=+fac*2d0*dble(
     .Dconjg(agamtree(1,2,4,3,5,zb,za,+1))*agamvirt(1,2,4,3,5,zb,za,+1))
     .    +fac*2d0*dble(
     .Dconjg(agamtree(1,2,4,3,5,zb,za,-1))*agamvirt(1,2,4,3,5,zb,za,-1))
       
      endif

      do j=-nf,nf
      do k=-nf,nf
c--- set msqv=0 to initalize
      msqv(j,k)=0d0
          if ((j .gt. 0) .and. (k .lt. 0)) then
            msqv(j,k)=Vsq(j,k)*qqb
          elseif ((j .lt. 0) .and. (k .gt. 0)) then
            msqv(j,k)=Vsq(j,k)*qbq
          endif
      enddo
      enddo

      return
      end


      

      double complex function agamvirt(p1,p2,p3,p4,p5,za,zb,hgamma)
      implicit none
C----Author R.K.Ellis August 2002
c     q(-p1)+qbar(-p2)-->e^-(p3)+nu(p4)+gamma(p5)
      include 'constants.f'
      include 'zprods_decl.f'
      include 'masses.f'
      double complex agamtree,vpole,fagamma,fbgamma,prop,vpl
      double precision s34,s12
      integer p1,p2,p3,p4,p5,hgamma

      
      s34=dble(za(p3,p4)*zb(p4,p3))
      s12=dble(za(p1,p2)*zb(p2,p1))
      prop=s34/dcmplx(s34-wmass**2,wmass*wwidth)
      vpl=vpole(s12)
      if    (hgamma .eq. +1) then
        agamvirt=prop*(
     .  +Qd*fagamma(p1,p2,p3,p4,p5,za,zb)
     .  +Qu*fbgamma(p1,p2,p3,p4,p5,za,zb))
     .          +vpl*agamtree(p1,p2,p3,p4,p5,za,zb,+1)
      elseif (hgamma .eq. -1) then
        agamvirt=prop*(
     .  +Qu*fagamma(p2,p1,p4,p3,p5,zb,za)
     .  +Qd*fbgamma(p2,p1,p4,p3,p5,zb,za))
     .          +vpl*agamtree(p1,p2,p3,p4,p5,za,zb,-1)
      endif
      
      return
      end
