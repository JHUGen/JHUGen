      subroutine qq_WZqqstrong(pin,msq)
      implicit none
      include 'types.f'
c--- Author: J.M. Campbell, January 2015
c--- q(-p1)+q(-p2)->W(p3,p4)+Z(p5,p6)+q(p7)+q(p8);
c--- at O(alpha_em^4 alpha_s^2)
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'cmplxmass.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'masses.f'
      include 'zprods_decl.f'
      include 'first.f'
      include 'nwz.f'
      integer:: nmax,jmax
      parameter(jmax=8,nmax=7)
      integer:: j,k,
     & uqsq_dqsq,uqcq_dqcq,uqdq_dqdq,cqdq_dqsq,
     & uqcq_uqsq,uquq_dquq,uquq_uqdq
      parameter(
     & uqsq_dqsq=1,uqcq_dqcq=2,uqdq_dqdq=3,cqdq_dqsq=4,
     & uqcq_uqsq=5,uquq_dquq=6,uquq_uqdq=7)
      integer:: h1,h2,h3,h5
      real(dp):: p(mxpart,4),msq(fn:nf,fn:nf),temp(fn:nf,fn:nf),
     & tempb(fn:nf,fn:nf),stat,spinavge,Colorfac,pin(mxpart,4),msqgg
      complex(dp)::zab(mxpart,4,mxpart),zba(mxpart,4,mxpart),
     & amp(nmax,2,2,2,2),ampa(nmax,2,2,2,2),ampb(nmax,2,2,2,2),cdotpr,
     & k7341(4),k7561(4),k7342(4),k7562(4),
     & tmpz1(mxpart,mxpart),tmpz2(mxpart,4,mxpart)
      complex(dp):: 
     & jw7_34_1g(2,4),jw7_34_2g(2,4),
     & jw8_34_1g(2,4),jw8_34_2g(2,4),
     & jz7_56_1g(4,2,2,2),jz7_56_2g(4,2,2,2),
     & jz8_56_1g(4,2,2,2),jz8_56_2g(4,2,2,2),

     & ampZ2_1728(2,2,2),ampZ2_1827(2,2,2),
     & ampZ2_2817(2,2,2),ampZ2_2718(2,2,2),
     & s7341,s7561,s7562,s7342
      parameter(spinavge=0.25_dp,stat=0.5_dp,Colorfac=V/four/xn**2)
      integer,parameter::j1(jmax)=(/1,2,8,7,7,2,7,1/)
      integer,parameter::j2(jmax)=(/2,1,7,8,2,7,1,7/)
      integer,parameter::j7(jmax)=(/7,7,1,1,1,8,2,8/)
      integer,parameter::j8(jmax)=(/8,8,2,2,8,1,8,2/)
      msq(:,:)=zip

c--- This calculation uses the complex-mass scheme (c.f. arXiv:hep-ph/0605312)
c--- and the following lines set up the appropriate masses and sin^2(theta_w)
      if (first) then
       cwmass2=cplx2(wmass**2,-wmass*wwidth)
       czmass2=cplx2(zmass**2,-zmass*zwidth)
       cxw=cone-cwmass2/czmass2
c       cxw=cplx2(xw,zip) ! DEBUG: Madgraph comparison
       write(6,*)
       write(6,*) '**************** Complex-mass scheme ***************'
       write(6,*) '*                                                  *'
       write(6,77) cwmass2
       write(6,78) czmass2
       write(6,79) cxw
       write(6,*) '*                                                  *'
       write(6,*) '****************************************************'
       write(6,*)
       first=.false.
       call flush(6)
      endif

      if (nwz == +1) then
        p(:,:)=pin(:,:)
      else
c--- W- amplitudes obtained by parity transformation on W+
        p(1:2,:)=pin(1:2,:)
        p(7,:)=pin(8,:)
        p(8,:)=pin(7,:)
        p(3,:)=pin(4,:)
        p(4,:)=pin(3,:)
        p(5,:)=pin(6,:)
        p(6,:)=pin(5,:)
      endif

C---setup spinors and spinorvector products
      call spinorcurr(8,p,za,zb,zab,zba)
      if (nwz == -1) then
c--- parity transformation corresponds to complex conjugation
        tmpz1(1:8,1:8)=za(1:8,1:8)
        za(1:8,1:8)=zb(1:8,1:8)
        zb(1:8,1:8)=tmpz1(1:8,1:8)
        tmpz2(1:8,:,1:8)=zab(1:8,:,1:8)
        zab(1:8,:,1:8)=zba(1:8,:,1:8)
        zba(1:8,:,1:8)=tmpz2(1:8,:,1:8)
      endif

      do j=1,jmax
      temp(:,:)=zip
      tempb(:,:)=zip
      amp(:,:,:,:,:)=czip
      ampa(:,:,:,:,:)=czip
      ampb(:,:,:,:,:)=czip

c--- jz7_34_1g:  current for gluon -> e^- e^+ qq
      call jonezstrong(j7(j),5,6,j1(j),za,zb,zab,zba,jz7_56_1g)
      call jonezstrong(j7(j),5,6,j2(j),za,zb,zab,zba,jz7_56_2g)
      call jonezstrong(j8(j),5,6,j1(j),za,zb,zab,zba,jz8_56_1g)
      call jonezstrong(j8(j),5,6,j2(j),za,zb,zab,zba,jz8_56_2g)

c--- jw7_34_1g:  current for gluon -> e^- nu~_e qq
      call jonewstrong(j7(j),3,4,j1(j),za,zb,zab,jw7_34_1g)
      call jonewstrong(j8(j),3,4,j1(j),za,zb,zab,jw8_34_1g)
      call jonewstrong(j8(j),3,4,j2(j),za,zb,zab,jw8_34_2g)
      call jonewstrong(j7(j),3,4,j2(j),za,zb,zab,jw7_34_2g)

c---- Z2 contributions
      call jZexch2strong
     & (j1(j),j2(j),3,4,5,6,j7(j),j8(j),za,zb,ampZ2_1728)
      call jZexch2strong
     & (j1(j),j2(j),3,4,5,6,j8(j),j7(j),za,zb,ampZ2_1827)
      call jZexch2strong
     & (j2(j),j1(j),3,4,5,6,j8(j),j7(j),za,zb,ampZ2_2817)
      call jZexch2strong
     & (j2(j),j1(j),3,4,5,6,j7(j),j8(j),za,zb,ampZ2_2718)
      
      k7341(:)=half*(zab(j1(j),:,j1(j))+zab(3,:,3)
     & +zab(4,:,4)+zab(j7(j),:,j7(j)))
      k7561(:)=half*(zab(j1(j),:,j1(j))+zab(5,:,5)
     & +zab(6,:,6)+zab(j7(j),:,j7(j)))
      k7342(:)=half*(zab(j2(j),:,j2(j))+zab(3,:,3)
     & +zab(4,:,4)+zab(j7(j),:,j7(j)))
      k7562(:)=half*(zab(j2(j),:,j2(j))+zab(5,:,5)
     & +zab(6,:,6)+zab(j7(j),:,j7(j)))
      s7341=cdotpr(k7341,k7341)
      s7561=cdotpr(k7561,k7561)
      s7562=cdotpr(k7562,k7562)
      s7342=cdotpr(k7342,k7342)
            
C-----setup for (uqsq_dqsq) (2,3)->(1,3)
      do h2=1,2
      do h5=1,2
      amp(uqsq_dqsq,1,h2,1,h5)=
     & +cdotpr(jw7_34_1g(2,:),jz8_56_2g(:,1,h2,h5))/s7341
      amp(uqsq_dqsq,1,h2,1,h5)=amp(uqsq_dqsq,1,h2,1,h5)
     & +ampZ2_1728(1,h2,h5)
      enddo
      enddo      

      do h1=1,2
      do h2=1,2
      do h3=1,2
      do h5=1,2
      temp(2,3)=temp(2,3)+esq**4*gsq**2*Colorfac*spinavge
     &   *real(amp(uqsq_dqsq,h1,h2,h3,h5)
     & *conjg(amp(uqsq_dqsq,h1,h2,h3,h5)))
      enddo
      enddo
      enddo
      enddo
      
c-----setup for (cqdq_dqsq)
      do h2=1,2
      do h5=1,2
      amp(cqdq_dqsq,1,h2,1,h5)=
     & +cdotpr(jw8_34_1g(2,:),jz7_56_2g(:,1,h2,h5))/s7562
      amp(cqdq_dqsq,1,h2,1,h5)=amp(cqdq_dqsq,1,h2,1,h5)
     & +ampZ2_1827(1,h2,h5)
      enddo
      enddo      

      do h1=1,2
      do h2=1,2
      do h3=1,2
      do h5=1,2
      temp(4,1)=temp(4,1)+esq**4*gsq**2*Colorfac*spinavge
     &   *real(amp(cqdq_dqsq,h1,h2,h3,h5)
     & *conjg(amp(cqdq_dqsq,h1,h2,h3,h5)))
      enddo
      enddo
      enddo
      enddo


C-----setup for (uqdq_dqdq) (2,1)->(1,1)
      ampa(uqdq_dqdq,:,:,:,:)=amp(uqsq_dqsq,:,:,:,:)
      ampb(uqdq_dqdq,:,:,:,:)=amp(cqdq_dqsq,:,:,:,:)
      
      do h1=1,2
      do h2=1,2
      do h3=1,2
      do h5=1,2
      temp(2,1)=temp(2,1)+esq**4*gsq**2*Colorfac*spinavge
     &   *real(ampa(uqdq_dqdq,h1,h2,h3,h5)
     & *conjg(ampa(uqdq_dqdq,h1,h2,h3,h5)))
      temp(2,1)=temp(2,1)+esq**4*gsq**2*Colorfac*spinavge
     &   *real(ampb(uqdq_dqdq,h1,h2,h3,h5)
     & *conjg(ampb(uqdq_dqdq,h1,h2,h3,h5)))
      if (h1 == h2) then
      temp(2,1)=temp(2,1)+two/xn*esq**4*gsq**2*Colorfac*spinavge
     &   *real(ampa(uqdq_dqdq,h1,h2,h3,h5)
     & *conjg(ampb(uqdq_dqdq,h1,h2,h3,h5)))
      endif
      enddo
      enddo
      enddo
      enddo


C-----setup for (uqcq_dqcq) (2,4)->(1,4)
      do h2=1,2
      do h5=1,2
      amp(uqcq_dqcq,1,h2,1,h5)=
     & +cdotpr(jw7_34_1g(2,:),jz8_56_2g(:,2,h2,h5))/s7341
      amp(uqcq_dqcq,1,h2,1,h5)=amp(uqcq_dqcq,1,h2,1,h5)
     & +ampZ2_1728(2,h2,h5)
      enddo
      enddo
      
      do h1=1,2
      do h2=1,2
      do h3=1,2
      do h5=1,2
      temp(2,4)=temp(2,4)+esq**4*gsq**2*Colorfac*spinavge
     &   *real(amp(uqcq_dqcq,h1,h2,h3,h5)
     & *conjg(amp(uqcq_dqcq,h1,h2,h3,h5)))
      enddo
      enddo
      enddo
      enddo
      
c      write(6,*) 'temp(2,4)',temp(2,4)
c      stop
      
C-----setup for (uqcq_uqsq) (2,4)->(2,3)
      do h2=1,2
      do h5=1,2
      amp(uqcq_uqsq,1,h2,1,h5)=
     & +cdotpr(jw8_34_2g(2,:),jz7_56_1g(:,2,h2,h5))/s7561
      amp(uqcq_uqsq,1,h2,1,h5)=amp(uqcq_uqsq,1,h2,1,h5)
     & +ampZ2_2817(2,h2,h5)
      enddo
      enddo
      
      do h1=1,2
      do h2=1,2
      do h3=1,2
      do h5=1,2
      tempb(2,4)=tempb(2,4)+esq**4*gsq**2*Colorfac*spinavge
     &   *real(amp(uqcq_uqsq,h1,h2,h3,h5)
     & *conjg(amp(uqcq_uqsq,h1,h2,h3,h5)))
      enddo
      enddo
      enddo
      enddo
            
C-----setup for (uquq_dquq) (2,4)->(1,4)
      ampa(uquq_dquq,:,:,:,:)=amp(uqcq_dqcq,:,:,:,:)
      do h2=1,2
      do h5=1,2
      ampb(uquq_dquq,1,h2,1,h5)=
     & +cdotpr(jw7_34_2g(2,:),jz8_56_1g(:,2,h2,h5))/s7342
      ampb(uquq_dquq,1,h2,1,h5)=ampb(uquq_dquq,1,h2,1,h5)
     & +ampZ2_2718(2,h2,h5)
      enddo
      enddo
      
      do h1=1,2
      do h2=1,2
      do h3=1,2
      do h5=1,2
      temp(2,2)=temp(2,2)+esq**4*gsq**2*Colorfac*spinavge
     &   *real(ampa(uquq_dquq,h1,h2,h3,h5)
     & *conjg(ampa(uquq_dquq,h1,h2,h3,h5)))
      temp(2,2)=temp(2,2)+esq**4*gsq**2*Colorfac*spinavge
     &   *real(ampb(uquq_dquq,h1,h2,h3,h5)
     & *conjg(ampb(uquq_dquq,h1,h2,h3,h5)))
      if (h1 == h2) then
      temp(2,2)=temp(2,2)+two/xn*esq**4*gsq**2*Colorfac*spinavge
     &   *real(ampa(uquq_dquq,h1,h2,h3,h5)
     & *conjg(ampb(uquq_dquq,h1,h2,h3,h5)))
      endif
      enddo
      enddo
      enddo
      enddo
      
C-----setup for (uquq_uqdq) (2,2)->(2,1) :: not canonical order, but useful for crossings
      ampa(uquq_uqdq,:,:,:,:)=amp(uqcq_uqsq,:,:,:,:)
      do h2=1,2
      do h5=1,2
      ampb(uquq_uqdq,1,h2,1,h5)=
     & +cdotpr(jw8_34_1g(2,:),jz7_56_2g(:,2,h2,h5))/s7562
      ampb(uquq_uqdq,1,h2,1,h5)=ampb(uquq_uqdq,1,h2,1,h5)
     & +ampZ2_1827(2,h2,h5)
      enddo
      enddo
      
      do h1=1,2
      do h2=1,2
      do h3=1,2
      do h5=1,2
      tempb(2,2)=tempb(2,2)+esq**4*gsq**2*Colorfac*spinavge
     &   *real(ampa(uquq_uqdq,h1,h2,h3,h5)
     & *conjg(ampa(uquq_uqdq,h1,h2,h3,h5)))
      tempb(2,2)=tempb(2,2)+esq**4*gsq**2*Colorfac*spinavge
     &   *real(ampb(uquq_uqdq,h1,h2,h3,h5)
     & *conjg(ampb(uquq_uqdq,h1,h2,h3,h5)))
      if (h1 == h2) then
      tempb(2,2)=tempb(2,2)+two/xn*esq**4*gsq**2*Colorfac*spinavge
     &   *real(ampa(uquq_uqdq,h1,h2,h3,h5)
     & *conjg(ampb(uquq_uqdq,h1,h2,h3,h5)))
      endif
      enddo
      enddo
      enddo
      enddo

c--- setup for (cquq_squq) (4,2)->(3,2) :: not canonical order, but useful for crossings
      do h1=1,2
      do h2=1,2
      do h3=1,2
      do h5=1,2
      tempb(4,2)=tempb(4,2)+esq**4*gsq**2*Colorfac*spinavge
     &   *real(ampb(uquq_uqdq,h1,h2,h3,h5)
     & *conjg(ampb(uquq_uqdq,h1,h2,h3,h5)))
      enddo
      enddo
      enddo
      enddo
      
      if (j==1) then
      msq(2,3)=temp(2,3)
      msq(4,1)=temp(4,1)
      msq(2,1)=temp(2,1)*half
      msq(4,3)=msq(2,1)
      msq(2,4)=temp(2,4)+tempb(2,4)
      msq(2,2)=temp(2,2)
      msq(4,4)=temp(2,2)

      elseif (j==2) then
      msq(3,2)=temp(2,3)
      msq(1,4)=temp(4,1)
      msq(1,2)=temp(2,1)*half
      msq(3,4)=msq(1,2)
      msq(4,2)=temp(2,4)+tempb(2,4)

c--- qbar-qbar
      elseif (j==3) then
      msq(-1,-4)=temp(2,4)
      msq(-1,-3)=temp(2,3)
      msq(-2,-3)=tempb(2,4)
      msq(-3,-1)=temp(4,1)
      msq(-2,-1)=tempb(2,2)*half
      msq(-4,-3)=tempb(2,2)*half

c--- qbar-qbar
      elseif (j==4) then
      msq(-1,-2)=temp(2,2)*half
      msq(-3,-4)=temp(2,2)*half
      msq(-3,-2)=temp(2,4)
      msq(-1,-1)=temp(2,1)
      msq(-3,-3)=temp(2,1)
      msq(-4,-1)=tempb(2,4)
      msq(-1,-3)=msq(-1,-3)+temp(4,1)
      msq(-3,-1)=msq(-3,-1)+temp(2,3)
      
c--- qbar-q
      elseif (j==5) then
      msq(-3,1)=temp(2,3)
      msq(-3,3)=temp(2,1)
      msq(-1,1)=temp(2,1)
      msq(-1,3)=temp(2,3)

c--- qbar-q
      elseif (j==6) then
      msq(-4,2)=temp(2,4)
      msq(-4,4)=temp(2,2)
      msq(-3,2)=temp(2,3)+tempb(2,4)
      msq(-1,4)=msq(-3,2)
      msq(-3,4)=temp(2,1)+tempb(2,2)+temp(4,1)+tempb(4,2)
      msq(-1,2)=msq(-3,4)
      msq(-2,2)=temp(2,2)
      msq(-2,4)=temp(2,4)

c--- q-qbar
      elseif (j==7) then
      msq(1,-3)=temp(2,3)
      msq(1,-1)=temp(2,1)
      msq(2,-4)=tempb(2,4)
      msq(3,-3)=temp(2,1)
      msq(3,-1)=temp(2,3)
      msq(4,-4)=tempb(2,2)

c--- q-qbar
      elseif (j==8) then
      msq(2,-3)=tempb(2,4)+temp(2,3)
      msq(2,-2)=temp(2,2)
      msq(2,-1)=tempb(4,2)+temp(4,1)+tempb(2,2)+temp(2,1)
      msq(4,-3)=msq(2,-1)
      msq(4,-2)=temp(2,4)
      msq(4,-1)=tempb(2,4)+temp(2,3)
      
c      write(6,*) 'temp(2,3)',temp(2,3)
c      write(6,*) 'temp(4,1)',temp(4,1)
c      write(6,*) 'temp(2,2)',temp(2,2)
c      write(6,*) 'temp(2,1)',temp(2,1)
c      write(6,*) 'temp(2,4)',temp(2,4)
c      write(6,*) 'tempb(2,4)',tempb(2,4)
c      write(6,*) 'tempb(2,2)',tempb(2,2)
c      write(6,*) 'tempb(4,2)',tempb(4,2)
c      write(6,*)

      endif

      enddo

c--- 2-gluon amplitudes
      call qqWZggampf(1,2,3,4,5,6,7,8,3,3,za,zb,msqgg)
      msq(2,-1)=msq(2,-1)+stat*aveqq*msqgg
      msq(4,-3)=msq(4,-3)+stat*aveqq*msqgg

      call qqWZggampf(2,1,3,4,5,6,7,8,3,4,za,zb,msqgg)
      msq(-1,2)=msq(-1,2)+stat*aveqq*msqgg
      msq(-3,4)=msq(-3,4)+stat*aveqq*msqgg

      call qqWZggampf(7,8,3,4,5,6,1,2,3,4,za,zb,msqgg)
      msq(0,0)=msq(0,0)+avegg*two*msqgg

      call qqWZggampf(7,2,3,4,5,6,1,8,3,4,za,zb,msqgg)
      msq(0,-1)=msq(0,-1)+aveqg*msqgg
      msq(0,-3)=msq(0,-3)+aveqg*msqgg

      call qqWZggampf(7,1,3,4,5,6,2,8,3,4,za,zb,msqgg)
      msq(-1,0)=msq(-1,0)+aveqg*msqgg
      msq(-3,0)=msq(-3,0)+aveqg*msqgg
      
      call qqWZggampf(2,8,3,4,5,6,1,7,3,4,za,zb,msqgg)
      msq(0,2)=msq(0,2)+aveqg*msqgg
      msq(0,4)=msq(0,4)+aveqg*msqgg

      call qqWZggampf(1,8,3,4,5,6,2,7,3,4,za,zb,msqgg)
      msq(2,0)=msq(2,0)+aveqg*msqgg
      msq(4,0)=msq(4,0)+aveqg*msqgg

c--- swap elements to correct places for nwz=-1
      if (nwz .eq .-1) then
        temp(:,:)=msq(:,:)
        do j=-nf,nf
        do k=-nf,nf
        msq(j,k)=temp(-j,-k)
        enddo
        enddo
      endif
      
      return

   77 format(' *      W-mass^2     (',f11.5,',',f11.5,')      *')
   78 format(' *      Z-mass^2     (',f11.5,',',f11.5,')      *')
   79 format(' *  sin^2(theta_w)   (',f11.5,',',f11.5,')      *')

      end
      
