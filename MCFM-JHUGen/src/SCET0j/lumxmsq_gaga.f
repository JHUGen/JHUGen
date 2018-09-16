      subroutine lumxmsq_gaga(p,xx,z1,z2,QB,order,xmsq)
      implicit none
      include 'types.f'
!====== C. Williams July 2015
c----Matrix element for Ga Ga production
C----averaged over initial colours and spins
c===== based upon similar routines for W and Z
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'scale.f'
      include 'facscale.f'
      include 'ewcharge.f'
      include 'scet_const.f'
      integer:: j,k,ih1,ih2,m,n,order
      real(dp),intent(out):: xmsq
      real(dp):: p(mxpart,4),s,fac,qqb,qbq,
     & xx(2),soft1(-1:1),soft2(-1:3),hard(2),
     & beama0(-5:5),beamb0(-5:5),
     & beama1(-5:5,-1:1),beamb1(-5:5,-1:1),
     & beama2(-5:5,-1:3),beamb2(-5:5,-1:3),
     & z1,z2,QB(2),lum0,lum1(-1:1),lum2(-1:3),bit,
     & msq(-nf:nf,-nf:nf),assemble,Haa(2)
      common/density/ih1,ih2
      real(dp)::statfac
      parameter(statfac=0.5_dp)
      real(dp):: msqt(-nf:nf,-nf:nf)
      real(dp):: qsum,tlrem,bitrem
      real(dp):: gggaga,msqgggaga,facgg
      
      s(j,k)=two*(p(j,4)*p(k,4)-p(j,1)*p(k,1)
     &           -p(j,2)*p(k,2)-p(j,3)*p(k,3))

      
      fac=xn*esq**2*statfac
      gggaga=zip
      qsum=zip
      do j=1,nf
         qsum=qsum+Q(j)**2
      enddo
!======= compute Matrix elements 
!      write(6,*) 'old routine'
!      call gamgamampsq(order,p,1,2,3,4,qqb,hard,tlrem)
!      write(6,*) 'hard(1),hard(2)',hard(1),hard(2)
!      write(6,*) 'new routine'
      call gamgamampsq_new(order,p,1,2,3,4,qqb,hard,tlrem)
!      write(6,*) 'hard(1),hard(2)',hard(1),hard(2)
!      pause
      
      qqb=fac*aveqq*qqb
      qbq=qqb

     
      call softqqbis(order,soft1,soft2)

      if (order >= 0) then
      call fdist(ih1,xx(1),facscale,beama0)
      call fdist(ih2,xx(2),facscale,beamb0)
      endif
      if (order >= 1) then
      call xbeam1bis(ih1,z1,xx(1),QB(1),beama1)
      call xbeam1bis(ih2,z2,xx(2),QB(2),beamb1)
      endif
      if (order >= 2) then
      call xbeam2bis(ih1,z1,xx(1),QB(1),beama2)
      call xbeam2bis(ih2,z2,xx(2),QB(2),beamb2)
      facgg=4._dp*esq*gsq/(16._dp*pisq)*Qsum
      gggaga=avegg*V*facgg**2*msqgggaga(s(1,2),s(1,3),s(2,3))*statfac
!      gggaga=zip
      endif

      xmsq=zip
      do j=-nf,nf
         k=-j
      
         if (j*k > 0) cycle    ! skip, qq, aa
      
         bit=assemble(order,
     &        beama0(j),beamb0(k),beama1(j,:),beamb1(k,:),
     &        beama2(j,:),beamb2(k,:),soft1,soft2,hard)
        

         bitrem=fac*aveqq*tlrem*beama0(j)*beamb0(k)*ason2pi**2
      
         if ((j > 0) .and. (k < 0)) then
            bit=bit*qqb*Q(j)**4
            if(order > 1) then
!======add on two-loop A functions which go like sum over quark charges
               bit=bit+bitrem*qsum*Q(j)**2
            endif
         elseif ((j < 0) .and. (k > 0)) then
            bit=bit*qbq*Q(k)**4
            if(order > 1) then
!======add on two-loop A functions which go like sum over quark charges
               bit=bit+bitrem*qsum*Q(k)**2
            endif
         elseif ((j==0).and.(k==0).and.(order >= 2 )) then
            bit=beama0(j)*beamb0(k)*gggaga
         else
            bit=zip
         endif
         
         xmsq=xmsq+bit
         
      enddo
      return
      end
