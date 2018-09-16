      subroutine phase8(r,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,wt,*)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'masses.f'
      include 'mxdim.f'
      include 'kprocess.f'
      include 'breit.f'
c******* generate phase space for 2-->4 process
c******* r(mxdim),p1(4),p2(4) are inputs reversed in sign from physical values 
c---- phase space for -p1-p2 --> p3+p4+p5+p6+p7+p8+p9+p10
c---- with all 2 pi's (ie 1/(2*pi)^20)
      integer:: nu,j
!      integer:: iflip
      real(dp):: r(mxdim)
      real(dp):: p1(4),p2(4),p3(4),p4(4),p5(4),p6(4),p7(4),p8(4)
     & ,p9(4),p10(4),p12(4),pa(4),pb(4),
     & p345(4),p678(4),p34(4),p78(4),
     & ph(4),smin,wt,wt0,wt12,wtxh,wt345,wt678,wt34,wt78,wth,dkmass
      parameter(wt0=1._dp/twopi**6)
!      data iflip/0/
!      save iflip

      wt=0._dp
      do j=1,4
      p12(j)=-p1(j)-p2(j)
      enddo
      smin=0._dp

      n2=1
      n3=0
      mass2=mt
      width2=twidth
      call phi1_2bis(r(1),r(2),r(3),r(4),p12,pb,pa,wt12,*99)

      n2=1
      n3=1

      if     (kcase==kqq_tth) then
        mass3=hmass
        width3=hwidth
        dkmass=mb
      elseif ((kcase==kqq_ttz) .or. (kcase==kqqtthz)) then
        mass3=zmass
        width3=zwidth
        dkmass=0._dp
      elseif ((kcase==kqq_ttw).or. (kcase==kttwldk))then
        mass3=wmass
        width3=wwidth
        dkmass=0._dp
      else
        write(6,*) 'Process not supported in phase8.f: ',kcase
        stop
      endif

!      if (iflip == 0) then
!        iflip=1
        call phi1_2(r(5),r(6),r(7),r(8),pa,p345,ph,wtxh,*99)
        do nu=1,4
        p678(nu)=pb(nu)
        enddo
!      elseif (iflip == 1) then 
!        iflip=0
!        call phi1_2(r(5),r(6),r(7),r(8),pa,p678,ph,wtxh,*99)
!        do nu=1,4
!        p345(nu)=pb(nu)
!        enddo
!      endif 

      mass3=wmass
      width3=wwidth
      call phi1_2m(mb,r(9),r(10),r(11),smin,p345,p5,p34,wt345,*99)
      call phi1_2m(mb,r(12),r(13),r(14),smin,p678,p6,p78,wt678,*99)
      if ((p5(4)<=0._dp).or.(p6(4)<=0._dp)) goto 99

      call phi3m0(r(15),r(16),p34,p3,p4,wt34,*99)
      if ((p3(4)<=0._dp).or.(p4(4)<=0._dp)) goto 99

      call phi3m0(r(17),r(18),p78,p7,p8,wt78,*99)
      if ((p7(4)<=0._dp).or.(p8(4)<=0._dp)) goto 99

      call phi3m(r(19),r(20),ph,p10,p9,dkmass,dkmass,wth,*99)
      if ((p9(4)<=0._dp).or.(p10(4)<=0._dp)) goto 99

      wt=wt0*wt12*wtxh*wt345*wt678*wt34*wt78*wth
      
      return
      
 99   wt=0._dp
      
      return 1
      end


