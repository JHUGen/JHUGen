      subroutine xzqqgg(mqqb)
      implicit none
      include 'types.f'

C     Author J.M.Campbell, February 2000
C     Returns the amplitudes squared for the process
C     0 ---> q(p1)+g(p2)+g(p3)+qbar(p4)+l(p5)+a(p6)
C     mqqb(2,2) has two indices;the first for the helicity quark line;
C     the second for helicity of lepton line.
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'zprods_com.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'mmsq_cs.f'
      real(dp):: mqqb(2,2),m0,fac
      complex(dp):: tempm0,m(2)
      complex(dp):: a6treeg1
      integer:: j,lh,h2,h3,hq,h(2:3)
      integer,parameter::i1(2)=(/1,4/),i2(2)=(/2,3/),i3(2)=(/3,2/),
     &                   i4(2)=(/4,1/),i5(2)=(/6,5/),i6(2)=(/5,6/)
      character*9,parameter:: st(2,2)=
     & reshape((/'q+g-g-qb-','q+g-g+qb-','q+g+g-qb-','q+g+g+qb-'/)
     & ,(/2,2/))
      include 'cplx.h'

c      compare=.false.

C ---final matrix element squared is needed as function of quark line helicity
C----and lepton line helicity
C----first argument is quark line helicity
C----second argument is lepton line helicity
      mqqb(1,1)=0._dp
      mqqb(1,2)=0._dp
      mqqb(2,1)=0._dp
      mqqb(2,2)=0._dp

      fac=avegg*2._dp*gsq**2*esq**2*cf*xn**2
c--- extra factor of 4 due to colour matrix normalization (rt2**4)
      fac=fac*4._dp

c--- USED TO COMPARE ONLY
c      if (compare) then
c      do hq=1,2
c      do lh=1,2
c      mqqb(hq,lh)=0._dp

c      do h2=1,2
c      do h3=1,2
c        h(2)=h2
c        h(3)=h3
C initialize loop sum to zero
c        tempm0=czip
c        do j=1,2
c        m(j)=amp_qqgg(1,hq,i2(j),h(i2(j)),i3(j),h(i3(j)),4,lh)
c        write(*,*) hq,h2,h3,lh,i2(j),i3(j),m(j)
c        tempm0=tempm0+m(j)
c        enddo
c      m1=Dble(conjg(m(1))*m(2))
c      m0=abs(tempm0)**2
c      mqqb(hq,lh)=mqqb(hq,lh)+fac*(m0-x*m1)
c      enddo
c      enddo
c      write(*,*) 'old mqqb(',hq,',',lh,')',mqqb(hq,lh)
c      enddo
c      enddo
c      endif
c--- USED TO COMPARE ONLY

      do hq=1,2
      do lh=1,2
      mqqb(hq,lh)=0._dp
      mmsq_cs(0,hq,lh)=0._dp
      mmsq_cs(1,hq,lh)=0._dp
      mmsq_cs(2,hq,lh)=0._dp

      do h2=1,2
      do h3=1,2
        h(2)=h2
        h(3)=h3
C initialize loop sum to zero
        tempm0=czip
        do j=1,2
        if (hq == 1) then
        m(j)=a6treeg1(st(3-h(i2(j)),3-h(i3(j))),
     &     i1(1),i2(j),i3(j),i4(1),i6(lh),i5(lh),zb,za)
        else
        m(j)=a6treeg1(st(h(i2(j)),h(i3(j))),
     &     i1(1),i2(j),i3(j),i4(1),i5(lh),i6(lh),za,zb)
c  if hq=2,lh=2,h2=2,h3=2,j=1 then st='q+g+g+qb-'
c  and we are calculating a6treeg1('q+g+g+qb-',1,2,3,4,5,6,za,zb)
C  which corresponds to
c  q(4)+l(5) --> q_R(1)+l_R(6)+g_R(2)+g_R(3)
c  hence, this routine thinks that particle 1 is an outgoing quark
        endif
        tempm0=tempm0+m(j)
        enddo
      m0=abs(tempm0)**2
      mmsq_cs(0,hq,lh)=mmsq_cs(0,hq,lh)-fac*m0/xn**2
      mmsq_cs(1,hq,lh)=mmsq_cs(1,hq,lh)+fac*abs(m(1))**2
      mmsq_cs(2,hq,lh)=mmsq_cs(2,hq,lh)+fac*abs(m(2))**2
      enddo
      enddo
      mqqb(hq,lh)=mmsq_cs(1,hq,lh)+mmsq_cs(2,hq,lh)+mmsq_cs(0,hq,lh)
c      write(*,*) 'tree ',hq,lh,mqqb(hq,lh)
c      if (compare) write(*,*) 'new mqqb(',hq,',',lh,')',mqqb(hq,lh)
      enddo
      enddo
c--- USED TO COMPARE ONLY
c      if (compare) then
c      tamp=a6treeg('q+g+g+qb-',4,3,2,1,6,5,zb,za)
c--- This amplitude with the choice of arguments we have made
c--- corresponds to
c    q(1)+l(6) --> q_R(4)+l_R(5)+g_R(2)+g_R(3)
c      write(*,*) 'q+g-g-qb-',tamp
c      tamp=a6treeg('q+g-g+qb-',1,2,3,4,5,6,za,zb)
c      write(*,*) 'q+g-g+qb-',tamp
c      tamp=a6treeg('q+g+g-qb-',1,2,3,4,5,6,za,zb)
c      write(*,*) 'q+g+g-qb-',tamp
c      tamp=a6treeg('q+g+g+qb-',1,2,3,4,5,6,za,zb)
c      write(*,*) 'q+g+g+qb-',tamp
c      tamp=a6treeg('q+g+qb-g-',1,2,3,4,5,6,za,zb)
c      write(*,*) 'q+g+qb-g-',tamp
c      tamp=a6treeg('q+g+qb-g+',1,2,3,4,5,6,za,zb)
c      write(*,*) 'q+g+qb-g+',tamp
c      tamp=a6treeg('q+qb-g-g+',1,2,3,4,5,6,za,zb)
c      write(*,*) 'q+qb-g-g+',tamp
c      tamp=a6treeg('q+qb-g+g-',1,2,3,4,5,6,za,zb)
c      write(*,*) 'q+qb-g+g-',tamp
c      tamp=a6treeg('q+qb-g+g+',1,2,3,4,5,6,za,zb)
c      write(*,*) 'q+qb-g+g+',tamp
c      pause
c      endif
c--- USED TO COMPARE ONLY

      return
      end

      function a6treeg1(st,j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: a6treeg1

c----wrapper to a6treeg that also includes config st='q+g-g-qb-'
      integer:: j1,j2,j3,j4,j5,j6
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      character*9 st
      complex(dp):: a6treeg

      if(st=='q+g-g-qb-') then
        a6treeg1=a6treeg('q+g+g+qb-',j4,j3,j2,j1,j6,j5,zb,za)
      else
        a6treeg1=a6treeg(st,j1,j2,j3,j4,j5,j6,za,zb)
      endif

      return
      end

