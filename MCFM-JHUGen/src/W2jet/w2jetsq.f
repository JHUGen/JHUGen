      subroutine w2jetsq(i1,i2,i3,i4,i5,i6,za,zb,msq)
      implicit none
      include 'types.f'

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'lc.f'
      include 'mmsq_cs.f'
      complex(dp):: qcd1(-1:1,-1:1),qcd2(-1:1,-1:1),qed(-1:1,-1:1)
      real(dp):: msq1,msq2,msqq,msq
      integer:: i1,i2,i3,i4,i5,i6

      call subqcd(i1,i2,i3,i4,i5,i6,za,zb,qcd1)
      call subqcd(i1,i2,i3,i4,i6,i5,za,zb,qcd2)

      qed(+1,+1)=qcd1(+1,+1)+qcd2(+1,+1)
      qed(+1,-1)=qcd1(+1,-1)+qcd2(-1,+1)
      qed(-1,+1)=qcd1(-1,+1)+qcd2(+1,-1)
      qed(-1,-1)=qcd1(-1,-1)+qcd2(-1,-1)

      msq1= abs(qcd1(+1,+1))**2+abs(qcd1(+1,-1))**2
     &     +abs(qcd1(-1,+1))**2+abs(qcd1(-1,-1))**2

      msq2= abs(qcd2(+1,+1))**2+abs(qcd2(+1,-1))**2
     &     +abs(qcd2(-1,+1))**2+abs(qcd2(-1,-1))**2

      msqq= abs( qed(+1,+1))**2+abs( qed(+1,-1))**2
     &     +abs( qed(-1,+1))**2+abs( qed(-1,-1))**2

      mmsq_cs(0,+1,+1)=0._dp
      mmsq_cs(1,+1,+1)=0._dp
      mmsq_cs(2,+1,+1)=0._dp

      if ((colourchoice == 1) .or. (colourchoice == 0)) then
      mmsq_cs(1,+1,+1)=msq1
      mmsq_cs(2,+1,+1)=msq2
      endif
      if ((colourchoice == 2) .or. (colourchoice == 0)) then
      mmsq_cs(0,+1,+1)=-ninth*msqq
      endif

      msq=mmsq_cs(0,+1,+1)+mmsq_cs(1,+1,+1)+mmsq_cs(2,+1,+1)

      return
      end
