      subroutine smalltau(p,npart,*)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'cutoff.f'
      include 'mxpart.f'
      include 'first.f'
      include 'plabel.f'
      integer npart,i5,i6,i7
      integer, save :: ipp
!$omp threadprivate(ipp)
      real(dp):: p(mxpart,4),t15,t25,t16,t26,t56,t17,t27,t57,t67,
     & n1(4),n2(4),n5(4),n6(4),n7(4)

c---  determine beginning of parton entries in plabel
      if (first) then
        first=.false.
        ipp=3
        do while ((ipp < mxpart ) .and. (plabel(ipp) .ne. 'pp'))
          ipp=ipp+1
        enddo
        if (ipp == mxpart) then
          write(6,*) 'Could not identify partons in smalltau.f'
          stop
        endif
c        write(6,*) 'found ipp=',ipp  
      endif
      i5=ipp
      i6=ipp+1
      i7=ipp+2
      
c--- do not apply any cuts if no partons in final state
      if (npart < i5-2) return
      
c--- reference directions for beam axes
      n1(:)=(/0._dp,0._dp, 1._dp,1._dp/)
      n2(:)=(/0._dp,0._dp,-1._dp,1._dp/)
c--- reference direction for particle 5
      n5(:)=p(i5,:)/p(i5,4)

c--- compute N-jettiness tau contributions
      t15=abs(n1(4)*p(i5,4)-n1(1)*p(i5,1)-n1(2)*p(i5,2)-n1(3)*p(i5,3))
      t25=abs(n2(4)*p(i5,4)-n2(1)*p(i5,1)-n2(2)*p(i5,2)-n2(3)*p(i5,3))

c--- apply small cut based on value of cutoff (set in input.DAT)
      if (t15 < cutoff) return 1
      if (t25 < cutoff) return 1

c--- no more cuts if less than 2 partons in final state
      if (npart < i6-2) return

c--- reference direction for particle 6 
      n6(:)=p(i6,:)/p(i6,4)

      t16=abs(n1(4)*p(i6,4)-n1(1)*p(i6,1)-n1(2)*p(i6,2)-n1(3)*p(i6,3))
      t26=abs(n2(4)*p(i6,4)-n2(1)*p(i6,1)-n2(2)*p(i6,2)-n2(3)*p(i6,3))

      if (p(i5,4) > p(i6,4)) then
        t56=n5(4)*p(i6,4)-n5(1)*p(i6,1)-n5(2)*p(i6,2)-n5(3)*p(i6,3)
      else
        t56=n6(4)*p(i5,4)-n6(1)*p(i5,1)-n6(2)*p(i5,2)-n6(3)*p(i5,3)
      endif

      if (t16 < cutoff) return 1
      if (t26 < cutoff) return 1
      if (t56 < cutoff) return 1
      
c--- no more cuts if less than 3 partons in final state
      if (npart < i7-2) return

c--- reference direction for particle 7
      n7(:)=p(i7,:)/p(i7,4)

      t17=abs(n1(4)*p(i7,4)-n1(1)*p(i7,1)-n1(2)*p(i7,2)-n1(3)*p(i7,3))
      t27=abs(n2(4)*p(i7,4)-n2(1)*p(i7,1)-n2(2)*p(i7,2)-n2(3)*p(i7,3))

      if (p(i5,4) > p(i7,4)) then
        t57=n5(4)*p(i7,4)-n5(1)*p(i7,1)-n5(2)*p(i7,2)-n5(3)*p(i7,3)
      else
        t57=n7(4)*p(i5,4)-n7(1)*p(i5,1)-n7(2)*p(i5,2)-n7(3)*p(i5,3)
      endif

      if (p(i6,4) > p(i7,4)) then
        t67=n6(4)*p(i7,4)-n6(1)*p(i7,1)-n6(2)*p(i7,2)-n6(3)*p(i7,3)
      else
        t67=n7(4)*p(i6,4)-n7(1)*p(i6,1)-n7(2)*p(i6,2)-n7(3)*p(i6,3)
      endif

      if (t17 < cutoff) return 1
      if (t27 < cutoff) return 1
      if (t57 < cutoff) return 1
      if (t67 < cutoff) return 1
      
      return
      end
      
