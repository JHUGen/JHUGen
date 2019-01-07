      subroutine Array2dim
      include 'Carraydef.f'
      integer n1,n2,n3,n4,n5,n6,n7,i
      integer nn2(2),nn3(3),nn4(4),nn5(5),nn6(6),nn7(7)
      logical,save:: first=.true.
!$omp threadprivate(first)
      
      if (first) then
        first=.false.
      else
        return
      endif 

      write(6,*) 'setting up the 2^n arrays'
      
C--fill ordered terms in array 
      do n1=1,2
      do n2=1,2
      if (n1.eq.n2) then
      delta(n1,n2)=1d0
      else
      delta(n1,n2)=0d0
      endif
      enddo
      enddo
      do n1=1,2
      z1(n1)=n1
      enddo

      i=0
      do n1=1,2
      do n2=n1,2
      i=i+1
      z2(n1,n2)=i
      enddo
      enddo
      do n1=1,2
      do n2=1,2
      nn2(1)=n1
      nn2(2)=n2
      call pvIpiksrt(2,nn2)
      z2(n1,n2)=z2(nn2(1),nn2(2))
      enddo
      enddo

      i=0
      do n1=1,2
      do n2=n1,2
      do n3=n2,2
      i=i+1
      z3(n1,n2,n3)=i
      enddo
      enddo
      enddo
      do n1=1,2
      do n2=1,2
      do n3=1,2
      nn3(1)=n1
      nn3(2)=n2
      nn3(3)=n3
      call pvIpiksrt(3,nn3)
      z3(n1,n2,n3)=z3(nn3(1),nn3(2),nn3(3))
      enddo
      enddo
      enddo


      i=0
      do n1=1,2
      do n2=n1,2
      do n3=n2,2
      do n4=n3,2
      i=i+1
      z4(n1,n2,n3,n4)=i
      enddo
      enddo
      enddo
      enddo
      do n1=1,2
      do n2=1,2
      do n3=1,2
      do n4=1,2
      nn4(1)=n1
      nn4(2)=n2
      nn4(3)=n3
      nn4(4)=n4
      call pvIpiksrt(4,nn4)
      z4(n1,n2,n3,n4)=z4(nn4(1),nn4(2),nn4(3),nn4(4))
      enddo
      enddo
      enddo
      enddo


      i=0
      do n1=1,2
      do n2=n1,2
      do n3=n2,2
      do n4=n3,2
      do n5=n4,2
      i=i+1
      z5(n1,n2,n3,n4,n5)=i
      enddo
      enddo
      enddo
      enddo
      enddo
      do n1=1,2
      do n2=1,2
      do n3=1,2
      do n4=1,2
      do n5=1,2
      nn5(1)=n1
      nn5(2)=n2
      nn5(3)=n3
      nn5(4)=n4
      nn5(5)=n5
      call pvIpiksrt(5,nn5)
      z5(n1,n2,n3,n4,n5)=z5(nn5(1),nn5(2),nn5(3),nn5(4),nn5(5))
      enddo
      enddo
      enddo
      enddo
      enddo


      i=0
      do n1=1,2
      do n2=n1,2
      do n3=n2,2
      do n4=n3,2
      do n5=n4,2
      do n6=n5,2
      i=i+1
      z6(n1,n2,n3,n4,n5,n6)=i
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      do n1=1,2
      do n2=1,2
      do n3=1,2
      do n4=1,2
      do n5=1,2
      do n6=1,2
      nn6(1)=n1
      nn6(2)=n2
      nn6(3)=n3
      nn6(4)=n4
      nn6(5)=n5
      nn6(6)=n6
      call pvIpiksrt(6,nn6)
      z6(n1,n2,n3,n4,n5,n6)
     . =z6(nn6(1),nn6(2),nn6(3),nn6(4),nn6(5),nn6(6))
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo

      i=0
      do n1=1,2
      do n2=n1,2
      do n3=n2,2
      do n4=n3,2
      do n5=n4,2
      do n6=n5,2
      do n7=n6,2
      i=i+1
      z7(n1,n2,n3,n4,n5,n6,n7)=i
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo

      do n1=1,2
      do n2=1,2
      do n3=1,2
      do n4=1,2
      do n5=1,2
      do n6=1,2
      do n7=1,2
      nn7(1)=n1
      nn7(2)=n2
      nn7(3)=n3
      nn7(4)=n4
      nn7(5)=n5
      nn7(6)=n6
      nn7(7)=n7
      call pvIpiksrt(7,nn7)
      z7(n1,n2,n3,n4,n5,n6,n7)=
     . z7(nn7(1),nn7(2),nn7(3),nn7(4),nn7(5),nn7(6),nn7(7))
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo

      return
      end
