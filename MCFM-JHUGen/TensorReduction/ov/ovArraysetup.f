      subroutine ovArraysetup
      include 'TRydef.f'
      integer n1,n2,n3,n4,n5,n6,n7,i
      integer nn2(2),nn3(3),nn4(4),nn5(5),nn6(6),nn7(7)
      logical,save:: first=.true.
!omp threadprivate(first)
      if (first) then
      first=.false.
      else
      return
      endif 
      write(6,*) 'setting up the arrays for the tensors'
C--fill ordered terms in array 
      do n1=1,4
      y1(n1)=n1
      enddo

      i=0
      do n1=1,4
      do n2=n1,4
      i=i+1
      y2(n1,n2)=i
      enddo
      enddo
      do n1=1,4
      do n2=1,4
      nn2(1)=n1
      nn2(2)=n2
      call ovIpiksrt(2,nn2)
      y2(n1,n2)=y2(nn2(1),nn2(2))
      enddo
      enddo

      i=0
      do n1=1,4
      do n2=n1,4
      do n3=n2,4
      i=i+1
      y3(n1,n2,n3)=i
      enddo
      enddo
      enddo
      do n1=1,4
      do n2=1,4
      do n3=1,4
      nn3(1)=n1
      nn3(2)=n2
      nn3(3)=n3
      call ovIpiksrt(3,nn3)
      y3(n1,n2,n3)=y3(nn3(1),nn3(2),nn3(3))
      enddo
      enddo
      enddo

      i=0
      do n1=1,4
      do n2=n1,4
      do n3=n2,4
      do n4=n3,4
      i=i+1
      y4(n1,n2,n3,n4)=i
      enddo
      enddo
      enddo
      enddo
      do n1=1,4
      do n2=1,4
      do n3=1,4
      do n4=1,4
      nn4(1)=n1
      nn4(2)=n2
      nn4(3)=n3
      nn4(4)=n4
      call ovIpiksrt(4,nn4)
      y4(n1,n2,n3,n4)=y4(nn4(1),nn4(2),nn4(3),nn4(4))
      enddo
      enddo
      enddo
      enddo


      i=0
      do n1=1,4
      do n2=n1,4
      do n3=n2,4
      do n4=n3,4
      do n5=n4,4
      i=i+1
      y5(n1,n2,n3,n4,n5)=i
      enddo
      enddo
      enddo
      enddo
      enddo
      do n1=1,4
      do n2=1,4
      do n3=1,4
      do n4=1,4
      do n5=1,4
      nn5(1)=n1
      nn5(2)=n2
      nn5(3)=n3
      nn5(4)=n4
      nn5(5)=n5
      call ovIpiksrt(5,nn5)
      y5(n1,n2,n3,n4,n5)=y5(nn5(1),nn5(2),nn5(3),nn5(4),nn5(5))
      enddo
      enddo
      enddo
      enddo
      enddo

      i=0
      do n1=1,4
      do n2=n1,4
      do n3=n2,4
      do n4=n3,4
      do n5=n4,4
      do n6=n5,4
      i=i+1
      y6(n1,n2,n3,n4,n5,n6)=i
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      do n1=1,4
      do n2=1,4
      do n3=1,4
      do n4=1,4
      do n5=1,4
      do n6=1,4
      nn6(1)=n1
      nn6(2)=n2
      nn6(3)=n3
      nn6(4)=n4
      nn6(5)=n5
      nn6(6)=n6
      call ovIpiksrt(6,nn6)
      y6(n1,n2,n3,n4,n5,n6)=
     . y6(nn6(1),nn6(2),nn6(3),nn6(4),nn6(5),nn6(6))
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo

      i=0
      do n1=1,4
      do n2=n1,4
      do n3=n2,4
      do n4=n3,4
      do n5=n4,4
      do n6=n5,4
      do n7=n6,4
      i=i+1
      y7(n1,n2,n3,n4,n5,n6,n7)=i
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      do n1=1,4
      do n2=1,4
      do n3=1,4
      do n4=1,4
      do n5=1,4
      do n6=1,4
      do n7=1,4
      nn7(1)=n1
      nn7(2)=n2
      nn7(3)=n3
      nn7(4)=n4
      nn7(5)=n5
      nn7(6)=n6
      nn7(7)=n7
      call ovIpiksrt(7,nn7)
      y7(n1,n2,n3,n4,n5,n6,n7)=
     . y7(nn7(1),nn7(2),nn7(3),nn7(4),nn7(5),nn7(6),nn7(7))
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo

      return
      end
