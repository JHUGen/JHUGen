      double precision function etdoublebin(pt1,pt2)
      implicit none
c--- returns a double precision number that indicates which bin
c--- the double-binned et's pt1,pt2 fall into
      double precision pt1,pt2,testpt(2),bound(20)
      integer i,j,bin(2),maxbin
      data bound/20d0,38.1d0,72.5d0,138d0,263d0,500d0,14*0d0/
      data maxbin/5/
      
c--- first order the pt's so that testpt(1) is largest      
      if (pt1 .gt. pt2) then
        testpt(1)=pt1
        testpt(2)=pt2
      else
        testpt(1)=pt2
        testpt(2)=pt1
      endif
      
c--- see which bin each pt falls into
      do i=1,2
        bin(i)=0
        do j=1,maxbin-1
        if ((testpt(i) .ge. bound(j)) .and. (testpt(i) .lt. bound(j+1)))
     .    bin(i)=j
        enddo
        if (bin(i) .eq. 0) bin(i)=maxbin
      enddo
      
      etdoublebin=dfloat((bin(1)-1)*maxbin+bin(2))
        
      return
      end
