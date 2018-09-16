      function etdoublebin(pt1,pt2)
      implicit none
      include 'types.f'
      real(dp):: etdoublebin
      
c--- returns a real(dp):: number that indicates which bin
c--- the double-binned et's pt1,pt2 fall into
      real(dp):: pt1,pt2,testpt(2),bound(20)
      integer:: i,j,bin(2),maxbin
      data bound/20._dp,38.1_dp,72.5_dp,138._dp,263._dp,500._dp,14*0._dp/
      data maxbin/5/
      
c--- first order the pt's so that testpt(1) is largest      
      if (pt1 > pt2) then
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
        if ((testpt(i) >= bound(j)) .and. (testpt(i) < bound(j+1)))
     &    bin(i)=j
        enddo
        if (bin(i) == 0) bin(i)=maxbin
      enddo
      
      etdoublebin=real((bin(1)-1)*maxbin+bin(2),dp)
        
      return
      end
