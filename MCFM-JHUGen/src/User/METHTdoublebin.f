      function METHTdoublebin(MET,HT)
      implicit none
      include 'types.f'
      real(dp):: METHTdoublebin
      
c--- returns a real(dp):: number that indicates which bin
c--- the double-binned MET and HT fall into
      real(dp):: MET,HT,boundMET(11),boundHT(9)
      integer:: j,bin(2),maxbinMET,maxbinHT
      data boundMET/30._dp,40._dp,50._dp,60._dp,70._dp,80._dp,90._dp,100._dp,
     &  110._dp,120._dp,130._dp/
      data maxbinMET/11/
      data boundHT/80._dp,140._dp,200._dp,260._dp,320._dp,380._dp,440._dp,500._dp,560._dp/
      data maxbinHT/9/
      save maxbinMET,maxbinHT,boundMET,boundHT
            
c--- see which bin MET falls into
      bin(1)=0
      do j=1,maxbinMET-1
        if ((MET >= boundMET(j)) .and. (MET < boundMET(j+1)))
     &    bin(1)=j
      enddo
      if (MET >= boundMET(maxbinMET)) bin(1)=maxbinMET
      
      bin(2)=0
      do j=1,maxbinHT-1
        if ((HT >= boundHT(j)) .and. (HT < boundHT(j+1)))
     &    bin(2)=j
      enddo
      if (HT >= boundHT(maxbinHT)) bin(2)=maxbinHT
      
      METHTdoublebin=real((bin(1)-1)*maxbinHT+bin(2),dp)

c--- underflow
      if ((bin(1) == 0) .or. (bin(2) == 0)) then
        METHTdoublebin=0._dp
      endif
      
      return
      end
