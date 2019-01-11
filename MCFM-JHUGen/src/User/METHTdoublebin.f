      double precision function METHTdoublebin(MET,HT)
      implicit none
c--- returns a double precision number that indicates which bin
c--- the double-binned MET and HT fall into
      double precision MET,HT,boundMET(11),boundHT(9)
      integer j,bin(2),maxbinMET,maxbinHT
      data boundMET/30d0,40d0,50d0,60d0,70d0,80d0,90d0,100d0,
     &  110d0,120d0,130d0/
      data maxbinMET/11/
      data boundHT/80d0,140d0,200d0,260d0,320d0,380d0,440d0,500d0,560d0/
      data maxbinHT/9/
      save maxbinMET,maxbinHT,boundMET,boundHT
            
c--- see which bin MET falls into
      bin(1)=0
      do j=1,maxbinMET-1
        if ((MET .ge. boundMET(j)) .and. (MET .lt. boundMET(j+1)))
     &    bin(1)=j
      enddo
      if (MET .ge. boundMET(maxbinMET)) bin(1)=maxbinMET
      
      bin(2)=0
      do j=1,maxbinHT-1
        if ((HT .ge. boundHT(j)) .and. (HT .lt. boundHT(j+1)))
     &    bin(2)=j
      enddo
      if (HT .ge. boundHT(maxbinHT)) bin(2)=maxbinHT
      
      METHTdoublebin=dfloat((bin(1)-1)*maxbinHT+bin(2))

c--- underflow
      if ((bin(1) .eq. 0) .or. (bin(2) .eq. 0)) then
        METHTdoublebin=0d0
      endif
      
      return
      end
