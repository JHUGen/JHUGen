      
      include 'ckmallowed.f'
      include 'montecarlorpp.f'
      integer:: i,k,nwz,awx,lwx
      integer:: i1,i2,i3,i4,i5,i6,ii1,ii2,ii3,ii4,ii5,ii6
      character*2 ch(-5:5),lb(-5:5),lb1,lb2,lb5,lb6
      ch(-5)='b~'
      ch(-4)='c~'
      ch(-3)='s~'
      ch(-2)='u~'
      ch(-1)='d~'
      ch(0)='g'
      ch(1)='d'
      ch(2)='u'
      ch(3)='s'
      ch(4)='c'
      ch(5)='b'

      lb(-5)='bb'
      lb(-4)='cb'
      lb(-3)='sb'
      lb(-2)='ub'
      lb(-1)='db'
      lb(0)='g0'
      lb(1)='dq'
      lb(2)='uq'
      lb(3)='sq'
      lb(4)='cq'
      lb(5)='bq'
c--- fill allowed CKM combinations
      do i=-5,5
      do k=-5,5
      ckmallowed(i,k)=.false.
      enddo
      enddo
      nwz=1      
c--- diagonal only
      if (nwz == 1) then
         ckmallowed(2,-1)=.true.
         ckmallowed(2,-3)=.true.
         ckmallowed(4,-1)=.true.
         ckmallowed(4,-3)=.true.
         ckmallowed(-1,2)=.true.
         ckmallowed(-3,2)=.true.
         ckmallowed(-1,4)=.true.
         ckmallowed(-3,4)=.true.
         awx=ma_pdg
         lwx=nml_pdg         
      elseif (nwz == -1) then
         ckmallowed(-2,1)=.true.
         ckmallowed(-2,3)=.true.
         ckmallowed(-4,1)=.true.
         ckmallowed(-4,3)=.true.
         ckmallowed(1,-2)=.true.
         ckmallowed(3,-2)=.true.
         ckmallowed(1,-4)=.true.
         ckmallowed(3,-4)=.true.
         awx=nea_pdg
         lwx=el_pdg         
      else
      write(6,*) 'Error in init_processes.f'
      write(6,*) 'Set jwnz to -1 or +1 in input, current value=',nwz
      stop
      endif

      nwz=1
         do i=-4,4
         do k=-4,4
           if(i+k .ne. nwz) ckmallowed(i,k)=.false.
         enddo
         enddo         

      ii1=1
      ii2=2
      ii3=3
      ii4=4
      ii5=5
      ii6=6

      open(unit=66,file='bit',status='unknown')
c      flst_nborn=0
C----annhilation topology
      do i1=-5,5
      do i2=-5,5
      do i5=-5,5
      do i6=-5,5

      if ((ckmallowed(i1,i2)) 
     & .and. (i5 == -i6) 
     & .and. (i5 >= 0) ) then
      lb1=lb(i1)
      lb2=lb(i2)
      lb5=lb(min(i5,i6))
      lb6=lb(max(i5,i6))
      write(66,*) ch(i1),' ',ch(i2),' -> ve e+ ',
     & ch(min(i5,i6)),' ',ch(max(i5,i6))
      write(66,*) '2'
      write(66,*) 'yes'
      write(66,*) lb1//lb2,'_veep',lb5//lb6


      elseif (ckmallowed(i1,-i5)
     &  .and. (i2 == i6)
     &  ) then
      lb1=lb(i1)
      lb2=lb(i2)
      lb5=lb(min(i5,i6))
      lb6=lb(max(i5,i6))
      write(66,*) ch(i1),' ',ch(i2),' -> ve e+ ',
     & ch(min(i5,i6)),' ',ch(max(i5,i6))
      write(66,*) '2'
      write(66,*) 'yes'
      write(66,*) lb1//lb2,'_veep',lb5//lb6

      elseif (ckmallowed(i2,-i5)
     &  .and. (i1 == i6)
     &  ) then
      lb1=lb(i1)
      lb2=lb(i2)
      lb5=lb(min(i5,i6))
      lb6=lb(max(i5,i6))
      write(66,*) ch(i1),' ',ch(i2),' -> ve e+ ',
     & ch(min(i5,i6)),' ',ch(max(i5,i6))
      write(66,*) '2'
      write(66,*) 'yes'
      write(66,*) lb1//lb2,'_veep',lb5//lb6


      elseif (ckmallowed(-i5,-i6)
     &  .and. (i1 == -i2)
     &  .and. (i5 < 0)
     &  ) then
      lb1=lb(i1)
      lb2=lb(i2)
      lb5=lb(min(i5,i6))
      lb6=lb(max(i5,i6))
      write(66,*) ch(i1),' ',ch(i2),' -> ve e+ ',
     & ch(min(i5,i6)),' ',ch(max(i5,i6))
      write(66,*) '2'
      write(66,*) 'yes'
      write(66,*) lb1//lb2,'_veep',lb5//lb6

      endif
      enddo
      enddo
      enddo
      enddo
      stop 
      end
