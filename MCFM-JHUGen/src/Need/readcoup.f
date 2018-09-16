      subroutine readcoup
      implicit none
      include 'types.f'
c--- reads in the anomalous couplings from the file anomcoup.DAT
      
      include 'anomcoup.f'
      include 'zerowidth.f'
     
c      if (newinput) goto 20     
      
c      open(unit=21,file='anomcoup.DAT',status='old',err=999)
c      call checkversion(21,'anomcoup.DAT')
c      read(21,*) delg1_z
c      read(21,*) delk_z
c      read(21,*) delk_g
c      read(21,*) lambda_z
c      read(21,*) lambda_g
c      read(21,*) tevscale
c      close(21)
c--- E-M gauge invariance requires that delg1_g=0
c      delg1_g=0._dp
      
c   20 continue   
      
      write(6,*)
      write(6,*)  '*************** Anomalous couplings ****************'
      write(6,*)  '*                                                  *'
      write(6,99) '*            Delta_g1(Z)  =  ',delg1_z,
     &                '                *'
      write(6,99) '*            Delta_g1(g)  =  ',0._dp,
     &                '                *'
      write(6,99) '*            Delta_K(Z)   =  ',delk_z,
     &                '                *'
      write(6,99) '*            Delta_K(g)   =  ',delk_g,
     &                '                *'
      write(6,99) '*            Lambda(Z)    =  ',lambda_z,
     &                '                *'
      write(6,99) '*            Lambda(g)    =  ',lambda_g,
     &                '                *'
      write(6,99) '*            TeV-scale    =  ',tevscale,
     &                ' TeV            *'
      write(6,*)  '****************************************************'
      
c--- Check to see whether anomalous couplings are being used
      if (max(abs(delg1_z),abs(lambda_g),abs(lambda_z),
     &        abs(delk_g),abs(delk_z)) > 1.e-8_dp) then
        anomtgc=.true.
      if (zerowidth .eqv. .false.) then
        write(6,*)
        write(6,*)'********************** WARNING *********************'
        write(6,*)'*                                                  *'
        write(6,*)'* No singly resonant diagrams will be included     *'
        write(6,*)'*  when anomalous TGCs are present.                *'
        write(6,*)'*                                                  *'
        write(6,*)'****************************************************'
        endif
      else
        anomtgc=.false.
      endif
      return
 
   99 format(1x,a29,f6.2,a17)
   
c  999 continue
      write(6,*) 'Error reading anomcoup.DAT'
      call flush(6)
      stop
 
      end
