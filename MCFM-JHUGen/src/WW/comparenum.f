      subroutine comparenum(nleg,cmax,coeff)
      implicit none
      include 'types.f'
      
      integer:: nleg,cmax,corder(cmax),h1,h2,i
      complex(dp):: coeff(2,2,cmax)
c--- writes out integral coefficients in the same order as in
c--- the output of the numerical code
      integer,parameter:: corder2(7)=(/6,1,2,4,3,5,7/),
     & corder3(12)=(/12,1,7,11,8,2,6,3,9,10,4,5/),
     & corder4(6)=(/5,1,6,2,3,4/)
      
      if (nleg == 2) then
        do i=1,cmax
      corder(i)=corder2(i)
      enddo
      endif
      
      if (nleg == 3) then
        do i=1,cmax
      corder(i)=corder3(i)
      enddo
      endif
      
      if (nleg == 4) then
        do i=1,cmax
      corder(i)=corder4(i)
      enddo
      endif
      
      do h1=1,2
      do h2=1,2
      write(6,*) 'h1=',h1,', h2=',h2
      
      do i=1,cmax
        if (abs(coeff(h1,h2,corder(i))) > 1.e-10_dp) then
        write(6,*) corder(i),
     &               coeff(h1,h2,corder(i)),
     &           abs(coeff(h1,h2,corder(i)))
        endif
      enddo
      write(6,*)
      
      enddo
      enddo
      
c      pause
      
      return
      end
      
