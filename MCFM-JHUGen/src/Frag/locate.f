c--- A subroutine for locating x in xx (an ordered array where 
c--- xx(n) > xx(1) has a meaning) 

      subroutine locate(xx,n,x,j) 
      implicit none
      integer n,inp,np,j,nav 
      double precision xx(n),x
       

      inp=0
      np=n+1
 10   if((np-inp) .gt. 1) then 
         nav=(np+inp)/2
         if((xx(n).gt.xx(1)).eqv.(x.gt.xx(nav))) then 
            inp=nav
         else
            np=nav
         endif
         GOTO 10
         endif
         j=inp
         return 
         end
      
