      subroutine DPOLIN2(X1A,X2A,YA,M,N,X1,X2,Y,DY) 
      implicit none
      integer nmax,mmax,j,k,M,N
      parameter(nmax=20,mmax=20)
      double precision X1A(M),X2A(N),YA(M,N),YNTMP(nmax),YMTMP(mmax) 
      double precision Y,DY,X1,X2

      do 12 j=1,M
         do 11 k=1,N
            YNTMP=YA(j,k)
 11      Continue
         call DPOLINT(X2A,YNTMP,N,X2,YMTMP(j),DY)
 12      Continue
         call DPOLINT(X1A,YMTMP,M,X1,Y,DY) 
         return 
         end


      subroutine DPOLINT(XA,YA,N,X,Y,DY) 
      implicit none
      integer nmax
      parameter(nmax=10)
      integer N,M
      double precision XA(N),YA(N),C(nmax),D(nmax) 
      double precision X,Y,DY,W,HO,HP,den
      integer ns,i
      double precision dif,difT
      

      ns=1
      dif=abs(X-XA(1))
      
      do 11 i=1,N
         difT = abs(X-XA(i))
         if (difT .lt. dif) then 
            ns=i
            dif=difT
            endif
            c(i)=YA(i) 
            d(i)=YA(i) 
 11      continue
         Y=YA(ns) 
         ns=ns-1
         do 13 m=1,N-1
            do 12 I=1,N-M
               HO=XA(i)-X
               HP=XA(i+m)-X
               W=C(i+1)-d(i)
               den=HO-HP
               if (abs(den) .lt. 1d-10) then
                  write(6,*) ' stop in polint ',den
                  stop
               endif
               den=W/den
               d(i)=HP*den
               c(i)=HO*den
 12            continue
               if (2*ns .lt. n-m) then 
                  DY=C(ns+1) 
               else
                  DY=D(ns) 
                  ns=ns-1
               endif
               Y=Y+DY
 13            continue
               return 
               end
