c      subroutine jeppelo1(nx,my,xx,yy,ff,cc)
c      implicit real*8(a-h,o-z)
c      dimension xx(nx),yy(my),ff(nx,my),ff1(nx,my),ff2(nx,my),
c     xff12(nx,my),yy0(4),yy1(4),yy2(4),yy12(4),z(16),wt(16,16),
c     xcl(16),cc(nx,my,4,4),iwt(16,16)

      subroutine jeppelo1(nx,my,xx,yy,ff,cc)
      implicit real*8(a-h,o-z)
      PARAMETER(NNX=49,MMY=37)
      dimension xx(nx),yy(my),ff(nx,my),ff1(NNX,MMY),ff2(NNX,MMY),
     xff12(NNX,MMY),yy0(4),yy1(4),yy2(4),yy12(4),z(16),wt(16,16),
     xcl(16),cc(nx,my,4,4),iwt(16,16)

      data iwt/1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     x		  0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,
     x		  -3,0,0,3,0,0,0,0,-2,0,0,-1,0,0,0,0,
     x		  2,0,0,-2,0,0,0,0,1,0,0,1,0,0,0,0,
     x		  0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,
     x		  0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,
     x		  0,0,0,0,-3,0,0,3,0,0,0,0,-2,0,0,-1,
     x		  0,0,0,0,2,0,0,-2,0,0,0,0,1,0,0,1,
     x		  -3,3,0,0,-2,-1,0,0,0,0,0,0,0,0,0,0,
     x		  0,0,0,0,0,0,0,0,-3,3,0,0,-2,-1,0,0,
     x		  9,-9,9,-9,6,3,-3,-6,6,-6,-3,3,4,2,1,2,
     x		  -6,6,-6,6,-4,-2,2,4,-3,3,3,-3,-2,-1,-1,-2,
     x		  2,-2,0,0,1,1,0,0,0,0,0,0,0,0,0,0,
     x		  0,0,0,0,0,0,0,0,2,-2,0,0,1,1,0,0,
     x		  -6,6,-6,6,-3,-3,3,3,-4,4,2,-2,-2,-2,-1,-1,
     x		  4,-4,4,-4,2,2,-2,-2,2,-2,-2,2,1,1,1,1/


      do 42 m=1,my
      dx=xx(2)-xx(1)
      ff1(1,m)=(ff(2,m)-ff(1,m))/dx
      dx=xx(nx)-xx(nx-1)
      ff1(nx,m)=(ff(nx,m)-ff(nx-1,m))/dx
      do 41 n=2,nx-1
      ff1(n,m)=polderiv(xx(n-1),xx(n),xx(n+1),ff(n-1,m),ff(n,m),
     xff(n+1,m))
   41 continue
   42 continue

      do 44 n=1,nx
      dy=yy(2)-yy(1)
      ff2(n,1)=(ff(n,2)-ff(n,1))/dy
      dy=yy(my)-yy(my-1)
      ff2(n,my)=(ff(n,my)-ff(n,my-1))/dy
      do 43 m=2,my-1
      ff2(n,m)=polderiv(yy(m-1),yy(m),yy(m+1),ff(n,m-1),ff(n,m),
     xff(n,m+1))
   43 continue
   44 continue

      do 46 m=1,my
      dx=xx(2)-xx(1)
      ff12(1,m)=(ff2(2,m)-ff2(1,m))/dx
      dx=xx(nx)-xx(nx-1)
      ff12(nx,m)=(ff2(nx,m)-ff2(nx-1,m))/dx
      do 45 n=2,nx-1
      ff12(n,m)=polderiv(xx(n-1),xx(n),xx(n+1),ff2(n-1,m),ff2(n,m),
     xff2(n+1,m))
   45 continue
   46 continue

      do 53 n=1,nx-1
      do 52 m=1,my-1
      d1=xx(n+1)-xx(n)
      d2=yy(m+1)-yy(m)
      d1d2=d1*d2

      yy0(1)=ff(n,m)
      yy0(2)=ff(n+1,m)
      yy0(3)=ff(n+1,m+1)
      yy0(4)=ff(n,m+1)

      yy1(1)=ff1(n,m)
      yy1(2)=ff1(n+1,m)
      yy1(3)=ff1(n+1,m+1)
      yy1(4)=ff1(n,m+1)

      yy2(1)=ff2(n,m)
      yy2(2)=ff2(n+1,m)
      yy2(3)=ff2(n+1,m+1)
      yy2(4)=ff2(n,m+1)

      yy12(1)=ff12(n,m)
      yy12(2)=ff12(n+1,m)
      yy12(3)=ff12(n+1,m+1)
      yy12(4)=ff12(n,m+1)

      do 47 k=1,4
      z(k)=yy0(k)
      z(k+4)=yy1(k)*d1
      z(k+8)=yy2(k)*d2
      z(k+12)=yy12(k)*d1d2
   47 continue

      do 49 l=1,16
      xxd=0.
      do 48 k=1,16
      xxd=xxd+iwt(k,l)*z(k)
   48 continue
      cl(l)=xxd
   49 continue
      l=0
      do 51 k=1,4
      do 50 j=1,4
      l=l+1
      cc(n,m,k,j)=cl(l)
   50 continue
   51 continue
   52 continue
   53 continue
      return
      end

      subroutine jeppelo2(x,y,nx,my,xx,yy,cc,z)
C--   G.W. 02/07/2007 Allow extrapolation to small x and large q.
      implicit real*8(a-h,o-z)
      dimension xx(nx),yy(my),cc(nx,my,4,4)      
      
      n=locx(xx,nx,x)
      m=locx(yy,my,y)
      
      if (n.gt.0.and.n.lt.nx.and.m.gt.0.and.m.lt.my) then
C--   Do usual interpolation.
         t=(x-xx(n))/(xx(n+1)-xx(n))
         u=(y-yy(m))/(yy(m+1)-yy(m))
         z=0.d0
         do l=4,1,-1
            z=t*z+((cc(n,m,l,4)*u+cc(n,m,l,3))*u
     &           +cc(n,m,l,2))*u+cc(n,m,l,1)
         enddo
         
      else if (n.eq.0.and.m.gt.0.and.m.lt.my) then
C--   Extrapolate to small x.
         call jeppelo3(xx(1),y,nx,my,xx,yy,cc,f0)
         call jeppelo3(xx(2),y,nx,my,xx,yy,cc,f1)
         if (f0.gt.0.d0.and.f1.gt.0.d0) then
            z = exp(log(f0)+(log(f1)-log(f0))/(xx(2)-xx(1))*(x-xx(1)))
         else
            z = f0+(f1-f0)/(xx(2)-xx(1))*(x-xx(1))
         end if
         
      else if (n.gt.0.and.m.eq.my) then
C--   Extrapolate to large q.
         call jeppelo3(x,yy(my),nx,my,xx,yy,cc,f0)
         call jeppelo3(x,yy(my-1),nx,my,xx,yy,cc,f1)
         if (f0.gt.0.d0.and.f1.gt.0.d0) then
            z = exp(log(f0)+(log(f0)-log(f1))/(yy(my)-yy(my-1))*
     &           (y-yy(my)))
         else
            z = f0+(f0-f1)/(yy(my)-yy(my-1))*(y-yy(my))
         end if
         
      else if (n.eq.0.and.m.eq.my) then
C--   Extrapolate to small x AND large q.
         call jeppelo3(xx(1),yy(my),nx,my,xx,yy,cc,f0)
         call jeppelo3(xx(1),yy(my-1),nx,my,xx,yy,cc,f1)
         if (f0.gt.0.d0.and.f1.gt.0.d0) then
            z0 = exp(log(f0)+(log(f0)-log(f1))/(yy(my)-yy(my-1))*
     &           (y-yy(my)))
         else
            z0 = f0+(f0-f1)/(yy(my)-yy(my-1))*(y-yy(my))
         end if
         call jeppelo3(xx(2),yy(my),nx,my,xx,yy,cc,f0)
         call jeppelo3(xx(2),yy(my-1),nx,my,xx,yy,cc,f1)
         if (f0.gt.0.d0.and.f1.gt.0.d0) then
            z1 = exp(log(f0)+(log(f0)-log(f1))/(yy(my)-yy(my-1))*
     &           (y-yy(my)))
         else
            z1 = f0+(f0-f1)/(yy(my)-yy(my-1))*(y-yy(my))
         end if
         if (z0.gt.0.d0.and.z1.gt.0.d0) then
            z = exp(log(z0)+(log(z1)-log(z0))/(xx(2)-xx(1))*(x-xx(1)))
         else
            z = z0+(z1-z0)/(xx(2)-xx(1))*(x-xx(1))
         end if

      else
C--   Set parton distribution to zero otherwise.
         z = 0.d0

      end if
      
      return
      end

C--   G.W. 02/07/2007 Copy of the original jeppelo2,
C--   only used for extrapolation.
      subroutine jeppelo3(x,y,nx,my,xx,yy,cc,z)
      implicit real*8(a-h,o-z)
      dimension xx(nx),yy(my),cc(nx,my,4,4)      
      n=locx(xx,nx,x)
      m=locx(yy,my,y)
      t=(x-xx(n))/(xx(n+1)-xx(n))
      u=(y-yy(m))/(yy(m+1)-yy(m))
      z=0.d0
      do l=4,1,-1
         z=t*z+((cc(n,m,l,4)*u+cc(n,m,l,3))*u
     &        +cc(n,m,l,2))*u+cc(n,m,l,1)
      enddo
      return
      end
