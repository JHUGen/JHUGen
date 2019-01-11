      subroutine twoloop(s,t,u,FL,FSL)
C-----Formula taken from hep-ph/0109078v1
      implicit none
      include 'constants.f'
      double complex FL(2,2,2,2),FSL(2,2,2,2)
      double complex FmpppL,FppmpL,FmpmpL,FmpppSL,FmmppSL,FmpmpSL
      double precision s,t,u,x,y,lX,lY,ddilog,Li3,Li4,zeta3,zeta4
      parameter(zeta3=1.202056903159594d0,zeta4=1.082323233711138d0)

c--- statement functions
      FmpppL(x,y,lX,lY)=
     & ((2d0+4d0*x/y**2-5d0*x**2/y**2)*((lX+im*pi)**2+pisq)
     & -(1d0-x*y)*((lX-lY)**2+pisq)
     & +2d0*(9d0/y-10d0*x)*(lX+im*pi))/8d0

      FppmpL(x,y,lX,lY)=
     & ((2d0+6d0*x/y**2-3d0*x**2/y**2)*((lX+im*pi)**2+pisq)
     & -(x-y)**2*((lX-lY)**2+pisq)
     & +2d0*(9d0/y-8d0*x)*(lX+im*pi))/8d0

      FmpmpL(x,y,lX,lY)=
     & -2d0*(x**2+1d0/y**2)*(Li4(-x)-zeta4
     & -0.5d0*(lX+im*pi)*(Li3(-x)-zeta3)

     & +pisqo6*(ddilog(-x)-pisqo6-0.5d0*lX**2)-(1d0/48d0)*lX**4

     & +1d0/24d0*(lX+im*pi)**2*((lX+im*pi)**2+pisq))
     & +2d0*(3d0*(1d0-x)**2-2d0)/y**2
     & *(Li4(-x)+Li4(-x/y)-Li4(-y)-(lY+im*pi)*(Li3(-x)-zeta3)
     & +pisqo6*(ddilog(-x)+0.5d0*lY**2)
     & -lX*lY**3/6d0+lY**4/24d0-7d0/360d0*pi**4)
     & -2d0/3d0*(8d0-x+30d0*x/y)*(Li3(-y)-zeta3
     & -(lY+im*pi)*(ddilog(-y)-pisqo6)-0.5d0*lX*((lY+im*pi)**2+pisq))
     & +(4d0*y+27d0+42d0/y+4d0/y**2)/6d0
     & *(Li3(-x)-zeta3-(lX+im*pi)*(ddilog(-x)-pisqo6)
     & +im*pi/2d0*lX**2-pisq*lX)
     & +(3d0-2d0/y-12d0*x/y**2)/12d0*(lX+im*pi)*((lX+im*pi)**2+pisq)
     & -y*(lX+im*pi)*((lY+im*pi)**2+pisq)/3d0
     & +2d0*(1d0+2d0/y)*(zeta3-pisqo6*(lY+im*pi))
     & +(y**2-24d0*y+44d0-8d0*x**3/y)/24d0*((lX-lY)**2+pisq)
     & -(15d0-14d0*x/y-48d0*x/y**2)/24d0*((lX+im*pi)**2+pisq)
     & +(8d0*x/y+60d0-24d0*y/x+27d0*y**2/x**2)/24d0
     & *((lY+im*pi)**2+pisq)
     & +(4d0/9d0)*pisq*(x/y)
     & +(2d0*x**2-54d0*x-27d0*y**2)/12d0*((lX+im*pi)/y+(lY+im*pi)/x)

      FmpppSL(x,y,lX,lY)=
     & ((x**2+1d0/y**2)*((lX+im*pi)**2+pisq)
     & +0.5d0*(x**2+y**2)*((lX-lY)**2+pisq)
     & -4d0*(1d0/y-x)*(lX+im*pi))/8d0

      FmmppSL(x,y,lX,lY)=
     & -2d0*x**2*(Li4(-x)+Li4(-y)-(lX+im*pi)*(Li3(-x)+Li3(-y))

     & +lX**4/12d0-lX**3*lY/3d0+pisq/12d0*lX*lY-4d0/90d0*pi**4
     & +im*pi/6d0*lX*(lX**2-3d0*lX*lY+pisq))
     & -(x-y)*(Li4(-x/y)-pisqo6*ddilog(-x))

     & -x*(2d0*Li3(-x)-Li3(-x/y)-3d0*zeta3-2d0*(lX+im*pi)*ddilog(-x)

     & +(lX-lY)*(ddilog(-x/y)+lX**2)
     & +(5d0*(lX-lY)+18d0*im*pi)*((lX-lY)**2+pisq)/12d0

     & -(2d0/3d0)*lX*(lX**2+pisq)-im*pi*(lY**2+pisq))

     & +(1d0-2d0*x**2)/(4d0*y**2)*((lX+im*pi)**2+pisq)
     & -0.125d0*(2d0*x*y+3d0)*((lX-lY)**2+pisq)+pisq/12d0

     & +(0.5d0/y+x)*(lX+im*pi)-0.25d0

      FmpmpSL(x,y,lX,lY)=-2d0*(x**2+1d0)/y**2
     & *(Li4(-x/y)-Li4(-y)+0.5d0*(lX-2d0*lY-im*pi)*(Li3(-x)-zeta3)
     & +1d0/24d0*(lX**4+2d0*im*pi*lX**3-4d0*lX*lY**3+lY**4
     & +2d0*pisq*lY**2)+7d0/360d0*pi**4)

     & -2d0*(x-1d0)/y*(Li4(-x)-zeta4-0.5d0*(lX+im*pi)*(Li3(-x)-zeta3)
     & +pisqo6*(ddilog(-x)-pisqo6-0.5d0*lX**2)-lX**4/48d0)

     & +(2d0*x/y-1d0)*(Li3(-x)-(lX+im*pi)*ddilog(-x)
     & +zeta3-lX**3/6d0-pisq/3d0*(lX+lY))

     & +2d0*(2d0*x/y+1d0)*(Li3(-y)+(lY+im*pi)*ddilog(-x)-zeta3
     & +0.25d0*lX*(2*lY**2+pisq)-0.125d0*lX**2*(lX+3d0*im*pi))

     & -0.25d0*(2d0*x**2-y**2)*((lX-lY)**2+pisq)

     & -0.25d0*(3d0+2d0*x/y**2)*((lX+im*pi)**2+pisq)
     & -(2d0-y**2)/(4d0*x**2)*((lY+im*pi)**2+pisq)+pisqo6

     & +0.5d0*(2d0*x+y**2)*((lX+im*pi)/y+(lY+im*pi)/x)-0.5d0
c--- end statement functions
      lX=log(-t/s)
      lY=log(-u/s)
      x=t/s
      y=u/s

c--- LEADING COLOR

C---Eq.(4.7)
      FL(2,2,2,2)=0.5d0
c--- by parity
      FL(1,1,1,1)=FL(2,2,2,2)
      
C---Eq.(4.8)
      FL(1,2,2,2)=FmpppL(x,y,lX,lY)+FmpppL(y,x,lY,lX)
c--- by 1<->2 exchange
      FL(2,1,2,2)=FL(1,2,2,2)

C---after Eq.(4.8)
      FL(2,2,1,2)=FppmpL(x,y,lX,lY)+FppmpL(y,x,lY,lX)
c--- by 3<->4 exchange
      FL(2,2,2,1)=FL(2,2,1,2)

c--- by parity
      FL(2,1,1,1)=FL(1,2,2,2)
      FL(1,2,1,1)=FL(2,1,2,2)
      FL(1,1,2,1)=FL(2,2,1,2)
      FL(1,1,1,2)=FL(2,2,2,1)

C---Eq.(4.9)
      FL(1,1,2,2)=
     & -(x**2+y**2)*(4d0*Li4(-x)+(lY-3d0*lX-2d0*im*pi)*Li3(-x)
     & +((lX+im*pi)**2+pisq)*ddilog(-x)+(lX+lY)**4/48d0
     & +im*pi/12d0*(lX+lY)**3+im*pi**3/2d0*lX
     & -pisq/12d0*lX**2-109d0/720d0*pi**4)

     & +0.5d0*x*(1d0-3d0*y)*(Li3(-x/y)-(lX-lY)*ddilog(-x/y)-zeta3
     & +0.5d0*lY*((lX-lY)**2+pisq))
     & +0.25d0*x**2*((lX-lY)**3+3d0*(lY+im*pi)*((lX-lY)**2+pisq))

     & +0.125d0*(14d0*(x-y)-8d0/y+9d0/y**2)*((lX+im*pi)**2+pisq)

     & +1d0/16d0*(38d0*x*y-13d0)*((lX-lY)**2+pisq)-pisq/12d0
     & -9d0/4d0*(1d0/y+2d0*x)*(lX+im*pi)+0.25d0

C---Eq.(4.10)
      FL(1,2,1,2)=FmpmpL(x,y,lX,lY)
c--- by 1<->2 exchange
      FL(2,1,1,2)=FmpmpL(y,x,lY,lX)

c--- by parity
      FL(2,2,1,1)=FL(1,1,2,2)
      FL(2,1,2,1)=FL(1,2,1,2)
      FL(1,2,2,1)=FL(2,1,1,2)


c--- SUBLEADING COLOR

C---Eq.(4.11)
      FSL(2,2,2,2)=-3d0/2d0
c--- by parity
      FSL(1,1,1,1)=FSL(2,2,2,2)

C---Eq.(4.12)
      FSL(1,2,2,2)=FmpppSL(x,y,lX,lY)+FmpppSL(y,x,lY,lX)
c--- by parity
      FSL(2,1,1,1)=FSL(1,2,2,2)



C---Eq.(4.13)
      FSL(2,2,1,2)=FSL(1,2,2,2)
      FSL(2,1,2,2)=FSL(1,2,2,2)
      FSL(2,2,2,1)=FSL(1,2,2,2)

c--- by parity
      FSL(1,1,2,1)=FSL(2,2,1,2)
      FSL(1,2,1,1)=FSL(2,1,2,2)
      FSL(1,1,1,2)=FSL(2,1,1,1)

C---Eq.(4.14)
      FSL(1,1,2,2)=FmmppSL(x,y,lX,lY)+FmmppSL(y,x,lY,lX)

C---Eq.(4.15)
      FSL(1,2,1,2)=FmpmpSL(x,y,lX,lY)

C---Eq.(4.16)
      FSL(2,1,1,2)=FmpmpSL(y,x,lY,lX)

c--- by parity
      FSL(2,2,1,1)=FSL(1,1,2,2)
      FSL(2,1,2,1)=FSL(1,2,1,2)
      FSL(1,2,2,1)=FSL(2,1,1,2)

      return
      end
