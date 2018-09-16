      subroutine twoloop(s,t,u,FL,FSL)
      implicit none
      include 'types.f'
C-----Formula taken from hep-ph/0109078v1
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      complex(dp):: FL(2,2,2,2),FSL(2,2,2,2)
      complex(dp):: FmpppL,FppmpL,FmpmpL,FmpppSL,FmmppSL,FmpmpSL
      real(dp):: s,t,u,x,y,lX,lY,ddilog,Li3,Li4,zeta3,zeta4
      parameter(zeta3=1.202056903159594_dp,zeta4=1.082323233711138_dp)

c--- statement functions
      FmpppL(x,y,lX,lY)=
     & ((2._dp+4._dp*x/y**2-5._dp*x**2/y**2)*((lX+im*pi)**2+pisq)
     & -(1._dp-x*y)*((lX-lY)**2+pisq)
     & +2._dp*(9._dp/y-10._dp*x)*(lX+im*pi))/8._dp

      FppmpL(x,y,lX,lY)=
     & ((2._dp+6._dp*x/y**2-3._dp*x**2/y**2)*((lX+im*pi)**2+pisq)
     & -(x-y)**2*((lX-lY)**2+pisq)
     & +2._dp*(9._dp/y-8._dp*x)*(lX+im*pi))/8._dp

      FmpmpL(x,y,lX,lY)=
     & -2._dp*(x**2+1._dp/y**2)*(Li4(-x)-zeta4
     & -0.5_dp*(lX+im*pi)*(Li3(-x)-zeta3)

     & +pisqo6*(ddilog(-x)-pisqo6-0.5_dp*lX**2)-(1._dp/48._dp)*lX**4

     & +1._dp/24._dp*(lX+im*pi)**2*((lX+im*pi)**2+pisq))
     & +2._dp*(3._dp*(1._dp-x)**2-2._dp)/y**2
     & *(Li4(-x)+Li4(-x/y)-Li4(-y)-(lY+im*pi)*(Li3(-x)-zeta3)
     & +pisqo6*(ddilog(-x)+0.5_dp*lY**2)
     & -lX*lY**3/6._dp+lY**4/24._dp-7._dp/360._dp*pi**4)
     & -2._dp/3._dp*(8._dp-x+30._dp*x/y)*(Li3(-y)-zeta3
     & -(lY+im*pi)*(ddilog(-y)-pisqo6)-0.5_dp*lX*((lY+im*pi)**2+pisq))
     & +(4._dp*y+27._dp+42._dp/y+4._dp/y**2)/6._dp
     & *(Li3(-x)-zeta3-(lX+im*pi)*(ddilog(-x)-pisqo6)
     & +im*pi/2._dp*lX**2-pisq*lX)
     & +(3._dp-2._dp/y-12._dp*x/y**2)/12._dp*(lX+im*pi)*((lX+im*pi)**2+pisq)
     & -y*(lX+im*pi)*((lY+im*pi)**2+pisq)/3._dp
     & +2._dp*(1._dp+2._dp/y)*(zeta3-pisqo6*(lY+im*pi))
     & +(y**2-24._dp*y+44._dp-8._dp*x**3/y)/24._dp*((lX-lY)**2+pisq)
     & -(15._dp-14._dp*x/y-48._dp*x/y**2)/24._dp*((lX+im*pi)**2+pisq)
     & +(8._dp*x/y+60._dp-24._dp*y/x+27._dp*y**2/x**2)/24._dp
     & *((lY+im*pi)**2+pisq)
     & +(4._dp/9._dp)*pisq*(x/y)
     & +(2._dp*x**2-54._dp*x-27._dp*y**2)/12._dp*((lX+im*pi)/y+(lY+im*pi)/x)

      FmpppSL(x,y,lX,lY)=
     & ((x**2+1._dp/y**2)*((lX+im*pi)**2+pisq)
     & +0.5_dp*(x**2+y**2)*((lX-lY)**2+pisq)
     & -4._dp*(1._dp/y-x)*(lX+im*pi))/8._dp

      FmmppSL(x,y,lX,lY)=
     & -2._dp*x**2*(Li4(-x)+Li4(-y)-(lX+im*pi)*(Li3(-x)+Li3(-y))

     & +lX**4/12._dp-lX**3*lY/3._dp+pisq/12._dp*lX*lY-4._dp/90._dp*pi**4
     & +im*pi/6._dp*lX*(lX**2-3._dp*lX*lY+pisq))
     & -(x-y)*(Li4(-x/y)-pisqo6*ddilog(-x))

     & -x*(2._dp*Li3(-x)-Li3(-x/y)-3._dp*zeta3-2._dp*(lX+im*pi)*ddilog(-x)

     & +(lX-lY)*(ddilog(-x/y)+lX**2)
     & +(5._dp*(lX-lY)+18._dp*im*pi)*((lX-lY)**2+pisq)/12._dp

     & -(2._dp/3._dp)*lX*(lX**2+pisq)-im*pi*(lY**2+pisq))

     & +(1._dp-2._dp*x**2)/(4._dp*y**2)*((lX+im*pi)**2+pisq)
     & -0.125_dp*(2._dp*x*y+3._dp)*((lX-lY)**2+pisq)+pisq/12._dp

     & +(0.5_dp/y+x)*(lX+im*pi)-0.25_dp

      FmpmpSL(x,y,lX,lY)=-2._dp*(x**2+1._dp)/y**2
     & *(Li4(-x/y)-Li4(-y)+0.5_dp*(lX-2._dp*lY-im*pi)*(Li3(-x)-zeta3)
     & +1._dp/24._dp*(lX**4+2._dp*im*pi*lX**3-4._dp*lX*lY**3+lY**4
     & +2._dp*pisq*lY**2)+7._dp/360._dp*pi**4)

     & -2._dp*(x-1._dp)/y*(Li4(-x)-zeta4-0.5_dp*(lX+im*pi)*(Li3(-x)-zeta3)
     & +pisqo6*(ddilog(-x)-pisqo6-0.5_dp*lX**2)-lX**4/48._dp)

     & +(2._dp*x/y-1._dp)*(Li3(-x)-(lX+im*pi)*ddilog(-x)
     & +zeta3-lX**3/6._dp-pisq/3._dp*(lX+lY))

     & +2._dp*(2._dp*x/y+1._dp)*(Li3(-y)+(lY+im*pi)*ddilog(-x)-zeta3
     & +0.25_dp*lX*(2*lY**2+pisq)-0.125_dp*lX**2*(lX+3._dp*im*pi))

     & -0.25_dp*(2._dp*x**2-y**2)*((lX-lY)**2+pisq)

     & -0.25_dp*(3._dp+2._dp*x/y**2)*((lX+im*pi)**2+pisq)
     & -(2._dp-y**2)/(4._dp*x**2)*((lY+im*pi)**2+pisq)+pisqo6

     & +0.5_dp*(2._dp*x+y**2)*((lX+im*pi)/y+(lY+im*pi)/x)-0.5_dp
c--- end statement functions
      lX=log(-t/s)
      lY=log(-u/s)
      x=t/s
      y=u/s

c--- LEADING COLOR

C---Eq.(4.7)
      FL(2,2,2,2)=0.5_dp
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
     & -(x**2+y**2)*(4._dp*Li4(-x)+(lY-3._dp*lX-2._dp*im*pi)*Li3(-x)
     & +((lX+im*pi)**2+pisq)*ddilog(-x)+(lX+lY)**4/48._dp
     & +im*pi/12._dp*(lX+lY)**3+im*pi**3/2._dp*lX
     & -pisq/12._dp*lX**2-109._dp/720._dp*pi**4)

     & +0.5_dp*x*(1._dp-3._dp*y)*(Li3(-x/y)-(lX-lY)*ddilog(-x/y)-zeta3
     & +0.5_dp*lY*((lX-lY)**2+pisq))
     & +0.25_dp*x**2*((lX-lY)**3+3._dp*(lY+im*pi)*((lX-lY)**2+pisq))

     & +0.125_dp*(14._dp*(x-y)-8._dp/y+9._dp/y**2)*((lX+im*pi)**2+pisq)

     & +1._dp/16._dp*(38._dp*x*y-13._dp)*((lX-lY)**2+pisq)-pisq/12._dp
     & -9._dp/4._dp*(1._dp/y+2._dp*x)*(lX+im*pi)+0.25_dp

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
      FSL(2,2,2,2)=-3._dp/2._dp
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
