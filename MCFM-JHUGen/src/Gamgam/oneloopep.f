      subroutine oneloopep(s,t,u,M1)
      implicit none
      include 'types.f'
c--- Formulae taken from hep-ph/0109078
c--- Routine written automatically by FORM program writef.frm
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'epinv.f'
      include 'nflav.f'
      include 'scale.f'
      integer:: h1,h2,h3,h4
      real(dp):: s,t,u,lX,lY,ddilog,Li3,Li4,
     & Li2mx,Li2my,Li3mx,Li3my,Li4mx,Li4my,Li4mxoy,Li4myox,Lns
      complex(dp):: Amp(2,2,2,2,0:2),M1(2,2,2,2)
      complex(dp):: Apppp0,Apppp1,Apppp2,Amppp0,Amppp1,Amppp2
      complex(dp):: Ammpp0,Ammpp1,Ammpp2,Ampmp0,Ampmp1,Ampmp2
      complex(dp):: xlog,lnrat,Igg_gaga0,Igg_gaga1,Igg_gaga2
      real(dp),parameter::zeta3=1.202056903159594_dp,
     &                    zeta4=1.082323233711138_dp

c--- statement functions for amplitudes
      Apppp0(t,u,lX,lY,Li2mx,Li2my,Li3mx,Li3my,Li4mx,Li4my,Li4mxoy)=
     & 1._dp

      Apppp1(t,u,lX,lY,Li2mx,Li2my,Li3mx,Li3my,Li4mx,Li4my,Li4mxoy)=
     & 8._dp/3._dp-1._dp/6._dp*s**(-2)*t*u*pisq
     & -1._dp/6._dp*s**(-2)*t*u*lY**2
     & +1._dp/3._dp*s**(-2)*t*u*lX*lY-1._dp/6._dp*s**(-2)*t*u*lX**2+1.
     & _dp/6._dp*s**(-1)*u*lY-1._dp/6._dp*s**(-1)*u*lX-1._dp/6._dp*
     & s**(-1)*t*lY+1._dp/6._dp*s**(-1)*t*lX+1._dp/3._dp*t**(-1)*u*lY
     & -Lns-1._dp/6._dp*lY-1._dp/6._dp*lX+1._dp/3._dp*t*u**(-1)*lX-
     & 1._dp/6._dp*s*t**(-2)*u*lY**2-1._dp/3._dp*s*t**(-2)*u*pi*lY*im
     & -1._dp/6._dp*s*t*u**(-2)*lX**2-1._dp/3._dp*s*t*u**(-2)*pi*lX*im
     & +1._dp/3._dp*s**2*t**(-1)*u**(-1)*pi*im

      Apppp2(t,u,lX,lY,Li2mx,Li2my,Li3mx,Li3my,Li4mx,Li4my,Li4mxoy)=
     & 52._dp/9._dp-1._dp/3._dp*s**(-2)*t*u*Li3mx-1._dp/3._dp*s**(-2)*t*
     & u*Li3my+1._dp/6._dp*s**(-2)*t*u*pisq*Lns+1._dp/3._dp*s**(-2)*t*u*
     & lY*Li2my+1._dp/6._dp*s**(-2)*t*u*lY*pisq+1._dp/6._dp*s**(-2)*t*u
     & *lY**2*Lns+1._dp/9._dp*s**(-2)*t*u*lY**3+1._dp/3._dp*s**(-2)*t*u
     & *lX*Li2mx+1._dp/6._dp*s**(-2)*t*u*lX*pisq+5._dp/9._dp*s**(-2)*t*
     & u*lX*lY-1._dp/3._dp*s**(-2)*t*u*lX*lY*Lns+1._dp/6._dp*s**(-2)*t*
     & u*lX**2*Lns+1._dp/9._dp*s**(-2)*t*u*lX**3+5._dp/18._dp*s**(-2)*
     & t**2*pisq+5._dp/18._dp*s**(-2)*t**2*lY**2+5._dp/18._dp*s**(-2)*
     & t**2*lX**2-8._dp/9._dp*s**(-1)*t**(-1)*u**2*lY-1._dp/6._dp*
     & s**(-1)*u*lY*Lns+1._dp/6._dp*s**(-1)*u*lX*Lns+5._dp/18._dp*
     & s**(-1)*t*pisq+1._dp/6._dp*s**(-1)*t*lY*Lns+4._dp/9._dp*s**(-1)*
     & t*lY**2+8._dp/9._dp*s**(-1)*t*lX-1._dp/6._dp*s**(-1)*t*lX*Lns-
     & 1._dp/9._dp*s**(-1)*t**2*u**(-1)*lX**2-1._dp/3._dp*t**(-1)*u*lY*
     & Lns-8._dp/3._dp*Lns+1._dp/2._dp*Lns**2-1._dp/6._dp*pisq+1._dp/6.
     & _dp*lY*Lns+1._dp/3._dp*lY**2+1._dp/6._dp*lX*Lns-1._dp/6._dp*t*
     & u**(-1)*pisq+8._dp/9._dp*t*u**(-1)*lX-1._dp/3._dp*t*u**(-1)*lX*
     & Lns+5._dp/18._dp*t**2*u**(-2)*lX**2+1._dp/3._dp*s*t**(-2)*u*
     & Li3my-1._dp/3._dp*s*t**(-2)*u*lY*Li2my+1._dp/6._dp*s*t**(-2)*u*
     & lY*pisq+1._dp/6._dp*s*t**(-2)*u*lY**2*Lns+1._dp/9._dp*s*t**(-2)*
     & u*lY**3-1._dp/6._dp*s*t**(-2)*u*lX*lY**2-1._dp/3._dp*s*t**(-2)*u
     & *zeta3-1._dp/3._dp*s*t**(-2)*u*pi*im*Li2my+1._dp/18._dp*s*
     & t**(-2)*u*pi*pisq*im-5._dp/9._dp*s*t**(-2)*u*pi*lY*im+1._dp/3._dp
     & *s*t**(-2)*u*pi*lY*im*Lns+1._dp/6._dp*s*t**(-2)*u*pi*lY**2*im-
     & 1._dp/3._dp*s*t**(-2)*u*pi*lX*lY*im+1._dp/6._dp*s*t**(-1)*pisq+4.
     & _dp/9._dp*s*t**(-1)*lY**2+1._dp/3._dp*s*t*u**(-2)*Li3mx-1._dp/3._dp
     & *s*t*u**(-2)*lX*Li2mx+1._dp/6._dp*s*t*u**(-2)*lX*pisq+1._dp/6._dp
     & *s*t*u**(-2)*lX**2*Lns-1._dp/6._dp*s*t*u**(-2)*lX**2*lY+1._dp/9.
     & _dp*s*t*u**(-2)*lX**3-1._dp/3._dp*s*t*u**(-2)*zeta3-1._dp/3._dp*s
     & *t*u**(-2)*pi*im*Li2mx+1._dp/18._dp*s*t*u**(-2)*pi*pisq*im-5._dp
     & /9._dp*s*t*u**(-2)*pi*lX*im+1._dp/3._dp*s*t*u**(-2)*pi*lX*im*Lns
     & -1._dp/3._dp*s*t*u**(-2)*pi*lX*lY*im+1._dp/6._dp*s*t*u**(-2)*pi*
     & lX**2*im+5._dp/18._dp*s**2*t**(-2)*lY**2+8._dp/9._dp*s**2*
     & t**(-1)*u**(-1)*pi*im-1._dp/3._dp*s**2*t**(-1)*u**(-1)*pi*im*Lns


      Amppp0(t,u,lX,lY,Li2mx,Li2my,Li3mx,Li3my,Li4mx,Li4my,Li4mxoy)=
     & 1._dp

      Amppp1(t,u,lX,lY,Li2mx,Li2my,Li3mx,Li3my,Li4mx,Li4my,Li4mxoy)=8._dp
     &  /3._dp+1._dp/3._dp*s**(-2)*t*u*pisq+1._dp/3._dp*s**(-2)*t*u*lY**2
     & -2._dp/3._dp*s**(-2)*t*u*lX*lY+1._dp/3._dp*s**(-2)*t*u*lX**2-1.
     & _dp/3._dp*s**(-1)*u*lY+1._dp/3._dp*s**(-1)*u*lX+1._dp/3._dp*
     & s**(-1)*t*lY-1._dp/3._dp*s**(-1)*t*lX-Lns+pi*im+1._dp/3._dp*
     & s*t**(-2)*u*lY**2+2._dp/3._dp*s*t**(-2)*u*pi*lY*im+2._dp/3._dp*s
     & *t**(-1)*lY+2._dp/3._dp*s*u**(-1)*lX+1._dp/3._dp*s*t*u**(-2)*
     & lX**2+2._dp/3._dp*s*t*u**(-2)*pi*lX*im-2._dp/3._dp*s**2*t**(-1)*
     & u**(-1)*pi*im

      Amppp2(t,u,lX,lY,Li2mx,Li2my,Li3mx,Li3my,Li4mx,Li4my,Li4mxoy)=52._dp
     & /9._dp+2._dp/3._dp*s**(-2)*t*u*Li3mx+2._dp/3._dp*s**(-2)*t*u*
     & Li3my+13._dp/18._dp*s**(-2)*t*u*pisq-1._dp/3._dp*s**(-2)*t*u*
     & pisq*Lns-2._dp/3._dp*s**(-2)*t*u*lY*Li2my-1._dp/3._dp*s**(-2)*t*
     & u*lY*pisq+13._dp/18._dp*s**(-2)*t*u*lY**2-1._dp/3._dp*s**(-2)*t*
     & u*lY**2*Lns-2._dp/9._dp*s**(-2)*t*u*lY**3-2._dp/3._dp*s**(-2)*t*
     & u*lX*Li2mx-1._dp/3._dp*s**(-2)*t*u*lX*pisq+2._dp/3._dp*s**(-2)*t
     & *u*lX*lY*Lns-1._dp/3._dp*s**(-2)*t*u*lX**2*Lns-2._dp/9._dp*
     & s**(-2)*t*u*lX**3+13._dp/9._dp*s**(-2)*t**2*lX*lY-13._dp/18._dp*
     & s**(-2)*t**2*lX**2+1._dp/3._dp*s**(-1)*u*lY*Lns-1._dp/3._dp*
     & s**(-1)*u*lX*Lns+19._dp/9._dp*s**(-1)*t*lY-1._dp/3._dp*s**(-1)*t
     & *lY*Lns-1._dp/3._dp*s**(-1)*t*lY**2-19._dp/9._dp*s**(-1)*t*lX+
     & 1._dp/3._dp*s**(-1)*t*lX*Lns+13._dp/9._dp*s**(-1)*t*lX*lY-7._dp/
     & 18._dp*s**(-1)*t*lX**2+1._dp/3._dp*t**(-1)*u*pisq-8._dp/3._dp*Lns
     & +1._dp/2._dp*Lns**2+1._dp/6._dp*pisq+11._dp/9._dp*lY-1._dp/6._dp
     & *lY**2-3._dp*lX+1._dp/2._dp*lX**2+8._dp/9._dp*pi*im-pi*im*Lns
     & +1._dp/3._dp*t*u**(-1)*pisq-19._dp/9._dp*t*u**(-1)*lX-7._dp/18._dp
     & *t*u**(-1)*lX**2-19._dp/9._dp*t*u**(-1)*pi*im+1._dp/3._dp*t*
     & u**(-1)*pi*lX**2*im-13._dp/18._dp*t**2*u**(-2)*lX**2+1._dp/3._dp
     & *t**2*u**(-2)*pi*lX**2*im-2._dp/3._dp*s*t**(-2)*u*Li3my+2._dp/3.
     & _dp*s*t**(-2)*u*lY*Li2my-1._dp/3._dp*s*t**(-2)*u*lY*pisq+13._dp/
     & 18._dp*s*t**(-2)*u*lY**2-1._dp/3._dp*s*t**(-2)*u*lY**2*Lns-2._dp/
     & 9._dp*s*t**(-2)*u*lY**3+1._dp/3._dp*s*t**(-2)*u*lX*lY**2+2._dp/3.
     & _dp*s*t**(-2)*u*zeta3-1._dp/9._dp*s*t**(-2)*u*pi*pisq*im+13._dp/
     & 9._dp*s*t**(-2)*u*pi*lY*im-2._dp/3._dp*s*t**(-2)*u*pi*lY*im*Lns
     & -1._dp/3._dp*s*t**(-2)*u*pi*lY**2*im+2._dp/3._dp*s*t**(-2)*u*pi*
     & lX*lY*im+19._dp/9._dp*s*t**(-1)*lY-2._dp/3._dp*s*t**(-1)*lY*Lns
     & -1._dp/3._dp*s*t**(-1)*lY**2+19._dp/9._dp*s*t**(-1)*pi*im-2._dp/
     & 3._dp*s*t**(-1)*pi*im*Li2my-2._dp/3._dp*s*u**(-1)*lX*Lns-2._dp/3.
     & _dp*s*t*u**(-2)*Li3mx+2._dp/3._dp*s*t*u**(-2)*lX*Li2mx-1._dp/3._dp
     & *s*t*u**(-2)*lX*pisq-1._dp/3._dp*s*t*u**(-2)*lX**2*Lns+1._dp/3._dp
     & *s*t*u**(-2)*lX**2*lY-2._dp/9._dp*s*t*u**(-2)*lX**3+2._dp/3._dp
     & *s*t*u**(-2)*zeta3+2._dp/3._dp*s*t*u**(-2)*pi*im*Li2mx-1._dp/9._dp
     & *s*t*u**(-2)*pi*pisq*im+13._dp/9._dp*s*t*u**(-2)*pi*lX*im-2._dp
     & /3._dp*s*t*u**(-2)*pi*lX*im*Lns+2._dp/3._dp*s*t*u**(-2)*pi*lX*lY*
     & im-2._dp/3._dp*s**2*t**(-2)*pi*im*Li2my+2._dp/3._dp*s**2*t**(-1)
     & *u**(-1)*pi*im*Lns


      Ammpp0(t,u,lX,lY,Li2mx,Li2my,Li3mx,Li3my,Li4mx,Li4my,Li4mxoy)=-
     & 1._dp-1._dp/2._dp*s**(-2)*u**2*pisq-1._dp/2._dp*s**(-2)*u**2*
     & lY**2+s**(-2)*u**2*lX*lY-1._dp/2._dp*s**(-2)*u**2*lX**2-1._dp/
     & 2._dp*s**(-2)*t**2*pisq-1._dp/2._dp*s**(-2)*t**2*lY**2+s**(-2)*
     & t**2*lX*lY-1._dp/2._dp*s**(-2)*t**2*lX**2-s**(-1)*u*lY+
     & s**(-1)*u*lX+s**(-1)*t*lY-s**(-1)*t*lX

      Ammpp1(t,u,lX,lY,Li2mx,Li2my,Li3mx,Li3my,Li4mx,Li4my,Li4mxoy)=-
     & 10._dp/3._dp-s**(-2)*u**2*Li3mx-s**(-2)*u**2*Li3my+1._dp/2._dp
     & *s**(-2)*u**2*pisq*Lns+s**(-2)*u**2*lY*Li2my+1._dp/2._dp*
     & s**(-2)*u**2*lY*pisq+1._dp/2._dp*s**(-2)*u**2*lY**2*Lns+1._dp/3.
     & _dp*s**(-2)*u**2*lY**3+s**(-2)*u**2*lX*Li2mx+1._dp/2._dp*
     & s**(-2)*u**2*lX*pisq-s**(-2)*u**2*lX*lY*Lns+1._dp/2._dp*
     & s**(-2)*u**2*lX**2*Lns+1._dp/3._dp*s**(-2)*u**2*lX**3+11._dp/6._dp
     & *s**(-2)*t*u*pisq+11._dp/6._dp*s**(-2)*t*u*lY**2-11._dp/3._dp*
     & s**(-2)*t*u*lX*lY+11._dp/6._dp*s**(-2)*t*u*lX**2-s**(-2)*t**2*
     & Li3mx-s**(-2)*t**2*Li3my+1._dp/2._dp*s**(-2)*t**2*pisq*Lns+
     & s**(-2)*t**2*lY*Li2my+1._dp/2._dp*s**(-2)*t**2*lY*pisq+1._dp/2._dp
     & *s**(-2)*t**2*lY**2*Lns+1._dp/3._dp*s**(-2)*t**2*lY**3+
     & s**(-2)*t**2*lX*Li2mx+1._dp/2._dp*s**(-2)*t**2*lX*pisq-s**(-2)
     & *t**2*lX*lY*Lns+1._dp/2._dp*s**(-2)*t**2*lX**2*Lns+1._dp/3._dp*
     & s**(-2)*t**2*lX**3-11._dp/3._dp*s**(-1)*u*lY+s**(-1)*u*lY*Lns
     & +s**(-1)*u*lY**2+2._dp*s**(-1)*u*lX-s**(-1)*u*lX*Lns+2._dp
     & *s**(-1)*t*lY-s**(-1)*t*lY*Lns-11._dp/3._dp*s**(-1)*t*lX+
     & s**(-1)*t*lX*Lns+s**(-1)*t*lX**2+1._dp/3._dp*t**(-1)*u*lY+
     & Lns-1._dp/2._dp*pisq+lX*lY+1._dp/3._dp*t*u**(-1)*lX-1._dp/3._dp
     & *s*t**(-1)*lY**2-2._dp/3._dp*s*t**(-1)*pi*lY*im-1._dp/3._dp*s*
     & u**(-1)*lX**2-2._dp/3._dp*s*u**(-1)*pi*lX*im+1._dp/6._dp*s**2*
     & t**(-2)*lY**2+1._dp/3._dp*s**2*t**(-2)*pi*lY*im+1._dp/3._dp*s**2
     & *t**(-1)*u**(-1)*pi*im+1._dp/6._dp*s**2*u**(-2)*lX**2+1._dp/3._dp
     & *s**2*u**(-2)*pi*lX*im

      Ammpp2(t,u,lX,lY,Li2mx,Li2my,Li3mx,Li3my,Li4mx,Li4my,Li4mxoy)=-
     & 56._dp/9._dp+s**(-2)*u**2*Li3mx*Lns+s**(-2)*u**2*Li3my*Lns-1.
     & _dp/4._dp*s**(-2)*u**2*pisq*Lns**2-s**(-2)*u**2*lY*Li2my*Lns-1.
     & _dp/2._dp*s**(-2)*u**2*lY*pisq*Lns-1._dp/4._dp*s**(-2)*u**2*lY**2*
     & Lns**2-1._dp/3._dp*s**(-2)*u**2*lY**3*Lns-s**(-2)*u**2*lX*
     & Li2mx*Lns-1._dp/2._dp*s**(-2)*u**2*lX*pisq*Lns+1._dp/2._dp*
     & s**(-2)*u**2*lX*lY*Lns**2-1._dp/4._dp*s**(-2)*u**2*lX**2*Lns**2
     & -1._dp/3._dp*s**(-2)*u**2*lX**3*Lns-11._dp/6._dp*s**(-2)*t*u*
     & pisq*Lns-11._dp/6._dp*s**(-2)*t*u*lY**2*Lns+11._dp/3._dp*s**(-2)
     & *t*u*lX*lY*Lns-11._dp/6._dp*s**(-2)*t*u*lX**2*Lns+2._dp*s**(-2)
     & *t**2*Li4mx+2._dp*s**(-2)*t**2*Li4my-11._dp/3._dp*s**(-2)*t**2*
     & Li3mx+s**(-2)*t**2*Li3mx*Lns-11._dp/3._dp*s**(-2)*t**2*Li3my
     & +s**(-2)*t**2*Li3my*Lns-67._dp/18._dp*s**(-2)*t**2*pisq-1._dp/
     & 4._dp*s**(-2)*t**2*pisq*Lns**2+2._dp*s**(-2)*t**2*lY*Li3mx+11._dp
     & /3._dp*s**(-2)*t**2*lY*Li2my-s**(-2)*t**2*lY*Li2my*Lns+11._dp/
     & 6._dp*s**(-2)*t**2*lY*pisq-1._dp/2._dp*s**(-2)*t**2*lY*pisq*Lns
     & -67._dp/18._dp*s**(-2)*t**2*lY**2-1._dp/4._dp*s**(-2)*t**2*lY**2
     & *Lns**2-s**(-2)*t**2*lY**2*Li2my-1._dp/2._dp*s**(-2)*t**2*
     & lY**2*pisq+11._dp/9._dp*s**(-2)*t**2*lY**3-1._dp/3._dp*s**(-2)*
     & t**2*lY**3*Lns-1._dp/4._dp*s**(-2)*t**2*lY**4+2._dp*s**(-2)*
     & t**2*lX*Li3my+11._dp/3._dp*s**(-2)*t**2*lX*Li2mx-s**(-2)*t**2*
     & lX*Li2mx*Lns+11._dp/6._dp*s**(-2)*t**2*lX*pisq-1._dp/2._dp*
     & s**(-2)*t**2*lX*pisq*Lns+67._dp/9._dp*s**(-2)*t**2*lX*lY+1._dp/
     & 2._dp*s**(-2)*t**2*lX*lY*Lns**2-4._dp/3._dp*s**(-2)*t**2*lX*lY*
     & pisq-1._dp/3._dp*s**(-2)*t**2*lX*lY**3-67._dp/18._dp*s**(-2)*
     & t**2*lX**2-1._dp/4._dp*s**(-2)*t**2*lX**2*Lns**2-s**(-2)*t**2*
     & lX**2*Li2mx-1._dp/2._dp*s**(-2)*t**2*lX**2*pisq+s**(-2)*t**2*
     & lX**2*lY**2+11._dp/9._dp*s**(-2)*t**2*lX**3-1._dp/3._dp*s**(-2)*
     & t**2*lX**3*Lns-1._dp/3._dp*s**(-2)*t**2*lX**3*lY-1._dp/4._dp*
     & s**(-2)*t**2*lX**4-8._dp*s**(-2)*t**2*zeta4+11._dp/3._dp*
     & s**(-1)*u*lY*Lns-1._dp/2._dp*s**(-1)*u*lY*Lns**2-s**(-1)*u*
     & lY**2*Lns-2._dp*s**(-1)*u*lX*Lns+1._dp/2._dp*s**(-1)*u*lX*
     & Lns**2+2._dp*s**(-1)*t*Li4mx+2._dp*s**(-1)*t*Li4my-11._dp/3._dp
     & *s**(-1)*t*Li3mx-11._dp/3._dp*s**(-1)*t*Li3my-67._dp/18._dp*
     & s**(-1)*t*pisq+118._dp/9._dp*s**(-1)*t*lY-2._dp*s**(-1)*t*lY*
     & Lns-1._dp/2._dp*s**(-1)*t*lY*Lns**2+2._dp*s**(-1)*t*lY*Li3mx+
     & 11._dp/3._dp*s**(-1)*t*lY*Li2my+11._dp/6._dp*s**(-1)*t*lY*pisq-
     & 59._dp/9._dp*s**(-1)*t*lY**2-s**(-1)*t*lY**2*Lns-s**(-1)*t*
     & lY**2*Li2my-1._dp/2._dp*s**(-1)*t*lY**2*pisq+11._dp/9._dp*
     & s**(-1)*t*lY**3-1._dp/4._dp*s**(-1)*t*lY**4-118._dp/9._dp*
     & s**(-1)*t*lX+11._dp/3._dp*s**(-1)*t*lX*Lns+1._dp/2._dp*s**(-1)*t
     & *lX*Lns**2+2._dp*s**(-1)*t*lX*Li3my+11._dp/3._dp*s**(-1)*t*lX*
     & Li2mx+11._dp/6._dp*s**(-1)*t*lX*pisq+67._dp/9._dp*s**(-1)*t*lX*
     & lY-4._dp/3._dp*s**(-1)*t*lX*lY*pisq-1._dp/3._dp*s**(-1)*t*lX*
     & lY**3-8._dp/9._dp*s**(-1)*t*lX**2-s**(-1)*t*lX**2*Li2mx-1._dp/
     & 2._dp*s**(-1)*t*lX**2*pisq+s**(-1)*t*lX**2*lY**2+11._dp/9._dp*
     & s**(-1)*t*lX**3-1._dp/3._dp*s**(-1)*t*lX**3*lY-1._dp/4._dp*
     & s**(-1)*t*lX**4-8._dp*s**(-1)*t*zeta4-1._dp/3._dp*t**(-1)*u*lY*
     & Lns+10._dp/3._dp*Lns-1._dp/2._dp*Lns**2+Li4mx+Li4my-2._dp*
     & Li3mx-Li3my-7._dp/6._dp*pisq+1._dp/2._dp*pisq*Lns+74._dp/9._dp
     & *lY-1._dp/2._dp*lY*Lns**2+lY*Li3mx+lY*Li2my+1._dp/2._dp*lY*
     & pisq-8._dp/3._dp*lY**2-1._dp/2._dp*lY**2*Lns-1._dp/2._dp*lY**2*
     & Li2my-1._dp/4._dp*lY**2*pisq+1._dp/3._dp*lY**3-1._dp/8._dp*lY**4
     & -4._dp*lX+1._dp/2._dp*lX*Lns**2+lX*Li3my+2._dp*lX*Li2mx+2._dp
     & *lX*lY-lX*lY*Lns-2._dp/3._dp*lX*lY*pisq-1._dp/6._dp*lX*lY**3
     & +lX**2+1._dp/2._dp*lX**2*Lns-1._dp/2._dp*lX**2*Li2mx-1._dp/4._dp
     & *lX**2*pisq+1._dp/2._dp*lX**2*lY+1._dp/2._dp*lX**2*lY**2-1._dp/
     & 6._dp*lX**3*lY-1._dp/8._dp*lX**4-4._dp*zeta4+zeta3+8._dp/9._dp
     & *pi*im+pi*im*Li2mx-1._dp/6._dp*pi*pisq*im+2._dp*pi*lX*im+pi
     & *lX*lY*im-1._dp/2._dp*pi*lX**2*im-4._dp/3._dp*t*u**(-1)*Li3mx-
     & 1._dp/6._dp*t*u**(-1)*pisq+8._dp/9._dp*t*u**(-1)*lX-1._dp/3._dp*t*
     & u**(-1)*lX*Lns+4._dp/3._dp*t*u**(-1)*lX*Li2mx-2._dp/3._dp*t*
     & u**(-1)*lX*pisq+10._dp/9._dp*t*u**(-1)*lX**2+2._dp/3._dp*t*
     & u**(-1)*lX**2*lY-4._dp/9._dp*t*u**(-1)*lX**3+4._dp/3._dp*t*
     & u**(-1)*zeta3+8._dp/9._dp*t*u**(-1)*pi*im+4._dp/3._dp*t*u**(-1)*
     & pi*im*Li2mx-2._dp/9._dp*t*u**(-1)*pi*pisq*im+23._dp/9._dp*t*
     & u**(-1)*pi*lX*im+4._dp/3._dp*t*u**(-1)*pi*lX*lY*im-2._dp/3._dp*t
     & *u**(-1)*pi*lX**2*im-1._dp/3._dp*t**2*u**(-2)*Li3mx+1._dp/3._dp*
     & t**2*u**(-2)*lX*Li2mx-1._dp/6._dp*t**2*u**(-2)*lX*pisq+5._dp/18.
     & _dp*t**2*u**(-2)*lX**2+1._dp/6._dp*t**2*u**(-2)*lX**2*lY-1._dp/9.
     & _dp*t**2*u**(-2)*lX**3+1._dp/3._dp*t**2*u**(-2)*zeta3+1._dp/3._dp
     & *t**2*u**(-2)*pi*im*Li2mx-1._dp/18._dp*t**2*u**(-2)*pi*pisq*im
     & +5._dp/9._dp*t**2*u**(-2)*pi*lX*im+1._dp/3._dp*t**2*u**(-2)*pi*
     & lX*lY*im-1._dp/6._dp*t**2*u**(-2)*pi*lX**2*im+2._dp/3._dp*s*
     & t**(-1)*Li3my+1._dp/6._dp*s*t**(-1)*pisq-8._dp/9._dp*s*t**(-1)*
     & lY-2._dp/3._dp*s*t**(-1)*lY*Li2my+1._dp/3._dp*s*t**(-1)*lY*pisq
     & -5._dp/9._dp*s*t**(-1)*lY**2+1._dp/3._dp*s*t**(-1)*lY**2*Lns+2.
     & _dp/9._dp*s*t**(-1)*lY**3-1._dp/3._dp*s*t**(-1)*lX*lY**2-2._dp/3._dp
     & *s*t**(-1)*zeta3-8._dp/9._dp*s*t**(-1)*pi*im-2._dp/3._dp*s*
     & t**(-1)*pi*im*Li2my+1._dp/9._dp*s*t**(-1)*pi*pisq*im-13._dp/9._dp
     & *s*t**(-1)*pi*lY*im+2._dp/3._dp*s*t**(-1)*pi*lY*im*Lns+1._dp/3._dp
     & *s*t**(-1)*pi*lY**2*im-2._dp/3._dp*s*t**(-1)*pi*lX*lY*im+1._dp/
     & 3._dp*s*u**(-1)*lX**2*Lns+2._dp/3._dp*s*u**(-1)*pi*lX*im*Lns-1._dp
     & /3._dp*s**2*t**(-2)*Li3my+1._dp/3._dp*s**2*t**(-2)*lY*Li2my-1._dp
     & /6._dp*s**2*t**(-2)*lY*pisq+5._dp/18._dp*s**2*t**(-2)*lY**2-1._dp
     & /6._dp*s**2*t**(-2)*lY**2*Lns-1._dp/9._dp*s**2*t**(-2)*lY**3+1.
     & _dp/6._dp*s**2*t**(-2)*lX*lY**2+1._dp/3._dp*s**2*t**(-2)*zeta3+1.
     & _dp/3._dp*s**2*t**(-2)*pi*im*Li2my-1._dp/18._dp*s**2*t**(-2)*pi*
     & pisq*im+5._dp/9._dp*s**2*t**(-2)*pi*lY*im-1._dp/3._dp*s**2*
     & t**(-2)*pi*lY*im*Lns-1._dp/6._dp*s**2*t**(-2)*pi*lY**2*im+1._dp/
     & 3._dp*s**2*t**(-2)*pi*lX*lY*im-1._dp/3._dp*s**2*t**(-1)*u**(-1)*
     & pi*im*Lns-1._dp/6._dp*s**2*u**(-2)*lX**2*Lns-1._dp/3._dp*s**2*
     & u**(-2)*pi*lX*im*Lns


      Ampmp0(t,u,lX,lY,Li2mx,Li2my,Li3mx,Li3my,Li4mx,Li4my,Li4mxoy)=-
     & 1._dp-t*u**(-1)*lX-t*u**(-1)*pi*im-1._dp/2._dp*t**2*u**(-2)*
     & lX**2-t**2*u**(-2)*pi*lX*im+s*u**(-1)*lX+s*u**(-1)*pi*im
     & -1._dp/2._dp*s**2*u**(-2)*lX**2-s**2*u**(-2)*pi*lX*im

      Ampmp1(t,u,lX,lY,Li2mx,Li2my,Li3mx,Li3my,Li4mx,Li4my,Li4mxoy)=-
     & 10._dp/3._dp-1._dp/6._dp*s**(-2)*t*u*pisq-1._dp/6._dp*s**(-2)*t*u*
     & lY**2+1._dp/3._dp*s**(-2)*t*u*lX*lY-1._dp/6._dp*s**(-2)*t*u*
     & lX**2-1._dp/3._dp*s**(-1)*t**(-1)*u**2*lY+1._dp/2._dp*s**(-1)*
     & t**(-1)*u**2*lY**2+s**(-1)*u*lX*lY-1._dp/2._dp*s**(-1)*u*lX**2
     & +1._dp/2._dp*s**(-1)*t*pisq+1._dp/3._dp*s**(-1)*t*lX-t**(-1)*u
     & *pi*lY*im+Lns-2._dp*lX-11._dp/3._dp*pi*im-pi*lX*im-17._dp/
     & 3._dp*t*u**(-1)*lX+t*u**(-1)*lX*Lns-t*u**(-1)*lX**2-17._dp/3.
     & _dp*t*u**(-1)*pi*im+t*u**(-1)*pi*im*Lns+t**2*u**(-2)*Li3mx-
     & t**2*u**(-2)*lX*Li2mx+1._dp/2._dp*t**2*u**(-2)*lX*pisq-2._dp*
     & t**2*u**(-2)*lX**2+1._dp/2._dp*t**2*u**(-2)*lX**2*Lns-1._dp/2._dp
     & *t**2*u**(-2)*lX**2*lY+1._dp/3._dp*t**2*u**(-2)*lX**3-t**2*
     & u**(-2)*zeta3-t**2*u**(-2)*pi*im*Li2mx+1._dp/6._dp*t**2*
     & u**(-2)*pi*pisq*im+t**2*u**(-2)*pi*lX*im*Lns-t**2*u**(-2)*pi
     & *lX*lY*im+1._dp/2._dp*t**2*u**(-2)*pi*lX**2*im-1._dp/6._dp*s*
     & t**(-2)*u*lY**2-1._dp/3._dp*s*t**(-2)*u*pi*lY*im-1._dp/3._dp*s*
     & t**(-1)*pi*im-s*u**(-1)*pisq-s*u**(-1)*lX*Lns-s*u**(-1)*pi
     & *im*Lns-1._dp/6._dp*s*t*u**(-2)*lX**2+11._dp/3._dp*s*t*u**(-2)*
     & pi*lX*im+s**2*u**(-2)*Li3mx-s**2*u**(-2)*lX*Li2mx+1._dp/2._dp
     & *s**2*u**(-2)*lX*pisq+1._dp/2._dp*s**2*u**(-2)*lX**2*Lns-1._dp/
     & 2._dp*s**2*u**(-2)*lX**2*lY+1._dp/3._dp*s**2*u**(-2)*lX**3-s**2
     & *u**(-2)*zeta3-s**2*u**(-2)*pi*im*Li2mx+1._dp/6._dp*s**2*
     & u**(-2)*pi*pisq*im+s**2*u**(-2)*pi*lX*im*Lns-s**2*u**(-2)*pi
     & *lX*lY*im+1._dp/2._dp*s**2*u**(-2)*pi*lX**2*im

      Ampmp2(t,u,lX,lY,Li2mx,Li2my,Li3mx,Li3my,Li4mx,Li4my,Li4mxoy)=-
     & 56._dp/9._dp+1._dp/6._dp*s**(-2)*t*u*pisq*Lns+1._dp/6._dp*s**(-2)*
     & t*u*lY**2*Lns-1._dp/3._dp*s**(-2)*t*u*lX*lY*Lns+1._dp/6._dp*
     & s**(-2)*t*u*lX**2*Lns+1._dp/3._dp*s**(-2)*t**2*Li3mx+1._dp/3._dp
     & *s**(-2)*t**2*Li3my+5._dp/18._dp*s**(-2)*t**2*pisq-1._dp/3._dp*
     & s**(-2)*t**2*lY*Li2my-1._dp/6._dp*s**(-2)*t**2*lY*pisq+5._dp/18.
     & _dp*s**(-2)*t**2*lY**2-1._dp/9._dp*s**(-2)*t**2*lY**3-1._dp/3._dp
     & *s**(-2)*t**2*lX*Li2mx-1._dp/6._dp*s**(-2)*t**2*lX*pisq-5._dp/9.
     & _dp*s**(-2)*t**2*lX*lY+5._dp/18._dp*s**(-2)*t**2*lX**2-1._dp/9._dp
     & *s**(-2)*t**2*lX**3+1._dp/3._dp*s**(-1)*t**(-1)*u**2*lY*Lns-1._dp
     & /2._dp*s**(-1)*t**(-1)*u**2*lY**2*Lns-s**(-1)*u*lX*lY*Lns+1._dp
     & /2._dp*s**(-1)*u*lX**2*Lns+4._dp/3._dp*s**(-1)*t*Li3mx+4._dp/3._dp
     & *s**(-1)*t*Li3my+23._dp/18._dp*s**(-1)*t*pisq-1._dp/2._dp*
     & s**(-1)*t*pisq*Lns-8._dp/9._dp*s**(-1)*t*lY-4._dp/3._dp*s**(-1)*
     & t*lY*Li2my-2._dp/3._dp*s**(-1)*t*lY*pisq+13._dp/9._dp*s**(-1)*t*
     & lY**2-4._dp/9._dp*s**(-1)*t*lY**3+8._dp/9._dp*s**(-1)*t*lX-1._dp
     & /3._dp*s**(-1)*t*lX*Lns-4._dp/3._dp*s**(-1)*t*lX*Li2mx-2._dp/3._dp
     & *s**(-1)*t*lX*pisq-23._dp/9._dp*s**(-1)*t*lX*lY+10._dp/9._dp*
     & s**(-1)*t*lX**2-4._dp/9._dp*s**(-1)*t*lX**3+t**(-1)*u*pi*lY*im
     & *Lns+10._dp/3._dp*Lns-1._dp/2._dp*Lns**2+Li4mxoy-Li4my+2._dp
     & *Li3mx+17._dp/6._dp*pisq+1._dp/2._dp*pisq*Lns+1._dp/2._dp*pisq*
     & Li2mx-23._dp/360._dp*pisq**2-16._dp/9._dp*lY-lY*pisq+7._dp/3._dp
     & *lY**2+1._dp/12._dp*lY**2*pisq-2._dp/3._dp*lY**3+1._dp/24._dp*
     & lY**4-4._dp*lX+2._dp*lX*Lns+1._dp/2._dp*lX*Lns**2+lX*Li3my
     & -2._dp*lX*Li2mx-2._dp*lX*lY+1._dp/3._dp*lX*lY*pisq+1._dp/2._dp
     & *lX*lY**2-1._dp/6._dp*lX*lY**3+lX**2+1._dp/2._dp*lX**2*Lns+1.
     & _dp/2._dp*lX**2*Li2mx-1._dp/4._dp*lX**2*pisq-1._dp/2._dp*lX**2*lY
     & +1._dp/4._dp*lX**2*lY**2+1._dp/3._dp*lX**3*lY-1._dp/8._dp*lX**4
     & -82._dp/9._dp*pi*im+11._dp/3._dp*pi*im*Lns+1._dp/2._dp*pi*im*
     & Lns**2+pi*im*Li3mx+pi*im*Li3my-pi*im*Li2mx+pi*im*Li2my
     & +2._dp*pi*lY*im-1._dp/2._dp*pi*lY**2*im-2._dp*pi*lX*im+pi*lX
     & *im*Lns+1._dp/2._dp*pi*lX**2*im+1._dp/2._dp*pi*lX**2*lY*im-1._dp
     & /6._dp*pi*lX**3*im-pi*zeta3*im+2._dp*t*u**(-1)*Li4mxoy-2._dp*
     & t*u**(-1)*Li4my+11._dp/3._dp*t*u**(-1)*Li3mx+17._dp/6._dp*t*
     & u**(-1)*pisq+t*u**(-1)*pisq*Lns+t*u**(-1)*pisq*Li2mx-23._dp/
     & 180._dp*t*u**(-1)*pisq**2+1._dp/6._dp*t*u**(-1)*lY**2*pisq+1._dp/
     & 12._dp*t*u**(-1)*lY**4-118._dp/9._dp*t*u**(-1)*lX+17._dp/3._dp*t*
     & u**(-1)*lX*Lns+1._dp/2._dp*t*u**(-1)*lX*Lns**2+2._dp*t*u**(-1)*
     & lX*Li3my-11._dp/3._dp*t*u**(-1)*lX*Li2mx+11._dp/6._dp*t*u**(-1)*
     & lX*pisq+2._dp/3._dp*t*u**(-1)*lX*lY*pisq-1._dp/3._dp*t*u**(-1)*
     & lX*lY**3-8._dp/9._dp*t*u**(-1)*lX**2+2._dp*t*u**(-1)*lX**2*Lns
     & +t*u**(-1)*lX**2*Li2mx-1._dp/2._dp*t*u**(-1)*lX**2*pisq-11._dp
     & /6._dp*t*u**(-1)*lX**2*lY+1._dp/2._dp*t*u**(-1)*lX**2*lY**2+11._dp
     & /9._dp*t*u**(-1)*lX**3+2._dp/3._dp*t*u**(-1)*lX**3*lY-1._dp/4._dp
     & *t*u**(-1)*lX**4-11._dp/3._dp*t*u**(-1)*zeta3-118._dp/9._dp*t*
     & u**(-1)*pi*im+17._dp/3._dp*t*u**(-1)*pi*im*Lns+1._dp/2._dp*t*
     & u**(-1)*pi*im*Lns**2+2._dp*t*u**(-1)*pi*im*Li3mx+2._dp*t*
     & u**(-1)*pi*im*Li3my-11._dp/3._dp*t*u**(-1)*pi*im*Li2mx+11._dp/
     & 18._dp*t*u**(-1)*pi*pisq*im-67._dp/9._dp*t*u**(-1)*pi*lX*im-11._dp
     & /3._dp*t*u**(-1)*pi*lX*lY*im+11._dp/6._dp*t*u**(-1)*pi*lX**2*im
     & +t*u**(-1)*pi*lX**2*lY*im-1._dp/3._dp*t*u**(-1)*pi*lX**3*im-
     & 2._dp*t*u**(-1)*pi*zeta3*im+2._dp*t**2*u**(-2)*Li4mxoy-2._dp*
     & t**2*u**(-2)*Li4my+11._dp/3._dp*t**2*u**(-2)*Li3mx-t**2*
     & u**(-2)*Li3mx*Lns+t**2*u**(-2)*pisq*Li2mx-23._dp/180._dp*t**2*
     & u**(-2)*pisq**2+1._dp/6._dp*t**2*u**(-2)*lY**2*pisq+1._dp/12._dp
     & *t**2*u**(-2)*lY**4+2._dp*t**2*u**(-2)*lX*Li3my-11._dp/3._dp*
     & t**2*u**(-2)*lX*Li2mx+t**2*u**(-2)*lX*Li2mx*Lns+11._dp/6._dp*
     & t**2*u**(-2)*lX*pisq-1._dp/2._dp*t**2*u**(-2)*lX*pisq*Lns+2._dp/
     & 3._dp*t**2*u**(-2)*lX*lY*pisq-1._dp/3._dp*t**2*u**(-2)*lX*lY**3
     & -67._dp/18._dp*t**2*u**(-2)*lX**2+2._dp*t**2*u**(-2)*lX**2*Lns
     & -1._dp/4._dp*t**2*u**(-2)*lX**2*Lns**2+t**2*u**(-2)*lX**2*
     & Li2mx-1._dp/2._dp*t**2*u**(-2)*lX**2*pisq-11._dp/6._dp*t**2*
     & u**(-2)*lX**2*lY+1._dp/2._dp*t**2*u**(-2)*lX**2*lY*Lns+1._dp/2._dp
     & *t**2*u**(-2)*lX**2*lY**2+11._dp/9._dp*t**2*u**(-2)*lX**3-1._dp
     & /3._dp*t**2*u**(-2)*lX**3*Lns+2._dp/3._dp*t**2*u**(-2)*lX**3*lY
     & -1._dp/4._dp*t**2*u**(-2)*lX**4-11._dp/3._dp*t**2*u**(-2)*zeta3
     & +t**2*u**(-2)*zeta3*Lns+2._dp*t**2*u**(-2)*pi*im*Li3mx+2._dp
     & *t**2*u**(-2)*pi*im*Li3my-11._dp/3._dp*t**2*u**(-2)*pi*im*Li2mx
     & +t**2*u**(-2)*pi*im*Li2mx*Lns+11._dp/18._dp*t**2*u**(-2)*pi*
     & pisq*im-1._dp/6._dp*t**2*u**(-2)*pi*pisq*im*Lns-67._dp/9._dp*
     & t**2*u**(-2)*pi*lX*im-1._dp/2._dp*t**2*u**(-2)*pi*lX*im*Lns**2
     & -11._dp/3._dp*t**2*u**(-2)*pi*lX*lY*im+t**2*u**(-2)*pi*lX*lY*
     & im*Lns+11._dp/6._dp*t**2*u**(-2)*pi*lX**2*im-1._dp/2._dp*t**2*
     & u**(-2)*pi*lX**2*im*Lns+t**2*u**(-2)*pi*lX**2*lY*im-1._dp/3._dp
     & *t**2*u**(-2)*pi*lX**3*im-2._dp*t**2*u**(-2)*pi*zeta3*im+1._dp/
     & 6._dp*s*t**(-2)*u*lY**2*Lns+1._dp/3._dp*s*t**(-2)*u*pi*lY*im*Lns
     & -4._dp/3._dp*s*t**(-1)*Li3my+1._dp/6._dp*s*t**(-1)*pisq-8._dp/9.
     & _dp*s*t**(-1)*lY+4._dp/3._dp*s*t**(-1)*lY*Li2my-2._dp/3._dp*s*
     & t**(-1)*lY*pisq+13._dp/9._dp*s*t**(-1)*lY**2-4._dp/9._dp*s*
     & t**(-1)*lY**3+2._dp/3._dp*s*t**(-1)*lX*lY**2+4._dp/3._dp*s*
     & t**(-1)*zeta3-8._dp/9._dp*s*t**(-1)*pi*im+1._dp/3._dp*s*t**(-1)*
     & pi*im*Lns+4._dp/3._dp*s*t**(-1)*pi*im*Li2my-2._dp/9._dp*s*
     & t**(-1)*pi*pisq*im+23._dp/9._dp*s*t**(-1)*pi*lY*im-2._dp/3._dp*s
     & *t**(-1)*pi*lY**2*im+4._dp/3._dp*s*t**(-1)*pi*lX*lY*im+s*
     & u**(-1)*pisq*Lns+1._dp/2._dp*s*u**(-1)*lX*Lns**2+1._dp/2._dp*s*
     & u**(-1)*pi*im*Lns**2+1._dp/6._dp*s*t*u**(-2)*lX**2*Lns-11._dp/3.
     & _dp*s*t*u**(-2)*pi*lX*im*Lns-1._dp/3._dp*s**2*t**(-2)*Li3my+1._dp
     & /3._dp*s**2*t**(-2)*lY*Li2my-1._dp/6._dp*s**2*t**(-2)*lY*pisq+5.
     & _dp/18._dp*s**2*t**(-2)*lY**2-1._dp/9._dp*s**2*t**(-2)*lY**3+1._dp
     & /6._dp*s**2*t**(-2)*lX*lY**2+1._dp/3._dp*s**2*t**(-2)*zeta3+1._dp
     & /3._dp*s**2*t**(-2)*pi*im*Li2my-1._dp/18._dp*s**2*t**(-2)*pi*pisq
     & *im+5._dp/9._dp*s**2*t**(-2)*pi*lY*im-1._dp/6._dp*s**2*t**(-2)*
     & pi*lY**2*im+1._dp/3._dp*s**2*t**(-2)*pi*lX*lY*im-s**2*u**(-2)*
     & Li3mx*Lns+s**2*u**(-2)*lX*Li2mx*Lns-1._dp/2._dp*s**2*u**(-2)*
     & lX*pisq*Lns-1._dp/4._dp*s**2*u**(-2)*lX**2*Lns**2+1._dp/2._dp*
     & s**2*u**(-2)*lX**2*lY*Lns-1._dp/3._dp*s**2*u**(-2)*lX**3*Lns+
     & s**2*u**(-2)*zeta3*Lns+s**2*u**(-2)*pi*im*Li2mx*Lns-1._dp/6._dp
     & *s**2*u**(-2)*pi*pisq*im*Lns-1._dp/2._dp*s**2*u**(-2)*pi*lX*im*
     & Lns**2+s**2*u**(-2)*pi*lX*lY*im*Lns-1._dp/2._dp*s**2*u**(-2)*
     & pi*lX**2*im*Lns

c---end statement functions for amplitudes

c--- special functions
      lX=log(-t/s)
      lY=log(-u/s)
      Li2mx=ddilog(-t/s)
      Li2my=ddilog(-u/s)
      Li3mx=Li3(-t/s)
      Li3my=Li3(-u/s)
      Li4mx=Li4(-t/s)
      Li4my=Li4(-u/s)
      Li4mxoy=Li4(-t/u)
      Li4myox=Li4(-u/t)
      Lns=log(s)

      Amp(2,2,2,2,0)=
     &  Apppp0(t,u,lX,lY,Li2mx,Li2my,Li3mx,Li3my,Li4mx,Li4my,Li4mxoy)
      Amp(2,2,2,2,1)=
     &  Apppp1(t,u,lX,lY,Li2mx,Li2my,Li3mx,Li3my,Li4mx,Li4my,Li4mxoy)
      Amp(2,2,2,2,2)=
     &  Apppp2(t,u,lX,lY,Li2mx,Li2my,Li3mx,Li3my,Li4mx,Li4my,Li4mxoy)

c--- by parity
      Amp(1,1,1,1,0)=Amp(2,2,2,2,0)
      Amp(1,1,1,1,1)=Amp(2,2,2,2,1)
      Amp(1,1,1,1,2)=Amp(2,2,2,2,2)

      Amp(1,2,2,2,0)=
     &  Amppp0(t,u,lX,lY,Li2mx,Li2my,Li3mx,Li3my,Li4mx,Li4my,Li4mxoy)
      Amp(1,2,2,2,1)=
     &  Amppp1(t,u,lX,lY,Li2mx,Li2my,Li3mx,Li3my,Li4mx,Li4my,Li4mxoy)
      Amp(1,2,2,2,2)=
     &  Amppp2(t,u,lX,lY,Li2mx,Li2my,Li3mx,Li3my,Li4mx,Li4my,Li4mxoy)

c--- by 1<->2 exchange (amplitude unchanged)
      Amp(2,1,2,2,0)=Amp(1,2,2,2,0)
      Amp(2,1,2,2,1)=Amp(1,2,2,2,1)
      Amp(2,1,2,2,2)=Amp(1,2,2,2,2)

c--- by permutation, 1<->3, 2<->4 (amplitude unchanged)
      Amp(2,2,1,2,0)=Amp(1,2,2,2,0)
      Amp(2,2,1,2,1)=Amp(1,2,2,2,1)
      Amp(2,2,1,2,2)=Amp(1,2,2,2,2)

c--- by exchange, 3<->4 (amplitude unchanged)
      Amp(2,2,2,1,0)=Amp(1,2,2,2,0)
      Amp(2,2,2,1,1)=Amp(1,2,2,2,1)
      Amp(2,2,2,1,2)=Amp(1,2,2,2,2)

c--- by parity
      Amp(2,1,1,1,0)=Amp(1,2,2,2,0)
      Amp(2,1,1,1,1)=Amp(1,2,2,2,1)
      Amp(2,1,1,1,2)=Amp(1,2,2,2,2)

c--- by parity
      Amp(1,2,1,1,0)=Amp(2,1,2,2,0)
      Amp(1,2,1,1,1)=Amp(2,1,2,2,1)
      Amp(1,2,1,1,2)=Amp(2,1,2,2,2)

c--- by parity
      Amp(1,1,2,1,0)=Amp(2,2,1,2,0)
      Amp(1,1,2,1,1)=Amp(2,2,1,2,1)
      Amp(1,1,2,1,2)=Amp(2,2,1,2,2)

c--- by parity
      Amp(1,1,1,2,0)=Amp(2,2,2,1,0)
      Amp(1,1,1,2,1)=Amp(2,2,2,1,1)
      Amp(1,1,1,2,2)=Amp(2,2,2,1,2)

      Amp(1,1,2,2,0)=
     &  Ammpp0(t,u,lX,lY,Li2mx,Li2my,Li3mx,Li3my,Li4mx,Li4my,Li4mxoy)
      Amp(1,1,2,2,1)=
     &  Ammpp1(t,u,lX,lY,Li2mx,Li2my,Li3mx,Li3my,Li4mx,Li4my,Li4mxoy)
      Amp(1,1,2,2,2)=
     &  Ammpp2(t,u,lX,lY,Li2mx,Li2my,Li3mx,Li3my,Li4mx,Li4my,Li4mxoy)

c--- by parity
      Amp(2,2,1,1,0)=Amp(1,1,2,2,0)
      Amp(2,2,1,1,1)=Amp(1,1,2,2,1)
      Amp(2,2,1,1,2)=Amp(1,1,2,2,2)

      Amp(1,2,1,2,0)=
     &  Ampmp0(t,u,lX,lY,Li2mx,Li2my,Li3mx,Li3my,Li4mx,Li4my,Li4mxoy)
      Amp(1,2,1,2,1)=
     &  Ampmp1(t,u,lX,lY,Li2mx,Li2my,Li3mx,Li3my,Li4mx,Li4my,Li4mxoy)
      Amp(1,2,1,2,2)=
     &  Ampmp2(t,u,lX,lY,Li2mx,Li2my,Li3mx,Li3my,Li4mx,Li4my,Li4mxoy)

c--- by 1<->2 exchange
      Amp(2,1,1,2,0)=
     &  Ampmp0(u,t,lY,lX,Li2my,Li2mx,Li3my,Li3mx,Li4my,Li4mx,Li4myox)
      Amp(2,1,1,2,1)=
     &  Ampmp1(u,t,lY,lX,Li2my,Li2mx,Li3my,Li3mx,Li4my,Li4mx,Li4myox)
      Amp(2,1,1,2,2)=
     &  Ampmp2(u,t,lY,lX,Li2my,Li2mx,Li3my,Li3mx,Li4my,Li4mx,Li4myox)

c--- by parity
      Amp(2,1,2,1,0)=Amp(1,2,1,2,0)
      Amp(2,1,2,1,1)=Amp(1,2,1,2,1)
      Amp(2,1,2,1,2)=Amp(1,2,1,2,2)

c--- by parity
      Amp(1,2,2,1,0)=Amp(2,1,1,2,0)
      Amp(1,2,2,1,1)=Amp(2,1,1,2,1)
      Amp(1,2,2,1,2)=Amp(2,1,1,2,2)

c--- Multiply this contribution by poles, c.f. Eqs.(2.14) and (2.11)
c--- This is taken from hep-ph/0109078 Eq.(2.11); note however that the
c--- log proportional to the beta-function coefficient is added in
c---  Eq.(2.11) and subtracted again in Eq.(4.5), so we omit it here.
      xlog=lnrat(musq,-s)
      Igg_gaga2=-xn
      Igg_gaga1=-xn*(xlog+(11._dp/6._dp-real(nflav,dp)/(3._dp*xn)))
      Igg_gaga0=-xn*(0.5_dp*xlog**2)
      do h1=1,2
      do h2=1,2
      do h3=1,2
      do h4=1,2
c      write(6,*) 'Amp 0',Amp(h1,h2,h3,h4,0)
c      write(6,*) 'Amp 1',Amp(h1,h2,h3,h4,1)
c      write(6,*) 'Amp 2',Amp(h1,h2,h3,h4,2)
      M1(h1,h2,h3,h4)=
     & +Igg_gaga2*Amp(h1,h2,h3,h4,0)*epinv**2
     & +(Igg_gaga2*Amp(h1,h2,h3,h4,1)+Igg_gaga1*Amp(h1,h2,h3,h4,0))*epinv
     & +(Igg_gaga2*Amp(h1,h2,h3,h4,2)+Igg_gaga1*Amp(h1,h2,h3,h4,1)
     &  +Igg_gaga0*Amp(h1,h2,h3,h4,0))
      enddo
      enddo
      enddo
      enddo
c      pause


      return
      end

