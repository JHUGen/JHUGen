      subroutine qqb_QQb_v(p,msq)
      implicit none
************************************************************************
*     Author: R.K. Ellis                                               *
*     March, 2002.                                                     *
*     Virtual corrections                                              *
*     Calculates the element squared for the process                   *
*                                                                      *
*     q(-p1)+qbar(-p2) -> Q(p3)+Qbar(p4)                               *
*                                                                      * 
*     The mass of the heavy quark is passed in the common block        *
*      via the variable mass2                                          *
************************************************************************
      include 'constants.f'
      include 'nflav.f'
      include 'qcdcouple.f'
      include 'sprods_com.f'
      include 'scale.f'
      include 'breit.f'
      integer j,k
      double precision msq(-nf:nf,-nf:nf),p(mxpart,4),t1,ro
      double precision qqsym,qqasy,ggsym,xm2

      xm2=mass2**2
      ro=4d0*xm2/s(1,2)
      t1=-s(1,3)/s(1,2)                                                         

      call dotem(4,p,s)
      call virteval(t1,ro,qqsym,qqasy,ggsym)  
      
      msq(:,:)=0d0
      do j=-nf,nf
      k=-j
      if     (j .eq. 0) then
      msq(0,0)=avegg*gsq**3*ggsym/8d0/pisq
      elseif (j .gt. 0) then
      msq(j,k)=aveqq*(qqsym+qqasy)/8d0/pisq
      elseif (j .lt. 0) then
      msq(j,k)=aveqq*(qqsym-qqasy)/8d0/pisq
      endif
      enddo

      return
      end


      subroutine virteval(t1,ro,qqsym,qqasy,ggsym)  
************************************************************************
*     The mass of the heavy quark is passed in the common block        *
*      via the variable mass2                                          *
************************************************************************
      implicit none
      include 'constants.f'
      include 'nflav.f'
      include 'qcdcouple.f'
      include 'epinv.f'
      include 'scheme.f'
      include 'scale.f'
      include 'breit.f'
      external ddilog
      double precision t1,t2,ro,tbar,ubar,b,xlp,xlm,vlpm,vlsm,vltm,vlwm
      double precision vlbl,vdmp,vdmb,f1,f2,f3,f4,f5t1,f5t2,
     . qqQQv_0,qqQQv_1,ddilog
      double precision vdt,vdw,xlf,rmuom2,epin2,epin
      double precision qqss,qqsv,qqas,qqav
      double precision ggs
      double precision ggv1,ggv2,ggv3,ggv4,ggv5,ggv6,
     . ggv7,ggv8,ggv9,ggv10,ggv11
      double precision ggQQv,ggQQov,ggQQsv,qqsym,qqasy,ggsym,
     . ggQQv_0,ggQQov_0,ggQQsv_0,ggQQv_1,ggQQov_1,ggQQsv_1,
     . ggQQv_2
      
C---we arrive at the dred scheme by taking a HV result and applying
C---a finite renormalization on the 
C---namely a factor of as/4/pi*cf/2 for every massles quark leg (amplitude)
C---namely a factor of as/4/pi*xn/6 for every massles gluon leg (amplitude)
C---namely a factor of as/2/pi*cf/2 for every massles quark leg (square)
C---namely a factor of as/2/pi*xn/6 for every massles gluon leg (square)
      scheme='dred'
     
C Notation for logarithms and dilogarithms.
C Integrals were in unphysical region, but contination has now 
C been performed; 
C In unphysical region,        In physical region after continuation.
C vltm=log(at/m^2),
C vlpm=log(-lp/lm),            log((1+b)/(1-b))
C vlsm=log(as/m^2),            log(s/m^2)
C vlwm=log(aw/m^2),
C vlbl=log(-b/lm),             log(2*b/(1-b))

C vdw=li[2]((aw-m^2)/aw)-1/2*vlwm^2,
C vdt=li[2]((at-m^2)/at)-1/2*vltm^2,
C vdmp=li[2](-lm/lp),
C vdmb=li[2](-lm/b)+1/2*vlbl^2,

C as=-s
C lp=(1+b)/2,
C lm=(1-b)/2,
C at=t1*s=m^2-t,
C aw=t2*s=m^2-u,
C b=sqrt(1-ro),
C m=sqrt(ro*s)/2
C TBAR=-T/S=t1-1/4*ro
C UBAR=-U/S=t2-1/4*ro

C VIRgg is the complete answer for virt diagrams in units of
C as/2/pi*g^4*Gamma(1-EP)/Gamma(1-2*EP)*(4*pi*mu^2/XM2)^EP

      xlf=dfloat(nflav)
      rmuom2=2d0*dlog(scale/mass2)
      t2=1d0-t1
      tbar=t1-0.25d0*ro
      ubar=t2-0.25d0*ro
      b=dsqrt(1d0-ro)   
      xlp=0.5d0*(1d0+b)
      xlm=0.5d0*(1d0-b)
      vlpm=dlog(xlp/xlm)     
      vlsm=dlog(4d0/ro)      
      vltm=dlog(4d0*t1/ro)      
      vlwm=dlog(4d0*t2/ro)      
      vlbl=dlog(b/xlm)       
      vdw=ddilog(1d0-ro/(4d0*t2))-0.5d0*vlwm**2      
      vdt=ddilog(1d0-ro/(4d0*t1))-0.5d0*vltm**2
      vdmp=ddilog(-xlm/xlp)    
      vdmb=ddilog(-xlm/b)+0.5d0*vlbl**2  

c--- Q-Qbar and Qbar-Q contributions

      f1=(vlpm**2/2d0-2d0*vdmb-pisq/3d0)/b
      f2=(-b*vlsm+vlpm**2/4d0+vdmp+pisq/12d0)/b**3
      f3=(-b**3*vlsm-3*b*vlsm+0.75d0*vlpm**2
     . +3d0*vdmp+pisq/4d0+2d0*b**3)/b**5
      f4 = (vlpm**2/4d0+vdmp+pisq/12d0)/b
      f5t1 = (vltm**2+vdt+pisq/6d0)/t1**3
      f5t2 = (vlwm**2+vdw+pisq/6d0)/t2**3
      
c---  Singular parts
c---  qqQQv is the lowest order matrix element squared without
c---  averaging over initial colors or spins:
c---  use _0 for O(1) piece, _1 for O(eps)
      qqQQv_0=gsq**2*V*(2d0*t1**2+2d0*t2**2+ro)
      if (scheme .eq. 'tH-V') then
      qqQQv_1=gsq**2*V*(-2d0)
      else
      qqQQv_1=0d0
      endif      
      
c--- These are the singular pieces, written in such a way that
c--- the limit EPINV -> 0 is smooth

c--- overall factor in the virtual terms is (fourpi*mass^2)^(epsilon)      
c--- to match with our usual definition, we should thus multiply by
c--- (musq/mass^2)^(epsilon)
c--- Note that the overall factor in this definition of the virtual
c--- functions also differs from ours in that we have 1/Gamma(1-ep).
c--- This gives an extra (-pisqo6) compared to Gamma(1-ep)/Gamma(1-2*ep)
      epin=epinv+rmuom2
      epin2=epinv**2+epinv*rmuom2+half*rmuom2**2-pisqo6

      qqss=
     .  gsq/V*qqQQv_0*( 
     . -2d0*V/2/xn*EPIN2-V/xn*EPIN
     . -3d0*V/2/xn*EPIN+xn*(vltm+vlwm)*EPIN
     . -(0.5d0*vlpm/b*(1d0+b**2)+vlsm)/xn*EPIN)
     . +gsq/V*qqQQv_1*(
     . -2d0*V/2/xn*EPIN-V/xn
     . -3d0*V/2/xn+xn*(vltm+vlwm)
     . -(0.5d0*vlpm/b*(1d0+b**2)+vlsm)/xn)
     
      qqas=gsq/V*qqQQv_0*(XN4/xn*(vltm-vlwm))*EPIN
     .    +gsq/V*qqQQv_1*(XN4/xn*(vltm-vlwm))

C---  Finite symmetric part
       qqsv=
     . +gsq/V*qqQQv_0*(xn*(22d0/9d0-0.5d0*vlsm*(vltm+vlwm)
     . +0.25d0*vlsm**2+11d0/3d0*rmuom2+11d0/6d0*vlsm-pisq/6d0)
     . +1d0/xn*(+6d0-pisq/3d0-1.5d0*vlsm+0.5d0*vlsm**2-b*vlpm
     . -0.5d0*(1d0+b**2)*(pisq/b+f1))
     . +TR*XLF*(-20d0/9d0+4d0/3d0*vlsm-4d0/3d0*rmuom2)
     . +TR*(-20d0/9d0+4d0/3d0*vlpm*b*(1d0+ro/2d0)-4d0/3d0*ro))
     . +gsq**3*0.5d0*xn*(t1-t2)**2*f3
     . +gsq**3*f2*(xn*(-5d0*(t1**2+t2**2)+2d0
     . +b**2*(0.5d0+6d0*t1**2+6d0*t2**2)-b**4))
     . -gsq**3*0.5d0/xn*vlpm/b*((t1-t2)**2+b**2)
     . +gsq**3*xn*((t2-t1)*(vdt-vdw)+0.5d0*ro*(vdt+vdw)
     . +0.5d0*vlsm**2*(t1**2 +t2**2)+1d0-pisq/12d0*ro
     . -vlsm*(vltm+vlwm)*(t1**2+t2**2)-vlsm*(vltm-vlwm)*(t1-t2)
     . -(vltm/TBAR+vlwm/UBAR)*(0.5d0*ro-t1*t2)
     . -vlsm*(1.5d0+2d0*ro))

c--- extra finite terms in DR scheme
      if (scheme .eq. 'dred') then
        qqsv=qqsv+0.5d0*(xn-1d0/xn)*gsq*qqQQv_0/V
      endif
        
C--- Finite antisymmetric part
      qqav=gsq/V*qqQQv_0*(-XN4/2d0/xn*vlsm*(vltm-vlwm))
     . +gsq**3*XN4/xn*(-(t1-t2)*(vdt+vdw)+0.5d0*ro*(vdt-vdw)
     . +(pisq/6d0+vlsm*(3d0-ro)+0.5d0*vlsm**2-vlsm*(vltm+vlwm))*(t1-t2)
     . -vlsm*(vltm-vlwm)*(t1**2+t2**2)
     . -(vltm/tbar-vlwm/ubar)*(0.5d0*ro-t1*t2)
     . +f2*(t1-t2)*(-1d0+2d0*b**2+b**4))

      qqsym=V*(qqss+qqsv)
      qqasy=V*(qqas+qqav)

c--- g-g contribution

c---  Singular parts
c---  ggQQv is the lowest order matrix element squared without
c---  averaging over initial colors or spins:
c---  use _0 for O(1) piece, _1 for O(eps) and _2 for O(eps**2)
      ggQQv_0=2d0/xn*(V/t1/t2-2d0*xnsq)
     . *((t1**2+t2**2)+ro-ro**2/4d0/t1/t2)
      ggQQov_0=2d0/xn*((xnsq-2d0)/t1/t2-2d0*xnsq)
     . *((t1**2+t2**2)+ro-ro**2/4d0/t1/t2)
      ggQQsv_0=2d0/xn*(1d0/t1/t2)
     . *((t1**2+t2**2)+ro-ro**2/4d0/t1/t2)
      if (scheme .eq. 'tH-V') then
      ggQQv_1=2d0/xn*(V/t1/t2-2d0*xnsq)
     . *(-(t1**2+t2**2)-1d0)
      ggQQv_2=2d0/xn*(V/t1/t2-2d0*xnsq)
      ggQQov_1=2d0/xn*((xnsq-2d0)/t1/t2-2d0*xnsq)
     . *(-(t1**2+t2**2)-1d0)
      ggQQsv_1=2d0/xn*(1d0/t1/t2)
     . *(-(t1**2+t2**2)-1d0)
      else
      ggQQv_1=0d0
      ggQQv_2=0d0
      ggQQov_1=0d0
      ggQQsv_1=0d0
      endif      

c--- These are the singular pieces, written in such a way that
c--- the limit EPINV -> 0 is smooth
      ggs=
     .  (-2d0*xn)*(EPIN2*ggQQv_0+EPIN*ggQQv_1+ggQQv_2)
     . +(2d0*(2d0*TR/3d0*XLF-11d0/6d0*xn)-V/xn)*(EPIN*ggQQv_0+ggQQv_1)
     . +(EPIN*ggQQsv_0+ggQQsv_1)*(1d0+b**2)/b*vlpm*V/2d0/xn
     . +(EPIN*ggQQov_0+ggQQov_1)*(vlsm*xn-(1d0+b**2)/b*vlpm/xn/2d0)
     . +2d0*(EPIN*ggQQsv_0+ggQQsv_1)*(vlsm-vltm-vlwm)*xn
     . +4d0*xnsq*EPIN*vltm
     . *(-ro**2/4d0/t1**2+t2/t1+ro*t2/t1-2d0*t2**2)
     . +4d0*xnsq*EPIN*vlwm
     . *(-2d0*t1**2-ro**2/4d0/t2**2+t1/t2+ro*t1/t2)
      if (scheme .eq. 'tH-V') then
c--- These are the extra finite pieces (akin to _1) for the last 4 lines
c--- and they appear to have already been added in the finite pieces
c--- below, so we subtract it later there
      ggs=ggs+(
     . +4d0*xnsq*vltm
     . *(-2d0*t2/t1+2d0*t2**2)
     . +4d0*xnsq*vlwm
     . *(-2d0*t1/t2+2d0*t1**2))
      endif
      
c--- extra finite terms in DR scheme
      if (scheme .eq. 'dred') then
        ggs=ggs+xn/3d0*ggQQv_0
      endif
        
C--- overall factor of V removed
      ggQQv=2d0/xn*(V/t1/t2-2d0*xnsq) 
     & *(t1**2+t2**2+ro-ro**2/(4d0*t1*t2))
      ggQQov=2d0/xn*((xnsq-2d0)/(t1*t2)-2d0*xnsq)
     & *(t1**2+t2**2+ro-ro**2/(4d0*t1*t2))
      ggQQsv=2d0/(xn*t1*t2)
     & *(t1**2+t2**2+ro-ro**2/(4d0*t1*t2))

      ggv1 = f3*(-2d0*t2**2+1d0/(t1*t2)/4d0-2d0*t1**2)*xnsq
     1   +f4*(24d0*t1*t2-ro**2/(t1*t2)/2d0+ro/(t1*t2)+11d0/(4d0*t1*t2)
     2  -4d0*ro-17d0)*xnsq+f2*(-20d0*t1*t2+(-3d0)/(2d0*t1*t2)+11d0)
     3  *xnsq-2*ggqqv*rmuom2*((-11d0)*xn/6d0+2d0*tr*xlf/3d0)+(b**2+1d0)
     4  *ggqqsv*pisq*v/(b*xn)/2d0-(b**2+1d0)*ggqqov*pisq/(b*xn)/2d0
     5  -f5t2*ro**2*v/xnsq/4d0-f5t1*ro**2*v/xnsq/4d0+f4*((1d0-ro/2d0)
     6  *ro**2*(t2**2+t1**2)/(t1**2*t2**2)+(-ro**3+5d0*ro**2-ro-6d0)
     7  /(t1*t2)+4d0)/xnsq+f1*(-(1d0-ro/2d0)*ro**2*(t2**2+t1**2)/
     8  (t1**2*t2**2)/2d0+(ro**3/2d0-2d0*ro**2+ro+2d0)/(t1*t2)+2d0*ro
     9  -4d0)/xnsq+f1*(4d0*t2**2+4d0*ro*t1*t2-(1d0-ro/2d0)*ro**2/(t1*t2)
     :  +4d0*t1**2-2d0*ro**2+2d0*ro)+f4*(-ro**3/(t1*t2)+4d0*ro**2
     ;  /(t1*t2)-4d0/(t1*t2)-2d0*ro**2-2d0*ro+6d0)
      ggv2 = (-ro**3/t2/2d0+2d0*ro**2/t2-ro/t2-2d0/t2-ro**3/t2**2/2d0
     1   +ro**2/t2**2-ro**3/t1/2d0+3d0*ro**2/t1-4d0/t1+2d0)*vlpm*vlwm/
     2   (b*xnsq)+(-ro**3/t2/2d0+3d0*ro**2/t2-4d0/t2-ro**3/t1/2d0
     3   +2d0*ro**2/t1-ro/t1-2/t1-ro**3/t1**2/2d0+ro**2/t1**2+2d0)
     4   *vlpm*vltm/(b*xnsq)+(ro**3/t2/2d0+(-5d0)*ro**2/(2d0*t2)
     5   +ro/t2/2d0+3d0/t2+ro**3/t2**2/4d0-ro**2/t2**2/2d0+ro**3/t1/2d0
     6   +(-5d0)*ro**2/(2d0*t1)+ro/t1/2d0+3/t1+ro**3/t1**2/4d0
     7   -ro**2/t1**2/2d0-2d0)*vlpm*vlsm/(b*xnsq)+(ro/(t1*t2)-4d0)*
     8   vlpm/(b*xnsq)
      ggv3 = (ro**3/t2-ro**2/t2+ro**3/t1-ro**2/t1-4d0*ro**3+4d0*ro**2)
     1   *tr*vlpm*xn/b+(-ro**3/t2/2d0+ro**2/t2-2d0*t1-ro**3/t1/2d0
     2   +3d0*ro**2/t1-4d0/t1-ro**2-ro+4d0)*vlpm*vlwm/b+(-ro**3/t2/2d0
     3   +3d0*ro**2/t2-4d0/t2+2d0*t1-ro**3/t1/2d0+ro**2/t1-ro**2-ro+2d0)
     4   *vlpm*vltm/b+(ro**3/t2/2d0-2*ro**2/t2+2/t2+ro**3/t1/2d0
     5   -2d0*ro**2/t1+2d0/t1+ro**2+ro-3d0)*vlpm*vlsm/b+(ro**2/t2/2d0
     6   -ro/t2/4d0-8d0*ro*t1**2+12d0*t1**2+8d0*ro*t1-12d0*t1
     7   +ro**2/t1/2d0-ro/t1/4d0-2d0*ro**2+ro+1d0)*vlpm/b
      ggv4 = (5d0*t1**2-2d0*ro*t1-8d0*t1+7d0*ro**2/8d0+3d0*ro-1d0)
     1   *vltm*xnsq/tbar+(-ro*t1**2/4d0+(-3d0)*ro**2*t1/8d0+ro*t1/4d0
     2   +3d0*ro**3/32d0-ro**2/8d0)*vltm*xnsq/tbar**2+(8d0*t1**2
     3   -16d0*t1+2d0*ro/t1-4d0/t1-ro**2/t1**2+16d0)*vltm*xnsq
     4   +(-14d0*ro/(3d0*t2)-7d0/(2d0*t2)+ro**2/t2**2+16d0*t1**2
     5   -16d0*t1+(-14d0)*ro/(3d0*t1)+(-7d0)/(2d0*t1)+ro**2/t1**2
     6   +26d0*ro/3d0+15d0)*xnsq-2d0*ro**2/t2+17d0*ro/(2d0*t2)
     7   +8d0/t2-2d0*ro**2/t2**2+ro/t2**2-16d0*t1**2+16d0*t1
     8   -2d0*ro**2/t1+17d0*ro/(2d0*t1)+8d0/t1-2d0*ro**2/t1**2
     9   +ro/t1**2-10d0*ro-19d0
      ggv5 = (-2d0*ro/t2-2d0/t2+ro**2/t2**2+8d0*t1**2-8d0*t1+2d0*ro+6d0)
     1  *vlsm*vlwm*xn**2+(2d0*ro/t2-4d0/t2-ro**2/t2**2+8d0*t1**2+8d0)
     2  *vlwm*xnsq+(8d0*t1**2-8d0*t1-2d0*ro/t1-2d0/t1+ro**2/t1**2+2d0*ro
     3  +6d0)*vlsm*vltm*xnsq+(-ro/t2/2d0-1d0/t2/2d0-ro/t1/2d0-1/t1/2d0
     4  +ro+1d0)*vlsm**2*xnsq+(ro/t2/2d0-1d0/t2/2d0-8d0*t1**2+8d0*t1
     5  +ro/t1/2d0-1d0/t1/2d0-2d0*ro+2d0)*vlsm*xnsq
      ggv6 = (5d0*t1**2+2d0*ro*t1-2d0*t1+7d0*ro**2/8d0+ro-4d0)*vlwm
     1  *xnsq/ubar+(-ro*t1**2/4d0+3d0*ro**2*t1/8d0+ro*t1/4d0
     2  +3d0*ro**3/32d0-ro**2/2d0)*vlwm*xnsq/ubar**2+(2d0*ro/t2+2d0/t2
     3  -8d0*t1-2d0*ro+2d0)*vdw*xnsq+(8d0*t1+2d0*ro/t1+2d0/t1-2d0*ro
     4  -6d0)*vdt*xnsq+(ro/2d0-ro*t1/4d0)*xnsq/ubar+(ro*t1/4d0+ro/4d0)
     5  *xnsq/tbar+pisq*(ro/t2/6d0+1d0/t2/6d0-ro**2/t2**2/12d0
     6  +(-4d0)*t1**2/3d0+4d0*t1/3d0+ro/t1/6d0+1d0/t1/6d0
     7  -ro**2/t1**2/12d0-ro/3d0-1d0)*xnsq+(-ro**2/t2+ro/t2-2d0/t2
     8  -ro**2/t1+3d0*ro/t1+4d0/t1-ro**2/t1**2+ro/t1**2)*vltm/xnsq
     9  +(2d0*ro**2/t2-5d0*ro/t2-5d0/t2+ro**2/t2**2-ro/t2**2
     :   +2d0*ro**2/t1-5d0*ro/t1-5d0/t1+ro**2/t1**2-ro/t1**2+6d0)/xnsq
      ggv7 = (t1+4d0/t1+ro-4d0)*vlwm/(ubar*xnsq)+(3d0*ro*t1/4d0-2d0/t1
     1  +ro**2/4d0-ro/4d0+2d0)*vlwm/(ubar**2*xnsq)+(-ro**2/t2
     2  +3d0*ro/t2+4/t2-ro**2/t2**2+ro/t2**2-ro**2/t1+ro/t1-2d0/t1)
     3  *vlwm/xnsq+(-3d0*ro**2/(4d0*t2)+ro/t2+4d0/t2
     4  -3d0*ro**2/(4d0*t1)-ro/t1+2d0/t1+ro**2/t1**2/4d0-2d0)*vltm**2
     5  /xnsq+(4d0/t2-t1+ro-3d0)*vltm/(tbar*xnsq)+(-2d0/t2-3d0*ro*t1/4d0
     6  +ro**2/4d0+ro/2d0+2d0)*vltm/(tbar**2*xnsq)+(-ro**2/t2/4d0
     7  +3d0/(2d0*t2)-ro**2/t1/4d0+3d0/(2d0*t1)-1d0)*vlpm**2/xn**2
      ggv8 = (-3d0*ro**2/(4d0*t2)-ro/t2+2d0/t2+ro**2/t2**2/4d0
     1  -3d0*ro**2/(4d0*t1)+ro/t1+4d0/t1-2d0)*vlwm**2/xnsq
     2  +(-3d0*ro**2/(4d0*t2)-ro/t2+2d0/t2+ro**2/t2**2/4d0
     3  +(-3d0)*ro**2/(4d0*t1)+ro/t1+4d0/t1-2d0)
     3  *vdw/xnsq+((-3d0)*ro**2/(4d0*t2)+ro/t2+4d0/t2
     4  -3d0*ro**2/(4d0*t1)-ro/t1+2d0/t1+ro**2/t1**2/4d0-2d0)*vdt/xnsq
     5  +(2d0/t1-ro/4d0-2d0)/(ubar*xnsq)+(2d0/t2-ro/4d0-2d0)/(tbar*xnsq)
     6  +pisq*(-1d0/t2/2d0+ro**2/t2**2/24d0-1d0/t1/2d0
     7  +ro**2/t1**2/24d0+1d0/3d0)/xnsq
      ggv9 = (-2d0*ro**2/t2+4d0*ro/t2+4d0/t2-2d0*ro**2/t2**2
     1   -2d0*ro**2/t1+4d0*ro/t1+4d0/t1-2d0*ro**2/t1**2-8d0)*vltm*vlwm
     2   +(-ro**2/t2+2d0*ro/t2+2d0/t2-2d0*t1+ro/t1-2d0/t1-2d0*ro+2d0)
     3   *vltm**2+(-12d0/t2-t1**2+3d0*ro*t1+9d0*t1-7d0*ro**2/8d0
     4   -4d0*ro+12d0)*vltm/tbar+(2d0/t2-3d0*ro*t1**2/4d0+5d0*ro**2*t1
     4   /8d0+5d0*ro*t1/2d0-3d0*ro**3/32d0-5d0*ro**2/8d0-ro/2d0-2d0)
     5   *vltm/tbar**2+(ro**2/t2+10d0/t2+ro**2/t1-4d0*ro/t1-8d0/t1
     6   +2d0*ro**2/t1**2-ro/t1**2)*vltm+(-ro**2/t2/8d0+ro/t2/4d0
     7   +1d0/t2-ro**2/t1/8d0+ro/t1/4d0+1d0/t1-ro+(-3d0)/2d0)*vlpm**2
      ggv10 = (ro/t2-2d0/t2+2d0*t1-ro**2/t1+2d0*ro/t1+2d0/t1-2d0*ro)
     1  *vlwm**2+(-t1**2-3d0*ro*t1-7*t1-12/t1-7d0*ro**2/8d0-ro+20d0)
     2  *vlwm/ubar+(-3d0*ro*t1**2/4d0-5d0*ro**2*t1/8d0-ro*t1+2d0/t1
     3  -3d0*ro**3/32d0+5d0*ro/4d0-2d0)*vlwm/ubar**2+(ro**2/t2-4d0*ro/t2
     4  -8d0/t2+2d0*ro**2/t2**2-ro/t2**2+ro**2/t1+10d0/t1)*vlwm
     5  +(ro**2/t2-ro/t2-4d0/t2+2d0*t1-2d0*ro+4d0)*vdw
     6  +(-2d0*t1+ro**2/t1-ro/t1-4d0/t1-2d0*ro+6d0)*vdt
      ggv11 = (ro/t2/3d0+ro/t1/3d0-4d0*ro/3d0)*tr*xlf*xn
     1 +(ro**2/t2/4d0+ro**2/t1/4d0-ro**2)*tr*vlpm**2*xn
     2 +(2d0*ro**2/t2+ro/t2/3d0+2d0*ro**2/t1+ro/t1/3d0-8d0*ro**2
     3 -4d0*ro/3d0)*tr*xn+pisq*(-ro**2/t2/4d0-ro**2/t1/4d0+ro**2)*tr*xn
     4 +(-3d0*ro*t1/4d0-2d0/t1-ro**2/4d0+3d0*ro/4d0+2d0)/ubar
     5 +(-2d0/t2+3d0*ro*t1/4d0-ro**2/4d0+2d0)/tbar
     6 +pisq*(-ro**2/t2/24d0+ro/t2/4d0-1d0/t2+ro**2/t2**2/3d0
     7 -ro**2/t1/24d0+ro/t1/4d0-1d0/t1+ro**2/t1**2/3d0+ro/3d0+11d0/6d0)


c--- This is subtracting the extra finite piece mentioned above
      ggv11=ggv11-(
     . +4d0*xnsq*vltm
     . *(-2d0*t2/t1+2d0*t2**2)
     . +4d0*xnsq*vlwm
     . *(-2d0*t1/t2+2d0*t1**2))

C---replace the overall factor of V which was removed
      ggsym=V*(ggs+ggv1+ggv2+ggv3+ggv4+ggv5
     &     +ggv6+ggv7+ggv8+ggv9+ggv10+ggv11)
      return        
      end 

