      subroutine qqb_QQb_v(p,msq)
      implicit none
      include 'types.f'

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
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'nflav.f'
      include 'qcdcouple.f'
      include 'sprods_com.f'
      include 'scale.f'
      include 'breit.f'
      integer:: j,k
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4),t1,ro
      real(dp):: qqsym,qqasy,ggsym,xm2

      xm2=mass2**2
      ro=4._dp*xm2/s(1,2)
      t1=-s(1,3)/s(1,2)

      call dotem(4,p,s)
      call virteval(t1,ro,qqsym,qqasy,ggsym)

      msq(:,:)=0._dp
      do j=-nf,nf
      k=-j
      if     (j == 0) then
      msq(0,0)=avegg*gsq**3*ggsym/8._dp/pisq
      elseif (j > 0) then
      msq(j,k)=aveqq*(qqsym+qqasy)/8._dp/pisq
      elseif (j < 0) then
      msq(j,k)=aveqq*(qqsym-qqasy)/8._dp/pisq
      endif
      enddo

      return
      end


      subroutine virteval(t1,ro,qqsym,qqasy,ggsym)
      implicit none
      include 'types.f'
************************************************************************
*     The mass of the heavy quark is passed in the common block        *
*      via the variable mass2                                          *
************************************************************************

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'nflav.f'
      include 'qcdcouple.f'
      include 'epinv.f'
      include 'scheme.f'
      include 'scale.f'
      include 'breit.f'
      external ddilog
      real(dp):: t1,t2,ro,tbar,ubar,b,xlp,xlm,vlpm,vlsm,vltm,vlwm
      real(dp):: vlbl,vdmp,vdmb,f1,f2,f3,f4,f5t1,f5t2,
     & qqQQv_0,qqQQv_1,ddilog
      real(dp):: vdt,vdw,xlf,rmuom2,epin2,epin
      real(dp):: qqss,qqsv,qqas,qqav
      real(dp):: ggs
      real(dp):: ggv1,ggv2,ggv3,ggv4,ggv5,ggv6,
     & ggv7,ggv8,ggv9,ggv10,ggv11
      real(dp):: ggQQv,ggQQov,ggQQsv,qqsym,qqasy,ggsym,
     & ggQQv_0,ggQQov_0,ggQQsv_0,ggQQv_1,ggQQov_1,ggQQsv_1,
     & ggQQv_2

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

      xlf=real(nflav,dp)
      rmuom2=2._dp*log(scale/mass2)
      t2=1._dp-t1
      tbar=t1-0.25_dp*ro
      ubar=t2-0.25_dp*ro
      b=sqrt(1._dp-ro)
      xlp=0.5_dp*(1._dp+b)
      xlm=0.5_dp*(1._dp-b)
      vlpm=log(xlp/xlm)
      vlsm=log(4._dp/ro)
      vltm=log(4._dp*t1/ro)
      vlwm=log(4._dp*t2/ro)
      vlbl=log(b/xlm)
      vdw=ddilog(1._dp-ro/(4._dp*t2))-0.5_dp*vlwm**2
      vdt=ddilog(1._dp-ro/(4._dp*t1))-0.5_dp*vltm**2
      vdmp=ddilog(-xlm/xlp)
      vdmb=ddilog(-xlm/b)+0.5_dp*vlbl**2

c--- Q-Qbar and Qbar-Q contributions

      f1=(vlpm**2/2._dp-2._dp*vdmb-pisq/3._dp)/b
      f2=(-b*vlsm+vlpm**2/4._dp+vdmp+pisq/12._dp)/b**3
      f3=(-b**3*vlsm-3*b*vlsm+0.75_dp*vlpm**2
     & +3._dp*vdmp+pisq/4._dp+2._dp*b**3)/b**5
      f4 = (vlpm**2/4._dp+vdmp+pisq/12._dp)/b
      f5t1 = (vltm**2+vdt+pisq/6._dp)/t1**3
      f5t2 = (vlwm**2+vdw+pisq/6._dp)/t2**3

c---  Singular parts
c---  qqQQv is the lowest order matrix element squared without
c---  averaging over initial colors or spins:
c---  use _0 for O(1) piece, _1 for O(eps)
      qqQQv_0=gsq**2*V*(2._dp*t1**2+2._dp*t2**2+ro)
      if (scheme == 'tH-V') then
      qqQQv_1=gsq**2*V*(-2._dp)
      else
      qqQQv_1=0._dp
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
     &  gsq/V*qqQQv_0*(
     & -2._dp*V/2/xn*EPIN2-V/xn*EPIN
     & -3._dp*V/2/xn*EPIN+xn*(vltm+vlwm)*EPIN
     & -(0.5_dp*vlpm/b*(1._dp+b**2)+vlsm)/xn*EPIN)
     & +gsq/V*qqQQv_1*(
     & -2._dp*V/2/xn*EPIN-V/xn
     & -3._dp*V/2/xn+xn*(vltm+vlwm)
     & -(0.5_dp*vlpm/b*(1._dp+b**2)+vlsm)/xn)

      qqas=gsq/V*qqQQv_0*(XN4/xn*(vltm-vlwm))*EPIN
     &    +gsq/V*qqQQv_1*(XN4/xn*(vltm-vlwm))

C---  Finite symmetric part
       qqsv=
     & +gsq/V*qqQQv_0*(xn*(22._dp/9._dp-0.5_dp*vlsm*(vltm+vlwm)
     & +0.25_dp*vlsm**2+11._dp/3._dp*rmuom2+11._dp/6._dp*vlsm-pisq/6._dp)
     & +1._dp/xn*(+6._dp-pisq/3._dp-1.5_dp*vlsm+0.5_dp*vlsm**2-b*vlpm
     & -0.5_dp*(1._dp+b**2)*(pisq/b+f1))
     & +TR*XLF*(-20._dp/9._dp+4._dp/3._dp*vlsm-4._dp/3._dp*rmuom2)
     & +TR*(-20._dp/9._dp+4._dp/3._dp*vlpm*b*(1._dp+ro/2._dp)-4._dp/3._dp*ro))
     & +gsq**3*0.5_dp*xn*(t1-t2)**2*f3
     & +gsq**3*f2*(xn*(-5._dp*(t1**2+t2**2)+2._dp
     & +b**2*(0.5_dp+6._dp*t1**2+6._dp*t2**2)-b**4))
     & -gsq**3*0.5_dp/xn*vlpm/b*((t1-t2)**2+b**2)
     & +gsq**3*xn*((t2-t1)*(vdt-vdw)+0.5_dp*ro*(vdt+vdw)
     & +0.5_dp*vlsm**2*(t1**2 +t2**2)+1._dp-pisq/12._dp*ro
     & -vlsm*(vltm+vlwm)*(t1**2+t2**2)-vlsm*(vltm-vlwm)*(t1-t2)
     & -(vltm/TBAR+vlwm/UBAR)*(0.5_dp*ro-t1*t2)
     & -vlsm*(1.5_dp+2._dp*ro))

c--- extra finite terms in DR scheme
      if (scheme == 'dred') then
        qqsv=qqsv+0.5_dp*(xn-1._dp/xn)*gsq*qqQQv_0/V
      endif

C--- Finite antisymmetric part
      qqav=gsq/V*qqQQv_0*(-XN4/2._dp/xn*vlsm*(vltm-vlwm))
     & +gsq**3*XN4/xn*(-(t1-t2)*(vdt+vdw)+0.5_dp*ro*(vdt-vdw)
     & +(pisq/6._dp+vlsm*(3._dp-ro)+0.5_dp*vlsm**2-vlsm*(vltm+vlwm))*(t1-t2)
     & -vlsm*(vltm-vlwm)*(t1**2+t2**2)
     & -(vltm/tbar-vlwm/ubar)*(0.5_dp*ro-t1*t2)
     & +f2*(t1-t2)*(-1._dp+2._dp*b**2+b**4))

      qqsym=V*(qqss+qqsv)
      qqasy=V*(qqas+qqav)

c--- g-g contribution

c---  Singular parts
c---  ggQQv is the lowest order matrix element squared without
c---  averaging over initial colors or spins:
c---  use _0 for O(1) piece, _1 for O(eps) and _2 for O(eps**2)
      ggQQv_0=2._dp/xn*(V/t1/t2-2._dp*xnsq)
     & *((t1**2+t2**2)+ro-ro**2/4._dp/t1/t2)
      ggQQov_0=2._dp/xn*((xnsq-2._dp)/t1/t2-2._dp*xnsq)
     & *((t1**2+t2**2)+ro-ro**2/4._dp/t1/t2)
      ggQQsv_0=2._dp/xn*(1._dp/t1/t2)
     & *((t1**2+t2**2)+ro-ro**2/4._dp/t1/t2)
      if (scheme == 'tH-V') then
      ggQQv_1=2._dp/xn*(V/t1/t2-2._dp*xnsq)
     & *(-(t1**2+t2**2)-1._dp)
      ggQQv_2=2._dp/xn*(V/t1/t2-2._dp*xnsq)
      ggQQov_1=2._dp/xn*((xnsq-2._dp)/t1/t2-2._dp*xnsq)
     & *(-(t1**2+t2**2)-1._dp)
      ggQQsv_1=2._dp/xn*(1._dp/t1/t2)
     & *(-(t1**2+t2**2)-1._dp)
      else
      ggQQv_1=0._dp
      ggQQv_2=0._dp
      ggQQov_1=0._dp
      ggQQsv_1=0._dp
      endif

c--- These are the singular pieces, written in such a way that
c--- the limit EPINV -> 0 is smooth
      ggs=
     &  (-2._dp*xn)*(EPIN2*ggQQv_0+EPIN*ggQQv_1+ggQQv_2)
     & +(2._dp*(2._dp*TR/3._dp*XLF-11._dp/6._dp*xn)-V/xn)*(EPIN*ggQQv_0+ggQQv_1)
     & +(EPIN*ggQQsv_0+ggQQsv_1)*(1._dp+b**2)/b*vlpm*V/2._dp/xn
     & +(EPIN*ggQQov_0+ggQQov_1)*(vlsm*xn-(1._dp+b**2)/b*vlpm/xn/2._dp)
     & +2._dp*(EPIN*ggQQsv_0+ggQQsv_1)*(vlsm-vltm-vlwm)*xn
     & +4._dp*xnsq*EPIN*vltm
     & *(-ro**2/4._dp/t1**2+t2/t1+ro*t2/t1-2._dp*t2**2)
     & +4._dp*xnsq*EPIN*vlwm
     & *(-2._dp*t1**2-ro**2/4._dp/t2**2+t1/t2+ro*t1/t2)
      if (scheme == 'tH-V') then
c--- These are the extra finite pieces (akin to _1) for the last 4 lines
c--- and they appear to have already been added in the finite pieces
c--- below, so we subtract it later there
      ggs=ggs+(
     & +4._dp*xnsq*vltm
     & *(-2._dp*t2/t1+2._dp*t2**2)
     & +4._dp*xnsq*vlwm
     & *(-2._dp*t1/t2+2._dp*t1**2))
      endif

c--- extra finite terms in DR scheme
      if (scheme == 'dred') then
        ggs=ggs+xn/3._dp*ggQQv_0
      endif

C--- overall factor of V removed
      ggQQv=2._dp/xn*(V/t1/t2-2._dp*xnsq)
     & *(t1**2+t2**2+ro-ro**2/(4._dp*t1*t2))
      ggQQov=2._dp/xn*((xnsq-2._dp)/(t1*t2)-2._dp*xnsq)
     & *(t1**2+t2**2+ro-ro**2/(4._dp*t1*t2))
      ggQQsv=2._dp/(xn*t1*t2)
     & *(t1**2+t2**2+ro-ro**2/(4._dp*t1*t2))

      ggv1 = f3*(-2._dp*t2**2+1._dp/(t1*t2)/4._dp-2._dp*t1**2)*xnsq
     1   +f4*(24._dp*t1*t2-ro**2/(t1*t2)/2._dp+ro/(t1*t2)+11._dp/(4._dp*t1*t2)
     2  -4._dp*ro-17._dp)*xnsq+f2*(-20._dp*t1*t2+(-3._dp)/(2._dp*t1*t2)+11._dp)
     3  *xnsq-2*ggqqv*rmuom2*((-11._dp)*xn/6._dp+2._dp*tr*xlf/3._dp)+(b**2+1._dp)
     4  *ggqqsv*pisq*v/(b*xn)/2._dp-(b**2+1._dp)*ggqqov*pisq/(b*xn)/2._dp
     5  -f5t2*ro**2*v/xnsq/4._dp-f5t1*ro**2*v/xnsq/4._dp+f4*((1._dp-ro/2._dp)
     6  *ro**2*(t2**2+t1**2)/(t1**2*t2**2)+(-ro**3+5._dp*ro**2-ro-6._dp)
     7  /(t1*t2)+4._dp)/xnsq+f1*(-(1._dp-ro/2._dp)*ro**2*(t2**2+t1**2)/
     8  (t1**2*t2**2)/2._dp+(ro**3/2._dp-2._dp*ro**2+ro+2._dp)/(t1*t2)+2._dp*ro
     9  -4._dp)/xnsq+f1*(4._dp*t2**2+4._dp*ro*t1*t2-(1._dp-ro/2._dp)*ro**2/(t1*t2)
     :  +4._dp*t1**2-2._dp*ro**2+2._dp*ro)+f4*(-ro**3/(t1*t2)+4._dp*ro**2
     ;  /(t1*t2)-4._dp/(t1*t2)-2._dp*ro**2-2._dp*ro+6._dp)
      ggv2 = (-ro**3/t2/2._dp+2._dp*ro**2/t2-ro/t2-2._dp/t2-ro**3/t2**2/2._dp
     1   +ro**2/t2**2-ro**3/t1/2._dp+3._dp*ro**2/t1-4._dp/t1+2._dp)*vlpm*vlwm/
     2   (b*xnsq)+(-ro**3/t2/2._dp+3._dp*ro**2/t2-4._dp/t2-ro**3/t1/2._dp
     3   +2._dp*ro**2/t1-ro/t1-2/t1-ro**3/t1**2/2._dp+ro**2/t1**2+2._dp)
     4   *vlpm*vltm/(b*xnsq)+(ro**3/t2/2._dp+(-5._dp)*ro**2/(2._dp*t2)
     5   +ro/t2/2._dp+3._dp/t2+ro**3/t2**2/4._dp-ro**2/t2**2/2._dp+ro**3/t1/2._dp
     6   +(-5._dp)*ro**2/(2._dp*t1)+ro/t1/2._dp+3/t1+ro**3/t1**2/4._dp
     7   -ro**2/t1**2/2._dp-2._dp)*vlpm*vlsm/(b*xnsq)+(ro/(t1*t2)-4._dp)*
     8   vlpm/(b*xnsq)
      ggv3 = (ro**3/t2-ro**2/t2+ro**3/t1-ro**2/t1-4._dp*ro**3+4._dp*ro**2)
     1   *tr*vlpm*xn/b+(-ro**3/t2/2._dp+ro**2/t2-2._dp*t1-ro**3/t1/2._dp
     2   +3._dp*ro**2/t1-4._dp/t1-ro**2-ro+4._dp)*vlpm*vlwm/b+(-ro**3/t2/2._dp
     3   +3._dp*ro**2/t2-4._dp/t2+2._dp*t1-ro**3/t1/2._dp+ro**2/t1-ro**2-ro+2._dp)
     4   *vlpm*vltm/b+(ro**3/t2/2._dp-2*ro**2/t2+2/t2+ro**3/t1/2._dp
     5   -2._dp*ro**2/t1+2._dp/t1+ro**2+ro-3._dp)*vlpm*vlsm/b+(ro**2/t2/2._dp
     6   -ro/t2/4._dp-8._dp*ro*t1**2+12._dp*t1**2+8._dp*ro*t1-12._dp*t1
     7   +ro**2/t1/2._dp-ro/t1/4._dp-2._dp*ro**2+ro+1._dp)*vlpm/b
      ggv4 = (5._dp*t1**2-2._dp*ro*t1-8._dp*t1+7._dp*ro**2/8._dp+3._dp*ro-1._dp)
     1   *vltm*xnsq/tbar+(-ro*t1**2/4._dp+(-3._dp)*ro**2*t1/8._dp+ro*t1/4._dp
     2   +3._dp*ro**3/32._dp-ro**2/8._dp)*vltm*xnsq/tbar**2+(8._dp*t1**2
     3   -16._dp*t1+2._dp*ro/t1-4._dp/t1-ro**2/t1**2+16._dp)*vltm*xnsq
     4   +(-14._dp*ro/(3._dp*t2)-7._dp/(2._dp*t2)+ro**2/t2**2+16._dp*t1**2
     5   -16._dp*t1+(-14._dp)*ro/(3._dp*t1)+(-7._dp)/(2._dp*t1)+ro**2/t1**2
     6   +26._dp*ro/3._dp+15._dp)*xnsq-2._dp*ro**2/t2+17._dp*ro/(2._dp*t2)
     7   +8._dp/t2-2._dp*ro**2/t2**2+ro/t2**2-16._dp*t1**2+16._dp*t1
     8   -2._dp*ro**2/t1+17._dp*ro/(2._dp*t1)+8._dp/t1-2._dp*ro**2/t1**2
     9   +ro/t1**2-10._dp*ro-19._dp
      ggv5 = (-2._dp*ro/t2-2._dp/t2+ro**2/t2**2+8._dp*t1**2-8._dp*t1+2._dp*ro+6._dp)
     1  *vlsm*vlwm*xn**2+(2._dp*ro/t2-4._dp/t2-ro**2/t2**2+8._dp*t1**2+8._dp)
     2  *vlwm*xnsq+(8._dp*t1**2-8._dp*t1-2._dp*ro/t1-2._dp/t1+ro**2/t1**2+2._dp*ro
     3  +6._dp)*vlsm*vltm*xnsq+(-ro/t2/2._dp-1._dp/t2/2._dp-ro/t1/2._dp-1/t1/2._dp
     4  +ro+1._dp)*vlsm**2*xnsq+(ro/t2/2._dp-1._dp/t2/2._dp-8._dp*t1**2+8._dp*t1
     5  +ro/t1/2._dp-1._dp/t1/2._dp-2._dp*ro+2._dp)*vlsm*xnsq
      ggv6 = (5._dp*t1**2+2._dp*ro*t1-2._dp*t1+7._dp*ro**2/8._dp+ro-4._dp)*vlwm
     1  *xnsq/ubar+(-ro*t1**2/4._dp+3._dp*ro**2*t1/8._dp+ro*t1/4._dp
     2  +3._dp*ro**3/32._dp-ro**2/2._dp)*vlwm*xnsq/ubar**2+(2._dp*ro/t2+2._dp/t2
     3  -8._dp*t1-2._dp*ro+2._dp)*vdw*xnsq+(8._dp*t1+2._dp*ro/t1+2._dp/t1-2._dp*ro
     4  -6._dp)*vdt*xnsq+(ro/2._dp-ro*t1/4._dp)*xnsq/ubar+(ro*t1/4._dp+ro/4._dp)
     5  *xnsq/tbar+pisq*(ro/t2/6._dp+1._dp/t2/6._dp-ro**2/t2**2/12._dp
     6  +(-4._dp)*t1**2/3._dp+4._dp*t1/3._dp+ro/t1/6._dp+1._dp/t1/6._dp
     7  -ro**2/t1**2/12._dp-ro/3._dp-1._dp)*xnsq+(-ro**2/t2+ro/t2-2._dp/t2
     8  -ro**2/t1+3._dp*ro/t1+4._dp/t1-ro**2/t1**2+ro/t1**2)*vltm/xnsq
     9  +(2._dp*ro**2/t2-5._dp*ro/t2-5._dp/t2+ro**2/t2**2-ro/t2**2
     :   +2._dp*ro**2/t1-5._dp*ro/t1-5._dp/t1+ro**2/t1**2-ro/t1**2+6._dp)/xnsq
      ggv7 = (t1+4._dp/t1+ro-4._dp)*vlwm/(ubar*xnsq)+(3._dp*ro*t1/4._dp-2._dp/t1
     1  +ro**2/4._dp-ro/4._dp+2._dp)*vlwm/(ubar**2*xnsq)+(-ro**2/t2
     2  +3._dp*ro/t2+4/t2-ro**2/t2**2+ro/t2**2-ro**2/t1+ro/t1-2._dp/t1)
     3  *vlwm/xnsq+(-3._dp*ro**2/(4._dp*t2)+ro/t2+4._dp/t2
     4  -3._dp*ro**2/(4._dp*t1)-ro/t1+2._dp/t1+ro**2/t1**2/4._dp-2._dp)*vltm**2
     5  /xnsq+(4._dp/t2-t1+ro-3._dp)*vltm/(tbar*xnsq)+(-2._dp/t2-3._dp*ro*t1/4._dp
     6  +ro**2/4._dp+ro/2._dp+2._dp)*vltm/(tbar**2*xnsq)+(-ro**2/t2/4._dp
     7  +3._dp/(2._dp*t2)-ro**2/t1/4._dp+3._dp/(2._dp*t1)-1._dp)*vlpm**2/xn**2
      ggv8 = (-3._dp*ro**2/(4._dp*t2)-ro/t2+2._dp/t2+ro**2/t2**2/4._dp
     1  -3._dp*ro**2/(4._dp*t1)+ro/t1+4._dp/t1-2._dp)*vlwm**2/xnsq
     2  +(-3._dp*ro**2/(4._dp*t2)-ro/t2+2._dp/t2+ro**2/t2**2/4._dp
     3  +(-3._dp)*ro**2/(4._dp*t1)+ro/t1+4._dp/t1-2._dp)
     3  *vdw/xnsq+((-3._dp)*ro**2/(4._dp*t2)+ro/t2+4._dp/t2
     4  -3._dp*ro**2/(4._dp*t1)-ro/t1+2._dp/t1+ro**2/t1**2/4._dp-2._dp)*vdt/xnsq
     5  +(2._dp/t1-ro/4._dp-2._dp)/(ubar*xnsq)+(2._dp/t2-ro/4._dp-2._dp)/(tbar*xnsq)
     6  +pisq*(-1._dp/t2/2._dp+ro**2/t2**2/24._dp-1._dp/t1/2._dp
     7  +ro**2/t1**2/24._dp+1._dp/3._dp)/xnsq
      ggv9 = (-2._dp*ro**2/t2+4._dp*ro/t2+4._dp/t2-2._dp*ro**2/t2**2
     1   -2._dp*ro**2/t1+4._dp*ro/t1+4._dp/t1-2._dp*ro**2/t1**2-8._dp)*vltm*vlwm
     2   +(-ro**2/t2+2._dp*ro/t2+2._dp/t2-2._dp*t1+ro/t1-2._dp/t1-2._dp*ro+2._dp)
     3   *vltm**2+(-12._dp/t2-t1**2+3._dp*ro*t1+9._dp*t1-7._dp*ro**2/8._dp
     4   -4._dp*ro+12._dp)*vltm/tbar+(2._dp/t2-3._dp*ro*t1**2/4._dp+5._dp*ro**2*t1
     4   /8._dp+5._dp*ro*t1/2._dp-3._dp*ro**3/32._dp-5._dp*ro**2/8._dp-ro/2._dp-2._dp)
     5   *vltm/tbar**2+(ro**2/t2+10._dp/t2+ro**2/t1-4._dp*ro/t1-8._dp/t1
     6   +2._dp*ro**2/t1**2-ro/t1**2)*vltm+(-ro**2/t2/8._dp+ro/t2/4._dp
     7   +1._dp/t2-ro**2/t1/8._dp+ro/t1/4._dp+1._dp/t1-ro+(-3._dp)/2._dp)*vlpm**2
      ggv10 = (ro/t2-2._dp/t2+2._dp*t1-ro**2/t1+2._dp*ro/t1+2._dp/t1-2._dp*ro)
     1  *vlwm**2+(-t1**2-3._dp*ro*t1-7*t1-12/t1-7._dp*ro**2/8._dp-ro+20._dp)
     2  *vlwm/ubar+(-3._dp*ro*t1**2/4._dp-5._dp*ro**2*t1/8._dp-ro*t1+2._dp/t1
     3  -3._dp*ro**3/32._dp+5._dp*ro/4._dp-2._dp)*vlwm/ubar**2+(ro**2/t2-4._dp*ro/t2
     4  -8._dp/t2+2._dp*ro**2/t2**2-ro/t2**2+ro**2/t1+10._dp/t1)*vlwm
     5  +(ro**2/t2-ro/t2-4._dp/t2+2._dp*t1-2._dp*ro+4._dp)*vdw
     6  +(-2._dp*t1+ro**2/t1-ro/t1-4._dp/t1-2._dp*ro+6._dp)*vdt
      ggv11 = (ro/t2/3._dp+ro/t1/3._dp-4._dp*ro/3._dp)*tr*xlf*xn
     1 +(ro**2/t2/4._dp+ro**2/t1/4._dp-ro**2)*tr*vlpm**2*xn
     2 +(2._dp*ro**2/t2+ro/t2/3._dp+2._dp*ro**2/t1+ro/t1/3._dp-8._dp*ro**2
     3 -4._dp*ro/3._dp)*tr*xn+pisq*(-ro**2/t2/4._dp-ro**2/t1/4._dp+ro**2)*tr*xn
     4 +(-3._dp*ro*t1/4._dp-2._dp/t1-ro**2/4._dp+3._dp*ro/4._dp+2._dp)/ubar
     5 +(-2._dp/t2+3._dp*ro*t1/4._dp-ro**2/4._dp+2._dp)/tbar
     6 +pisq*(-ro**2/t2/24._dp+ro/t2/4._dp-1._dp/t2+ro**2/t2**2/3._dp
     7 -ro**2/t1/24._dp+ro/t1/4._dp-1._dp/t1+ro**2/t1**2/3._dp+ro/3._dp+11._dp/6._dp)


c--- This is subtracting the extra finite piece mentioned above
      ggv11=ggv11-(
     & +4._dp*xnsq*vltm
     & *(-2._dp*t2/t1+2._dp*t2**2)
     & +4._dp*xnsq*vlwm
     & *(-2._dp*t1/t2+2._dp*t1**2))

C---replace the overall factor of V which was removed
      ggsym=V*(ggs+ggv1+ggv2+ggv3+ggv4+ggv5
     &     +ggv6+ggv7+ggv8+ggv9+ggv10+ggv11)
      return
      end

