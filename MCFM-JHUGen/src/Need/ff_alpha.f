      double precision function ff_alpha(muj,muk)
c--- This is an implementation of the modifications to the final-final
c--- dipoles ff_qq, according to the formulae in Appendix A of:
c---
c---  G.~Bevilacqua, M.~Czakon, C.~G.~Papadopoulos, R.~Pittau and M.~Worek,
c---  ``Assault on the NLO Wishlist: pp -> tt bb,''
c---  JHEP {\bf 0909}, 109 (2009) [arXiv:0907.4723 [hep-ph]].
c---
c---  Note: according to the MCFM implementation of the final-final dipoles,
c---  alpha in the paper is replaced by (aff*yp) here
      implicit none
      include 'alfacut.f'
      double precision muj,muk,mujsq,muksq,fac,rt
      double precision Deik,Dcoll,yp,ddilog,aa,bb,cc,dd,xx,xp,xm,vjk

c--- some definitions, for speed      
      mujsq=muj**2
      muksq=muk**2
      fac=1d0-mujsq-muksq
      rt=dsqrt(1d0+mujsq**2+muksq**2-2d0*(mujsq+muksq+mujsq*muksq))

c--- basic definitions      
      yp=1d0-2d0*muk*(1d0-muk)/fac

      aa=2d0*muk/fac
      bb=2d0*(1d0-muk)/fac
      cc=bb*muk
      dd=fac/2d0
      
      xp=((1d0-muk)**2-mujsq+rt)/fac
      xm=((1d0-muk)**2-mujsq-rt)/fac
      vjk=rt/fac

c--- Equation (A.8) (remember alpha -> alpha*yp)
      xx=yp*(1d0-aff)+dsqrt(yp*(1d0-aff)*
     & (1d0/yp-yp*aff+4d0*mujsq*muksq/(mujsq-(1d0-muk)**2)/fac))

c--- Equation (A.9)      
      Deik=(
     & -ddilog((aa+xx)/(aa+xp))+ddilog(aa/(aa+xp))
     & +ddilog((xp-xx)/(xp-bb))-ddilog(xp/(xp-bb))
     & +ddilog((cc+xx)/(cc+xp))-ddilog(cc/(cc+xp))
     & +ddilog((xm-xx)/(xm+aa))-ddilog(xm/(xm+aa))
     & -ddilog((bb-xx)/(bb-xm))+ddilog(bb/(bb-xm))
     & -ddilog((xm-xx)/(xm+cc))+ddilog(xm/(xm+cc))
     & +ddilog((bb-xx)/(bb+aa))-ddilog(bb/(bb+aa))
     & -ddilog((cc+xx)/(cc-aa))+ddilog(cc/(cc-aa))
     & +dlog(cc+xx)*dlog((aa-cc)*(xp-xx)/(aa+xx)/(cc+xp))
     & -dlog(cc)*dlog((aa-cc)*xp/aa/(cc+xp))
     & +dlog(bb-xx)*dlog((aa+xx)*(xm-bb)/(aa+bb)/(xm-xx))
     & -dlog(bb)*dlog(aa*(xm-bb)/(aa+bb)/xm)
     & -dlog((aa+xx)*(bb-xp))*dlog(xp-xx)
     & +dlog(aa*(bb-xp))*dlog(xp)
     & +dlog(dd)*dlog((aa+xx)*xp*xm/aa/(xp-xx)/(xm-xx))
     & +dlog((xm-xx)/xm)*dlog((cc+xm)/(aa+xm))
     & +0.5d0*dlog((aa+xx)/aa)*dlog(aa*(aa+xx)*(aa+xp)**2)
     &  )/vjk
      
c--- Equation (A.20)  (remember alpha -> alpha*yp)
      Dcoll=1.5d0*(1d0+aff*yp)+1d0/(1d0-muk)-2d0*(2d0-2d0*mujsq-muk)/fac
     & +(1d0-aff*yp)*mujsq/2d0/(mujsq+aff*yp*fac)
     & -2d0*dlog(aff*yp*fac/((1d0-muk)**2-mujsq))
     & +(1d0+mujsq-muksq)/2d0/fac*dlog((mujsq+aff*yp*fac)/(1d0-muk)**2)
      
c--- the final dipole is equal to this combination (c.f. dipoles_mass.f)
      ff_alpha=2d0*Deik+Dcoll        
      
c      write(6,*) 'ff_alpha: 2d0*Deik   ',2d0*Deik
c      write(6,*) 'ff_alpha:    Dcoll   ',Dcoll
c      write(6,*) 'ff_alpha: ff_alpha   ',ff_alpha
c      pause
      
      return
      end
      
