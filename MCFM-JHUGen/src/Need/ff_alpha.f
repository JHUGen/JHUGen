      function ff_alpha(muj,muk)
      implicit none
      include 'types.f'
      real(dp):: ff_alpha
c--- This is an implementation of the modifications to the final-final
c--- dipoles ff_qq, according to the formulae in Appendix A of:
c---
c---  G.~Bevilacqua, M.~Czakon, C.~G.~Papadopoulos, R.~Pittau and M.~Worek,
c---  ``Assault on the NLO Wishlist: pp -> tt bb,''
c---  JHEP {\bf 0909}, 109 (2009) [arXiv:0907.4723 [hep-ph]].
c---
c---  Note: according to the MCFM implementation of the final-final dipoles,
c---  alpha in the paper is replaced by (aff*yp) here
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'alfacut.f'
      real(dp):: muj,muk,mujsq,muksq,fac,rt
      real(dp):: Deik,Dcoll,yp,ddilog,aa,bb,cc,dd,xx,xp,xm,vjk

c--- some definitions, for speed      
      mujsq=muj**2
      muksq=muk**2
      fac=1._dp-mujsq-muksq
      rt=sqrt(1._dp+mujsq**2+muksq**2-2._dp*(mujsq+muksq+mujsq*muksq))

c--- basic definitions      
      yp=1._dp-2._dp*muk*(1._dp-muk)/fac

      aa=2._dp*muk/fac
      bb=2._dp*(1._dp-muk)/fac
      cc=bb*muk
      dd=fac/2._dp
      
      xp=((1._dp-muk)**2-mujsq+rt)/fac
      xm=((1._dp-muk)**2-mujsq-rt)/fac
      vjk=rt/fac

c--- Equation (A.8) (remember alpha -> alpha*yp)
      xx=yp*(1._dp-aff)+sqrt(yp*(1._dp-aff)*
     & (1._dp/yp-yp*aff+4._dp*mujsq*muksq/(mujsq-(1._dp-muk)**2)/fac))

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
     & +log(cc+xx)*log((aa-cc)*(xp-xx)/(aa+xx)/(cc+xp))
     & -log(cc)*log((aa-cc)*xp/aa/(cc+xp))
     & +log(bb-xx)*log((aa+xx)*(xm-bb)/(aa+bb)/(xm-xx))
     & -log(bb)*log(aa*(xm-bb)/(aa+bb)/xm)
     & -log((aa+xx)*(bb-xp))*log(xp-xx)
     & +log(aa*(bb-xp))*log(xp)
     & +log(dd)*log((aa+xx)*xp*xm/aa/(xp-xx)/(xm-xx))
     & +log((xm-xx)/xm)*log((cc+xm)/(aa+xm))
     & +0.5_dp*log((aa+xx)/aa)*log(aa*(aa+xx)*(aa+xp)**2)
     &  )/vjk
      
c--- Equation (A.20)  (remember alpha -> alpha*yp)
      Dcoll=1.5_dp*(1._dp+aff*yp)+1._dp/(1._dp-muk)-2._dp*(2._dp-2._dp*mujsq-muk)/fac
     & +(1._dp-aff*yp)*mujsq/2._dp/(mujsq+aff*yp*fac)
     & -2._dp*log(aff*yp*fac/((1._dp-muk)**2-mujsq))
     & +(1._dp+mujsq-muksq)/2._dp/fac*log((mujsq+aff*yp*fac)/(1._dp-muk)**2)
      
c--- the final dipole is equal to this combination (c.f. dipoles_mass.f)
      ff_alpha=2._dp*Deik+Dcoll        
      
c      write(6,*) 'ff_alpha: 2._dp*Deik   ',2._dp*Deik
c      write(6,*) 'ff_alpha:    Dcoll   ',Dcoll
c      write(6,*) 'ff_alpha: ff_alpha   ',ff_alpha
c      pause
      
      return
      end
      
