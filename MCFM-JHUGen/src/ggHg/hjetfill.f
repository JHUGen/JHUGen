      subroutine hjetfill(s,t,u,virtgg,virtqa,virtaq,virtqg,virtgq)
      implicit none
      include 'constants.f'
      include 'epinv.f'
      include 'scale.f'
      include 'b0.f'
      double precision virtgg,virtqa,virtaq,virtqg,virtgq,
     . logg,loqa,loaq,loqg,logq,ddilog,Li2s,Li2t,Li2u,
     . lnm,lns,lnt,lnu,ln2t,ln2u,mhsq,s,t,u,xlf,subuv,Delta
      
      mhsq=s+t+u

      xlf=dfloat(nf)
      Li2t=ddilog(t/mhsq)
      Li2u=ddilog(u/mhsq)
      Li2s=ddilog((s-mhsq)/s)
      lns=dlog(s/mhsq)
      lnt=dlog(-t/mhsq)
      lnu=dlog(-u/mhsq)
      lnm=dlog(musq/mhsq)
      ln2t=dlog((mhsq-t)/mhsq)
      ln2u=dlog((mhsq-u)/mhsq)

      logg=+V*xn*(mhsq**4+s**4+t**4+u**4)/(s*t*u)
      loqa=xn*cf/s*(t**2+u**2)
      loaq=loqa
      logq=-xn*cf/u*(s**2+t**2)
      loqg=-xn*cf/t*(s**2+u**2)

c--- UV counterterm in MS bar scheme. 
      subuv=-epinv*b0
      Delta=11d0
c--- See C.R.Schmidt, PLB (413) 391, eq. (16),(17)
c--- Factor of ason2pi included in gg_hg_v.f
C--- Three powers of as in Born --> 3      
      subuv=3d0*subuv+Delta
  
      virtgg=-3d0*epinv**2*xn*logg 
     . +epinv*xn*logg*(lns+lnt+lnu-3d0*lnm )
     . +xn*logg
     . *(2d0*(Li2t+Li2u+Li2s)
     . +lnm*(lns+lnt+lnu)-lns*lnt-lns*lnu-lnt*lnu
     . +0.5d0*(lns**2-lnt**2-lnu**2)-1.5d0*lnm**2
     . +2d0*(lnu*ln2u+lnt*ln2t)+4d0/3d0*pisq)
     . +V*xn*(xn-xlf)/3d0*mhsq*(1d0+mhsq/s+mhsq/t+mhsq/u)
     . +subuv*logg
 
      virtqa=+(-2d0*xn+1d0/xn)*loqa*epinv**2
     . -2d0/3d0*xlf*epinv*loqa
     . +epinv*xn*loqa*(13d0/6d0-2d0*lnm+lnt+lnu)
     . +epinv/xn*loqa*(1.5d0-lns+lnm)
     . +loqa*xlf*(-10d0/9d0+2d0/3d0*lns-2d0/3d0*lnm)
     . +xn*loqa* (40d0/9d0+Li2t+Li2u+2d0*Li2s-13d0/6d0*(lns-lnm)
     . +(lnm-lns)*(lnt+lnu)+lns**2-lnm**2-0.5d0*lnt**2-0.5d0*lnu**2
     . +lnt*ln2t+lnu*ln2u)
     . +loqa/xn
     . *(4d0-Li2t-Li2u-1.5d0*(lns-lnm)+0.5d0*(lns-lnm)**2
     . +lnt*lnu-lnt*ln2t-lnu*ln2u)
     . -4d0/3d0*pi**2/xn*loqa
     . -0.25d0*(xn**3-1d0/xn)*(t+u)
     . +subuv*loqa

      virtaq=(-2d0*xn+1d0/xn)*loaq*epinv**2
     . -2d0/3d0*xlf*epinv*loaq
     . +epinv*xn*loaq*(13d0/6d0-2d0*lnm+lnu+lnt)
     . +epinv/xn*loaq*(1.5d0-lns+lnm)
     . +loaq*xlf*(-10d0/9d0+2d0/3d0*lns-2d0/3d0*lnm)
     . +xn*loaq* (40d0/9d0+Li2u+Li2t+2d0*Li2s-13d0/6d0*(lns-lnm)
     . +(lnm-lns)*(lnu+lnt)+lns**2-lnm**2-0.5d0*lnu**2-0.5d0*lnt**2
     . +lnu*ln2u+lnt*ln2t)
     . +loaq/xn
     . *(4d0-Li2u-Li2t-1.5d0*(lns-lnm)+0.5d0*(lns-lnm)**2
     . +lnu*lnt-lnu*ln2u-lnt*ln2t)
     . -4d0/3d0*pi**2/xn*loaq
     . -0.25d0*(xn**3-1d0/xn)*(u+t)
     . +subuv*loaq
 

      virtgq=(-2d0*xn+1d0/xn)*epinv**2*logq
     . -2d0/3d0*xlf*epinv*logq
     . +epinv*xn*logq*(13d0/6d0+lns-2d0*lnm+lnt)
     . +epinv/xn*logq*(3d0/2d0+lnm-lnu)
     . +logq*xlf*(-10d0/9d0-2d0/3d0*lnm+2d0/3d0*lnu)
     . +xn*logq*(40d0/9d0+Li2t+2d0*Li2u+Li2s
     . +lns*lnm-lns*lnu-13d0/6d0*(lnu-lnm)
     . +lnm*lnt-lnm**2-lnt*lnu-0.5d0*lnt**2
     . +2d0*lnu*ln2u+lnt*ln2t)
     . +logq/xn*(4d0-Li2t-Li2s+lns*lnt+0.5d0*lnu**2-0.5d0*lns**2
     . -lnm*lnu+0.5d0*lnm**2-lnt*ln2t-1.5d0*(lnu-lnm))
     . +4d0/3d0*pi**2*xn*logq
     . +0.25d0*(xn**3-1d0/xn)*(t+s)
     . +subuv*logq

      virtqg=(-2d0*xn+1d0/xn)*epinv**2*loqg
     . -2d0/3d0*xlf*epinv*loqg
     . +epinv*xn*loqg*(13d0/6d0+lns-2d0*lnm+lnu)
     . +epinv/xn*loqg*(3d0/2d0+lnm-lnt)
     . +loqg*xlf*(-10d0/9d0-2d0/3d0*lnm+2d0/3d0*lnt)
     . +xn*loqg*(40d0/9d0+Li2u+2d0*Li2t+Li2s
     . +lns*lnm-lns*lnt-13d0/6d0*(lnt-lnm)
     . +lnm*lnu-lnm**2-lnu*lnt-0.5d0*lnu**2
     . +2d0*lnt*ln2t+lnu*ln2u)
     . +loqg/xn*(4d0-Li2u-Li2s+lns*lnu+0.5d0*lnt**2-0.5d0*lns**2
     . -lnm*lnt+0.5d0*lnm**2-lnu*ln2u-1.5d0*(lnt-lnm))
     . +4d0/3d0*pi**2*xn*loqg
     . +0.25d0*(xn**3-1d0/xn)*(u+s)
     . +subuv*loqg

      return
      end
