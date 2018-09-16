      subroutine hjetfill(s,t,u,virtgg,virtqa,virtaq,virtqg,virtgq)
      implicit none
      include 'types.f'

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'epinv.f'
      include 'scale.f'
      include 'b0.f'
      real(dp):: virtgg,virtqa,virtaq,virtqg,virtgq,
     & logg,loqa,loaq,loqg,logq,ddilog,Li2s,Li2t,Li2u,
     & lnm,lns,lnt,lnu,ln2t,ln2u,mhsq,s,t,u,xlf,subuv,Delta

      mhsq=s+t+u

      xlf=real(nf,dp)
      Li2t=ddilog(t/mhsq)
      Li2u=ddilog(u/mhsq)
      Li2s=ddilog((s-mhsq)/s)
      lns=log(s/mhsq)
      lnt=log(-t/mhsq)
      lnu=log(-u/mhsq)
      lnm=log(musq/mhsq)
      ln2t=log((mhsq-t)/mhsq)
      ln2u=log((mhsq-u)/mhsq)

      logg=+V*xn*(mhsq**4+s**4+t**4+u**4)/(s*t*u)
      loqa=xn*cf/s*(t**2+u**2)
      loaq=loqa
      logq=-xn*cf/u*(s**2+t**2)
      loqg=-xn*cf/t*(s**2+u**2)

c--- UV counterterm in MS bar scheme.
      subuv=-epinv*b0
      Delta=11._dp
c--- See C.R.Schmidt, PLB (413) 391, eq. (16),(17)
c--- Factor of ason2pi included in gg_hg_v.f
C--- Three powers of as in Born --> 3
      subuv=3._dp*subuv+Delta

      virtgg=-3._dp*epinv**2*xn*logg
     & +epinv*xn*logg*(lns+lnt+lnu-3._dp*lnm )
     & +xn*logg
     & *(2._dp*(Li2t+Li2u+Li2s)
     & +lnm*(lns+lnt+lnu)-lns*lnt-lns*lnu-lnt*lnu
     & +0.5_dp*(lns**2-lnt**2-lnu**2)-1.5_dp*lnm**2
     & +2._dp*(lnu*ln2u+lnt*ln2t)+4._dp/3._dp*pisq)
     & +V*xn*(xn-xlf)/3._dp*mhsq*(1._dp+mhsq/s+mhsq/t+mhsq/u)
     & +subuv*logg

      virtqa=+(-2._dp*xn+1._dp/xn)*loqa*epinv**2
     & -2._dp/3._dp*xlf*epinv*loqa
     & +epinv*xn*loqa*(13._dp/6._dp-2._dp*lnm+lnt+lnu)
     & +epinv/xn*loqa*(1.5_dp-lns+lnm)
     & +loqa*xlf*(-10._dp/9._dp+2._dp/3._dp*lns-2._dp/3._dp*lnm)
     & +xn*loqa* (40._dp/9._dp+Li2t+Li2u+2._dp*Li2s-13._dp/6._dp*(lns-lnm)
     & +(lnm-lns)*(lnt+lnu)+lns**2-lnm**2-0.5_dp*lnt**2-0.5_dp*lnu**2
     & +lnt*ln2t+lnu*ln2u)
     & +loqa/xn
     & *(4._dp-Li2t-Li2u-1.5_dp*(lns-lnm)+0.5_dp*(lns-lnm)**2
     & +lnt*lnu-lnt*ln2t-lnu*ln2u)
     & -4._dp/3._dp*pi**2/xn*loqa
     & -0.25_dp*(xn**3-1._dp/xn)*(t+u)
     & +subuv*loqa

      virtaq=(-2._dp*xn+1._dp/xn)*loaq*epinv**2
     & -2._dp/3._dp*xlf*epinv*loaq
     & +epinv*xn*loaq*(13._dp/6._dp-2._dp*lnm+lnu+lnt)
     & +epinv/xn*loaq*(1.5_dp-lns+lnm)
     & +loaq*xlf*(-10._dp/9._dp+2._dp/3._dp*lns-2._dp/3._dp*lnm)
     & +xn*loaq* (40._dp/9._dp+Li2u+Li2t+2._dp*Li2s-13._dp/6._dp*(lns-lnm)
     & +(lnm-lns)*(lnu+lnt)+lns**2-lnm**2-0.5_dp*lnu**2-0.5_dp*lnt**2
     & +lnu*ln2u+lnt*ln2t)
     & +loaq/xn
     & *(4._dp-Li2u-Li2t-1.5_dp*(lns-lnm)+0.5_dp*(lns-lnm)**2
     & +lnu*lnt-lnu*ln2u-lnt*ln2t)
     & -4._dp/3._dp*pi**2/xn*loaq
     & -0.25_dp*(xn**3-1._dp/xn)*(u+t)
     & +subuv*loaq


      virtgq=(-2._dp*xn+1._dp/xn)*epinv**2*logq
     & -2._dp/3._dp*xlf*epinv*logq
     & +epinv*xn*logq*(13._dp/6._dp+lns-2._dp*lnm+lnt)
     & +epinv/xn*logq*(3._dp/2._dp+lnm-lnu)
     & +logq*xlf*(-10._dp/9._dp-2._dp/3._dp*lnm+2._dp/3._dp*lnu)
     & +xn*logq*(40._dp/9._dp+Li2t+2._dp*Li2u+Li2s
     & +lns*lnm-lns*lnu-13._dp/6._dp*(lnu-lnm)
     & +lnm*lnt-lnm**2-lnt*lnu-0.5_dp*lnt**2
     & +2._dp*lnu*ln2u+lnt*ln2t)
     & +logq/xn*(4._dp-Li2t-Li2s+lns*lnt+0.5_dp*lnu**2-0.5_dp*lns**2
     & -lnm*lnu+0.5_dp*lnm**2-lnt*ln2t-1.5_dp*(lnu-lnm))
     & +4._dp/3._dp*pi**2*xn*logq
     & +0.25_dp*(xn**3-1._dp/xn)*(t+s)
     & +subuv*logq

      virtqg=(-2._dp*xn+1._dp/xn)*epinv**2*loqg
     & -2._dp/3._dp*xlf*epinv*loqg
     & +epinv*xn*loqg*(13._dp/6._dp+lns-2._dp*lnm+lnu)
     & +epinv/xn*loqg*(3._dp/2._dp+lnm-lnt)
     & +loqg*xlf*(-10._dp/9._dp-2._dp/3._dp*lnm+2._dp/3._dp*lnt)
     & +xn*loqg*(40._dp/9._dp+Li2u+2._dp*Li2t+Li2s
     & +lns*lnm-lns*lnt-13._dp/6._dp*(lnt-lnm)
     & +lnm*lnu-lnm**2-lnu*lnt-0.5_dp*lnu**2
     & +2._dp*lnt*ln2t+lnu*ln2u)
     & +loqg/xn*(4._dp-Li2u-Li2s+lns*lnu+0.5_dp*lnt**2-0.5_dp*lns**2
     & -lnm*lnt+0.5_dp*lnm**2-lnu*ln2u-1.5_dp*(lnt-lnm))
     & +4._dp/3._dp*pi**2*xn*loqg
     & +0.25_dp*(xn**3-1._dp/xn)*(u+s)
     & +subuv*loqg

      return
      end
