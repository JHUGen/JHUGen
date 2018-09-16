      subroutine branch(brwen,brzee,brznn,brtau,brtop,brcharm)
      implicit none
      include 'types.f'
      
C     Returns the lowest order branching ratios for 
C     1) W   --> e nu
C     2) Z   --> e e
C     3) Z   --> (nu nubar) x 3
C     4) tau --> e nu nubar
C     5) t   --> b W
C     6) c   --> s W (with Vcs omitted)
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'ewcouple.f'
      include 'zcouple.f'
      
      real(dp):: facz,facw,factau,factop,faccharm,
     & pwidth_e,pwidth_n
c      real(dp):: pwidth_u,pwidth_d,width
      real(dp):: brwen,brzee,brznn,brtau,brtop,brcharm,
     & xwsq,xbsq,root

      facz=esq/4._dp*zmass/(6._dp*pi)
      facw=gwsq/8._dp*wmass/(6._dp*pi)
      factau=gwsq**2/32._dp/wmass**4*mtau**5/192._dp/pi**3
c      xwsq=(wmass/mt)**2
c      xbsq=(mb/mt)**2
c      write(6,*) '(mb/mt)**2',xbsq
c      root=sqrt((1._dp+xbsq-xwsq)**2-4._dp*xbsq)
c      factop=(gw/wmass)**2*mt**3/(64._dp*pi)(1._dp-xwsq)**2*(1._dp+2._dp*xwsq)
c      factop=(gw/wmass)**2*mt**3/(64._dp*pi)*root
c     & *(1._dp+xwsq-2._dp*xwsq**2-2._dp*xbsq+xbsq*xwsq+xbsq**2)   
      faccharm=(gwsq/wmass**2)**2/32._dp*mc**5/192._dp/pi**3

      pwidth_e=facz*(le**2+re**2)
      pwidth_n=facz*(ln**2)*3._dp
c      pwidth_d=3*facz*(L(1)**2+R(1)**2)
c      pwidth_u=3*facz*(L(2)**2+R(2)**2)
c calculated zwidth=3*pwidth_.e+2_dp*pwidth_u+3*pwidth_e+3*pwidth_n
      brzee=pwidth_e/zwidth
      brznn=pwidth_n/zwidth
      brwen=facw/wwidth
      brtau=factau/tauwidth
c      brtop=factop/twidth
      brcharm=faccharm
      
      brtop=1._dp
      
      return
      end
