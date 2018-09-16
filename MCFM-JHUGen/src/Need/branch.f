      subroutine branch(brwen,brzee,brznn,brtau,brtop,brcharm)
      implicit none
C     Returns the lowest order branching ratios for 
C     1) W   --> e nu
C     2) Z   --> e e
C     3) Z   --> (nu nubar) x 3
C     4) tau --> e nu nubar
C     5) t   --> b W
C     6) c   --> s W (with Vcs omitted)
      include 'constants.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'zcouple.f'
      
      double precision facz,facw,factau,factop,faccharm,
     . pwidth_e,pwidth_n
c      double precision pwidth_u,pwidth_d,width
      double precision brwen,brzee,brznn,brtau,brtop,brcharm,
     & xwsq,xbsq,root

      facz=esq/4d0*zmass/(6d0*pi)
      facw=gwsq/8d0*wmass/(6d0*pi)
      factau=gwsq**2/32d0/wmass**4*mtau**5/192d0/pi**3
c      xwsq=(wmass/mt)**2
c      xbsq=(mb/mt)**2
c      write(6,*) '(mb/mt)**2',xbsq
c      root=sqrt((1d0+xbsq-xwsq)**2-4d0*xbsq)
c      factop=(gw/wmass)**2*mt**3/(64d0*pi)(1d0-xwsq)**2*(1d0+2d0*xwsq)
c      factop=(gw/wmass)**2*mt**3/(64d0*pi)*root
c     & *(1d0+xwsq-2d0*xwsq**2-2d0*xbsq+xbsq*xwsq+xbsq**2)   
      faccharm=(gwsq/wmass**2)**2/32d0*mc**5/192d0/pi**3

      pwidth_e=facz*(le**2+re**2)
      pwidth_n=facz*(ln**2)*3d0
c      pwidth_d=3*facz*(L(1)**2+R(1)**2)
c      pwidth_u=3*facz*(L(2)**2+R(2)**2)
c calculated zwidth=3*pwidth_d+2*pwidth_u+3*pwidth_e+3*pwidth_n
      brzee=pwidth_e/zwidth
      brznn=pwidth_n/zwidth
      brwen=facw/wwidth
      brtau=factau/tauwidth
c      brtop=factop/twidth
      brcharm=faccharm
      
      brtop=1d0
      
      return
      end
