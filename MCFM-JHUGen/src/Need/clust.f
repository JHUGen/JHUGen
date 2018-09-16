      subroutine clust(jets,npar)
      implicit none
      include 'types.f'
      implicit real(dp) (a-h,o-z)
      parameter(pi=3.141592653589793238_dp)
      common /jetdef/ etminj,etmaxj,delrjj,rapmaxj,rapminj
      common /clusdef/ rsep,jalg1,jalg2
      common /jetcom/ icol,ji,jj,jk
      common /parmom/ ppar(4,10)
      common /jetmom/ pjet(8,10),jp(10)
      common /phypar/ w,ipp1,ipp2,rmw,rgw,rmz,rgz,sw2,qcdl
      data init/0/
      if(init==0) then
        init=1
        write(*,*)' jet clustering as of 12/8/95 '
      endif
*
* pjet(5,j) = ET
* pjet(6,j) = pseudorapidity
* pjet(7,j) = azimuthal angle
* pjet(8,j) = 0 (possible mass entry)
*


      do i=1,4+npar
         do j=1,4
            pjet(j,i)=ppar(j,i)
         enddo
      enddo
c--added by RKE
      pjet(8,5)=+1._dp
      pjet(8,6)=+1._dp
      pjet(8,7)=-1._dp
c--added by RKE
      if (npar==0) then
         jets=0
         return
      endif
      icol=0
      ji=-1
      jj=-1
      dij=w**2
      dib=w**2
      do i=1,npar
         jp(i)=i+4
         j=jp(i)
         pjet(5,j)=sqrt(pjet(1,j)**2+pjet(2,j)**2)
         theta=atan2(pjet(5,j),pjet(3,j))
         pjet(6,j)=-log(abs(tan(theta/2._dp)))
         pjet(7,j)=atan2(pjet(1,j),pjet(2,j))
      enddo
* cluster the partons
      if(npar>1)then
      do  i=1,npar-1
        do j=i+1,npar
          j1=jp(i)
          j2=jp(j)
          if ((pjet(4,j1)>0._dp).and.(pjet(4,j2)>0._dp)) then
*
* clustering criterion
*       jalg1 = 1 ; deltaR(i,j)   < delrjj
*       jalg1 = 2 ; deltaR(i,jet) < delrjj and deltaR(j,jet) < delrjj
*       jalg1 = 3 ; kt algorithm; R = delrjj
*       jalg1 = 4 ; deltaR(i,jet) < delrjj and deltaR(j,jet) < delrjj
*                      but deltaR(i,j) < Rsep
*
               if (jalg1==1) then
                  dely=pjet(6,j1)-pjet(6,j2)
                  rar=(pjet(1,j1)*pjet(1,j2)+pjet(2,j1)*pjet(2,j2))
     &                 /pjet(5,j1)/pjet(5,j2)
                  if (rar<-1._dp) then
                      delfi=pi
                  elseif (rar>1._dp) then
                      delfi=0._dp
                  else
                      delfi=acos(rar)
                  endif
                  delr=sqrt(dely**2+delfi**2)
                  if (delr<delrjj) then
                     icol=1
                     ji=j1
                     jj=j2
                  endif
               endif
               if (jalg1==2.or.jalg1==4) then
*
                 if(jalg1==4)then
                   dely=pjet(6,j1)-pjet(6,j2)
                   rar=(pjet(1,j1)*pjet(1,j2)+pjet(2,j1)*pjet(2,j2))
     &                 /pjet(5,j1)/pjet(5,j2)
                   if (rar<-1._dp) then
                      delfi=pi
                   elseif (rar>1._dp) then
                      delfi=0._dp
                   else
                      delfi=acos(rar)
                   endif
                   delr=sqrt(dely**2+delfi**2)
                 endif
                 if(jalg2==1.or.jalg2==3.or.jalg2==4)then
                   pt1=pjet(5,j1)
                   pt2=pjet(5,j2)
                   px=pjet(1,j1)+pjet(1,j2)
                   py=pjet(2,j1)+pjet(2,j2)
                   pz=pjet(3,j1)+pjet(3,j2)
                   ee=pjet(4,j1)+pjet(4,j2)
                   pt=sqrt(px**2+py**2)
                   theta=atan2(pt,pz)
                   if(jalg2==1)then
                     etjet=pt1+pt2
                   elseif(jalg2==3)then
                     etjet=pt
                   elseif(jalg2==4)then
                     etjet=ee*sin(theta)
                   endif
                   etajet=-log(abs(tan(theta/2.0_dp)))
                   phijet=atan2(px,py)
                   rar=(pjet(1,j1)*px+pjet(2,j1)*py)
     &                 /pt1/pt
                   if (rar<-1._dp) then
                      delfi1=pi
                   elseif (rar>1._dp) then
                      delfi1=0._dp
                   else
                      delfi1=acos(rar)
                   endif
                   rar=(pjet(1,j2)*px+pjet(2,j2)*py)
     &                 /pt2/pt
                   if (rar<-1._dp) then
                      delfi2=pi
                   elseif (rar>1._dp) then
                      delfi2=0._dp
                   else
                      delfi2=acos(rar)
                   endif
                 endif
                 if(jalg2==2)then
                   etjet=pjet(5,j1)+pjet(5,j2)
                   etajet=(pjet(6,j1)*pjet(5,j1)
     &                    +pjet(6,j2)*pjet(5,j2))/etjet
                   phijet=(pjet(7,j1)*pjet(5,j1)
     &                    +pjet(7,j2)*pjet(5,j2))/etjet
                   rar=(pjet(1,j1)*pjet(1,j2)+pjet(2,j1)*pjet(2,j2))
     &                 /pjet(5,j1)/pjet(5,j2)
                   if (rar<-1._dp) then
                      delfi=pi
                   elseif (rar>1._dp) then
                      delfi=0._dp
                   else
                      delfi=acos(rar)
                   endif
                   delfi1= pjet(5,j2)*delfi/etjet ! phi_1 - phijet
                   delfi2=-pjet(5,j1)*delfi/etjet ! phi_2 - phijet
                 endif
                 dely1=pjet(6,j1)-etajet
                 dely2=pjet(6,j2)-etajet
                 delr1=sqrt(dely1**2+delfi1**2)
                 delr2=sqrt(dely2**2+delfi2**2)
                 if ((delr1<delrjj).and.(delr2<delrjj)) then
                   if(jalg1==2)then
                     icol=1
                     ji=j1
                     jj=j2
                   elseif((jalg1==4).and.(delr<rsep)) then
                     icol=1
                     ji=j1
                     jj=j2
                   endif
                 endif
               endif
               if (jalg1==3) then
                  dely=pjet(6,j1)-pjet(6,j2)
                  rar=(pjet(1,j1)*pjet(1,j2)+pjet(2,j1)*pjet(2,j2))
     &                /pjet(5,j1)/pjet(5,j2)
                  if (rar<-1._dp) then
                     delfi=pi
                  elseif (rar>1._dp) then
                     delfi=0._dp
                  else
                     delfi=acos(rar)
                  endif
                  delr=sqrt(dely**2+delfi**2)
                  et=dmin1(pjet(5,j1),pjet(5,j2))
                  if(et**2*delr**2<dij)then
                    ji=j1
                    jj=j2
                    dij=et**2*delr**2
                  endif
                  if(et**2*delrjj**2<dib)dib=et**2*delrjj**2
               endif
            endif
         enddo
      enddo
      if(jalg1==3)then
         if(dij<dib)then
             icol=1
          endif
      endif
      endif
*
*
      if(icol==1)then
c----Added by RKE
             if ((ji == 5).or.(ji == 6)
     &       .or.(jj == 5).or.(jj == 6)) then
             pjet(8,ji)=+1._dp
             else
             pjet(8,ji)=-1._dp
             endif
        jk=ji
*    pjet(.,jk) is made of ppar(.,ji)+ppar(.,jj)
        if(jalg2==1.or.jalg2==3.or.jalg2==4)then
             do  k=1,4
               pjet(k,ji)=pjet(k,ji)+pjet(k,jj)
             enddo
             pp=sqrt(pjet(1,ji)**2+pjet(2,ji)**2)
             theta=atan2(pp,pjet(3,ji))
             if(jalg2==1)then
               pjet(5,ji)=pjet(5,ji)+pjet(5,jj)
               pjet(4,ji)=pjet(5,ji)/sin(theta)
             elseif(jalg2==3)then
               pjet(5,ji)=sqrt(pjet(1,ji)**2+pjet(2,ji)**2)
             elseif(jalg2==4)then
               pjet(5,ji)=pjet(4,ji)*sin(theta)
             endif
             pjet(6,ji)=-log(abs(tan(theta/2.0_dp)))
             pjet(7,ji)=atan2(pjet(1,ji),pjet(2,ji))
             pjet(4,jj)=-1._dp
        endif
        if(jalg2==2)then

             write(6,*) 'ji',ji
             write(6,*) 'jj',jj

             pause
             etjet=pjet(5,ji)+pjet(5,jj)
             etajet=
     &          (pjet(6,ji)*pjet(5,ji)+pjet(6,jj)*pjet(5,jj))/etjet
             phijet=
     &          (pjet(7,ji)*pjet(5,ji)+pjet(7,jj)*pjet(5,jj))/etjet
             eejet=exp(etajet)
             rar=(pjet(1,ji)*pjet(1,jj)+pjet(2,ji)*pjet(2,jj))
     &           /pjet(5,ji)/pjet(5,jj)
             if (rar<-1._dp) then
                  delfi=pi
             elseif (rar>1._dp) then
                  delfi=0._dp
             else
                  delfi=acos(rar)
             endif
             delfi1= pjet(5,jj)*delfi/etjet ! phi_1 - phijet
             cde=cos(delfi1)
             sde=sin(delfi1)
             if (pjet(2,ji)*pjet(1,jj)>pjet(2,jj)*pjet(1,ji))
     &       sde=-sde
             cphi=(pjet(1,ji)*cde-pjet(2,ji)*sde)/pjet(5,ji)
             sphi=(pjet(2,ji)*cde+pjet(1,ji)*sde)/pjet(5,ji)
             pjet(1,ji)=etjet*cphi
             pjet(2,ji)=etjet*sphi
             pjet(3,ji)=etjet/2._dp*(eejet-1._dp/eejet)
             pjet(4,ji)=etjet/2._dp*(eejet+1._dp/eejet)
             pjet(5,ji)=etjet
             pjet(6,ji)=etajet
             pjet(7,ji)=phijet
             pjet(4,jj)=-1._dp
        endif
      endif
*
* sort jets according to decreasing transverse energy
* j1=most energetic,...,j4=least energetic inside rapidity region
* then sort non-jets according to Et outside rap region
*
      do i=1,npar-1
         do j=i+1,npar
           j1=jp(i)
           j2=jp(j)
           if(pjet(4,j1)<0._dp)pjet(5,j1)=0._dp
           if(pjet(4,j2)<0._dp)pjet(5,j2)=0._dp
           if (pjet(5,j1)<pjet(5,j2)) then
               jt   =jp(i)
               jp(i)=jp(j)
               jp(j)=jt
           endif
         enddo
      enddo
* count number of observed jets
*      etminj  <  et < etmaxj
*      rapminj < eta < rapmaxj
      jets=0
      do i=1,npar
         j=jp(i)
         if (pjet(4,j)>0._dp) then
           if (pjet(5,j)>=etminj
     &        .and.pjet(5,j)<=etmaxj
     &        .and.abs(pjet(6,j))<=rapmaxj
     &        .and.abs(pjet(6,j))>=rapminj) then
              jets=jets+1
           else
              pjet(4,j)=-1._dp
           endif
         endif
      enddo
*
             write(6,*) 'end of clust:pjet(4,5)',pjet(4,5)
             write(6,*) 'end of clust:pjet(4,6)',pjet(4,6)
             write(6,*) 'end of clust:pjet(4,7)',pjet(4,7)
      end
*
************************************************************************
