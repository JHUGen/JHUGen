      subroutine qqb_dm_monojet_g(p,msq) 
      implicit none
      include 'types.f'
       
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h' 
      include 'dm_params.f' 
      include 'qcdcouple.f' 
      real(dp):: p(mxpart,4),msq(-nf:nf,-nf:nf) 
c      real(dp):: msqo
      real(dp):: qgXqbg1(2),qgXqbg2(2),qgXqbg3(2)
      real(dp):: qbgXqg1(2),qbgXqg2(2),qbgXqg3(2) 
     &     ,ggXqqb1(2),ggXqqb2(2),ggXqqb3(2) 
     &     ,gqXqbg1(2),gqXqbg2(2),gqXqbg3(2) 
     &     ,gqbXqg1(2),gqbXqg2(2),gqbXqg3(2) 
     &     ,qqbXgg1(2),qqbXgg2(2),qqbXgg3(2) 
     &     ,qbqXgg1(2),qbqXgg2(2),qbqXgg3(2) 
      real(dp):: fac
      complex(dp):: qR_a(2,2,2,2),qR_b(2,2,2,2),qbRb_a(2,2,2,2)
     &     ,qbRb_b(2,2,2,2),qq_a(2,2,2,2),qq_b(2,2,2,2)

      complex(dp):: qbqb_a(2,2,2,2),qbqb_b(2,2,2,2),qqb_a(2,2,2,2)
     &    ,qqb_b(2,2,2,2),qRb_a(2,2,2,2),qRb_b(2,2,2,2) 
      complex(dp):: qbR_a(2,2,2,2),qbR_b(2,2,2,2)
      complex(dp):: qbq_a(2,2,2,2),qbq_b(2,2,2,2)
      integer:: j,k,nq
      complex(dp):: a1111,a1211,a2111,a2211
      complex(dp):: a1112,a1212,a2112,a2212
      complex(dp):: a1121,a1221,a2121,a2221
      complex(dp):: a1122,a1222,a2122,a2222
      complex(dp):: b1111,b1211,b2111,b2211
      complex(dp):: b1112,b1212,b2112,b2212
      complex(dp):: b1121,b1221,b2121,b2221
      complex(dp):: b1122,b1222,b2122,b2222
      real(dp):: faclo,s34
      complex(dp):: cprop
      real(dp):: propsq

      if(dm_mediator=='gluonO') then 
         call gg_dm_monojet_g(p,msq) 
         return 
      endif
      

!      check_QED=.false. 
!      s34=(p(3,4)*p(4,4)-p(3,3)*p(4,3)-p(3,2)*p(4,2)-p(3,1)*p(4,1))*2d0
!      if(check_QED) then 
!         dm_lam=sqrt(s34/esq)
!         do j=1,nf 
!            dmL(j)=Q(j) 
!            dmR(j)=Q(j)
!         enddo
!      endif
         
!------ Fix Med-width 
!      if((medwidth==1d0).and.(first)) then 
!         medwidth=medmass/8d0/pi 
!      elseif((medwidth==0d0).and.(first)) then
!         medwidth=medmass/3d0
!      endif
!      if(first.and.(effective_th.eqv..false.)) then 
!         write(6,*) 'automatically set medwidth ',medwidth 
!      endif

      if(effective_th) then 
!--------- effective theory pre-factor
         fac=four*V*xn*gsq**2/dm_lam**4
         faclo=4d0*V*gsq**2*aveqq/dm_lam**4
       else
!-------- full theory => for V,A and PS is simply g_dmx**2*g_dmq**2/prop(34)**2 
!-------- for S and G need to be more careful, but G happens elsewhere and S can be done 
!-------- in special routine
         s34=(p(3,4)+p(4,4))**2
     &        -(p(3,1)+p(4,1))**2-(p(3,2)+p(4,2))**2-(p(3,3)+p(4,3))**2
         cprop=cone/cplx2((s34-medmass**2),medmass*medwidth)
         propsq=abs(cprop)**2 
         fac=four*V*xn*gsq**2*propsq*g_dmq**2*g_dmx**2
         faclo=four*V*gsq**2*aveqq*propsq*g_dmq**2*g_dmx**2
      endif

!-------two gluon amplitudes 
c      msqo=0d0

      if(dm_mediator=='vector') then 
         call qqb_dm_gg_Vamps(p,1,2,6,5,3,4,qgXqbg1,qgXqbg2,qgXqbg3) 
         call qqb_dm_gg_Vamps(p,5,2,6,1,3,4,qbgXqg1,qbgXqg2,qbgXqg3) 
         call qqb_dm_gg_Vamps(p,5,1,2,6,3,4,ggXqqb1,ggXqqb2,ggXqqb3) 
         call qqb_dm_gg_Vamps(p,2,1,6,5,3,4,gqXqbg1,gqXqbg2,gqXqbg3) 
         call qqb_dm_gg_Vamps(p,5,1,6,2,3,4,gqbXqg1,gqbXqg2,gqbXqg3) 
         call qqb_dm_gg_Vamps(p,1,5,6,2,3,4,qqbXgg1,qqbXgg2,qqbXgg3) 
         call qqb_dm_gg_Vamps(p,2,5,6,1,3,4,qbqXgg1,qbqXgg2,qbqXgg3) 
      elseif(dm_mediator=='axvect') then 
         call qqb_dm_gg_Axamps(p,1,2,6,5,3,4,qgXqbg1,qgXqbg2,qgXqbg3) 
         call qqb_dm_gg_Axamps(p,5,2,6,1,3,4,qbgXqg1,qbgXqg2,qbgXqg3) 
         call qqb_dm_gg_Axamps(p,5,1,2,6,3,4,ggXqqb1,ggXqqb2,ggXqqb3) 
         call qqb_dm_gg_Axamps(p,2,1,6,5,3,4,gqXqbg1,gqXqbg2,gqXqbg3) 
         call qqb_dm_gg_Axamps(p,5,1,6,2,3,4,gqbXqg1,gqbXqg2,gqbXqg3) 
         call qqb_dm_gg_Axamps(p,1,5,6,2,3,4,qqbXgg1,qqbXgg2,qqbXgg3) 
         call qqb_dm_gg_Axamps(p,2,5,6,1,3,4,qbqXgg1,qbqXgg2,qbqXgg3) 
      elseif(dm_mediator=='scalar') then 
         call qqb_dm_gg_Samps(p,1,2,6,5,3,4,qgXqbg1,qgXqbg2,qgXqbg3) 
         call qqb_dm_gg_Samps(p,5,2,6,1,3,4,qbgXqg1,qbgXqg2,qbgXqg3) 
         call qqb_dm_gg_Samps(p,5,1,2,6,3,4,ggXqqb1,ggXqqb2,ggXqqb3) 
         call qqb_dm_gg_Samps(p,2,1,6,5,3,4,gqXqbg1,gqXqbg2,gqXqbg3) 
         call qqb_dm_gg_Samps(p,5,1,6,2,3,4,gqbXqg1,gqbXqg2,gqbXqg3) 
         call qqb_dm_gg_Samps(p,1,5,6,2,3,4,qqbXgg1,qqbXgg2,qqbXgg3) 
         call qqb_dm_gg_Samps(p,2,5,6,1,3,4,qbqXgg1,qbqXgg2,qbqXgg3)  
         fac=fac/4d0
         faclo=faclo/4d0
      elseif(dm_mediator=='pseudo') then 
         call qqb_dm_gg_PSamps(p,1,2,6,5,3,4,qgXqbg1,qgXqbg2,qgXqbg3) 
         call qqb_dm_gg_PSamps(p,5,2,6,1,3,4,qbgXqg1,qbgXqg2,qbgXqg3) 
         call qqb_dm_gg_PSamps(p,5,1,2,6,3,4,ggXqqb1,ggXqqb2,ggXqqb3) 
         call qqb_dm_gg_PSamps(p,2,1,6,5,3,4,gqXqbg1,gqXqbg2,gqXqbg3) 
         call qqb_dm_gg_PSamps(p,5,1,6,2,3,4,gqbXqg1,gqbXqg2,gqbXqg3) 
         call qqb_dm_gg_PSamps(p,1,5,6,2,3,4,qqbXgg1,qqbXgg2,qqbXgg3) 
         call qqb_dm_gg_PSamps(p,2,5,6,1,3,4,qbqXgg1,qbqXgg2,qbqXgg3) 
         fac=fac/4d0
         faclo=faclo/4d0
      endif
    
!------- FOUR quark amplitudes 
!----- order is qqbRRb for my amplitude    
      
!-------- qRb-> qRb
      call qqb_dm_qqb(p,1,5,2,6,qRb_a,qRb_b) 
!-------- qR-> qR
      call qqb_dm_qqb(p,1,5,6,2,qR_a,qR_b) 
!-------- qbR-> qbR
      call qqb_dm_qqb(p,6,1,5,2,qbR_a,qbR_b) 
!-------- qbRb-> qRb
      call qqb_dm_qqb(p,5,1,2,6,qbRb_a,qbRb_b) 
!-------- qqb-> qqbq
      call qqb_dm_qqb(p,1,2,5,6,qqb_a,qqb_b) 
!-------- qbq-> qbq
      call qqb_dm_qqb(p,2,1,5,6,qbq_a,qbq_b) 
!-------- qq-> qq 
      call qqb_dm_qqb(p,1,6,5,2,qq_a,qq_b) 
!--------qbqb -> qbqb
      call qqb_dm_qqb(p,6,1,2,5,qbqb_a,qbqb_b) 
            
      
!------ assemble final matrix element 
       
!------- gluon pieces (plus initialization) 
!------- at the moment do not seperate by col struc 
      do j=-nf,nf 
         do k=-nf,nf 
            msq(j,k)=0d0 
            if((j.ne.0).and.(k.ne.0).and.(j.ne.-k)) goto 22 

!-------- gluon only pieces             
            if((j==0).and.(k==0)) then 
               do nq=1,nf 
                  msq(j,k)=msq(j,k)
     &    +avegg*fac*(ggXqqb1(1)+ggXqqb2(1)+ggXqqb3(1))*abs(dmL(nq))**2
     &    +avegg*fac*(ggXqqb1(2)+ggXqqb2(2)+ggXqqb3(2))*abs(dmR(nq))**2
               enddo
            elseif((j==0).and.(k<0)) then 
           msq(j,k)=aveqg*fac*(
     &           (gqbXqg1(1)+gqbXqg2(1)+gqbXqg3(1))*abs(dmL(abs(k)))**2
     &          +(gqbXqg1(2)+gqbXqg2(2)+gqbXqg3(2))*abs(dmR(abs(k)))**2)  
            elseif((j==0).and.(k>0)) then 
               msq(j,k)=aveqg*fac*(
     &            (gqXqbg1(1)+gqXqbg2(1)+gqXqbg3(1))*abs(dmL(k))**2
     &           +(gqXqbg1(2)+gqXqbg2(2)+gqXqbg3(2))*abs(dmR(k))**2)            
            elseif((j<0).and.(k==0)) then 
               msq(j,k)=aveqg*fac*(
     &         (qbgXqg1(1)+qbgXqg2(1)+qbgXqg3(1))*abs(dmL(abs(j)))**2
     &        +(qbgXqg1(2)+qbgXqg2(2)+qbgXqg3(2))*abs(dmR(abs(j)))**2)
            elseif((j>0).and.(k==0)) then 
               msq(j,k)=aveqg*fac*(
     &              (qgXqbg1(1)+qgXqbg2(1)+qgXqbg3(1))*abs(dmL(j))**2
     &             +(qgXqbg1(2)+qgXqbg2(2)+qgXqbg3(2))*abs(dmR(j))**2)
            elseif((j<0).and.(k>0).and.(j==-k)) then 
               msq(j,k)=half*aveqq*fac*(
     &              (qbqXgg1(1)+qbqXgg2(1)+qbqXgg3(1))*abs(dmL(k))**2
     &             +(qbqXgg1(2)+qbqXgg2(2)+qbqXgg3(2))*abs(dmR(k))**2)          
            elseif((j>0).and.(k<0).and.(j==-k)) then 
               msq(j,k)=half*aveqq*fac*(
     &              (qqbXgg1(1)+qqbXgg2(1)+qqbXgg3(1))*abs(dmL(j))**2
     &             +(qqbXgg1(2)+qqbXgg2(2)+qqbXgg3(2))*abs(dmR(j))**2)
            endif

 22         continue 
         enddo
      enddo
       
  
!========= four quark pieces, 
      
      do j=-nf,nf
         do k=-nf,nf 
            if((j>0).and.(k>0)) then 
               
               if(j.ne.k) then 
!--------------- DM = LL                  
                  a1111=dmL(j)*qR_a(1,1,1,1)+dmL(k)*qR_b(1,1,1,1)               
                  a1211=dmL(j)*qR_a(1,2,1,1)+dmR(k)*qR_b(1,2,1,1)
                  a2111=dmR(j)*qR_a(2,1,1,1)+dmL(k)*qR_b(2,1,1,1)
                  a2211=dmR(j)*qR_a(2,2,1,1)+dmR(k)*qR_b(2,2,1,1)
!--------------  DM = RR 
                  a1122=dmL(j)*qR_a(1,1,2,2)+dmL(k)*qR_b(1,1,2,2)
                  a1222=dmL(j)*qR_a(1,2,2,2)+dmR(k)*qR_b(1,2,2,2)
                  a2122=dmR(j)*qR_a(2,1,2,2)+dmL(k)*qR_b(2,1,2,2)
                  a2222=dmR(j)*qR_a(2,2,2,2)+dmR(k)*qR_b(2,2,2,2)
!--------------  DM = LR 
                  a1112=dmL(j)*qR_a(1,1,1,2)+dmL(k)*qR_b(1,1,1,2)
                  a1212=dmL(j)*qR_a(1,2,1,2)+dmR(k)*qR_b(1,2,1,2)
                  a2112=dmR(j)*qR_a(2,1,1,2)+dmL(k)*qR_b(2,1,1,2)
                  a2212=dmR(j)*qR_a(2,2,1,2)+dmR(k)*qR_b(2,2,1,2)
!--------------  DM = RL 
                  a1121=dmL(j)*qR_a(1,1,2,1)+dmL(k)*qR_b(1,1,2,1)
                  a1221=dmL(j)*qR_a(1,2,2,1)+dmR(k)*qR_b(1,2,2,1)
                  a2121=dmR(j)*qR_a(2,1,2,1)+dmL(k)*qR_b(2,1,2,1)
                  a2221=dmR(j)*qR_a(2,2,2,1)+dmR(k)*qR_b(2,2,2,1)
!------------------- LL                   
                  msq(j,k)=msq(j,k)
     &                 +faclo*(abs(a1111)**2+abs(a1211)**2
     &                 +abs(a2111)**2+abs(a2211)**2)
!------------------- RR 
                  msq(j,k)=msq(j,k)
     &                 +faclo*(abs(a1122)**2+abs(a1222)**2
     &                 +abs(a2122)**2+abs(a2222)**2)
!------------------- RL
                  msq(j,k)=msq(j,k)
     &                 +faclo*(abs(a1121)**2+abs(a1221)**2
     &                 +abs(a2121)**2+abs(a2221)**2)
!------------------- LR 
                  msq(j,k)=msq(j,k)
     &                 +faclo*(abs(a1112)**2+abs(a1212)**2
     &                 +abs(a2112)**2+abs(a2212)**2)
!=================== prefac 
                  msq(j,k)=msq(j,k) 
               
                  
               elseif( j==k) then 
!=========== diagrams which intefere 
!------------------DM LL 
                  a1111=dmL(j)*qR_a(1,1,1,1)+dmL(k)*qR_b(1,1,1,1)               
                  b1111=dmL(j)*qq_a(1,1,1,1)+dmL(k)*qq_b(1,1,1,1)               
                  a2211=dmR(j)*qR_a(2,2,1,1)+dmR(k)*qR_b(2,2,1,1)               
                  b2211=dmR(j)*qq_a(2,2,1,1)+dmR(k)*qq_b(2,2,1,1)               
!------------------DM RR
                  a1122=dmL(j)*qR_a(1,1,2,2)+dmL(k)*qR_b(1,1,2,2)               
                  b1122=dmL(j)*qq_a(1,1,2,2)+dmL(k)*qq_b(1,1,2,2)               
                  a2222=dmR(j)*qR_a(2,2,2,2)+dmR(k)*qR_b(2,2,2,2)               
                  b2222=dmR(j)*qq_a(2,2,2,2)+dmR(k)*qq_b(2,2,2,2)               
!------------------DM LR 
                  a1112=dmL(j)*qR_a(1,1,1,2)+dmL(k)*qR_b(1,1,1,2)               
                  b1112=dmL(j)*qq_a(1,1,1,2)+dmL(k)*qq_b(1,1,1,2)               
                  a2212=dmR(j)*qR_a(2,2,1,2)+dmR(k)*qR_b(2,2,1,2)               
                  b2212=dmR(j)*qq_a(2,2,1,2)+dmR(k)*qq_b(2,2,1,2)               
!------------------DM RL 
                  a1121=dmL(j)*qR_a(1,1,2,1)+dmL(k)*qR_b(1,1,2,1)               
                  b1121=dmL(j)*qq_a(1,1,2,1)+dmL(k)*qq_b(1,1,2,1)               
                  a2221=dmR(j)*qR_a(2,2,2,1)+dmR(k)*qR_b(2,2,2,1)               
                  b2221=dmR(j)*qq_a(2,2,2,1)+dmR(k)*qq_b(2,2,2,1)               
                  
!==================|sq| only pieces. 
!------------------ DM LL 
                  a1211=dmL(j)*qR_a(1,2,1,1)+dmR(k)*qR_b(1,2,1,1)               
                  b1211=dmL(j)*qq_a(1,2,1,1)+dmR(k)*qq_b(1,2,1,1)               
                  a2111=dmR(j)*qR_a(2,1,1,1)+dmL(k)*qR_b(2,1,1,1)               
                  b2111=dmR(j)*qq_a(2,1,1,1)+dmL(k)*qq_b(2,1,1,1)    
!------------------DM RR 
                  a1222=dmL(j)*qR_a(1,2,2,2)+dmR(k)*qR_b(1,2,2,2)               
                  b1222=dmL(j)*qq_a(1,2,2,2)+dmR(k)*qq_b(1,2,2,2)               
                  a2122=dmR(j)*qR_a(2,1,2,2)+dmL(k)*qR_b(2,1,2,2)               
                  b2122=dmR(j)*qq_a(2,1,2,2)+dmL(k)*qq_b(2,1,2,2)    
!------------------DM LR 
                  a1212=dmL(j)*qR_a(1,2,1,2)+dmR(k)*qR_b(1,2,1,2)               
                  b1212=dmL(j)*qq_a(1,2,1,2)+dmR(k)*qq_b(1,2,1,2)               
                  a2112=dmR(j)*qR_a(2,1,1,2)+dmL(k)*qR_b(2,1,1,2)               
                  b2112=dmR(j)*qq_a(2,1,1,2)+dmL(k)*qq_b(2,1,1,2)    
!------------------DM RL 
                  a1221=dmL(j)*qR_a(1,2,2,1)+dmR(k)*qR_b(1,2,2,1)               
                  b1221=dmL(j)*qq_a(1,2,2,1)+dmR(k)*qq_b(1,2,2,1)               
                  a2121=dmR(j)*qR_a(2,1,2,1)+dmL(k)*qR_b(2,1,2,1)               
                  b2121=dmR(j)*qq_a(2,1,2,1)+dmL(k)*qq_b(2,1,2,1)    
!----------inteference (suppressed by xn) 
                  msq(j,k)=msq(j,k)+half*faclo*(
     &              +Dble(a1111*conjg(b1111))+Dble(a2211*conjg(b2211))
     &              +Dble(a1112*conjg(b1112))+Dble(a2212*conjg(b2212))
     &              +Dble(a1121*conjg(b1121))+Dble(a2221*conjg(b2221))
     &      +Dble(a1122*conjg(b1122))+Dble(a2222*conjg(b2222)))*two/xn
!----------- |a^2| pieces 
               msq(j,k)=msq(j,k)+half*faclo*(
     &       abs(a1111)**2+abs(a2211)**2+abs(a1211)**2+abs(a2111)**2
     &      +abs(a1112)**2+abs(a2212)**2+abs(a1212)**2+abs(a2112)**2
     &      +abs(a1121)**2+abs(a2221)**2+abs(a1221)**2+abs(a2121)**2
     &      +abs(a1122)**2+abs(a2222)**2+abs(a1222)**2+abs(a2122)**2)
!----------- |b^2| pieces 
               msq(j,k)=msq(j,k)+half*faclo*(
     &       abs(b1111)**2+abs(b2211)**2+abs(b1211)**2+abs(b2111)**2
     &      +abs(b1112)**2+abs(b2212)**2+abs(b1212)**2+abs(b2112)**2
     &      +abs(b1121)**2+abs(b2221)**2+abs(b1221)**2+abs(b2121)**2
     &      +abs(b1122)**2+abs(b2222)**2+abs(b1222)**2+abs(b2122)**2)

            endif
         elseif((j<0).and.(k<0)) then 
            if(j.ne.k) then 
!--------------- DM = LL                  
                  a1111=dmL(-j)*qbRb_a(1,1,1,1)+dmL(-k)*qbRb_b(1,1,1,1)               
                  a1211=dmL(-j)*qbRb_a(1,2,1,1)+dmR(-k)*qbRb_b(1,2,1,1)
                  a2111=dmR(-j)*qbRb_a(2,1,1,1)+dmL(-k)*qbRb_b(2,1,1,1)
                  a2211=dmR(-j)*qbRb_a(2,2,1,1)+dmR(-k)*qbRb_b(2,2,1,1)
!--------------  DM = RR 
                  a1122=dmL(-j)*qbRb_a(1,1,2,2)+dmL(-k)*qbRb_b(1,1,2,2)
                  a1222=dmL(-j)*qbRb_a(1,2,2,2)+dmR(-k)*qbRb_b(1,2,2,2)
                  a2122=dmR(-j)*qbRb_a(2,1,2,2)+dmL(-k)*qbRb_b(2,1,2,2)
                  a2222=dmR(-j)*qbRb_a(2,2,2,2)+dmR(-k)*qbRb_b(2,2,2,2)
!--------------  DM = LR 
                  a1112=dmL(-j)*qbRb_a(1,1,1,2)+dmL(-k)*qbRb_b(1,1,1,2)
                  a1212=dmL(-j)*qbRb_a(1,2,1,2)+dmR(-k)*qbRb_b(1,2,1,2)
                  a2112=dmR(-j)*qbRb_a(2,1,1,2)+dmL(-k)*qbRb_b(2,1,1,2)
                  a2212=dmR(-j)*qbRb_a(2,2,1,2)+dmR(-k)*qbRb_b(2,2,1,2)
!--------------  DM = RL 
                  a1121=dmL(-j)*qbRb_a(1,1,2,1)+dmL(-k)*qbRb_b(1,1,2,1)
                  a1221=dmL(-j)*qbRb_a(1,2,2,1)+dmR(-k)*qbRb_b(1,2,2,1)
                  a2121=dmR(-j)*qbRb_a(2,1,2,1)+dmL(-k)*qbRb_b(2,1,2,1)
                  a2221=dmR(-j)*qbRb_a(2,2,2,1)+dmR(-k)*qbRb_b(2,2,2,1)
!------------------- LL                   
                  msq(j,k)=msq(j,k)
     &                 +faclo*(abs(a1111)**2+abs(a1211)**2
     &                 +abs(a2111)**2+abs(a2211)**2)
!------------------- RR 
                  msq(j,k)=msq(j,k)
     &                 +faclo*(abs(a1122)**2+abs(a1222)**2
     &                 +abs(a2122)**2+abs(a2222)**2)
!------------------- RL
                  msq(j,k)=msq(j,k)
     &                 +faclo*(abs(a1121)**2+abs(a1221)**2
     &                 +abs(a2121)**2+abs(a2221)**2)
!------------------- LR 
                  msq(j,k)=msq(j,k)
     &                 +faclo*(abs(a1112)**2+abs(a1212)**2
     &                 +abs(a2112)**2+abs(a2212)**2)
!=================== prefac 
                  msq(j,k)=msq(j,k) 
               
                  
               elseif( j==k) then 
!=========== diagrams which intefere 
!------------------DM LL 
                  a1111=dmL(-j)*qbRb_a(1,1,1,1)+dmL(-k)*qbRb_b(1,1,1,1)               
                  b1111=dmL(-j)*qbqb_a(1,1,1,1)+dmL(-k)*qbqb_b(1,1,1,1)               
                  a2211=dmR(-j)*qbRb_a(2,2,1,1)+dmR(-k)*qbRb_b(2,2,1,1)               
                  b2211=dmR(-j)*qbqb_a(2,2,1,1)+dmR(-k)*qbqb_b(2,2,1,1)               
!------------------DM RR
                  a1122=dmL(-j)*qbRb_a(1,1,2,2)+dmL(-k)*qbRb_b(1,1,2,2)               
                  b1122=dmL(-j)*qbqb_a(1,1,2,2)+dmL(-k)*qbqb_b(1,1,2,2)               
                  a2222=dmR(-j)*qbRb_a(2,2,2,2)+dmR(-k)*qbRb_b(2,2,2,2)               
                  b2222=dmR(-j)*qbqb_a(2,2,2,2)+dmR(-k)*qbqb_b(2,2,2,2)               
!------------------DM LR 
                  a1112=dmL(-j)*qbRb_a(1,1,1,2)+dmL(-k)*qbRb_b(1,1,1,2)               
                  b1112=dmL(-j)*qbqb_a(1,1,1,2)+dmL(-k)*qbqb_b(1,1,1,2)               
                  a2212=dmR(-j)*qbRb_a(2,2,1,2)+dmR(-k)*qbRb_b(2,2,1,2)               
                  b2212=dmR(-j)*qbqb_a(2,2,1,2)+dmR(-k)*qbqb_b(2,2,1,2)               
!------------------DM RL 
                  a1121=dmL(-j)*qbRb_a(1,1,2,1)+dmL(-k)*qbRb_b(1,1,2,1)               
                  b1121=dmL(-j)*qbqb_a(1,1,2,1)+dmL(-k)*qbqb_b(1,1,2,1)               
                  a2221=dmR(-j)*qbRb_a(2,2,2,1)+dmR(-k)*qbRb_b(2,2,2,1)               
                  b2221=dmR(-j)*qbqb_a(2,2,2,1)+dmR(-k)*qbqb_b(2,2,2,1)               
                  
!==================|sq| only pieces. 
!------------------ DM LL 
                  a1211=dmL(-j)*qbRb_a(1,2,1,1)+dmR(-k)*qbRb_b(1,2,1,1)               
                  b1211=dmL(-j)*qbqb_a(1,2,1,1)+dmR(-k)*qbqb_b(1,2,1,1)               
                  a2111=dmR(-j)*qbRb_a(2,1,1,1)+dmL(-k)*qbRb_b(2,1,1,1)               
                  b2111=dmR(-j)*qbqb_a(2,1,1,1)+dmL(-k)*qbqb_b(2,1,1,1)    
!------------------DM RR 
                  a1222=dmL(-j)*qbRb_a(1,2,2,2)+dmR(-k)*qbRb_b(1,2,2,2)               
                  b1222=dmL(-j)*qbqb_a(1,2,2,2)+dmR(-k)*qbqb_b(1,2,2,2)               
                  a2122=dmR(-j)*qbRb_a(2,1,2,2)+dmL(-k)*qbRb_b(2,1,2,2)               
                  b2122=dmR(-j)*qbqb_a(2,1,2,2)+dmL(-k)*qbqb_b(2,1,2,2)    
!------------------DM LR 
                  a1212=dmL(-j)*qbRb_a(1,2,1,2)+dmR(-k)*qbRb_b(1,2,1,2)               
                  b1212=dmL(-j)*qbqb_a(1,2,1,2)+dmR(-k)*qbqb_b(1,2,1,2)               
                  a2112=dmR(-j)*qbRb_a(2,1,1,2)+dmL(-k)*qbRb_b(2,1,1,2)               
                  b2112=dmR(-j)*qbqb_a(2,1,1,2)+dmL(-k)*qbqb_b(2,1,1,2)    
!------------------DM RL 
                  a1221=dmL(-j)*qbRb_a(1,2,2,1)+dmR(-k)*qbRb_b(1,2,2,1)               
                  b1221=dmL(-j)*qbqb_a(1,2,2,1)+dmR(-k)*qbqb_b(1,2,2,1)               
                  a2121=dmR(-j)*qbRb_a(2,1,2,1)+dmL(-k)*qbRb_b(2,1,2,1)               
                  b2121=dmR(-j)*qbqb_a(2,1,2,1)+dmL(-k)*qbqb_b(2,1,2,1)    
!----------inteference (suppressed by xn) 
                  msq(j,k)=msq(j,k)+half*faclo*(
     &              +Dble(a1111*conjg(b1111))+Dble(a2211*conjg(b2211))
     &              +Dble(a1112*conjg(b1112))+Dble(a2212*conjg(b2212))
     &              +Dble(a1121*conjg(b1121))+Dble(a2221*conjg(b2221))
     &      +Dble(a1122*conjg(b1122))+Dble(a2222*conjg(b2222)))*two/xn
!----------- |a^2| pieces 
               msq(j,k)=msq(j,k)+half*faclo*(
     &       abs(a1111)**2+abs(a2211)**2+abs(a1211)**2+abs(a2111)**2
     &      +abs(a1112)**2+abs(a2212)**2+abs(a1212)**2+abs(a2112)**2
     &      +abs(a1121)**2+abs(a2221)**2+abs(a1221)**2+abs(a2121)**2
     &      +abs(a1122)**2+abs(a2222)**2+abs(a1222)**2+abs(a2122)**2)
!----------- |b^2| pieces 
               msq(j,k)=msq(j,k)+half*faclo*(
     &       abs(b1111)**2+abs(b2211)**2+abs(b1211)**2+abs(b2111)**2
     &      +abs(b1112)**2+abs(b2212)**2+abs(b1212)**2+abs(b2112)**2
     &      +abs(b1121)**2+abs(b2221)**2+abs(b1221)**2+abs(b2121)**2
     &      +abs(b1122)**2+abs(b2222)**2+abs(b1222)**2+abs(b2122)**2)

            endif

!======= q qb case 
         elseif((j>0).and.(k<0)) then 
            if(j.ne.-k) then 
               !--------------- DM = LL                  
                  a1111=dmL(j)*qRb_a(1,1,1,1)+dmL(-k)*qRb_b(1,1,1,1)               
                  a1211=dmL(j)*qRb_a(1,2,1,1)+dmR(-k)*qRb_b(1,2,1,1)
                  a2111=dmR(j)*qRb_a(2,1,1,1)+dmL(-k)*qRb_b(2,1,1,1)
                  a2211=dmR(j)*qRb_a(2,2,1,1)+dmR(-k)*qRb_b(2,2,1,1)
!--------------  DM = RR 
                  a1122=dmL(j)*qRb_a(1,1,2,2)+dmL(-k)*qRb_b(1,1,2,2)
                  a1222=dmL(j)*qRb_a(1,2,2,2)+dmR(-k)*qRb_b(1,2,2,2)
                  a2122=dmR(j)*qRb_a(2,1,2,2)+dmL(-k)*qRb_b(2,1,2,2)
                  a2222=dmR(j)*qRb_a(2,2,2,2)+dmR(-k)*qRb_b(2,2,2,2)
!--------------  DM = LR 
                  a1112=dmL(j)*qRb_a(1,1,1,2)+dmL(-k)*qRb_b(1,1,1,2)
                  a1212=dmL(j)*qRb_a(1,2,1,2)+dmR(-k)*qRb_b(1,2,1,2)
                  a2112=dmR(j)*qRb_a(2,1,1,2)+dmL(-k)*qRb_b(2,1,1,2)
                  a2212=dmR(j)*qRb_a(2,2,1,2)+dmR(-k)*qRb_b(2,2,1,2)
!--------------  DM = RL 
                  a1121=dmL(j)*qRb_a(1,1,2,1)+dmL(-k)*qRb_b(1,1,2,1)
                  a1221=dmL(j)*qRb_a(1,2,2,1)+dmR(-k)*qRb_b(1,2,2,1)
                  a2121=dmR(j)*qRb_a(2,1,2,1)+dmL(-k)*qRb_b(2,1,2,1)
                  a2221=dmR(j)*qRb_a(2,2,2,1)+dmR(-k)*qRb_b(2,2,2,1)
                  
!------------------- LL                   
                  msq(j,k)=msq(j,k)
     &                 +faclo*(abs(a1111)**2+abs(a1211)**2
     &                 +abs(a2111)**2+abs(a2211)**2)
!------------------- RR 
                  msq(j,k)=msq(j,k)
     &                 +faclo*(abs(a1122)**2+abs(a1222)**2
     &                 +abs(a2122)**2+abs(a2222)**2)
!------------------- RL
                  msq(j,k)=msq(j,k)
     &                 +faclo*(abs(a1121)**2+abs(a1221)**2
     &                 +abs(a2121)**2+abs(a2221)**2)
!------------------- LR 
                  msq(j,k)=msq(j,k)
     &                 +faclo*(abs(a1112)**2+abs(a1212)**2
     &                 +abs(a2112)**2+abs(a2212)**2)
!=================== prefac 
                  msq(j,k)=msq(j,k) 

               elseif(j==-k) then 
!=========== diagrams which intefere 
!------------------DM LL 
                  a1111=dmL(j)*qRb_a(1,1,1,1)+dmL(-k)*qRb_b(1,1,1,1)               
                  b1111=dmL(j)*qqb_a(1,1,1,1)+dmL(-k)*qqb_b(1,1,1,1)               
                  a2211=dmR(j)*qRb_a(2,2,1,1)+dmR(-k)*qRb_b(2,2,1,1)               
                  b2211=dmR(j)*qqb_a(2,2,1,1)+dmR(-k)*qqb_b(2,2,1,1)               
!------------------DM RR
                  a1122=dmL(j)*qRb_a(1,1,2,2)+dmL(-k)*qRb_b(1,1,2,2)               
                  b1122=dmL(j)*qqb_a(1,1,2,2)+dmL(-k)*qqb_b(1,1,2,2)               
                  a2222=dmR(j)*qRb_a(2,2,2,2)+dmR(-k)*qRb_b(2,2,2,2)               
                  b2222=dmR(j)*qqb_a(2,2,2,2)+dmR(-k)*qqb_b(2,2,2,2)               
!------------------DM LR 
                  a1112=dmL(j)*qRb_a(1,1,1,2)+dmL(-k)*qRb_b(1,1,1,2)               
                  b1112=dmL(j)*qqb_a(1,1,1,2)+dmL(-k)*qqb_b(1,1,1,2)               
                  a2212=dmR(j)*qRb_a(2,2,1,2)+dmR(-k)*qRb_b(2,2,1,2)               
                  b2212=dmR(j)*qqb_a(2,2,1,2)+dmR(-k)*qqb_b(2,2,1,2)               
!------------------DM RL 
                  a1121=dmL(j)*qRb_a(1,1,2,1)+dmL(-k)*qRb_b(1,1,2,1)               
                  b1121=dmL(j)*qqb_a(1,1,2,1)+dmL(-k)*qqb_b(1,1,2,1)               
                  a2221=dmR(j)*qRb_a(2,2,2,1)+dmR(-k)*qRb_b(2,2,2,1)               
                  b2221=dmR(j)*qqb_a(2,2,2,1)+dmR(-k)*qqb_b(2,2,2,1)               
                  
!==================|sq| only pieces. 
!------------------ DM LL 
                  a1211=dmL(j)*qRb_a(1,2,1,1)+dmR(-k)*qRb_b(1,2,1,1)               
                  b1211=dmL(j)*qqb_a(1,2,1,1)+dmR(-k)*qqb_b(1,2,1,1)               
                  a2111=dmR(j)*qRb_a(2,1,1,1)+dmL(-k)*qRb_b(2,1,1,1)               
                  b2111=dmR(j)*qqb_a(2,1,1,1)+dmL(-k)*qqb_b(2,1,1,1)    
!------------------DM RR 
                  a1222=dmL(j)*qRb_a(1,2,2,2)+dmR(-k)*qRb_b(1,2,2,2)               
                  b1222=dmL(j)*qqb_a(1,2,2,2)+dmR(-k)*qqb_b(1,2,2,2)               
                  a2122=dmR(j)*qRb_a(2,1,2,2)+dmL(-k)*qRb_b(2,1,2,2)               
                  b2122=dmR(j)*qqb_a(2,1,2,2)+dmL(-k)*qqb_b(2,1,2,2)    
!------------------DM LR 
                  a1212=dmL(j)*qRb_a(1,2,1,2)+dmR(-k)*qRb_b(1,2,1,2)               
                  b1212=dmL(j)*qqb_a(1,2,1,2)+dmR(-k)*qqb_b(1,2,1,2)               
                  a2112=dmR(j)*qRb_a(2,1,1,2)+dmL(-k)*qRb_b(2,1,1,2)               
                  b2112=dmR(j)*qqb_a(2,1,1,2)+dmL(-k)*qqb_b(2,1,1,2)    
!------------------DM RL 
                  a1221=dmL(j)*qRb_a(1,2,2,1)+dmR(-k)*qRb_b(1,2,2,1)               
                  b1221=dmL(j)*qqb_a(1,2,2,1)+dmR(-k)*qqb_b(1,2,2,1)               
                  a2121=dmR(j)*qRb_a(2,1,2,1)+dmL(-k)*qRb_b(2,1,2,1)               
                  b2121=dmR(j)*qqb_a(2,1,2,1)+dmL(-k)*qqb_b(2,1,2,1)    
                  
!----------inteference (suppressed by xn) 
                  msq(j,k)=msq(j,k)+faclo*(
     &              Dble(a1111*conjg(b1111))+Dble(a2211*conjg(b2211))
     &              +Dble(a1112*conjg(b1112))+Dble(a2212*conjg(b2212))
     &              +Dble(a1121*conjg(b1121))+Dble(a2221*conjg(b2221))
     &      +Dble(a1122*conjg(b1122))+Dble(a2222*conjg(b2222)))*two/xn
!-----------|a^2| pieces 
                  msq(j,k)=msq(j,k)+faclo*(
     &       abs(a1111)**2+abs(a2211)**2+abs(a1211)**2+abs(a2111)**2
     &      +abs(a1112)**2+abs(a2212)**2+abs(a1212)**2+abs(a2112)**2
     &      +abs(a1121)**2+abs(a2221)**2+abs(a1221)**2+abs(a2121)**2
     &      +abs(a1122)**2+abs(a2222)**2+abs(a1222)**2+abs(a2122)**2)
!----------- |b^2| pieces 
               msq(j,k)=msq(j,k)+faclo*(
     &       abs(b1111)**2+abs(b2211)**2+abs(b1211)**2+abs(b2111)**2
     &      +abs(b1112)**2+abs(b2212)**2+abs(b1212)**2+abs(b2112)**2
     &      +abs(b1121)**2+abs(b2221)**2+abs(b1221)**2+abs(b2121)**2
     &      +abs(b1122)**2+abs(b2222)**2+abs(b1222)**2+abs(b2122)**2)

!-------------- s-channel anhilation into non-identical quark pair 
!------ note that our dm could be flavour dependent so we simply do a sum 
!------ over nq ne j 
                  b1111=0d0
                  b1211=0d0
                  b2111=0d0
                  b2211=0d0
                  b1122=0d0
                  b1222=0d0
                  b2122=0d0
                  b2222=0d0
                  b1112=0d0
                  b1212=0d0
                  b2112=0d0
                  b2212=0d0
                  b1121=0d0
                  b1221=0d0
                  b2121=0d0
                  b2221=0d0
               do nq=1,nf
                  if(nq.ne.j) then 
                     
!-------------------- LL 
                b1111=dmL(j)*qqb_a(1,1,1,1)+dmL(nq)*qqb_b(1,1,1,1)
                b1211=dmL(j)*qqb_a(1,2,1,1)+dmR(nq)*qqb_b(1,2,1,1)
                b2111=dmR(j)*qqb_a(2,1,1,1)+dmL(nq)*qqb_b(2,1,1,1)
                b2211=dmR(j)*qqb_a(2,2,1,1)+dmR(nq)*qqb_b(2,2,1,1)
!------------------- LR 
                b1112=+(dmL(j))*qqb_a(1,1,1,2)+dmL(nq)*qqb_b(1,1,1,2)
                b1212=+(dmL(j))*qqb_a(1,2,1,2)+dmR(nq)*qqb_b(1,2,1,2)
                b2112=+(dmR(j))*qqb_a(2,1,1,2)+dmL(nq)*qqb_b(2,1,1,2)
                b2212=+(dmR(j))*qqb_a(2,2,1,2)+dmR(nq)*qqb_b(2,2,1,2)
!-------------------- RL 
                b1121=+(dmL(j))*qqb_a(1,1,2,1)+dmL(nq)*qqb_b(1,1,2,1)
                b1221=+(dmL(j))*qqb_a(1,2,2,1)+dmR(nq)*qqb_b(1,2,2,1)
                b2121=+(dmR(j))*qqb_a(2,1,2,1)+dmL(nq)*qqb_b(2,1,2,1)
                b2221=+(dmR(j))*qqb_a(2,2,2,1)+dmR(nq)*qqb_b(2,2,2,1)
!---------------- RR 
                b1122=+(dmL(j))*qqb_a(1,1,2,2)+dmL(nq)*qqb_b(1,1,2,2)
                b1222=+(dmL(j))*qqb_a(1,2,2,2)+dmR(nq)*qqb_b(1,2,2,2)
                b2122=+(dmR(j))*qqb_a(2,1,2,2)+dmL(nq)*qqb_b(2,1,2,2)
                b2222=+(dmR(j))*qqb_a(2,2,2,2)+dmR(nq)*qqb_b(2,2,2,2)
                
                msq(j,k)=msq(j,k)+faclo*(
     &       abs(b1111)**2+abs(b2211)**2+abs(b1211)**2+abs(b2111)**2
     &      +abs(b1112)**2+abs(b2212)**2+abs(b1212)**2+abs(b2112)**2
     &      +abs(b1121)**2+abs(b2221)**2+abs(b1221)**2+abs(b2121)**2
     &      +abs(b1122)**2+abs(b2222)**2+abs(b1222)**2+abs(b2122)**2)
             endif
              enddo
               endif

!======= qb q case 
         elseif((j<0).and.(k>0)) then 
            if( j.ne.-k) then 
!--------------- DM = LL                  
                  a1111=dmL(-j)*qbR_a(1,1,1,1)+dmL(k)*qbR_b(1,1,1,1)               
                  a1211=dmL(-j)*qbR_a(1,2,1,1)+dmR(k)*qbR_b(1,2,1,1)
                  a2111=dmR(-j)*qbR_a(2,1,1,1)+dmL(k)*qbR_b(2,1,1,1)
                  a2211=dmR(-j)*qbR_a(2,2,1,1)+dmR(k)*qbR_b(2,2,1,1)
!--------------  DM = RR 
                  a1122=dmL(-j)*qbR_a(1,1,2,2)+dmL(k)*qbR_b(1,1,2,2)
                  a1222=dmL(-j)*qbR_a(1,2,2,2)+dmR(k)*qbR_b(1,2,2,2)
                  a2122=dmR(-j)*qbR_a(2,1,2,2)+dmL(k)*qbR_b(2,1,2,2)
                  a2222=dmR(-j)*qbR_a(2,2,2,2)+dmR(k)*qbR_b(2,2,2,2)
!--------------  DM = LR 
                  a1112=dmL(-j)*qbR_a(1,1,1,2)+dmL(k)*qbR_b(1,1,1,2)
                  a1212=dmL(-j)*qbR_a(1,2,1,2)+dmR(k)*qbR_b(1,2,1,2)
                  a2112=dmR(-j)*qbR_a(2,1,1,2)+dmL(k)*qbR_b(2,1,1,2)
                  a2212=dmR(-j)*qbR_a(2,2,1,2)+dmR(k)*qbR_b(2,2,1,2)
!--------------  DM = RL 
                  a1121=dmL(-j)*qbR_a(1,1,2,1)+dmL(k)*qbR_b(1,1,2,1)
                  a1221=dmL(-j)*qbR_a(1,2,2,1)+dmR(k)*qbR_b(1,2,2,1)
                  a2121=dmR(-j)*qbR_a(2,1,2,1)+dmL(k)*qbR_b(2,1,2,1)
                  a2221=dmR(-j)*qbR_a(2,2,2,1)+dmR(k)*qbR_b(2,2,2,1)
!------------------- LL                   
                  msq(j,k)=msq(j,k)
     &                 +faclo*(abs(a1111)**2+abs(a1211)**2
     &                 +abs(a2111)**2+abs(a2211)**2)
!------------------- RR 
                  msq(j,k)=msq(j,k)
     &                 +faclo*(abs(a1122)**2+abs(a1222)**2
     &                 +abs(a2122)**2+abs(a2222)**2)
!------------------- RL
                  msq(j,k)=msq(j,k)
     &                 +faclo*(abs(a1121)**2+abs(a1221)**2
     &                 +abs(a2121)**2+abs(a2221)**2)
!------------------- LR 
                  msq(j,k)=msq(j,k)
     &                 +faclo*(abs(a1112)**2+abs(a1212)**2
     &                 +abs(a2112)**2+abs(a2212)**2)
!=================== prefac 
                  msq(j,k)=msq(j,k) 
               elseif(j==-k) then 
!=========== diagrams which intefere 
!------------------DM LL 
                  a1111=dmL(-j)*qbR_a(1,1,1,1)+dmL(k)*qbR_b(1,1,1,1)               
                  b1111=dmL(-j)*qbq_a(1,1,1,1)+dmL(k)*qbq_b(1,1,1,1)               
                  a2211=dmR(-j)*qbR_a(2,2,1,1)+dmR(k)*qbR_b(2,2,1,1)               
                  b2211=dmR(-j)*qbq_a(2,2,1,1)+dmR(k)*qbq_b(2,2,1,1)               
!------------------DM RR
                  a1122=dmL(-j)*qbR_a(1,1,2,2)+dmL(k)*qbR_b(1,1,2,2)               
                  b1122=dmL(-j)*qbq_a(1,1,2,2)+dmL(k)*qbq_b(1,1,2,2)               
                  a2222=dmR(-j)*qbR_a(2,2,2,2)+dmR(k)*qbR_b(2,2,2,2)               
                  b2222=dmR(-j)*qbq_a(2,2,2,2)+dmR(k)*qbq_b(2,2,2,2)               
!------------------DM LR 
                  a1112=dmL(-j)*qbR_a(1,1,1,2)+dmL(k)*qbR_b(1,1,1,2)               
                  b1112=dmL(-j)*qbq_a(1,1,1,2)+dmL(k)*qbq_b(1,1,1,2)               
                  a2212=dmR(-j)*qbR_a(2,2,1,2)+dmR(k)*qbR_b(2,2,1,2)               
                  b2212=dmR(-j)*qbq_a(2,2,1,2)+dmR(k)*qbq_b(2,2,1,2)               
!------------------DM RL 
                  a1121=dmL(-j)*qbR_a(1,1,2,1)+dmL(k)*qbR_b(1,1,2,1)               
                  b1121=dmL(-j)*qbq_a(1,1,2,1)+dmL(k)*qbq_b(1,1,2,1)               
                  a2221=dmR(-j)*qbR_a(2,2,2,1)+dmR(k)*qbR_b(2,2,2,1)               
                  b2221=dmR(-j)*qbq_a(2,2,2,1)+dmR(k)*qbq_b(2,2,2,1)               
                  
!==================|sq| only pieces. 
!------------------ DM LL 
                  a1211=dmL(-j)*qbR_a(1,2,1,1)+dmR(k)*qbR_b(1,2,1,1)               
                  b1211=dmL(-j)*qbq_a(1,2,1,1)+dmR(k)*qbq_b(1,2,1,1)               
                  a2111=dmR(-j)*qbR_a(2,1,1,1)+dmL(k)*qbR_b(2,1,1,1)               
                  b2111=dmR(-j)*qbq_a(2,1,1,1)+dmL(k)*qbq_b(2,1,1,1)    

!------------------DM RR 
                  a1222=dmL(-j)*qbR_a(1,2,2,2)+dmR(k)*qbR_b(1,2,2,2)               
                  b1222=dmL(-j)*qbq_a(1,2,2,2)+dmR(k)*qbq_b(1,2,2,2)               
                  a2122=dmR(-j)*qbR_a(2,1,2,2)+dmL(k)*qbR_b(2,1,2,2)               
                  b2122=dmR(-j)*qbq_a(2,1,2,2)+dmL(k)*qbq_b(2,1,2,2)    
!------------------DM LR 
                  a1212=dmL(-j)*qbR_a(1,2,1,2)+dmR(k)*qbR_b(1,2,1,2)               
                  b1212=dmL(-j)*qbq_a(1,2,1,2)+dmR(k)*qbq_b(1,2,1,2)               
                  a2112=dmR(-j)*qbR_a(2,1,1,2)+dmL(k)*qbR_b(2,1,1,2)               
                  b2112=dmR(-j)*qbq_a(2,1,1,2)+dmL(k)*qbq_b(2,1,1,2)    
!------------------DM RL 
                  a1221=dmL(-j)*qbR_a(1,2,2,1)+dmR(k)*qbR_b(1,2,2,1)               
                  b1221=dmL(-j)*qbq_a(1,2,2,1)+dmR(k)*qbq_b(1,2,2,1)               
                  a2121=dmR(-j)*qbR_a(2,1,2,1)+dmL(k)*qbR_b(2,1,2,1)               
                  b2121=dmR(-j)*qbq_a(2,1,2,1)+dmL(k)*qbq_b(2,1,2,1)    
                  
!----------inteference (suppressed by xn) 
                  msq(j,k)=msq(j,k)+faclo*(
     &              Dble(a1111*conjg(b1111))+Dble(a2211*conjg(b2211))
     &              +Dble(a1112*conjg(b1112))+Dble(a2212*conjg(b2212))
     &              +Dble(a1121*conjg(b1121))+Dble(a2221*conjg(b2221))
     &      +Dble(a1122*conjg(b1122))+Dble(a2222*conjg(b2222)))*two/xn
!-----------|a^2| pieces 
                  msq(j,k)=msq(j,k)+faclo*(
     &       abs(a1111)**2+abs(a2211)**2+abs(a1211)**2+abs(a2111)**2
     &      +abs(a1112)**2+abs(a2212)**2+abs(a1212)**2+abs(a2112)**2
     &      +abs(a1121)**2+abs(a2221)**2+abs(a1221)**2+abs(a2121)**2
     &      +abs(a1122)**2+abs(a2222)**2+abs(a1222)**2+abs(a2122)**2)
!----------- |b^2| pieces 
               msq(j,k)=msq(j,k)+faclo*(
     &       abs(b1111)**2+abs(b2211)**2+abs(b1211)**2+abs(b2111)**2
     &      +abs(b1112)**2+abs(b2212)**2+abs(b1212)**2+abs(b2112)**2
     &      +abs(b1121)**2+abs(b2221)**2+abs(b1221)**2+abs(b2121)**2
     &      +abs(b1122)**2+abs(b2222)**2+abs(b1222)**2+abs(b2122)**2)

!-------------- s-channel anhilation into non-identical quark pair 
!------ note that our dm could be flavour dependent so we simply do a sum 
!------ over nq ne j 
                  b1111=0d0
                  b1211=0d0
                  b2111=0d0
                  b2211=0d0
                  b1122=0d0
                  b1222=0d0
                  b2122=0d0
                  b2222=0d0
                  b1112=0d0
                  b1212=0d0
                  b2112=0d0
                  b2212=0d0
                  b1121=0d0
                  b1221=0d0
                  b2121=0d0
                  b2221=0d0
               do nq=1,nf
                  if(nq.ne.k) then 
!-------------------- LL 
                b1111=+(dmL(k))*qbq_a(1,1,1,1)+dmL(nq)*qbq_b(1,1,1,1)
                b1211=+(dmL(k))*qbq_a(1,2,1,1)+dmR(nq)*qbq_b(1,2,1,1)
                b2111=+(dmR(k))*qbq_a(2,1,1,1)+dmL(nq)*qbq_b(2,1,1,1)
                b2211=+(dmR(k))*qbq_a(2,2,1,1)+dmR(nq)*qbq_b(2,2,1,1)
!------------------- LR 
                b1112=+(dmL(k))*qbq_a(1,1,1,2)+dmL(nq)*qbq_b(1,1,1,2)
                b1212=+(dmL(k))*qbq_a(1,2,1,2)+dmR(nq)*qbq_b(1,2,1,2)
                b2112=+(dmR(k))*qbq_a(2,1,1,2)+dmL(nq)*qbq_b(2,1,1,2)
                b2212=+(dmR(k))*qbq_a(2,2,1,2)+dmR(nq)*qbq_b(2,2,1,2)
!-------------------- RL 
                b1121=+(dmL(k))*qbq_a(1,1,2,1)+dmL(nq)*qbq_b(1,1,2,1)
                b1221=+(dmL(k))*qbq_a(1,2,2,1)+dmR(nq)*qbq_b(1,2,2,1)
                b2121=+(dmR(k))*qbq_a(2,1,2,1)+dmL(nq)*qbq_b(2,1,2,1)
                b2221=+(dmR(k))*qbq_a(2,2,2,1)+dmR(nq)*qbq_b(2,2,2,1)
!---------------- RR 
                b1122=+(dmL(k))*qbq_a(1,1,2,2)+dmL(nq)*qbq_b(1,1,2,2)
                b1222=+(dmL(k))*qbq_a(1,2,2,2)+dmR(nq)*qbq_b(1,2,2,2)
                b2122=+(dmR(k))*qbq_a(2,1,2,2)+dmL(nq)*qbq_b(2,1,2,2)
                b2222=+(dmR(k))*qbq_a(2,2,2,2)+dmR(nq)*qbq_b(2,2,2,2)
                
                msq(j,k)=msq(j,k)+faclo*(
     &       abs(b1111)**2+abs(b2211)**2+abs(b1211)**2+abs(b2111)**2
     &      +abs(b1112)**2+abs(b2212)**2+abs(b1212)**2+abs(b2112)**2
     &      +abs(b1121)**2+abs(b2221)**2+abs(b1221)**2+abs(b2121)**2
     &      +abs(b1122)**2+abs(b2222)**2+abs(b1222)**2+abs(b2122)**2)
             ENDIF                
             enddo
             endif
             
                
            endif

            
         enddo
      enddo



!      call spinoru(6,p,za,zb)
!      cor=4d0/(za(5,6)*zb(6,5))
!      cord=abs(cor)**2
      
!      call z2jetsq(1,4,5,6,2,3,za,zb,msq_test)
!      do h1=1,2
!         do h2=1,2
!            msqo=msqo+msq_test(h1,h2)
!         enddo
!      enddo
!      write(6,*) msqo
!      write(6,*)qqbgg1*cord,qqbgg2*cord,cord*qqbgg3,
!     & qqbgg1*cord+qqbgg2*cord+cord*qqbgg3
!      pause
      return 
      end 
