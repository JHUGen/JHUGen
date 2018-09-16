      subroutine qg_Hint_ZZ(p,msq)
      implicit none
      include 'types.f'

!==== C.W Oct 2013
!==== routine for calculating q(qb) g initiated inteferences in
!==== H=>ZZ
!===== returns only the interference
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'qlfirst.f'
      real(dp):: p(mxpart,4),msq(-nf:nf,-nf:nf)
      complex(dp):: qg_tH(2,2,2,2),qbg_tH(2,2,2,2),
     & gq_tH(2,2,2,2),gqb_tH(2,2,2,2)
      complex(dp):: qg_bH(2,2,2,2),qbg_bH(2,2,2,2),
     & gq_bH(2,2,2,2),gqb_bH(2,2,2,2)
      complex(dp):: qg_contu(2,2,2,2),qbg_contu(2,2,2,2),
     & gq_contu(2,2,2,2),gqb_contu(2,2,2,2)
      complex(dp):: qg_contd(2,2,2,2),qbg_contd(2,2,2,2),
     & gq_contd(2,2,2,2),gqb_contd(2,2,2,2)
      complex(dp):: A_higgs_XG,A_higgs_GX,A_cont_XG,A_cont_GX
      integer:: j,h1,h2,h34,h56,k
      real(dp):: fac

      if(qlfirst) then
         qlfirst=.false.
         call qlinit
      endif

      msq(:,:)=zip
      call spinoru(7,p,za,zb)

C-----Ordering of call is outgoing quark,outgoing antiquark,gluon
      call qg_HZZjet_amp(7,1,2,za,zb,qg_tH,qg_bH)
      call qg_HZZjet_amp(1,7,2,za,zb,qbg_tH,qbg_bH)
      call qg_HZZjet_amp(7,2,1,za,zb,gq_tH,gq_bH)
      call qg_HZZjet_amp(2,7,1,za,zb,gqb_tH,gqb_bH)

!===== continuum  amplitudes

C-----Ordering of call is outgoing quark,outgoing antiquark,gluon
      call qg_Cont_ZZj_amp(7,1,2,za,zb,qg_contu,qg_contd)
      call qg_Cont_ZZj_amp(1,7,2,za,zb,qbg_contu,qbg_contd)
      call qg_Cont_ZZj_amp(7,2,1,za,zb,gq_contu,gq_contd)
      call qg_Cont_ZZj_amp(2,7,1,za,zb,gqb_contu,gqb_contd)

!======== QB G and G QB summation :
      do j=-nf,-1


!======= down type pieces
         if(mod(abs(j),2)==1) then
            do h1=1,2
               do h2=1,2
                  do h34=1,2
                     do h56=1,2

                        A_higgs_XG=
     &                       qbg_tH(h1,h2,h34,h56)
     &                +      qbg_bH(h1,h2,h34,h56)
                        A_higgs_GX=
     &                       gqb_tH(h1,h2,h34,h56)
     &                +      gqb_bH(h1,h2,h34,h56)

                        A_cont_XG=
     &                       qbg_contd(h1,h2,h34,h56)
                        A_cont_GX=
     &                       gqb_contd(h1,h2,h34,h56)

                        msq(j,0)=msq(j,0)+abs(A_cont_XG+A_higgs_XG)**2
                        msq(0,j)=msq(0,j)+abs(A_cont_GX+A_higgs_GX)**2
!======== subtract S**2 and B**2
                        msq(j,0)=msq(j,0)
     &                       -abs(A_cont_XG)**2-abs(A_higgs_XG)**2
                        msq(0,j)=msq(0,j)
     &                       -abs(A_cont_GX)**2-abs(A_higgs_GX)**2
                     enddo
                  enddo
               enddo
            enddo
         else
!======= up type pieces
            do h1=1,2
               do h2=1,2
                  do h34=1,2
                     do h56=1,2

                        A_higgs_XG=
     &                       qbg_tH(h1,h2,h34,h56)
     &                +      qbg_bH(h1,h2,h34,h56)
                        A_higgs_GX=
     &                       gqb_tH(h1,h2,h34,h56)
     &                +      gqb_bH(h1,h2,h34,h56)

                        A_cont_XG=
     &                       qbg_contu(h1,h2,h34,h56)
                        A_cont_GX=
     &                       gqb_contu(h1,h2,h34,h56)

                        msq(j,0)=msq(j,0)+abs(A_cont_XG+A_higgs_XG)**2
                        msq(0,j)=msq(0,j)+abs(A_cont_GX+A_higgs_GX)**2
!======== subtract S**2 and B**2
                        msq(j,0)=msq(j,0)
     &                       -abs(A_cont_XG)**2-abs(A_higgs_XG)**2
                        msq(0,j)=msq(0,j)
     &                       -abs(A_cont_GX)**2-abs(A_higgs_GX)**2
                     enddo
                  enddo
               enddo
            enddo

         endif
      enddo
!==================================================================



!=============== Q G and G Q pieces :
      do j=1,nf

!=======down type pieces
         if(mod(j,2)==1) then

            do h1=1,2
               do h2=1,2
                  do h34=1,2
                     do h56=1,2

                        A_higgs_XG=
     &                       qg_tH(h1,h2,h34,h56)
     &                +      qg_bH(h1,h2,h34,h56)
                        A_higgs_GX=
     &                       gq_tH(h1,h2,h34,h56)
     &                +      gq_bH(h1,h2,h34,h56)

                        A_cont_XG=
     &                       qg_contd(h1,h2,h34,h56)
                        A_cont_GX=
     &                       gq_contd(h1,h2,h34,h56)

                        msq(j,0)=msq(j,0)+abs(A_cont_XG+A_higgs_XG)**2
                        msq(0,j)=msq(0,j)+abs(A_cont_GX+A_higgs_GX)**2
!======== subtract S**2 and B**2
                        msq(j,0)=msq(j,0)
     &                       -abs(A_cont_XG)**2-abs(A_higgs_XG)**2
                        msq(0,j)=msq(0,j)
     &                       -abs(A_cont_GX)**2-abs(A_higgs_GX)**2
                     enddo
                  enddo
               enddo
            enddo
         else
!=======up type pieces
            do h1=1,2
               do h2=1,2
                  do h34=1,2
                     do h56=1,2

                        A_higgs_XG=
     &                       qg_tH(h1,h2,h34,h56)
     &                +      qg_bH(h1,h2,h34,h56)
                        A_higgs_GX=
     &                       gq_tH(h1,h2,h34,h56)
     &                +      gq_bH(h1,h2,h34,h56)

                        A_cont_XG=
     &                       qg_contu(h1,h2,h34,h56)
                        A_cont_GX=
     &                       gq_contu(h1,h2,h34,h56)

                        msq(j,0)=msq(j,0)+abs(A_cont_XG+A_higgs_XG)**2
                        msq(0,j)=msq(0,j)+abs(A_cont_GX+A_higgs_GX)**2
!======== subtract S**2 and B**2
                        msq(j,0)=msq(j,0)
     &                       -abs(A_cont_XG)**2-abs(A_higgs_XG)**2
                        msq(0,j)=msq(0,j)
     &                       -abs(A_cont_GX)**2-abs(A_higgs_GX)**2
                     enddo
                  enddo
               enddo
            enddo
         endif
      enddo

!===== colour factor
      fac=aveqg*V/two

!==== rescale matrix elements by appropriate factor
      do j=-nf,nf
         do k=-nf,nf
            msq(j,k)=fac*msq(j,k)
         enddo
      enddo



      return
      end

