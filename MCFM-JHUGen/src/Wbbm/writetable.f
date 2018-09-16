      subroutine writetable(
     & A6treemp,A6treepm,A6treemm,A6treepp,
     & ALCmp,ALCpm,ALCmm,ALCpp,
     & ASLmp,ASLpm,ASLmm,ASLpp,
     & Afmp,Afpm,Afmm,Afpp)
      implicit none
      include 'types.f'
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      integer:: j
      complex(dp):: 
     & A6treemp,A6treepm,A6treemm,A6treepp,
     & ALCmp(-2:0),ALCpm(-2:0),ALCmm(-2:0),ALCpp(-2:0),
     & ASLmp(-2:0),ASLpm(-2:0),ASLmm(-2:0),ASLpp(-2:0),
     & Afmp(-2:0),Afpm(-2:0),Afmm(-2:0),Afpp(-2:0)
      real(dp):: p(mxpart,4),wmass,mb
      
      do j=-2,0
      ALCmp(j)=ALCmp(j)/A6treemp
      ALCpm(j)=ALCpm(j)/A6treepm
      ALCmm(j)=ALCmm(j)/A6treemm
      ALCpp(j)=ALCpp(j)/A6treepp
      if (abs(aimag(ALCmp(j)))<1e-15_dp) 
     & ALCmp(j)=cplx2(real(ALCmp(j)),zip)
      if (abs(aimag(ALCpm(j)))<1e-15_dp) 
     & ALCpm(j)=cplx2(real(ALCpm(j)),zip)
      if (abs(aimag(ALCmm(j)))<1e-15_dp) 
     & ALCmm(j)=cplx2(real(ALCmm(j)),zip)
      if (abs(aimag(ALCpp(j)))<1e-15_dp) 
     & ALCpp(j)=cplx2(real(ALCpp(j)),zip)

      ASLmp(j)=ASLmp(j)/A6treemp
      ASLpm(j)=ASLpm(j)/A6treepm
      ASLmm(j)=ASLmm(j)/A6treemm
      ASLpp(j)=ASLpp(j)/A6treepp
      if (abs(aimag(ASLmp(j)))<1e-15_dp) 
     & ASLmp(j)=cplx2(real(ASLmp(j)),zip)
      if (abs(aimag(ASLpm(j)))<1e-15_dp) 
     & ASLpm(j)=cplx2(real(ASLpm(j)),zip)
      if (abs(aimag(ASLmm(j)))<1e-15_dp) 
     & ASLmm(j)=cplx2(real(ASLmm(j)),zip)
      if (abs(aimag(ASLpp(j)))<1e-15_dp) 
     & ASLpp(j)=cplx2(real(ASLpp(j)),zip)

      Afmp(j)=Afmp(j)/A6treemp
      Afpm(j)=Afpm(j)/A6treepm
      Afmm(j)=Afmm(j)/A6treemm
      Afpp(j)=Afpp(j)/A6treepp
      if (abs(aimag(Afmp(j)))<1e-15_dp) 
     & Afmp(j)=cplx2(real(Afmp(j)),zip)
      if (abs(aimag(Afpm(j)))<1e-15_dp) 
     & Afpm(j)=cplx2(real(Afpm(j)),zip)
      if (abs(aimag(Afmm(j)))<1e-15_dp) 
     & Afmm(j)=cplx2(real(Afmm(j)),zip)
      if (abs(aimag(Afpp(j)))<1e-15_dp) 
     & Afpp(j)=cplx2(real(Afpp(j)),zip)
      enddo
      open(unit=66,file='table1.tex',status='unknown')
      open(unit=67,file='table2.tex',status='unknown')

      write(66,*) '\\begin{table}'
      write(66,*) '\\begin{tabular}{|c||c|c|c||c|c|c|}'
      write(66,*) '\\hline'
      write(66,*) '&\\multicolumn{3}{c||}{',
     & '$A_6(1_q^-,2_\\Qb^+,3_Q^-,4_',
     & '\\qb^+,5_{\\bar{\\ell}}^+,6_{\\ell}^-)$}'
      write(66,*) '&\\multicolumn{3}{c||}{',
     & '$A_6(1_q^-,2_\Qb^-,3_Q^+,4_\qb^+,',
     & '5_{\\bar{\\ell}}^+,6_{\\ell}^-)$}' 
      write(66,*) '\\\\'
      write(66,*) '\\hline'
      write(66,*) '\\hline'
      write(66,*) '& $1/\\e^2$ & $1/\\e$ & $\\e^0$'
      write(66,*) '&$1/\\e^2$ & $1/\e$ & $\\e^0$ \\\\'
      write(66,*) '\\hline'
      write(66,*) '\\hline $A_6^{\\rm tree}$ '
      write(66,88) zip,zip,real(A6treemp),zip,zip,real(A6treepm),
     &             zip,zip,aimag(A6treemp),zip,zip,aimag(A6treepm)
      write(66,*) '\\hline $A/A_6^{\\rm tree}$'
      write(66,88) 
     & real(ALCmp(-2)),
     & real(ALCmp(-1)),
     & real(ALCmp( 0)),
     & real(ALCpm(-2)),
     & real(ALCpm(-1)),
     & real(ALCpm( 0)),
     & aimag(ALCmp(-2)),
     & aimag(ALCmp(-1)),
     & aimag(ALCmp( 0)),
     & aimag(ALCpm(-2)),
     & aimag(ALCpm(-1)),
     & aimag(ALCpm( 0))
      write(66,*) '\\hline $A^{sl}/A_6^{\\rm tree}$'
      write(66,88) 
     & real(ASLmp(-2)),
     & real(ASLmp(-1)),
     & real(ASLmp( 0)),
     & real(ASLpm(-2)),
     & real(ASLpm(-1)),
     & real(ASLpm( 0)),
     & aimag(ASLmp(-2)),
     & aimag(ASLmp(-1)),
     & aimag(ASLmp( 0)),
     & aimag(ASLpm(-2)),
     & aimag(ASLpm(-1)),
     & aimag(ASLpm( 0))

      write(66,*) '\\hline $A^f/A_6^{\\rm tree}$ '
      write(66,88) 
     & real(Afmp(-2)),
     & real(Afmp(-1)),
     & real(Afmp( 0)),
     & real(Afpm(-2)),
     & real(Afpm(-1)),
     & real(Afpm( 0)),
     & aimag(Afmp(-2)),
     & aimag(Afmp(-1)),
     & aimag(Afmp( 0)),
     & aimag(Afpm(-2)),
     & aimag(Afpm(-1)),
     & aimag(Afpm( 0))

      write(66,*) '\\hline'

   88 format(6('& $',f19.11,'   $ '),'\\\\ \n',
     &                   6('& $',f19.11,'\\,i$ '),'\\\\ \n')   
      write(66,*) '\\hline'
      write(66,*) '\\end{tabular}'
      write(66,*) '\\caption{Numerical values',
     & ' of primitive amplitudes'
      write(66,*) 'at the kinematic point defined in'
      write(66,*) 'Eq.~(\\ref{kinpoint}).'
      write(66,*) '\\label{numres}}'
      write(66,*) '\\end{table}'


      write(67,*) '\\begin{table}'
      write(67,*) '\\begin{tabular}{|c||c|c|c||c|c|c|}'
      write(67,*) '\\hline'
      write(67,*) '&\\multicolumn{3}{c||}{',
     & '$A_6(1_q^-,2_\\Qb^-,3_Q^-,4_',
     & '\\qb^+,5_{\\bar{\\ell}}^+,6_{\\ell}^-)$}'
      write(67,*) '&\\multicolumn{3}{c||}{',
     & '$A_6(1_q^-,2_\Qb^+,3_Q^+,4_\qb^+,',
     & '5_{\\bar{\\ell}}^+,6_{\\ell}^-)$}' 
      write(67,*) '\\\\'
      write(67,*) '\\hline'
      write(67,*) '\\hline'
      write(67,*) '& $1/\\e^2$ & $1/\\e$ & $\\e^0$'
      write(67,*) '&$1/\\e^2$ & $1/\e$ & $\\e^0$ \\\\'
      write(67,*) '\\hline'
      write(67,*) '\\hline $A_6^{\\rm tree}$ '
      write(67,88) zip,zip,real(A6treemm),zip,zip,real(A6treepp),
     &             zip,zip,aimag(A6treemm),zip,zip,aimag(A6treepp)
      write(67,*) '\\hline $A/A_6^{\\rm tree}$'
      write(67,88) 
     & real(ALCmm(-2)),
     & real(ALCmm(-1)),
     & real(ALCmm( 0)),
     & real(ALCpp(-2)),
     & real(ALCpp(-1)),
     & real(ALCpp( 0)),
     & aimag(ALCmm(-2)),
     & aimag(ALCmm(-1)),
     & aimag(ALCmm( 0)),
     & aimag(ALCpp(-2)),
     & aimag(ALCpp(-1)),
     & aimag(ALCpp( 0))
      write(67,*) '\\hline $A^{sl}/A_6^{\\rm tree}$'
      write(67,88) 
     & real(ASLmm(-2)),
     & real(ASLmm(-1)),
     & real(ASLmm( 0)),
     & real(ASLpp(-2)),
     & real(ASLpp(-1)),
     & real(ASLpp( 0)),
     & aimag(ASLmm(-2)),
     & aimag(ASLmm(-1)),
     & aimag(ASLmm( 0)),
     & aimag(ASLpp(-2)),
     & aimag(ASLpp(-1)),
     & aimag(ASLpp( 0))

      write(67,*) '\\hline $A^f/A_6^{\\rm tree}$ '
      write(67,88) 
     & real(Afmm(-2)),
     & real(Afmm(-1)),
     & real(Afmm( 0)),
     & real(Afpp(-2)),
     & real(Afpp(-1)),
     & real(Afpp( 0)),
     & aimag(Afmm(-2)),
     & aimag(Afmm(-1)),
     & aimag(Afmm( 0)),
     & aimag(Afpp(-2)),
     & aimag(Afpp(-1)),
     & aimag(Afpp( 0))

      write(67,*) '\\hline'

      write(67,*) '\\hline'
      write(67,*) '\\end{tabular}'
      write(67,*) '\\caption{Numerical values',
     & ' of primitive amplitudes'
      write(67,*) 'at the kinematic point defined in'
      write(67,*) 'Eq.~(\\ref{kinpoint}).'
      write(67,*) '\\label{numres}}'
      write(67,*) '\\end{table}'

      include 'MCFMpoint.f'
      
      write(6,*) 'mb   =',mb
      write(6,*) 'wmass=',wmass
      do j=1,6
      write(6,77) j,p(j,4),p(j,1),p(j,2),p(j,3)
      enddo
      
      stop
      return
      
   77 format('p_',i1,' &=& (',3(f16.11,','),f19.11,'), \\nn \\\\')   
      
      end
