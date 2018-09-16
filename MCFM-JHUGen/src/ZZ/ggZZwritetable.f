      subroutine ggZZwritetable
      implicit none
      include 'types.f'
c--- routine to write out contents of res common block
      
      include 'blabels.f'
      include 'ZZclabels.f'
      include 'ZZdlabels.f'
      integer:: imt0,imt2,imt4,imp,ipp,iunit
      complex(dp):: res(2,4,10,3)
      parameter(imt0=1,imt2=2,imt4=3)
      parameter(ipp=1,imp=2)
      common/ggZZcaptureres/res
!$omp threadprivate(/ggZZcaptureres/)
      iunit=11

      if (iunit .ne. 6) then
        open(unit=iunit,file='coefftable.tex',status='unknown')
      endif
    
      write(iunit,*) '\begin{math}'
      write(iunit,*) '\begin{array}{| l || l || l |}'
      write(iunit,*) '\hline'
      write(iunit,*) ' & (1^+,2^+) & (1^-,2^+) \\'
      write(iunit,*) '\hline'
      
c--- 6-d boxes      
      write(iunit,90) 'd_2^{d=6}',
     & res(ipp,4,d6_2_1_34,imt0),res(imp,4,d6_2_1_34,imt0)
      write(iunit,*) '\hline' 
      write(iunit,90) 'd_3^{d=6}',
     & res(ipp,4,d6_1_2_34,imt0),res(imp,4,d6_1_2_34,imt0)
      write(iunit,*) '\hline' 

c--- 4-d boxes      
      write(iunit,94) 'd_1^{(0)}',
     & res(ipp,4,d1_34_2,imt0),res(imp,4,d1_34_2,imt0),
     & 'd_1^{(2)}',
     & res(ipp,4,d1_34_2,imt2),res(imp,4,d1_34_2,imt2),
     & 'd_1^{(4)}',
     & res(ipp,4,d1_34_2,imt4),res(imp,4,d1_34_2,imt4)
      write(iunit,*) '\hline' 
      write(iunit,94) 'd_2^{(0)}',
     & res(ipp,4,d2_1_34,imt0),res(imp,4,d2_1_34,imt0),
     & 'd_2^{(2)}',
     & res(ipp,4,d2_1_34,imt2),res(imp,4,d2_1_34,imt2),
     & 'd_2^{(4)}',
     & res(ipp,4,d2_1_34,imt4),res(imp,4,d2_1_34,imt4)
      write(iunit,*) '\hline' 
      write(iunit,94) 'd_3^{(0)}',
     & res(ipp,4,d1_2_34,imt0),res(imp,4,d1_2_34,imt0),
     & 'd_3^{(2)}',
     & res(ipp,4,d1_2_34,imt2),res(imp,4,d1_2_34,imt2),
     & 'd_3^{(4)}',
     & res(ipp,4,d1_2_34,imt4),res(imp,4,d1_2_34,imt4)
      write(iunit,*) '\hline' 

      write(iunit,*) '\hline' 

c--- triangles    
      write(iunit,92) 'c_1^{(0)}',
     & res(ipp,3,c1_2,imt0),res(imp,3,c1_2,imt0),
     & 'c_1^{(2)}',
     & res(ipp,3,c1_2,imt2),res(imp,3,c1_2,imt2)
      write(iunit,*) '\hline' 

      write(iunit,92) 'c_2^{(0)}',
     & res(ipp,3,c12_34,imt0),res(imp,3,c12_34,imt0),
     & 'c_2^{(2)}',
     & res(ipp,3,c12_34,imt2),res(imp,3,c12_34,imt2)
      write(iunit,*) '\hline' 

      write(iunit,92) 'c_3^{(0)}',
     & res(ipp,3,c1_34,imt0),res(imp,3,c1_34,imt0),
     & 'c_3^{(2)}',
     & res(ipp,3,c1_34,imt2),res(imp,3,c1_34,imt2)
      write(iunit,*) '\hline' 

      write(iunit,92) 'c_4^{(0)}',
     & res(ipp,3,c2_34,imt0),res(imp,3,c2_34,imt0),
     & 'c_4^{(2)}',
     & res(ipp,3,c2_34,imt2),res(imp,3,c2_34,imt2)
      write(iunit,*) '\hline' 

      write(iunit,92) 'c_5^{(0)}',
     & res(ipp,3,c1_56,imt0),res(imp,3,c1_56,imt0),
     & 'c_5^{(2)}',
     & res(ipp,3,c1_56,imt2),res(imp,3,c1_56,imt2)
      write(iunit,*) '\hline' 

      write(iunit,92) 'c_6^{(0)}',
     & res(ipp,3,c2_56,imt0),res(imp,3,c2_56,imt0),
     & 'c_6^{(2)}',
     & res(ipp,3,c2_56,imt2),res(imp,3,c2_56,imt2)
      write(iunit,*) '\hline' 

      write(iunit,*) '\hline' 

c--- bubbles    
      write(iunit,90) 'b_1',
     & res(ipp,2,b12,imt0),res(imp,2,b12,imt0)
      write(iunit,*) '\hline' 
      write(iunit,90) 'b_2',
     & res(ipp,2,b34,imt0),res(imp,2,b34,imt0)
      write(iunit,*) '\hline' 
      write(iunit,90) 'b_3',
     & res(ipp,2,b56,imt0),res(imp,2,b56,imt0)
      write(iunit,*) '\hline' 
      write(iunit,90) 'b_4',
     & res(ipp,2,b134,imt0),res(imp,2,b134,imt0)
      write(iunit,*) '\hline' 
      write(iunit,90) 'b_5',
     & res(ipp,2,b234,imt0),res(imp,2,b234,imt0)
      write(iunit,*) '\hline' 

      write(iunit,*) '\hline' 

c--- rational    
      write(iunit,90) 'R',
     & res(ipp,2,rat,imt0),res(imp,2,rat,imt0)
      write(iunit,*) '\hline' 
   
   90 format(a20,' & ',SP,2e20.10,'i','&',2e20.10,'i \\')
   92 format(a20,' & ',SP,2e20.10,'i','&',2e20.10,'i \\',/,
     &       a20,' & ',SP,2e20.10,'i','&',2e20.10,'i \\')
   94 format(a20,' & ',SP,2e20.10,'i','&',2e20.10,'i \\',/,
     &       a20,' & ',SP,2e20.10,'i','&',2e20.10,'i \\',/,
     &       a20,' & ',SP,2e20.10,'i','&',2e20.10,'i \\')
      
      write(iunit,*) '\end{array}'
      write(iunit,*) '\end{math}'

      if (iunit .ne. 6) close(iunit)

      return
      end
      
