#ifndef MCFM_UTILITIES
#define MCFM_UTILITIES
#include <string>
#include <vector>
// PDGHelpers
#include "PDGHelpers.h"
// ROOT includes
#include "TLorentzVector.h"


namespace TMCFMUtils{

  template<typename T> struct quad{
    T value[4];
    quad(T i1, T i2, T i3, T i4){
      value[0]=i1;
      value[1]=i2;
      value[2]=i3;
      value[3]=i4;
    }
    quad(T i1){
      value[0]=i1;
      value[1]=i1;
      value[2]=i1;
      value[3]=i1;
    }
    quad(){}
    T& operator[](std::size_t ipos){ return value[ipos]; } // Return by reference
    const T& operator[](std::size_t ipos)const{ return value[ipos]; } // Return by value
  };
  typedef quad<int> intQuad_t;


  void AssociatedParticleOrdering_QQVVQQAny(int iSel, int jSel, int rSel, int sSel, int order[2]);
  std::vector<intQuad_t> Hash_QQVVQQAny();
  std::vector<intQuad_t> Hash_QQVVQQ();
  std::vector<intQuad_t> Hash_QQVVQQStrong();
}

#endif

