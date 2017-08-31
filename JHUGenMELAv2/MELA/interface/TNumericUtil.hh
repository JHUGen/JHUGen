#ifndef NUMERIC_UTILITIES
#define NUMERIC_UTILITIES

#include <vector>


namespace TNumericUtil{

  template<typename T> struct triplet{
    T value[3];
    triplet(T i1, T i2, T i3){
      value[0]=i1;
      value[1]=i2;
      value[2]=i3;
    }
    triplet(T i1){
      value[0]=i1;
      value[1]=i1;
      value[2]=i1;
    }
    triplet(){}
    T& operator[](std::size_t ipos){ return value[ipos]; } // Return by reference
    const T& operator[](std::size_t ipos)const{ return value[ipos]; } // Return by const reference
  };
  template<typename T> struct quadruplet{
    T value[4];
    quadruplet(T i1, T i2, T i3, T i4){
      value[0]=i1;
      value[1]=i2;
      value[2]=i3;
      value[3]=i4;
    }
    quadruplet(T i1){ for (unsigned int idim=0; idim<4; idim++) value[idim] = i1; }
    quadruplet(){}
    T& operator[](std::size_t ipos){ return value[ipos]; } // Return by reference
    const T& operator[](std::size_t ipos)const{ return value[ipos]; } // Return by const reference
  };

  typedef triplet<int> intTriplet_t;
  typedef triplet<float> floatTriplet_t;
  typedef triplet<double> doubleTriplet_t;

  typedef quadruplet<int> intQuad_t;
  typedef quadruplet<float> floatQuad_t;
  typedef quadruplet<double> doubleQuad_t;

}

#endif

