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
    bool operator == (const triplet<T>& other)const{ bool res = true; for (std::size_t i=0; i<3; i++) res &= (*this)[i]==other[i]; return res; }
    bool operator != (const triplet<T>& other)const{ return !(*this==other);  }
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
    bool operator == (const quadruplet<T>& other)const{ bool res = true; for (std::size_t i=0; i<4; i++) res &= (*this)[i]==other[i]; return res; }
    bool operator != (const quadruplet<T>& other)const{ return !(*this==other); }
  };

  typedef triplet<unsigned int> uintTriplet_t;
  typedef triplet<int> intTriplet_t;
  typedef triplet<float> floatTriplet_t;
  typedef triplet<double> doubleTriplet_t;

  typedef quadruplet<unsigned int> uintQuad_t;
  typedef quadruplet<int> intQuad_t;
  typedef quadruplet<float> floatQuad_t;
  typedef quadruplet<double> doubleQuad_t;

  void PermutationGenerator(int n, int k, std::vector<std::vector<int>>& perms, int valbegin=1, int increment=1); // n!/(n-k)! different permutations of {i_1,...,i_k}, i+j in {valbegin,...,valbegin+increment*(n-1)}
  void CombinationGenerator(int n, int k, std::vector<std::vector<int>>& perms, int valbegin=1, int increment=1); // n!/(k!(n-k)!) different combinations (unordered permutations) of {i_1,...,i_k}, i+j in {valbegin,...,valbegin+increment*(n-1)}

}

#endif

