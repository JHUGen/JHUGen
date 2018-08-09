#ifndef MELAHXSWIDTH_H
#define MELAHXSWIDTH_H

#include <string>
#include "TGraph.h"
#include "TSpline.h"


class MELAHXSWidth{
public:
  MELAHXSWidth(std::string fileLoc = "../txtFiles", std::string strAppend="YR3");
  ~MELAHXSWidth();
  double HiggsWidth(double mH);
protected:
  std::string fileName;
  std::vector<double> mass_BR;
  std::vector<double> BR;
  double* xmhW;
  double* sigW;
  TGraph* graphW;
  TSpline3* gsW;
};

#endif

