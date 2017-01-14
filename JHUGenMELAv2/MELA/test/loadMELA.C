{
  TString LIBMCFMPATH = "../data/slc6_amd64_gcc530/";
  TString LIBMCFM="libmcfm_703.so";
  TString LIBJHUGENMELA="libjhugenmela.so";
  TString LIBMELA="libMELA.so";

  gSystem->AddIncludePath("-I$ROOFITSYS/include/");
  gSystem->AddIncludePath("-I../interface/");
  gSystem->Load(LIBMCFMPATH + LIBMCFM);
  gSystem->Load(LIBMCFMPATH + LIBJHUGENMELA);
  gSystem->Load(LIBMCFMPATH + LIBMELA);
}
