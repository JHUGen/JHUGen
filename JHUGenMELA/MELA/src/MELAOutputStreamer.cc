#include "MELAOutputStreamer.h"


using namespace std;


MELAOutputStreamer::MELAOutputStreamer(const char* fname, std::ios_base::openmode fmode, bool printError) :
theFile(fname, fmode)
{
  if (!printError) stdout_ptr = &std::cout;
  else stdout_ptr = &std::cerr;
}
MELAOutputStreamer::~MELAOutputStreamer(){ this->close(); }

MELAOutputStreamer& MELAOutputStreamer::operator<<(std::ostream& (*fcn)(std::ostream&)){
  fcn(theFile);
  if (stdout_ptr) fcn(*stdout_ptr);
  return *this;
}
MELAOutputStreamer& MELAOutputStreamer::operator<<(std::ios& (*fcn)(std::ios&)){
  fcn(theFile);
  if (stdout_ptr) fcn(*stdout_ptr);
  return *this;
}
MELAOutputStreamer& MELAOutputStreamer::operator<<(std::ios_base& (*fcn)(std::ios_base&)){
  fcn(theFile);
  if (stdout_ptr) fcn(*stdout_ptr);
  return *this;
}
template<> MELAOutputStreamer& MELAOutputStreamer::operator<< <TLorentzVector>(TLorentzVector const& val){
  theFile << "{" << val.X() << ", " << val.Y() << ", " << val.Z() << ", " << val.T() << ", " << val.M() << "}";
  if (stdout_ptr) *stdout_ptr << "{" << val.X() << ", " << val.Y() << ", " << val.Z() << ", " << val.T() << ", " << val.M() << "}";
  return *this;
}
template<> MELAOutputStreamer& MELAOutputStreamer::operator<< <MELAParticle>(MELAParticle const& val){
  *this << "(" << val.id << ") (X,Y,Z,T,M)=( "
    << val.x() << " , "
    << val.y() << " , "
    << val.z() << " , "
    << val.t() << " , "
    << val.m()
    << " )";
  return *this;
}
template<> MELAOutputStreamer& MELAOutputStreamer::operator<< <int, TLorentzVector>(SimpleParticle_t const& val){
  *this << "(" << val.first << ") (X,Y,Z,T,M)=( "
    << val.second.X() << " , "
    << val.second.Y() << " , "
    << val.second.Z() << " , "
    << val.second.T() << " , "
    << val.second.M() << " )";
  return *this;
}
template<> MELAOutputStreamer& MELAOutputStreamer::operator<< <MELAThreeBodyDecayCandidate>(MELAThreeBodyDecayCandidate const& val){
  *this << "(" << val.id << ") (X,Y,Z,T,M)=( "
    << val.x() << " , "
    << val.y() << " , "
    << val.z() << " , "
    << val.t() << " , "
    << val.m()
    << " )"
    << endl;

  int ip=0;
  *this << "\tHas " << val.getNMothers() << " mothers" << endl;
  for (auto const& part:val.getMothers()){
    *this << "\t\tV" << ip << ' ' << *part << endl;
    ip++;
  }

  *this << "\tHas " << val.getNDaughters() << " daughters" << endl;
  *this << "\t\tPartner particle "; if (val.getPartnerParticle()) *this << *(val.getPartnerParticle()); else *this << "N/A"; *this << endl;
  *this << "\t\tW fermion "; if (val.getWFermion()) *this << *(val.getWFermion()); else *this << "N/A"; *this << endl;
  *this << "\t\tW antifermion "; if (val.getWAntifermion()) *this << *(val.getWAntifermion()); else *this << "N/A"; *this << endl;

  return *this;
}
template<> MELAOutputStreamer& MELAOutputStreamer::operator<< <MELACandidate>(MELACandidate const& val){
  int ip=0;
  *this << "\tHas " << val.getNMothers() << " mothers" << endl;
  for (auto const& part:val.getMothers()){
    *this << "\t\tV" << ip << ' ' << *part << endl;
    ip++;
  }

  *this << "\tHas " << val.getNSortedVs() << " sorted Vs" << endl;
  ip=0;
  for (auto const& part:val.getSortedVs()){
    *this << "\t\tV" << ip << ' ' << *part << endl;
    unsigned int ivd=0;
    for (auto const& dau:part->getDaughters()){
      *this << "\t\t- V" << ip << ivd << ' ' << *dau << endl;
      ivd++;
    }
    ip++;
  }

  *this << "\tHas " << val.getNAssociatedLeptons() << " leptons or neutrinos" << endl;
  ip=0;
  for (auto const& part:val.getAssociatedLeptons()){
    *this << "\t\tV" << ip << ' ' << *part << endl;
    ip++;
  }

  *this << "\tHas " << val.getNAssociatedPhotons() << " photons" << endl;
  ip=0;
  for (auto const& part:val.getAssociatedPhotons()){
    *this << "\t\tV" << ip << ' ' << *part << endl;
    ip++;
  }

  *this << "\tHas " << val.getNAssociatedJets() << " jets" << endl;
  ip=0;
  for (auto const& part:val.getAssociatedJets()){
    *this << "\t\tV" << ip << ' ' << *part << endl;
    ip++;
  }

  *this << "\tHas " << val.getNAssociatedTops() << " tops" << endl;
  ip=0;
  for (auto const& part:val.getAssociatedTops()){
    *this << "\t\tTop" << ip << ' ' << static_cast<MELAParticle const>(*part) << endl;
    { MELAParticle* bottom=part->getPartnerParticle(); if (bottom) *this << "\t\t- Top" << ip << " b " <<  ' ' << *bottom << endl; }
    { MELAParticle* Wf=part->getWFermion(); if (Wf) *this << "\t\t- Top" << ip << " Wf " << ' ' << *Wf << endl; }
    { MELAParticle* Wfb=part->getWAntifermion(); if (Wfb) *this << "\t\t- Top" << ip << " Wfb " <<  ' ' << *Wfb << endl; }
    ip++;
  }

  return *this;
}


std::streamsize MELAOutputStreamer::width() const{ return theFile.width(); }
std::streamsize MELAOutputStreamer::width(std::streamsize wide){
  if (stdout_ptr) stdout_ptr->width(wide);
  return theFile.width(wide);
}

char MELAOutputStreamer::fill() const{ return theFile.fill(); }
char MELAOutputStreamer::fill(char fillch){
  if (stdout_ptr) stdout_ptr->fill(fillch);
  return theFile.fill(fillch);
}

void MELAOutputStreamer::close(){
  theFile.close();
}
void MELAOutputStreamer::open(const char* fname, std::ios_base::openmode fmode){
  this->close();
  theFile.open(fname, fmode);
}
