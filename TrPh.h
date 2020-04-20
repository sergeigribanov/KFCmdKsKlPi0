#ifndef __KsKlPi0_TrPh__
#define __KsKlPi0_TrPh__

#include <KFCmd/TrPh.hpp>

class TrPh : public KFCmd::TrPh {
 public:
  TrPh(TTree *tree = 0);
  virtual ~TrPh();
  virtual Int_t Cut(Long64_t entry) override final;
  virtual void Loop(const std::string&, double magneticField = 1.3) override final;
 private:
  bool cutTracks();
  bool cutPhotons();
  static const double _dZ;
  static const double _dRho;
  static const double _mindEdX;
  static const double _maxdEdX;
  static const double _minTPtot;
  static const double _maxTPtot;
  std::vector<std::size_t> _trackIndices;
  std::vector<std::size_t> _photonIndices;
};

#endif
