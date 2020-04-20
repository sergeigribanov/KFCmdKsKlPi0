#include <cmath>
#include <iostream>
#include <algorithm>
#include <limits>
#include <set>

#include <TCanvas.h>
#include <TH2.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TStopwatch.h>
#include <TMath.h>

#include "TrPh.h"

#include <KFCmd/HypoKsKlPi0.hpp>

const double TrPh::_dZ = 30;
const double TrPh::_dRho = 30;
const double TrPh::_mindEdX = 0;
const double TrPh::_maxdEdX = 15000;
const double TrPh::_minTPtot = 5;
const double TrPh::_maxTPtot = 1000;

TrPh::TrPh(TTree *tree) :
  KFCmd::TrPh(tree) {}

TrPh::~TrPh() {}

bool TrPh::cutTracks() {
  _trackIndices.clear();
  for (int i = 0; i < nt; i++) {
    bool point = (std::fabs(tz[i]) < _dZ) && (std::fabs(trho[i]) < _dRho);
    bool dedx = (tdedx[i] > _mindEdX) && (tdedx[i] < _maxdEdX);
    bool ptot = (tptot[i] > _minTPtot) && (tptot[i] < _maxTPtot);
    if (point && dedx && ptot) _trackIndices.push_back(i);
  }
  if (_trackIndices.size() == 2) {
    int totalCharge = 0;
    for (int i = 0; i < 4; ++i) totalCharge += tcharge[_trackIndices[i]];
    return (totalCharge == 0);
  }
  return false;
}

bool TrPh::cutPhotons() {
  _photonIndices.clear();
  for (int i = 0; i < nph; ++i) {
    if (phen[i] > 20)
      _photonIndices.push_back(i);
  }
  if (_photonIndices.size() > 1) return true;
  return false;
}

Int_t TrPh::Cut(Long64_t) {
  if (nt < 2) return -1;
  if (!cutTracks()) return -1;
  if (!cutPhotons()) return -1;
  // if (!cutPhotons()) return -1;
  std::vector<Int_t> charges(nt);
  std::copy(tcharge, tcharge + nt, charges.begin());
  std::sort(_trackIndices.begin(), _trackIndices.end(),
            [&charges](int i, int j) { return charges[i] < charges[j]; });
  return 1;
}

void TrPh::Loop(const std::string& outpath, double magneticField) {
  if (fChain == 0) return;
  std::set<std::string> sKs = {"pi+_1", "pi-_1"};
  std::set<std::string> sgg = {"g0", "g1"};
  TStopwatch timeLoop;
  TStopwatch timeCommon;
  auto outfl = TFile::Open(outpath.c_str(), "recreate");
  TH1F h_kf_mks("h_kf_mks", "", 500, 0, 1000);
  TH1F h_in_mks("h_in_mks", "", 500, 0, 1000);
  TH1F h_kf_mgg("h_kf_mgg", "", 500, 0, 1000);
  TH1F h_in_mgg("h_in_mgg", "", 500, 0, 1000);
  TH1F h_kf_chi2("h_kf_chi2", "", 150, 0, 300);
  TH1F h_vtx_dr("h_vtx_dr", "", 120, 0, 30);
  TH1F h_dirang("h_dirang", "", 128, 0, TMath::Pi());
  TH1F h_vtx0_z("h_vtx0_z", "", 1000, -50, 50);
  TH1F h_vtx1_z("h_vtx1_z", "", 1000, -50, 50);
  TH2F h_vtx0_xy("h_vtx0_xy", "", 1000, -10, 10, 1000, -10, 10);
  TH2F h_vtx1_xy("h_vtx1_xy", "", 1000, -10, 10, 1000, -10, 10);
  fChain->GetEntry(0);
  KFCmd::HypoKsKlPi0 hypo(2 * emeas, magneticField);
  hypo.setBeamXY(xbeam, ybeam);
  hypo.fixVertexComponent("vtx0", xbeam, KFBase::VERTEX_X);
  hypo.fixVertexComponent("vtx0", ybeam, KFBase::VERTEX_Y);
  double kf_chi2;
  double in_mks = 0;
  double kf_mks = 0;
  double in_mgg = 0;
  double kf_mgg = 0;
  double vtx_dr = 0;
  double vtx0_x = 0;
  double vtx0_y = 0;
  double vtx0_z = 0;
  double vtx1_x = 0;
  double vtx1_y = 0;
  double vtx1_z = 0;
  double dirang = 0;
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  timeCommon.Start();
  timeLoop.Start();
  for (Long64_t jentry = 0; jentry < nentries; jentry++) {
    if (jentry % 1000 == 0) {
      std::cout << jentry << std::endl;
      timeLoop.Stop();
      std::cout << "CPU TIME : " << timeLoop.CpuTime() << std::endl;
      timeLoop.Start();
    }
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);
    nbytes += nb;
    if (Cut(ientry) < 0) continue;
    if (nks != 1) continue;
    hypo.setBeamXY(xbeam, ybeam);
    hypo.fixVertexComponent("vtx0", xbeam, KFBase::VERTEX_X);
    hypo.fixVertexComponent("vtx0", ybeam, KFBase::VERTEX_Y);
    if (!hypo.fillTrack("pi-_1", _trackIndices[0], *this)) continue;
    if (!hypo.fillTrack("pi+_1", _trackIndices[1], *this)) continue;
    kf_chi2 = std::numeric_limits<double>::infinity();
    for (int iph = 0; iph + 1 < (int) _photonIndices.size(); ++iph) {
      if (!hypo.fillPhoton("g0", _photonIndices[iph], *this)) continue;
      for (int jph = iph + 1; jph < (int) _photonIndices.size(); ++jph) {
	if (!hypo.fillPhoton("g1", _photonIndices[jph], *this)) continue;
	hypo.optimize();
	if (hypo.getErrorCode() != 0) continue;
	if (hypo.getChiSquare() > kf_chi2) continue;
	kf_chi2 = hypo.getChiSquare();
	in_mks = hypo.getInitialMomentum(sKs).M();
	kf_mks = hypo.getFinalMomentum(sKs).M();
	in_mgg = hypo.getInitialMomentum(sgg).M();
	kf_mgg = hypo.getFinalMomentum(sgg).M();
	TVector3 vtx1 = hypo.getFinalVertex("vtx1");
	TVector3 vtx0 = hypo.getFinalVertex("vtx0");
	vtx_dr = (vtx1 - vtx0).Mag();
	vtx0_z = vtx0.Z();
	vtx1_z = vtx1.Z();
	vtx0_x = vtx0.X();
	vtx0_y = vtx0.Y();
	vtx0_z = vtx0.Z();
	vtx1_x = vtx1.X();
	vtx1_y = vtx1.Y();
	vtx1_z = vtx1.Z();
	TVector3 ta = vtx1 - vtx0;
	TVector3 tb = hypo.getFinalMomentum(sKs).Vect();
	dirang = TMath::ACos((ta * tb) / ta.Mag() / tb.Mag());
      }
    }

    h_kf_chi2.Fill(kf_chi2);
    h_in_mks.Fill(in_mks);
    h_kf_mks.Fill(kf_mks);
    h_in_mgg.Fill(in_mgg);
    h_kf_mgg.Fill(kf_mgg);
    h_vtx_dr.Fill(vtx_dr);
    h_dirang.Fill(dirang);
    h_vtx0_z.Fill(vtx0_z);
    h_vtx1_z.Fill(vtx1_z);
    h_vtx0_xy.Fill(vtx0_x, vtx0_y);
    h_vtx1_xy.Fill(vtx1_x, vtx1_y);
  }
  outfl->cd();
  h_kf_chi2.Write();
  h_in_mks.Write();
  h_kf_mks.Write();
  h_in_mgg.Write();
  h_kf_mgg.Write();
  h_vtx_dr.Write();
  h_dirang.Write();
  h_vtx0_z.Write();
  h_vtx1_z.Write();
  h_vtx0_xy.Write();
  h_vtx1_xy.Write();
  outfl->Close();
  delete outfl;
}
