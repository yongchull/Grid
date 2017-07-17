#ifndef QEDFVOL_WILSONLOOPS_H
#define QEDFVOL_WILSONLOOPS_H

#include <Global.hpp>

BEGIN_QEDFVOL_NAMESPACE

template <class Gimpl> class NewWilsonLoops : public Gimpl {
public:
  INHERIT_GIMPL_TYPES(Gimpl);

  typedef typename Gimpl::GaugeLinkField GaugeMat;
  typedef typename Gimpl::GaugeField GaugeLorentz;

  //////////////////////////////////////////////////
  // directed plaquette oriented in mu,nu plane
  //////////////////////////////////////////////////
  static void dirPlaquette(GaugeMat &plaq, const std::vector<GaugeMat> &U,
                           const int mu, const int nu) {
    // Annoyingly, must use either scope resolution to find dependent base
    // class,
    // or this-> ; there is no "this" in a static method. This forces explicit
    // Gimpl scope
    // resolution throughout the usage in this file, and rather defeats the
    // purpose of deriving
    // from Gimpl.
    plaq = Gimpl::CovShiftBackward(
        U[mu], mu, Gimpl::CovShiftBackward(
                       U[nu], nu, Gimpl::CovShiftForward(U[mu], mu, U[nu])));
  }
  //////////////////////////////////////////////////
  // trace of directed plaquette oriented in mu,nu plane
  //////////////////////////////////////////////////
  static void traceDirPlaquette(LatticeComplex &plaq,
                                const std::vector<GaugeMat> &U, const int mu,
                                const int nu) {
    GaugeMat sp(U[0]._grid);
    dirPlaquette(sp, U, mu, nu);
    plaq = trace(sp);
  }
  //////////////////////////////////////////////////
  // sum over all planes of plaquette
  //////////////////////////////////////////////////
  static void sitePlaquette(LatticeComplex &Plaq,
                            const std::vector<GaugeMat> &U) {
    LatticeComplex sitePlaq(U[0]._grid);
    Plaq = zero;
    for (int mu = 1; mu < U[0]._grid->_ndimension; mu++) {
      for (int nu = 0; nu < mu; nu++) {
        traceDirPlaquette(sitePlaq, U, mu, nu);
        Plaq = Plaq + sitePlaq;
      }
    }
  }
  //////////////////////////////////////////////////
  // sum over all x,y,z,t and over all planes of plaquette
  //////////////////////////////////////////////////
  static Real sumPlaquette(const GaugeLorentz &Umu) {
    std::vector<GaugeMat> U(4, Umu._grid);

    for (int mu = 0; mu < Umu._grid->_ndimension; mu++) {
      U[mu] = PeekIndex<LorentzIndex>(Umu, mu);
    }

    LatticeComplex Plaq(Umu._grid);

    sitePlaquette(Plaq, U);

    TComplex Tp = sum(Plaq);
    Complex p = TensorRemove(Tp);
    return p.real();
  }
  //////////////////////////////////////////////////
  // average over all x,y,z,t and over all planes of plaquette
  //////////////////////////////////////////////////
  static Real avgPlaquette(const GaugeLorentz &Umu) {
    int ndim = Umu._grid->_ndimension;
    Real sumplaq = sumPlaquette(Umu);
    Real vol = Umu._grid->gSites();
    Real faces = (1.0 * ndim * (ndim - 1)) / 2.0;
    return sumplaq / vol / faces / Nc; // Nc dependent... FIXME
  }

  //////////////////////////////////////////////////
  // Wilson loop of size (R1, R2), oriented in mu,nu plane
  //////////////////////////////////////////////////
  static void wilsonLoop(GaugeMat &wl, const std::vector<GaugeMat> &U,
                           const int Rmu, const int Rnu,
                           const int mu, const int nu) {
    wl = U[nu];

    for(int i = 0; i < Rnu-1; i++){
      wl = Gimpl::CovShiftForward(U[nu], nu, wl);
    }

    for(int i = 0; i < Rmu; i++){
      wl = Gimpl::CovShiftForward(U[mu], mu, wl);
    }

    for(int i = 0; i < Rnu; i++){
      wl = Gimpl::CovShiftBackward(U[nu], nu, wl);
    }

    for(int i = 0; i < Rmu; i++){
      wl = Gimpl::CovShiftBackward(U[mu], mu, wl);
    }
  }
  //////////////////////////////////////////////////
  // trace of Wilson Loop oriented in mu,nu plane
  //////////////////////////////////////////////////
  static void traceWilsonLoop(LatticeComplex &wl,
                                const std::vector<GaugeMat> &U,
                                const int Rmu, const int Rnu,
                                const int mu, const int nu) {
    GaugeMat sp(U[0]._grid);
    wilsonLoop(sp, U, Rmu, Rnu, mu, nu);
    wl = trace(sp);
  }
  //////////////////////////////////////////////////
  // sum over all planes of Wilson loop
  //////////////////////////////////////////////////
  static void siteWilsonLoop(LatticeComplex &Wl,
                            const std::vector<GaugeMat> &U,
                            const int R1, const int R2) {
    LatticeComplex siteWl(U[0]._grid);
    Wl = zero;
    for (int mu = 1; mu < U[0]._grid->_ndimension; mu++) {
      for (int nu = 0; nu < mu; nu++) {
        traceWilsonLoop(siteWl, U, R1, R2, mu, nu);
        Wl = Wl + siteWl;
        traceWilsonLoop(siteWl, U, R2, R1, mu, nu);
        Wl = Wl + siteWl;
      }
    }
  }
  //////////////////////////////////////////////////
  // sum over planes of Wilson loop with length R1
  // in the time direction
  //////////////////////////////////////////////////
  static void siteTimelikeWilsonLoop(LatticeComplex &Wl,
                            const std::vector<GaugeMat> &U,
                            const int R1, const int R2) {
    LatticeComplex siteWl(U[0]._grid);

    int ndim = U[0]._grid->_ndimension;

    Wl = zero;
    for (int nu = 0; nu < ndim - 1; nu++) {
      traceWilsonLoop(siteWl, U, R1, R2, ndim-1, nu);
      Wl = Wl + siteWl;
    }
  }
  //////////////////////////////////////////////////
  // sum Wilson loop over all planes orthogonal to the time direction
  //////////////////////////////////////////////////
  static void siteSpatialWilsonLoop(LatticeComplex &Wl,
                            const std::vector<GaugeMat> &U,
                            const int R1, const int R2) {
    LatticeComplex siteWl(U[0]._grid);

    Wl = zero;
    for (int mu = 1; mu < U[0]._grid->_ndimension - 1; mu++) {
      for (int nu = 0; nu < mu; nu++) {
        traceWilsonLoop(siteWl, U, R1, R2, mu, nu);
        Wl = Wl + siteWl;
        traceWilsonLoop(siteWl, U, R2, R1, mu, nu);
        Wl = Wl + siteWl;
      }
    }
  }
  //////////////////////////////////////////////////
  // sum over all x,y,z,t and over all planes of Wilson loop
  //////////////////////////////////////////////////
  static Real sumWilsonLoop(const GaugeLorentz &Umu,
                            const int R1, const int R2) {
    std::vector<GaugeMat> U(4, Umu._grid);

    for (int mu = 0; mu < Umu._grid->_ndimension; mu++) {
      U[mu] = PeekIndex<LorentzIndex>(Umu, mu);
    }

    LatticeComplex Wl(Umu._grid);

    siteWilsonLoop(Wl, U, R1, R2);

    TComplex Tp = sum(Wl);
    Complex p = TensorRemove(Tp);
    return p.real();
  }
  //////////////////////////////////////////////////
  // sum over all x,y,z,t and over all planes of timelike Wilson loop
  //////////////////////////////////////////////////
  static Real sumTimelikeWilsonLoop(const GaugeLorentz &Umu,
                            const int R1, const int R2) {
    std::vector<GaugeMat> U(4, Umu._grid);

    for (int mu = 0; mu < Umu._grid->_ndimension; mu++) {
      U[mu] = PeekIndex<LorentzIndex>(Umu, mu);
    }

    LatticeComplex Wl(Umu._grid);

    siteTimelikeWilsonLoop(Wl, U, R1, R2);

    TComplex Tp = sum(Wl);
    Complex p = TensorRemove(Tp);
    return p.real();
  }
  //////////////////////////////////////////////////
  // sum over all x,y,z,t and over all planes of spatial Wilson loop
  //////////////////////////////////////////////////
  static Real sumSpatialWilsonLoop(const GaugeLorentz &Umu,
                            const int R1, const int R2) {
    std::vector<GaugeMat> U(4, Umu._grid);

    for (int mu = 0; mu < Umu._grid->_ndimension; mu++) {
      U[mu] = PeekIndex<LorentzIndex>(Umu, mu);
    }

    LatticeComplex Wl(Umu._grid);

    siteSpatialWilsonLoop(Wl, U, R1, R2);

    TComplex Tp = sum(Wl);
    Complex p = TensorRemove(Tp);
    return p.real();
  }
  //////////////////////////////////////////////////
  // average over all x,y,z,t and over all planes of Wilson loop
  //////////////////////////////////////////////////
  static Real avgWilsonLoop(const GaugeLorentz &Umu,
                            const int R1, const int R2) {
    int ndim = Umu._grid->_ndimension;
    Real sumWl = sumWilsonLoop(Umu, R1, R2);
    Real vol = Umu._grid->gSites();
    Real faces = 1.0 * ndim * (ndim - 1);
    return sumWl / vol / faces / Nc; // Nc dependent... FIXME
  }
  //////////////////////////////////////////////////
  // average over all x,y,z,t and over all planes of timelike Wilson loop
  //////////////////////////////////////////////////
  static Real avgTimelikeWilsonLoop(const GaugeLorentz &Umu,
                            const int R1, const int R2) {
    int ndim = Umu._grid->_ndimension;
    Real sumWl = sumTimelikeWilsonLoop(Umu, R1, R2);
    Real vol = Umu._grid->gSites();
    Real faces = 1.0 * (ndim - 1);
    return sumWl / vol / faces / Nc; // Nc dependent... FIXME
  }
  //////////////////////////////////////////////////
  // average over all x,y,z,t and over all planes of spatial Wilson loop
  //////////////////////////////////////////////////
  static Real avgSpatialWilsonLoop(const GaugeLorentz &Umu,
                            const int R1, const int R2) {
    int ndim = Umu._grid->_ndimension;
    Real sumWl = sumSpatialWilsonLoop(Umu, R1, R2);
    Real vol = Umu._grid->gSites();
    Real faces = 1.0 * (ndim - 1) * (ndim - 2);
    return sumWl / vol / faces / Nc; // Nc dependent... FIXME
  }
};

END_QEDFVOL_NAMESPACE

#endif // QEDFVOL_WILSONLOOPS_H