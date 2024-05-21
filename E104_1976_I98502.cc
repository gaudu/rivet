// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class E104_1976_I98502 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(E104_1976_I98502);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      declare(ChargedFinalState(), "CFS");

      // beam.first.pid() can be a proton, antiproton, K+, K-, pi+ or pi- for d01-x01-y01 to d16-x01-y01
      // beam.second.pid() can be a proton, deuteron or neutron for d01-x01-y01 to d16-x01-y01
      // d17-x01-y01 to d24-x01-y01 are ratios
      
      book(_h["p_p"],   1,  1, 1);
      book(_h["p_d"],   2,  1, 1);
      book(_h["p_n"],   3,  1, 1);
      book(_h["ap_p"],  4,  1, 1);
      book(_h["ap_d"],  5,  1, 1);
      book(_h["ap_n"],  6,  1, 1);
      book(_h["kp_p"],  7,  1, 1);
      book(_h["kp_d"],  8,  1, 1);
      book(_h["kp_n"],  9,  1, 1);
      book(_h["km_p"],  10, 1, 1);
      book(_h["km_d"],  11, 1, 1);
      book(_h["km_n"],  12, 1, 1);
      book(_h["pip_p"], 13, 1, 1);
      book(_h["pip_d"], 14, 1, 1);
      book(_h["pim_p"], 15, 1, 1);
      book(_h["pim_d"], 16, 1, 1);

      book(_h["ap_p-p_p"],    17, 1, 1);
      book(_h["ap_d-p_d"],    18, 1, 1);
      book(_h["ap_n-p_n"],    19, 1, 1);
      book(_h["km_p-kp_p"],   20, 1, 1);
      book(_h["km_d-kp_d"],   21, 1, 1);
      book(_h["km_n-kp_n"],   22, 1, 1);
      book(_h["pim_p-pip_p"], 23, 1, 1);
      book(_h["pim_d-pip_d"], 24, 1, 1);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const ParticlePair& beam = beams();

      if (beam.first.pid() == PID::PROTON && beam.second.pid() == PID::PROTON) {
        _h["p_p"]->fill();
      }
      else if (beam.first.pid() == PID::PIMINUS && beam.second.pid() == PID::PROTON) {
        _h["pim_p"]->fill();
      }

      // Fill histogram 
      _h["XXXX"]->fill();

    }


    /// Normalise histograms etc., after the run
    void finalize() {

      normalize(_h["YYYY"], crossSection()/millibarn); // normalize to generated cross-section in mb (no cuts)
      scale(_h["ZZZZ"], crossSection()/millibarn/sumW()); // norm to generated cross-section in mb (after cuts)

    }

    /// @}


    /// @name Histograms
    /// @{
    map<string, Histo1DPtr> _h;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(E104_1976_I98502);

}