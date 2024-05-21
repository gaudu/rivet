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
      const ParticlePair& beam = beams();
      
      if (beam.first.pid() == PID::PIPLUS && beam.second.pid() == PID::PROTON) {
        btype = 1;
      }
      else if (beam.first.pid() == PID::PIMINUS && beam.second.pid() == PID::PROTON) {
        btype = 2;
      }
      
      book(_h["AAAA"], 1, 1, 1);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Fill histogram 
      _h["XXXX"]->fill();

    }


    /// Normalise histograms etc., after the run
    void finalize() {

      normalize(_h["XXXX"]); // normalize to unity
      normalize(_h["YYYY"], crossSection()/microbarn); // normalize to generated cross-section in pb (no cuts)
      scale(_h["ZZZZ"], crossSection()/microbarn/sumW()); // norm to generated cross-section in pb (after cuts)

    }

    /// @}


    /// @name Histograms
    /// @{
    map<string, Histo1DPtr> _h;
    map<string, Profile1DPtr> _p;
    map<string, CounterPtr> _c;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(E104_1976_I98502);

}
