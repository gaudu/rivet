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

      const ChargedFinalState& cfs = ChargedFinalState();
      declare(ChargedFinalState(), "CFS");

      // beam.first.pid() can be a proton, antiproton, K+, K-, pi+ or pi- for d01-x01-y01 to d16-x01-y01
      // beam.second.pid() can be a proton, deuteron or neutron for d01-x01-y01 to d16-x01-y01
      // d17-x01-y01 to d24-x01-y01 are ratios

      const ParticlePair& beam = beams();
      
      std::unordered_map<std::pair<PID, PID>, std::function<void()>> booking_beam_pid_map = {
        {{PID::PROTON,      PID::PROTON},   []() {book(_h["p_p"],    1,  1, 1);}},
        {{PID::PROTON,      PID::DEUTERON}, []() {book(_h["p_d"],    2,  1, 1);}},
        {{PID::PROTON,      PID::NEUTRON},  []() {book(_h["p_n"],    3,  1, 1);}},
        {{PID::ANTIPROTON,  PID::PROTON},   []() {book(_h["ap_p"],   4,  1, 1);}},
        {{PID::ANTIPROTON,  PID::DEUTERON}, []() {book(_h["ap_d"],   5,  1, 1);}},
        {{PID::ANTIPROTON,  PID::NEUTRON},  []() {book(_h["ap_n"],   6,  1, 1);}},
        {{PID::KPLUS,       PID::PROTON},   []() {book(_h["kp_p"],   7,  1, 1);}},
        {{PID::KPLUS,       PID::DEUTERON}, []() {book(_h["kp_d"],   8,  1, 1);}},
        {{PID::KPLUS,       PID::NEUTRON},  []() {book(_h["kp_n"],   9,  1, 1);}},
        {{PID::KMINUS,      PID::PROTON},   []() {book(_h["km_p"],  10,  1, 1);}},
        {{PID::KMINUS,      PID::DEUTERON}, []() {book(_h["km_d"],  11,  1, 1);}},
        {{PID::KMINUS,      PID::NEUTRON},  []() {book(_h["km_n"],  12,  1, 1);}},
        {{PID::PIPLUS,      PID::PROTON},   []() {book(_h["pip_p"], 13,  1, 1);}},
        {{PID::PIPLUS,      PID::DEUTERON}, []() {book(_h["pip_d"], 14,  1, 1);}},
        {{PID::PIMINUS,     PID::PROTON},   []() {book(_h["pim_p"], 15,  1, 1);}},
        {{PID::PIMINUS,     PID::DEUTERON}, []() {book(_h["pim_d"], 16,  1, 1);}},
      };

      std::pair<PID, PID> key = {beam.first.pid(), beam.second.pid()};
      if (booking_beam_pid_map.find(key) != booking_beam_pid_map.end()) {
        booking_beam_pid_map[key]();
      } else {
        MSG_WARNING("Beam type not compatible with this analysis.");
      }
      
      /*
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
      */

      book(_h["ap_p-p_p"],    17, 1, 1);
      book(_h["ap_d-p_d"],    18, 1, 1);
      book(_h["ap_n-p_n"],    19, 1, 1);
      book(_h["km_p-kp_p"],   20, 1, 1);
      book(_h["km_d-kp_d"],   21, 1, 1);
      book(_h["km_n-kp_n"],   22, 1, 1);
      book(_h["pim_p-pip_p"], 23, 1, 1);
      book(_h["pim_d-pip_d"], 24, 1, 1);

    }

    void analyze(const Event& event) {

      const ParticlePair& beam = beams();

      // beam.first.pid(): PID::PROTON, PID::ANTIPROTON, PID::KPLUS, PID::KMINUS, PID::PIPLUS or PID::PIMINUS
      // beam.second.pid(): PID::PROTON, PID::DEUTERON or PID::NEUTRON
      
      //switch/case structure
      switch ( beam.second.pid() ) {

      case PID::PROTON:
        switch ( beam.first.pid() ) {
        case PID::PROTON:
          _h["p_p"]->fill();
          break;
        case PID::ANTIPROTON:
          _h["ap_p"]->fill();
          break;
        }
      break;
    
      case PID::DEUTERON:
      switch ( beam.first.pid() ) {
        case PID::PROTON:
          _h["p_d"]->fill();
          break;
        case PID::ANTIPROTON:
          _h["ap_d"]->fill();
          break;
        }
      break;

      case PID::NEUTRON:
      switch ( beam.first.pid() ) {
        case PID::PROTON:
          _h["p_n"]->fill();
          break;
        case PID::ANTIPROTON:
          _h["ap_n"]->fill();
          break;
        }
      break;

      }

      // if/else structure
      if (beam.first.pid() == PID::PROTON && beam.second.pid() == PID::PROTON) {
        _h["p_p"]->fill();
      }
      else if (beam.first.pid() == PID::PIMINUS && beam.second.pid() == PID::PROTON) {
        _h["pim_p"]->fill();
      }

      // unordered_map structure?

    }


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
