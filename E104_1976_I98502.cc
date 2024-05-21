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

      const ParticlePair& beam = beams();
      
      std::unordered_map<std::pair<PID, PID>, std::function<void()>> book_histo_map = {
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
      if (book_histo_map.find(key) != book_histo_map.end()) {
        book_histo_map[key]();
      } else {
        MSG_WARNING("Beam type not compatible with this analysis.");
      }
      
      /* total cross-section
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

      /* total cross-section ratios
      book(_h["ap_p-p_p"],    17, 1, 1);
      book(_h["ap_d-p_d"],    18, 1, 1);
      book(_h["ap_n-p_n"],    19, 1, 1);
      book(_h["km_p-kp_p"],   20, 1, 1);
      book(_h["km_d-kp_d"],   21, 1, 1);
      book(_h["km_n-kp_n"],   22, 1, 1);
      book(_h["pim_p-pip_p"], 23, 1, 1);
      book(_h["pim_d-pip_d"], 24, 1, 1);
      */
    }

    void analyze(const Event& event) {

      const ChargedFinalState& cfs = ChargedFinalState();
      const ParticlePair& beam = beams();

      for (const Particle& p : cfs.particles(Cuts::eta > 0)) {
      
        std::unordered_map<std::pair<PID, PID>, std::function<void()>> fill_histo_map = {
          {{PID::PROTON,      PID::PROTON},   []() {_h["p_p"]->fill();}},
          {{PID::PROTON,      PID::DEUTERON}, []() {_h["p_d"]->fill();}},
          {{PID::PROTON,      PID::NEUTRON},  []() {_h["p_n"]->fill();}},
          {{PID::ANTIPROTON,  PID::PROTON},   []() {_h["ap_p"]->fill();}},
          {{PID::ANTIPROTON,  PID::DEUTERON}, []() {_h["ap_d"]->fill();}},
          {{PID::ANTIPROTON,  PID::NEUTRON},  []() {_h["ap_n"]->fill();}},
          {{PID::KPLUS,       PID::PROTON},   []() {_h["kp_p"]->fill();}},
          {{PID::KPLUS,       PID::DEUTERON}, []() {_h["kp_d"]->fill();}},
          {{PID::KPLUS,       PID::NEUTRON},  []() {_h["kp_n"]->fill();}},
          {{PID::KMINUS,      PID::PROTON},   []() {_h["km_p"]->fill();}},
          {{PID::KMINUS,      PID::DEUTERON}, []() {_h["km_d"]->fill();}},
          {{PID::KMINUS,      PID::NEUTRON},  []() {_h["km_n"]->fill();}},
          {{PID::PIPLUS,      PID::PROTON},   []() {_h["pip_p"]->fill();}},
          {{PID::PIPLUS,      PID::DEUTERON}, []() {_h["pip_d"]->fill();}},
          {{PID::PIMINUS,     PID::PROTON},   []() {_h["pim_p"]->fill();}},
          {{PID::PIMINUS,     PID::DEUTERON}, []() {_h["pim_d"]->fill();}},
        };

        std::pair<PID, PID> key = {beam.first.pid(), beam.second.pid()};
        if (fill_histo_map.find(key) != fill_histo_map.end()) {
          fill_histo_map[key]();
        }

      }
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
