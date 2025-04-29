// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
<<<<<<< HEAD
#include "Rivet/Projections/ChargedFinalState.hh"
=======
#include "Rivet/Tools/BinnedHistogram.hh"
>>>>>>> 673eabb (restart plugin for E104_1976_I98502 on 4.1.0 with insights from Lund)

namespace Rivet {


  /// @brief Add a short analysis description here
  class E104_1976_I98502 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(E104_1976_I98502);


    /// @name Analysis methods
    /// @{

<<<<<<< HEAD
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
=======
    void init() {
      
      declare(Beam(), "Beam");
      const std::map<int, int> proj = {
        {2212, 0},
        {-2212, 3},
        {321, 6},
        {-321, 9},
        {211, 12},
        {-211, 14}
      };
      const std::map<int, int> targ = {
        {2212, 1},
        {1000010020, 2},
        {2112, 3}
      }; 

      auto proj_it = proj.find(beamIDs().first);
      auto targ_it = targ.find(beamIDs().second);
      if (proj_it != proj.end() && targ_it != targ.end()) {
          int i = proj_it->second + targ_it->second;
          book(_sig_tot, i, 1, 1);
      } else {
          std::cerr << "Error: Invalid beam PID. Did not properly book counter." << std::endl;
      }

      _plab_edges = {23, 35, 50, 70, 100, 120, 150, 170, 200, 240, 280};

>>>>>>> 673eabb (restart plugin for E104_1976_I98502 on 4.1.0 with insights from Lund)
    }

    void analyze(const Event& event) {

<<<<<<< HEAD
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

=======
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
            
      if (isZero(beams.second.momentum().p3().mod())) { // fixed-target mode
        const double plab = beams.first.momentum().p3().mod() / GeV; 

        bool found = false;
        for (size_t i = 0; i < plab_paper.size(); ++i) {
          if (std::fabs(plab - plab_paper[i]) < tolerance) {
            _sig_tot[i]->fill(sqrtS()/GeV);  
            found = true;
            break;
          }
        }

        if (!found) {
          MSG_WARNING("Warning: Beam plab = " << plab << " GeV/c does not match any reference plab.");
        }

      } else {
        std::cerr << "Error: Not using fixed-target mode." << std::endl;
        vetoEvent; // reject event
      }
    }


    void finalize() {

      _plab_edges = {23, 35, 50, 70, 100, 120, 150, 170, 200, 240, 280};
      for (size_t i = 0; i < _plab_edges.size(); ++i) {
        scale(_sig_tot[i], crossSection()/millibarn/sumW()); 
      }
>>>>>>> 673eabb (restart plugin for E104_1976_I98502 on 4.1.0 with insights from Lund)
    }

    /// @}


    /// @name Histograms
    /// @{
<<<<<<< HEAD
    map<string, Histo1DPtr> _h;
=======
    BinnedHistoPtr<Counter> _sig_tot[11];
    std::vector<double> _plab_edges;
    static const double tolerance = 0.5;
>>>>>>> 673eabb (restart plugin for E104_1976_I98502 on 4.1.0 with insights from Lund)
    /// @}


  };


  RIVET_DECLARE_PLUGIN(E104_1976_I98502);

}
