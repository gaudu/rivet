// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Tools/BinnedHistogram.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class E104_1976_I98502 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(E104_1976_I98502);


    /// @name Analysis methods
    /// @{

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

    }

    void analyze(const Event& event) {

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
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<Counter> _sig_tot[11];
    std::vector<double> _plab_edges;
    static const double tolerance = 0.5;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(E104_1976_I98502);

}
