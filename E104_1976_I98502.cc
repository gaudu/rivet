// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class E104_1976_I98502 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(E104_1976_I98502);


    /// @name Analysis methods
    /// @{

    void init() {
      
      _proj = {{2212, 0}, {-2212, 3}, {321, 6}, {-321, 9}, {211, 12}, {-211, 14}};
      _targ = {{2212, 1}, {1000010020, 2}, {2112, 3}}; 
      _plab_edges = {23, 35, 50, 70, 100, 120, 150, 170, 200, 240, 280};

      declare(Beam(), "Beam");
      
      auto proj_it = _proj.find(beamIDs().first);
      auto targ_it = _targ.find(beamIDs().second);
      if (proj_it != _proj.end() && targ_it != _targ.end()) {
          int i = proj_it->second + targ_it->second;
          cout << "Booking d0" << i << "-x01-y01." << std::endl;
          book(_h_sig_tot, i, 1, 1);
      } else {
          std::cerr << "Error: Invalid beam PID. Did not properly book." << std::endl;
      }

    }

    void analyze(const Event& event) {

      const ParticlePair& beams = apply<Beam>(event, "Beam").beams();
            
      if (isZero(beams.second.momentum().p3().mod())) { // fixed-target mode
        const double plab = beams.first.momentum().p3().mod();

        bool found = false;
        for (size_t i = 0; i < _plab_edges.size(); ++i) {
          if (std::fabs(plab - _plab_edges[i]) < tolerance) {
            cout << "_h_sig_tot->fill(" << _plab_edges[i] << ");" << std::endl;
            _h_sig_tot->fill(_plab_edges[i]);  
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

      scale(_h_sig_tot, crossSection()/millibarn/sumOfWeights());
    
    }   

    /// @}


    /// @name Histograms
    /// @{
    std::map<int, int> _proj; 
    std::map<int, int> _targ;
    BinnedHistoPtr<int> _h_sig_tot;
    std::vector<double> _plab_edges;
    double tolerance = 0.5;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(E104_1976_I98502);

}
