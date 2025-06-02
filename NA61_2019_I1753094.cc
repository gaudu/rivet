// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Event.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class NA61_2019_I1753094 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(NA61_2019_I1753094);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      
      declare(FinalState(), "FS");
      declare(Beam(),"Beam");
      const ParticlePair& beam = beams();

      if (beam.first.pid() == PID::PROTON && beam.second.pid() == 1000060120) {
        book(_h_sig_prod, 1, 1, 1);
        book(_h_sig_inel, 2, 1, 1);
      }    
      else {
        MSG_ERROR("Beam error: Not compatible!");
        return;
      }       

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      if (_plab_edges.empty()) {
        _plab_edges = _h_sig_prod->xEdges(); 
        _plab_edges = _h_sig_inel->xEdges();
      }

      const ParticlePair& beams = apply<Beam>(event, "Beam").beams();
      const FinalState& fs = apply<FinalState>(event, "FS");
            
      if (isZero(beams.second.momentum().p3().mod())) { 
        const double plab = beams.first.momentum().p3().mod();

        bool plab_found = false;
        for (size_t i = 0; i < _plab_edges.size(); ++i) {
          if (std::fabs(plab - _plab_edges[i]) < tolerance) {
            
            /* veto elastic events -> inelastic cross-section */

            for (const Particle& p : fs.particles(Cuts::abspid == 1000060120)) {
              vetoEvent; // veto elastic event
            }
            
            int n_charged = 0;
            for (const Particle& p : fs.particles()) {
                if (p.charge() != 0 && p.pT() > 0.1*GeV) ++n_charged;
            }
            if (n_charged < 2) {
              vetoEvent; // veto elastic event
            }
            _h_sig_inel->fill(_plab_edges[i]); 

            /* tag event producing non-nucleon particles -> production cross-section */

            bool isProduction = false;
            for (const Particle& p : fs.particles()) {
              const int pdgid = abs(p.pid());
              // ignore nucleons (proton, neutron) and leptons (e, mu, nu)
              if (pdgid == 2212 || pdgid == 2112) continue;
              if (pdgid == 11 || pdgid == 13 || (pdgid >= 12 && pdgid <= 16)) continue; 
              if (p.isHadron()) {
                isProduction = true;
                break;
              }
            }
            if (isProduction) {
              _h_sig_prod->fill(_plab_edges[i]);
            }

            plab_found = true;
            break;
          }
        }

        if (!plab_found) {
          MSG_WARNING("Warning: Beam plab = " << plab << " GeV/c does not match any reference plab.");
        }

      } else {
        std::cerr << "Error: Not using fixed-target mode." << std::endl;
        vetoEvent;
      }      
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      scale(_h_sig_prod, crossSection()/millibarn/sumOfWeights());
      scale(_h_sig_inel, crossSection()/millibarn/sumOfWeights());

    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<int> _h_sig_prod;
    BinnedHistoPtr<int> _h_sig_inel;
    vector<int> _plab_edges;
    double tolerance = 0.5;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(NA61_2019_I1753094);

}
