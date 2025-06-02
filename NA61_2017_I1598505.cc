// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Event.hh"
#include "Rivet/Projections/FinalState.hh"
//#include "Rivet/Analyses/AliceCommon.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class NA61_2017_I1598505 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(NA61_2017_I1598505);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      declare(Beam(),"Beam");
      const ParticlePair& beam = beams();
      declare(FinalState(Cuts::rap >= 0 && Cuts::rap <= 3), "FS");
      //declare(ALICE::PrimaryParticles(Cuts::rap >= 0 && Cuts::rap <= 3), "FS");

      book(_c_inel, "n_inelastic");

      if (beam.first.pid() != PID::PROTON && beam.second.pid() != PID::PROTON) {
        MSG_ERROR("Beam error: Not compatible! (projectile id, target id)" );
        return;
      }

      _plab_map = {{158., 0}, {80., 1}, {40., 2}, {30.9, 3}, {20., 4}};
      int i = -1;
      for (const auto& [plab, index] : _plab_map) {
        if (std::fabs(beam.first.momentum().p3().mod() - plab*GeV) < tolerance) {
          std::cout << "beam.first.momentum().p3().mod() = " << beam.first.momentum().p3().mod() << std::endl;
          i = index;
          break;
        }
      }
      if (i >= 0) {
        book(_h_dndy_km,   59+i*12, 1, 1);
        book(_h_dndy_kp,   61+i*12, 1, 1);
        book(_h_dndy_pim,  63+i*12, 1, 1);
        book(_h_dndy_pip,  65+i*12, 1, 1);
        book(_h_dndy_p,    67+i*12, 1, 1);
        book(_h_dndy_pbar, 69+i*12, 1, 1);
      }
      else {
        MSG_ERROR("Beam error: Not compatible! (plab)");
        return;
      }

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      //const ALICE::PrimaryParticles& fs =  apply<ALICE::PrimaryParticles>(event,"FS");
      const FinalState& fs = apply<FinalState>(event, "FS");
      for (const Particle& p : fs.particles(Cuts::abspid == 321 || Cuts::abspid == 211 || Cuts::abspid == 2212 )) {

        const size_t n_particles = fs.particles().size();
        if (n_particles <= 2) {
          MSG_DEBUG("Elastic event found!");
          vetoEvent;
        }
        _c_inel -> fill();

        if(p.pid() == PID::KMINUS){
          _h_dndy_km->fill(p.rapidity());
        }
        else if (p.pid() == PID::KPLUS){
          _h_dndy_kp->fill(p.rapidity());
        }
        else if (p.pid() == PID::PIMINUS){
          _h_dndy_pim->fill(p.rapidity());
        }
        else if (p.pid() == PID::PIPLUS){
          _h_dndy_pip->fill(p.rapidity());
        }
        else if (p.pid() == PID::PROTON){
          _h_dndy_p->fill(p.rapidity());
        }
        else if (p.pid() == PID::ANTIPROTON){
          _h_dndy_pbar->fill(p.rapidity());
        }
      
        /*switch (p.pid()){
          case PID::KMINUS:
            _h_dndy_km->fill(p.rapidity());
            break;
          case PID::KPLUS:
            _h_dndy_kp->fill(p.rapidity());
            break;
          case PID::PIMINUS:
            _h_dndy_pim->fill(p.rapidity());
            break;
          case PID::PIPLUS:
            _h_dndy_pip->fill(p.rapidity());
            break;
          case PID::PROTON:
            _h_dndy_p->fill(p.rapidity());
            break;
          case PID::ANTIPROTON:
            _h_dndy_pbar->fill(p.rapidity());
            break;
        }*/
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      
      /* red line
      scale(_h_dndy_km, 1./ *_c_inel);
      for (auto& by : _h_dndy_km->bins()) by.scaleW(1./by.xWidth());
      scale(_h_dndy_kp, 1./ *_c_inel);
      for (auto& by : _h_dndy_kp->bins()) by.scaleW(1./by.xWidth());
      
      scale(_h_dndy_pim, 1./ *_c_inel);
      for (auto& by : _h_dndy_pim->bins()) by.scaleW(1./by.xWidth());
      scale(_h_dndy_pip, 1./ *_c_inel);
      for (auto& by : _h_dndy_pip->bins()) by.scaleW(1./by.xWidth());
      scale(_h_dndy_p, 1./ *_c_inel);

      for (auto& by : _h_dndy_p->bins()) by.scaleW(1./by.xWidth());
      scale(_h_dndy_pbar, 1./ *_c_inel);
      for (auto& by : _h_dndy_pbar->bins()) by.scaleW(1./by.xWidth());*/
      
      /* blue line      
      normalize(_h_dndy_km);
      normalize(_h_dndy_kp);
      normalize(_h_dndy_pim);
      normalize(_h_dndy_pip);
      normalize(_h_dndy_p);
      normalize(_h_dndy_pbar);*/

      /* green line (or orange line with ALICE::PrimaryParticles) */
      const double scale_factor = 1.0 / sumOfWeights();
      scale(_h_dndy_km,   scale_factor);
      scale(_h_dndy_kp,   scale_factor);
      scale(_h_dndy_pim,  scale_factor);
      scale(_h_dndy_pip,  scale_factor);
      scale(_h_dndy_p,    scale_factor);
      scale(_h_dndy_pbar, scale_factor);

    }

    /// @}


    /// @name Histograms
    /// @{
    std::vector<std::pair<double, int>> _plab_map;
    Histo1DPtr _h_dndy_km, _h_dndy_kp, _h_dndy_pim, _h_dndy_pip, _h_dndy_p, _h_dndy_pbar;
    double tolerance = 0.5*GeV;
    CounterPtr _c_inel;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(NA61_2017_I1598505);

}
