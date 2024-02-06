// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Tools/BinnedHistogram.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/AnalysisHandler.hh"

namespace Rivet {


  /// Double-diferential cross-sections for prompt J/psi in pp collisions at 7 TeV 
  class LHCB_2011_I891233 : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(LHCB_2011_I891233);  	

    /// Book histograms and initialise projections before the run
    void init() {
	
	// Projection for J/psi    
    	declare(UnstableParticles(Cuts::abspid == PID::JPSI && Cuts::absrapIn(2.0, 4.5) && Cuts::ptIn(0.0*GeV, 14.0*GeV)), "JPSI");
	
	// Book histograms for different rapidity ranges
	{Histo1DPtr tmp; _h_jpsi_pT_y.add(2.0, 2.5, book(tmp, 6, 1, 1));}
	{Histo1DPtr tmp; _h_jpsi_pT_y.add(2.5, 3.0, book(tmp, 6, 1, 2));}
	{Histo1DPtr tmp; _h_jpsi_pT_y.add(3.0, 3.5, book(tmp, 6, 1, 3));}
	{Histo1DPtr tmp; _h_jpsi_pT_y.add(3.5, 4.0, book(tmp, 6, 1, 4));}
	{Histo1DPtr tmp; _h_jpsi_pT_y.add(4.0, 4.5, book(tmp, 6, 1, 5));}

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
	const double weight = 1.0;
	const UnstableParticles& jpsi = apply<UnstableParticles> (event, "JPSI");
	
	for (const Particle& p : jpsi.particles()) {
        	if (p.abspid() != 443) continue;
        	ConstGenVertexPtr gv = p.genParticle()->production_vertex();
        	bool nonPrompt = false;
        	if (gv) {
          		for (ConstGenParticlePtr pi: HepMCUtils::particles(gv, Relatives::ANCESTORS)) {
            			const PdgId pid2 = pi->pdg_id();
            			if (PID::isHadron(pid2) && PID::hasBottom(pid2)) {
              				nonPrompt = true;
              				break;
            			}	
          		}
        	}
		double pT = p.pt();
		double y = p.absrap();
		if (!nonPrompt) _h_jpsi_pT_y.fill(y, pT/GeV, weight);
 	}	
    }


    /// Normalise histogram, Scale the acceptance histogram etc., after the run
    void finalize() {
	// Avoid the implicit division by the bin width in the BinnedHistogram::scale method
	double scale_factor = crossSection()/nanobarn/2.0/sumOfWeights();//*0.017; //0.017 factor for EPOS3447 data
	for (Histo1DPtr h : _h_jpsi_pT_y.histos()) h->scaleW(scale_factor);
    }

  private:

    /// Histogram
    BinnedHistogram _h_jpsi_pT_y;

  };

  // Hook for the plugin system
  RIVET_DECLARE_PLUGIN(LHCB_2011_I891233);

}
