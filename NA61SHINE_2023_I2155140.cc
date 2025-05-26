// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/AnalysisHandler.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/AliceCommon.hh"

namespace Rivet {

  class NA61SHINE_2023_I2155140 : public Analysis {
  public:

    //RIVET_DEFAULT_ANALYSIS_CTOR(NA61SHINE_2023_I2155140);
	NA61SHINE_2023_I2155140()
      : Analysis("NA61SHINE_2023_I2155140")
    {    }

    void init() {
    
    	// Selection of charged final states used in analyze()
		declare(ALICE::PrimaryParticles(Cuts::abscharge > 0), "APRIM");
    	
    	// Beam object used in histogram booking
		const ParticlePair& beam = beams();
    	int i = 0;
    	
    	// Selection of pi^-C @ 158 GeV/c
    	if (beam.first.pid() == PID::PIMINUS 
    	    && beam.first.pz() > 157 && beam.first.pz() < 159 
    	    && beam.second.pid() == PID::CARBON) i = 1;
    	// Selection of pi^-C @ 350 GeV/c
    	else if (beam.first.pid() == PID::PIMINUS
    	    && beam.first.pz() > 349 && beam.first.pz() < 351
    	    && beam.second.pid() == PID::CARBON) i = 2;
        
        // Histogram booking using HepData format (i=1 @158, i=2 @350): d0a-x0b-y0c
        book(_h_piP, 	i, 1, 1);
		book(_h_piM, 	i, 1, 2);
		book(_h_kP, 	i, 1, 3);
		book(_h_kM, 	i, 1, 4);
		book(_h_pP, 	i, 1, 5);
		book(_h_pM, 	i, 1, 6); 

    }

    void analyze(const Event& event) {
	
	// Working with charged particles only excluding the weak decay feed-down
	const Particles cfs = apply<ALICE::PrimaryParticles>(event, "APRIM").particles();
	for (const Particle& p : cfs) {
		
		// Definition of the weight to fulfill p(dn/dp) distributions
		double weight = p.momentum().p3().mod()/GeV;
		
		// Histogram filling of the different identitified charged final states
        	switch (p.pid()) {	
				case 211: { // positive pion s
               		_h_piP->fill(p.momentum().p3().mod()/GeV, weight);
               		break;
	    		}
	    		case -211: { // negative pions
					_h_piM->fill(p.momentum().p3().mod()/GeV, weight);
               		break;
	    		}
	    		case 321: { // positive kaons
               		_h_kP->fill(p.momentum().p3().mod()/GeV, weight);
               		break;
	    		}
	    		case -321: { // negative kaons
					_h_kM->fill(p.momentum().p3().mod()/GeV, weight);
               		break;
	    		}
	    		case 2212: { // protons
               			_h_pP->fill(p.momentum().p3().mod()/GeV, weight);
               			break;
	    		}
	    		case -2212: { // anti-protons
					_h_pM->fill(p.momentum().p3().mod()/GeV, weight);
               			break;
	    		}
				// TODO add contributions from K-Shorts, lambdas and anti-lambdas
	  	}
 	}

    }

    void finalize() {
	
	// Normalisation factor
	double scale_factor = 1./sumOfWeights();
	
	// Histogram normalisation
	scale(_h_piP,  	scale_factor);
	scale(_h_piM,  	scale_factor);
	scale(_h_kP,  	scale_factor);
	scale(_h_kM,  	scale_factor);
	scale(_h_pP,  	scale_factor);
	scale(_h_pM,  	scale_factor);

   }
   
   private:

    Histo1DPtr _h_piP, _h_piM, _h_kP, _h_kM, _h_pP, _h_pM;

  };

  RIVET_DECLARE_PLUGIN(NA61SHINE_2023_I2155140);

}
