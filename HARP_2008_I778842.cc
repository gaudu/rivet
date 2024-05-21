// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Tools/BinnedHistogram.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {

  class HARP_2008_I778842 : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(HARP_2008_I778842);

    void init() {
    
    	/* Could this function be used instead?
    	Rivet::Analysis::isCompatible (PdgId beam1, PdgId beam2, double e1, double e2) const
    	Not if the IgnoreBeam flag is turned on! */
    
    	const ParticlePair& beam = beams();
    	int i = 0;
      	if (beam.first.pid() == PID::PROTON) { // && beam.second.pid() == PID::CARBON) if added to ParticleName.hh  	
        	i = 1;
		}
      	else if (beam.first.pid() == PID::PIPLUS) {
        	i = 8;
		}
      	else if (beam.first.pid() == PID::PIMINUS) {
        	i = 14;
        }	
      	
    	/* 	
    	p C -> pi X:   i = 1	 1 -  7
    	pi+ C -> pi X: i = 8	 8 - 13
    	pi- C -> pi X: i = 14	14 - 19
    	*/
    	
   	// -> pi+ X  		  theta ranges	        d0a-x0b-y01
     	{Histo1DPtr tmp; _h_pip.add(.03, .06, 	book(tmp, i,   1, 1));}
		{Histo1DPtr tmp; _h_pip.add(.06, .09, 	book(tmp, i+1, 1, 1));}
		{Histo1DPtr tmp; _h_pip.add(.09, .12, 	book(tmp, i+2, 1, 1));}
		{Histo1DPtr tmp; _h_pip.add(.12, .15, 	book(tmp, i+3, 1, 1));}
		{Histo1DPtr tmp; _h_pip.add(.15, .18, 	book(tmp, i+4, 1, 1));}
		{Histo1DPtr tmp; _h_pip.add(.18, .21, 	book(tmp, i+5, 1, 1));}
		{Histo1DPtr tmp; _h_pip.add(.21, .24, 	book(tmp, i+6, 1, 1));}
	// -> pi- X  		  theta ranges			d0a-x0b-y02
    	{Histo1DPtr tmp; _h_pim.add(.03, .06, 	book(tmp, i,   1, 2));}    	
    	{Histo1DPtr tmp; _h_pim.add(.06, .09, 	book(tmp, i+1, 1, 2));}    	
    	{Histo1DPtr tmp; _h_pim.add(.09, .12, 	book(tmp, i+2, 1, 2));}    	
    	{Histo1DPtr tmp; _h_pim.add(.12, .15, 	book(tmp, i+3, 1, 2));}    	
    	{Histo1DPtr tmp; _h_pim.add(.15, .18, 	book(tmp, i+4, 1, 2));}    	
    	{Histo1DPtr tmp; _h_pim.add(.18, .21, 	book(tmp, i+5, 1, 2));}    	
    	{Histo1DPtr tmp; _h_pim.add(.21, .24, 	book(tmp, i+6, 1, 2));} 
    	
    	declare(UnstableParticles(), "UFS");
    }

    void analyze(const Event& event) {

	const double weight = 1.0;
	const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");

      	for (const Particle& p : ufs.particles()) {
      	
      		const Vector3& v = p.momentum().p3();
      		//cout << "p.pid() = " << p.pid() << ", v = " << v << ", v.theta() = " << v.theta() << "\n" << endl;
      	
      		if (p.pid() == 211) {
      			//cout << "entered PIPLUS, v.theta():" << std::setprecision(4) << v.theta() << ", v.mod():" << v.mod()/GeV << "\n" << endl;
      			_h_pip.fill(v.theta(), v.mod()/GeV, weight);
      		}
      		else if (p.pid() == -211) { 
      			//cout << "entered PIMINUS, v.theta():" << std::setprecision(4) << v.theta() << ", v.mod():" << v.mod()/GeV << "\n" << endl;
      			_h_pim.fill(v.theta(), v.mod()/GeV, weight);
      		}
      	  	else {
	    		continue;
	    	}
      	}
    }

    void finalize() {
    	
    	double scale_factor = crossSection()/millibarn/sumOfWeights();
    	
	for (Histo1DPtr h : _h_pip.histos()) h->scaleW(scale_factor);
	for (Histo1DPtr h : _h_pim.histos()) h->scaleW(scale_factor);

    }

    private:

    BinnedHistogram _h_pip, _h_pim;

  };


  RIVET_DECLARE_PLUGIN(HARP_2008_I778842);

}
