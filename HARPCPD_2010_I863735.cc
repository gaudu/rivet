// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Tools/BinnedHistogram.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {

  class HARPCPD_2010_I863735 : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(HARPCPD_2010_I863735);

    void init() {

	declare(UnstableParticles(Cuts::ptIn(0*GeV, 1.25*GeV)), "UFS");
	
/*		 		a) pi+ C			b) pi- C
	3 GeV/c	   	 30 - 53	i = 30		 54 - 77	i = 54
	5 GeV/c    	102 - 125	i = 102		126 - 149	i = 126
	8 GeV/c	   	174 - 197	i = 174		198 - 221	i = 198
	12 GeV/c   	246 - 269 	i = 246		270 - 293	i = 270
	15 GeV/c   	318 - 341	i = 318		342 - 365	i = 342
*/	
	const int i=246;
	cout << "i:" << i << "\n" << endl;
	
	//BinnedHistogram over polar angle (theta()) in degrees
	// -> p X
	{Histo1DPtr tmp; _h_p.add(	20.,  	30., 	book(tmp, i, 	1, 1));}
    	{Histo1DPtr tmp; _h_p.add(	30.,  	40., 	book(tmp, i+1, 	1, 1));}
    	{Histo1DPtr tmp; _h_p.add(	40.,  	50., 	book(tmp, i+2, 	1, 1));}
    	{Histo1DPtr tmp; _h_p.add(	50.,  	60., 	book(tmp, i+3, 	1, 1));}
    	{Histo1DPtr tmp; _h_p.add(	60.,  	75., 	book(tmp, i+4, 	1, 1));}
    	{Histo1DPtr tmp; _h_p.add(	75.,  	90., 	book(tmp, i+5, 	1, 1));}
    	{Histo1DPtr tmp; _h_p.add(	90.,  	105., 	book(tmp, i+6, 	1, 1));}
    	{Histo1DPtr tmp; _h_p.add(	105., 	125., 	book(tmp, i+7, 	1, 1));}
    	
    	// -> pi+ X
    	{Histo1DPtr tmp; _h_pip.add(	20.,  	30., 	book(tmp, i+8, 	1, 1));}
    	{Histo1DPtr tmp; _h_pip.add(	30.,  	40., 	book(tmp, i+9, 	1, 1));}
    	{Histo1DPtr tmp; _h_pip.add(	40.,  	50., 	book(tmp, i+10,	1, 1));}
    	{Histo1DPtr tmp; _h_pip.add(	50.,  	60., 	book(tmp, i+11,	1, 1));}
    	{Histo1DPtr tmp; _h_pip.add(	60.,  	75., 	book(tmp, i+12, 1, 1));}
    	{Histo1DPtr tmp; _h_pip.add(	75.,  	90., 	book(tmp, i+13, 1, 1));}
    	{Histo1DPtr tmp; _h_pip.add(	90.,  	105., 	book(tmp, i+14, 1, 1));}
    	{Histo1DPtr tmp; _h_pip.add(	105., 	125., 	book(tmp, i+15, 1, 1));}
    	
    	// -> pi- X
    	{Histo1DPtr tmp; _h_pim.add(	20.,  	30., 	book(tmp, i+16,	1, 1));}
    	{Histo1DPtr tmp; _h_pim.add(	30.,  	40., 	book(tmp, i+17, 1, 1));}
    	{Histo1DPtr tmp; _h_pim.add(	40.,  	50., 	book(tmp, i+18, 1, 1));}
    	{Histo1DPtr tmp; _h_pim.add(	50.,  	60., 	book(tmp, i+19, 1, 1));}
    	{Histo1DPtr tmp; _h_pim.add(	60.,  	75., 	book(tmp, i+20, 1, 1));}
    	{Histo1DPtr tmp; _h_pim.add(	75.,  	90., 	book(tmp, i+21, 1, 1));}
    	{Histo1DPtr tmp; _h_pim.add(	90.,  	105., 	book(tmp, i+22, 1, 1));}
    	{Histo1DPtr tmp; _h_pim.add(	105., 	125., 	book(tmp, i+23, 1, 1));}
    }

    void analyze(const Event& event) {
    	//cout << "void analyze()" << "\n" << endl;
	
	const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      	for (const Particle& p : ufs.particles()) {
      		//cout << "particle()" << "\n" << endl;
      		
      		const Vector3& v = p.momentum().p3();
      		double weight = 1./(2.*pi*(p.pT()/GeV));
      		double thetadeg = (v.theta())*180./pi;
      		
      		if (p.pid() == PID::PROTON) {
      			//cout << "entered PROTON, v.theta()*180./pi:" << v.theta()*180./pi << ", p.pT()/GeV:" << p.pT()/GeV << "\n" << endl;
      			_h_p.fill(thetadeg, p.pT()/GeV, weight);
      		}
      		
      		if (p.pid() == PID::PIPLUS) {
      			//cout << "entered PIPLUS, v.theta()*180./pi:" << v.theta()*180./pi << ", p.pT()/GeV:" << p.pT()/GeV << "\n" << endl;
      			_h_pip.fill(thetadeg, p.pT()/GeV, weight);
      		}
      		
      		if (p.pid() == PID::PIMINUS) {
      			//cout << "entered PIMINUS v.theta()*180./pi:" << v.theta()*180./pi << ", p.pT()/GeV:" << p.pT()/GeV << "\n" << endl;
      			_h_pim.fill(thetadeg, p.pT()/GeV, weight);
      		}
      	}
      	
    }

    void finalize() {

	double factor = crossSection()/millibarn/sumOfWeights();
	
	for (Histo1DPtr h : _h_p.histos()) 	h->scaleW(factor);
	for (Histo1DPtr h : _h_pip.histos()) 	h->scaleW(factor);
	for (Histo1DPtr h : _h_pim.histos()) 	h->scaleW(factor);
	
    }
    
    private:

    BinnedHistogram _h_p, _h_pip, _h_pim;


  };


  RIVET_DECLARE_PLUGIN(HARPCPD_2010_I863735);

}
