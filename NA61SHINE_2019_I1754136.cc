// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Tools/BinnedHistogram.hh"

namespace Rivet {

  class NA61SHINE_2019_I1754136 : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(NA61SHINE_2019_I1754136);

    void init() {

      declare(UnstableParticles(), "UFS");
      
      // -> pi+X  		 theta ranges	       d01-x01-y0a
      {Histo1DPtr tmp; _h_pip.add(   0, .01, 	book(tmp, 1, 1, 1));}
      {Histo1DPtr tmp; _h_pip.add(.003, .01, 	book(tmp, 1, 1, 2));}
      {Histo1DPtr tmp; _h_pip.add( .01, .02, 	book(tmp, 1, 1, 3));}
      {Histo1DPtr tmp; _h_pip.add( .02, .04, 	book(tmp, 1, 1, 4));}
      {Histo1DPtr tmp; _h_pip.add( .04, .06, 	book(tmp, 1, 1, 5));}
      {Histo1DPtr tmp; _h_pip.add( .06, .1 , 	book(tmp, 1, 1, 6));}
      {Histo1DPtr tmp; _h_pip.add( .1 , .14, 	book(tmp, 1, 1, 7));}
      {Histo1DPtr tmp; _h_pip.add( .14, .18, 	book(tmp, 1, 1, 8));}
      {Histo1DPtr tmp; _h_pip.add( .18, .24, 	book(tmp, 1, 1, 9));}
      {Histo1DPtr tmp; _h_pip.add( .24, .3 ,  book(tmp, 1, 1, 10));}
      {Histo1DPtr tmp; _h_pip.add( .3 , .36, 	book(tmp, 1, 1, 11));}
      {Histo1DPtr tmp; _h_pip.add( .36, .42, 	book(tmp, 1, 1, 12));}
      
      // -> pi-X  		 theta ranges	       d02-x01-y0a
      {Histo1DPtr tmp; _h_pim.add(   0, .01, 	book(tmp, 2, 1, 1));}
      {Histo1DPtr tmp; _h_pim.add(.003, .01, 	book(tmp, 2, 1, 2));}
      {Histo1DPtr tmp; _h_pim.add( .01, .02, 	book(tmp, 2, 1, 3));}
      {Histo1DPtr tmp; _h_pim.add( .02, .04, 	book(tmp, 2, 1, 4));}
      {Histo1DPtr tmp; _h_pim.add( .04, .06, 	book(tmp, 2, 1, 5));}
      {Histo1DPtr tmp; _h_pim.add( .06, .1 , 	book(tmp, 2, 1, 6));}
      {Histo1DPtr tmp; _h_pim.add( .1 , .14, 	book(tmp, 2, 1, 7));}
      {Histo1DPtr tmp; _h_pim.add( .14, .18, 	book(tmp, 2, 1, 8));}
      {Histo1DPtr tmp; _h_pim.add( .18, .24, 	book(tmp, 2, 1, 9));}
      {Histo1DPtr tmp; _h_pim.add( .24, .3 ,  book(tmp, 2, 1, 10));}
      {Histo1DPtr tmp; _h_pim.add( .3 , .36, 	book(tmp, 2, 1, 11));}
      {Histo1DPtr tmp; _h_pim.add( .36, .42, 	book(tmp, 2, 1, 12));}
      
      // -> K+X  		theta ranges	       d03-x01-y0a
      {Histo1DPtr tmp; _h_kp.add(   0, .02, 	book(tmp, 3, 1, 1));}
      {Histo1DPtr tmp; _h_kp.add( .02, .04, 	book(tmp, 3, 1, 2));}
      {Histo1DPtr tmp; _h_kp.add( .04, .06, 	book(tmp, 3, 1, 3));}
      {Histo1DPtr tmp; _h_kp.add( .06, .1 , 	book(tmp, 3, 1, 4));}
      {Histo1DPtr tmp; _h_kp.add( .1 , .14, 	book(tmp, 3, 1, 5));}
      {Histo1DPtr tmp; _h_kp.add( .14, .18, 	book(tmp, 3, 1, 6));}
      {Histo1DPtr tmp; _h_kp.add( .18, .24, 	book(tmp, 3, 1, 7));}
      {Histo1DPtr tmp; _h_kp.add( .24, .3 ,  	book(tmp, 3, 1, 8));}
      {Histo1DPtr tmp; _h_kp.add( .3 , .36, 	book(tmp, 3, 1, 9));}
      
      // -> K-X  		theta ranges	       d04-x01-y0a
      {Histo1DPtr tmp; _h_km.add(   0, .02, 	book(tmp, 4, 1, 1));}
      {Histo1DPtr tmp; _h_km.add( .02, .04, 	book(tmp, 4, 1, 2));}
      {Histo1DPtr tmp; _h_km.add( .04, .06, 	book(tmp, 4, 1, 3));}
      {Histo1DPtr tmp; _h_km.add( .06, .1 , 	book(tmp, 4, 1, 4));}
      {Histo1DPtr tmp; _h_km.add( .1 , .14, 	book(tmp, 4, 1, 5));}
      {Histo1DPtr tmp; _h_km.add( .14, .18, 	book(tmp, 4, 1, 6));}
      {Histo1DPtr tmp; _h_km.add( .18, .24, 	book(tmp, 4, 1, 7));}
      {Histo1DPtr tmp; _h_km.add( .3 , .36, 	book(tmp, 4, 1, 8));}
      
      // -> pX  	       theta ranges	       d05-x01-y0a
      {Histo1DPtr tmp; _h_p.add(   0, .02, 	book(tmp, 5, 1, 1));}
      {Histo1DPtr tmp; _h_p.add( .02, .04, 	book(tmp, 5, 1, 2));}
      {Histo1DPtr tmp; _h_p.add( .04, .06, 	book(tmp, 5, 1, 3));}
      {Histo1DPtr tmp; _h_p.add( .06, .1 , 	book(tmp, 5, 1, 4));}
      {Histo1DPtr tmp; _h_p.add( .1 , .14, 	book(tmp, 5, 1, 5));}
      {Histo1DPtr tmp; _h_p.add( .14, .18, 	book(tmp, 5, 1, 6));}
      {Histo1DPtr tmp; _h_p.add( .18, .24, 	book(tmp, 5, 1, 7));}
      {Histo1DPtr tmp; _h_p.add( .24, .3 ,  book(tmp, 5, 1, 8));}
      {Histo1DPtr tmp; _h_p.add( .3 , .36, 	book(tmp, 5, 1, 9));}
      {Histo1DPtr tmp; _h_p.add( .36, .42, 	book(tmp, 5, 1, 10));}
      
      // -> K0sX  		 theta ranges	       d06-x01-y0a
      {Histo1DPtr tmp; _h_k0s.add(   0, .02, 	book(tmp, 6, 1, 1));}
      {Histo1DPtr tmp; _h_k0s.add( .02, .04, 	book(tmp, 6, 1, 2));}
      {Histo1DPtr tmp; _h_k0s.add( .04, .06, 	book(tmp, 6, 1, 3));}
      {Histo1DPtr tmp; _h_k0s.add( .06, .1 , 	book(tmp, 6, 1, 4));}
      {Histo1DPtr tmp; _h_k0s.add( .1 , .14, 	book(tmp, 6, 1, 5));}
      {Histo1DPtr tmp; _h_k0s.add( .14, .18, 	book(tmp, 6, 1, 6));}
      {Histo1DPtr tmp; _h_k0s.add( .18, .24, 	book(tmp, 6, 1, 7));}
      {Histo1DPtr tmp; _h_k0s.add( .24, .3 ,  book(tmp, 6, 1, 8));}
      {Histo1DPtr tmp; _h_k0s.add( .3 , .36, 	book(tmp, 6, 1, 9));}
      
      // -> lambdaX  		    theta ranges	       d07-x01-y0a
      {Histo1DPtr tmp; _h_lambda.add(   0, .02, 	book(tmp, 7, 1, 1));}
      {Histo1DPtr tmp; _h_lambda.add( .02, .04, 	book(tmp, 7, 1, 2));}
      {Histo1DPtr tmp; _h_lambda.add( .04, .06, 	book(tmp, 7, 1, 3));}
      {Histo1DPtr tmp; _h_lambda.add( .06, .1 , 	book(tmp, 7, 1, 4));}
      {Histo1DPtr tmp; _h_lambda.add( .1 , .14, 	book(tmp, 7, 1, 5));}
      {Histo1DPtr tmp; _h_lambda.add( .14, .18, 	book(tmp, 7, 1, 6));}
      {Histo1DPtr tmp; _h_lambda.add( .18, .24, 	book(tmp, 7, 1, 7));}
      {Histo1DPtr tmp; _h_lambda.add( .24, .3 ,  	book(tmp, 7, 1, 8));}
      {Histo1DPtr tmp; _h_lambda.add( .3 , .36, 	book(tmp, 7, 1, 9));}
      
      // -> antilambdaX 	        theta ranges	       d08-x01-y0a
      {Histo1DPtr tmp; _h_antilambda.add(   0, .02, 	book(tmp, 8, 1, 1));}
      {Histo1DPtr tmp; _h_antilambda.add( .02, .04, 	book(tmp, 8, 1, 2));}
      {Histo1DPtr tmp; _h_antilambda.add( .04, .06, 	book(tmp, 8, 1, 3));}
      {Histo1DPtr tmp; _h_antilambda.add( .06, .1 , 	book(tmp, 8, 1, 4));}
      {Histo1DPtr tmp; _h_antilambda.add( .1 , .14, 	book(tmp, 8, 1, 5));}
      {Histo1DPtr tmp; _h_antilambda.add( .14, .18, 	book(tmp, 8, 1, 6));}
      {Histo1DPtr tmp; _h_antilambda.add( .18, .24, 	book(tmp, 8, 1, 7));}

    }

    void analyze(const Event& event) {
	
      const UnstableParticles& ufs = apply<UnstableParticles> (event, "UFS");
	    for (const Particle& p : ufs.particles()) {
	
        const Vector3& v = p.momentum().p3();	
        double weight = 1.0;
      
        //cout << "p.pid(): " << p.pid() << " v.theta(): " << v.theta() << " v.mod()/GeV: " << v.mod()/GeV << endl;
    
        switch (p.pid()) {	
        case 211: {
          /*if (v.theta()<0.01):
            if (v.theta()>0.003 && v.mod()/GeV>33):
              _h_pip.fill(	*/			
          _h_pip.fill(v.theta(), v.mod()/GeV, weight);
          break;
        }
            
        case -211: {
          _h_pim.fill(v.theta(), v.mod()/GeV, weight);
          break;
        }

        case 321: {
          _h_kp.fill(v.theta(), v.mod()/GeV, weight);
          break;
        }
            
        case -321: {
          _h_km.fill(v.theta(), v.mod()/GeV, weight);
          break;
        }
            
        case 2212: {
          _h_p.fill(v.theta(), v.mod()/GeV, weight);
          break;
        }
            
        case 310: {
          _h_k0s.fill(v.theta(), v.mod()/GeV, weight);
          break;
        }
            
        case 3122: {
          _h_lambda.fill(v.theta(), v.mod()/GeV, weight);
          break;
        }
            
        case -3122: {
          _h_antilambda.fill(v.theta(), v.mod()/GeV, weight);
          break;
        }
        }
 	    }
    }

    void finalize() {

      //double scale_factor = 1./sumOfWeights(); 
    
      for (Histo1DPtr h :        _h_pip.histos()) 	normalize(h);
      for (Histo1DPtr h :        _h_pim.histos()) 	normalize(h);
      for (Histo1DPtr h :         _h_kp.histos()) 	normalize(h);
      for (Histo1DPtr h :         _h_km.histos())	  normalize(h);
      for (Histo1DPtr h :          _h_p.histos()) 	normalize(h);
      for (Histo1DPtr h :        _h_k0s.histos()) 	normalize(h);
      for (Histo1DPtr h :     _h_lambda.histos()) 	normalize(h);
      for (Histo1DPtr h : _h_antilambda.histos()) 	normalize(h);

    }

    private:
    
     BinnedHistogram _h_pip, _h_pim, _h_kp, _h_km, _h_p, _h_k0s, _h_lambda, _h_antilambda;

  };

  RIVET_DECLARE_PLUGIN(NA61SHINE_2019_I1754136);

}
