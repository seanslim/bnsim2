/**
 * BNSim is an open-source, parallel, stochastic, and multi-scale modeling
 * platform which integrates various simulation algorithms, together with
 * chemotaixs, quorum sensing, and biofilm models in a 3D environment.
 *
 * BNSim is developed by CMU SLD Group, and released under GPLv2 license
 * Please check http://www.ece.cmu.edu/~sld/ for more information
 *
 * BNSim is developed using C++ and pthread under Linux and Mac OS
 *
 * If you use it or part of it for your work, please cite
 * Wei, Guopeng, Paul Bogdan, and Radu Marculescu. "Efficient Modeling and
 * Simulation of Bacteria-Based Nanonetworks with BNSim." Selected Areas in
 * Communications, IEEE Journal on 31.12 (2013): 868-878.
 *
 * @author      Guopeng (Daniel) Wei  1@weiguopeng.com
 * @version     2.0
 * @since       1.0
 */

#include "regulatoryNet.h"

namespace BNSim {

SimpleMetabolism::SimpleMetabolism(Agent* host) :
		RegulatoryNet(host) {
            name = "SimpleMetabolism";
}

SimpleMetabolism::~SimpleMetabolism() {
}

double SimpleMetabolism::randomUniform() {
	return ((double) rand() / (RAND_MAX));
}

void SimpleMetabolism::update() {

	int substrate = -1;
	substrate =
			CONFIG::universe->getMoleculeInfo("substrate") == NULL ?
					-1 :
					CONFIG::universe->getMoleculeInfo("substrate")->getIndex();

	intVector3d GPos = getHost()->getGridPos();
	Grid* myGrid = CONFIG::universe->getGrid(GPos.x, GPos.y, GPos.z);
	double timestep = CONFIG::timestep;

	// input to the system

    _s = myGrid->getConc(substrate);
    double C = _qs->getC1();
    double X = getHost()->getMass("X");
    //	double eps = getHost()->getMass("EPS");
 
    _u = ((_u_max - _u_PG) * _s / (_s + _Ks) - _m)*timestep;
    double deltaX = _u * X+0.1*_u * X*(randomUniform()-0.5);
    
    double consume =( (_u_max + _u_PG) * _s / (_s + _Ks) + _m)*X*timestep;
    
    X += deltaX;
    // consume substrate
    myGrid->consumeChemical(substrate, consume);
    
    // let's set the mass
    getHost()->updateMass("X",X);
    //	getHost()->updateMass("EPS",eps);
}

} /* namespace BNSim */
