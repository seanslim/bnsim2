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

#include "headers/EPS.h"

namespace BNSim {

EPS::EPS(const myVector3d& absPosition, EPSType type, double radius,
		double shovek) :
		Agent(absPosition, radius, shovek), _type(type) {

	// add cell mass
	addMass("EPS", 10, 75);   // density g.L-1

	// initialize EPS mass volume
	_volume = _total_volume = 4 * 3.141592653 * radius * radius * radius / 3;

	// EPS has only one type of mass
	massInfos[0]->mass = _volume * massInfos[0]->density;

	updateRadius();
}

EPS::~EPS() {

}

double EPS::getEPSMass() {
	return massInfos[0]->mass;
}

void EPS::update() {
	Agent::update();
}

} /* namespace BNSim */
