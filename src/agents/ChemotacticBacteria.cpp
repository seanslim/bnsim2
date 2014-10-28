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

#include "headers/ChemotacticBacteria.h"

namespace BNSim {

ChemotacticBacteria::ChemotacticBacteria(const myVector3d& absPosition, double radius, double shovek, double velocity) :
		Agent(absPosition, radius, shovek), _velocity(velocity), _active(true) {

	// add a simple chemotaxis network
    ChemotaxisSystem * chemopathway = new ChemotaxisSystem(this);

    _nets.push_back(chemopathway);

	_volume = _total_volume = 4 * 3.141592653 * _cell_radius * _cell_radius
			* _cell_radius / 3;
    _total_radius =_cell_radius;
}

ChemotacticBacteria::~ChemotacticBacteria() {

}

void ChemotacticBacteria::update() {
    
	if (_active) {
        myVector3d velocity1;
        velocity1.scale(_velocity*CONFIG::timestep, ((ChemotaxisSystem * )_nets[0])->getDirection()); // 30 micrometers/sec
        _absPosition.add(velocity1);
		Agent::update();
	}
}


} /* namespace BNSim */
