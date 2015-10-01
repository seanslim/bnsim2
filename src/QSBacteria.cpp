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

#include "QSBacteria.h"

namespace BNSim {

QSBacteria::QSBacteria(const myVector3d& absPosition, double radius, double shovek, double T_div,
		double T_eps, double T_die) :
		Agent(absPosition, radius, shovek), _T_div(T_div), _T_eps(T_eps), _T_die(
				T_die), _active(true) {

	// add a Lux-type QS network
	QSLux * qs = new QSLux(this);
	_nets.push_back(qs);

	// add a simple metabolism network
            SimpleMetabolism * meta = new SimpleMetabolism(this);
            meta->setQs(qs);
            meta->setKc(0.001);
            meta->setKs(1);
            meta->setMaintenance(1e-4);
            meta->setuMax(0.002);
            meta->setPgMax(0.0001);
            meta->setY(0.85);
            _nets.push_back(meta);

	// add cell mass
	addMass("X", 10, 150);   // density g.L-1

	_volume = _total_volume = 4 * 3.141592653 * _cell_radius * _cell_radius
			* _cell_radius / 3;

	for (std::vector<massInfo *>::iterator itr = massInfos.begin();
			itr != massInfos.end(); itr++) {
		if ((*itr)->name == "X")
			(*itr)->mass = _volume * (*itr)->density;
	}

	// add EPS mass, not adding anything now for the ICC paper
//	addMass("EPS", 1, 75);
}

QSBacteria::~QSBacteria() {

}

void QSBacteria::update() {
	if (_active) {
		Agent::update();
		divide();
	}
}

void QSBacteria::divide() {

	if (_cell_radius < _T_die) {
		// Remain in the space, dead mass
		_active = false;
	}

	if (_cell_radius > _T_div) {

		// place my daughter near me
		myVector3d Position(_absPosition);
		Position.pos.x = Position.pos.x + 0.01 * ((double) rand() / (RAND_MAX))
				- 0.005;
		Position.pos.y = Position.pos.y + 0.01 * ((double) rand() / (RAND_MAX))
				- 0.005;
		Position.pos.z = Position.pos.z + 0.01 * ((double) rand() / (RAND_MAX))
				- 0.005;

		// give birth to my daughter
		Agent* newBac = new QSBacteria(Position, 1.7, 1.2, 2, 2, 0.2);
		newBac->updateMass("X", getMass("X") * 0.45); // 0.05 for the loss in the division process
		newBac->updateVolume();
		newBac->updateRadius();
		CONFIG::universe->addAgent(newBac);

		// cut my mass by half
		updateMass("X", getMass("X") * 0.45);
		updateVolume();
		updateRadius();
	}

	// release an EPS agent from the capsule volume

	if (_total_radius - _cell_radius > _T_eps) {

		myVector3d Position(_absPosition);
		Position.pos.x = Position.pos.x + 0.01 * ((double) rand() / (RAND_MAX))
				- 0.005;
		Position.pos.y = Position.pos.y + 0.01 * ((double) rand() / (RAND_MAX))
				- 0.005;
		Position.pos.z = Position.pos.z + 0.01 * ((double) rand() / (RAND_MAX))
				- 0.005;

		Agent* newEPS = new EPS(Position, protein, _total_radius - _cell_radius,
				1.2);
		CONFIG::universe->addAgent(newEPS);

		updateMass("EPS", 0);
		updateVolume();
		updateRadius();
	}
}

} /* namespace BNSim */
