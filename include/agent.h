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

#ifndef AGENT_H_
#define AGENT_H_

#include "common.h"
#include "regulatoryNet.h"
#include "configuration.h"
#include "spacegrid.h"
#include <vector>
#include <cmath>

namespace BNSim {

class RegulatoryNet;
class Grid;

/*
 * Agent represents all sorts of particles in the universe,
 * such as bacteria and EPS capsules
 */

class Agent {

public:
	Agent(const myVector3d& absPosition, double radius, double shovek);
	virtual ~Agent();
	virtual void update();
	virtual void pos_update();
	intVector3d getGridPos() {
		return _gridPosition;
	}

	double getShovek() {
		return _shovek;
	}
	myVector3d& getabsPosition() {
		return _absPosition;
	}
	double getDistance(Agent& agentB);
	bool updateMass(const std::string& name, double mass);
	double getMass(const std::string& name);
	double getCellRadius() const {
		return _cell_radius;
	}
	double getTotalRadius() const {
		return _total_radius;
	}
	void updateVolume();
	void updateRadius() {
		_total_radius = pow(_total_volume * 0.75 / 3.141592653, 0.33);
		_cell_radius = pow(_volume * 0.75 / 3.141592653, 0.33);
	}
	unsigned long getID() const {
		return ID;
	}
	RegulatoryNet * getRegulatoryNet(const std::string& name);
protected:
	myVector3d _absPosition, _deltaMovement;
	intVector3d _gridPosition;
	std::vector<RegulatoryNet *> _nets;
	std::vector<massInfo *> massInfos;
	void updateGridPos();
	void registerGridPos();
	void deRegisterGridPos();
	void shove();
	void addMass(const std::string& name, double mass, double density);
	double _total_radius, _cell_radius, _shovek, _total_volume, _volume; // volume: vol without EPS
	unsigned long ID;
private:
	Grid* myGrid;
};

}

#endif /* AGENT_H_ */
