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

#include "headers/agent.h"

namespace BNSim {

Agent::Agent(const myVector3d& absPosition, double radius, double shovek) :
		_absPosition(absPosition), _cell_radius(radius), _total_radius(radius), _shovek(
				shovek) {
	registerGridPos();

	_deltaMovement.pos.x = 0;
	_deltaMovement.pos.y = 0;
	_deltaMovement.pos.z = 0;

	ID = CONFIG::universe->retriveID();
}

Agent::~Agent() {
	for (std::vector<RegulatoryNet *>::iterator itr = _nets.begin();
			itr != _nets.end(); itr++)
		delete (*itr);
	_nets.clear();

	for (std::vector<massInfo *>::iterator itr = massInfos.begin();
			itr != massInfos.end(); itr++)
		delete (*itr);
	massInfos.clear();
}

RegulatoryNet* Agent::getRegulatoryNet(const std::string& name) {
    for (std::vector<RegulatoryNet *>::iterator itr = _nets.begin();
         itr != _nets.end(); itr++) {
        if((*itr)->getName()==name)
            return(*itr);
    }
    return NULL;
}
    
void Agent::update() {
	// update regulatory networks
	for (std::vector<RegulatoryNet *>::iterator itr = _nets.begin();
			itr != _nets.end(); itr++)
		(*itr)->update();

	shove();
	updateVolume();
	updateRadius();
	updateGridPos();
}

void Agent::pos_update() {
	_absPosition.add(_deltaMovement);
	_absPosition.add(_deltaMovement);

	// check boundaries
	if (_absPosition.pos.x > CONFIG::worldSizeX) {
		_absPosition.pos.x = CONFIG::worldSizeX - 1;
	}
	if (_absPosition.pos.x < 0) {
		_absPosition.pos.x = 1;
	}
	if (_absPosition.pos.y > CONFIG::worldSizeY) {
		_absPosition.pos.y = CONFIG::worldSizeY - 1;
	}
	if (_absPosition.pos.y < 0) {
		_absPosition.pos.y = 1;
	}
	if (_absPosition.pos.z > CONFIG::worldSizeZ) {
		_absPosition.pos.z = CONFIG::worldSizeZ - 1;
	}
	if (_absPosition.pos.z < 0) {
		_absPosition.pos.z = 1;
	}

	updateGridPos();

	_deltaMovement.pos.x = 0;
	_deltaMovement.pos.y = 0;
	_deltaMovement.pos.z = 0;
}

void Agent::updateVolume() {
	_total_volume = 0;
	_volume = 0;

	for (std::vector<massInfo *>::iterator itr = massInfos.begin();
			itr != massInfos.end(); itr++) {
		_total_volume += (*itr)->mass / (*itr)->density;
		if ((*itr)->name != "EPS")
			_volume += (*itr)->mass / (*itr)->density;
	}
}

double Agent::getDistance(Agent& agentB) {

	double dist = sqrt(
			pow(_absPosition.pos.x - agentB.getabsPosition().pos.x, 2)
					+ pow(_absPosition.pos.y - agentB.getabsPosition().pos.y, 2)
					+ pow(_absPosition.pos.z - agentB.getabsPosition().pos.z,
							2));

	return dist;
}

void Agent::deRegisterGridPos() {
	updateGridPos();
	myGrid->deleteAgent(this);
}

void Agent::registerGridPos() {
	intVector3d newPos;

	double x = (_absPosition.pos.x + _deltaMovement.pos.x)
			/ (double) (CONFIG::gridSizeX);
	double y = (_absPosition.pos.y + _deltaMovement.pos.y)
			/ (double) (CONFIG::gridSizeY);
	double z = (_absPosition.pos.z + _deltaMovement.pos.z)
			/ (double) (CONFIG::gridSizeZ);

	double remainderx = (int) (x * 100) % 100;
	double remaindery = (int) (y * 100) % 100;
	double remainderz = (int) (z * 100) % 100;

	if (remainderx >= 50)
		newPos.x = floor(x) + 1;
	else
		newPos.x = floor(x);
	if (remaindery >= 50)
		newPos.y = floor(y) + 1;
	else
		newPos.y = floor(y);
	if (remainderz >= 50)
		newPos.z = floor(z) + 1;
	else
		newPos.z = floor(z);

	if (newPos.x >= CONFIG::gridNumberX)
		newPos.x = CONFIG::gridNumberX - 1;
	if (newPos.y >= CONFIG::gridNumberY)
		newPos.y = CONFIG::gridNumberY - 1;
	if (newPos.z >= CONFIG::gridNumberZ)
		newPos.z = CONFIG::gridNumberZ - 1;

	_gridPosition = newPos;
	myGrid = CONFIG::universe->getGrid(newPos.x, newPos.y, newPos.z);
	myGrid->addAgent(this);
}

void Agent::updateGridPos() {

	intVector3d newPos;

	double x = (_absPosition.pos.x + _deltaMovement.pos.x)
			/ (double) (CONFIG::gridSizeX);
	double y = (_absPosition.pos.y + _deltaMovement.pos.y)
			/ (double) (CONFIG::gridSizeY);
	double z = (_absPosition.pos.z + _deltaMovement.pos.z)
			/ (double) (CONFIG::gridSizeZ);

	double remainderx = (int) (x * 100) % 100;
	double remaindery = (int) (y * 100) % 100;
	double remainderz = (int) (z * 100) % 100;

	if (remainderx >= 50)
		newPos.x = floor(x) + 1;
	else
		newPos.x = floor(x);
	if (remaindery >= 50)
		newPos.y = floor(y) + 1;
	else
		newPos.y = floor(y);
	if (remainderz >= 50)
		newPos.z = floor(z) + 1;
	else
		newPos.z = floor(z);

	if (newPos.x >= CONFIG::gridNumberX)
		newPos.x = CONFIG::gridNumberX - 1;
	if (newPos.y >= CONFIG::gridNumberY)
		newPos.y = CONFIG::gridNumberY - 1;
	if (newPos.z >= CONFIG::gridNumberZ)
		newPos.z = CONFIG::gridNumberZ - 1;

	if (!(newPos.x == _gridPosition.x && newPos.y == _gridPosition.y
			&& newPos.z == _gridPosition.z)) {
		CONFIG::universe->getGrid(_gridPosition.x, _gridPosition.y,
				_gridPosition.z)->deleteAgent(this);

		myGrid = CONFIG::universe->getGrid(newPos.x, newPos.y, newPos.z);

		myGrid->addAgent(this);
	}
	_gridPosition = newPos;
}

void Agent::shove() {
	// Get all the grids
	BNSimVector<Agent*> * gridsArray[7];
	gridsArray[0] = myGrid->_agents;
	gridsArray[1] =
			CONFIG::universe->getEastGrid(_gridPosition.x, _gridPosition.y,
					_gridPosition.z) == NULL ?
					NULL :
					CONFIG::universe->getEastGrid(_gridPosition.x,
							_gridPosition.y, _gridPosition.z)->_agents;
	gridsArray[2] =
			CONFIG::universe->getWestGrid(_gridPosition.x, _gridPosition.y,
					_gridPosition.z) == NULL ?
					NULL :
					CONFIG::universe->getWestGrid(_gridPosition.x,
							_gridPosition.y, _gridPosition.z)->_agents;
	gridsArray[3] =
			CONFIG::universe->getNorthGrid(_gridPosition.x, _gridPosition.y,
					_gridPosition.z) == NULL ?
					NULL :
					CONFIG::universe->getNorthGrid(_gridPosition.x,
							_gridPosition.y, _gridPosition.z)->_agents;
	gridsArray[4] =
			CONFIG::universe->getSouthGrid(_gridPosition.x, _gridPosition.y,
					_gridPosition.z) == NULL ?
					NULL :
					CONFIG::universe->getSouthGrid(_gridPosition.x,
							_gridPosition.y, _gridPosition.z)->_agents;
	gridsArray[5] =
			CONFIG::universe->getUpGrid(_gridPosition.x, _gridPosition.y,
					_gridPosition.z) == NULL ?
					NULL :
					CONFIG::universe->getUpGrid(_gridPosition.x,
							_gridPosition.y, _gridPosition.z)->_agents;
	gridsArray[6] =
			CONFIG::universe->getDownGrid(_gridPosition.x, _gridPosition.y,
					_gridPosition.z) == NULL ?
					NULL :
					CONFIG::universe->getDownGrid(_gridPosition.x,
							_gridPosition.y, _gridPosition.z)->_agents;

	// interactions in my grid and neighboring grids
	for (unsigned int i = 0; i != 7; ++i) {
		if (gridsArray[i] == NULL)
			continue;

		for (unsigned int j = 0; j < (*gridsArray[i]).getSize(); ++j) {
			Agent * agentB = (*gridsArray[i])[j];
            try {
                // make sure it's not me
                if (agentB == this || agentB == NULL || agentB==nullptr)
                    continue;
                
                double radiusSum = agentB->getTotalRadius() + _total_radius;
                
                double distance = getDistance(*agentB);
                
                // packed
                if (distance < radiusSum) {
                    myVector3d delta(_absPosition);
                    delta.sub(agentB->getabsPosition());
                    delta.normalize();
                    delta.scale(0.5 * (radiusSum - distance));
                    
                    _deltaMovement.add(delta);
                }
            } catch (...) {
                ;
            }
			
		}
	}
}

/*
 * Update biomass
 * return false if no record found
 */

bool Agent::updateMass(const std::string& name, double mass) {
	for (std::vector<massInfo *>::iterator itr = massInfos.begin();
			itr != massInfos.end(); itr++) {
		if (name == (*itr)->name) {
			(*itr)->mass = mass;
			return true;
		}
	}
	return false;
}

double Agent::getMass(const std::string& name) {
	for (std::vector<massInfo *>::iterator itr = massInfos.begin();
			itr != massInfos.end(); itr++) {
		if (name == (*itr)->name) {
			return (*itr)->mass;
		}
	}
	return -1;
}

void Agent::addMass(const std::string& name, double mass, double density) {
	massInfo* newMass = new massInfo;
	newMass->name = name;
	newMass->mass = mass;
	newMass->density = density;
	massInfos.push_back(newMass);
}

}

