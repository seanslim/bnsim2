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

#ifndef SPACEGRID_H_
#define SPACEGRID_H_

#include <vector>
#include "../../agents/headers/agent.h"
#include "../../core/configuration.h"
#include "../../utils/mutexlock.h"
#include "common.h"
#include <mutex>

namespace BNSim {

class Agent;

class Grid {

public:
	Grid(const int gridIndex);
	~Grid();
	const double getConc(const unsigned int moleculeSpeciesIndex) { return _moleculeSpeciesConc[moleculeSpeciesIndex];}
    void setConc(const unsigned int moleculeSpeciesIndex, const double conc);
    void consumeChemical(const unsigned int moleculeSpeciesIndex, const double conc);
	void updateParticles();
	const std::size_t getAgentNumber() { return _agents->getCount(); }
	void addAgent(Agent* p) ;
	void deleteAgent(Agent* p);
	void deltaChemical(const unsigned int moleculeSpeciesIndex, const double mass);
	layerType getLayerType() {return _type;}
	void setLayerType(layerType type) {_type = type; }
	BNSimVector<Agent*> *_agents;
private:
	unsigned int _gridIndex;
	double *_moleculeSpeciesConc;
	std::mutex locker;
	layerType _type;
	double volume;
};

}


#endif /* SPACEGRID_H_ */
