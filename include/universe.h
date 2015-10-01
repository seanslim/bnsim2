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

#ifndef UNIVERSE_H_
#define UNIVERSE_H_

#include<vector>
#include<string>
#include<map>
#include<cmath>
#include<iostream>
#include <stdlib.h>
#include"spacegrid.h"
#include"moleculeInfo.h"
#include"agent.h"
#include"configuration.h"

namespace BNSim {

// Forward declaration
class Grid;
class Agent;

struct thread_data_t {
	unsigned int start;
	unsigned int end;
};

class Universe {
public:
	Universe();
	virtual ~Universe();
	void addMoleculeSpecies(MoleculeInfo* moleculeInfo) { _moleculeMAP[moleculeInfo->getName()]=moleculeInfo; _moleculeMAPIndexed[moleculeInfo->getIndex()]=moleculeInfo;}
	const MoleculeInfo* getMoleculeInfo (const std::string& name) ;
	const MoleculeInfo* getMoleculeInfo (const unsigned int index) ;
	void diffuse();
	Grid* getGrid(unsigned int GridIndex) { return _Grids[GridIndex]; }
	Grid* getGrid(unsigned int x, unsigned int y, unsigned int z);
	void addAgent(Agent* agent) { _Agents.add(agent);}
    Agent* getAgent(unsigned int AgentIndex) { if(AgentIndex>=_Agents.getSize()) return NULL; else return _Agents[AgentIndex]; }
	std::size_t getTotalAgentNumber() { return _Agents.getSize();}
	void evolute();
	Grid* getEastGrid(unsigned int x, unsigned int y, unsigned int z);
	Grid* getWestGrid(unsigned int x, unsigned int y, unsigned int z);
	Grid* getUpGrid(unsigned int x, unsigned int y, unsigned int z);
	Grid* getDownGrid(unsigned int x, unsigned int y, unsigned int z);
	Grid* getNorthGrid(unsigned int x, unsigned int y, unsigned int z);
	Grid* getSouthGrid(unsigned int x, unsigned int y, unsigned int z);
    unsigned long retriveID() {return IDcounts++;}
private:
	std::vector<Grid*> _Grids;
	BNSimVector<Agent*> _Agents;
	std::map<std::string,MoleculeInfo*> _moleculeMAP;
	std::map<unsigned int,MoleculeInfo*> _moleculeMAPIndexed;
	pthread_t *thr;
	thread_data_t *thr_data;
	thread_data_t *evn_thr_data;
	void prepare_multithreading();
	void update_agent_parallel();
	void update_environment_p();
    unsigned long IDcounts;
};

}

#endif /* UNIVERSE_H_ */
