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

#include "spacegrid.h"

using namespace BNSim;



Grid::Grid(const int gridIndex):_gridIndex(gridIndex)
{
	_moleculeSpeciesConc = new double[CONFIG::numberMoleculeSpecies];

	for(unsigned int i = 0; i != CONFIG::numberMoleculeSpecies; ++i)
		_moleculeSpeciesConc[i] = 0;

	_agents = new BNSimVector<Agent*>();
	_type = bulk;   // all grids initialized to bulk type

	volume = CONFIG::gridSizeZ*CONFIG::gridSizeY*CONFIG::gridSizeZ;
}


Grid::~Grid()
{
    for(unsigned int i = 0; i!= _agents->getSize();++i) {
		delete (*_agents)[i];
	}

	delete _agents;
	delete [] _moleculeSpeciesConc;
}

void Grid::updateParticles()
{
	// todo
}

void Grid::setConc(const unsigned int moleculeSpeciesIndex, const double conc){
	mutexLock lock(locker);
     _moleculeSpeciesConc[moleculeSpeciesIndex] = conc>0?conc:0;
}

void Grid::deltaChemical(const unsigned int moleculeSpeciesIndex, const double mass) {

	mutexLock lock(locker);
	double orgmass = _moleculeSpeciesConc[moleculeSpeciesIndex]*volume;

	_moleculeSpeciesConc[moleculeSpeciesIndex] = (orgmass+mass)/volume;

	if(_moleculeSpeciesConc[moleculeSpeciesIndex]<0)
        _moleculeSpeciesConc[moleculeSpeciesIndex] = 0;
}

void Grid::consumeChemical(const unsigned int moleculeSpeciesIndex, const double conc) {
    mutexLock lock(locker);
    _moleculeSpeciesConc[moleculeSpeciesIndex] -= conc;
    
    if(_moleculeSpeciesConc[moleculeSpeciesIndex]<0)
        _moleculeSpeciesConc[moleculeSpeciesIndex] = 0;
}

void Grid::addAgent(Agent* p)
{
	mutexLock lock(locker);

	_agents->add(p);
}

void Grid::deleteAgent(Agent* p)
{
    mutexLock lock(locker);
	(*_agents).remove(p);
}
