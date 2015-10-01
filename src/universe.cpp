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

#include "universe.h"

namespace BNSim {

unsigned int CONFIG::worldSizeX = 0, CONFIG::worldSizeY = 0,
		CONFIG::worldSizeZ = 0;
unsigned int CONFIG::gridNumberX = 0, CONFIG::gridNumberY = 0,
		CONFIG::gridNumberZ = 0;
unsigned int CONFIG::gridSizeX = 0, CONFIG::gridSizeY = 0,
		CONFIG::gridSizeZ = 0;
double CONFIG::timestep = 1, CONFIG::simLength = 0, CONFIG::time = 0;
std::size_t CONFIG::numberMoleculeSpecies = 0;
std::size_t CONFIG::threadNumber = 1;
std::string CONFIG::workdir = "";
Universe* CONFIG::universe = NULL;
unsigned int CONFIG::boundaryLayerThickness = 0;
bool CONFIG::diffusion = true;

Universe::Universe() {

	// Calculate some constants that are useful

	CONFIG::gridSizeX = CONFIG::worldSizeX / CONFIG::gridNumberX;
	CONFIG::gridSizeY = CONFIG::worldSizeY / CONFIG::gridNumberY;
	CONFIG::gridSizeZ = CONFIG::worldSizeZ / CONFIG::gridNumberZ;
	std::size_t gridNum = CONFIG::gridNumberX * CONFIG::gridNumberY
			* CONFIG::gridNumberZ;

	// Add the grids

	unsigned int i = 0;
	while (i++ < gridNum) {
		_Grids.push_back(new Grid(i - 1));
	}

	// Prepare multi-threading
	thr = new pthread_t[CONFIG::threadNumber];
	thr_data = new thread_data_t[CONFIG::threadNumber];
	evn_thr_data = new thread_data_t[CONFIG::threadNumber];

	CONFIG::universe = this;
	IDcounts = 0;
	_Agents.setCapacity(10000000);
}

Universe::~Universe() {

	// Release maps
	for (std::map<std::string, MoleculeInfo*>::iterator itr =
			_moleculeMAP.begin(); itr != _moleculeMAP.end(); ++itr) {
		delete (*itr).second;
	}
	_moleculeMAP.clear();
	_moleculeMAPIndexed.clear();

	// Release grids
	for (std::vector<Grid*>::iterator itr = _Grids.begin(); itr != _Grids.end();
			++itr) {
		delete (*itr);
	}

	_Grids.clear();
}

void * environment_thread(void *arg) {
	thread_data_t *data = (thread_data_t *) arg;
	unsigned int startIndex = data->start;
	unsigned int endIndex = data->end;

	int xEast, xWest, yNorth, ySouth, zUp, zDown;
	// Sim::univ->Sim::cubesList of chemical leaving the grid in the positive (negative) .. direction
	double movingQuant;

	/**
	 * Flux of molecules crossing in the positive x-direction (Fick's law)
	 * 	J = -D(dC/dx) = -D*(C(x+dx)-C(x))/dx =  -D*(N(x+dx)-N(x))/((dx)^2*dy*dz)  molecules/(micron)^2/sec
	 * Number of molecules transferred in the positive x-direction over dt
	 * 	xAbove = J*(dy*dz)*dt = -((D*dt)/(dx)^2)*(N(x+dx)-N(x)) = -kX*(N(x+dx)-N(x))
	 * where kX = (D*dt)/(dx)^2 is a dimensionless constant
	 */
	double normX = CONFIG::timestep / pow(CONFIG::gridSizeX, 2);
	double normY = CONFIG::timestep / pow(CONFIG::gridSizeY, 2);
	double normZ = CONFIG::timestep / pow(CONFIG::gridSizeZ, 2);

	double temp = 0;

	for (unsigned int p = 0; p != CONFIG::numberMoleculeSpecies; p++) {

		const MoleculeInfo* pInfo = CONFIG::universe->getMoleculeInfo(p);

		// Shuffle order first

		unsigned int orderX[endIndex - startIndex];
		unsigned int orderY[CONFIG::gridNumberY];
		unsigned int orderZ[CONFIG::gridNumberZ];

		for (unsigned int i = startIndex; i != endIndex; ++i)
			orderX[i - startIndex] = i;

		for (unsigned int i = 0; i != CONFIG::gridNumberY; ++i)
			orderY[i] = i;

		for (unsigned int i = 0; i != CONFIG::gridNumberZ; ++i)
			orderZ[i] = i;

		for (unsigned int i = 0; i != (endIndex - startIndex); ++i) {
			size_t j = i
					+ rand() / (RAND_MAX / ((endIndex - startIndex) - i) + 1);
			int t = orderX[j];
			orderX[j] = orderX[i];
			orderX[i] = t;
		}

		for (unsigned int i = 0; i != CONFIG::gridNumberY; ++i) {
			size_t j = i + rand() / (RAND_MAX / (CONFIG::gridNumberY - i) + 1);
			int t = orderY[j];
			orderY[j] = orderY[i];
			orderY[i] = t;
		}

		for (unsigned int i = 0; i != CONFIG::gridNumberZ; ++i) {
			size_t j = i + rand() / (RAND_MAX / (CONFIG::gridNumberZ - i) + 1);
			int t = orderZ[j];
			orderZ[j] = orderZ[i];
			orderZ[i] = t;
		}

		// Diffusion

		for (unsigned int x = startIndex; x != endIndex; ++x)
			for (unsigned int y = 0; y != CONFIG::gridNumberY; ++y)
				for (unsigned int z = 0; z != CONFIG::gridNumberZ; ++z) {

					unsigned int i = orderX[x - startIndex];
					unsigned int j = orderY[y];
					unsigned int k = orderZ[z];

					xEast = (i == CONFIG::gridNumberX - 1 ? -1 : i + 1);
					xWest = (i == 0 ? -1 : i - 1);
					yNorth = (j == CONFIG::gridNumberY - 1 ? -1 : j + 1);
					ySouth = (j == 0 ? -1 : j - 1);
					zUp = (k == CONFIG::gridNumberZ - 1 ? -1 : k + 1);
					zDown = (k == 0 ? -1 : k - 1);

					Grid* thisGrid = CONFIG::universe->getGrid(i, j, k);

					if (thisGrid->getLayerType() == bulk)
						continue;

					Grid* neighbor = NULL;

					if (xEast != -1) {

						neighbor = CONFIG::universe->getGrid(xEast, j, k);

						double kX = normX * 0.5
								* (pInfo->getDiffusionCoefficient(
										thisGrid->getLayerType())
										+ pInfo->getDiffusionCoefficient(
												neighbor->getLayerType()));

						if (kX > 1) {
							std::cout
									<< "Incorrect numerical setup detected, may cause numerical instability"
									<< std::endl;
						}

						/* Calculate the Sim::cubesList of chemical leaving the grid in this direction */
						movingQuant = -kX
								* (neighbor->getConc(p) - thisGrid->getConc(p));

						/* Adjust the amount to avoid numerical instability*/
						double minQuant = 0;

						if (movingQuant > 0)
							minQuant =
									movingQuant > thisGrid->getConc(p) ?
											thisGrid->getConc(p) : movingQuant;
						else
							minQuant =
									movingQuant > neighbor->getConc(p) ?
											neighbor->getConc(p) : movingQuant;

						if (neighbor->getLayerType() != bulk) {
							temp = neighbor->getConc(p) + minQuant;
							/* Add that Sim::cubesList to the grid above */
							neighbor->setConc(p, (temp > 0 ? temp : 0));
						}

						/* Remove it from this grid */
						temp = thisGrid->getConc(p) - minQuant;
						thisGrid->setConc(p, (temp > 0 ? temp : 0));
					}
					if (xWest != -1) {

						neighbor = CONFIG::universe->getGrid(xWest, j, k);

						double kX = normX * 0.5
								* (pInfo->getDiffusionCoefficient(
										thisGrid->getLayerType())
										+ pInfo->getDiffusionCoefficient(
												neighbor->getLayerType()));

						if (kX > 1) {
							std::cout
									<< "Incorrect numerical setup detected, may cause numerical instability"
									<< std::endl;
						}

						movingQuant = -kX
								* (neighbor->getConc(p) - thisGrid->getConc(p));

						double minQuant = 0;

						if (movingQuant > 0)
							minQuant =
									movingQuant > thisGrid->getConc(p) ?
											thisGrid->getConc(p) : movingQuant;
						else
							minQuant =
									movingQuant > neighbor->getConc(p) ?
											neighbor->getConc(p) : movingQuant;

						if (neighbor->getLayerType() != bulk) {
							temp = neighbor->getConc(p) + minQuant;
							neighbor->setConc(p, (temp > 0 ? temp : 0));
						}
						temp = thisGrid->getConc(p) - minQuant;
						thisGrid->setConc(p, (temp > 0 ? temp : 0));
					}
					if (yNorth != -1) {

						neighbor = CONFIG::universe->getGrid(i, yNorth, k);

						double kY = normY * 0.5
								* (pInfo->getDiffusionCoefficient(
										thisGrid->getLayerType())
										+ pInfo->getDiffusionCoefficient(
												neighbor->getLayerType()));

						if (kY > 1) {
							std::cout
									<< "Incorrect numerical setup detected, may cause numerical instability"
									<< std::endl;
						}

						movingQuant = -kY
								* (neighbor->getConc(p) - thisGrid->getConc(p));

						double minQuant = 0;

						if (movingQuant > 0)
							minQuant =
									movingQuant > thisGrid->getConc(p) ?
											thisGrid->getConc(p) : movingQuant;
						else
							minQuant =
									movingQuant > neighbor->getConc(p) ?
											neighbor->getConc(p) : movingQuant;

						if (neighbor->getLayerType() != bulk) {
							temp = neighbor->getConc(p) + minQuant;
							neighbor->setConc(p, (temp > 0 ? temp : 0));
						}

						temp = thisGrid->getConc(p) - minQuant;
						thisGrid->setConc(p, (temp > 0 ? temp : 0));
					}
					if (ySouth != -1) {

						neighbor = CONFIG::universe->getGrid(i, ySouth, k);

						double kY = normY * 0.5
								* (pInfo->getDiffusionCoefficient(
										thisGrid->getLayerType())
										+ pInfo->getDiffusionCoefficient(
												neighbor->getLayerType()));

						if (kY > 1) {
							std::cout
									<< "Incorrect numerical setup detected, may cause numerical instability"
									<< std::endl;
						}

						movingQuant = -kY
								* (neighbor->getConc(p) - thisGrid->getConc(p));

						double minQuant = 0;

						if (movingQuant > 0)
							minQuant =
									movingQuant > thisGrid->getConc(p) ?
											thisGrid->getConc(p) : movingQuant;
						else
							minQuant =
									movingQuant > neighbor->getConc(p) ?
											neighbor->getConc(p) : movingQuant;

						if (neighbor->getLayerType() != bulk) {
							temp = neighbor->getConc(p) + minQuant;
							neighbor->setConc(p, (temp > 0 ? temp : 0));
						}

						temp = thisGrid->getConc(p) - minQuant;
						thisGrid->setConc(p, (temp > 0 ? temp : 0));
					}
					if (zUp != -1) {

						neighbor = CONFIG::universe->getGrid(i, j, zUp);

						double kZ = normZ * 0.5
								* (pInfo->getDiffusionCoefficient(
										thisGrid->getLayerType())
										+ pInfo->getDiffusionCoefficient(
												neighbor->getLayerType()));

						if (kZ > 1) {
							std::cout
									<< "Incorrect numerical setup detected, may cause numerical instability"
									<< std::endl;
						}

						movingQuant = -kZ
								* (neighbor->getConc(p) - thisGrid->getConc(p));

						double minQuant = 0;

						if (movingQuant > 0)
							minQuant =
									movingQuant > thisGrid->getConc(p) ?
											thisGrid->getConc(p) : movingQuant;
						else
							minQuant =
									movingQuant > neighbor->getConc(p) ?
											neighbor->getConc(p) : movingQuant;

						if (neighbor->getLayerType() != bulk) {
							temp = neighbor->getConc(p) + minQuant;
							neighbor->setConc(p, (temp > 0 ? temp : 0));
						}

						temp = thisGrid->getConc(p) - minQuant;
						thisGrid->setConc(p, (temp > 0 ? temp : 0));
					}
					if (zDown != -1) {
						neighbor = CONFIG::universe->getGrid(i, j, zDown);

						double kZ = normZ * 0.5
								* (pInfo->getDiffusionCoefficient(
										thisGrid->getLayerType())
										+ pInfo->getDiffusionCoefficient(
												neighbor->getLayerType()));

						if (kZ > 1) {
							std::cout
									<< "Incorrect numerical setup detected, may cause numerical instability"
									<< std::endl;
						}

						movingQuant = -kZ
								* (neighbor->getConc(p) - thisGrid->getConc(p));

						double minQuant = 0;

						if (movingQuant > 0)
							minQuant =
									movingQuant > thisGrid->getConc(p) ?
											thisGrid->getConc(p) : movingQuant;
						else
							minQuant =
									movingQuant > neighbor->getConc(p) ?
											neighbor->getConc(p) : movingQuant;

						if (neighbor->getLayerType() != bulk) {
							temp = neighbor->getConc(p) + minQuant;
							neighbor->setConc(p, (temp > 0 ? temp : 0));
						}

						temp = thisGrid->getConc(p) - minQuant;
						thisGrid->setConc(p, (temp > 0 ? temp : 0));
					}

					// Molecular decay
					double newValue = ((1
							- pInfo->getDecayRate(thisGrid->getLayerType()))
							* thisGrid->getConc(p)) * CONFIG::timestep;
					thisGrid->setConc(p, newValue > 0 ? newValue : 0);
				}

	}

	pthread_exit(NULL);
}

void * agent_thread(void *arg) {
	thread_data_t *data = (thread_data_t *) arg;

	unsigned int order[data->end - data->start];

	for (unsigned int i = data->start; i != data->end; ++i)
		order[i - data->start] = i;

	for (unsigned int i = 0; i != (data->end - data->start); ++i) {
		size_t j = i
				+ rand() / (RAND_MAX / ((data->end - data->start) - i) + 1);
		int t = order[j];
		order[j] = order[i];
		order[i] = t;
	}

	for (unsigned int i = data->start; i != data->end; ++i) {
		unsigned int id = order[i - data->start];
		Agent * age = CONFIG::universe->getAgent(id);
		if (age != NULL)
			age->update();
	}
	pthread_exit(NULL);
}

// Handle agent delta movement

void * post_agent_thread(void *arg) {
	thread_data_t *data = (thread_data_t *) arg;

	for (unsigned int i = data->start; i != data->end; ++i) {
		if (CONFIG::universe->getAgent(i) != NULL)
			CONFIG::universe->getAgent(i)->pos_update();
	}

	pthread_exit(NULL);
}

void Universe::update_agent_parallel() {

	unsigned int i;
	int rc;

	unsigned int order[CONFIG::threadNumber];

	for (unsigned int i = 0; i != CONFIG::threadNumber; ++i)
		order[i] = i;

	for (unsigned int i = 0; i != CONFIG::threadNumber; ++i) {
		size_t j = i + rand() / (RAND_MAX / (CONFIG::threadNumber - i) + 1);
		int t = order[j];
		order[j] = order[i];
		order[i] = t;
	}

	for (i = 0; i != CONFIG::threadNumber; ++i) {
		if ((rc = pthread_create(&thr[i], NULL, agent_thread,
				&thr_data[order[i]]))) {
			return;
		}
	}
	/* block until all threads complete */
	for (i = 0; i != CONFIG::threadNumber; ++i) {
		pthread_join(thr[i], NULL);
	}

	for (i = 0; i != CONFIG::threadNumber; ++i) {
		if ((rc = pthread_create(&thr[i], NULL, post_agent_thread,
				&thr_data[order[i]]))) {
			return;
		}
	}
	/* block until all threads complete */
	for (i = 0; i != CONFIG::threadNumber; ++i) {
		pthread_join(thr[i], NULL);
	}
}

void Universe::update_environment_p() {
	unsigned int i;
	int rc;

	unsigned int order[CONFIG::threadNumber];

	for (unsigned int i = 0; i != CONFIG::threadNumber; ++i)
		order[i] = i;

	for (unsigned int i = 0; i != CONFIG::threadNumber; ++i) {
		size_t j = i + rand() / (RAND_MAX / (CONFIG::threadNumber - i) + 1);
		int t = order[j];
		order[j] = order[i];
		order[i] = t;
	}

	for (unsigned int i = 0; i != CONFIG::threadNumber; ++i) {
		if ((rc = pthread_create(&thr[i], NULL, environment_thread,
				&evn_thr_data[order[i]]))) {
			return;
		}
	}

	/* block until all threads complete */
	for (i = 0; i != CONFIG::threadNumber; ++i) {
		pthread_join(thr[i], NULL);
	}
}

void Universe::prepare_multithreading() {
	unsigned int threadstep = CONFIG::universe->getTotalAgentNumber()
			/ CONFIG::threadNumber;
	unsigned int envthreadstep = CONFIG::gridNumberX / CONFIG::threadNumber;

	// prepare data-structure

	for (unsigned int i = 0; i < CONFIG::threadNumber - 1; ++i) {
		thr_data[i].start = i * threadstep;
		thr_data[i].end = (i + 1) * threadstep;

		evn_thr_data[i].start = i * envthreadstep;
		evn_thr_data[i].end = (i + 1) * envthreadstep;

	}
	thr_data[CONFIG::threadNumber - 1].start =
			thr_data[CONFIG::threadNumber - 2].end;
	thr_data[CONFIG::threadNumber - 1].end =
			CONFIG::universe->getTotalAgentNumber();

	evn_thr_data[CONFIG::threadNumber - 1].start =
			evn_thr_data[CONFIG::threadNumber - 2].end;
	evn_thr_data[CONFIG::threadNumber - 1].end = CONFIG::gridNumberX;
}

void Universe::evolute() {
	prepare_multithreading();
	update_agent_parallel();

	if (CONFIG::diffusion)
		update_environment_p();

	CONFIG::time += CONFIG::timestep;
}

Grid* Universe::getGrid(unsigned int x, unsigned int y, unsigned int z) {
	return _Grids[x * CONFIG::gridNumberY * CONFIG::gridNumberZ
			+ y * CONFIG::gridNumberZ + z];
}

const MoleculeInfo* Universe::getMoleculeInfo(const std::string& name) {
	std::map<std::string, MoleculeInfo*>::iterator it = _moleculeMAP.find(name);

	if (it != _moleculeMAP.end()) {
		//element found;
		return (*it).second;
	} else
		return NULL;
}

const MoleculeInfo* Universe::getMoleculeInfo(const unsigned int index) {
	std::map<unsigned int, MoleculeInfo*>::iterator it =
			_moleculeMAPIndexed.find(index);

	if (it != _moleculeMAPIndexed.end()) {
		//element found;
		return (*it).second;
	} else
		return NULL;
}

Grid* Universe::getEastGrid(unsigned int x, unsigned int y, unsigned int z) {
	int xEast = (x == CONFIG::gridNumberX - 1 ? -1 : x + 1);

	if (xEast == -1)
		return NULL;
	else
		return getGrid(xEast, y, z);
}

Grid* Universe::getWestGrid(unsigned int x, unsigned int y, unsigned int z) {
	int xWest = (x == 0 ? -1 : x - 1);

	if (xWest == -1)
		return NULL;
	else
		return getGrid(xWest, y, z);
}

Grid* Universe::getUpGrid(unsigned int x, unsigned int y, unsigned int z) {
	int yEast = (y == CONFIG::gridNumberY - 1 ? -1 : y + 1);

	if (yEast == -1)
		return NULL;
	else
		return getGrid(x, yEast, z);
}

Grid* Universe::getDownGrid(unsigned int x, unsigned int y, unsigned int z) {
	int yWest = (y == 0 ? -1 : y - 1);

	if (yWest == -1)
		return NULL;
	else
		return getGrid(x, yWest, z);
}

Grid* Universe::getNorthGrid(unsigned int x, unsigned int y, unsigned int z) {
	int zEast = (z == CONFIG::gridNumberZ - 1 ? -1 : z + 1);

	if (zEast == -1)
		return NULL;
	else
		return getGrid(x, y, zEast);
}

Grid* Universe::getSouthGrid(unsigned int x, unsigned int y, unsigned int z) {
	int zWest = (z == 0 ? -1 : z - 1);

	if (zWest == -1)
		return NULL;
	else
		return getGrid(x, y, zWest);
}

}
