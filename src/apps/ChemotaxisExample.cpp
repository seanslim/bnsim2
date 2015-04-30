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

#include <iostream>
#include "../core/configuration.h"
#include "../core/headers/universe.h"
#include "../datastructures/headers/moleculeInfo.h"
#include "../exporters/regularExporters.h"
#include "../agents/headers/ChemotacticBacteria.h"
#include <ctime>
#include <cstring>

using namespace std;
using namespace BNSim;

// Quorum sensing example

int main(int argn, char **argv) {

	time_t start, end;
	time(&start);

	cout << "BNSim 2.0 -- Chamotaxis Example" << endl;

	//---------------------------------------Setup Critical Parameters-----------------------------------

	CONFIG::workdir = "/home/daniel/test/";

	CONFIG::worldSizeX = 1000; // 0.5 cm
	CONFIG::worldSizeY = 500;
	CONFIG::worldSizeZ = 150;

	CONFIG::gridNumberX = 10;
	CONFIG::gridNumberY = 250;
	CONFIG::gridNumberZ = 10;

	CONFIG::numberMoleculeSpecies = 1;
	CONFIG::threadNumber = 8;
	CONFIG::timestep = 0.1;  // 0.4 seconds
	CONFIG::simLength = 10 * 60;  // 60 mins

	for (int i = 1; i < argn; i++) {
		if (strcmp(argv[i], "-dir") == 0) {
			CONFIG::workdir = argv[++i];
			continue;
		}
	}


	//---------------------------------------Build the Universe-------------------------------------------

	Universe universe;

	universe.addMoleculeSpecies(
			new MoleculeInfo("Aspartate", 0, 890, 890 * 0.6, 0, 0));

	//---------------------------Create a stable one dimensional chemoattractant gradient -----------------

	double chemotaxis_gradient = 1e-4;
	CONFIG::diffusion = false; // constant linear gradient

	for (unsigned int x = 0; x != CONFIG::gridNumberX; ++x)
		for (unsigned int y = 0; y != CONFIG::gridNumberY; ++y)
			for (unsigned int z = 0; z != CONFIG::gridNumberZ; ++z) {
				Grid* thisGrid = CONFIG::universe->getGrid(x, y, z);
				thisGrid->setConc(0,
						(CONFIG::worldSizeY - y) * chemotaxis_gradient
								* CONFIG::gridSizeY);
			}

	//------------------------------------------Add Some Bacteria-----------------------------------------------

	std::size_t bacNum = 10000;
	unsigned int i = 0;
	while (i++ < bacNum) {
		myVector3d Position;
		Position.pos.x = 0.1 * CONFIG::worldSizeX
				+ 0.8 * CONFIG::worldSizeX * ((double) rand() / (RAND_MAX));
		Position.pos.y = 0.1 * CONFIG::worldSizeY
				+ 0.8 * CONFIG::worldSizeY * ((double) rand() / (RAND_MAX));
		Position.pos.z = 0.05 * CONFIG::worldSizeZ
				+ 0.1 * CONFIG::worldSizeZ * ((double) rand() / (RAND_MAX));
		Agent* newBac = new ChemotacticBacteria(Position,
				1.5 + 0.49 * ((double) rand() / (RAND_MAX)), 1.2, 15);   // velocity changed to 15 um/s
		universe.addAgent(newBac);
	}

	//-----------------------------------------Main Simulation Loop----------------------------------------------
	regularExporters::dump_Agent();
	regularExporters::dump_Con(0);

	while (CONFIG::time < CONFIG::simLength) {
		universe.evolute();
		cout << "Simulation time: " << CONFIG::time
				<< " seconds Agents number: " << universe.getTotalAgentNumber()
				<< endl;

		if (((int) CONFIG::time) % 1 == 0
				&& (CONFIG::time - (int) CONFIG::time) < CONFIG::timestep) {
			regularExporters::dump_Agent();
		}
	}

	time(&end);
	double dif = difftime(end, start);

	cout << "Simulation done. Elapsed time: " << dif << " secs." << endl;

	return 0;
}
