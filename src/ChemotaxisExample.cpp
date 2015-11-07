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
#include "configuration.h"
#include "universe.h"
#include "moleculeInfo.h"
#include "regularExporters.h"
#include "ChemotacticBacteria.h"
#include <ctime>
#include <cstring>
#include <cstdint>

using namespace std;
using namespace BNSim;

struct Parameters {
    uint32_t agent_count    = 3600; 

    /* world size in micrometers */
    uint32_t sizex          = 1000; 
    uint32_t sizey          = 3000;
    uint32_t sizez          = 1000;
    
    /* world resolution */
    uint32_t gridx          = 100;
    uint32_t gridy          = 300;
    uint32_t gridz          = 100;
    
    uint32_t species        = 1;
    uint32_t threads        = 8;
    double   timestep       = 0.1;
    double   simtime        = 3600;
};

Parameters params;

// Quorum sensing example

int main(int argn, char **argv) {

	time_t start, end;
	time(&start);

    // Init random number generator
    srand48(time (0));

	cout << "BNSim 2.0 -- Chamotaxis Example" << endl;

	CONFIG::workdir = "./";

	CONFIG::worldSizeX  = params.sizex; 
	CONFIG::worldSizeY  = params.sizey;
	CONFIG::worldSizeZ  = params.sizez;

	CONFIG::gridNumberX = params.gridx;
	CONFIG::gridNumberY = params.gridy;
	CONFIG::gridNumberZ = params.gridz;

	CONFIG::numberMoleculeSpecies = params.species;
	CONFIG::threadNumber = params.threads;
	CONFIG::timestep    = params.timestep;  
	CONFIG::simLength   = params.simtime;  

	Universe universe;

	universe.addMoleculeSpecies(
			new MoleculeInfo("Aspartate", 0, 890, 890 * 0.6, 0, 0));

	/* Add gradient */
	CONFIG::diffusion = false; // constant linear gradient

	for (uint32_t x = 0; x < params.gridx; x++) {
		for (uint32_t y = 0; y < params.gridy; y++) {
			for (uint32_t z = 0; z < params.gridz; z++) {
				
				Grid* thisGrid = CONFIG::universe->getGrid(x, y, z);
				
				if (y <= 1.0 / 3.0 * params.gridy) {    
				    /* No food */
				    thisGrid->setConc(0, 0.0);
				    				
				} else if (y <= 2.0 / 3.0 * params.gridy) {
				    /* High gradient 0 to 1 */
				    thisGrid->setConc(0, (double(y) / params.gridy - 1.0 / 3.0) *
				                         3.3 / params.gridz);
				    
				} else {
				    /* Food without gradient */
				    thisGrid->setConc(0, 1.0 / params.gridz);
				}
			}
		}
	}

	/* Add bacteria */
    for (uint32_t i = 0; i < params.agent_count; i++) {
        double x, y, z, bias, radius;
        
        radius  =  1.5 + 0.49 * drand48();
        bias    = (i * 3 / params.agent_count) / 3.0 + 0.1;
        x       = (0.25 + 0.5 * drand48()) * params.sizex;
        y       = (bias + 0.1 * drand48()) * params.sizey;
        z       = (0.25 + 0.5 * drand48()) * params.sizez;
          
		universe.addAgent(
		        new ChemotacticBacteria(myVector3d(x, y, z), radius, 1.2, 15));
    }

	/* Simulate */
	ofstream dump_file("bacteria.txt", ofstream::trunc);
	regularExporters::dump_Con(0, "concentration.txt");

    dump_file.precision(5);

	while (CONFIG::time < CONFIG::simLength) {
		universe.evolute();
		
		if (((int) CONFIG::time) % 1 == 0
				&& (CONFIG::time - (int) CONFIG::time) < CONFIG::timestep) {
			regularExporters::export_agent_position(dump_file);
			cout << "Simulation time: " << CONFIG::time << "s" << endl;
		}
	}

	dump_file.close();
	
	time(&end);
	cout << "Simulation done." << endl;
	cout << "Elapsed time: " << difftime(end, start) << " secs." << endl;

	return 0;
}
