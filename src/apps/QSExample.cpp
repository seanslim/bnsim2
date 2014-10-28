///**
// * BNSim is an open-source, parallel, stochastic, and multi-scale modeling
// * platform which integrates various simulation algorithms, together with
// * chemotaixs, quorum sensing, and biofilm models in a 3D environment.
// *
// * BNSim is developed by CMU SLD Group, and released under GPLv2 license
// * Please check http://www.ece.cmu.edu/~sld/ for more information
// *
// * BNSim is developed using C++ and pthread under Linux and Mac OS
// *
// * If you use it or part of it for your work, please cite
// * Wei, Guopeng, Paul Bogdan, and Radu Marculescu. "Efficient Modeling and
// * Simulation of Bacteria-Based Nanonetworks with BNSim." Selected Areas in
// * Communications, IEEE Journal on 31.12 (2013): 868-878.
// *
// * @author      Guopeng (Daniel) Wei  1@weiguopeng.com
// * @version     2.0
// * @since       1.0
// */
//
//#include <iostream>
//#include "../core/configuration.h"
//#include "../core/headers/universe.h"
//#include "../datastructures/headers/moleculeInfo.h"
//#include "../exporters/regularExporters.h"
//#include "../agents/headers/QSBacteria.h"
//#include <ctime>
//
//using namespace std;
//using namespace BNSim;
//
//// Quorum sensing example
//
//int main(int argn, char **argv) {
//
//	time_t start, end;
//	time(&start);
//
//	cout << "BNSim 2.0 -- Quorum Sensing Example"
//			<< endl;
//
//	//---------------------------------------Setup Critical Parameters-----------------------------------
//
//	CONFIG::workdir = "/home/guopeng/data/";
//
//	CONFIG::worldSizeX = 500; // 0.5 cm
//	CONFIG::worldSizeY = 500;
//	CONFIG::worldSizeZ = 500;
//	CONFIG::boundaryLayerThickness = 100;
//
//	CONFIG::gridNumberX = 25;
//	CONFIG::gridNumberY = 25;
//	CONFIG::gridNumberZ = 25;
//
//	CONFIG::numberMoleculeSpecies = 2;
//	CONFIG::threadNumber = 8;
//	CONFIG::timestep = 0.4;  // 0.4 seconds
//	CONFIG::simLength = 100000;  // in seconds
//
//    	for (int i = 1; i < argn; i++) {
//    		if (strcmp(argv[i], "-dir") == 0) {
//    			CONFIG::workdir = argv[++i];
//    			continue;
//    		}
//    	}
//
//	CONFIG::boundaryLayerThickness = 100; // this will update the height of layers automatically
//
//	//---------------------------------------Build the Universe-------------------------------------------
//
//	Universe universe;
//
//	universe.addMoleculeSpecies(
//			new MoleculeInfo("substrate", 0, 890, 890 * 0.6, 0, 0));
//	universe.addMoleculeSpecies(
//			new MoleculeInfo("AI", 1, 890, 890 * 0.6, 0.01, 0.01));
//
//	//---------------------------Create a constant substrate concentration in the bulk layer -----------------
//
//    double constantSubstrate = 1;
//    
//	for (unsigned int x = 0; x != CONFIG::gridNumberX; ++x)
//		for (unsigned int y = 0; y != CONFIG::gridNumberY; ++y)
//			for (unsigned int z = 0; z != CONFIG::gridNumberZ; ++z) {
//				Grid* thisGrid = CONFIG::universe->getGrid(x, y, z);
//                if(thisGrid->getLayerType()==bulk)
//                    thisGrid->setConc(0, 1);
//			}
//
//	//------------------------------------------Add Some Bacteria-----------------------------------------------
//
//	std::size_t bacNum = 100;
//	unsigned int i = 0;
//	while (i++ < bacNum) {
//		myVector3d Position;
//		Position.pos.x = 0.1 * CONFIG::worldSizeX
//				+ 0.8 * CONFIG::worldSizeX * ((double) rand() / (RAND_MAX));
//		Position.pos.y = 0.1 * CONFIG::worldSizeY
//				+ 0.8 * CONFIG::worldSizeY * ((double) rand() / (RAND_MAX));
//		Position.pos.z = 0.05 * CONFIG::worldSizeZ
//				+ 0.1 * CONFIG::worldSizeZ * ((double) rand() / (RAND_MAX));
//		Agent* newBac = new QSBacteria(Position,
//				1.5 + 0.49 * ((double) rand() / (RAND_MAX)), 1.2, 2, 2, 0.2);
//		universe.addAgent(newBac);
//	}
//
//	//-----------------------------------------Main Simulation Loop----------------------------------------------
//	regularExporters::dump_Agent();
//	regularExporters::dump_Con(0);
//
//	string path = CONFIG::workdir + "/agent_QS.txt";
//	ofstream QSStatus(path.c_str());
//
//	path = CONFIG::workdir + "/agent_META.txt";
//	ofstream METAStatus(path.c_str());
//
//	while (CONFIG::time < CONFIG::simLength) {
//		universe.evolute();
//		cout << "Simulation time: " << CONFIG::time
//				<< " seconds Agents number: " << universe.getTotalAgentNumber()
//				<< endl;
//
//		if (((int) CONFIG::time) % 20 == 0
//				&& (CONFIG::time - (int) CONFIG::time) < CONFIG::timestep) {
//			regularExporters::dump_Con(1);
//			regularExporters::dump_Con(0);
//			regularExporters::dump_Agent();
////			regularExporters::dump_CUDA_Biofilm();
//			regularExporters::dump_QS_Status(QSStatus);
//			regularExporters::dump_Metabolism_Status(METAStatus);
//		}
//	}
//
//	time(&end);
//	double dif = difftime(end, start);
//
//	cout << "Simulation done. Elapsed time: " << dif << " secs." << endl;
//
//	return 0;
//}
