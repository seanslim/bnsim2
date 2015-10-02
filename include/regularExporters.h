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

#ifndef REGULAREXPORTERS_H_
#define REGULAREXPORTERS_H_

#include<iostream>
#include<fstream>
#include<string>
#include"configuration.h"

namespace BNSim {

class regularExporters {
public:
	static void dump_Con(unsigned int chemicalIndex);
	static void dump_Agent();
	static void dump_CUDA_Biofilm();
	static void dump_QS_Status(std::ofstream& QSStatus);
	static void dump_Metabolism_Status(std::ofstream& MetaStatus);
	static void export_agent_position(std::ofstream& output_file);
};

} /* namespace BNSim */

#endif /* REGULAREXPORTERS_H_ */
