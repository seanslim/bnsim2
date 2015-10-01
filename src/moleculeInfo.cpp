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

#include "moleculeInfo.h"

namespace BNSim {

MoleculeInfo::MoleculeInfo(const std::string& name, int index, double Dc_boundaryLayer, double Dc_biofilmLayer, double decay_boundaryLayer,
		double decay_biofilmLayer):_name(name),_index(index),_Dc_boundaryLayer(Dc_boundaryLayer),_Dc_biofilmLayer(Dc_biofilmLayer),
		_decay_boundaryLayer(decay_boundaryLayer),_decay_biofilmLayer(decay_biofilmLayer)
{

}

MoleculeInfo::~MoleculeInfo() {
	// TODO Auto-generated destructor stub
}

} /* namespace BNSim */
