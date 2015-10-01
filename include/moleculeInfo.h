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

#ifndef MOLECULEINFO_H_
#define MOLECULEINFO_H_

#include<string>
#include "common.h"

namespace BNSim {

class MoleculeInfo {
public:
	MoleculeInfo(const std::string& name, int index, double Dc_boundaryLayer,
			double Dc_biofilmLayer, double decay_boundaryLayer,
			double decay_biofilmLayer);
	~MoleculeInfo();
	const std::string& getName() const {
		return _name;
	}
	unsigned int getIndex() const {
		return _index;
	}
	double getDiffusionCoefficient(layerType type) const {
		if (type == biofilm)
			return _Dc_biofilmLayer;
		else
			return _Dc_boundaryLayer;   // bulk also use the boundary diffusion coefficient for averaging purpose
	}
	double getDecayRate(layerType type) const {
		if (type == biofilm)
			return _decay_biofilmLayer;
		else if (type == boundary)
			return _decay_boundaryLayer;
	}
private:
	std::string _name;
	unsigned int _index;
	double _Dc_boundaryLayer, _Dc_biofilmLayer; // diffusion coefficient in diffusion and biofilm layer, respectively
	double _decay_boundaryLayer, _decay_biofilmLayer;  // molecular degradation
};

} /* namespace BNSim */

#endif /* MOLECULEINFO_H_ */
