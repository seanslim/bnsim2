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

#ifndef EPS_H_
#define EPS_H_

#include "agent.h"


namespace BNSim {

class EPS: public Agent {
public:
	EPS(const myVector3d& absPosition,EPSType type, double radius, double shovek);
	virtual ~EPS();
	virtual void update();
	EPSType getType();
	double getEPSMass();
private:
	EPSType _type;
};

} /* namespace BNSim */

#endif /* PAQS_H_ */
