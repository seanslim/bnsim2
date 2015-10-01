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

#ifndef REGULATORYNET_H_
#define REGULATORYNET_H_


#include "configuration.h"
#include <string>

namespace BNSim {

class Agent;
class SimpleMetabolism;

class RegulatoryNet {
public:
	RegulatoryNet(Agent* host);
	virtual ~RegulatoryNet();
	virtual void update() = 0;
	Agent* getHost() { return _host;}
    std::string& getName() {return name;}
protected:
	Agent* _host;
    std::string name;
};

class QSLux: public RegulatoryNet {
public:
	QSLux(Agent* host);
	virtual ~QSLux();
	virtual void update();
	bool isActivated() { return activated; }
	double getA1() { return A1;}
	double getR1() { return R1;}
	double getA2() { return A2;}
	double getR2() { return R2;}
    double getC1() { return C1;}
	double getC2() { return C2;}

private:
	// las
	double A1;
	double C1;
	double S;
	double R1;
	double R1A1;
	double R1QSI1;

	// rhl
	double A2;
	double C2;
	double R2;
	double R2A2;

	// QSI
	double qsi1;
	double qsi2;
	double qsi3;

	// define rate constants
	double kA1;
	double kA2;
	double kR1;
	double kR2;
	double kS;
	double CR1;
	double CR2;
	double CA1;
	double CA2;
	double CS;

	// define binding constants
	double KL1;
	double KL2;
	double KL3;
	double KL4;
	double KS1;
	double KS;
	double KR1;
	double k[13];

	bool activated;

    double randomUniform();
};

class SimpleMetabolism: public RegulatoryNet {
private:
	double _m,_u_PG,_u_PGMax;
	double _Y;
	double _u,_u_max; // substrate
	double _Ks,_Kc;
	QSLux *_qs;
	double _s;
	void grow();
    double randomUniform();
public:
	SimpleMetabolism(Agent* host);
	virtual ~SimpleMetabolism();
	virtual void update();
	double getS() const {	return _s;}
	double getU() const {	return _u;	}
	void setU(double u) {	_u = u;	}
	double getMaintenance() const {	return _m;	}
	void setMaintenance(double maintenance) {	_m = maintenance;	}
	double getuMax() const {	return _u_max;	}
	void setuMax(double max) {	_u_max = max;	}
	double getPg() const {	return _u_PG;	}
	void setPg(double pg) {	_u_PG = pg;	}
	double getY() const {	return _Y;	}
	void setY(double y) {	_Y = y;	}
	const QSLux* getQs() const { return _qs;	}
	void setQs(QSLux* qs) {	_qs = qs;	}
	double getKc() const {	return _Kc;	}
	void setKc(double kc) {	_Kc = kc;	}
	double getKs() const {	return _Ks;	}
	void setKs(double ks) {	_Ks = ks;}
	double getPgMax() const {	return _u_PGMax;	}
	void setPgMax(double pgMax) {	_u_PGMax = pgMax;	}
};

class ChemotaxisSystem : public RegulatoryNet{
public:
    ChemotaxisSystem(Agent* host);
    virtual void update();
    void rotationalDiffusion();
    double tumbleAngle();
    void tumble();
    double sampleGamma(double k, double theta);
    void rotatePerp(myVector3d& v, double theta);
    void rotate(myVector3d& v, myVector3d axis, double theta);
    bool isCCW() {return CCW;}
    void setDirection(const myVector3d& dir) { direction = dir;};
    const myVector3d& getDirection() { return direction;}
private:
    double Y;
    double Ya;
    double Ye;
    double H;
    double inverseTau;
    double kr;
    double kb;
    double m;
    double m0;
    double aspcon;
    double activity;
    double alpha;
    double KA;
    double KI;
    int N_tar;
    double Kd;
    double g0;
    double g1;
    double w;
    bool CCW;
    myVector3d direction;
    double randomUniform();

    // motor
    double tau;
    double runtime;
    double tumbletime;
    double kPlus();
    double kMinus();
};

} /* namespace BNSim */

#endif /* REGULATORYNET_H_ */
