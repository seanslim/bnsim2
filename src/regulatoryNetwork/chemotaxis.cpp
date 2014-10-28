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

#include "headers/regulatoryNet.h"

namespace BNSim {

ChemotaxisSystem::ChemotaxisSystem(Agent * host) :
		RegulatoryNet(host) {
	CCW = false;
	Y = 0;
	Ya = 5;
	Ye = 3; //um
	H = 10;
	inverseTau = 5;
	kr = 0.2;
	kb = 0.2;
	m = 2;
	m0 = 1;
	alpha = 1.7;
	N_tar = 6;
	aspcon = 0;
	activity = 0.5;
	KA = 3;
	KI = 0.0182;
	Kd = 3.06;
	g0 = 40;
	g1 = 40;
	w = 1.3;
	tau = 0;
	runtime = 0;
	tumbletime = 1;
    myVector3d random(0.5-randomUniform(),0.5-randomUniform(),0.5-randomUniform());
            random.normalize();
    direction = random;
}

double ChemotaxisSystem::randomUniform() {
	return ((double) rand() / (RAND_MAX));
}

void ChemotaxisSystem::update() {

	// Read chemoattract concentration from environment

	int Aspartate = -1;
	Aspartate =
			CONFIG::universe->getMoleculeInfo("Aspartate") == NULL ?
					-1 :
					CONFIG::universe->getMoleculeInfo("Aspartate")->getIndex();

	intVector3d GPos = getHost()->getGridPos();
	Grid* myGrid = CONFIG::universe->getGrid(GPos.x, GPos.y, GPos.z);

	if (Aspartate != -1)
		aspcon = myGrid->getConc(Aspartate);

	m = m
			+ CONFIG::timestep
					* (kr * CONFIG::timestep * (1 - activity)
							- kb * CONFIG::timestep * activity);

	double fm = alpha * (m0 - m);
	double epsilon = fm - log((1 + aspcon / KA) / (1 + aspcon / KI));
	activity = 1 / (1 + exp(N_tar * epsilon));

	Y = Ya * activity;

	if (CCW && randomUniform() < kMinus()) {
		CCW = false;
	} else if (!CCW && randomUniform() < kPlus()) {
		CCW = true;
	}
    
    if(!CCW)
        tumble();
    else
        rotationalDiffusion();
}

/*
 Sneddon, Michael W., William Pontius, and Thierry Emonet. "Stochastic coordination of multiple actuators reduces latency and improves chemotactic response in bacteria." Proceedings of the National Academy of Sciences 109.3 (2012): 805-810.
 */

// rate from CW to CCW
double ChemotaxisSystem::kPlus() {
	double r = w * exp((g0 / 4.0) - (g0 / 2.0) * (Y / (Kd + Y)));
	return r;
}

// rate from CCW to CW
double ChemotaxisSystem::kMinus() {
	double r = w * exp(-(g1 / 4.0) + (g1 / 2.0) * (Y / (Kd + Y)));
	return r;
}

/**
 * Causes the cell to rotate such that Var(theta(dt)) = 4*D*dt.
 */
void ChemotaxisSystem::rotationalDiffusion() {
	double Dr = 0.28;  // rotational diffusion constant rad^2/s
	double dTheta = randomUniform()*sqrt(2*Dr*CONFIG::timestep);
	rotatePerp(direction, dTheta);
}

/**
 * Return a tumble angle in radians distributed according to Fig. 3, 'Chemotaxis
 * in Escherichia Coli', Berg et al. (claim from 'AgentCell: a digital single-cell
 * assay for bacterial chemotaxis', Emonet et al.).
 */
double ChemotaxisSystem::tumbleAngle() {
	double tumbleShape = 4;
	double tumbleScale = 18.32;
	double tumbleLocation = -4.60;

	double tumbleAngle;
	do {
		tumbleAngle = sampleGamma(tumbleShape, tumbleScale) + tumbleLocation;
	} while (tumbleAngle > 180);

	return (tumbleAngle*3.14/360.0);
}

void ChemotaxisSystem::tumble() {
	rotatePerp(direction, tumbleAngle());
}

double ChemotaxisSystem::sampleGamma(double k, double theta) {
	bool accept = false;
	if (k < 1) {
		// Weibull algorithm
		double c = (1 / k);
		double d = ((1 - k) * pow(k, (k / (1 - k))));
		double u, v, z, e, x;
		do {
			u = randomUniform();	v = randomUniform();
			z = -log(u);		e = -log(v);
			x = pow(z, c);
			if ((z + e) >= (d + x)) {
				accept = true;
			}
		} while (!accept);
		return (x * theta);
	} else {
		// Cheng's algorithm
		double b = (k - log(4));
		double c = (k + sqrt(2 * k - 1));
		double lam = sqrt(2 * k - 1);
		double cheng = (1 + log(4.5));
		double u, v, x, y, z, r;
		do {
			u = randomUniform();		v = randomUniform();
			y = ((1 / lam) * log(v / (1 - v)));
			x = (k * exp(y));
			z = (u * v * v);
			r = (b + (c * y) - x);
			if ((r >= ((4.5 * z) - cheng)) ||
					(r >= log(z))) {
				accept = true;
			}
		} while (!accept);
		return (x * theta);
	}
}

/**
 * Rotates the vector v by an angle theta in a random direction perpendicular to v.
 * From BSim, Gorochowski, Thomas E., et al. "BSim: an agent-based tool for modeling
 * bacterial populations in systems and synthetic biology." PloS one 7.8 (2012): e42790.
 */
void ChemotaxisSystem::rotatePerp(myVector3d& v, double theta) {
	/* Obtain a random direction perpendicular to v */
	myVector3d random(0.5-randomUniform(),0.5-randomUniform(),0.5-randomUniform());
	myVector3d randomPerp;
	randomPerp.cross(v, random);
	rotate(v, randomPerp, theta);
	// Ensure the vector is unit
	v.normalize();
}

/**
 * Rotates the vector v towards the specified axis by an angle theta.
 * From BSim, Gorochowski, Thomas E., et al. "BSim: an agent-based tool for modeling
 * bacterial populations in systems and synthetic biology." PloS one 7.8 (2012): e42790.
 */
void ChemotaxisSystem::rotate(myVector3d& v, myVector3d axis, double theta) {
	/* Generate the rotation matrix for rotating about the axis by an angle theta */
	myMatrix3d  r;
	r.setFromAxisAngle(axis.pos.x, axis.pos.y, axis.pos.z, theta);
	/* Apply the rotation */
	r.transform(v);
}

} /* namespace BNSim */
