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

#include "regulatoryNet.h"

namespace BNSim {

QSLux::QSLux(Agent* host):RegulatoryNet(host){
	A1 = 0;
	C1 = 0;
	S = 0;
	R1 = 0;
	R1A1 = 0;

	R1QSI1 = 0;

	// rhl
	A2 = 0;
	C2 = 0;
	R2 = 0;
	R2A2 = 0;

	// QSI

	qsi1 = 0;
	qsi2 = 0;
	qsi3 = 0;

	// define rate constants
	kA1 = 2e-3;
	kA2 = 2e-3;
	kR1 = 2e-3;
	kR2 = 2e-3;
	kS = 2e-3;
	CA1 = 1e-4;
	CA2 = 1e-4;
	CR1 =  1e-4;
	CR2 =  1e-4;
	CS =  1e-4;

	KL1 = 1e-9;
	KL2 = 1e-3;
	KL3 = 1e-6;
	KL4 = 1e-6;
	KS1 = 1e-6;
	KS = 1e-6;
	KR1 = 1e-9;

	int i = 0;
	while(i++!=12) {
		k[i] = 0.1;
	}

	k[0] = k[7] = 0.01;
	k[3] = k[10] = k[6] = 0.01;

	activated = false;
    name = "QSLux";
}

QSLux::~QSLux() {
}

double QSLux::randomUniform() {
	return ((double) rand() / (RAND_MAX));
}

void QSLux::update() {

	//double deltaA1 = (CA1 + kA1*C1*KS1/((KL1+C1)*(KS1+S)) - k[0]*A1 - k[1]*A1*R1 + k[2]*R1A1)*timestep;

	int QSI1 = -1, QSI2 = -1, QSI3 = -1;
	QSI1 = CONFIG::universe->getMoleculeInfo("QSI1")==NULL?-1:CONFIG::universe->getMoleculeInfo("QSI1")->getIndex();
	QSI2 = CONFIG::universe->getMoleculeInfo("QSI2")==NULL?-1:CONFIG::universe->getMoleculeInfo("QSI2")->getIndex();
	QSI3 = CONFIG::universe->getMoleculeInfo("QSI3")==NULL?-1:CONFIG::universe->getMoleculeInfo("QSI3")->getIndex();
	// now we disable the S feedback

	intVector3d GPos = getHost()->getGridPos();

	Grid* myGrid = CONFIG::universe->getGrid(GPos.x,GPos.y,GPos.z);
	double timestep = CONFIG::timestep;

	// sense AI

	double myqsi1 = 0, myqsi2 = 0, myqsi3 = 0;
	double kh = 0.3, g= 2;

	double deltaA1 = (CA1 + kA1*C1/(KL1+C1)- k[0]*A1 - k[1]*A1*R1 + k[2]*R1A1 - A1*myqsi2*kh/(1+g*myqsi2))*timestep;
	double deltaR1 = ((CR1 + kR1*C1/(KR1+C1))/(1+myqsi3/KL1) - k[3]*R1 - k[1]*A1*R1 - k[1]*myqsi1*R1 + k[2]*R1A1)*timestep;
	double deltaR1A1 = (k[1]*A1*R1 - k[2]*R1A1 -2*k[4]*R1A1*R1A1 + 2*k[5]*C1)*timestep;
	double deltaR1QSI1 =  (k[1]*myqsi1*R1 - k[2]*R1QSI1)*timestep;
	double deltaQSI1 =  (- k[1]*myqsi1*R1 + k[2]*R1QSI1)*timestep;
	double deltaC1 = (k[4]*R1A1*R1A1 - k[5]*C1)*timestep;
	//double deltaS = (CS + kS*C1/(KS +C1) - k[6]*S)*timestep;

	A1 += (deltaA1 + 0.1*deltaA1*(randomUniform()-0.5));
	C1 += (deltaC1 + 0.1*deltaC1*(randomUniform()-0.5));
	R1 += (deltaR1 + 0.1*deltaR1*(randomUniform()-0.5));
	R1A1 += (deltaR1A1 + 0.1*deltaR1A1*(randomUniform()-0.5));

	if(QSI1!=-1)  {
		qsi1 += deltaQSI1;
		R1QSI1 += deltaR1QSI1;
	}
	//S += (deltaS + 0.1*deltaS*(randomUniform()-0.5));;

//	double deltaA2 = (CA2 + kA2*C2/((KL2+C2)) - k[7]*A2 - k[8]*A2*R2 + k[9]*R2A2)*timestep;
//	double deltaR2 = (CR2 + kR2*C2*C1/((KL3+C2)*(KL4+C1)) - k[10]*R2 - k[8]*A2*R2 + k[9]*R2A2)*timestep;
//	double deltaR2A2 = (k[8]*A2*R2 - k[9]*R2A2 -2*k[11]*R2A2*R2A2 + 2*k[12]*C2)*timestep;
//	double deltaC2 = (k[11]*R2A2*R2A2 - k[12]*C2)*timestep;
//
//	A2 += (deltaA2 + 0.5*deltaA2*(randomUniform()-0.5));
//	C2 += (deltaC2 + 0.5*deltaC2*(randomUniform()-0.5));
//	R2 += (deltaR2 + 0.5*deltaR2*(randomUniform()-0.5));
//	R2A2 += (deltaR2A2 + 0.5*deltaR2A2*(randomUniform()-0.5));

	//	if(Sim::param.networkExport && C1>Sim::param.QST && !this->activated) {
	//		//   Sim::univ->BacNetwork->addNode(this->host->getID());
	//		activated = true;
	//	}
	//	else if(Sim::param.networkExport && C1<Sim::param.QST && this->activated) {
	//		//    Sim::univ->BacNetwork->addNode(this->host->getID());
	//		activated = false;
	//	}
	double agentVolume = getHost()->getTotalRadius()
			* getHost()->getTotalRadius() * getHost()->getTotalRadius();

	int AIIndex = -1;
	AIIndex = CONFIG::universe->getMoleculeInfo("AI")==NULL?-1:CONFIG::universe->getMoleculeInfo("AI")->getIndex();

	// Now let's exchange molecules with extracellular space
	double D = 0.5, ex = 0, delta = 0;

	if(AIIndex!=-1) {
		ex = myGrid->getConc(AIIndex);
		double deltaCon = (A1 - ex) * D * timestep;
		double deltaMass = deltaCon*agentVolume;
		A1 = (A1 * agentVolume-deltaMass)/agentVolume;
		myGrid->deltaChemical(AIIndex, deltaMass);
	}
}

} /* namespace BNSim */
