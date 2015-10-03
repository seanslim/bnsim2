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

#include "regularExporters.h"
#include <stdio.h>

using namespace std;

namespace BNSim {

void regularExporters::dump_Con(unsigned int chemicalIndex, const char file_name[]) {
	int clk = (int) (CONFIG::time);
	ofstream gradient(file_name);

	for (unsigned int i = 0; i != CONFIG::gridNumberX; ++i) {
		for (unsigned int j = 0; j != CONFIG::gridNumberY; ++j) {
			double total = 0;
			for (unsigned int k = 0; k != CONFIG::gridNumberZ; ++k) {

				Grid * check = CONFIG::universe->getGrid(i, j, k);
				total += check->getConc(chemicalIndex);
			}
			gradient << total << " ";
		}
		gradient << endl;
	}
	gradient.close();
}

void regularExporters::dump_QS_Status(ofstream& QSStatus) {

	int count = 0;
	double LuxR = 0, LuxC = 0, AI = 0;
	int clk = (int) (CONFIG::time);

	for (unsigned int i = 0; i != CONFIG::universe->getTotalAgentNumber();
			++i,++count) {
		Agent* agent = CONFIG::universe->getAgent(i);
		QSLux* qs = (QSLux*) agent->getRegulatoryNet("QSLux");
		LuxR += qs->getR1();
		LuxC += qs->getC1();
		AI += qs->getA1();
	}

	QSStatus<<clk<<","<<count<<","<<LuxR/count<<","<<LuxC/count<<","<<AI/count<<endl;
}

void regularExporters::dump_Metabolism_Status(ofstream& MetaStatus) {

	int count = 0;
	double substrate = 0, u=0;
	int clk = (int) (CONFIG::time);

	for (unsigned int i = 0; i != CONFIG::universe->getTotalAgentNumber();
			++i,++count) {
		Agent* agent = CONFIG::universe->getAgent(i);
		SimpleMetabolism* meta = (SimpleMetabolism*) agent->getRegulatoryNet("SimpleMetabolism");
		substrate += meta->getS();
		u += meta->getU();
	}

	MetaStatus<<clk<<","<<count<<","<<substrate/count<<","<<u/count<<endl;
}

void regularExporters::dump_Agent() {
	char file_name[100];
	int clk = (int) (CONFIG::time);
	string path = CONFIG::workdir + "/agent_%d.txt";
	sprintf(file_name, path.c_str(), clk);
	ofstream agents(file_name);

	cout << "dump bacteria spatial distribution at:" << file_name << endl;

	for (unsigned int i = 0; i != CONFIG::gridNumberX; ++i) {
		for (unsigned int j = 0; j != CONFIG::gridNumberY; ++j) {
			unsigned int total = 0;
			for (unsigned int k = 0; k != CONFIG::gridNumberZ; ++k) {
				total += CONFIG::universe->getGrid(i, j, k)->getAgentNumber();
			}
			agents << total << " ";
		}
		agents << endl;
	}
	agents.close();
}

void regularExporters::export_agent_position(ofstream& output_file) {
	
	for (int i = 0; i < CONFIG::universe->getTotalAgentNumber(); i++) {
		myVector3d abs_pos = CONFIG::universe->getAgent(i)->getabsPosition();
		
		output_file << abs_pos.pos.x << " " << abs_pos.pos.y << " " 
					<< abs_pos.pos.z << " ";
	}
	
	output_file << endl;
}

void regularExporters::dump_CUDA_Biofilm() {
	char file_name[100];
	int clk = (int) (CONFIG::time);
	string path = CONFIG::workdir + "/CUDA_Biofilm_%d.txt";
	sprintf(file_name, path.c_str(), clk);
	ofstream agents(file_name);

	cout << "dump bacteria spatial positions at:" << file_name << endl;

	for (unsigned int i = 0; i != CONFIG::universe->getTotalAgentNumber();
			++i) {
		Agent* agent = CONFIG::universe->getAgent(i);

		if (agent != NULL)
			agents << agent->getID() << " " << agent->getabsPosition().pos.x
					<< " " << agent->getabsPosition().pos.y << " "
					<< agent->getabsPosition().pos.z << " "
					<< agent->getTotalRadius() << endl;
	}
	agents.close();
}

} /* namespace BNSim */
