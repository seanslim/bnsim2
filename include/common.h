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

#ifndef COMMON_H_
#define COMMON_H_

#include<cmath>
#include<string>

namespace BNSim {

struct intVector3d {
	int x, y, z;
};

struct doubleVector3d {
	double x, y, z;
};

struct massInfo {
	std::string name;
	double mass;
	double density;
};

enum layerType {
	bulk, boundary, biofilm
};

enum EPSType {
	DNA, protein
};

class myVector3d {
public:
	myVector3d();
	myVector3d(double x, double y, double z);
	myVector3d(const myVector3d& v);
	~myVector3d();
	myVector3d& operator=(myVector3d rhs) {
		pos.x = rhs.pos.x;
		pos.y = rhs.pos.y;
		pos.z = rhs.pos.z;
		return *this;
	}
	;
	void normalize();
	double length();
	double lengthSquared();
	void add(const myVector3d& v);
	void sub(const myVector3d& v);
	void negate(const myVector3d& v);
	void scale(double s);
	void scale(double s, const myVector3d& v);
	void cross(const myVector3d& v1, const myVector3d& v2);
	doubleVector3d pos;
};

class myMatrix3d {
public:
	myMatrix3d();
	~myMatrix3d();
	void transform(myVector3d& t);
	void setFromAxisAngle(double x, double y, double z, double angle);

private:
	double m00, m01, m02, m10, m11, m12, m20, m21, m22;
};

template<typename T>
class BNSimVector {
private:
	T* elem;
	unsigned long size;
	unsigned long count;
	unsigned long capacity;
public:
	unsigned long getSize() {
		return size;
	}
	unsigned long getCount() {
		return count;
	}
	BNSimVector() {
		elem = new T[10];
		count = size = 0;
		capacity = 10;
	}
	BNSimVector(unsigned long s) {
		elem = new T[s];
		count = size = 0;
		capacity = s;
	}
	~BNSimVector() {
		delete[] elem;
	}
	void setCapacity(unsigned long s) {
		elem = new T[s];
		count = size = 0;
		capacity = s;
	}
	T& operator[](int i) {
		return elem[i];
	}
	void add(T& newElem) {
		elem[size++] = newElem;
		count++;
		if (size / capacity > 0.8) {
			capacity *= 2;
			T* temp = new T[capacity];
			int j = 0;
			for (int i = 0; i != size; ++i) {
				if (elem[i] != NULL)
					temp[j++] = elem[i];
			}
			delete[] elem;
			elem = temp;
			count = size = j;
		}
	}
	;
	void remove(T& newElem) {
		for (int i = 0; i != size; ++i)
			if (newElem == elem[i]) {
				elem[i] = NULL;
				count--;
				return;
			}
	}
	;
};

} /* namespace BNSim */

#endif /* COMMON_H_ */
