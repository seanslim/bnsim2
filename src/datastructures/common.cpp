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

#include "headers/common.h"

namespace BNSim {

myMatrix3d::myMatrix3d() {
	m00 = m01= m02= m10= m11= m12= m20= m21= m22=0;
}

myMatrix3d::~myMatrix3d() {}

void myMatrix3d::transform(myVector3d& t) {

	double newx = m00*t.pos.x + m01*t.pos.y + m02*t.pos.z;
	double newy = m10*t.pos.x + m11*t.pos.y + m12*t.pos.z;
	double newz = m20*t.pos.x + m21*t.pos.y + m22*t.pos.z;

	t.pos.x = newx;	t.pos.y = newy; t.pos.z = newz;
}

void myMatrix3d::setFromAxisAngle(double x, double y, double z, double angle) {
	// Taken from Rick's which is taken from Wertz. pg. 412
	// Bug Fixed and changed into right-handed by hiranabe
	double n = sqrt(x*x + y*y + z*z);
	// zero-div may occur
	n = 1/n;
	x *= n;
	y *= n;
	z *= n;
	double c = cos(angle);
	double s = sin(angle);
	double omc = 1.0 - c;
	m00 = c + x*x*omc;
	m11 = c + y*y*omc;
	m22 = c + z*z*omc;

	double tmp1 = x*y*omc;
	double tmp2 = z*s;
	m01 = tmp1 - tmp2;
	m10 = tmp1 + tmp2;

	tmp1 = x*z*omc;
	tmp2 = y*s;
	m02 = tmp1 + tmp2;
	m20 = tmp1 - tmp2;

	tmp1 = y*z*omc;
	tmp2 = x*s;
	m12 = tmp1 - tmp2;
	m21 = tmp1 + tmp2;
}

myVector3d::myVector3d () {
	this->pos.x = 0;
	this->pos.y = 0;
	this->pos.z = 0;
}

myVector3d::myVector3d (double x, double y, double z) {
	this->pos.x = x;
	this->pos.y = y;
	this->pos.z = z;
}

myVector3d::myVector3d(const myVector3d& v) {
	this->pos.x = v.pos.x;
	this->pos.y = v.pos.y;
	this->pos.z = v.pos.z;
}

myVector3d::~myVector3d() {
}

void myVector3d::normalize() {
	double d = length();
	// zero-div may occur.
	pos.x /= d;
	pos.y /= d;
	pos.z /= d;

	if(length()<0.9)
		pos.x=2;
}

double myVector3d::lengthSquared() {
	return pos.x*pos.x + pos.y*pos.y + pos.z*pos.z;
}

double myVector3d::length() {
	return sqrt(lengthSquared());
}

void myVector3d::add(const myVector3d& v) {
	pos.x += v.pos.x;
	pos.y += v.pos.y;
	pos.z += v.pos.z;
}

void myVector3d::sub(const myVector3d& v) {
	pos.x -= v.pos.x;
	pos.y -= v.pos.y;
	pos.z -= v.pos.z;
}

void myVector3d::negate(const myVector3d& v) {
	pos.x = -v.pos.x;
	pos.y = -v.pos.y;
	pos.z = -v.pos.z;
}

void myVector3d::scale(double s, const myVector3d& v) {
	pos.x = s*v.pos.x;
	pos.y = s*v.pos.y;
	pos.z = s*v.pos.z;
}

void myVector3d::scale(double s) {
	pos.x *= s;
	pos.y *= s;
	pos.z *= s;
}

void myVector3d::cross(const myVector3d& v1, const myVector3d& v2) {
	// store on stack once for aliasing-safty
	// i.e. safe when a.cross(a, b)

	double newx = v1.pos.y*v2.pos.z - v1.pos.z*v2.pos.y;
	double newy = v1.pos.z*v2.pos.x - v1.pos.x*v2.pos.z;
	double newz = v1.pos.x*v2.pos.y - v1.pos.y*v2.pos.x;
	pos.x = newx; pos.y = newy; pos.z = newz;
}
    
} /* namespace BNSim */
