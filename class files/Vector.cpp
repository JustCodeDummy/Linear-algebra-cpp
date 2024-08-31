//
// Created by fsociety on 18/09/2023.
//

#include "../Header files/Vector.h"
#include <cmath>
#include <vector>
#include "../Header files/matrix.h"

//TODO implement Vector exceptions
using u_int = unsigned int;
using u_long = unsigned long;

Vector::Vector(const std::vector<double>& v) {
	this-> _vector = v;
}

Vector::Vector(int n, double element) {
	for (u_int i = 0; i<n; i++){
		this->_vector.push_back(element);
	}
}

void Vector::normalize() {
	double scalar = 1/length();
	for (double & i: this->_vector){
		i*=scalar;
	}
}

Vector Vector::crossProduct(const Vector &rhs) {
	std::vector<double> a(3, 0.0);
	a[0] = rhs[2] * this->_vector[1] - rhs[1] * this->_vector[2];
	a[1] = this->_vector[2] * rhs[0] - this->_vector[0] * rhs[2];
	a[3] = this->_vector[0] * rhs[1] - this->_vector[1] * rhs[0];
	return Vector(a);
}

double Vector::scalarTripleProduct(Vector v1, const Vector& v2) const {
	return *this * v1.crossProduct(v2);
}

Vector Vector::vectorTripleProduct(Vector& v1, const Vector &v2) {
	Vector a(3, 0.0);
	return a.crossProduct(v1.crossProduct(v2));
}

Vector Vector::projection( Vector v, double angle) const {
	v.normalize();
	v.scale(std::cos(angle) * this->length());
	return v;
}

Vector Vector::projection(const Vector &v1) const {
	double a = v1 * (*this);
	Vector _b(v1.size());
	_b.scale(a/v1.length());
	return _b;
}

double Vector::angle(const Vector &rhs, bool degrees) const {
	double p = (*this * rhs) / ((rhs.length()) * (*this).length());
	double r = std::acos(p);
	if (degrees){
		return r * (180 / M_PI);
	}
	return r;
}

void Vector::scale(double scalar) {
	for (double & i : this->_vector){
		i *= scalar;
	}
}

double Vector::length() const {
	double val = 0.0;
	for (double v : this->_vector){
		val += pow(v,2);
	}
	return sqrt(val);
}

bool Vector::isNull() const {
	return this->size() == 0.0;
}

Vector Vector::operator+(const Vector &rhs) {
	Vector v = *this;
	v += rhs;
	return v;
}

void Vector::operator*=(double scalar) {
	this->scale(scalar);
}

Vector Vector::operator-(const Vector &rhs) {
	Vector v = *this;
	v -= rhs;
	return v;
}

double Vector::operator*(const Vector &rhs) const {
	double value = 0.0;
	for (int i = 0; i<rhs.size(); i++){
		value += this->_vector[i] * rhs[i];
	}
	return value;
}

double& Vector::operator[](int index) {
	return this->_vector[index];
}
const double& Vector::operator[](int index) const{
	return this->_vector[index];
}

Vector Vector::operator*( double scalar) {
	Vector v = Vector(this->_vector);
	v.scale(scalar);
	return v;
}

void Vector::operator-=(const Vector& rhs){
	for (int i=0; i<rhs.size(); i++){
		this->_vector[i] -= rhs[i];
	}
}
void Vector::operator+=(const Vector& rhs){
	for (int i=0; i<rhs.size(); i++){
		this->_vector[i] += rhs[i];
	}
}

bool Vector::isOrthogonal(const Vector &v) const {
	return v * (*this) == 0;
}

bool Vector::isPerpendicular(const Vector &v) const {
	return (*this).isOrthogonal(v);
}

bool Vector::isCoplanar(const Vector &v) const {
	return (*this).copy().crossProduct(v).length() == 0;
}

bool Vector::isParallel(const Vector &v) const {
	return (*this).isCoplanar(v);
}

Matrix Vector::asMatrix() const {
	std::vector<Vector> v;
	v.push_back(*this);
	return Matrix(v);
}
void Vector::rotateX(double angle) {
	if (this->size() == 2) {
		std::vector<double> _t1 = {{std::cos(angle), -std::sin(angle)}};
		Vector t1(_t1);
		std::vector<double> _t2 = {std::sin(angle), std::cos(angle)};
		Vector t2(_t2);
		Matrix rotation({t1, t2});
		this->_vector = rotation.vectorMatrixMultiplication(Vector(this->_vector)).asVector();
	}
	if (this->size() == 3){
		Vector t3({1, 0, 0});
		Vector t1({0, std::cos(angle), -std::sin(angle)});
		Vector t2({0, std::sin(angle), std::cos(angle)});

		Matrix rotation({t3, t1, t2 });
		this->_vector = rotation.vectorMatrixMultiplication(Vector(this->_vector)).asVector();
	}

}

void Vector::rotateY(double angle) {
	if (this->size() == 3){
		Vector t3({0, 1, 0});
		Vector t1({std::cos(angle), 0 ,std::sin(angle)});
		Vector t2({ -std::sin(angle),0 ,std::cos(angle)});
		Matrix rotation({ t1, t3, t2 });
		this->_vector = rotation.vectorMatrixMultiplication(Vector(this->_vector)).asVector();
	}
}

void Vector::rotateZ(double angle) {
	if (this->size() == 3){
		Vector t3({1, 0, 0});
		Vector t1({std::cos(angle), -std::sin(angle), 0});
		Vector t2({std::sin(angle), std::cos(angle), 0});
		Matrix rotation({t1, t2, t3 });
		this->_vector = rotation.vectorMatrixMultiplication(Vector(this->_vector)).asVector();
	}
}

std::string Vector::toString() {
	std::string s = "( ";
	for (double value: this->_vector){
		s+= value;
		s+= " ";
	}
	s+= ")";

	return s;
}












