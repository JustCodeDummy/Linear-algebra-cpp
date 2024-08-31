// Vector.h
// Created by fsociety on 18/09/2023.
//
#include "matrix.h"
class Matrix;
#ifndef VECTOR_H
#define VECTOR_H
#include "ForwardDeclarations.h"
#include <vector>
#include <string>


using u_int = unsigned int;
using u_long = unsigned long;



class Vector{
	private:
		std::vector<double> _vector;
	public:
		// ----------------------------constructors------------------------------
		explicit Vector(const std::vector<double>& v);

		explicit Vector(int n, double element = 0.0);
		
		// -----------------------------------main fucntions -------------------------------------
		/**
		 * @description normaliation of a vector ~ length is 1
		 */
		void normalize();

		[[nodiscard]] Vector vectorTripleProduct(Vector& v1, const Vector& v2);

		[[nodiscard]] Vector projection(const Vector& v1) const;

		[[nodiscard]] Vector projection(Vector v, double angle) const;

		[[nodiscard]] double scalarTripleProduct(Vector v1, const Vector& v2) const;

		/**
		 * @return dimension of the vector
		 */
		[[nodiscard]] int size() const{
			return (int) _vector.size();
		}
		/**
		 * @desc currently only for 2d and 3d vectors
		 */
		void rotateX(double angle);
		/**
		 * @desc currently only for 2d and 3d vectors
		 */
		void rotateY(double angle);
		/**
		 * @desc currently only for 3d vectors
		 */
		void rotateZ(double angle);

		/**
		 * @return are vectors perpendicular to each other
		 * @param comparing vector
		 */

		[[nodiscard]] bool isPerpendicular(const Vector& v) const;

		/**
		 * @return are vectors coplanar to eachother
		 * @param v comparing vector
		 */

		[[nodiscard]] bool isCoplanar(const Vector& v) const;

		/**
		 * @return are two vectors parallel to each-other
		 * @param v comparing vector
		 */

		[[nodiscard]] bool isParallel(const Vector& v) const;

		/**
		 * @return Lenght of a vector
		 */

		[[nodiscard]] double length() const;

		/**
		 * @return are two vectors orthogonal
		 * @param v vector comparing
		 */

		[[nodiscard]] bool isOrthogonal(const Vector& v) const;

		/**
		 * multiply vector by a scalar
		 */

		void scale(double scalar);

		/**
		 * Cross product
		 * @param rhs second vector in cross product
		 * @return cross product of two vectors
		 */
		Vector crossProduct(const Vector& rhs);

		/**
		 *
		 * @return angle between two vectors
		 */
		 [[nodiscard]] double angle(const Vector& rhs, bool degrees = false) const;

		 /**
		  *@return is a null vector
		  * */
		[[nodiscard]] bool isNull() const;
		//------------------------------------------operators--------------------------------------------
		/**
		 * @desc vector addition
		 * @return vector addition
		 * @param rhs Vector that is being added
		 */
		Vector operator+(const Vector &rhs);

		/**
		 * vector subtraction
		 */
		Vector operator-(const Vector& rhs);

		/**
		 * vector subtraction
		 */
		void operator-=(const Vector& rhs);

		void operator+=(const Vector& rhs);

		/**
		 * dot/inner product
		 * @return
		 */
		double operator*(const Vector& rhs) const;

		/**
		 * @return Scaled vector
		 */
		Vector operator*(double scalar) ;

		void operator*=(double scalar);

		/**
		 *@return element on index position. Can be used for setting/getting
		 */
		double& operator[](int index) ;

		/**
		 *@return element on index position. Read-only
		 */

		const double& operator[](int index) const;

		// ----------------------------------------------------helpers-----------------------------------------------
		/**
		 * @return Vector as a Matrix object
		 */
		[[nodiscard]] Matrix asMatrix() const;

		/**
		 *
		 */
		std::vector<double> asVector(){
			return _vector;
		}

		std::string toString();

		[[nodiscard]] Vector copy() const {
			return Vector(_vector);
		}

};



#endif //VECTOR_H
