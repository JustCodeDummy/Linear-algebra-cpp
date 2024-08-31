//matrix.h
#ifndef ALGEBRAFORSOMEREASON_MATRIX_H
#define ALGEBRAFORSOMEREASON_MATRIX_H

#include <vector>
#include "Vector.h"
#include <string>
using u_long = unsigned long;
using u_int = unsigned int;
class Matrix{

    private:
        std::vector<Vector> matrix;

	public:

    // Constructors
		explicit Matrix(const std::vector<Vector> & _matrix){
			matrix = _matrix;
		};

		explicit Matrix(){
			matrix = std::vector<Vector>();
		}

		explicit Matrix( const std::vector<std::vector<double>>& v);

		Matrix(int rows, int columns, double element=0.0);

		// Functions
		/**
		 * @return Returns determinant of the matrix
		 */

		double determinant();
		/**
		 * @return the inverse of a matrix
		 */
		Matrix inverse();

		/**
		 * Performs transposition of matrix
		 */
		[[nodiscard]] Matrix transpose() const;


		/**
		 * @return Return a cofactor of a matrix
		 */
		 Matrix cofactor(int row, int column);

		/**
		 * Multiply matrix by scalar
		*/

		void scale(double scalar);

		/**
		 *
		 * @param row1 First row to be swapped
		 * @param row2 Second row to be swapped
		 */

		void swapRows(int row1, int row2);

		/**
		 *
		 * @param row Index of row to be scaled
		 * @param scalar Scalar
		 */
		void scaleRow(int row, double scalar);

		/**
		 *
		 * @param row1 Index to row that will be modified
		 * @param row2 Index to row that will be added
		 * @param scalar Scalar to multiply row before adding it
		 */

		void addRows(int row1, int row2, double scalar=1);
		/**
		 *
		 * @return Matrix in reduced row echelon form
		 */
		Matrix rref();

		 int rowCount(){
			 return (int) matrix.size();
		 }
		 int columnCount(){
			 return (int) matrix[0].size();
		 }

		 void add(const Vector& v);

		 void add(const Vector& v, int index);

		 Matrix orthonormalize() const;

		Vector vectorMatrixMultiplication(const Vector& v);

		Vector toVector();

		// Operators

		Matrix operator*(Matrix& _matrix);

		Vector& operator[](int index){
			return matrix[index];
		}

		void operator*=(Matrix& rhs);

		bool operator==(Matrix &rhs);

		bool operator !=( Matrix &rhs);

		Matrix operator+(Matrix& rhs);

		Matrix operator-(Matrix& rhs);

		void operator+=(Matrix& rhs);

		void operator-=(Matrix& rhs);
		// Helpers
		std::string toString();

		Vector asVector(){
			return matrix[0];
		}

		std::vector<Vector> asVectorList(){
			return matrix;
		}

		std::vector<std::vector<double>> asStdVectorList();
};

#endif //ALGEBRAFORSOMEREASON_MATRIX_H

