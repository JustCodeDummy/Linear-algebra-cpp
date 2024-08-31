//
// Created by fsociety on 10/09/2023.
//

#include "../Header files/matrix.h"
#include <cmath>
#include <iostream>
#include <string>
using u_long = unsigned long;
using u_int = unsigned int;

Matrix::Matrix( const std::vector<std::vector<double>>& v) {
	for (const std::vector<double>& u : v){
		this->matrix.emplace_back(u);
	}
}

Matrix::Matrix(int rows, int columns, double element) {
	this->matrix = std::vector<Vector>(rows, Vector(columns, element));
}

double Matrix::determinant() {

	//TODO exception if not square
	if(this->matrix.size() == 1){
		return this->matrix[0][0];
	}
	if (this->matrix.size() == 2){
		return this->matrix[0][0] * this->matrix[1][1] - this->matrix[0][1] * this->matrix[1][0];
	}

	std::vector<double> determinants;

	determinants.reserve(this->matrix.size());
	for (int i = 0; i<this->matrix.size(); i++){
		determinants.push_back(Matrix::cofactor(0, i).determinant());
	}

	double det = 0.0;

	for (int i =0; i<determinants.size(); i++){
		determinants[i]*=this->matrix[0][i];
	}

	for (int i=0; i<determinants.size();i++){
		det += pow(-1, i) * determinants[i];
	}

	return det;
}

Matrix Matrix::cofactor(int row, int column) {

	std::vector<std::vector<double>> temp = this->asStdVectorList();
	temp.erase(temp.begin() + row);

	for (int i=0; i< temp[0].size(); i++){
		temp[i].erase(temp[i].begin()+column);
	}
	return Matrix(temp);

}

void Matrix::scale(double scalar){
	for ( Vector & i : this->matrix){
		i.scale(scalar);
	}
}

Matrix Matrix::inverse(){

	// TODO det=0 & !square exceptions

	double det = 1/this->determinant();

	std::vector<std::vector<double>> inverse;
	inverse.reserve(this->matrix.size());
	for (int i = 0; i<this->matrix.size(); i++){
		std::vector<double> row;
		row.reserve(this->matrix[0].size());
		for (int j=0; j<this->matrix[0].size(); j++){
			row.push_back(Matrix::cofactor(i, j).determinant());
		}
		inverse.push_back(row);
	}

	Matrix in(inverse);
	in.scale(det);
	return in.transpose();

}

Matrix Matrix::transpose() const{
	unsigned long rows = this->matrix.size();
	unsigned long cols = this->matrix[0].size();
	std::vector<std::vector<double>> transposed(cols, std::vector<double>(rows));
	for (int i = 0; i<rows; i++){
		for (int j=0; j<cols; j++){
			transposed[j][i] = this->matrix[i][j];
		}
	}
	return Matrix(transposed);
}

double rounded(double number, u_int decimal){
	return (round(number * pow(10, decimal))/pow(10, decimal));
}

Vector Matrix::vectorMatrixMultiplication(const Vector& v) {
	Matrix a = v.asMatrix();
	Matrix b(this->matrix);
	return (b*a).asVector();
}

Vector Matrix::toVector() {
	return this->matrix[0];
}

Matrix Matrix::operator*( Matrix & _matrix) {
	double value;
	Matrix m(this->rowCount(), _matrix.columnCount());
	for (int k=0; k<this->rowCount(); k++){
		for (int i = 0; i<_matrix.columnCount(); i++){
			value = 0.0;
			for (int j = 0; j< this->columnCount();j++){
				value+= this->matrix[k][j] * _matrix[j][i];

			}
			m[k][i] = value;
		}
	}
	return m;
}

Matrix Matrix::operator+( Matrix &rhs) {
	Matrix m(this->matrix);
	for (int i = 0; i<rhs.rowCount();i++){
		for (int j=0; j<rhs.columnCount(); j++){
			m[i][j] += rhs[i][j];
		}
	}
	return m;
}

Matrix Matrix::operator-(Matrix &rhs) {
	Matrix m(this->matrix);
	for (int i = 0; i<rhs.rowCount();i++){
		for (int j=0; j<rhs.columnCount(); j++){
			m[i][j] -= rhs[i][j];
		}
	}
	return m;
}

void Matrix::operator+=(Matrix &rhs) {
	for (int i = 0; i<rhs.rowCount();i++){
		for (int j=0; j<rhs.columnCount(); j++){
			this->matrix[i][j] += rhs[i][j];
		}
	}
}

void Matrix::operator-=(Matrix &rhs) {
	for (int i = 0; i<rhs.rowCount();i++){
		for (int j=0; j<rhs.columnCount(); j++){
			this->matrix[i][j] -= rhs[i][j];
		}
	}
}

bool Matrix::operator==(Matrix &rhs) {
	if (this->columnCount() != rhs.columnCount() || this->rowCount() != rhs.rowCount()){
		return false;
	}
	for (int i = 0; i<this->rowCount(); i++){
		for (int j = 0; j<rowCount(); j++){
			if (this->matrix[i][j] != rhs[i][j]){
				return false;
			}
		}
	}
	return true;
}

void Matrix::operator*=(Matrix &rhs) {
	this->matrix = (Matrix(this->matrix) * rhs).asVectorList();
}

std::string Matrix::toString() {
	std::string s;
	for (Vector& v: this->matrix){
		s+= "| ";
		for (int i = 0;i<columnCount(); i++){
			s+= std::to_string(v[i]) + " ";

		}
		s+="|\n";
	}
	return s;

}

std::vector<std::vector<double>> Matrix::asStdVectorList() {
	std::vector<std::vector<double>> v;
	for (Vector u : this->matrix){
		v.push_back(u.asVector());
	}
	return v;
}

Matrix Matrix::orthonormalize() const {
	Matrix orthonormalized;
	Matrix basis = this->transpose();

	for (Vector v : basis.matrix){
		std::cout<<v.toString();
		for (int i=0; i<orthonormalized.rowCount(); i++){
			v -= orthonormalized[i] * (v * orthonormalized[i]);
		}
		v.normalize();
		orthonormalized.add(v);
	}
	return orthonormalized.transpose();
}

void Matrix::add(const Vector& v) {
	this->matrix.push_back(v);
}

void Matrix::add(const Vector &v, int index) {
	this->matrix.insert(this->matrix.begin() + index, v);
}

bool Matrix::operator!=( Matrix &rhs) {
	if (this->columnCount() != rhs.columnCount() || this->rowCount() != rhs.rowCount()){
		return true;
	}
	for (int i = 0; i<this->rowCount(); i++){
		for (int j = 0; j<rowCount(); j++){
			if (this->matrix[i][j] != rhs[i][j]){
				return true;
			}
		}
	}
	return false;
}

void Matrix::swapRows(int row1, int row2) {
	Vector t1 = this->matrix[row1];
	Vector t2 = this->matrix[row2];
	this->matrix[row1] = t2;
	this->matrix[row2] = t1;
}

void Matrix::scaleRow(int row, double scalar) {
	this->matrix[row] *= scalar;
}

void Matrix::addRows(int row1, int row2, double scalar) {
	this->matrix[row1]+=this->matrix[row2] * scalar;
}
// TODO testerej če dela
Matrix Matrix::rref() {
	Matrix rrefMatrix(this->matrix);
	const int cols = rrefMatrix.columnCount();
	const int rows = rrefMatrix.rowCount();
	int lead = 0;
	for (int r =0; r<rows;r++){
		if (lead>=cols){
			return rrefMatrix; // TODO break or return?
		}
		int i = r;
		while (rrefMatrix[i][lead] == 0){
			i++;
			if (i == rows){
				i = r;
				lead++;
				if (lead == cols){
					return rrefMatrix;
				}
			}
		}
		rrefMatrix.swapRows(i, r);
		rrefMatrix.scaleRow(r, 1/rrefMatrix[r][lead]);
		double lv = 1/rrefMatrix[r][lead];
		for (int j = 0; j<cols; j++){
			rrefMatrix[r][j] = rrefMatrix[r][j] * lv;
		}

		for (int j = 0; j<rows; j++){
			if (j!=r){
				rrefMatrix.addRows(r, j, -rrefMatrix[j][lead]);
			}
		}
		lead++;

	}
	return rrefMatrix;
}

