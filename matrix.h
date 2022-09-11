#pragma once
#include <iostream>
#include <math.h>
#include <vector>

using namespace std;

class Vector {
	vector<double> columns;
public:
	Vector(vector<double> v) {
		columns = v;
	}

	Vector(int n) {
		columns.resize(n);
	}

	Vector() {}

	double& operator [] (int j) {
		return columns[j];
	}

	Vector operator = (Vector temp) {
		if (columns.size() != temp.get().size())
			resize(temp.get().size());
		for (int i = 0; i < columns.size(); i++)
			columns[i] = temp.get()[i];
		return (*this);
	}

	Vector operator = (vector< double> temp) {
		if (columns.size() != temp.size())
			resize(temp.size());
		for (int i = 0; i < columns.size(); i++)
			columns[i] = temp[i];
		return (*this);
	}

	vector<double> get() {
		return columns;
	}

	int size() const {
		return columns.size();
	}

	//double get(int i) { return columns[i]; }

	void set(vector<double>v) {
		columns = v;
	}

	void set(int i, double value) { columns[i] = value; }

	void resize(int i) {
		columns.resize(i);
	}

	Vector operator + (Vector v) const {
		Vector temp(columns.size());
		for (int i = 0; i < columns.size(); i++) {
			temp.set(i, columns[i] + v[i]);
		}
		return temp;
	}

	Vector operator - (Vector v) const {
		Vector result(columns.size());
		for (int i = 0; i < columns.size(); i++) {
			result[i] = columns[i] - v[i];
		}
		return result;
	}

	Vector operator *(double v) const {
		Vector result(columns.size());
		for (int i = 0; i < columns.size(); i++)
			result[i] = columns[i] * v;
		return result;
	}

	Vector operator / (double v)const {
		Vector result(columns.size());
		for (int i = 0; i < columns.size(); i++)
			result[i] = columns[i] / v;
		return result;
	}

	friend Vector operator * (double tmp, Vector v) {
		Vector temp(v.get().size());
		for (int i = 0; i < v.get().size(); i++)
			temp.set(i, v[i] * tmp);
		return temp;
	}

	double operator * (Vector v) {
		double sum = 0;
		for (int i = 0; i < columns.size(); i++)
			sum += columns[i] * v[i];
		return sum;
	}

	double cosPhi(Vector v) {
		auto a = (*this);
		auto b = v;
		double scalar = a * b;
		double denominator = a.modul() * b.modul();
		return scalar / denominator;
	}


	double modul() {
		double sum = 0;
		for (int i = 0; i < columns.size(); i++)
			sum = sum + columns[i] * columns[i];
		return sqrt(sum);
	}

	~Vector() {}

	void print() {
		int n = columns.size();
		for (int i = 0; i < n; i++)
			cout << columns[i] << " ";
		cout << endl;
	}
};



class Matrix {
	vector<vector<double>> p;
	int rows; // i
	int columns; // j

public:

	Matrix() : Matrix(0) { }

	Matrix(vector<vector<double>> p) {
		this->p = p;
		rows = p.size();
		columns = p[0].size();
	}

	Matrix(int ni, int nj) {
		this->rows = ni;
		this->columns = nj;
		this->p = vector<vector<double>>(ni);
		for (int i = 0; i < ni; i++) {
			p.push_back(vector<double>(nj));
		}
	}

	Matrix(int n) {
		this->rows = n;
		this->columns = n;
		for (int i = 0; i < n; i++) {
			p.push_back(vector<double>(n));
		}
	}

	int NI() const {
		return rows;
	}

	int NJ() const {
		return columns;
	}




	double& operator()(int i, int j)
	{
		return p[i][j];
	}

	Vector operator * (Vector v) const {
		return (*this) * v.get();
	}

	Vector operator * (vector<double> vec) const {
		auto columns = p;
		Vector production(vec.size());
		auto ni = NI();
		auto nj = NJ();
		for (int i = 0; i < ni; i++)
		{
			double sum = 0;
			for (int j = 0; j < nj; j++)
				sum += p[i][j] * vec[j];
			production.set(i, sum);
		}
		return production;
	}
};