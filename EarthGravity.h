#pragma once

#include "integrator.h"
#include "matrix.h"

#include <cmath>

#include <iostream>
#include <fstream>


const double mu = 398600.4418e9; // m3/s2

class GravityModel {
public:
	virtual Vector getAcceleration(Vector x) = 0;
};


class EarthGravity : public GravityModel {
public:
	const double
		OMEGA = 7.292115E-5, // rad/sec
		a = 6378.136e3; // m
};


Matrix A(double W, double Theta, double Omega1, double I) {
	Matrix A(3);
	A(0, 0) = cos(W + Theta) * cos(Omega1) - sin(W + Theta) * sin(Omega1) * cos(I);
	A(0, 1) = -sin(W + Theta) * cos(Omega1) - cos(W + Theta) * sin(Omega1) * cos(I);
	A(0, 2) = sin(I) * sin(Omega1);
	A(1, 0) = cos(W + Theta) * sin(Omega1) + sin(W + Theta) * cos(Omega1) * cos(I);
	A(1, 1) = -sin(W + Theta) * sin(Omega1) + cos(W + Theta) * cos(Omega1) * cos(I);
	A(1, 2) = -sin(I) * cos(Omega1);
	A(2, 0) = sin(W + Theta) * sin(I);
	A(2, 1) = cos(W + Theta) * sin(I);
	A(2, 2) = cos(I);
	return A;
}

double getfocalparam(double Alpha, double E) {
	return Alpha * (1 - E * E);
}

Vector getRosculating(double Alpha, double E, double Theta) {
	Vector Rosculating = Vector(3);
	Rosculating.set(0, getfocalparam(Alpha, E) / (1 + E * cos(Theta)));
	Rosculating.set(1, 0);
	Rosculating.set(2, 0);
	return Rosculating;
}

Vector getVosculating(double Alpha, double E, double Theta) {
	Vector Vosculating = Vector(3);
	auto fp = getfocalparam(Alpha, E);
	Vosculating[0] = sqrt(mu / fp) * E * sin(Theta);
	Vosculating[1] = sqrt(mu / fp) * (1 + E * cos(Theta));
	Vosculating[2] = 0;
	return Vosculating;
}




class CentralGravitation : public EarthGravity {
public:

	Vector getAcceleration(Vector x) override {
		return -mu * x / pow(x.modul(), 3);
	}
};









class NormalGravitation : public EarthGravity {


public:
	double delta(int k) { // ok
		return k == 0 ? 0.5 : 1;
	}

	double J(int n) { // totally ok
		if (n == 2) { return  1082.62575E-6; }
		if (n == 4) { return  -2.37089E-6; }
		if (n == 6) { return  6.08E-9; }
		return  -1.40E-11;
	}

	double C0(int n) { // totally ok
		return -J(n) / sqrt(2 * n + 1);
	}

	double Precderivative(int n, double phi) {
		return sqrt(0.5 * n * (n + 1)) * Prec(n, 1, phi);
	}

	double Prec(int n, int m, double phi) {
		if (n < m)
			return 0;
		if (n == m && n == 0 && m == 0)
			return 1;
		if (n > m)
		{
			return
				Prec(n - 1, m, phi) *
				sin(phi) *
				sqrt((4 * n * n - 1) / (double)(n * n - m * m)) -
				Prec(n - 2, m, phi) *
				sqrt((pow(n - 1, 2) - m * m) * (2 * n + 1) / (double)((n * n - m * m) * (2 * n - 3)));
		}
		// n == m != 0
		return Prec(n - 1, m - 1, phi) *
			cos(phi) *
			sqrt((2 * n + 1) / (2 * n * delta(m - 1)));
	}


	double Sum1(double ro, double phi) {
		double sum1 = 0;
		for (int n = 2; n <= 8; n = n + 2)
			sum1 +=
			(n + 1) * // ok
			pow(a / ro, n) *
			C0(n) * // ok
			Prec(n, 0, phi);
		return sum1;
	}

	double Sum2(double ro, double phi) {
		double sum = 0;
		for (int n = 2; n <= 8; n = n + 2)
			sum += pow(a / ro, n) * C0(n) * Precderivative(n, phi);
		return sum;
	}

	Vector gosphere(Vector v) {
		Vector sph(3);
		double x = v[0];
		double y = v[1];
		double z = v[2];
		sph[0] = atan2(y, x);
		sph[1] = atan2(z, sqrt(x * x + y * y));
		sph[2] = sqrt(x * x + y * y + z * z);
		return sph;
	}


	Matrix gorectangularISK(Vector vec, double r0, double ro) {
		Matrix Res(3);
		double x = vec[0];
		double y = vec[1];
		double z = vec[2];
		bool r0Is0 = r0 == 0; // todo compare to epsilon
		Res(0, 0) = x / ro;
		if (!r0Is0)
		{
			Res(0, 1) = -(x * z) / (ro * r0);
			Res(0, 2) = -y / r0;
		}

		Res(1, 0) = y / ro;
		if (!r0Is0)
		{
			Res(1, 1) = -(y * z) / (ro * r0);
			Res(1, 2) = x / r0;
		}

		Res(2, 0) = z / ro;
		Res(2, 1) = r0 / ro;
		Res(2, 2) = 0;

		return Res;
	}


	Vector getAcceleration(Vector v) override
	{
		Vector lpr = gosphere(v);

		double lambda = lpr[0];
		double phi = lpr[1];
		double ro = lpr[2];

		Vector gPPL(3);
		gPPL[0] = -mu / pow(ro, 2) * (1 + Sum1(ro, phi)); // todo invalid
		//cout << endl << gPPL[0] << endl;
		gPPL[1] = mu / pow(ro, 2) * Sum2(ro, phi);
		gPPL[2] = 0;

		double r0 = sqrt(v[0] * v[0] + v[1] * v[1]);
		Matrix A = gorectangularISK(v, r0, ro);
		Vector gXYZ = A * gPPL;

		return gXYZ.get();
	}
};






void readTM60(vector<double>& es, vector<Vector>& cs) {
	ifstream fin;
	fin.open("tm60.txt");
	if (!fin)
		return;

	for (int j = 0; j < 60; j++)
	{
		vector<double> v;
		for (int i = 0; i < 5 && !fin.eof(); i++)
		{
			double x;
			fin >> x;
			v.push_back(x);
		}
		double
			e = v[1] / 1e10,
			x = v[2] * 1e3,
			y = v[3] * 1e3,
			z = v[4] * 1e3;
		es.push_back(e);
		cs.push_back(Vector(vector<double> { x, y, z }));
	}
	fin.close();
}



class AnomalGravitation : public NormalGravitation {
	vector<double> es;
	vector<Vector> cs;
	const double N = 60;

public:

	AnomalGravitation()
	{
		readTM60(es, cs);
	}


	Vector getAcceleration(Vector x) override {
		Vector gn = NormalGravitation::getAcceleration(x);
		Vector dg(3);
		for (int i = 0; i < N; i++)
		{
			Vector r = x - cs[i];
			dg = dg + es[i] * (x - cs[i]) / pow(r.modul(), 3);
		}
		dg = -mu * dg;
		Vector res = gn + dg;
		return res;
	}
};






struct OrbitParameters
{
	double Alpha = 0;
	double Theta = 0;
	double E = 0;
	double I = 0;
	double W = 0;
	double w = 0;
};


class Satellite : public IMathModel {
	OrbitParameters op;
	GravityModel* gm;

public:
	Satellite(OrbitParameters op, GravityModel* gm) : op(op), gm(gm)
	{
		IMathModel::SetX0(getX0());
	}

protected:
	vector<double> getX0() {
		double Alpha = op.Alpha;
		double E = op.E;
		double Theta = op.Theta;
		double Omega1 = op.w;
		double I = op.I;
		double W = op.W;

		Vector x0(6);
		Vector r1 = getRosculating(Alpha, E, Theta);
		Matrix A1 = A(W, Theta, Omega1, I);
		Vector r = A1 * r1;
		for (int i = 0; i < r.size(); i++)
			x0[i] = r[i];
		Vector v1 = getVosculating(Alpha, E, Theta);
		Vector v = A1 * v1;
		for (int i = v1.size(); i < 2 * v1.size(); i++)
			x0[i] = v[i - v1.size()];
		return x0.get();
	}


	void addResult(vector<double> vec, double time) override {
		for (int i = 0; i < 6; i++)
			vec[i] /= 1000;

		vec.push_back(time);

		Results.push_back(vec);
	}


	vector<double> f(double t, vector<double> x) override {
		vector<double> y(6);
		y[0] = x[3];
		y[1] = x[4];
		y[2] = x[5];
		auto g = gm->getAcceleration(Vector(vector<double> { x[0], x[1], x[2] }));
		y[3] = g[0];
		y[4] = g[1];
		y[5] = g[2];
		return y;
	}
};
