#pragma once

#include <iostream>
#include <vector>
#include <stdlib.h>
#include <fstream>
#include <random>
//#include <graphics.h>
#include <conio.h>
#include "attractor_ui.h"

using namespace std;
const double g = 9.8;
const int m = 2;
const int a = 1;
const int n = 5;
const double Cx = 0.47;
const double S = 0.00785;
const double mass = 1;


class gauss {
private:
	double mean;
	double stddev;

	double uniformRand()
	{
		return static_cast<double>(rand()) / RAND_MAX;
	}

public:

	double normalRand()
	{
		return sqrt(-2.0 * log(uniformRand())) * cos(2.0 * 3.14 * uniformRand()) * stddev + mean;
	}

	gauss(double mean, double stddev) {
		this->mean = mean;
		this->stddev = stddev;
	}
};

class IMathModel {
public:
	bool doPrintAll = false;

protected:
	vector<double> X;
	vector<vector<double>> Vectors;
	vector<vector<double>> Results;
public:


	virtual void SetX0(vector<double> x0) {
		X = x0;
		Vectors.push_back(x0);
		if (doPrintAll)
			printLast();
	}

	virtual double VGX() { return 0; }
	virtual double VGY() { return 0; }

	vector<double> X0() { return X; }
	vector<vector<double>> X_all() { return Vectors; }
	vector<vector<double>> GetResults() { return Results; }

	int size() { return X.size(); }

	void print(vector<vector<double>> vectors)
	{
		int ni = vectors.size();
		for (int i = 0; i < ni; i++) {
			auto v = vectors[i];
			int nj = v.size();
			for (int j = 0; j < nj; j++) {
				cout << v[j] << " ";
			}
			cout << endl;
		}
		cout << endl;
	}

	void printAll() {
		print(Vectors);
	}

	void print() {
		int n = X.size();
		for (int i = 0; i < n; i++) {
			cout << X[i] << " ";
		}
		cout << endl;
	}

	void printResults() {
		print(Results);
	}

	void printLast() {
		auto v = Vectors[Vectors.size() - 1];
		int n = v.size();
		for (int i = 0; i < n; i++) {
			cout << v[i] << " ";
		}
		cout << endl;
	}

	virtual vector<double> RigthParts() { return X; }

	virtual vector<double> f(double t, vector<double> r) { return r; }

	virtual void addResult(vector<double> vec, double time) { }
};



class MotionModel : public IMathModel {
protected:
	double VGx;
	double VGy;
	gauss* Gauss;
public:
	MotionModel(vector<double>m, gauss* Gauss) {
		this->Gauss = Gauss;
		SetX0(m);
		VGx = Gauss->normalRand();
		VGy = Gauss->normalRand();
	}

	double VGX() override { return VGx; }
	double VGY() override { return VGy; }

	vector<double> RigthParts() override {
		vector<double> vTmp(n);
		double x = X[0];
		double y = X[1];
		double Vx = X[2];
		double Vy = X[3];
		vTmp[0] = Vx;
		vTmp[1] = Vy;
		vTmp[2] = -sqrt(Vx * Vx + VGx * VGx) * (Vx - VGx) * Cx * S * 0.5 / mass;
		vTmp[3] = -sqrt(Vy * Vy + VGy * VGy) * (Vy - VGy) * Cx * S * 0.5 / mass - g;
		return vTmp;
	}

	void SetX0(vector<double> x0) override {
		X = x0;
	}
};


class attractor : public IMathModel {
public:
	attractor(vector<double> m) {
		SetX0(m);
	}

	vector<double> RigthParts() override {
		vector<double> vTmp(n);
		double x = X[0];
		double y = X[1];
		double z = X[2];
		double Vy = X[3];

		vTmp[0] = 10 * (y - x);
		vTmp[1] = x * (28 - z) - y;
		vTmp[2] = x * y - 8.0 / 3.0 * z;
		vTmp[3] = 0;
		return vTmp;
	}

	//void SetX0(vector<double> x0) override {
	//	X = x0;
	//}
};

struct integratorArgs {
	integratorArgs(IMathModel* m, double tk) {
		model = m;
		this->tk = tk;
	}
	IMathModel* model = nullptr;
	double tk;
};

struct eylerArgs : integratorArgs {
	eylerArgs(IMathModel* m, double tk, double e, double dt, double t0) :integratorArgs(m, tk), e(e), dt(dt), t0(t0) {};
	double e;
	double dt;
	double t0;
};

struct rungeArgs : integratorArgs {
	rungeArgs(IMathModel* m, double tk, double h) : integratorArgs(m, tk), h(h) {};
	double h;
};

struct eylerCorrectionArgs : integratorArgs {
	eylerCorrectionArgs(IMathModel* m, double tk, double h, double dt, double t0) :integratorArgs(m, tk), h(h), dt(dt), t0(t0) {};
	double dt;
	double h;
	double t0;
};

class integrator {
protected:
	IMathModel* model;
	double tk;
public:

	integrator(integratorArgs args) : tk(args.tk), model(args.model) { }

	~integrator() {
		cout << "class is destroyed" << endl;
	}

	vector<double> RigthParts() {
		return model->RigthParts();
	}

	void SetX0(vector<double> x0) {
		model->SetX0(x0);
		//cout << "X = " << x0[0] << " Y = " << x0[1] << " Vx = " << x0[2] << endl;
	}

	void print() {
		model->print();
	}
	void printLast() {
		model->printLast();
	}
	virtual void integration() = 0;
};

class eyler : public integrator {
protected:
	ofstream* fout;
	const double e;
	const double dt;
	const double tk;
	const double t0;
public:

	void printF(vector<double> x0, double t) {
		ofstream& pFout = *fout;
		pFout << "   t=" << t << "    x=" << x0[0] << "    y=" << x0[1] << endl;
	}

	eyler(eylerArgs args) : integrator(args), dt(args.dt), e(args.e), tk(args.tk), t0(args.t0) {
		fout = new ofstream("f.txt");
	}

	~eyler() {
		fout->close();
		delete fout;
	}

	double t0_() { return t0; }

	void integration() override {
		double t = t0_();
		int n = model->size();
		vector<double> vPrev = model->X0();
		double y = vPrev[1];
		double x = vPrev[0];
		while (t < tk) {
			vector<double> f = RigthParts();
			vector<double> vNext(n);
			t = t + dt;
			for (int i = 0; i < n; i++)
				vNext[i] = vPrev[i] + f[i] * dt;

			/// ATTRACTOR
			SetX0(vNext);
			printF(vNext, t);
			vPrev = vNext;

			/// VECTOR
			 /* if (vNext[1] >= 0)
			{
				SetX0(vNext);
				printF(vNext, t);
				vPrev = vNext;
				y = vNext[1];
			}
			else
				break; */
		}
		//cout << "S=" << model->X0()[0] - x << endl; 
	}
};

class runge : public integrator {
protected:
	ofstream* fout;
	const double h;
public:
	runge(rungeArgs args) : integrator(args), h(args.h) {
		fout = new ofstream("f.txt");
	}

	~runge() {
		fout->close();
		delete fout;
	}

	void printF(vector<double> x0, double t) {
		ofstream& pFout = *fout;
		pFout << "   t=" << t << "    x=" << x0[0] << "    y=" << x0[1] << endl;
	}

	struct Point {
		double x, y, vx, vy;

		Point() {}
		Point(double x, double y, double vx, double vy) : x(x), y(y), vx(vx), vy(vy) { }
		Point add(double d) {
			return Point(x + d, y + d, vx + d, vy + d);
		}
	};

	double f(int i, double t, Point point) {
		switch (i) {
		case 0:
			/// ATTRACTOR
			return 10 * (point.y - point.x);

			/// VECTOR
			//return point.vx;
		case 1:
			/// ATTRACTOR
			return point.x * (28 - point.vx) - point.y;

			/// VECTOR
			//return point.vy;
		case 2: {
			/// ATTRACTOR
			return point.x * point.y - 8.0 / 3.0 * point.vx;

			/// VECTOR
			//double Vx = point.vx;
			//double VGx = this->model->VGX();
			//return -sqrt(Vx * Vx + VGx * VGx) * (Vx - VGx) * Cx * S * 0.5 / mass;
		}
		case 3:
			double Vy = point.vy;
			double VGy = this->model->VGY();
			return -sqrt(Vy * Vy + VGy * VGy) * (Vy - VGy) * Cx * S * 0.5 / mass - g;
		}
	}

	double oneStep(int i, vector<double> vPrev, double t) {
		Point columns;
		columns.x = vPrev[0];
		columns.y = vPrev[1];
		columns.vx = vPrev[2];
		columns.vy = vPrev[3];

		double k1 = h * f(i, t, columns);
		double k2 = h * f(i, t, columns.add(k1 / 2));
		double k3 = h * f(i, t, columns.add(k2 / 2));
		double k4 = h * f(i, t, columns.add(k3));
		double res = (vPrev[i] + 1.0 / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4));
		return res;
	}

	void integration() override {
		int n_count = 1000;
		int n = model->size();
		double x = model->X0()[0];
		double t = 0;
		vector<double> vPrev = model->X0();
		for (int i = 0; i < n_count; i++) {
			t = t + 1;
			vector<double> vNext(n);
			for (int j = 0; j < n; j++) {
				vNext[j] = oneStep(j, vPrev, t);
			}
			vPrev = vNext;
			/// ATTRACTOR
			SetX0(vNext);
			printF(vNext, t);
			vPrev = vNext;

			/// VECTOR
			 /* if (vNext[1] >= 0)
			{
				SetX0(vNext);
				printF(vNext, t);
				vPrev = vNext;
				y = vNext[1];
			}
			else
				break; */
		}
	}
};

class eylerCorrection : public integrator {
protected:
	ofstream* fout;
	const double h;
	const double dt;
	const double tk;
	const double t0;
public:

	eylerCorrection(eylerCorrectionArgs args) : integrator(args), dt(args.dt), h(args.h), tk(args.tk), t0(args.t0) { fout = new ofstream("f.txt"); }
	~eylerCorrection() {
		fout->close();
		delete fout;
	}

	void printF(vector<double> x0, double t) {
		ofstream& pFout = *fout;
		pFout << "   t=" << t << "    x=" << x0[0] << "    y=" << x0[1] << endl;
	}

	double t0_() { return t0; }

	struct Point {
		double x, y, vx, vy;

		Point() {}
		Point(double x, double y, double vx, double vy) : x(x), y(y), vx(vx), vy(vy) { }

		Point addX(double d) {
			return Point(x + d, y, vx, vy);
		}
		Point addY(double d) {
			return Point(x, y + d, vx, vy);
		}
		Point addVx(double d) {
			return Point(x, y, vx + d, vy);
		}
		Point addVy(double d) {
			return Point(x, y, vx, vy + d);
		}
	};



	double f(int i, Point point) {
		switch (i) {
		case 0:
			return point.vx;
		case 1:
			return point.vy;
		case 2: {
			double Vx = point.vx;
			double VGx = this->model->VGX();
			return -sqrt(Vx * Vx + VGx * VGx) * (Vx - VGx) * Cx * S * 0.5 / mass;
		}
		case 3:
			double Vy = point.vy;
			double VGy = this->model->VGY();
			return -sqrt(Vy * Vy + VGy * VGy) * (Vy - VGy) * Cx * S * 0.5 / mass - g;
		}
	}


	double oneStep(int i, vector<double> vPrev) {
		Point columns;
		columns.x = vPrev[0];
		columns.y = vPrev[1];
		columns.vx = vPrev[2];
		columns.vy = vPrev[3];
		switch (i) {
		case 0:
			return (columns.x + h / 2 * (f(i, columns) + f(i, columns.addX(h * f(i, columns)))));
		case 1:
			return (columns.y + h / 2 * (f(i, columns) + f(i, columns.addY(h * f(i, columns)))));
		case 2:
			return (columns.vx + h / 2 * (f(i, columns) + f(i, columns.addVx(h * f(i, columns)))));
		case 3:
			return (columns.vy + h / 2 * (f(i, columns) + f(i, columns.addVy(h * f(i, columns)))));
		}
	}


	void integration() override {
		int n = model->size();
		double t = t0_();
		double x = model->X0()[0];
		vector<double> vPrev = model->X0();
		while (t < tk) {
			t = t + dt;
			vector<double> vNext(n);
			for (int j = 0; j < n; j++) {
				vNext[j] = oneStep(j, vPrev);
			}
			if (vNext[1] < 0)
				break;
			else
			{
				SetX0(vNext);
				printF(vNext, t);
				vPrev = vNext;
			}
		}
		cout << "S=" << (model->X0()[0] - x) << endl;
	}
};
