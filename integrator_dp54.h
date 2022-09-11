#pragma once
#define _CRT_SECURE_NO_WARNINGS

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <vector>

#include "integrator.h"
#include "matrix.h"
#include <ctime>

#define M_PIl 3.1415926535897932384626433832795028841968L

class dp54_Model : public IMathModel {
public:
	dp54_Model(vector<double> m) {
		SetX0(m);
	}

	vector<double> RigthParts() override {
		return X;
	}

	vector<double> f(double t, vector<double> r) override {
		const double mu = 0.012277471;

		double D1 = pow(pow(r[0] + mu, 2) + pow(r[1], 2), (3.0 / 2));
		double D2 = pow(pow(r[0] - (1 - mu), 2) + pow(r[1], 2), (3.0 / 2));

		vector<double> rd(4);
		rd[0] = r[2];
		rd[1] = r[3];
		rd[2] = r[0] + 2 * r[3] - (1 - mu) * (r[0] + mu) / D1 - mu * (r[0] - (1 - mu)) / D2;
		rd[3] = r[1] - 2 * r[2] - (1 - mu) * r[1] / D1 - mu * r[1] / D2;
		return rd;
	}
};

class EarthMovementModel : public IMathModel {
public:
	const long double mu = 132712.43994e+6;
	// x y z Vx Vy Vz
	EarthMovementModel(double t0, double x, double y, double z, double vx, double vy, double vz) {
		vector<double>  temp(7);
		temp[0] = t0;
		temp[1] = x;
		temp[2] = y;
		temp[3] = z;
		temp[4] = vx;
		temp[5] = vy;
		temp[6] = vz;
		SetX0(temp);
	}

	vector<double> RigthParts() override {
		return X;
	}

	vector<double> f(double t, vector<double> r) override {
		double tmp = 0;
		for (int i = 1; i < 7; i++) {
			tmp += r[i] * r[i];
		}

		vector< double> temp(7);
		temp[0] = 1;
		temp[1] = r[4]; //x
		temp[2] = r[5]; //y
		temp[3] = r[6]; //z

		temp[4] = -mu * r[1] / pow(sqrt(tmp), 3); //vx 	
		temp[5] = -mu * r[2] / pow(sqrt(tmp), 3); //vy
		temp[6] = -mu * r[3] / pow(sqrt(tmp), 3); //vz

		return temp;
	}
};

class myTime {
private:
	tm date = { 0 };

	void datePlusSeconds(struct tm* date, int seconds) {
		time_t date_seconds = mktime(date) + seconds;
		*date = *localtime(&date_seconds); ;
	}

public:

	myTime() { }

	myTime(int day, int month, int year) {
		date.tm_year = year - 1900;
		date.tm_mon = month - 1;
		date.tm_mday = day;
		addSeconds(0);
	}

	void addDays(int days) { 
		addSeconds(days * 86400);
	}

	void addSeconds(int seconds) {
		datePlusSeconds(&date, seconds);
	}

	int year() {
		return date.tm_year + 1900;
	}

	int month() {
		return date.tm_mon + 1;
	}

	int day() {
		return date.tm_mday;
	}

	int hour() {
		return date.tm_hour;
	}

	int minute() {
		return date.tm_min;
	}

	int seconds() {
		return date.tm_sec;
	}

	void print() {
		printf("%s", asctime(&date));
	}

	int dayOfWeek() {
		return date.tm_wday;
	}
};


class SunModel : public  EarthMovementModel {
private:
	const double
		R = 6371.3,
		w = 7.292115e-5;
	double
		phi = 0,
		lambda = 0,
		s0 = 0,
		l = 1;
	
	double tLastDay;
protected:
	myTime date;

public:

	SunModel(
		double phi, double lambda,
		int dd, int mm, int yyyy,
		double t0, double tf,
		double x, double y, double z, double vx, double vy, double vz) :
		EarthMovementModel(t0, x, y, z, vx, vy, vz)
	{
		this->phi = phi;
		this->lambda = lambda;
		this->tLastDay = tf - 86400.0;
		this->date = myTime(dd, mm, yyyy);
	}

	double siderialTime(
		double year, double month, double day, 
		double hour, double minute, double second)
	{
		int a = (14 - month) / 12;
		int M = month + 12 * a - 3;
		int Y = year + 4800 - a;
		double
			JDN = day + ((153 * M + 2) / 5) + 365 * Y + (Y / 4) - (Y / 100) + (Y / 400) - 32045,
			JD = JDN + (hour - 12) / 24 + minute / 1440 + second / 86400,
			tStar = int((JD - 2451544.5) / 36525),
			sG = 24110.54841 + 8640184.812866 * tStar + 0.093104 * tStar * tStar - 6.2e-6 * tStar * tStar * tStar,
			t = hour * 3600 + minute * 60 + second,
			lambdaRad = lambda * M_PIl / 180,
			S = sG * (M_PIl / 180) + w * t + lambdaRad;
		return S;
	}

	virtual bool mustProcessShadow(double time) {
		return time >= tLastDay;
	}

	virtual void processShadow(double time, Vector re, Vector rsh) {
		Results.push_back(rsh.get());
	}

	int dayOf(double time) {
		return time / 86400;
	}

	void addResult(vector<double> vec, double time) override {
		if (!mustProcessShadow(time))
			return;

		// CALC SHADOW begin
		Vector
			Re_(3),
			Re0(3), Rsh(4), Rg(3),
			Re({ vec[1], vec[2], vec[3] }),
			r(3),
			r0(3);

		myTime date = this->date;
		date.addSeconds(time);

		int Y = date.year();
		int M = date.month();
		int D = date.day();
		int h = date.hour();
		int m = date.minute();
		int sec = date.seconds();

		double sRad = siderialTime(
			Y, M, D,
			h, m, sec);
		double phiRad = phi * M_PIl / 180;

		const Matrix A({
				{	-cos(sRad) * sin(phiRad),	-sin(sRad) * sin(phiRad),	cos(phiRad)},
				{	cos(phiRad) * cos(sRad),	cos(phiRad) * sin(sRad),	sin(phiRad)},
				{	-sin(sRad),			cos(sRad),				0}
			});

		r.set(0, R * cos(phiRad) * cos(sRad));
		r.set(1, R * cos(phiRad) * sin(sRad));
		r.set(2, R * sin(phiRad));

		r0 = (1 / r.modul()) * r;
		Re0 = (1 / Re.modul()) * Re;
		Re_ = (-1 / (Re0 * r0)) * Re0;
		Rg = l * r0;
		Rsh = Rg + Re_;
		Rsh = A * Rsh.get();
		// CALC SHADOW end

		processShadow(time, Re, Rsh);
	}
};

class DayLengthModel : public SunModel {
	enum class State {
		LightBefore,
		Light,
		LightAfter
	};

private:
	int day_last;
	int t_sun_up = 0;
	int t_sun_down = 0;
	State state;
	
	bool enableEightLimit = true;
	//bool eightLimit = false;

	/*bool enableDaylight = false;*/
	bool enableDaylight = true;

	double utcSec;

public:
	DayLengthModel(
		double phi, double lambda, 
		int dd, int mm, int yyyy,
		double t0, double tf, double utc,
		double x, double y, double z, double vx, double vy, double vz) :
		SunModel(
			phi, lambda, 
			dd, mm, yyyy,
			t0,
			tf,
			x, y, z, vx, vy, vz)
	{
		this->day_last = dayOf(t0);
		this->utcSec = utc * 3600;
		dayReset();
	}

	bool mustProcessShadow(double time) override {
		return true;
	}


	bool isSummerTime(myTime date) {
		int m = date.month();
		int d = date.day();
		if (4 <= m && m <= 9)
			return true;
		if (m == 3 || m == 10) { // 31
			if (d > 20) {
				myTime dateSundayLast(31, m, date.year());
				while (dateSundayLast.dayOfWeek() != 0) {
					dateSundayLast.addDays(-1);
				}
				int dSundayLast = dateSundayLast.day();
				return dSundayLast <= d;
			}
		}
		return false;
	}


	void dayReset() {
		t_sun_up = t_sun_down = 0;
		state = State::LightBefore;
	}


	void processShadow(double timeUtc, Vector re, Vector rsh) override {
		int time = timeUtc + this->utcSec;
		double cosPhi = rsh.cosPhi(re);
		double phi = (acos(cosPhi));
		bool isLight = phi > M_PIl / 2;
		double angle = 180 / M_PIl * phi;

		int day = int(time / 86400);
		int h = int((time - day * 86400) / 3600);

		if (day != day_last) {
			// save previous day results
			double time_delta = t_sun_down - t_sun_up;
			vector<double> result{ (double)day_last, time_delta / 3600 };
			Results.push_back(result);
			cout << result[0] << " " << result[1] << endl << endl;
			dayReset();
		}
		day_last = day;

		bool conditionUp = true;
		if (enableDaylight) {
			myTime date = this->date;
			date.addSeconds(time);
			if (isSummerTime(date)) {
				h += 1;
			}
		}
		if (enableEightLimit) {
			conditionUp = 8 <= h && h <= 20;
		}

		if (state == State::LightBefore && isLight && conditionUp) {
			if (Results.size() == 0) {
				cout << "pi/2 = " << M_PIl / 2 << endl << endl;
			}
			cout << "up:   d = " << day << " h = " << h << " t = " << time << " phi = " << phi << " angle = " << angle << endl;
			t_sun_up = time;
			state = State::Light;
		}

		bool conditionDown = false;
		if (enableEightLimit) {
			conditionDown = h >= 20;
			// TODO store previous time
		}
		if (state == State::Light && (!isLight || conditionDown)) {
			cout << "down: d = " << day << " h = " << h << " t = " << time << " phi = " << phi << " angle = " << angle << endl;
			t_sun_down = time;
			state = State::LightAfter;
		}
	}
};





class pendulum_base_Model : public IMathModel {
public:
	pendulum_base_Model(double t0) {
		//	SetX0({ t0, 0, M_PIl });
		SetX0({ t0, 0, 10 });
	}

	vector<double> RigthParts() override {
		return X;
	}
};


class pendulum_ideal_spring_Model : public pendulum_base_Model {
public:
	pendulum_ideal_spring_Model(double t0) : pendulum_base_Model(t0) { }

	vector<double> f(double t, vector<double> r) override {
		vector<double> rd(3);
		const double k = 10;
		const double m = 1;
		rd[0] = 1;
		rd[1] = -k / m * r[2];
		rd[2] = r[1];
		return rd;
	}
};

class pendulum_ideal_Model : public pendulum_base_Model {
public:
	pendulum_ideal_Model(double t0) : pendulum_base_Model(t0) { }

	vector<double> f(double t, vector<double> r) override {
		const double l = 2;
		vector<double> rd(3);
		rd[0] = 1;
		rd[1] = -9.81 * sin(r[2]);
		rd[2] = 1.0 / l * r[1];
		return rd;
	}
};

class pendulum_real_Model : public pendulum_ideal_Model {
public:
	pendulum_real_Model(double t0) : pendulum_ideal_Model(t0) { }

	vector<double> f(double t, vector<double> r) override {
		const double b = 0.2;
		const double l = 0.2;
		vector<double> rd(3);
		rd[0] = 1;
		rd[1] = -9.81 * sin(r[2]) - b * r[1];
		rd[2] = 1.0 / l * r[1];
		return rd;
	}
};

class pendulum_real_spring_slide_Model : public pendulum_base_Model {
public:
	pendulum_real_spring_slide_Model(double t0) : pendulum_base_Model(t0) { }

	vector<double> f(double t, vector<double> r) override {
		vector<double> rd(3);
		const double k = 1;
		const double m = 2;
		const double C = 0.2;
		rd[0] = 1;
		rd[1] = -k / m * r[2] - C / m * r[1];
		rd[2] = r[1];
		return rd;
	}
};

class pendulum_real_viscous_Model : public IMathModel {
public:
	pendulum_real_viscous_Model(double t0) {
		SetX0({ t0, 0, 20 });
	}

	vector<double> RigthParts() override {
		return X;
	}

	vector<double> f(double t, vector<double> r) override {
		vector<double> rd(3);
		const double k = 5.5;
		const double m = 0.1;
		//const double C = 0.015; big
		double C;
		if (r[1] > 0)
		{
			C = -0.15;
		}
		else
			if (r[1] < 0)
			{
				C = 0.15;
			}
			else
			{
				C = 0.;
			}

		rd[0] = 1;
		rd[1] = -k / m * r[2] + C * 9.81;
		rd[2] = r[1];
		return rd;
	}
};

struct dp54_args : integratorArgs {
	double t0, hmin, hmax, h0, eps;
	double k;

public:
	dp54_args(
		IMathModel* m,
		double t0,
		double tf,
		double hmin,
		double hmax,
		double h0,
		double eps,
		double k)
		: integratorArgs(m, tf), t0(t0), hmin(hmin), hmax(hmax), h0(h0), eps(eps), k(k)
	{ }


	double getSamplingIncrement() {
		return k;
	}
};

double calcU() {
	double v = 1.0;
	double u;
	while (1.0 + v > 1.0) {
		u = v;
		v = v / 2;
	}
	return u;
}

class dp54_integrator : public integrator {
	const int n;
	double u;
	dp54_args args;

public:
	dp54_integrator(dp54_args args) : integrator(args), n(args.model->size()), args(args)
	{
		this->u = calcU();
	}

	virtual bool canContinue(int t) { return t < this->tk; }

	void integration() override {
		vector<double>
			v = this->model->X0(),
			x = v,
			x0 = v;
		double
			epsMax = args.eps,
			t0 = args.t0,
			t = args.t0,
			t_out = args.t0,
			h = args.h0,
			tf = this->tk,
			k1 = v[2];

		vector<double> k[6];

		double
			t1 = tf,
			h_new = args.h0,
			eps;
		vector<double> x1(n);

		int i;
		int j = 0;
		while (canContinue(t))
		{
			//print();
			h = h_new;
			if (t + h >= tf)
				h = tf - t;
			oneStep(h, t, x, epsMax, t1, x1, k, eps);

			h_new = h / (max(0.1, min(5.0, pow(eps / epsMax, 0.2) / 0.9)));

			if (eps > epsMax)
				continue;

			while ((t_out < t + h) && (t_out <= t1))
			{
				double
					theta = (t_out - t) / h,
					b[6];

				double T = theta;
				b[0] = T * (1 + T * (-1337 / 480.0 + T * (1039 / 360.0 + T * (-1163 / 1152.0))));
				b[1] = 0;
				b[2] = 100.0 * T * T * (1054 / 9275.0 + T * (-4682 / 27825.0 + T * (379 / 5565.0))) / 3.0;
				b[3] = -5.0 * T * T * (27 / 40.0 + T * (-9 / 5.0 + T * (83 / 96.0))) / 2.0;
				b[4] = 18225.0 * T * T * (-3 / 250.0 + T * (22 / 375.0 + T * (-37 / 600.0))) / 848.0;
				b[5] = -22.0 * T * T * (-3 / 10.0 + T * (29 / 30.0 + T * (-17 / 24.0))) / 7.0;

				vector<double> xOut(n);
				for (int mi = 0; mi < n; mi++) {
					xOut[mi] = x[mi];
					for (int mj = 0; mj < 6; mj++) {
						xOut[mi] += b[mj] * k[mj][mi];
					}
				}

				SetX0(xOut);
				j++;

				model->addResult(xOut, t_out);

				t_out += this->args.getSamplingIncrement();
			}
			for (i = 0; i < n; i++)
				x[i] = x1[i];
			t += h;
		}
		//cout << "t        X/Fi      " << endl;
		//cout <<t << "          " << k1-0.1353 << endl;
	//	printLast();
	}


	double getSamplingIncrement() {
		return 1; // TODO move to model
		//return 1e-1; // TODO move to model
	}

private:
	void oneStep(
		double h,
		double t0,
		vector<double> r0,
		double epsMax,
		double& t1,
		vector<double>& r1,
		vector<double> k[6],
		double& erra);

	void f(double t, vector<double> r, vector<double>& rd);
};



void dp54_integrator::f(double t, vector<double> x0, vector<double>& x1) {
	x1 = model->f(t, x0);
}

double _max(double d1, double d2, double d3, double d4) {
	return max(max(d1, d2), max(d3, d4));
}

void dp54_integrator::oneStep(
	double h,
	double t0,
	vector<double> x0,
	double epsMax,
	double& t1,
	vector<double>& x1,
	vector<double> k[6],
	double& eps)
{
	vector<double> k1(n), k2(n), k3(n), k4(n), k5(n), k6(n), k7(n);
	double
		a2 = 1.0L / 5,
		b21 = 1.0L / 5;
	double
		a3 = 3.0L / 10,
		b31 = 3.0L / 40,
		b32 = 9.0L / 40;
	double
		a4 = 4.0L / 5,
		b41 = 44.0L / 45,
		b42 = -56.0L / 15,
		b43 = 32.0L / 9;
	double
		a5 = 8.0L / 9,
		b51 = 19372.0L / 6561,
		b52 = -25360.0L / 2187,
		b53 = 64448.0L / 6561,
		b54 = -212.0L / 729;
	double
		a6 = 1.0L,
		b61 = 9017.0L / 3168,
		b62 = -355.0L / 33,
		b63 = 46732.0L / 5247,
		b64 = 49.0L / 176,
		b65 = -5103.0L / 18656;
	double
		a7 = 1.0L,
		b71 = 35.0L / 384,
		b73 = 500.0L / 1113,
		b74 = 125.0L / 192,
		b75 = -2187.0L / 6784,
		b76 = 11.0L / 84;
	double
		c1 = 35.0L / 384,
		c3 = 500.0L / 1113,
		c4 = 125.0L / 192,
		c5 = -2187.0L / 6784,
		c6 = 11.0L / 84;
	double tk0, tk1, tk2, tk3, tk4, tk5, tk6;
	vector<double> yk0(n), yk1(n), yk2(n), yk3(n), yk4(n), yk5(n), yk6(n);
	vector<double> rz(n);
	double s;
	int i;
	tk0 = t0;
	for (i = 0; i < n; i++)
		yk0[i] = x0[i];
	f(tk0, yk0, k1);
	for (i = 0; i < n; i++)
		k1[i] *= h;
	tk1 = t0 + a2 * h;
	for (i = 0; i < n; i++)
		yk1[i] = x0[i] + b21 * k1[i];
	f(tk1, yk1, k2);
	for (i = 0; i < n; i++)
		k2[i] *= h;
	tk2 = t0 + a3 * h;
	for (i = 0; i < n; i++)
		yk2[i] = x0[i] + b31 * k1[i] + b32 * k2[i];
	f(tk2, yk2, k3);
	for (i = 0; i < n; i++)
		k3[i] *= h;
	tk3 = t0 + a4 * h;
	for (i = 0; i < n; i++)
		yk3[i] = x0[i] + b41 * k1[i] + b42 * k2[i] + b43 * k3[i];
	f(tk3, yk3, k4);
	for (i = 0; i < n; i++)
		k4[i] *= h;
	tk4 = t0 + a5 * h;
	for (i = 0; i < n; i++)
		yk4[i] = x0[i] + b51 * k1[i] + b52 * k2[i] + b53 * k3[i] + b54 * k4[i];
	f(tk4, yk4, k5);
	for (i = 0; i < n; i++)
		k5[i] *= h;
	tk5 = t0 + a6 * h;
	for (i = 0; i < n; i++)
		yk5[i] = x0[i] + b61 * k1[i] + b62 * k2[i] + b63 * k3[i] + b64 * k4[i] + b65 * k5[i];
	f(tk5, yk5, k6);
	for (i = 0; i < n; i++)
		k6[i] *= h;
	tk6 = t0 + a7 * h;
	for (i = 0; i < n; i++)
		yk6[i] = x0[i] + b71 * k1[i] + b73 * k3[i] + b74 * k4[i] + b75 * k5[i] + b76 * k6[i];
	f(tk6, yk6, k7);
	for (i = 0; i < n; i++)
		k7[i] *= h;
	for (i = 0; i < n; i++)
		x1[i] = x0[i] + c1 * k1[i] + c3 * k3[i] + c4 * k4[i] + c5 * k5[i] + c6 * k6[i];

	// e calculated as difference between x1^ and x1
	double
		e1 = 71.0L / 57600,
		e3 = -71.0L / 16695,
		e4 = 71.0L / 1920,
		e5 = -17253.0L / 339200,
		e6 = 22.0L / 525,
		e7 = -1.0L / 40;
	vector<double> e = { e1, 0, e3, e4, e5, e6, e7 };
	k[0] = k1;
	k[1] = k2;
	k[2] = k3;
	k[3] = k4;
	k[4] = k5;
	k[5] = k6;
	{
		// calculate erra
		eps = 0.0;
		for (i = 0; i < n; i++)
		{
			//cout << " " << h << endl;
			double d = 0;
			for (int j = 0; j < 6; j++) {
				d += e[j] * k[j][i];
			}
			//d = yk4[i] - yk3[i];
			eps += pow((h * d) / _max(pow(10.0, -5.0), abs(x1[i]), abs(x0[i]), 2 * u / epsMax), 2.0);
		}
		eps = sqrt(eps / n);
	}
}