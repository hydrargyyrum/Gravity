#pragma once

#define _CRT_SECURE_NO_WARNINGS

#include "integrator.h"
#include "integrator_dp54.h"
#include <functional>

#include <string>
#include <iostream>
#include <vector>

#include "EarthGravity.h"

using namespace std;




IMathModel* test_runge() {
	//gauss G(-1, 4);
	vector<double> m = { 10, 10, 10, 3 };

	//MotionModel* model = new MotionModel(m, &G);

	attractor* modelA = new attractor(m);
	/*eylerArgs args1(modelA, 10, 0.01, 0.01, 0);
	eyler c1(args1); */

	/*eylerCorrectionArgs args2(model, 100, 0.01, 1, 0);
	eylerCorrection c1(args2); */

	rungeArgs args(modelA, 1000, 0.01);
	runge c1(args);

	c1.integration();
	//cin.get();

	return modelA;
}

IMathModel* test_dp54() {
	double
		t0 = 0,
		tf = 100,
		h0 = 0.001,
		epsMax = 0.0000001,
		k = 1;
	vector<double> m = { 0.994, 0, 0, -2.00158510637908252240537862224 };
	IMathModel* modelA = new dp54_Model(m);
	dp54_args args(modelA, t0, tf, 0, 0, h0, epsMax, k);
	dp54_integrator c1(args);
	c1.integration();
	return modelA;
}

IMathModel* test_dp54_pendulum() {
	double
		t0 = 0,
		tf = 160,
		h0 = 0.001,
		epsMax = 0.00001,
		k = 0.01;
	//IMathModel* modelA = new pendulum_ideal_Model(t0);    //1
	//IMathModel* modelA = new pendulum_ideal_spring_Model(t0);    //2
	//IMathModel* modelA = new pendulum_real_Model(t0);     //3
	//IMathModel* modelA = new pendulum_real_spring_slide_Model(t0);    //4
	IMathModel* modelA = new pendulum_real_viscous_Model(t0);   //5
	dp54_args args(modelA, t0, tf, 0, 0, h0, epsMax, k);
	dp54_integrator c1(args); // 1 2  for ideal
	//dp54_integrator_2 c1(args); // 3 4 5 for real
	c1.integration();
	return modelA;
}


IMathModel* test_dp54_EarthMovementModel() {
	const double secInDay = 86400;
	int
		dd = 1,
		mm = 1,
		yyyy = 2022;

	double
		lat = 54.43,
		lng = 20.30,
		utc = +2,
		//utc = +0,

		/*x0 = 55.45,
		y0 = 37.37,
		utc = +3,*/


		t0 = 0 - utc * 3600,
		tf = secInDay * 100 - utc * 3600,

		//tf = 31536000.0,
		h0 = tf,

		/*h0 = 31536000.0,*/

		/*tf = 86400,
		h0 = 86400,*/
		/*epsMax = 1e-3,*/
		/*epsMax = 1,*/
		epsMax = 0.0001,
		k = 60,
		x = -2.741248720132488e7,
		y = 1.332710216571607e8,
		z = 5.780072873943681e7,
		vx = -2.980512290685886e1,
		vy = -4.971761443577496,
		vz = -2.154845308279276;


	//IMathModel* model1 = new EarthMovementModel(t0, x, y, z, vx, vy, vz);
	//dp54_args args(model1, t0, tf, 0, 0, h0, epsMax, k);
	//dp54_integrator c1(args);
	//c1.integration();
	//auto v = model1->X0();
	/*SunModel* model2 = new SunModel(
		43.07, 1337.37, 2459659.5,
		v[0], v[1], v[2], v[3], v[4], v[5], v[6]);*/

	if (false)
	{
		// task 1. shadow
		SunModel* model2 = new SunModel(
			lat, lng,
			dd, mm, yyyy,
			t0, tf,
			-1.503101968223290 * 1e+10,
			2.388181224842258 * 1e+5,
			1.344286818174785 * 1e+5,
			-4.633867856978692 * 1e-1,
			-2.744939467460866 * 1e+1,
			-1.189889265710981 * 1e+1
		);
		dp54_args args(model2, t0, tf, 0, 0, h0, epsMax, k);
		dp54_integrator c2(args);
		c2.integration();
		return model2;
	}
	else
	{
		// task 2. day langth
		IMathModel* model2 = new DayLengthModel(
			lat, lng,
			dd, mm, yyyy,
			t0, tf, utc,
			-1.503101968223290 * 1e+10,
			2.388181224842258 * 1e+5,
			1.344286818174785 * 1e+5,
			-4.633867856978692 * 1e-1,
			-2.744939467460866 * 1e+1,
			-1.189889265710981 * 1e+1
		);
		dp54_args args(model2, t0, tf, 0, 0, h0, epsMax, k);
		dp54_integrator c2(args);
		c2.integration();
		return model2;
	}
}

void attractor_show(IMathModel* model, function<vector<double>(vector<double>)> map) {
	attractor_show_args args_attr;
	args_attr.example_to_show = false;
	args_attr.example_continiously_draw = true;
	auto v = model->X_all();
	cout << "size : " << v.size() << endl;
	vector<vector<double>> points(v.size());
	std::transform(v.begin(), v.end(), points.begin(), map);
	args_attr.points = points;
	args_attr.print_dots = false;
	attractor_show(args_attr);
}

void attractor_show_shadow(IMathModel* model, function<vector<double>(vector<double>)> map) {
	attractor_show_args args_attr;
	args_attr.example_to_show = false;
	args_attr.example_continiously_draw = true;
	auto v = model->GetResults();
	cout << "size : " << v.size() << endl;
	//return;


	vector<vector<double>> points(v.size());
	std::transform(v.begin(), v.end(), points.begin(), map);
	args_attr.points = points;
	args_attr.print_dots = false;
	attractor_show(args_attr);
}



void mainOld() {
	srand(time(nullptr));
	//attractor_show(test_runge()); // todo coefficiens * 40 to remove

	// for DP
	/*attractor_show(test_dp54(), [](vector<double> v) -> vector<double> {
		return vector<double>{ v[0] * 40, 0, v[1] * 40 };
	});*/


	// for pendulum

	//attractor_show(test_dp54_pendulum(), [](vector<double> v) -> vector<double> {
	//attractor_show(test_dp54_EarthMovementModel(), [](vector<double> v) -> vector<double> {
	//	//return vector<double>{ v[0] * 3, 0, v[1] * 3 };                             // 1 2 3 4 
	//	return vector<double>{ v[1] * 1e-7, 0, v[2] * 1e-7 };                         // for pendulum_real_viscous_Model
	//	});

	attractor_show_shadow(test_dp54_EarthMovementModel(), [](vector<double> v) -> vector<double> {
		//	//return vector<double>{ v[0] * 3, 0, v[1] * 3 };                             // 1 2 3 4 
			//return vector<double>{ v[0] * 1e-7, 0, v[2] * 1e-7 };                         // for pendulum_real_viscous_Model

		return vector<double>{ v[0] * 20, v[1] * 20, v[2] * 20 }; // task 1
		//return vector<double>{ v[0] * 20, 0, v[1] * 20  }; // task 2
		});

	std::cin.get();
}