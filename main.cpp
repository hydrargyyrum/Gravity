#define _CRT_SECURE_NO_WARNINGS

#include "integrator.h"
#include "integrator_dp54.h"
#include <functional>

#include <string>
#include <fstream>
#include <iostream>
#include <vector>

#include "EarthGravity.h"
#include "utils.h"

using namespace std;


void writeToCsvFile(vector<vector<double>> vecs) {
	ofstream f("f.csv");
	const int n = vecs.size();
	for (int i = 0; i < n; i++)
	{
		std::vector<double> vec = vecs[i];
		const int m = vec.size();
		if (m == 0)
			return;
		f << vec[0];
		for (int j = 1; j < m; j++) {
			f << "," << vec[j];
		}
		f << std::endl;
	}
	f.close();
}

void attractor_show(IMathModel* model, function<vector<double>(vector<double>)> map) {
	attractor_show_args args_attr;
	args_attr.example_to_show = false;
	args_attr.example_continiously_draw = true;
	auto v = model->GetResults();
	cout << "size : " << v.size() << endl;
	vector<vector<double>> points(v.size());
	std::transform(v.begin(), v.end(), points.begin(), map);
	args_attr.points = points;
	args_attr.print_dots = false;
	attractor_show(args_attr);
}

enum class Test
{
	Central,
	Normal,
	Anomal
};

GravityModel* getGravityModel(Test test)
{
	switch (test) {
	case Test::Central:
		return new CentralGravitation();
	case Test::Normal:
		return new NormalGravitation();
	case Test::Anomal:
		return new AnomalGravitation();
	}
	return nullptr;
}

IMathModel* getModel(Test test)
{
	auto x = getGravityModel(test);
	OrbitParameters op;
	op.Alpha = 6500e3;
	/*op.E = 0.5;*/
	op.E = 0.1;
	return new Satellite(op, x);
}

void main() {
	cout << "!! ";
	Vector v(vector<double> {
		1.09987728084142e-5,
			1.37974458075151e-6,
			3.98075410778898
	});
	v.print();


	srand(time(nullptr));
	//Test test = Test::Central;
	//Test test = Test::Normal;
	Test test = Test::Anomal;
	auto gm = getGravityModel(test);
	auto result = gm->getAcceleration(Vector(vector<double>
	{
		0, 0, 1e7,
			//1e7, 0, 0
	}));


	cout << endl << "test: ";
	result.print();
	cout << endl << endl;

	/*cin.get();
	return;*/

	IMathModel* model = getModel(test);

	double t0 = 0;
	double tf = 25000;
	/*double hmin;
	double hmax;*/
	double h0 = tf;
	double eps = 0.1;
	double k = 60;

	//model->doPrintAll = true;
	dp54_args args(
		model,
		t0,
		tf,
		0,
		0,
		h0,
		eps,
		k);
	dp54_integrator integrator(args);
	integrator.integration();

	model->printResults();

	auto results = model->GetResults();
	writeToCsvFile(results);

	cin.get();

	return;
	

	attractor_show(model, [test](vector<double> v) -> vector<double> {
		if (test == Test::Normal)
			return vector<double>
		{

			v[6] / 1000,
				0,
				v[1] / 1000,
		};
		if (test == Test::Anomal)
		{
			return vector<double>
			{

				v[6] / 500,
					0,
					v[0] / 500,
			};
		}
		return vector<double>
		{
			v[6] / 1000,
				0,
				v[0] / 1000,
		};
		});


	std::cin.get();
}