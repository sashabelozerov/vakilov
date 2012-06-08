#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>

#include "mersenne.h"
#include "swendsen-wang.h"
#include "csvwriter.h"

using namespace std;

int main(int argc, char** argv) {
	int L;                // lattice linear size
	int d;                // number of dimenstions
	double T_min;         // temperature range start
	double T_max;         // temperature range end
	double T_step;        // temperature step
	int MCSteps;          // number of Monte Carlo steps
	std::string outFile;  // output file name

	map<string, string> opts;
	for (int i = 1; i < argc; i++) {
		string arg(argv[i]);
		string::size_type delim_pos = arg.find('=');
		if (delim_pos != string::npos) {
			string name(arg.begin(), arg.begin() + delim_pos);
			string value(arg.begin() + delim_pos + 1, arg.end());
			opts.insert(make_pair(name, value));
		}
	}

	map<string, string>::const_iterator opt;

	if ((opt = opts.find("L")) != opts.end()) {
		L = L = atoi(opt->second.c_str());
	} else {
		cerr << "L is not set" << endl;
		exit(EXIT_FAILURE);
	}
	if ((opt = opts.find("Tmin")) != opts.end()) {
		T_min = atof(opt->second.c_str());
	} else {
		cerr << "Tmin is not set" << endl;
		exit(EXIT_FAILURE);
	}
	if ((opt = opts.find("Tmax")) != opts.end()) {
		T_max = atof(opt->second.c_str());
	} else {
		cerr << "Tmax is not set" << endl;
		exit(EXIT_FAILURE);
	}
	if ((opt = opts.find("Tstep")) != opts.end()) {
		T_step = atof(opt->second.c_str());
	} else {
		cerr << "Tstep is not set" << endl;
		exit(EXIT_FAILURE);
	}
	if ((opt = opts.find("mcs")) != opts.end()) {
		MCSteps = atoi(opt->second.c_str());
	} else {
		cerr << "mcs is not set" << endl;
		exit(EXIT_FAILURE);
	}

	SwendsenWang *sw;
	if ((opt = opts.find("d")) != opts.end()) {
		d = atoi(opt->second.c_str());
		if (d == 2) {
			sw = new SwendsenWang2D;
		} else if (d == 3) {
			sw = new SwendsenWang3D;
		} else {
			cerr << "d can be either 2 or 3" << endl;
			exit(EXIT_FAILURE);
		}
	} else {
		cerr << "d is not set" << endl;
		exit(EXIT_FAILURE);
	}

	cout << d << "-dimensional Ising Model - Swendsen-Wang Algorithm\n"
	     << "-----------------------------------------------------"
		 << endl;

	if ((opt = opts.find("outfile")) != opts.end()) {
		outFile = opt->second;
	} else {
		stringstream name;
		name << d << "d_" << L << "_" << T_min << "-" << T_max << "-" << T_step << "_" << MCSteps << ".txt";
		outFile = name.str();
	}

	CSVWriter writer(outFile);
	sw->setLatticeSize(L);

	for (double T = T_min; T <= T_max; T += T_step) {
		sw->setTemperature(T);

		sw->initialize();
		sw->initializeClusterVariables();

		int thermSteps = MCSteps / 5;
		for (int i = 0; i < thermSteps; i++)
			sw->oneMonteCarloStep();

		sw->initializeObservables();
		for (int i = 0; i < MCSteps; i++) {
			sw->oneMonteCarloStep();
			sw->measureObservables();
		}

		sw->computeAverages();

		double c = (sw->energySquare().first - pow(sw->energy().first, 2)) / (T * T);
		double cError = (pow(sw->energy().second, 2) + sw->energySquare().second) / (T * T);

		double x = (sw->magnetSquare().first - pow(sw->magnet().first, 2)) / T;
		double xError = (2 * pow(sw->magnet().second, 2)) / T; 

		double percent_complete = (T - T_min) / (T_max - T_min) * 100;
		cout << "\r" << setprecision(3) << percent_complete << "% complete   ";
		cout.flush();
		
		writer.write(T,
			sw->energy().first,
			sw->energy().second,
			sw->energySquare().first,
			sw->energySquare().second,
			sw->magnet().first,
			sw->magnet().second,
			sw->magnetSquare().first,
			sw->magnetSquare().second,
			c,
			cError,
			x,
			xError
		);
	}

	cout << endl;

	exit(EXIT_SUCCESS);
}


