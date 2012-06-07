/* 
 * File:   main.cpp
 * Author: wert
 *
 * Created on June 6, 2012, 10:23 PM
 */

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
	int L;
	double T_min, T_max, T_step;
	int MCSteps;
	std::string outFile;

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

	cout << " Two-dimensional Ising Model - Swendsen-Wang Algorithm\n"
		 << " -----------------------------------------------------\n";

	if (opts["cli"] == "yes") {
		map<string, string>::const_iterator iterator;
		if ((iterator = opts.find("L")) != opts.end()) {
			L = L = atoi(iterator->second.c_str());
		} else {
			cerr << "L is not set" << endl;
			exit(1);
		}
		if ((iterator = opts.find("Tmin")) != opts.end()) {
			T_min = atof(iterator->second.c_str());
		} else {
			cerr << "Tmin is not set" << endl;
			exit(1);
		}
		if ((iterator = opts.find("Tmax")) != opts.end()) {
			T_max = atof(iterator->second.c_str());
		} else {
			cerr << "Tmax is not set" << endl;
			exit(1);
		}
		if ((iterator = opts.find("Tstep")) != opts.end()) {
			T_step = atof(iterator->second.c_str());
		} else {
			cerr << "Tstep is not set" << endl;
			exit(1);
		}
		if ((iterator = opts.find("mcs")) != opts.end()) {
			MCSteps = atoi(iterator->second.c_str());
		} else {
			cerr << "mcs is not set" << endl;
			exit(1);
		}
		if ((iterator = opts.find("outfile")) != opts.end()) {
			outFile = iterator->second;
		} else {
			stringstream name;
			name << L << "_" << T_min << "-" << T_max << "-" << T_step << "_" << MCSteps << ".txt";
			outFile = name.str();
		}
	} else {
		cout << " Enter number of spins L in each direction: ";
		cin >> L;
		cout << " Enter min T: ";
		cin >> T_min;
		cout << " Enter max T: ";
		cin >> T_max;
		cout << " Enter T step: ";
		cin >> T_step;
		cout << " Enter number of Monte Carlo steps: ";
		cin >> MCSteps;
		cout << " Enter output filename: ";
		cin >> outFile;
	}

	CSVWriter writer(outFile);
	
	SwendsenWang sw;

	sw.setSizeX(L);
	sw.setSizeY(L);

	for (double T = T_min; T <= T_max; T += T_step) {
		sw.setTemperature(T);

		sw.initialize();
		sw.initializeClusterVariables();

		int thermSteps = MCSteps / 5;
		for (int i = 0; i < thermSteps; i++)
			sw.oneMonteCarloStep();

		sw.initializeObservables();
		for (int i = 0; i < MCSteps; i++) {
			sw.oneMonteCarloStep();
			sw.measureObservables();
		}

		sw.computeAverages();

		double c = (sw.energySquare().first - pow(sw.energy().first, 2)) / (T * T);
		double cError = (pow(sw.energy().second, 2) + sw.energySquare().second) / (T * T);

		double x = (sw.magnetSquare().first - pow(sw.magnet().first, 2)) / T;
		double xError = (2 * pow(sw.magnet().second, 2)) / T; 

		double percent_complete = (T - T_min) / (T_max - T_min) * 100;
		cout << "\r" << setprecision(3) << percent_complete << "% complete  ";
		cout.flush();
		
		writer.write(T,
			sw.energy().first,
			sw.energy().second,
			sw.energySquare().first,
			sw.energySquare().second,
			sw.magnet().first,
			sw.magnet().second,
			sw.magnetSquare().first,
			sw.magnetSquare().second,
			c,
			cError,
			x,
			xError
		);
	}

	cout << endl;
}


