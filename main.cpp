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
			Lx = Ly = atoi(iterator->second.c_str());
			N = Lx * Ly; 
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
			name << Lx << "_" << T_min << "-" << T_max << "-" << T_step << "_" << MCSteps << ".txt";
			outFile = name.str();
		}
	} else {
		cout << " Enter number of spins L in each direction: ";
		cin >> Lx;
		Ly = Lx;
		N = Lx * Ly;
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
	
    T = T_min;
    while(T <= T_max) {
        initialize();
        initializeClusterVariables();

        int thermSteps = MCSteps / 5;
        for (int i = 0; i < thermSteps; i++)
            oneMonteCarloStep();

        initializeObservables();
        for (int i = 0; i < MCSteps; i++) {
            oneMonteCarloStep();
            measureObservables();
        }

        computeAverages();

		double c = (e2Ave - eAve*eAve) / T*T;
		double cError = (eError*eError + e2Error) / T*T;

		double x = (m2Ave - mAve*mAve) / T;
		double xError = (2 * mError*mError) / T; 

		double percent_complete = (T - T_min) / (T_max - T_min) * 100;
		cout << "\r" << setprecision(3) << percent_complete << "% complete ";
		cout.flush();
		
		writer.write(T, eAve, eError, e2Ave, e2Error, mAve, mError, m2Ave, mError, c, cError, x, xError);
        T += T_step;
    }

	cout << endl;
}


