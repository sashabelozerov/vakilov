/* 
 * File:   main.cpp
 * Author: wert
 *
 * Created on June 6, 2012, 10:23 PM
 */

#include <cstdlib>
#include "mersenne.h"
#include "swendsen-wang.h"

using namespace std;

/*
 * 
 */

int main(int argc, char** argv) {

    float T_min, T_max, T_step;
    cout << " Two-dimensional Ising Model - Swendsen-Wang Algorithm\n"
         << " -----------------------------------------------------\n"
         << " Enter number of spins L in each direction: ";
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
    int MCSteps;
    cin >> MCSteps;

    T = T_min;
    while(T <= T_max) {
        initialize();
        initializeClusterVariables();

        int thermSteps = MCSteps / 5;
        cout << " Performing " << thermSteps 
            << " thermalization steps ..." << flush;
        for (int i = 0; i < thermSteps; i++)
            oneMonteCarloStep();
        cout << " done\n Performing production steps ..." << flush;

        initializeObservables();
        for (int i = 0; i < MCSteps; i++) {
            oneMonteCarloStep();
            measureObservables();
        }
        cout << " done" << endl;
        computeAverages();
        cout << "T = " << T << " | Energy per spin = " << eAve << " +- " << eError << endl;
        T += T_step;
    }
}


