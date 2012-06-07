#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <list>                 // to save values for autocorrelations

#include "mersenne.h"
#include "swendsen-wang.h"

using namespace std;

SwendsenWang::SwendsenWang()
	: J(+1)
	, H(0)
	, steps(0)
{
}

inline double qadran() {
	return mt_random_d();
}

void SwendsenWang::initialize() {
    mt_random_init();

    s = new int* [Lx];
    for (int i = 0; i < Lx; i++)
        s[i] = new int [Ly];
    for (int i = 0; i < Lx; i++)
        for (int j = 0; j < Ly; j++)
            s[i][j] = qadran() < 0.5 ? +1 : -1;   // hot start

    steps = 0;
}

void SwendsenWang::initializeClusterVariables() {
    // allocate 2-D arrays for bonds in x and y directions    
    iBondFrozen = new bool* [Lx];
    jBondFrozen = new bool* [Lx];
    for (int i = 0; i < Lx; i++) {
        iBondFrozen[i] = new bool [Ly];
        jBondFrozen[i] = new bool [Ly];
    }

    // compute the bond freezing probability
    freezeProbability = 1 - exp(-2*J/T);

    // allocate 2-D array for spin cluster labels
    cluster = new int* [Lx];
    for (int i = 0; i < Lx; i++)
        cluster[i] = new int [Ly];

    // allocate arrays of size = number of spins for
    labelLabel = new int [N];        // proper label pointers
    sNewChosen = new bool [N];       // setting new cluster spin values
    sNew = new int [N];              // new cluster spin values
}

void SwendsenWang::oneMonteCarloStep() {
    // first construct a bond lattice with frozen bonds
    freezeOrMeltBonds();

    // use the Hoshen-Kopelman algorithm to identify and label clusters
    labelClusters();

    // re-set cluster spins randomly up or down
    flipClusterSpins();

    ++steps;
}

void SwendsenWang::freezeOrMeltBonds() {
    // visit all the spins in the lattice
    for (int i = 0; i < Lx; i++)
    for (int j = 0; j < Ly; j++) {

        // freeze or melt the two bonds connected to this spin
        // using a criterion which depends on the Boltzmann factor
        iBondFrozen[i][j] = jBondFrozen[i][j] = false;

        // bond in the i direction
        int iNext = i == Lx-1 ? 0 : i+1;
        if (s[i][j] == s[iNext][j] && qadran() < freezeProbability)
            iBondFrozen[i][j] = true;

        // bond in the j direction
        int jNext = j == Ly-1 ? 0 : j+1;
        if (s[i][j] == s[i][jNext] && qadran() < freezeProbability)
            jBondFrozen[i][j] = true;
    }
}

int SwendsenWang::properLabel(int label) {
    while (labelLabel[label] != label)
        label = labelLabel[label];
    return label;
}

void SwendsenWang::labelClusters() {
    int label = 0;

    // visit all lattice sites
    for (int i = 0; i < Lx; i++)
    for (int j = 0; j < Ly; j++) {

        // find previously visited sites connected to i,j by frozen bonds
        int bonds = 0;
        int iBond[4], jBond[4];

        // check bond to i-1,j
        if (i > 0 && iBondFrozen[i - 1][j]) {
            iBond[bonds] = i - 1;
            jBond[bonds++] = j;
        }

        // apply periodic conditions at the boundary:
        // if i,j is the last site, check bond to i+1,j
        if (i == Lx - 1 && iBondFrozen[i][j]) {
            iBond[bonds] = 0;
            jBond[bonds++] = j;
        }

        // check bond to i,j-1
        if (j > 0 && jBondFrozen[i][j - 1]) {
            iBond[bonds] = i;
            jBond[bonds++] = j - 1;
        }

        // periodic boundary conditions at the last site
        if (j == Ly - 1 && jBondFrozen[i][j]) {
            iBond[bonds] = i;
            jBond[bonds++] = 0;
        }

        // check number of bonds to previously visited sites
        if (bonds == 0) { // need to start a new cluster
            cluster[i][j] = label;
            labelLabel[label] = label;
            ++label;
        } else {          // re-label bonded spins with smallest proper label
            int minLabel = label;
            for (int b = 0; b < bonds; b++) {
                int pLabel = properLabel(cluster[iBond[b]][jBond[b]]);
                if (minLabel > pLabel)
                    minLabel = pLabel;
            }

            // set current site label to smallest proper label
            cluster[i][j] = minLabel;

            // re-set the proper label links on the previous labels
            for (int b = 0; b < bonds; b++) {
                int pLabel = cluster[iBond[b]][jBond[b]];
                labelLabel[pLabel] = minLabel;

                // re-set label on connected sites
                cluster[iBond[b]][jBond[b]] = minLabel;
            }
        }
    }
}

void SwendsenWang::flipClusterSpins() {
    for (int i = 0; i < Lx; i++)
    for (int j = 0; j < Ly; j++) {

        // random new cluster spins values have not been set
        int n = i * Lx + j;
        sNewChosen[n] = false;

        // replace all labels by their proper values
        cluster[i][j] = properLabel(cluster[i][j]);
    }    

    int flips = 0;    // to count number of spins that are flipped
    for (int i = 0; i < Lx; i++)
    for (int j = 0; j < Ly; j++) {

        // find the now proper label of the cluster
        int label = cluster[i][j];

        // choose a random new spin value for cluster
        // only if this has not already been done
        if (!sNewChosen[label]) {    
            sNew[label] = qadran() < 0.5 ? +1 : -1;
            sNewChosen[label] = true;
        }

        // re-set the spin value and count number of flips
        if (s[i][j] != sNew[label]) {
            s[i][j] = sNew[label];
            ++flips;
        }
    }
}

void SwendsenWang::initializeObservables() {
    eSum = eSqdSum = 0;     // zero energy accumulators
    nSum = 0;               // no terms so far
    mSum = mSqdSum = 0;
}

void SwendsenWang::measureObservables() {
    int sSum = 0, ssSum = 0;

    for (int i = 0; i < Lx; i++)
    for (int j = 0; j < Ly; j++) {
        sSum += s[i][j];
        int iNext = i == Lx-1 ? 0 : i+1;
        int jNext = j == Ly-1 ? 0 : j+1;
        ssSum += s[i][j]*(s[iNext][j] + s[i][jNext]);
    }

    double e = -(J*ssSum + H*sSum)/ (double)N;
    double m = fabs((double)sSum / (double)N);
	
    eSum += e;
    eSqdSum += e * e;
	eQuadSum += e * e * e * e;
	
    mSum += m;
    mSqdSum += m * m;
	mQuadSum += m * m * m * m;

    ++nSum;
}

void SwendsenWang::computeAverages() {
    eAve = eSum / nSum;
    eError = eSqdSum / nSum;
    eError = sqrt(eError - eAve*eAve);
    eError /= sqrt(double(nSum));
    
    e2Ave = eSqdSum / nSum;
	e2Error = eQuadSum / nSum;
	e2Error = sqrt(e2Error - e2Ave*e2Ave);
	e2Error /= sqrt(double(nSum));
	
    mAve = mSum / nSum;
    mError = mSqdSum / nSum;
    mError = sqrt(mError - mAve*mAve);
    mError /= sqrt(double(nSum));
    
    m2Ave = mSqdSum / nSum;
	m2Error = mQuadSum / nSum;
	m2Error = sqrt(m2Error - m2Ave*m2Ave);
	m2Error /= sqrt(double(nSum));
}