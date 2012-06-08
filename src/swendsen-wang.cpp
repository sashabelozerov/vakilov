#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <list>                 // to save values for autocorrelations

#include "mersenne.h"
#include "swendsen-wang.h"

using namespace std;

// ====================================================================================
//                                     SwendsenWang
// ====================================================================================

SwendsenWang::SwendsenWang()
	: J(+1)
	, H(0)
	, steps(0)
{
}

SwendsenWang::~SwendsenWang() {
}

inline double qadran() {
	return mt_random_d();
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

int SwendsenWang::properLabel(int label) {
    while (labelLabel[label] != label)
        label = labelLabel[label];
    return label;
}

void SwendsenWang::initializeObservables() {
    eSum = eSqdSum = 0;     // zero energy accumulators
    nSum = 0;               // no terms so far
    mSum = mSqdSum = 0;
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

// ====================================================================================
//                                     SwendsenWang2D
// ====================================================================================

SwendsenWang2D::SwendsenWang2D()
	: SwendsenWang()
{
}

SwendsenWang2D::~SwendsenWang2D() {
	free();
}

void SwendsenWang2D::initialize() {
    mt_random_init();

    s = new int* [L];
    for (int i = 0; i < L; i++)
        s[i] = new int [L];
    for (int i = 0; i < L; i++)
        for (int j = 0; j < L; j++)
            s[i][j] = qadran() < 0.5 ? +1 : -1;   // hot start

    steps = 0;
}

void SwendsenWang2D::initializeClusterVariables() {
    // allocate 2-D arrays for bonds in x and y directions    
    iBondFrozen = new bool* [L];
    jBondFrozen = new bool* [L];
    for (int i = 0; i < L; i++) {
        iBondFrozen[i] = new bool [L];
        jBondFrozen[i] = new bool [L];
    }

    // compute the bond freezing probability
    freezeProbability = 1 - exp(-2*J/T);

    // allocate 2-D array for spin cluster labels
    cluster = new int* [L];
    for (int i = 0; i < L; i++)
        cluster[i] = new int [L];

    // allocate arrays of size = number of spins for
    labelLabel = new int [N];        // proper label pointers
    sNewChosen = new bool [N];       // setting new cluster spin values
    sNew = new int [N];              // new cluster spin values
}

void SwendsenWang2D::freezeOrMeltBonds() {
    // visit all the spins in the lattice
    for (int i = 0; i < L; i++)
    for (int j = 0; j < L; j++) {

        // freeze or melt the two bonds connected to this spin
        // using a criterion which depends on the Boltzmann factor
        iBondFrozen[i][j] = jBondFrozen[i][j] = false;

        // bond in the i direction
        int iNext = i == L-1 ? 0 : i+1;
        if (s[i][j] == s[iNext][j] && qadran() < freezeProbability)
            iBondFrozen[i][j] = true;

        // bond in the j direction
        int jNext = j == L-1 ? 0 : j+1;
        if (s[i][j] == s[i][jNext] && qadran() < freezeProbability)
            jBondFrozen[i][j] = true;
    }
}

void SwendsenWang2D::labelClusters() {
    int label = 0;

    // visit all lattice sites
    for (int i = 0; i < L; i++)
    for (int j = 0; j < L; j++) {

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
        if (i == L - 1 && iBondFrozen[i][j]) {
            iBond[bonds] = 0;
            jBond[bonds++] = j;
        }

        // check bond to i,j-1
        if (j > 0 && jBondFrozen[i][j - 1]) {
            iBond[bonds] = i;
            jBond[bonds++] = j - 1;
        }

        // periodic boundary conditions at the last site
        if (j == L - 1 && jBondFrozen[i][j]) {
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

void SwendsenWang2D::flipClusterSpins() {
    for (int i = 0; i < L; i++)
    for (int j = 0; j < L; j++) {

        // random new cluster spins values have not been set
        int n = i * L + j;
        sNewChosen[n] = false;

        // replace all labels by their proper values
        cluster[i][j] = properLabel(cluster[i][j]);
    }    

    int flips = 0;    // to count number of spins that are flipped
    for (int i = 0; i < L; i++)
    for (int j = 0; j < L; j++) {

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

void SwendsenWang2D::measureObservables() {
    int sSum = 0, ssSum = 0;

    for (int i = 0; i < L; i++)
    for (int j = 0; j < L; j++) {
        sSum += s[i][j];
        int iNext = i == L-1 ? 0 : i+1;
        int jNext = j == L-1 ? 0 : j+1;
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

void SwendsenWang2D::free() {
	for (int i = 0; i < L; i++) {
		delete[] s[i];
	}
	delete[] s;

    for (int i = 0; i < L; i++) {
        delete[] iBondFrozen[i];
        delete[] jBondFrozen[i];
    }
    delete[] iBondFrozen;
    delete[] jBondFrozen;

    for (int i = 0; i < L; i++) {
		delete[] cluster[i];
	}
	delete[] cluster;

    delete[] labelLabel;
    delete[] sNewChosen;
    delete[] sNew;
}

// ====================================================================================
//                                     SwendsenWang3D
// ====================================================================================

SwendsenWang3D::SwendsenWang3D()
	: SwendsenWang()
{
}

SwendsenWang3D::~SwendsenWang3D() {
	free();
}

void SwendsenWang3D::initialize() {
    mt_random_init();

    s = new int** [L];
    for (int i = 0; i < L; i++) {
        s[i] = new int* [L];
        for (int j = 0; j < L; j++) {
			s[i][j] = new int [L];
			for (int k = 0; k < L; k++) {
				s[i][j][k] = qadran() < 0.5 ? +1 : -1;   // hot start
			}
		}
	}
    steps = 0;
}

void SwendsenWang3D::initializeClusterVariables() {
    // allocate 2-D arrays for bonds in x and y directions    
    iBondFrozen = new bool** [L];
    jBondFrozen = new bool** [L];
	kBondFrozen = new bool** [L];
    for (int i = 0; i < L; i++) {
        iBondFrozen[i] = new bool* [L];
        jBondFrozen[i] = new bool* [L];
		kBondFrozen[i] = new bool* [L];
		for (int j = 0; j < L; j++) {
			iBondFrozen[i][j] = new bool [L];
			jBondFrozen[i][j] = new bool [L];
			kBondFrozen[i][j] = new bool [L];
		}
    }

    // compute the bond freezing probability
    freezeProbability = 1 - exp(-2*J/T);

    // allocate 2-D array for spin cluster labels
    cluster = new int** [L];
    for (int i = 0; i < L; i++) {
        cluster[i] = new int* [L];
		for (int j = 0; j < L; j++) {
			cluster[i][j] = new int [L];
		}
	}

    // allocate arrays of size = number of spins for
    labelLabel = new int [N];        // proper label pointers
    sNewChosen = new bool [N];       // setting new cluster spin values
    sNew = new int [N];              // new cluster spin values
}

void SwendsenWang3D::freezeOrMeltBonds() {
    // visit all the spins in the lattice
    for (int i = 0; i < L; i++)
    for (int j = 0; j < L; j++)
	for (int k = 0; k < L; k++) {

        // freeze or melt the two bonds connected to this spin
        // using a criterion which depends on the Boltzmann factor
        iBondFrozen[i][j][k] = jBondFrozen[i][j][k] = kBondFrozen[i][j][k] = false;

        // bond in the i direction
        int iNext = i == L-1 ? 0 : i+1;
        if (s[i][j][k] == s[iNext][j][k] && qadran() < freezeProbability)
            iBondFrozen[i][j][k] = true;

        // bond in the j direction
        int jNext = j == L-1 ? 0 : j+1;
        if (s[i][j][k] == s[i][jNext][k] && qadran() < freezeProbability)
            jBondFrozen[i][j][k] = true;

        // bond in the k direction
        int kNext = k == L-1 ? 0 : k+1;
        if (s[i][j][k] == s[i][j][kNext] && qadran() < freezeProbability)
            kBondFrozen[i][j][k] = true;
    }
}

void SwendsenWang3D::labelClusters() {
    int label = 0;

    // visit all lattice sites
    for (int i = 0; i < L; i++)
    for (int j = 0; j < L; j++)
	for (int k = 0; k < L; k++) {

        // find previously visited sites connected to i,j,k by frozen bonds
        int bonds = 0;
        int iBond[6], jBond[6], kBond[6];

        // check bond to i-1,j,k
        if (i > 0 && iBondFrozen[i - 1][j][k]) {
            iBond[bonds] = i - 1;
            jBond[bonds] = j;
			kBond[bonds] = k;
			++bonds;
        }

        // apply periodic conditions at the boundary:
        // if i,j,k is the last site, check bond to i+1,j,k
        if (i == L - 1 && iBondFrozen[i][j][k]) {
            iBond[bonds] = 0;
            jBond[bonds] = j;
			kBond[bonds] = k;
			++bonds;
        }

        // check bond to i,j-1,k
        if (j > 0 && jBondFrozen[i][j - 1][k]) {
            iBond[bonds] = i;
            jBond[bonds] = j - 1;
			kBond[bonds] = k;
			++bonds;
        }

        // periodic boundary conditions at the last site
        if (j == L - 1 && jBondFrozen[i][j][k]) {
            iBond[bonds] = i;
            jBond[bonds] = 0;
			kBond[bonds] = k;
			++bonds;
        }

        // check bond to i,j,k-1
        if (k > 0 && kBondFrozen[i][j][k - 1]) {
            iBond[bonds] = i;
			jBond[bonds] = j;
            kBond[bonds] = k - 1;
			++bonds;
        }

        // periodic boundary conditions at the last site
        if (k == L - 1 && kBondFrozen[i][j][k]) {
            iBond[bonds] = i;
			jBond[bonds] = j;
            kBond[bonds] = 0;
			++bonds;
        }

        // check number of bonds to previously visited sites
        if (bonds == 0) { // need to start a new cluster
            cluster[i][j][k] = label;
            labelLabel[label] = label;
            ++label;
        } else {          // re-label bonded spins with smallest proper label
            int minLabel = label;
            for (int b = 0; b < bonds; b++) {
                int pLabel = properLabel(cluster[iBond[b]][jBond[b]][kBond[b]]);
                if (minLabel > pLabel)
                    minLabel = pLabel;
            }

            // set current site label to smallest proper label
            cluster[i][j][k] = minLabel;

            // re-set the proper label links on the previous labels
            for (int b = 0; b < bonds; b++) {
                int pLabel = cluster[iBond[b]][jBond[b]][kBond[b]];
                labelLabel[pLabel] = minLabel;

                // re-set label on connected sites
                cluster[iBond[b]][jBond[b]][kBond[b]] = minLabel;
            }
        }
    }
}

void SwendsenWang3D::flipClusterSpins() {
    for (int i = 0; i < L; i++)
    for (int j = 0; j < L; j++)
	for (int k = 0; k < L; k++) {
        // random new cluster spins values have not been set
        int n = (i * L * L) + (j * L) + k;
        sNewChosen[n] = false;

        // replace all labels by their proper values
        cluster[i][j][k] = properLabel(cluster[i][j][k]);
    }    

    int flips = 0;    // to count number of spins that are flipped
    for (int i = 0; i < L; i++)
    for (int j = 0; j < L; j++)
	for (int k = 0; k < L; k++) {

        // find the now proper label of the cluster
        int label = cluster[i][j][k];

        // choose a random new spin value for cluster
        // only if this has not already been done
        if (!sNewChosen[label]) {    
            sNew[label] = qadran() < 0.5 ? +1 : -1;
            sNewChosen[label] = true;
        }

        // re-set the spin value and count number of flips
        if (s[i][j][k] != sNew[label]) {
            s[i][j][k] = sNew[label];
            ++flips;
        }
    }
}

void SwendsenWang3D::measureObservables() {
    int sSum = 0, ssSum = 0;

    for (int i = 0; i < L; i++)
    for (int j = 0; j < L; j++)
	for (int k = 0; k < L; k++) {
        sSum += s[i][j][k];
        int iNext = i == L-1 ? 0 : i+1;
        int jNext = j == L-1 ? 0 : j+1;
		int kNext = k == L-1 ? 0 : k+1;
        ssSum += s[i][j][k]*(s[iNext][j][k] + s[i][jNext][k] + s[i][j][kNext]);
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

void SwendsenWang3D::free() {
	for (int i = 0; i < L; i++) {
		for (int j = 0; j < L; j++) {
			delete[] s[i][j];
		}
		delete[] s[i];
	}
	delete[] s;

    for (int i = 0; i < L; i++) {
		for (int j = 0; j < L; j++) {
			delete[] iBondFrozen[i][j];
			delete[] jBondFrozen[i][j];
			delete[] kBondFrozen[i][j];
		}
        delete[] iBondFrozen[i];
        delete[] jBondFrozen[i];
		delete[] kBondFrozen[i];
    }
    delete[] iBondFrozen;
    delete[] jBondFrozen;
	delete[] kBondFrozen;

    for (int i = 0; i < L; i++) {
		for (int j = 0; j < L; j++) {
			delete[] cluster[i][j];
		}
		delete[] cluster[i];
	}
	delete[] cluster;

    delete[] labelLabel;
    delete[] sNewChosen;
    delete[] sNew;
}