// Swendsen-Wang cluster algorithm for the 2-D Ising Model

#ifndef SWENDSEN_WANG_H
#define SWENDSEN_WANG_H

#include <utility>

class SwendsenWang {
public:
	SwendsenWang();

	void setCoupling(double J) { this->J = J; }
	double coupling() const { return J; }

	void setMagneticField(double H) { this->H = H; }
	double magneticField() const { return H; }

	void setTemperature(double T) { this->T = T; }
	double temperature() const { return T; }

	void setSizeX(int Lx) { this->Lx = Lx; N = Lx * Ly; }
	void setSizeY(int Ly) { this->Ly = Ly; N = Lx * Ly; }

	int sizeX() const { return Lx; }
	int sizeY() const { return Ly; }

	void initialize();
	void initializeClusterVariables();

	void freezeOrMeltBonds();
	int properLabel(int label);
	void labelClusters();
	void flipClusterSpins();

	void oneMonteCarloStep();

	void initializeObservables();
	void measureObservables();
	void computeAverages();

	std::pair<double, double> energy() const {
		return std::make_pair(this->eAve, this->eError);
	}
	std::pair<double, double> energySquare() const {
		return std::make_pair(this->e2Ave, this->e2Error);
	}

	std::pair<double, double> magnet() const {
		return std::make_pair(this->mAve, this->mError); 
	}
	std::pair<double, double> magnetSquare() const {
		return std::make_pair(this->m2Ave, this->mError); 
	}

private:
	int Lx, Ly; // number of spins in x and y
	int N;      // number of spins
	int **s;    // the spins
	double J;   // ferromagnetic coupling
	double T;   // temperature
	double H;   // magnetic field
	int steps;  // steps so far

private:
	bool **iBondFrozen, **jBondFrozen;  // bond lattice - two bonds per spin
	double freezeProbability;           // 1 - e^(-2J/kT)
	int **cluster;                      // cluster labels for spins
	int *labelLabel;                    // to determine proper labels
	bool *sNewChosen;                   // has the new spin value been chosen?
	int *sNew;                          // random new spin values in each cluster

private:
	double eSum;                // accumulator for energy per spin
	double eSqdSum;             // accumulator for square of energy per spin
	double eQuadSum;			// accumulator for quad of energy per spin

	double mSum;
	double mSqdSum;
	double mQuadSum;

	int nSum;                   // number of terms in sum

private:
	double eAve;                // average energy per spin
	double eError;              // Monte Carlo error estimate
	double e2Ave;
	double e2Error;
	double mAve;
	double m2Ave;
	double mError;
	double m2Error;
};

#endif // SWENDSEN_WANG_H