// Swendsen-Wang cluster algorithm for the 2-D Ising Model

#ifndef SWENDSEN_WANG_H
#define SWENDSEN_WANG_H

#include <utility>

class SwendsenWang {
public:
	SwendsenWang();
	virtual ~SwendsenWang();

	void setCoupling(double J) { this->J = J; }
	double coupling() const { return J; }

	void setMagneticField(double H) { this->H = H; }
	double magneticField() const { return H; }

	void setTemperature(double T) { this->T = T; }
	double temperature() const { return T; }

	void setLatticeSize(int L) { this->L = L; N = L * L * L; }
	int latticeSize() const { return L; }

	virtual void initialize() = 0;
	virtual void initializeClusterVariables() = 0;

	virtual void freezeOrMeltBonds() = 0;
	int properLabel(int label);
	virtual void labelClusters() = 0;
	virtual void flipClusterSpins() = 0;

	void oneMonteCarloStep();

	void initializeObservables();
	virtual void measureObservables() = 0;
	void computeAverages();

	std::pair<double, double> energy() const {
		return std::make_pair(eAve, eError);
	}
	std::pair<double, double> energySquare() const {
		return std::make_pair(e2Ave, e2Error);
	}

	std::pair<double, double> magnet() const {
		return std::make_pair(mAve, mError); 
	}
	std::pair<double, double> magnetSquare() const {
		return std::make_pair(m2Ave, mError); 
	}

	virtual void free() = 0;

protected:
	int L;     // number of spins in x and y
	int N;     // number of spins
	double J;  // ferromagnetic coupling
	double T;  // temperature
	double H;  // magnetic field
	int steps; // steps so far

protected:
	int *labelLabel;            // to determine proper labels
	bool *sNewChosen;           // has the new spin value been chosen?
	int *sNew;                  // random new spin values in each cluster

protected:
	double eSum;                // accumulator for energy per spin
	double eSqdSum;             // accumulator for square of energy per spin
	double eQuadSum;			// accumulator for quad of energy per spin

	double mSum;
	double mSqdSum;
	double mQuadSum;

	int nSum;                   // number of terms in sum

protected:
	double eAve;                // average energy per spin
	double eError;              // Monte Carlo error estimate
	double e2Ave;
	double e2Error;
	double mAve;
	double m2Ave;
	double mError;
	double m2Error;

};

class SwendsenWang2D : public SwendsenWang {
public:
	SwendsenWang2D();

	virtual ~SwendsenWang2D();

	virtual void initialize();
	virtual void initializeClusterVariables();

	virtual void freezeOrMeltBonds();
	virtual void labelClusters();
	virtual void flipClusterSpins();

	virtual void measureObservables();

	virtual void free();

private:
	SwendsenWang2D(const SwendsenWang2D &);
	SwendsenWang2D &operator=(const SwendsenWang2D &);

private:
	int **s;                            // the spins
	bool **iBondFrozen;                 // bond lattice - three bonds per spin
	bool **jBondFrozen;
	double freezeProbability;           // 1 - e^(-2J/kT)
	int **cluster;                      // cluster labels for spins
};

class SwendsenWang3D : public SwendsenWang {
public:
	SwendsenWang3D();
	virtual ~SwendsenWang3D();

	virtual void initialize();
	virtual void initializeClusterVariables();

	virtual void freezeOrMeltBonds();
	virtual void labelClusters();
	virtual void flipClusterSpins();

	virtual void measureObservables();

	virtual void free();

private:
	SwendsenWang3D(const SwendsenWang3D &);
	SwendsenWang3D &operator=(const SwendsenWang3D &);

private:
	int ***s;                           // the spins
	bool ***iBondFrozen;                // bond lattice - three bonds per spin
	bool ***jBondFrozen;
	bool ***kBondFrozen;
	double freezeProbability;           // 1 - e^(-2J/kT)
	int ***cluster;                     // cluster labels for spins
};

#endif // SWENDSEN_WANG_H