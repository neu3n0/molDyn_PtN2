#pragma once
#include "Atom.hpp"
#include "Outer/Out.hpp"

const int MAX_CELLS = 30;
const double MASS_FOR_N = 2.33e-6;
const double MASS_FOR_PT = 3.24e-5;

class Cell {
public:
	~Cell();
	std::vector<Atom*> atoms{};
};

class Space {
public:
	void setConfig(std::istream& inp);
	void setConfig(const std::string& fname);
	void setConfig(const std::vector<double>&);
	void printConfig() const;

	void init(std::istream& inp);
	void init(const std::string& fname);
	void init(const std::vector<std::vector<double>>&);
	void printInit() const;

	bool initFromEnergy(std::istream& inp);
	bool initFromEnergy(const std::string& fname);
	bool initFromParams(const double E_tr, const double E_rot, const double E_vib, const double alpha);
	std::vector<double> initFromCoordsAndVel(
        const std::vector<double>& r1, 
            const std::vector<double>& r2, 
                const std::vector<double>& vAbs1, 
                    const std::vector<double>& vAbs2);
	bool prepareMolecule(const double E_tr, const double E_rot, const double E_vib, const double alpha);

	void printInfo() const;
	void printMol() const;

	std::vector<double> getStartMol();

public:
	Cell cells[MAX_CELLS][MAX_CELLS][MAX_CELLS];
	std::vector<MoleculeN2> molsN2{};
	bool isZ_periodic{false};

public:
	size_t countMolN2{0}, countAtPt{0};
	double lengthCell{0};
	double spaceLength[3]{0};  
	size_t numberCellsX{0}, numberCellsY{0}, numberCellsZ{0};
	double hMax{startHeight};

public:
	int MDStep();
	void saveAndStepsForAvg();
	void changeCell(size_t i, size_t j, size_t k, size_t iNew, size_t kNew, size_t jNew, int indAt);
	void SetNullMacro();
	void resetChecker();

public:
	int writeVTK(const std::string& name);
	inline static int vtkNum{0};

public:
	double temp{0}, p{0}, dt{0};
	double energy{0}, kinEn{0}, potEn{0}, vibEn{0};

public:
	double tempr();
	double rmsVel();
	std::vector<double> averVel();
	void getEnergy(Outer& out, int step);
	void saveStruct();
	double calcRotEn();
	double calcVibEn();
	double calcKinVib();

public:
	std::vector<double> calcShiftForMol(const MoleculeN2& mol);
	void calcEnergiesForN2();
	void saveAvgEn();
	double avgVibEn{0};
	double avgRotEn{0};
	double avgTrEn{0};
	size_t indAvg{0};
	inline static size_t period = 29;
	void calcAvgEn();
	double eTr{0};
	std::vector<double> vMol{};

public:
	static void saveInfo(Outer& out, Outer& out2, int i, const double E_tr, const double E_rot, const double E_vib,
		const double alpha, std::vector<double>& vel, const double eTr, const double eVib, 
			const double eRot, size_t n, double time, size_t steps, const std::vector<double>& startInf, bool flag);

	static bool getMoleculesParams(const std::string& filename, std::vector<std::vector<double>>& coord1, 
						std::vector<std::vector<double>>& coord2, 
						std::vector<std::vector<double>>& vel1, 
						std::vector<std::vector<double>>& vel2);

	static std::vector<double> getCfg(const std::string& fname);
	static std::vector<std::vector<double>> getPlatinumSurf(const std::string& fname);

public:
	void print_avel();
	void print_energy();
	void print_temp_wall();
};


