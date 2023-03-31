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
	void init(std::istream& inp);
	void init(const std::string& fname);
	bool initFromEnergy(std::istream& inp);
	bool initFromEnergy(const std::string& fname);
	bool initFromParams(const double E_tr, const double E_rot, const double E_vib, const double alpha);
	std::vector<double> initFromCoordsAndVel(
        const std::vector<double>& r1, 
            const std::vector<double>& r2, 
                const std::vector<double>& vAbs1, 
                    const std::vector<double>& vAbs2);
public:
	void printConfig() const;
	void printInit() const;
	void printInfo() const;
	void printMol() const;
public:
	size_t countMolN2{0}, countAtPt{0};
	double temp{0}, p{0}, dt{0};
	// double spaceLength{0}, lengthCell{0};
	double lengthCell{0};
	double spaceLength[3]{0};  
	// size_t numberCells{0};
	size_t numberCellsX{0}, numberCellsY{0}, numberCellsZ{0};
	double hMax{startHeight};
	double energy{0}, kinEn{0}, potEn{0}, vibEn{0};
public:
	Cell cells[MAX_CELLS][MAX_CELLS][MAX_CELLS];
	std::vector<MoleculeN2> molsN2{};
	bool isZ_periodic{false};
public:
	int MDStep();
	void changeCell(size_t i, size_t j, size_t k, size_t iNew, size_t kNew, size_t jNew, int indAt);
	void SetNullMacro();
	void resetChecker();
public:
	int writeVTK(const std::string& name);
	inline static int vtkNum{0};
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
	bool saveAvg{false};
	void saveAvgEn();
	double avgVibEn{0};
	double avgRotEn{0};
	size_t indAvg{0};

public:
	std::vector<double> calcShiftForMol(const MoleculeN2& mol);
};

void saveInfo(Outer& out, const double E_tr, const double E_rot, const double E_vib,
    const double alpha, std::vector<double>& vel, const double eTr, const double eVib, 
		const double eRot, size_t n, double time, size_t steps);