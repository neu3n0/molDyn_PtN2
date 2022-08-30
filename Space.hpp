#pragma once
#include "Atom.hpp"

const int MAX_CELLS = 30;

class Cell {
public:
	~Cell();
	std::vector<Atom*> atoms{};
};

class Space {
public:
	size_t countMolN2{0}, countAtPt{0};
	double temp{0}, p{0}, dt{0};
	double lengthCell{0};
	double spaceLength[3]{0};  
	size_t numberCellsX{0}, numberCellsY{0}, numberCellsZ{0};
	double hMax{33.5};
	double energy{0}, kinEn{0}, potEn{0}, vibEn{0};
public:
	Cell cells[MAX_CELLS][MAX_CELLS][MAX_CELLS];
	std::vector<MoleculeN2> molsN2{};
	bool isZ_periodic{true};
};