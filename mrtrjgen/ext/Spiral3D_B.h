#pragma once

#include "Trajectory.h"

class Spiral3D_B: public Trajectory
{
public:
    Spiral3D_B();
    ~Spiral3D_B();
	Spiral3D_B& operator=(Spiral3D_B &o);
    void Update(std::vector<double> *pvdPara);
private:
    double SovD1Tht(double dS, double dNp, double dD0Tht, double dTht0, double dPhi0, double dUTht, double dUPhi);
};