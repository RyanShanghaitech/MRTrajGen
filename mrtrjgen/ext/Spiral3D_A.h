#pragma once

#include <cmath>
#include <vector>
#include <list>
#include "Trajectory.h"

class Spiral3D_A: public Trajectory
{
public:
    Spiral3D_A();
    ~Spiral3D_A();
	Spiral3D_A& operator=(Spiral3D_A &o);
    void Update(std::vector<double> vdPara);
private:
    double SovD1Phi(double dS, double dNp, double dD0Phi, double dTht0, double dPhi0, double dUTht, double dUPhi);
};