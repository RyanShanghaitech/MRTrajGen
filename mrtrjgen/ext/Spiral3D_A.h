#pragma once

#include "Trajectory.h"

class Spiral3D_A: public Trajectory
{
public:
    Spiral3D_A();
    ~Spiral3D_A();
	Spiral3D_A& operator=(Spiral3D_A &o);
    void Update(double dNp, double dUTht, double dUPhi, double dTht0, double dPhi0, double dKmax, double dS, double dDt);
private:
    double SovD1Phi(double dD0Phi, double dNp, double dUTht, double dUPhi, double dTht0, double dPhi0, double dS);
};