#pragma once

#include "Trajectory.h"

class Spiral2D: public Trajectory
{
public:
    Spiral2D();
    ~Spiral2D();
    Spiral2D &operator=(Spiral2D &o);
    void Update(double dKtht1, double dKtht2, double dTht0, double dKmax, double dS, double dDt, double dOv);
private:
    double SovD2Tht(double dD0Tht, double dD1Tht, double dKtht1, double dKtht2, double dS);
};
