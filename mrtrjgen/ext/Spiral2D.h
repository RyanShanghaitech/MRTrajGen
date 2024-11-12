#pragma once

#include "Trajectory.h"

class Spiral2D: public Trajectory
{
public:
    Spiral2D();
    ~Spiral2D();
    Spiral2D &operator=(Spiral2D &o);
    void Update(std::vector<double> *pvdPara);
private:
    double SovD2Tht(double dD0Tht, double dD1Tht, double dNp, double dU, double dS);
};
