#pragma once

#include "Trajectory.h"

class Spiral2D_Inv: public Trajectory
{
public:
    Spiral2D_Inv();
    ~Spiral2D_Inv();
    Spiral2D_Inv &operator=(Spiral2D_Inv &o);
    void Update(std::vector<double> *pvdPara);
private:
    double SovD2Tht(double dD0Tht, double dD1Tht, double dNp, double dU, double dS, double dSign);
    double SovD1Tht(double dD0Tht, double dNp, double dU, double dS);
};
