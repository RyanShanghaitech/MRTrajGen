#pragma once

#include <cmath>
#include <vector>
#include <list>
#include <cstdio>

class Spiral3D
{
public:
    Spiral3D();
    ~Spiral3D();
	Spiral3D& operator=(Spiral3D &o);
    void Update(double dS, double dNp, double dTht0, double dPhi0, double dUTht, double dUPhi, double dDt=10e-6, double dKmax=0.5e0);
    void GetTraj(std::vector<double>* pvdKx, std::vector<double>* pvdKy, std::vector<double>* pvdKz, double dScale=1e0) const;
    void GetGrad(std::vector<float>* pvfGx, std::vector<float>* pvfGy, std::vector<float>* pvfGz, double dScale=1e0) const;
    int64_t GetNpt() const { return m_lNpt; }
    void SaveGrad(FILE* pFile) const;
private:
    std::vector<double> m_vdKx, m_vdKy, m_vdKz; // trajectory
    std::vector<double> m_vdGx, m_vdGy, m_vdGz; // gradient
    int64_t m_lNpt;

    const double m_dPi;
    const double m_dSqrt2;

    double SovD1Phi(double dS, double dNp, double dD0Phi, double dTht0, double dPhi0, double dUTht, double dUPhi);
    
    static double SovQDE(double dA, double dB, double dC);
    static void AppendDown(std::vector<float> *pvfGrad, int64_t lNptDown);
    static void AppendZero(std::vector<float> *pvfGrad, int64_t lNptZero);
};