#pragma once

#include <vector>
#include <cstdio> // FILE
#include <cstdint> // int64_t
#include <cmath> // acos
#include <string>

#ifndef M_SQRT2
#define M_SQRT2		1.41421356237309504880e0
#endif
#ifndef M_PI
#define M_PI		3.14159265358979323846e0
#endif

class Trajectory
{
public:
    Trajectory();
    virtual ~Trajectory();
	virtual Trajectory& operator=(Trajectory &o);
    virtual void Update(std::vector<double> vdPara) {};
    void GetTraj(std::vector<double>* pvdKx, std::vector<double>* pvdKy, std::vector<double>* pvdKz, double dScale=1e0) const;
    void GetGrad(std::vector<float>* pvfGx, std::vector<float>* pvfGy, std::vector<float>* pvfGz, double dScale=1e0) const;
    int64_t GetNpt() const;
    void SaveTraj(FILE* pFile) const;
    void SaveGrad(FILE* pFile) const;
    static void AppendDown(std::vector<float> *pvfGrad, int64_t lNptDown);
    static void AppendZero(std::vector<float> *pvfGrad, int64_t lNptZero);
    std::vector<double> m_vdPara;
protected:
    std::vector<double> m_vdKx, m_vdKy, m_vdKz; // trajectory
    std::vector<double> m_vdGx, m_vdGy, m_vdGz; // gradient
    int64_t m_lNpt;

    static double SovQDE(double dA, double dB, double dC);
    static int64_t ExtractPara(std::vector<double> vdSrc, std::vector<double*> vpdDst);
};