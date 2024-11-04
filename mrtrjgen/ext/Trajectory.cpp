#include "Trajectory.h"

Trajectory::Trajectory() {}

Trajectory::~Trajectory() {}

Trajectory& Trajectory::operator=(Trajectory &o)
{
	return *this;
}

void Trajectory::GetTraj(std::vector<double>* pvdKx, std::vector<double>* pvdKy, std::vector<double>* pvdKz, double dScale) const
{
    if (pvdKx->size() < (size_t)m_lNpt) pvdKx->resize(m_lNpt);
    if (pvdKy->size() < (size_t)m_lNpt) pvdKy->resize(m_lNpt);
    if (pvdKz->size() < (size_t)m_lNpt) pvdKz->resize(m_lNpt);
    for (int64_t lIdxPt = 0; lIdxPt < m_lNpt; ++lIdxPt)
    {
        pvdKx->at(lIdxPt) = m_vdKx[lIdxPt]*dScale;
        pvdKy->at(lIdxPt) = m_vdKy[lIdxPt]*dScale;
        pvdKz->at(lIdxPt) = m_vdKz[lIdxPt]*dScale;
    }
}

void Trajectory::GetGrad(std::vector<float>* pvfGx, std::vector<float>* pvfGy, std::vector<float>* pvfGz, double dScale) const
{
    if (pvfGx->size() < (size_t)m_lNpt) pvfGx->resize(m_lNpt);
    if (pvfGy->size() < (size_t)m_lNpt) pvfGy->resize(m_lNpt);
    if (pvfGz->size() < (size_t)m_lNpt) pvfGz->resize(m_lNpt);
    for (int64_t lIdxPt = 0; lIdxPt < m_lNpt; ++lIdxPt)
    {
        pvfGx->at(lIdxPt) = m_vdGx[lIdxPt]*dScale;
        pvfGy->at(lIdxPt) = m_vdGy[lIdxPt]*dScale;
        pvfGz->at(lIdxPt) = m_vdGz[lIdxPt]*dScale;
    }
}

int64_t Trajectory::GetNpt() const
{
    return m_lNpt;
}

void Trajectory::SaveTraj(FILE* pFile) const
{
    if (pFile == NULL) return;
    for (int64_t lIdxPt = 0; lIdxPt < m_lNpt; ++lIdxPt)
    {
        fprintf(pFile, "%.8e, %.8e, %.8e\n", m_vdKx[lIdxPt], m_vdKy[lIdxPt], m_vdKz[lIdxPt]);
    }
}

void Trajectory::SaveGrad(FILE* pFile) const
{
    if (pFile == NULL) return;
    for (int64_t lIdxPt = 0; lIdxPt < m_lNpt; ++lIdxPt)
    {
        fprintf(pFile, "%.8e, %.8e, %.8e\n", m_vdGx[lIdxPt], m_vdGy[lIdxPt], m_vdGz[lIdxPt]);
    }
}

void Trajectory::AppendDown(std::vector<float> *pvfGrad, int64_t lNptDown)
{
    int64_t lNptOld = pvfGrad->size();
    int64_t lNptNew = lNptOld + lNptDown;
    pvfGrad->resize(lNptNew);

    double dGradMax = pvfGrad->at(lNptOld - 1);
    double dScale;
    for (int64_t lIdxGrad = lNptOld; lIdxGrad < lNptNew; ++lIdxGrad)
    {
        dScale = (double)(lNptNew - lIdxGrad - 1)/lNptDown;
        pvfGrad->at(lIdxGrad) = dGradMax*dScale;
    }
}

void Trajectory::AppendZero(std::vector<float> *pvfGrad, int64_t lNptZero)
{
    int64_t lNptOld = pvfGrad->size();
    int64_t lNptNew = lNptOld + lNptZero;

    pvfGrad->resize(lNptNew);
    std::fill(pvfGrad->begin()+lNptOld, pvfGrad->end(), 0e0f);
}

double Trajectory::SovQDE(double dA, double dB, double dC)
{
    double dDelta;

    dDelta = pow(dB, 2) - 4*dA*dC;
    dDelta = dDelta < 0 ? 0 : dDelta;

    return (-dB + sqrt(dDelta))/(2*dA);
}

int64_t Trajectory::ExtractPara(std::vector<double> vdSrc, std::vector<double*> vpdDst)
{
    if (vdSrc.size() != vpdDst.size()) return 1;

    int64_t lNPara = vpdDst.size();
    for (int64_t lIdxPara = 0; lIdxPara < lNPara; ++lIdxPara)
    {
        *vpdDst[lIdxPara] = vdSrc[lIdxPara];
    }

    return 0;
}