#include "Trajectory.h"

Trajectory::Trajectory() {}

Trajectory::~Trajectory() {}

Trajectory& Trajectory::operator=(Trajectory &o)
{
	return *this;
}

void Trajectory::GetTraj(std::vector<double>* pvdKx, std::vector<double>* pvdKy, std::vector<double>* pvdKz, double dScale) const
{
    if (pvdKx != NULL && pvdKx->size() < (size_t)m_lNpt) pvdKx->resize(m_lNpt);
    if (pvdKy != NULL && pvdKy->size() < (size_t)m_lNpt) pvdKy->resize(m_lNpt);
    if (pvdKz != NULL && pvdKz->size() < (size_t)m_lNpt) pvdKz->resize(m_lNpt);
    for (int64_t lIdxPt = 0; lIdxPt < m_lNpt; ++lIdxPt)
    {
        if (pvdKx != NULL) pvdKx->at(lIdxPt) = m_vdKx[lIdxPt]*dScale;
        if (pvdKy != NULL) pvdKy->at(lIdxPt) = m_vdKy[lIdxPt]*dScale;
        if (pvdKz != NULL) pvdKz->at(lIdxPt) = m_vdKz[lIdxPt]*dScale;
    }
}

void Trajectory::GetGrad(std::vector<float>* pvfGx, std::vector<float>* pvfGy, std::vector<float>* pvfGz, double dScale) const
{
    if (pvfGx != NULL && pvfGx->size() < (size_t)m_lNpt) pvfGx->resize(m_lNpt);
    if (pvfGy != NULL && pvfGy->size() < (size_t)m_lNpt) pvfGy->resize(m_lNpt);
    if (pvfGz != NULL && pvfGz->size() < (size_t)m_lNpt) pvfGz->resize(m_lNpt);
    for (int64_t lIdxPt = 0; lIdxPt < m_lNpt; ++lIdxPt)
    {
        if (pvfGx != NULL) pvfGx->at(lIdxPt) = m_vdGx[lIdxPt]*dScale;
        if (pvfGy != NULL) pvfGy->at(lIdxPt) = m_vdGy[lIdxPt]*dScale;
        if (pvfGz != NULL) pvfGz->at(lIdxPt) = m_vdGz[lIdxPt]*dScale;
    }
}

int64_t Trajectory::GetNpt() const
{
    return m_lNpt;
}

void Trajectory::SaveTraj(FILE* pfBin, FILE* pfHdr) const
{
    // save binary file
    if (pfBin != NULL)
    {
        fwrite(m_vdKx.data(), sizeof(double), m_vdKx.size(), pfBin);
        fwrite(m_vdKy.data(), sizeof(double), m_vdKy.size(), pfBin);
        fwrite(m_vdKz.data(), sizeof(double), m_vdKz.size(), pfBin);
    }
    // save header file
    if (pfHdr != NULL)
    {
        fprintf(pfHdr, "double m_vdKx[%ld]\n", m_vdKx.size());
        fprintf(pfHdr, "double m_vdKy[%ld]\n", m_vdKy.size());
        fprintf(pfHdr, "double m_vdKz[%ld]\n", m_vdKz.size());
    }
}

void Trajectory::SaveGrad(FILE* pfBin, FILE* pfHdr) const
{
    // save binary file
    if (pfBin != NULL)
    {
        fwrite(m_vdGx.data(), sizeof(double), m_vdGx.size(), pfBin);
        fwrite(m_vdGy.data(), sizeof(double), m_vdGy.size(), pfBin);
        fwrite(m_vdGz.data(), sizeof(double), m_vdGz.size(), pfBin);
    }
    // save header file
    if (pfHdr != NULL)
    {
        fprintf(pfHdr, "double m_vdGx[%ld]\n", m_vdGx.size());
        fprintf(pfHdr, "double m_vdGy[%ld]\n", m_vdGy.size());
        fprintf(pfHdr, "double m_vdGz[%ld]\n", m_vdGz.size());
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

double Trajectory::SovQuadEq(double dA, double dB, double dC, double dSign)
{
    double dDelta;

    dDelta = pow(dB, 2) - 4*dA*dC;
    dDelta = dDelta < 0 ? 0 : dDelta;

    return (-dB + dSign*sqrt(dDelta))/(2*dA);
}

void Trajectory::ExtractPara(std::vector<double> *pvdSrc, std::vector<double*> *pvpdDst)
{
    if (pvdSrc->size() != pvpdDst->size()) throw std::runtime_error("");;

    int64_t lNPara = pvpdDst->size();
    for (int64_t lIdxPara = 0; lIdxPara < lNPara; ++lIdxPara)
    {
        *pvpdDst->at(lIdxPara) = pvdSrc->at(lIdxPara);
    }
}