#include "Spiral3D.h"

Spiral3D::Spiral3D(): m_dPi(acos(-1e0)), m_dSqrt2(sqrt(2e0))
{
}

Spiral3D::~Spiral3D()
{
}

Spiral3D& Spiral3D::operator=(Spiral3D &o)
{
    m_vdKx = o.m_vdKx;
    m_vdKy = o.m_vdKy;
    m_vdKz = o.m_vdKz; // trajectory
    m_vdGx = o.m_vdGx;
    m_vdGy = o.m_vdGy;
    m_vdGz = o.m_vdGz; // gradient
    m_lNpt = o.m_lNpt;
    return *this;
}

double Spiral3D::SovQDE(double dA, double dB, double dC)
{
    double dDelta;

    dDelta = pow(dB, 2) - 4*dA*dC;
    dDelta = dDelta < 0 ? 0 : dDelta;

    return (-dB + sqrt(dDelta))/(2*dA);
}

double Spiral3D::SovD1Phi(double dS, double dNp, double dD0Phi, double dTht0, double dPhi0, double dUTht, double dUPhi)
{
    double dA, dB, dC;

    dA = (1.0/32.0)*dUPhi*(16*pow(dD0Phi, 23.0/2.0)*dUTht*(dD0Phi*dUTht*pow(sin(dPhi0 + dD0Phi), 2)*pow(sin(m_dSqrt2*sqrt(dD0Phi)*sqrt(dUPhi/dUTht) + dTht0), 2) + dD0Phi*dUTht*pow(sin(m_dSqrt2*sqrt(dD0Phi)*sqrt(dUPhi/dUTht) + dTht0), 2)*pow(cos(dPhi0 + dD0Phi), 2) + dUPhi*pow(sin(dPhi0 + dD0Phi), 2)*pow(sin(m_dSqrt2*sqrt(dD0Phi)*sqrt(dUPhi/dUTht) + dTht0), 2) + 2*dUPhi*pow(sin(dPhi0 + dD0Phi), 2)*pow(cos(m_dSqrt2*sqrt(dD0Phi)*sqrt(dUPhi/dUTht) + dTht0), 2) + dUPhi*pow(sin(m_dSqrt2*sqrt(dD0Phi)*sqrt(dUPhi/dUTht) + dTht0), 2)*pow(cos(dPhi0 + dD0Phi), 2) + 2*dUPhi*pow(cos(dPhi0 + dD0Phi), 2)*pow(cos(m_dSqrt2*sqrt(dD0Phi)*sqrt(dUPhi/dUTht) + dTht0), 2)) + 4*pow(dD0Phi, 21.0/2.0)*pow(dUPhi, 2) + 24*pow(dD0Phi, 21.0/2.0)*pow(dUTht, 2)*pow(sin(m_dSqrt2*sqrt(dD0Phi)*sqrt(dUPhi/dUTht) + dTht0), 2) + 6*pow(dD0Phi, 19.0/2.0)*dUPhi*dUTht + pow(dD0Phi, 17.0/2.0)*pow(dUTht, 2) + 12*m_dSqrt2*pow(dD0Phi, 11)*pow(dUTht, 2)*sqrt(dUPhi/dUTht)*sin(2*m_dSqrt2*sqrt(dD0Phi)*sqrt(dUPhi/dUTht) + 2*dTht0))/(pow(dNp, 2)*pow(dD0Phi, 23.0/2.0)*pow(m_dPi, 2)*dUTht);
    dB = 0;
    dC = -pow(dS, 2);

    return sqrt(fabs(SovQDE(dA, dB, dC)));
}

void Spiral3D::Update(double dS, double dNp, double dTht0, double dPhi0, double dUTht, double dUPhi, double dDt, double dKmax)
{
    std::list<double> ldKx, ldKy, ldKz;

    // k[0]
    double dRho = 0;
    double dTht = 0;
    double dPhi = 0;
    ldKx.push_back(dRho*sin(dTht + dTht0)*cos(dPhi + dPhi0));
    ldKy.push_back(dRho*sin(dTht + dTht0)*sin(dPhi + dPhi0));
    ldKz.push_back(dRho*cos(dTht + dTht0));

    // k[1], since SovD1Phi() doesn't allow phi=0, we have to calculate k[1]
    dRho = (0.01*dS)*dDt*dDt/2;
    dTht = dRho*2*m_dPi*dNp/dUTht;
    dPhi = (2*dRho*dRho*m_dPi*m_dPi*dNp*dNp)/(dUTht*dUPhi);
    ldKx.push_back(dRho*sin(dTht + dTht0)*cos(dPhi + dPhi0));
    ldKy.push_back(dRho*sin(dTht + dTht0)*sin(dPhi + dPhi0));
    ldKz.push_back(dRho*cos(dTht + dTht0));

    // k[2] ~ k[inf]
    int64_t lIdxPt = 2;
    double dD1Phi = 0;
    while (dRho <= dKmax)
    {
        dD1Phi = SovD1Phi(dS, dNp, dPhi, dTht0, dPhi0, dUTht, dUPhi);
        dPhi += dD1Phi*dDt;

        dTht = sqrt((2*dUPhi)/(dUTht))*sqrt(dPhi);
        dRho = sqrt(dUTht*dUPhi/2)/(m_dPi*dNp)*sqrt(dPhi);

        ldKx.push_back(dRho*sin(dTht + dTht0)*cos(dPhi + dPhi0));
        ldKy.push_back(dRho*sin(dTht + dTht0)*sin(dPhi + dPhi0));
        ldKz.push_back(dRho*cos(dTht + dTht0));

        ++lIdxPt;
    }
    m_lNpt = lIdxPt;
    
    // copy to member variables
    m_vdKx.resize(m_lNpt); m_vdKy.resize(m_lNpt); m_vdKz.resize(m_lNpt);
    m_vdKx.assign(ldKx.begin(), ldKx.end());
    m_vdKy.assign(ldKy.begin(), ldKy.end());
    m_vdKz.assign(ldKz.begin(), ldKz.end());

    // update gradient vector
    m_vdGx.resize(m_lNpt); m_vdGy.resize(m_lNpt); m_vdGz.resize(m_lNpt);
    m_vdGx[0] = 0e0; m_vdGy[0] = 0e0; m_vdGz[0] = 0e0;
    for (int64_t lIdxPt = 1; lIdxPt < m_lNpt; ++lIdxPt)
    {
        m_vdGx[lIdxPt] = (m_vdKx[lIdxPt] - m_vdKx[lIdxPt-1])/dDt;
        m_vdGy[lIdxPt] = (m_vdKy[lIdxPt] - m_vdKy[lIdxPt-1])/dDt;
        m_vdGz[lIdxPt] = (m_vdKz[lIdxPt] - m_vdKz[lIdxPt-1])/dDt;
    }
}

void Spiral3D::GetTraj(std::vector<double>* pvdKx, std::vector<double>* pvdKy, std::vector<double>* pvdKz, double dScale) const
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

void Spiral3D::GetGrad(std::vector<float>* pvfGx, std::vector<float>* pvfGy, std::vector<float>* pvfGz, double dScale) const
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

void Spiral3D::SaveGrad(FILE* pFile) const
{
    // open file
    if (pFile == NULL) return;

    // write gradient waveform into log
    for (int64_t idxPt = 0; idxPt < m_lNpt; ++idxPt)
    {
        fprintf(pFile, "%.8e, %.8e, %.8e\n", m_vdGx[idxPt], m_vdGy[idxPt], m_vdGz[idxPt]);
    }
}

void Spiral3D::AppendDown(std::vector<float> *pvfGrad, int64_t lNptDown)
{
    int64_t lNptRead = pvfGrad->size();
    int64_t lNptTotal = lNptRead + lNptDown;
    pvfGrad->resize(lNptTotal);

    double dGradMax = *pvfGrad->end();
    double dScale;
    for (int64_t iIdxGrad = lNptRead; iIdxGrad < lNptTotal; ++iIdxGrad)
    {
        dScale = (double)(lNptTotal - iIdxGrad - 1)/lNptDown;
        pvfGrad->at(iIdxGrad) = dGradMax*dScale;
    }
}

void Spiral3D::AppendZero(std::vector<float> *pvfGrad, int64_t lNptZero)
{
    int64_t lNptOld = pvfGrad->size();
    int64_t lNptNew = lNptOld + lNptZero;

    pvfGrad->resize(lNptNew);
    std::fill(pvfGrad->begin()+lNptOld, pvfGrad->end(), 0e0f);
}