#include "Spiral3D_A.h"

Spiral3D_A::Spiral3D_A() {}

Spiral3D_A::~Spiral3D_A() {}

Spiral3D_A& Spiral3D_A::operator=(Spiral3D_A &o)
{
	return *this;
}

double Spiral3D_A::SovD1Phi(double dS, double dNp, double dD0Phi, double dTht0, double dPhi0, double dUTht, double dUPhi)
{
    double dA = (1.0/32.0)*dUPhi*(16*std::pow(dD0Phi, 23.0/2.0)*dUTht*(dD0Phi*dUTht*std::pow(std::sin(dD0Phi + dPhi0), 2)*std::pow(std::sin(M_SQRT2*std::sqrt(dD0Phi)*std::sqrt(dUPhi/dUTht) + dTht0), 2) + dD0Phi*dUTht*std::pow(std::sin(M_SQRT2*std::sqrt(dD0Phi)*std::sqrt(dUPhi/dUTht) + dTht0), 2)*std::pow(std::cos(dD0Phi + dPhi0), 2) + dUPhi*std::pow(std::sin(dD0Phi + dPhi0), 2)*std::pow(std::sin(M_SQRT2*std::sqrt(dD0Phi)*std::sqrt(dUPhi/dUTht) + dTht0), 2) + 2*dUPhi*std::pow(std::sin(dD0Phi + dPhi0), 2)*std::pow(std::cos(M_SQRT2*std::sqrt(dD0Phi)*std::sqrt(dUPhi/dUTht) + dTht0), 2) + dUPhi*std::pow(std::sin(M_SQRT2*std::sqrt(dD0Phi)*std::sqrt(dUPhi/dUTht) + dTht0), 2)*std::pow(std::cos(dD0Phi + dPhi0), 2) + 2*dUPhi*std::pow(std::cos(dD0Phi + dPhi0), 2)*std::pow(std::cos(M_SQRT2*std::sqrt(dD0Phi)*std::sqrt(dUPhi/dUTht) + dTht0), 2)) + 4*std::pow(dD0Phi, 21.0/2.0)*std::pow(dUPhi, 2) + 24*std::pow(dD0Phi, 21.0/2.0)*std::pow(dUTht, 2)*std::pow(std::sin(M_SQRT2*std::sqrt(dD0Phi)*std::sqrt(dUPhi/dUTht) + dTht0), 2) + 6*std::pow(dD0Phi, 19.0/2.0)*dUPhi*dUTht + std::pow(dD0Phi, 17.0/2.0)*std::pow(dUTht, 2) + 12*M_SQRT2*std::pow(dD0Phi, 11)*std::pow(dUTht, 2)*std::sqrt(dUPhi/dUTht)*std::sin(2*M_SQRT2*std::sqrt(dD0Phi)*std::sqrt(dUPhi/dUTht) + 2*dTht0))/(std::pow(dD0Phi, 23.0/2.0)*std::pow(dNp, 2)*dUTht*std::pow(M_PI, 2));
    double dB = 0;
    double dC = -std::pow(dS, 2);

    return sqrt(fabs(SovQuadEq(dA, dB, dC)));
}

void Spiral3D_A::Update(std::vector<double> vdPara)
{
    m_vdPara = vdPara;
    double dS, dNp, dTht0, dPhi0, dUTht, dUPhi, dDt, dKmax;
    double *adDst[] = {&dS, &dNp, &dTht0, &dPhi0, &dUTht, &dUPhi, &dDt, &dKmax};
    ExtractPara(vdPara, std::vector<double*>(adDst, adDst + sizeof(adDst)/sizeof(double*)));
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
    dTht = dRho*2*M_PI*dNp/dUTht;
    dPhi = (2*dRho*dRho*M_PI*M_PI*dNp*dNp)/(dUTht*dUPhi);
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
        dRho = sqrt(dUTht*dUPhi/2)/(M_PI*dNp)*sqrt(dPhi);

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
