#include "Spiral3D_B.h"

Spiral3D_B::Spiral3D_B() {}

Spiral3D_B::~Spiral3D_B() {}

Spiral3D_B& Spiral3D_B::operator=(Spiral3D_B &o)
{
	return *this;
}

double Spiral3D_B::SovD1Tht(double dS, double dNp, double dD0Tht, double dTht0, double dPhi0, double dUTht, double dUPhi)
{
    double dA = (1.0/32.0)*dUTht*(8*std::pow(dD0Tht, 3)*dUPhi*(2*dD0Tht*dUPhi + dUTht*std::cos(2*dD0Tht + 2*dTht0) + 3*dUTht) + 24*std::pow(dD0Tht, 2)*std::pow(dUPhi, 2) + 4*std::pow(dD0Tht, 2)*std::pow(dUTht, 2)*std::pow(std::sin(dD0Tht + dTht0), 2) + 6*dD0Tht*dUPhi*dUTht*std::pow(std::sin(dD0Tht + dTht0), 2) + std::pow(dUPhi, 2))/(std::pow(dD0Tht, 3)*std::pow(dNp, 2)*dUPhi*std::pow(M_PI, 2));
    double dB = 0;
    double dC = -std::pow(dS, 2);

    return sqrt(fabs(SovQuadEq(dA, dB, dC)));
}

void Spiral3D_B::Update(std::vector<double> vdPara)
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
    dPhi = (2*M_PI*dNp/dUPhi)*dRho;
    dTht = dUPhi/(2*dUTht)*(dPhi*dPhi);
    ldKx.push_back(dRho*sin(dTht + dTht0)*cos(dPhi + dPhi0));
    ldKy.push_back(dRho*sin(dTht + dTht0)*sin(dPhi + dPhi0));
    ldKz.push_back(dRho*cos(dTht + dTht0));

    // k[2] ~ k[inf]
    int64_t lIdxPt = 2;
    double dD1Tht = 0;
    while (dRho <= dKmax)
    {
        dD1Tht = SovD1Tht(dS, dNp, dTht, dTht0, dPhi0, dUTht, dUPhi);
        dTht += dD1Tht*dDt;

        dPhi = sqrt((2*dUTht)/(dUPhi))*sqrt(dTht);
        dRho = dUPhi/(2*M_PI*dNp)*dPhi;

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
